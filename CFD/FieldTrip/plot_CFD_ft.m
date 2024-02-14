function [fig_cfc,comodulogram]=plot_CFD_ft(comodulogram,cnfg)

%% Plot CFD with some smooth
%
% USE:
%   plot_CFD_ft(comodulogram,cfg);
%
% INPUT:
%   comodulogram = output of CFD_parallel
%   cfg.stats = 'none'    - Do not compute stats
%               'mean'    - Substract the mean of the surrogates
%               'single'  - Pval of each value with surrogates. No multiple comparisons
%               'cluster' - [Default] Correct pvalues with clustering
%               'pixel'   - Pixel-based correction for multiple comparisons
%
%   cfg.pval       - [Def = 0.01] pval for single point correction
%   cfg.pval_clus  - [Def = 0.01] pval for cluster correction
%
% OUTPUT:
%
%
% See also: computeCFC comodulogram_ft computeCFC_trial CFD_parallel

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Apr. 2016; Last revision: 21-Feb-2023

% Change log:
% 17/10/2023: Changed to contour plot
% 16/08/2023: Modified to store statistical mask too
% 21/02/2023: Corrected bug when ploting zero values
% 21/02/2023: Corrected bug with negative pval in pixel-based stats
% 15/02/2023: Clims of plots between [-max +max]. 
% 15/02/2023: Names modified to refer to phase and amplitude. 
% 17/08/2022: Can be applied to CFD
% 25/08/2021: Cluster and pixel-based stats included

%% PARAMETERS

if nargin == 1
    cnfg.stats     = 'cluster';
    cnfg.pval      = 0.01;
    cnfg.pval_clus = 0.01;
end

if ~isfield(cnfg,'stats'), cnfg.stats='cluster'; end
if ~isfield(cnfg,'pval'), cnfg.pval=0.01; end
if ~isfield(cnfg,'pval_clus'), cnfg.pval_clus=0.01; end

if isempty(comodulogram.pval)
    warning('No surrogates were computed to do stats on the comodulogram')
    cnfg.stats='none';
end
    
pval = cnfg.pval;
pval_clus = cnfg.pval_clus;

%Extract the information from the struct
f_min_phase = comodulogram.f_phase.f_min;
f_max_phase = comodulogram.f_phase.f_max;
step_phase  = comodulogram.f_phase.step;
f_min_amp = comodulogram.f_amp.f_min;
f_max_amp = comodulogram.f_amp.f_max;
MI          = comodulogram.PSI;
MI_pval     = comodulogram.pval;

%% STATISTICS

%Different options of pval can be applyed here:
%%% MEAN
if strcmp(cnfg.stats,'mean')
    MI=MI - mean(MI_pval,3);
end


%%% SINGLE PIXEL STATS
if strcmp(cnfg.stats,'single') || strcmp(cnfg.stats,'cluster')
    pval=1-pval;
    MIpval=MI*0;
    MIpval_op=MI*0;
    for i=1:size(comodulogram.PSI,1)
        for j=1:size(comodulogram.PSI,2)
            MU=mean(comodulogram.pval(i,j,:));
            SIGMA=std(comodulogram.pval(i,j,:));
            MIpval(i,j)=norminv(pval,MU,SIGMA);
            MIpval_op(i,j)=norminv(1-pval,MU,SIGMA);
            if MIpval(i,j)>MI(i,j) && MIpval_op(i,j)<MI(i,j)
                MI(i,j)=0;
            end
        end
    end
end

%%% Cluster-based statistics
% Threshold is based on the sum of t-stats in each cluster
if strcmp(cnfg.stats,'cluster')
    pval_clus = 1-pval_clus;
    MImu = mean(comodulogram.pval,3);
    MIsigma = std(comodulogram.pval,[],3);
    
    % Cluster surrogate
    surro_cluster = zeros(1,size(comodulogram.pval,3));
    for s=1:size(comodulogram.pval,3)
        %z-values (t-stats)
        MIz_surro    = (comodulogram.pval(:,:,s) - MImu) ./ MIsigma;
        MIsig_surro  = comodulogram.pval(:,:,s)>=MIpval + comodulogram.pval(:,:,s)<=MIpval_op;
        MIsigz_surro = MIz_surro.*MIsig_surro;
        
        CC = bwconncomp(MIsigz_surro);
        SS = regionprops(CC,'Area','PixelIdxList');
        
        if ~isempty(SS)
            cluster_aux = zeros(1,length(SS));
            for nreg = 1:length(SS)
                cluster_aux(nreg) = sum(MIz_surro(SS(nreg).PixelIdxList));
            end
            surro_cluster(s) = max(cluster_aux);
        end
    end
    
    % Cluster Original
    MIz = (comodulogram.PSI - MImu) ./ MIsigma; %z-values (t-stats)
    MIsig  = comodulogram.PSI>=MIpval + comodulogram.PSI<=MIpval_op;
    MIsigz = MIz.*MIsig;
    
    % Put those nonsignificant clusters to zero
    MU_clus    = mean(surro_cluster);
    SIGMA_clus = std(surro_cluster);
    TH_clus    = norminv(pval_clus,MU_clus,SIGMA_clus);
    CC = bwconncomp(MIsigz);
    SS = regionprops(CC,'Area','PixelIdxList');
    
    if ~isempty(SS)
        for nreg = 1:length(SS)
            cluster_tval = sum(MIz(SS(nreg).PixelIdxList));
            if cluster_tval < TH_clus
                MI(SS(nreg).PixelIdxList) = 0;
            end
        end
    end
end


%%% PIXEL-BASED STATS
if strcmp(cnfg.stats,'pixel')
    pval=1-pval;
    %z-values (t-stats) surrogates
    MIz_surro = comodulogram.pval*0;
    MImu = mean(comodulogram.pval,3);
    MIsigma = std(comodulogram.pval,[],3);
    for s=1:size(comodulogram.pval,3)
        MIz_surro(:,:,s) = (comodulogram.pval(:,:,s) - MImu) ./ MIsigma;
    end
    % Cluster Original
    MIz = (comodulogram.PSI - MImu) ./ MIsigma; %z-values (t-stats)
    
    pixel_surro = squeeze(max(max(MIz_surro,[],1),[],2));
    MU_pixel    = mean(pixel_surro);
    SIGMA_pixel = std(pixel_surro);
    TH_pixel    = norminv(pval,MU_pixel,SIGMA_pixel); 
    pixel_surro_op = squeeze(min(min(MIz_surro,[],1),[],2));
    MU_pixel_op    = mean(pixel_surro_op);
    SIGMA_pixel_op = std(pixel_surro_op);
    TH_pixel_op    = norminv(1-pval,MU_pixel_op,SIGMA_pixel_op); 
    
    % Apply threshold
    MI   = comodulogram.PSI;
    MIth = (MIz>=TH_pixel) + (MIz<=TH_pixel_op);
    MI   = MI.*MIth;
end

mask = MI*0;
mask(MI~=0) = 1;
comodulogram.mask = mask;

%% MI plot
[ly,~]=size(MI);
x=[f_min_phase f_max_phase];
y=[f_min_amp f_max_amp];

%figure,
%fig_cfc=imagesc_filter(MI,4,x,y);
%fig_cfc=imagesc(x,y,MI); axis('xy');
xx=linspace(x(1),x(2),size(MI,2));
yy=linspace(y(1),y(2),size(MI,1));
[~,fig_cfc]=contourf(xx,yy,MI); axis('xy');

xlabel('Phase frequency (Hz)'),
ylabel('Amplitude frequency (Hz)'),
colormap('jet')
h = colorbar;
set(get(h,'title'),'string','Phase Slope Index');
if isprop(fig_cfc,'ZData')
    maxval = max(abs(fig_cfc.ZData(:)));
else
    maxval = max(abs(fig_cfc.CData(:)));
end
if maxval>0, clim([-maxval maxval]); end

end

function fig_cfc=imagesc_filter(im_in,smooth,x_axis,y_axis)

%
% Plot a smooth image increasing the number of pixels and filtering with a 
% low-pass.
%
% USE:
%    imagesc_filter(im_in,smooth,x,y)
%
% INPUT:
%    im_in: matrix (MxN) with the image to plot.
%    smooth: Integer with the intensity of the smooth (Default = 0 -no smooth).
%    x,y (Optional): picture axis with format x=[x_min x_max].

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Jun. 2016; Last revision: 14-Jul-2020

if nargin ==1
    n = 0;
elseif nargin==2
    n = smooth;
    [X,Y]=size(im_in);
    y_axis=[0 X];
    x_axis=[0 Y];
elseif nargin>2
    n = smooth;
end

if n==0
    imagesc(im_in)
    return
end

h=1/9*[1 1 1; 1 1 1; 1 1 1];
M_pre=max(max(im_in));
im_in=imfilter(imfilter(im_in,h),h);
M_post=max(max(im_in));
im_in=im_in*M_pre/M_post;

[X,Y]=size(im_in);
im_out=zeros(X*n,Y*n);

for x=1:X
    for y=1:Y
        im_out((x-1)*n+1:x*n, (y-1)*n+1:y*n) = im_in(x,y);
    end
end

IM_f=imfilter(imfilter(im_out,h),h);

xx=linspace(x_axis(1),x_axis(2),X*4);
yy=linspace(y_axis(1),y_axis(2),Y*4);
fig_cfc=imagesc(xx,yy,IM_f); axis('xy');

end
