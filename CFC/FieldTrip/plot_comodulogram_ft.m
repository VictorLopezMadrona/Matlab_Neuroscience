function [fig_cfc,comodulogram]=plot_comodulogram_ft(comodulogram,cnfg)

%% Plot comodulogram with some smooth
%
% USE:
%   plot_comodulogram_ft(comodulogram,cfg);
%
% INPUT:
%   comodulogram = output of computeCFC
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
% See also: computeCFC comodulogram_ft computeCFC_trial comodulogram_trial

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Apr. 2016; Last revision: 16-Aug-2023

% Change log:
% 17/10/2023: Changed to contour plot
% 16/08/2023: Added output mask with stats in comodulogram
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

if isempty(comodulogram.MI_pval)
    warning('No surrogates were computed to do stats on the comodulogram')
    cnfg.stats='none';
end
    
pval = cnfg.pval;
pval_clus = cnfg.pval_clus;

%Extract the information from the struct
f_min_theta = comodulogram.f_theta.f_min;
f_max_theta = comodulogram.f_theta.f_max;
step_theta  = comodulogram.f_theta.step;
f_min_gamma = comodulogram.f_gamma.f_min;
f_max_gamma = comodulogram.f_gamma.f_max;
CFC         = comodulogram.CFC;
MI          = comodulogram.MI;

%% STATISTICS

%Different options of pval can be applyed here:
%%% MEAN
if strcmp(cnfg.stats,'mean')
    MI=comodulogram.MI - mean(comodulogram.MI_pval,3);
    MI(MI<0)=0;
end


%%% SINGLE PIXEL STATS
if strcmp(cnfg.stats,'single') || strcmp(cnfg.stats,'cluster')
    pval=1-pval;
    MI=comodulogram.MI;
    MIpval=MI*0;
    for i=1:size(comodulogram.MI,1)
        for j=1:size(comodulogram.MI,2)
            MU=mean(comodulogram.MI_pval(i,j,:));
            SIGMA=std(comodulogram.MI_pval(i,j,:));
            MIpval(i,j)=norminv(pval,MU,SIGMA);
            if MIpval(i,j)>MI(i,j)
                MI(i,j)=0;
            end
        end
    end
end

%%% Cluster-based statistics
% Threshold is based on the sum of t-stats in each cluster
if strcmp(cnfg.stats,'cluster')
    pval_clus = 1-pval_clus;
    MImu = mean(comodulogram.MI_pval,3);
    MIsigma = std(comodulogram.MI_pval,[],3);
    
    % Cluster surrogate
    surro_cluster = zeros(1,size(comodulogram.MI_pval,3));
    for s=1:size(comodulogram.MI_pval,3)
        %z-values (t-stats)
        MIz_surro    = (comodulogram.MI_pval(:,:,s) - MImu) ./ MIsigma;
        MIsig_surro  = comodulogram.MI_pval(:,:,s)>=MIpval;
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
    MIz = (comodulogram.MI - MImu) ./ MIsigma; %z-values (t-stats)
    MIsig  = comodulogram.MI>=MIpval;
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
    MIz_surro = comodulogram.MI_pval*0;
    MImu = mean(comodulogram.MI_pval,3);
    MIsigma = std(comodulogram.MI_pval,[],3);
    for s=1:size(comodulogram.MI_pval,3)
        MIz_surro(:,:,s) = (comodulogram.MI_pval(:,:,s) - MImu) ./ MIsigma;
    end
    % Cluster Original
    MIz = (comodulogram.MI - MImu) ./ MIsigma; %z-values (t-stats)
    
    pixel_surro = squeeze(max(max(MIz_surro,[],1),[],2));
    MU_pixel    = mean(pixel_surro);
    SIGMA_pixel = std(pixel_surro);
    TH_pixel    = norminv(pval,MU_pixel,SIGMA_pixel);    
    % Apply threshold
    MI   = comodulogram.MI;
    MIth = (MIz>=TH_pixel);
    MI   = MI.*MIth;
end

mask = MI*0;
mask(MI>0) = 1;
comodulogram.mask = mask;

%% MI plot
[ly,~]=size(MI);
x=[f_min_theta f_max_theta];
y=[f_min_gamma f_max_gamma];

MI(MI<0)=0; %Remove negative values
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
set(get(h,'title'),'string','Modulation Index');

%% Phase distribution plot
% 
% [~,~,bins]=size(CFC);
% xx=f_min_theta:step_theta:f_max_theta;
% [~,x_phase]=max(xx>=phase);
% 
% CFC_plot=zeros(ly,bins);
% CFC_plot(:,:)=CFC(:,x_phase,:);
%     
% CFC_180(:,1:bins/2)=CFC_plot(:,bins/2+1:end);
% CFC_180(:,bins/2+1:bins)=CFC_plot(:,1:bins/2);
% CFC_180(:,bins+1:bins*2)=CFC_180;
% 
% %Make a z-score
% for i=1:ly
%     CFC_180(i,:)=(CFC_180(i,:)-mean(CFC_180(i,:)));
% end
% 
% figure,
% imagesc_filter(CFC_180,4,[0 720],y);
% 
% %Seno
% x_sin=0:0.01:4*pi;
% x_frec_sin=(1:1257)*720/1257;
% A_sin=(max(y)-min(y))/8;
% V=max(y)-A_sin;
% 
% colormap('jet')
% h=colorbar;
% set(get(h,'title'),'string','Norm. amplitude');
% hold on
% plot(x_frec_sin,V+A_sin*sin(x_sin+pi/2),'k')
% xlabel('Phase (degrees)')
% ylabel('Amplitude frequency (Hz)')

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
