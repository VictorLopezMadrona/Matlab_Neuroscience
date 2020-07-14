function plot_comodulogram(comodulogram,phase)

%% Plot comodulogram and phase distribution with some smooth
%
% USE:
%   plot_comodulogram(comodulogram,phase);
%
% INPUT:
%   comodulogram = output of comodulogram_parallel
%   phase [Optional] = freq. of the main phase (Def = 8 Hz)
%
% OUTPUT:
%
%
% See also: comodulogram_parallel

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Apr. 2016; Last revision: 14-Jul-2020

if nargin == 1
    phase = 8;
end

%Extract the information from the struct
f_min_theta=comodulogram.f_theta.f_min;
f_max_theta=comodulogram.f_theta.f_max;
step_theta=comodulogram.f_theta.step;
f_min_gamma=comodulogram.f_gamma.f_min;
f_max_gamma=comodulogram.f_gamma.f_max;

%Different options of pval can be applyed here:
MI=comodulogram.MI - mean(comodulogram.MI_pval,3);
CFC=comodulogram.CFC;

[ly,lx]=size(MI);
x=[f_min_theta f_max_theta];
y=[f_min_gamma f_max_gamma];

%% MI plot
MI(MI<0)=0; %Remove negative values
figure,
imagesc_filter(MI,4,x,y);

xlabel('Phase frequency (Hz)'),
ylabel('Amplitude frequency (Hz)'),
colormap('jet')
h = colorbar;
set(get(h,'title'),'string','Modulation Index');

%% Phase distribution plot

[~,~,bins]=size(CFC);
xx=f_min_theta:step_theta:f_max_theta;
[~,x_phase]=max(xx>=phase);

CFC_plot=zeros(ly,bins);
CFC_plot(:,:)=CFC(:,x_phase,:);
    
CFC_180(:,1:bins/2)=CFC_plot(:,bins/2+1:end);
CFC_180(:,bins/2+1:bins)=CFC_plot(:,1:bins/2);
CFC_180(:,bins+1:bins*2)=CFC_180;

%Make a z-score
for i=1:ly
    CFC_180(i,:)=(CFC_180(i,:)-mean(CFC_180(i,:)));
end

figure,
imagesc_filter(CFC_180,4,[0 720],y);

%Seno
x_sin=0:0.01:4*pi;
x_frec_sin=(1:1257)*720/1257;
A_sin=(max(y)-min(y))/8;
V=max(y)-A_sin;

colormap('jet')
h=colorbar;
set(get(h,'title'),'string','Norm. amplitude');
hold on
plot(x_frec_sin,V+A_sin*sin(x_sin+pi/2),'k')
xlabel('Phase (degrees)')
ylabel('Amplitude frequency (Hz)')
