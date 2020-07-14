function plot_CFD(CFD)

%% Plot CFD
%
% See also: CFD_parallel

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Jul. 2020; Last revision: 14-Jul-2020

%Extract the information from the struct
f_min_theta=CFD.f_theta.f_min;
f_max_theta=CFD.f_theta.f_max;
f_min_gamma=CFD.f_gamma.f_min;
f_max_gamma=CFD.f_gamma.f_max;

x=[f_min_theta f_max_theta];
y=[f_min_gamma f_max_gamma];

%% CFD plot
figure,
imagesc_filter(CFD.PSI,4,x,y);

xlabel('Phase frequency (Hz)'),
ylabel('Amplitude frequency (Hz)'),
colormap('parula')
h = colorbar;
set(get(h,'title'),'string','CFD');
