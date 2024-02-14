function comodulogram=computeCFC(cnfg,data_phase,data_amplitude)

%% Compute phase-amplitude CFC using Modulation Index 
% It computes the phase and amplitude with filters and Hilbert transform.
%
% USE:
%   comodulogram=computeCFC(cfg,data_phase,data_amplitude);
%
% INPUT:
%   data_phase (1,samples): Signal which will be used as phase reference.
%   data_amplitude(1,samples): (Optional) Signal to analyze amplitude coupling.
%                               
%   (Optional. Default: data_amplitude = data_phase)
%
% PARAMETERS: 
%   cfg.Fs      - Sampling frequency
%   cfg.bins    - Number of bins to divide each theta cycle (Def=16)
%   cfg.f_phase - Information to make the frequency phase vector
%                 [f_min, f_max, f_step, BW] 
%   cfg.f_amp   - Information to make the frequency amplitude vector
%                 [f_min, f_max, f_step, BW] 
%   cfg.Nsurro  - Number of surrogates to statistical significance (Def=0)
%
% OUTPUT:
%   comodulogram =
%       Fs: Sampling Frequency (Hz).
%       bins: Number of bins.
%       f_theta: Information about the frequency theta vector.
%       f_gamma: Information about the frequency gamma vector.
%       MI: (N_gamma_pixels, N_theta_pixels) Matrix with the Modulation.
%           Index of each gamma-theta frequency.
%       CFC: (N_gamma_pixels, N_theta_pixels, 2*bins) Matrix with the
%            average of gamma amplitue during a cycle for each 
%            gamma-theta frequency.
%
% See also: modulation_index plot_comodulogram

% This code measures the Cross-Frequency Coupling (CFC) between the phase
% and he amplitude of two signals (can be the same) following the method
% in [1], and a statistical analysis based on block-resampling at single
% locations [2]. For more information about common errors in CFC see [3].
%
% [1] Tort, A. B., Komorowski, R., Eichenbaum, H., & Kopell, N. (2010). Measuring phase-amplitude coupling between neuronal oscillations of different frequencies. Journal of neurophysiology, 104(2), 1195-1210.
% [2] Canolty, R. T., Edwards, E., Dalal, S. S., Soltani, M., Nagarajan, S. S., Kirsch, H. E., ... & Knight, R. T. (2006). High gamma power is phase-locked to theta oscillations in human neocortex. science, 313(5793), 1626-1628.
% [3] Aru, J., Aru, J., Priesemann, V., Wibral, M., Lana, L., Pipa, G., ... & Vicente, R. (2015). Untangling cross-frequency coupling in neuroscience. Current opinion in neurobiology, 31, 51-61.

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Dec. 2016; Last revision: 16-Aug-2021

% Change log:
% 16/8/2021: Updated to enter the parameters in a cfg struct

%% Initial parameters. 

if ~isfield(cnfg,'Fs')
    error('Fs parameter should be specified'); end
if ~isfield(cnfg,'bins')
    cnfg.bins = 16; end
if ~isfield(cnfg,'f_phase')
    cnfg.f_phase = [2 16 1 2]; end
if ~isfield(cnfg,'f_amp')
    cnfg.f_amp = [20 100 5 20]; end
if ~isfield(cnfg,'Nsurro')
    cnfg.Nsurro = 0; end

Fs      = cnfg.Fs;
bins    = cnfg.bins;
Nsurro  = cnfg.Nsurro;
f_theta = cnfg.f_phase;
f_gamma = cnfg.f_amp;

if nargin == 2 
    data_amplitude = data_phase;
end

%Move the frequency information to different variables
f_min_theta=f_theta(1);
f_max_theta=f_theta(2);
f_step_theta=f_theta(3);
BW_theta=f_theta(4);

f_min_gamma=f_gamma(1);
f_max_gamma=f_gamma(2);
f_step_gamma=f_gamma(3);
BW_gamma=f_gamma(4);

%Check if the number of bins is too big for any theta frequency
if round(Fs/f_max_theta)+2 < bins
    disp('Solutions: -Increase de sampling frequency')
    disp('           -Reduce the number of bins')
    disp('           -Reduce f_max_theta')
    error('Number of bins too big to divide the theta cycles in all frequencies');
end

%Make the frequency vectors 
x_theta=(f_min_theta:f_step_theta:f_max_theta);
y_gamma=(f_min_gamma:f_step_gamma:f_max_gamma);

MI=zeros(length(y_gamma),length(x_theta));
MI_pval = zeros(length(y_gamma),length(x_theta),Nsurro);
CFC=zeros(length(y_gamma),length(x_theta),bins);    

K=length(x_theta)*length(y_gamma); %Total iterations
nCores=feature('numCores'); %Number of cores 

%First iteration, just to make a rude estimation of time.
tic
  x=1;
    x_min=x_theta(x)-BW_theta/2;
    x_max=x_theta(x)+BW_theta/2;
    data_theta=eegfilt(data_phase,Fs,x_min,x_max);
  y=1;
    y_min=y_gamma(y)-BW_gamma/2;
    y_max=y_gamma(y)+BW_gamma/2;
    data_gamma=eegfilt(data_amplitude,Fs,y_min,y_max);
    x_phase = angle(hilbert(data_theta));
    modulation_index(x_phase,data_gamma,bins);
    if Nsurro > 0
        compute_surrogate(x_phase,data_gamma,bins,Nsurro);
    end
toc
clc
tiempo_toc=toc;
tiempo_res=K*tiempo_toc/nCores;
TIME=sec2hms(tiempo_res);
barra2=['Estimated remaining time: ',TIME];
barra=waitbar(0,barra2);    

%Loop to get the Modulation Index for each pixel
lx=length(x_theta);
ly=length(y_gamma);
for x=1:lx
    x_min=x_theta(x)-BW_theta/2;
    x_max=x_theta(x)+BW_theta/2;
    data_theta=eegfilt(data_phase,Fs,x_min,x_max);    
    parfor y=1:ly
        y_min=y_gamma(y)-BW_gamma/2;
        y_max=y_gamma(y)+BW_gamma/2;
        data_gamma=eegfilt(data_amplitude,Fs,y_min,y_max);
        
        %Compute phase. Other methods can be applyed here:
        x_phase = angle(hilbert(data_theta));
        [MI(y,x),CFC(y,x,:)]=modulation_index(x_phase,data_gamma,bins);
        if Nsurro > 0
            MI_pval(y,x,:)=compute_surrogate(x_phase,data_gamma,bins,Nsurro);
        end
        
        remaining_iterations = (lx-x+1)*ly-y;
        disp(['Remaining iterations: ' num2str(remaining_iterations)])
    end
end
close(barra)
    
%Prepare the ouput struct
comodulogram.Fs=Fs;
comodulogram.Nsurro=Nsurro;
comodulogram.bins=bins;
comodulogram.MI=MI;
comodulogram.MI_pval=MI_pval;
comodulogram.CFC=CFC;
comodulogram.f_theta.f_min=f_min_theta;
comodulogram.f_theta.f_max=f_max_theta;
comodulogram.f_theta.BW=BW_theta;
comodulogram.f_theta.step=f_step_theta;
comodulogram.f_gamma.f_min=f_min_gamma;
comodulogram.f_gamma.f_max=f_max_gamma;
comodulogram.f_gamma.BW=BW_gamma;
comodulogram.f_gamma.step=f_step_gamma;

function MI_surro=compute_surrogate(x_phase,data_gamma,bins,Nsurro)

%Create surrogates and estimate the normalized value of Modulation Index.
[~,peaks]=findpeaks(x_phase);
MI_surro=zeros(1,Nsurro);

for s=1:Nsurro
    break_point = round( (length(peaks)-1)*rand+1);
    theta_surro = [x_phase(peaks(break_point)+1:end), x_phase(1:peaks(break_point))];
    MI_surro(s) = modulation_index(theta_surro,data_gamma,bins);
end











