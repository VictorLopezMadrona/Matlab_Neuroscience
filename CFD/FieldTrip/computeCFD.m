function CFD=computeCFD(cnfg,data_phase,data_amplitude)

%% Compute phase-amplitude CFD using Phase Slope Index 
%
% USE:
%   CFD=computeCFD(cfg,data_phase,data_amplitude);
%
% INPUT:
%   data_phase (1,samples): Signal which will be used as phase reference.
%   data_amplitude(1,samples): (Optional) Signal to analyze amplitude coupling.
%                               
%   (Optional. Default: data_amplitude = data_phase)
%
% PARAMETERS: 
%   cfg.Fs      - Sampling frequency
%   cfg.f_phase - Information to make the frequency phase vector
%                 [f_min, f_max, f_step, BW] 
%   cfg.f_amp   - Information to make the frequency amplitude vector
%                 [f_min, f_max, f_step, BW] 
%   cfg.Nsurro  - Number of surrogates to statistical significance (Def=0)
%   cfg.seglen  - [Def: Fs*2]
%   cfg.B       - [Def: 4*(Fs/seglen)]
%
% OUTPUT:
%   comodulogram =
%
% See also: CFD_ft plot_CFD_ft

% This function is based or uses code from:
% [1] Jiang, H., Bahramisharif, A., van Gerven, M. A., & Jensen, O. (2015). Measuring directionality between neuronal oscillations of different frequencies. Neuroimage, 118, 359-367.
% [2] Guido Nolte, Andreas Ziehe, Vadim Nikulin, Alois Schlögl, Nicole Krämer, Tom Brismar, Klaus-Robert Müller; Robustly estimating the flow direction of information in complex physical systems; Physical Review Letters 100, 234101, 2008
% [3] Niso, G., Bruña, R., Pereda, E., Gutiérrez, R., Bajo, R., Maestú, F., & del-Pozo, F. (2013). HERMES: towards an integrated toolbox to characterize functional and effective brain connectivity. Neuroinformatics, 11(4), 405-434. DOI: 10.1007/s12021-013-9186-1.

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Feb 23; Last revision: 09-Jun-2023

% Change log:
% 09/06/2023: Corrected error in freqbins

%% Initial parameters. 

if ~isfield(cnfg,'Fs')
    error('Fs parameter should be specified'); 
else, cnfg.Fs = round(cnfg.Fs); end
if ~isfield(cnfg,'f_phase')
    cnfg.f_phase = [2 16 1 2]; end
if ~isfield(cnfg,'f_amp')
    cnfg.f_amp = [20 100 5 20]; end
if ~isfield(cnfg,'Nsurro')
    cnfg.Nsurro = 0; end
if ~isfield(cnfg,'seglen')
    cnfg.seglen = cnfg.Fs*2; end
if ~isfield(cnfg,'B')
    cnfg.B = 4*(cnfg.Fs/cnfg.seglen); end

Fs      = cnfg.Fs;
Nsurro  = cnfg.Nsurro;
f_theta = cnfg.f_phase;
f_gamma = cnfg.f_amp;
seglen  = cnfg.seglen;
B       = cnfg.B;

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

%Make the frequency vectors 
x_theta=(f_min_theta:f_step_theta:f_max_theta);
y_gamma=(f_min_gamma:f_step_gamma:f_max_gamma);


%Loop to get the PSI for each pixel
lx=length(x_theta);
ly=length(y_gamma);
PSI=zeros(ly,lx);
PSI_pval=zeros(ly,lx,Nsurro);
for x=1:lx
    parfor y=1:ly %parfor
        y_min=y_gamma(y)-BW_gamma/2;
        y_max=y_gamma(y)+BW_gamma/2;
        data_gamma=abs(hilbert(eegfilt(data_amplitude,Fs,y_min,y_max)));
        
        data_psi = [data_phase' data_gamma'];
        freqbins = floor((x_theta(x)-B/2)*seglen/Fs) : ceil((x_theta(x)+B/2)*seglen/Fs) + 1;   %I must add one, because f(1) is frequency 0, not 1 Hz.        
        [psi_p]=data2psi(data_psi,seglen,[],freqbins); %faster
        %[psi_p, stdpsi]=data2psi(data_psi,seglen,seglen*2,freqbins);        
        %psi_p = psi_p./(stdpsi+eps);
        PSI(y,x) = mean(psi_p(1,2,:));
        if Nsurro>0
            PSI_pval(y,x,:)=compute_surrogate(data_psi,seglen,seglen*2,freqbins,Nsurro);
        end        

        remaining_iterations = (lx-x+1)*ly-y;
        disp(['Remaining iterations: ' num2str(remaining_iterations)])
    end
end
    
%Prepare the ouput struct
CFD.Fs=Fs;
CFD.Nsurro=Nsurro;
CFD.PSI=PSI;
CFD.pval=PSI_pval;
CFD.f_phase.f_min=f_min_theta;
CFD.f_phase.f_max=f_max_theta;
CFD.f_phase.BW=BW_theta;
CFD.f_phase.step=f_step_theta;
CFD.f_amp.f_min=f_min_gamma;
CFD.f_amp.f_max=f_max_gamma;
CFD.f_amp.BW=BW_gamma;
CFD.f_amp.step=f_step_gamma;

function PSI_pval=compute_surrogate(data_psi,seglen,seglen2,freqbins,Nsurro)

n=randperm(length(data_psi));
PSI_pval = zeros(1,Nsurro);
for s=1:Nsurro
    data_surro(:,1)=data_psi(:,1);
    data_surro(:,2)=vertcat(data_psi(n(s):end,2),data_psi(1:n(s)-1,2));
    
    psi_p=data2psi(data_surro,seglen,[],freqbins);     
    %[psi_p, stdpsi]=data2psi(data_surro,seglen,seglen2,freqbins); 
    %psi_p = psi_p./(stdpsi+eps);
    PSI_pval(s) = mean(psi_p(1,2,:));
end
