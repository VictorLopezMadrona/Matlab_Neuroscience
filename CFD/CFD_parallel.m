function CFD=CFD_parallel(data_phase,data_amplitude,varargin)

%% Computes Cross-Frequency Directionality
%
% USE:
%   CFD=CFD_parallel(data_phase,data_amplitude,'param1',value1,...);
%
% INPUT:
%   data_phase (1,samples): Original signal which will be used as theta reference.
%   data_amplitude (1,samples): Original signal which will be used as modulated 
%                               amplitude by the data_phase.
%                               (Optional. Default: data_amplitude = data_phase)
%
% PARAMETERS: 
%   'Fs' - Sampling Frequency 
%   'f_theta' - Information to make the frequency theta vector
%               [f_min, f_max, f_step, BW] 
%   'f_gamma' - Information to make the frequency gamma vector
%               [f_min, f_max, f_step, BW] 
%   'Nsurro' - Number of surrogates to statistical significance
%
% Note: Both data inputs are neccesary before parameters.
%
% See also: main_psi, plot_CFD

% This function is based or uses code from:
% [1] Jiang, H., Bahramisharif, A., van Gerven, M. A., & Jensen, O. (2015). Measuring directionality between neuronal oscillations of different frequencies. Neuroimage, 118, 359-367.
% [2] Guido Nolte, Andreas Ziehe, Vadim Nikulin, Alois Schlögl, Nicole Krämer, Tom Brismar, Klaus-Robert Müller; Robustly estimating the flow direction of information in complex physical systems; Physical Review Letters 100, 234101, 2008
% [3] Niso, G., Bruña, R., Pereda, E., Gutiérrez, R., Bajo, R., Maestú, F., & del-Pozo, F. (2013). HERMES: towards an integrated toolbox to characterize functional and effective brain connectivity. Neuroinformatics, 11(4), 405-434. DOI: 10.1007/s12021-013-9186-1.

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Dec. 2018; Last revision: 14-Jul-2020

%% Initial parameters. 
[errormsg,Fs,f_theta,f_gamma,Nsurro,seglen,B,winlen,overlap]=getparams(varargin{:});
% Error during "getparams"
if ~isempty(errormsg) && nargin > 1
    error(errormsg);
end
if nargin == 1 
    data_amplitude = data_phase;
end

PSI_pval=[];

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

%Loop to get the CFD for each pixel
lx=length(x_theta);
ly=length(y_gamma);
for x=1:lx
    for y=1:ly
        y_min=y_gamma(y)-BW_gamma/2;
        y_max=y_gamma(y)+BW_gamma/2;
        data_gamma=abs(hilbert(eegfilt(data_amplitude,Fs,y_min,y_max)));
        
        data_psi = [data_phase' data_gamma'];
        freqbins = floor((x_theta(x)-B/2)*seglen/Fs) : ceil((x_theta(x)+B/2)*seglen/Fs);     
        [psi_p, stdpsi]=data2psi(data_psi,seglen,seglen*2,freqbins);        
        PSI(y,x) = mean(psi_p(1,2,:));
        if Nsurro>0
            PSI_pval(y,x,:)=compute_surrogate(data_psi,seglen,seglen*2,freqbins,Nsurro);
        end        
        psi_p = psi_p./(stdpsi+eps);

        remaining_iterations = (lx-x+1)*ly-y;
        disp(['Remaining iterations: ' num2str(remaining_iterations)])
    end
end
    
%Prepare the ouput struct
CFD.Fs=Fs;
CFD.Nsurro=Nsurro;
CFD.PSI=PSI;
CFD.pval=PSI_pval;
CFD.f_theta.f_min=f_min_theta;
CFD.f_theta.f_max=f_max_theta;
CFD.f_theta.BW=BW_theta;
CFD.f_theta.step=f_step_theta;
CFD.f_gamma.f_min=f_min_gamma;
CFD.f_gamma.f_max=f_max_gamma;
CFD.f_gamma.BW=BW_gamma;
CFD.f_gamma.step=f_step_gamma;

function PSI_pval=compute_surrogate(data_psi,seglen,seglen2,freqbins,Nsurro)

n=randperm(length(data_psi));
for s=1:Nsurro
    data_surro(:,1)=data_psi(:,1);
    data_surro(:,2)=vertcat(data_psi(n(s):end,2),data_psi(1:n(s)-1,2));
    
    [psi_p, ~]=data2psi(data_surro,seglen,seglen2,freqbins);        
    PSI_pval(s) = mean(psi_p(1,2,:));
end


function [error,Fs,f_theta,f_gamma,Nsurro,seglen,B,winlen,overlap]=getparams(varargin)

%GETPARAMS Process input parameters for comodulograma.
%% Default params. 
Fs=1000;
seglen=Fs*2; % 2s -> 0.5Hz
f_theta=[1.5 12 0.5 1];
f_gamma=[20 100 5 10];
B=4*(Fs/seglen); % 4 x seglen
Nsurro=100;
winlen=Fs*10; % 10 seconds
overlap = round(winlen/2); %overlap 50%

error = [];
%%

while ~isempty(varargin)
   if length(varargin)==1
      error = ('Parameters need a name and its value');
      return
   end
   pname = varargin{1};
   if ~ischar(pname)
      error = ('Invalid name for parameter');
      return
   end
   pvalue = varargin{2};
   j = find(strcmp(pname,{'SegLen','Fs','f_theta','f_gamma','Nsurro','B','WinLen','overlap'}));
   if isempty(j)
      error = ('Parameter does not exist');
      return
   end
   if j==1
       seglen = pvalue;
   elseif j==2
       Fs = pvalue;
   elseif j==3
       f_theta = pvalue;
   elseif j==4
       f_gamma = pvalue; 
   elseif j==5
       Nsurro = pvalue;
   elseif j==6
       B = pvalue;
   elseif j==7
       winlen = pvalue;
   elseif j==8
       overlap = pvalue;
   end
   varargin(1:2) = [];
end
