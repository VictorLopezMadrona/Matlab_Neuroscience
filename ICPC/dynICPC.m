function [ICPC,sample_cycle] = dynICPC(x_phase)

%% Computes Inter Cycle Phase Clustering (based on ITPC) across cycles
%
% USE:
%   [ICPC,sample_cycle] = dynICPC(x_phase);
%   
% INPUT:
%   x_phase - [Nch, samples] Phase of the signals in radians.
% 
% OUTPUT:
%   ICPC - ICPC value for each cycle and pair of signals.
%   sample_cycle - Sample point that corresponds to each ICPC value
%
% See also: 

% For more info:
%   Lopez-Madrona, V. J., Perez-Montoyo, E., Alvarez-Salvado, E., Moratal,
%   D., Herreras, O., Pereda, E., ... & Canals, S. (2020). Different theta
%   frameworks coexist in the rat hippocampus and are coordinated during
%   memory-guided exploration and novelty detection
%

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Jul. 2020; Last revision: 15-Jul-2020

Nch = size(x_phase,1);
if Nch <= 2
    error('At least two signals are required.')
end

peaks = cell(1,Nch);
for ch=1:Nch
    [~,peaks{ch}] = findpeaks(x_phase(ch,:),'MinPeakWidth',2.5);
end

ICPC = cell(1,Nch);
for chi=1:Nch
    for chj=1:Nch
        if chi~=chj
        for p=2:length(peaks{chi})-1
            ICPC{chi}(chj,p) = abs(mean(exp(1i* x_phase(chj,peaks{chi}(p-1:p+1)))));
        end
        end
    end
end

sample_cycle = peaks;

