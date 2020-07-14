function [MI,CFC]=modulation_index(thetaphase,data_gamma,nbins)

%% Function to estimare the Modulation Index from two filtered signals. 
%
% USE:
%   [MI,CFC]=modulation_index(thetaphase,data_gamma,nbins);
%
% INPUT:
%   thetaphase: Phase component (Phase must be already computed)
%   data_gamma: Signal with amplitude component (it will be computed with Hilbert)
%
% OUTPUT:
%    MI:  Modulation Index.           
%    CFC: Average of gamma amplitue during through the cycle phase
%
% See also: comodulogram_parallel

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Nov. 2018; Last revision: 14-Jul-2020

%%
hil_gamma=hilbert(data_gamma);
amplitude=abs(hil_gamma);

% Traditional bins
% thetaphase_bin = ceil( tiedrank( thetaphase ) / (length(thetaphase) / nbins) );
    
% Equalized bins. It divided the cycle into four epochs (0-90; 90-180;
% 180-270; 270-360) and the bin length is different in each epoch.
thetaphase_bin1 = ceil( tiedrank( thetaphase(thetaphase<-pi/2) ) / (length(thetaphase(thetaphase<-pi/2)) / (nbins/4)) );
thetaphase_bin2 = ceil( tiedrank( thetaphase(thetaphase>=-pi/2 & thetaphase<0) ) / (length(thetaphase(thetaphase>=-pi/2 & thetaphase<0)) / (nbins/4)) );
thetaphase_bin3 = ceil( tiedrank( thetaphase(thetaphase>=0  & thetaphase<pi/2) ) / (length(thetaphase(thetaphase>=0 & thetaphase<pi/2))  / (nbins/4)) );
thetaphase_bin4 = ceil( tiedrank( thetaphase(thetaphase>=pi/2) ) / (length(thetaphase(thetaphase>=pi/2)) / (nbins/4)) );

thetaphase_bin = ceil( tiedrank( thetaphase ) / (length(thetaphase) / nbins) );

thetaphase_bin(thetaphase<-pi/2) = thetaphase_bin1;
thetaphase_bin(thetaphase>=-pi/2 & thetaphase<0) = thetaphase_bin2 + nbins/4;
thetaphase_bin(thetaphase>=0 & thetaphase<pi/2) = thetaphase_bin3 + 2*nbins/4;
thetaphase_bin(thetaphase>=pi/2) = thetaphase_bin4 + 3*nbins/4;

%%
gammapow_bin = zeros(1,nbins);
for k=1:nbins
    gammapow_bin(k) = squeeze(mean(amplitude(thetaphase_bin==k)));
end
gammapow_bin = gammapow_bin ./ sum(gammapow_bin);
tmpmi = (log(nbins) + sum(gammapow_bin.*log(gammapow_bin)) ) ./ log(nbins);
MI = tmpmi;
CFC = gammapow_bin;
