function [sGC, Amax] = sGC_sGC(A,A_TD,A_BU,GC,Amax)

%Actualiar con Amax

% Function to estimate synaptic Granger Causality (sGC) from A coefficients.
% This function is part of the sGC_toolbox, and it is called during the
% main function "sGC_computesGC()".
%
% USE:
%   sGC = sGC_sGC(A,A_TD,A_BU,GC);
%
% INPUT:
%   A: matrix with the coefficients of the Autoregressive Model.
%   A_TD: mask with zero constraints for A matrix after Top-Down strategy.
%   A_BU: mask with zero constraints for A matrix after Bottom-Up strategy.
%
% OUTPUT:
%   sGC: matrix with sGC measurement for each pair of signals.
%
%
% Developed by V.J. López Madrona - 15/02/2018
% Doubts and comments: v.lopez@umh.es

[nch,~,p] = size(A);
if nargin == 1
    A_TD = (A*0)+1;
    A_BU = (A*0)+1;
    GC = (A(:,:,1)*0)+1;
    Amax = zeros(nch);
elseif nargin == 2
    A_BU = (A*0)+1;
    GC = (A(:,:,1)*0)+1;
    Amax = zeros(nch);
elseif nargin == 3
    GC = (A(:,:,1)*0)+1;
    Amax = zeros(nch);
elseif nargin == 4
    Amax = zeros(nch);
elseif nargin == 5
    A_TD = (A*0)+1;
    A_BU = (A*0)+1;
    GC = (A(:,:,1)*0)+1;
end


%% synaptic Granger Causality Result

A_sGC = A.*A_BU.*A_TD.*(repmat(GC,1,1,p)>0);
sGC = zeros(nch);
for ni=1:nch
    for nj=1:nch
        if Amax(ni,nj)==0
            Amax(ni,nj) = max(sum(A_sGC(ni,nj,A_sGC(ni,nj,:)>0).^2) , sum(A_sGC(ni,nj,A_sGC(ni,nj,:)<0).^2));
        end
        if Amax(ni,nj)~=0
            sGC(ni,nj) = sum(A_sGC(ni,nj,A_sGC(ni,nj,:)>0).^2)/Amax(ni,nj) - sum(A_sGC(ni,nj,A_sGC(ni,nj,:)<0).^2)/Amax(ni,nj);
        end
    end
end

