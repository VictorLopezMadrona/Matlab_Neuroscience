function [GC,A] = sGC_GC(X,p)

% Function to estimate Granger Causality (GC) between signals in a given dataset X.
% This function is part of the sGC_toolbox, and it is called during the
% main function "sGC_computesGC()".
% The functions called in this code are part of MVGC_toolbox, developed by
% Lionel Barnett and Anil K. Seth. 
% See "http://users.sussex.ac.uk/~lionelb/MVGC/" for more information.
%
% USE:
%   [GC,A] = sGC_GC(X,p);
%
% INPUT:
%   X [nch x samples]: Time-series. Each channel should be in a different row.
%   p: Model order to compute the AR model
%
% OUTPUT:
%   GC [nch x nch]: Matrix with the Granger-Causality results.
%   A [nch x nch x p]: Coefficients of the AR model.
%
% Developed by V.J. López-Madrona - 15/02/2018
% Doubts and comments: v.lopez@umh.es


[nch, samples, trials]=size(X);
alpha = 0.05;
mhtc='FDR';
[A,SIG] = tsdata_to_var(X,p);

[G,~] = var_to_autocov(A,SIG);
F = autocov_to_pwcgc(G);
pval = mvgc_pval(F,p,samples,trials,1,1,nch-2); % take careful note of arguments!
sig  = significance(pval,alpha,mhtc);
GC=F.*sig;



