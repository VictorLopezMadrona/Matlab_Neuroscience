function [pval] = sGC_statistics(sGC,sGC_s)

% Function to compute p-value of sGC using estimated surrogates.
% This function is part of the sGC_toolbox, and it is called during the
% main function "sGC_computesGC()".
%
% USE:
%   [pval] = sGC_statistics(sGC,sGC_s);
%
% INPUT:
%   sGC: synaptic Granger Causality matrix.
%   sGC_s: surrogate sGC matrices.
%
% OUTPUT:
%   pval: p-value associated to each sGC value.
%
% Developed by V.J. López Madrona - 15/02/2018
% Doubts and comments: v.lopez@umh.es

[Nch,~,Nw]=size(sGC);
aver = mean(sGC_s,3);
desvest = zeros(Nch);
pval = zeros(Nch,Nch,Nw);
val=-10:0.001:10;
for i=1:Nch
    for j=1:Nch
        %figure, histogram(sGC_s(i,j,:)), title([num2str(i) num2str(j)])
        desvest(i,j) = std(sGC_s(i,j,:));
        disp(['Link ' num2str(i) num2str(j) ', Mean = ' num2str(aver(i,j)) ', SD = ' num2str(desvest(i,j))])
        normdist = normpdf(val,aver(i,j),desvest(i,j));
        for w=1:Nw
            if sGC(i,j,w)>=0
                pval(i,j,w) = sum(normdist(val>sGC(i,j,w)))/sum(normdist);
            else
                pval(i,j,w) = sum(normdist(val<sGC(i,j,w)))/sum(normdist);
            end
        end
    end
end







