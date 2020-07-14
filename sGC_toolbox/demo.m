
%% sGC-DEMO
% Example of how to use sGC
% 
% More more info and citation:
% López-Madrona, V. J., Matias, F. S., Mirasso, C. R., Canals, S., & Pereda, E. (2019). Inferring correlations associated to causal interactions in brain signals using autoregressive models. Scientific reports, 9(1), 1-15.

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Jul. 2020; Last revision: 14-Jul-2020

%Add sGC to you path
addpath(genpath(cd))

disp('Loading data')

%Load data
load([cd '/Examples/example_motif.mat']);
disp('Data loaded')
disp(Args.desc)

%Setting parameters
Fs=Args.Fs;
param.window = Fs*5; % 2 seconds in samples
param.overlap = 50; % [0% - 100%[
param.Nsurro = 1000;

disp('Computing sGC')
%Computing sGC
[sGC,GC,pval] = sGC_computesGC(data,param);
sGC = sGC.*(pval<0.05);

%Plotting the results
figure,

f1=subplot(1,2,1); %GC
imagesc(mean(GC,3))
axis('xy','square')
title('GC')
xlabel('FROM')
ylabel('TO')
colormap(f1,'gray')
colorbar

f2=subplot(1,2,2); %sGC 
imagesc(mean(sGC,3))
axis('xy','square')
title('sGC')
xlabel('FROM')
ylabel('TO')
caxis([-1 1])
colormap(f2,flipud(jet))
colorbar






