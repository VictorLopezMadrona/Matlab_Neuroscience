function [sGC,GC,pval,jaja] = sGC_computesGC(X,param)

% Main function to compute synaptic Granger Causality
% This function is part of the sGC_toolbox.
%
% USE:
%   [sGC,GC,pval] = sGC_computesGC(X,param);
%
% INPUT:
%   X [nch x samples]: Time-series. Each channel should be in a different row.
%   param: struct with parameters:
%        .window: length of the sliding window. Default = length(data);
%        .overlap: [0% - 100%[ Overlap of the sliding windows. Default = 0%.
%        .p: Model order. If p=0, estimates model order. Default = 0.
%        .Nsurro: Number of surrogates. Default = 0.
%        .criterion: Information criterion ('AIC' 'SC' 'HQ'). Default 'AIC'.
%
% OUTPUT:
%   sGC: matrix with sGC measurement for each pair of signals.
%   GC: Matrix with the Granger-Causality results.
%   pval: p-value associated to each sGC value.
%
% Developed by V.J. López Madrona - 15/02/2018
% Doubts and comments: v.lopez@umh.es


%% PARAMETERS

if nargin == 1 %Default values
    [~,param.window] = size(X);
    param.overlap = 0;
    param.p = 0;
    param.Nsurro = 0;
    param.criterion = 'AIC';
end

if ~isfield(param,'window')
    [~,param.window] = size(X);
end
if ~isfield(param,'overlap')
    param.overlap = 0;
end
if ~isfield(param,'p') 
    param.p = 0; 
end
if ~isfield(param,'Nsurro')
    param.Nsurro = 0;
end
if ~isfield(param,'criterion')
    param.criterion = 'AIC';
end

window=param.window;
overlap=param.overlap;
p=param.p;
Nsurro=param.Nsurro;
criterion=param.criterion;

%% CODE
[Nch,samples] = size(X);

if overlap >= 100
    error('Overlap cannot be equal or superior 100%')
end

slidw = round(window*(100-overlap)/100); 
Nw = floor((samples-window)/slidw)+1;

GC = zeros(Nch,Nch,Nw);
sGC = zeros(Nch,Nch,Nw);

%% Estimate model order
if p==0
wb = waitbar(0,'Estimating model order. Estimating remaining time...','Name','sGC toolbox');
t_pre = 0;
pw=zeros(1,Nw);
for w=1:Nw
    tic
    [~,~,pw(w),~] = tsdata_to_infocrit(X(:,(w-1)*slidw+1:(w-1)*slidw+window),50);

    %Remaining time
    t=toc;
    t_pre = (t_pre*(w-1)+t)/w;
    t_rem = (Nw-w)*t_pre;
    waitbar(w/Nw,wb,['Estimating model order. Remaining time: ' sec2hms(t_rem)])
end

param.p = round(mean(pw));
p = round(mean(pw));
disp(['Estimated model order = ' num2str(p)]);
close(wb)
end

%% Compute sGC
wb = waitbar(0,'Computing sGC. Estimating remaining time...','Name','sGC toolbox');
t_pre = 0;
for w=1:Nw
    tic
    p=15
    [A_BU,A_TD]=sGC_strateg(X(:,(w-1)*slidw+1:(w-1)*slidw+window),p,criterion);
    [GC(:,:,w),A] = sGC_GC(X(:,(w-1)*slidw+1:(w-1)*slidw+window),p);
    miau=(A*0)+1;
    miau2=sum(A_BU.*A_TD);
    miau3=A.*A_BU.*A_TD;
    %sum(abs(A(:)))
    %sum(abs(miau3(:)))
    jaja=[sum(miau(:))   sum(miau(:))-sum(A_BU(:))  sum(A_BU(:))-sum(miau2(:)) 100-sum(abs(miau3(:)))/sum(abs(A(:)))*100] 
    %sum(A_BU(:))
    %sum(A_TD(:))
    miau=sum(A_BU.*A_TD);
    %sum(miau(:))
    %sGC(:,:,w) = sGC_sGC(A,A_TD,A_BU,GC(:,:,w));
    %A_TD = (A*0)+1;
    %A_BU = (A*0)+1;
    [sGC(:,:,w), Amax(:,:,w)] = sGC_sGC(A,A_TD,A_BU,GC(:,:,w));
    
    %Remaining time
    t=toc;
    t_pre = (t_pre*(w-1)+t)/w;
    t_rem = (Nw-w)*t_pre;
    waitbar(w/Nw,wb,['Computing sGC. Remaining time: ' sec2hms(t_rem)])
end
close(wb)

%% Surrogates
if Nsurro>0
    
if Nw < Nch
    error('Number of windows is not long enough to perform permuting surrogates. Try to reduce window length or increase the overlap')
elseif Nw < (Nch*3)
    warning('Number of windows too low for good surrogates by permuting. Try to reduce window length or increase the overlap')
end

wb = waitbar(0,'Computing statistics. Estimating remaining time...','Name','sGC toolbox');
t_pre = 0;
X_surro = zeros(Nch,window);
sGC_s = zeros(Nch,Nch,Nsurro);
for s=1:Nsurro
    tic
    %Creating surrogate
    sperm=randperm(Nw,Nch);
    for n=1:Nch
        X_surro(n,:) = X(n,(sperm(n)-1)*slidw+1:(sperm(n)-1)*slidw+window);
    end
    [A] = tsdata_to_var(X_surro,p);
    sGC_s(:,:,s) = sGC_sGC(A,1,1,1,Amax);
    %sGC_s(:,:,s) = sGC_sGC(A);
    
    %Remaining time
    t=toc;
    t_pre = (t_pre*(s-1)+t)/s;
    t_rem = (Nsurro-s)*t_pre;
    waitbar(s/Nsurro,wb,['Computing statistics. Remaining time: ' sec2hms(t_rem)])
end

pval = sGC_statistics(sGC,sGC_s);

close(wb)

else
pval=sGC*0;
end








