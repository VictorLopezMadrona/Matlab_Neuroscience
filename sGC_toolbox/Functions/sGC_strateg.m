function [A_BU,A_TD]=sGC_strateg(X,p,criterion)

% Function to impose zero constraints to A coefficients.
% This function is part of the sGC_toolbox, and it is called during the
% main function "sGC_computesGC()".
%
% USE:
%   [A_BU,A_TD]=sGC_strateg(X,p,criterion);
%
% INPUT:
%   X [nch x samples]: Time-series. Each channel should be in a different row.
%   p: Model order to compute the AR model
%   criterion: information criterion ('AIC' 'SC' 'HQ'). Default 'AIC'.
%
% OUTPUT:
%   A_TD: mask with zero constraints for A matrix after Top-Down strategy.
%   A_BU: mask with zero constraints for A matrix after Bottom-Up strategy.
%
% Developed by V.J. López Madrona - 15/02/2018
% Doubts and comments: v.lopez@umh.es



if nargin == 2
    criterion = 'AIC'; % Akaike Information Criterion
end

if strcmp(criterion,'AIC')
    criterion = 1;
elseif strcmp(criterion,'SC')
    criterion = 2;
elseif strcmp(criterion,'HQ')
    criterion = 3;    
else
    error(['Selected criterion (' criterion ') is not an available option. Use "AIC", "SC" or "HQ".'])
end
    
[nch,samples]=size(X);


%%

%Estimation of matrix Z
Z=NaN*ones(p*nch+1,samples-p);
for i=p:samples-1 
    Z_aux=1;
    for pp=i:-1:i-p+1
        Z_aux=horzcat(Z_aux,X(:,pp)');
    end
    Z(:,i-p+1)=Z_aux';
end

X=X(:,p+1:end); %We have sampled and presampled values
[~,samples]=size(X);

%% Bottom-Up strategy

Rk = zeros(nch*p+1,nch*p+1,nch);
for k = 1:nch % Each equation
    Rk(1,1,k) = 1; % The "mean" component v is always 1
    yk = X(k,:)';
    
    for kk = 1:nch % The contribution of the kk-channel in the k equation
        
        % First iteration without setting any kk-coefficient to 1 (No
        % relation between channels kk and k)
        
        R_aux = Rk(:,:,k);
        R_aux( :, ~any(R_aux,1) ) = [];  %Delete zero-columns
        ck = inv(R_aux'*(Z*Z')*R_aux) * R_aux'*Z*yk;
        bk = R_aux*ck;
        sig = (yk-Z'*bk)' * (yk-Z'*bk) / samples;
            
        if criterion == 1
            AIC(1) = log(sig) + 2/samples*rank(R_aux);
        elseif criterion == 2
            SC(1) = log(sig) + log(samples)/samples*rank(R_aux);
        elseif criterion == 3
            HQ(1) = log(sig) + 2*log(log(samples))/samples*rank(R_aux);
        end
        
        for pi = 1:p % Compute the criterion for each order p
            Rk( (kk-1)*p+pi+1, (kk-1)*p+pi+1, k) = 1; % We put to 1 the pi coefficient (p=pi) of channel kk.
            
            R_aux = Rk(:,:,k);
            R_aux( :, ~any(R_aux,1) ) = [];  %Delete zero-columns
            ck = inv(R_aux'*(Z*Z')*R_aux) * R_aux'*Z*yk;
            bk = R_aux*ck;
            sig = (yk-Z'*bk)' * (yk-Z'*bk) / samples;
            
            if criterion == 1
                AIC(pi+1) = log(sig) + 2/samples*rank(R_aux);
            elseif criterion == 2
                SC(pi+1) = log(sig) + log(samples)/samples*rank(R_aux);
            elseif criterion == 3
                HQ(pi+1) = log(sig) + 2*log(log(samples))/samples*rank(R_aux);
            end
        end
        
        % Find the order pi which minimizes each criterion
        if criterion == 1
            [~,AICp] = min(AIC);
            AICp_TD(k,kk) = AICp;
            for pzero = AICp:p
                Rk( (kk-1)*p+pzero+1, (kk-1)*p+pzero+1, k) = 0; %Set to zero the orders higher than the minimum
            end
        elseif criterion == 2
            [~,SCp] = min(SC);
            SCp_TD(k,kk) = SCp;
            for pzero = SCp:p
                Rk( (kk-1)*p+pzero+1, (kk-1)*p+pzero+1, k) = 0; %Set to zero the orders higher than the minimum
            end
        elseif criterion == 3
            [~,HQp] = min(HQ);
            HQp_TD(k,kk) = HQp;
            for pzero = HQp:p
                Rk( (kk-1)*p+pzero+1, (kk-1)*p+pzero+1, k) = 0; %Set to zero the orders higher than the minimum
            end
        end
    end
end

% Create mask for Bottom-Up Strategy
R = Rk(:,:,1);
for i=2:nch
    R = blkdiag(R,Rk(:,:,i));
end
Rd = diag(R);
A_BU = zeros(nch,nch,p);
% Create Am mask
for i=1:nch
    B_BU(i,:) = Rd((i-1)*(nch*p+1)+1:i*(nch*p+1));
end 
for i=1:p
    A_BU(:,:,i)=B_BU(:,(i-1)*nch+2:i*nch+1);
end
v_BU = B_BU(:,1);


%% Top-Down strategy
for k=1:nch 
    Rk(:,:,k) = eye(nch*p+1); %Identity matrix for the restriction matrix
    yk = X(k,:)';
    
    R_aux = Rk(:,:,k);    
    ck = (R_aux'*Z*Z'*(inv(R_aux))) * R_aux'*Z*yk; 
    bk = R_aux*ck;
    sig = (yk-Z'*bk)' * (yk-Z'*bk) / samples;
    AIC = log(sig) + 2/samples*rank(R_aux);
    SC = log(sig) + log(samples)/samples*rank(R_aux);
    HQ = log(sig) + 2*log(log(samples))/samples*rank(R_aux);
   
    Rk((nch*p+1),(nch*p+1),k) = 0; %First coefficient to test
    for kk = (nch*p+1):-1:1
        
        R_aux = Rk(:,:,k);
        R_aux( :, ~any(R_aux,1) ) = [];  %Delete zero-columns
        ck = inv(R_aux'*(Z*Z')*R_aux) * R_aux'*Z*yk;
        bk = R_aux*ck;
        sig = (yk-Z'*bk)' * (yk-Z'*bk) / samples;
        
        if criterion == 1
            AIC_new = log(sig) + 2/samples*rank(R_aux);
            if AIC_new <= AIC
                AIC = AIC_new; %Update the value of the criterion
            else
                Rk(kk,kk,k) = 1; %This coefficient is necessary. 
            end
            
        elseif criterion == 2
            SC_new = log(sig) + log(samples)/samples*rank(R_aux);
            if SC_new <= SC
                SC = SC_new; %Update the value of the criterion
            else
                Rk(kk,kk,k) = 1; %This coefficient is necessary. 
            end
        
        elseif criterion == 3
            HQ_new = log(sig) + 2*log(log(samples))/samples*rank(R_aux);
            if HQ_new <= HQ
                HQ = HQ_new; %Update the value of the criterion
            else
                Rk(kk,kk,k) = 1; %This coefficient is necessary. 
            end
        end
        
        if kk>1
            Rk((kk-1),(kk-1),k) = 0; %Next coefficient to test
        end
    end
end

% Create mask for Top-Down Strategy
R = Rk(:,:,1);
for i=2:nch
    R = blkdiag(R,Rk(:,:,i));
end
Rd = diag(R);
A_TD = zeros(nch,nch,p);
% Create Am mask
for i=1:nch
    B_TD(i,:) = Rd((i-1)*(nch*p+1)+1:i*(nch*p+1));
end 
for i=1:p
    A_TD(:,:,i)=B_TD(:,(i-1)*nch+2:i*nch+1);
end
v_TD = B_TD(:,1);


