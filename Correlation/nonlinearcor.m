
function [h2,lag_max,h2_max] = nonlinearcor(x,y,maxlag,bins)

%% Compute nonlinear correlation from x to y using h2.
% 
% Syntax:  
%   [h2,lag_max,h2_max] = nonlinearcor(x,y,maxlag,bins);
%
% Inputs:
%   x -
%   y -
%   maxlag - [Def=0] maximum delay to analyze in samples
%   bins   - [Def=8] number of bins for 
%
% Outputs:
%   h2
%   lag_max - lag for the maximum h2 value
%   h2_max  - h2 value at lag lag_max
%
% See also: 

% References: 
% da Silva, F. L., Pijn, J. P., & Boeijinga, P. (1989). 
% Interdependence of EEG signals: linear vs. nonlinear associations and the
% significance of time delays and phase shifts. Brain topography, 2(1-2), 9-18.
%
% Wendling, F., Bartolomei, F., Bellanger, J. J., & Chauvel, P. (2001). 
% Interpretation of interdependencies in epileptic signals using a macroscopic 
% physiological model of the EEG. Clinical neurophysiology, 112(7), 1201-1218.

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Feb. 2020; Last revision: 27-Apr-2020

if nargin == 2 
    maxlag = 0;
    bins = 8;
end
if nargin == 3, bins = 8; end
h2 = zeros(1,length(-maxlag:maxlag));

c=1;
for lag=-maxlag:maxlag
    
    if lag<0
        x_aux = x(1:end+lag);
        y_aux = y(-lag+1:end);
    elseif lag == 0
        x_aux = x;
        y_aux = y;
    elseif lag>0
        x_aux = x(lag+1:end);
        y_aux = y(1:end-lag);
    end
        

    x_bin = ceil( tiedrank( x_aux ) / (length(x_aux) / bins) );
    %binsize = (max(x)-min(x))/bins;
    %pi = (min(x)+binsize/2) : binsize : (max(x)-binsize/2);

    qi = zeros(1,bins);
    pi = zeros(1,bins);
    for k=1:bins
        pi(k) = squeeze(mean(x_aux(x_bin==k))); % -- compute mean of x in each bin
        qi(k) = squeeze(mean(y_aux(x_bin==k))); % -- compute mean of y in each bin
    end

    %We want the intervals between bins, so we divide the signal x in 2*bins, the
    %first three bins will correspond to the first segment, the last three to
    %the segment (bins-1) and each two bin will be a different segment.

    x_bin2 = ceil( tiedrank( x_aux ) / (length(x_aux) / (bins*2)) ); 
    x_bin2(x_bin2==1 | x_bin2==2 | x_bin2==3) = 1;
    for k=2:bins-2
        x_bin2(x_bin2==2*k | x_bin2==(2*k+1)) = k;
    end
    x_bin2(x_bin2==(2*bins-2) | x_bin2==(2*bins-1) | x_bin2==2*bins) = bins-1;

    gi = [diff(qi)./diff(pi); -((diff(qi)./diff(pi)).*pi(1:bins-1))+qi(1:bins-1)];  


    %We compute y' or h(x), which is the estimation of y based on x
    hx = y_aux*0;
    for k=1:bins-1
        hx(x_bin2==k) = x_aux(x_bin2==k)*gi(1,k) + gi(2,k);
    end

    % figure, 
    % plot(x_aux,y_aux,'*')
    % hold on,
    % plot(x_aux,hx,'*')
    % plot(pi,qi)

    h2(c) = 1-(var(y_aux-hx))/var(y_aux);
    %h2(c) = (var(y)-var(y-hx))/var(y);
    c=c+1;
end

lags = -maxlag:maxlag;
[h2_max,p]=max(abs(h2));
lag_max = lags(p);






