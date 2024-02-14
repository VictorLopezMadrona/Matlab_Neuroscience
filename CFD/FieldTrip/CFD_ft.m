function CFD=CFD_ft(cnfg,ftdata)

%% Compute phase-amplitude CFD with fieldtrip structures 
%
% USE:
%   CFD=CFD_ft(cfg,ft_data);
%
% PARAMETERS (cfg.): 
%   
%   latency - [t_start t_end]
%   channel - (Def='all') Channel to compute intraCFD
%   
% If you want to compute CFD between phases and amplitudes of different
% signals, it is possible to define a vector of channels with the
% specific pairs (ch_interCFD). For example: ch_interCFD = [1 1; 1 2; 1 3]
% will compute the CFD between the phase of ch1 and the amplitudes of
% channels 1, 2 and 3. If defined, cfg.channel won't be used.
%
%   ch_interCFD (N x 2)- Two-column vector with pairs of channels for CFD.                       
%
% For CFD analysis:
%
%   'f_phase' - Information to make the frequency theta vector
%               [f_min, f_max, f_step, BW] 
%   'f_amp'   - Information to make the frequency gamma vector
%               [f_min, f_max, f_step, BW] 
%   'Nsurro'  - Number of surrogates to statistical significance
%
% To save the data:
%
%   doplot   - true/false (Def: true)
%   dosave   - true/false (Def: false)
%   outpath  - string with path 
%   infosave - string with filename
%
% OUTPUT:
%   CFD = vector of cells:
%
% See also: computeCFD plot_CFD_ft comodulogram_ft

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Feb 23; Last revision: 20-Feb-2023

% Change log:
% 16/08/23 - Modified to store statistical mask too
% 20/02/23 - Corrected bug during inter-regional CFD

%% Initial parameters:

if ~isfield(cnfg,'latency')
    cnfg.latency=[ftdata.time{1}(1) ftdata.time{1}(end)]; end
if isfield(cnfg,'ch_interCFD') && isfield(cnfg,'channel')
    if size(cnfg.ch_interCFD,2) ~= 2 %Should have two columns
        error('cnfg.ch_interCFD must have two columns'), end
    warning('ch_interCFD was defined. Information in cfg.channel won''t be used.') 
    cnfg=rmfield(cnfg,'channel'); 
end
if ~isfield(cnfg,'ch_interCFD') && ~isfield(cnfg,'channel')
    cnfg.channel = 'all'; end
if ~isfield(cnfg,'f_phase')
    cnfg.f_phase = [3 16 1 2]; end
if ~isfield(cnfg,'f_amp')
    cnfg.f_amp = [20 120 5 20]; end
if ~isfield(cnfg,'Nsurro')
    cnfg.Nsurro = 0; end
if ~isfield(cnfg,'dosave'), cnfg.dosave=false; end
if ~isfield(cnfg,'doplot'), cnfg.doplot=true; end
if isfield(cnfg,'outpath')
    if ~strcmp(cnfg.outpath(end),'\')
        cnfg.outpath = [cnfg.outpath '\'];
    end
end
if ~isfield(cnfg,'infosave'), cnfg.infosave=''; end

%% COMPUTE CFD

% Case "intraCFD"
if isfield(cnfg,'channel')
    cfg = [];
    cfg.latency = cnfg.latency;
    cfg.channel = cnfg.channel;
    ftdata      = ft_selectdata(cfg,ftdata); 
    conc_data   = concatenate_fttrials(ftdata);
    
    cfg = [];
    cfg.Fs      = conc_data.fsample;
    cfg.f_phase = cnfg.f_phase;
    cfg.f_amp   = cnfg.f_amp;
    cfg.Nsurro  = cnfg.Nsurro;
    
    CFD = cell(1,size(conc_data.trial{1},1));
    for ch = 1:size(conc_data.trial{1},1) 
        data_phase = conc_data.trial{1}(ch,:);
        CFD{ch} = computeCFD(cfg,data_phase);
        CFD{ch}.label_phase = conc_data.label{ch};
        CFD{ch}.label_amp   = conc_data.label{ch};
        % Include labels in the structure. Something else?
    end
end


% Case "interCFD"
if isfield(cnfg,'ch_interCFD') 
    cfg = [];
    cfg.latency = cnfg.latency;
    ftdata      = ft_selectdata(cfg,ftdata); 
    conc_data   = concatenate_fttrials(ftdata);
    
    cfg = [];
    cfg.Fs      = conc_data.fsample;
    cfg.f_phase = cnfg.f_phase;
    cfg.f_amp   = cnfg.f_amp;
    cfg.Nsurro = cnfg.Nsurro;
    
    CFD = cell(1,size(cnfg.ch_interCFD,1));
    for ch = 1:size(cnfg.ch_interCFD,1) 
        data_phase     = conc_data.trial{1}(cnfg.ch_interCFD(ch,1),:);
        data_amplitude = conc_data.trial{1}(cnfg.ch_interCFD(ch,2),:);
        CFD{ch} = computeCFD(cfg,data_phase,data_amplitude);
        CFD{ch}.label_phase = conc_data.label{cnfg.ch_interCFD(ch,1)};
        CFD{ch}.label_amp   = conc_data.label{cnfg.ch_interCFD(ch,2)};
    end    
end

%% Plot and save results

if cnfg.dosave
    if ~strcmp(cnfg.outpath(end),'\'), cnfg.outpath = [cnfg.outpath '\']; end
    if ~exist(cnfg.outpath,'dir'), mkdir(cnfg.outpath); end
end

if cnfg.doplot || cnfg.dosave    
    h=figure;
    nrows = ceil(sqrt(length(CFD)));
    ncols = ceil(length(CFD)/nrows);
    for iter=1:length(CFD)
        G(iter) = subplot(nrows,ncols,iter);
        G(iter).ButtonDownFcn = @newFigure1;
        [~,CFD{iter}]=plot_CFD_ft(CFD{iter});
        if strcmp(CFD{iter}.label_phase,CFD{iter}.label_amp)
            title(['CFD: ' CFD{iter}.label_phase]);            
        else
            title(['CFD: ' CFD{iter}.label_phase ' (ph) - ' CFD{iter}.label_amp ' (amp)']);
        end
    end
    
    if cnfg.dosave
        savefig(h,[cnfg.outpath 'CFD_' cnfg.infosave ])
        saveas(h,[cnfg.outpath 'CFD_' cnfg.infosave '.png'])
    end
    if ~cnfg.doplot, close(h), end
    
end

% Save data
if cnfg.dosave
    save([cnfg.outpath 'CFDval_' cnfg.infosave],'CFD')
end

end

 function newFigure1(h1,~)
%% Function to act on subplot click
%         Mouse click: Plots the selected subplot to a new figure
%  Ctrl + Mouse click: Delete subplot

    switch get(gcf,'SelectionType')
        case 'normal'
            F = figure();
            copyobj(h1.Children,gca(F));
            % Copy the selected subplot title
            tmp = get(h1,'title'); tmp = tmp.String;
            % Set title to the new figure
            title(gca(F), tmp);
            
            % Set axis to the new figure
            tmp = get(h1,'XLim');
            xlim(gca(F), tmp)
            tmp = get(h1,'YLim');
            ylim(gca(F), tmp)
        case 'alt'
            delete(h1);
    end
 end


