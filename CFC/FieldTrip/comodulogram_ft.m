function CFC=comodulogram_ft(cnfg,ftdata)

%% Compute phase-amplitude CFC with fieldtrip structures 
% It computes the phase and amplitude with filters and Hilbert transform.
%
% USE:
%   CFC=comodulogram_ft(cfg,ft_data);
%
% PARAMETERS (cfg.): 
%   
%   latency - [t_start t_end]
%   channel - (Def='all') Channel to compute intraCFC
%   
% If you want to compute CFC between phases and amplitudes of different
% signals, it is possible to define a vector of channels with the
% specific pairs (ch_interCFC). For example: ch_interCFC = [1 1; 1 2; 1 3]
% will compute the CFC between the phase of ch1 and the amplitudes of
% channels 1, 2 and 3. If defined, cfg.channel won't be used.
%
%   ch_interCFC (N x 2)- Two-column vector with pairs of channels for CFC.                       
%
% For CFC analysis:
%
%   bins     - Num  ber of bins to divide each theta cycle (Def=16)
%   f_phase  - Information to make the frequency phase vector
%                 [f_min, f_max, f_step, BW] 
%   f_amp    - Information to make the frequency amplitude vector
%                 [f_min, f_max, f_step, BW] 
%   Nsurro   - Number of surrogates to statistical significance (Def=0)
%
% To save the data:
%
%   doplot   - true/false (Def: true)
%   dosave   - true/false (Def: false)
%   outpath  - string with path 
%   infosave - string with filename
%
% OUTPUT:
%   comodulogram = vector of cells:
%       Fs: Sampling Frequency (Hz).
%       bins: Number of bins.
%       f_theta: Information about the frequency theta vector.
%       f_gamma: Information about the frequency gamma vector.
%       MI: (N_gamma_pixels, N_theta_pixels) Matrix with the Modulation.
%           Index of each gamma-theta frequency.
%       CFC: (N_gamma_pixels, N_theta_pixels, 2*bins) Matrix with the
%            average of gamma amplitue during a cycle for each 
%            gamma-theta frequency.
%
% See also: computeCFC modulation_index plot_comodulogram_ft

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Aug 21; Last revision: 16-Aug-2023

% Change log:
% 06/09/23 - Changed pval as double instead of char
% 16/08/23 - Modified to store statistical mask too
% 01/06/23 - Included plot_comodulogram stats as parameters

%% Initial parameters:

if ~isfield(cnfg,'latency')
    cnfg.latency=[ftdata.time{1}(1) ftdata.time{1}(end)]; end
if isfield(cnfg,'ch_interCFC') && isfield(cnfg,'channel')
    if size(cnfg.ch_interCFC,2) ~= 2 %Should have two columns
        error('cnfg.ch_interCFC must have two columns'), end
    warning('ch_interCFC was defined. Information in cfg.channel won''t be used.') 
    cnfg=rmfield(cnfg,'channel'); 
end
if ~isfield(cnfg,'ch_interCFC') && ~isfield(cnfg,'channel')
    cnfg.channel = 'all'; end
if ~isfield(cnfg,'bins')
    cnfg.bins = 16; end
if ~isfield(cnfg,'f_phase')
    cnfg.f_phase = [3 16 1 2]; end
if ~isfield(cnfg,'f_amp')
    cnfg.f_amp = [20 120 5 20]; end
if ~isfield(cnfg,'stats'), cnfg.stats = 'cluster'; end
if ~isfield(cnfg,'pval'), cnfg.pval = 0.01; end
if ~isfield(cnfg,'pval_clus'), cnfg.pval_clus = 0.01; end
if ~isfield(cnfg,'Nsurro'), cnfg.Nsurro = 0; end
if ~isfield(cnfg,'dosave'), cnfg.dosave=false; end
if ~isfield(cnfg,'doplot'), cnfg.doplot=true; end
if isfield(cnfg,'outpath')
    if ~strcmp(cnfg.outpath(end),'\')
        cnfg.outpath = [cnfg.outpath '\'];
    end
end
if ~isfield(cnfg,'infosave'), cnfg.infosave=''; end

%% COMPUTE CFC

% Case "intraCFC"
if isfield(cnfg,'channel')
    cfg = [];
    cfg.latency = cnfg.latency;
    cfg.channel = cnfg.channel;
    ftdata      = ft_selectdata(cfg,ftdata); 
    conc_data   = concatenate_fttrials(ftdata);
    
    cfg = [];
    cfg.Fs      = conc_data.fsample;
    cfg.bins    = cnfg.bins;
    cfg.f_phase = cnfg.f_phase;
    cfg.f_amp   = cnfg.f_amp;
    cfg.Nsurro  = cnfg.Nsurro;
    
    CFC = cell(1,size(conc_data.trial{1},1));
    for ch = 1:size(conc_data.trial{1},1) 
        data_phase = conc_data.trial{1}(ch,:);
        CFC{ch} = computeCFC(cfg,data_phase);
        CFC{ch}.label_phase = conc_data.label{ch};
        CFC{ch}.label_amp   = conc_data.label{ch};
        % Include labels in the structure. Something else?
    end
end


% Case "interCFC"
if isfield(cnfg,'ch_interCFC') 
    cfg = [];
    cfg.latency = cnfg.latency;
    ftdata      = ft_selectdata(cfg,ftdata); 
    conc_data   = concatenate_fttrials(ftdata);
    
    cfg = [];
    cfg.Fs      = conc_data.fsample;
    cfg.bins    = cnfg.bins;
    cfg.f_phase = cnfg.f_phase;
    cfg.f_amp   = cnfg.f_amp;
    cfg.Nsurro = cnfg.Nsurro;
    
    CFC = cell(1,size(cnfg.ch_interCFC,1));
    for ch = 1:size(cnfg.ch_interCFC,1) 
        data_phase     = conc_data.trial{1}(cnfg.ch_interCFC(ch,1),:);
        data_amplitude = conc_data.trial{1}(cnfg.ch_interCFC(ch,2),:);
        CFC{ch} = computeCFC(cfg,data_phase,data_amplitude);
        CFC{ch}.label_phase = conc_data.label{cnfg.ch_interCFC(ch,1)};
        CFC{ch}.label_amp   = conc_data.label{cnfg.ch_interCFC(ch,2)};
    end    
end

%% Plot and save results

if cnfg.dosave
    if ~strcmp(cnfg.outpath(end),'\'), cnfg.outpath = [cnfg.outpath '\']; end
    if ~exist(cnfg.outpath,'dir'), mkdir(cnfg.outpath); end
end

if cnfg.doplot || cnfg.dosave    
    h=figure;
    nrows = ceil(sqrt(length(CFC)));
    ncols = ceil(length(CFC)/nrows);
    for iter=1:length(CFC)
        G(iter) = subplot(nrows,ncols,iter);
        G(iter).ButtonDownFcn = @newFigure1;
        cfg = [];
        cfg.stats     = cnfg.stats;
        cfg.pval      = cnfg.pval;
        cfg.pval_clus = cnfg.pval_clus;
        [~,CFC{iter}]=plot_comodulogram_ft(CFC{iter},cfg);
        if strcmp(CFC{iter}.label_phase,CFC{iter}.label_amp)
            title(['CFC: ' CFC{iter}.label_phase]);            
        else
            title(['CFC: ' CFC{iter}.label_phase ' (ph) - ' CFC{iter}.label_amp ' (amp)']);
        end
    end
    
    if cnfg.dosave
        savefig(h,[cnfg.outpath 'CFC_' cnfg.infosave ])
        saveas(h,[cnfg.outpath 'CFC_' cnfg.infosave '.png'])
    end
    if ~cnfg.doplot, close(h), end
    
end

% Save data
if cnfg.dosave
    save([cnfg.outpath 'CFCval_' cnfg.infosave],'CFC')
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


