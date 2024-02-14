function conc_data = concatenate_fttrials(ftdata)
%% Combines all the trials into a single (long) trial
% 
% Syntax:  
%    conc_data = concatenate_fttrials(ftdata) concatenates the trials
%       of the Fildtrip data producing a single trial.
%
% Inputs:
%    ftdata - Fildtrip data obtained by ft_preprocessing.
%
% Outputs:
%    conc_data - concatenated Fildtrip data with a single (long) trial.
%
% Tips:
%    Could be used for concatenating the trials to compute a covariance
%    matrix based on raw (not average) data.
%
% Example: 
%    ftdata = ft_preprocessing( ... );
%    conc_data = concatenate_fttrials(ftdata);
%    
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ft_preprocessing,  ft_redefinetrial,  ft_selectdata

% Author: Christos Papageorgakis
% email: christos.papageorgakis@gmail.com
% Jan. 2019; Last revision: 16-Jan-2019

% Make a FieldTrip data template for data with one trials
cfg = {};
cfg.trials = 1;
conc_data = ft_preprocessing(cfg, ftdata);

% Concatenate all the trials to a single trials
conc_data.trial = {horzcat(ftdata.trial{1:end})};

% keep the sample-info of the 1st and last trial
if isfield(ftdata, 'sampleinfo')
    warning('Sampleinfo correspondence to raw is lost. DO NOT use conc_data.sampleinfo');
    conc_data.sampleinfo = [ftdata.sampleinfo(1,1), ftdata.sampleinfo(end,end)];
end

% Make a plausible time axis for the single trial
warning('Time correspondence to raw is lost. DO NOT use conc_data.time');
conc_data.time = {(1:size(conc_data.trial{1},2)) / conc_data.fsample};

end