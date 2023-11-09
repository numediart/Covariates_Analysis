function [trl, event] = SenSem_trialfun_trial(cfg)

% SenSem_trialfun_trial creates a trial definition that corresponds to the
% "Sensorimotor-Semantic ARC (SenSem)" project.
% The dataset consists of a priming task where the primer and the target
% are semantically linked or not. The type of the link between both stimuli
% defines the condition we want to express in the "trl" variable
% Conditions range from 0 to 8 as follow:
% 0: target = non-existing entity
% 1: thematic link between manufactured entities
% 2: thematic link between natural entities
% 3: taxonomic link between manufactured entities
% 4: taxonomic link between natural entities
% 5: No semantic link between manufactured entities
% 6: No semantic link between natural entities
% 7: primer = non-existing (neutral) entity, target = manufactured entities
% 8: primer = non-existing (neutral) entity, target = natural entities
% Note: if the subject answers too early/late to the requested task, the
% condition is put to NaN
%
% On top of that, some continuous variables describing each trial are
% extracted and referenced in trialinfo

if isfield(cfg, 'cov_description')
    var_file = cfg.cov_description;
else
    var_file = 'continuous_variables.csv'; % file contining variable information
end

non_exist_val = 12;
man_val       = 10;
nat_val       = 11;
thema_val     = 1;
taxo_val      = 2;
no_sem_val    = 3;
neutral_val   = 4;


% default file formats and chanindx for trigger detection
cfg.eventformat   = ft_getopt(cfg, 'eventformat');
cfg.headerformat  = ft_getopt(cfg, 'headerformat');
cfg.dataformat    = ft_getopt(cfg, 'dataformat');

% get the header, among others for the sampling frequency
if isfield(cfg, 'hdr')
  ft_info('using the header from the configuration structure\n');
  hdr = cfg.hdr;
else
  % read the header, contains the sampling frequency
  ft_info('reading the header from ''%s''\n', cfg.headerfile);
  hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
end

if isfield(cfg, 'event')
% for BCI applications events should be specified in the cfg
% to prevent reading the same events many times
    event = cfg.event;
else
    event = ft_read_event(cfg.dataset);
end

% start by selecting all events
sel = true(1, length(event)); % this should be a row vector

% select all events of the specified type
if isfield(cfg.trialdef, 'eventtype') && ~isempty(cfg.trialdef.eventtype)
  for i=1:numel(event)
    sel(i) = sel(i) && strcmp(event(i).type, cfg.trialdef.eventtype);
  end
elseif ~isfield(cfg.trialdef, 'eventtype') || isempty(cfg.trialdef.eventtype)
  % search for STATUS events
  for i=1:numel(event)
    sel(i) = sel(i) && strcmp(event(i).type, 'STATUS');
  end
end

% convert from boolean vector into a list of indices
sel = find(sel);

if isfield(cfg.trialdef, 'eventvalue') && ~isempty(cfg.trialdef.eventvalue)
    sel = sel(ismember([event(sel).value],cfg.trialdef.eventvalue));
end

trl = [];
cond = [];
for i = sel
    if ~isfield(cfg.trialdef, 'prestim')
        trloff = event(i).offset;
        trlbeg = event(i).sample;
    else
        % override the offset of the event
        trloff = round(-cfg.trialdef.prestim*hdr.Fs);
        % also shift the begin sample with the specified amount
        trlbeg = event(i).sample + trloff;
    end
    % determine the number of samples that has to be read (excluding the begin sample)
    if ~isfield(cfg.trialdef, 'poststim')
        trldur = max(event(i).duration - 1, 0);
    else
        % this will not work if prestim was not defined, the code will then crash
        try
            trldur = round((cfg.trialdef.poststim+cfg.trialdef.prestim)*hdr.Fs) - 1;
        catch
            error('a trial prestimulus has to be defined as cfg.trialdef.prestim');
        end
    end
    trlend = trlbeg + trldur;
    
    % add the beginsample, endsample and offset of this trial to the list
    % if all samples are in the dataset
    if trlbeg>0 && trlend<=hdr.nSamples*hdr.nTrials
        trl = [trl; [trlbeg trlend trloff]];
        if isnumeric(event(i).value)
            % Detect bad trials
            min_delay = 0.3; % seconds
            max_delay = 1; % seconds
            diff = (event(i+1).sample - event(i).sample)/hdr.Fs;
            if diff<min_delay || diff>max_delay
                cond = [cond; nan];
            elseif event(i).value == non_exist_val
                cond = [cond; 0];
            elseif event(i).value == man_val && event(i-1).value == thema_val
                cond = [cond; 1];
            elseif event(i).value == nat_val && event(i-1).value == thema_val
                cond = [cond; 2];
            elseif event(i).value == man_val && event(i-1).value == taxo_val
                cond = [cond; 3];
            elseif event(i).value == nat_val && event(i-1).value == taxo_val
                cond = [cond; 4];
            elseif event(i).value == man_val && event(i-1).value == no_sem_val
                cond = [cond; 5];
            elseif event(i).value == nat_val && event(i-1).value == no_sem_val
                cond = [cond; 6];
            elseif event(i).value == man_val && event(i-1).value == neutral_val
                cond = [cond; 7];
            elseif event(i).value == nat_val && event(i-1).value == neutral_val
                cond = [cond; 8];
            else
                sprintf("unknown combination %d-%d \n",event(i).value,event(i-1).value)
                cond = [cond; nan];
            end
        else
            sprintf("Non nummeric value\n")
            cond = [cond; nan];   
        end
    end
end

% read txt file
txt_file = fileread([cfg.dataset(1:end-3) 'txt']);
% myfile = 'D:\__EEG-data\BIDS_files\sub-030\eeg\sub-030_task-semantic-priming_eeg.txt';
% txt_file = fileread(myfile);
txt_file = double(txt_file);
txt_file(txt_file==0) = [];
txt_file = char(txt_file);
first_trial = strfind(txt_file,'TrialProc');
first_trial = first_trial(1);
txt_file = txt_file(first_trial:end);

% trial_idx_pat = 'TrialList: ';
primer_pat = 'amorce: ';
% primer_phoneme_pat = 'nbPhonemeAmorce: ';
% primer_lexical_pat = 'freqLivreAmorce: ';
% primer_visual_pat = 'freqFilmAmorce: ';
target_pat = 'cible: ';
% target_phoneme_pat = 'nbPhoneme: ';
% target_lexical_pat = 'freqLivre: ';
% target_visual_pat = 'freqFilm: ';

if isfield(cfg, "cov_description")
    T = readtable(var_file);
else
    T = readtable(fullfile(cfg.root,var_file));
end
% writetable(T,'continuous_variables.csv')
var_names = T.Properties.VariableNames;

% trial_idx = str2double(extractBetween(txt_file,trial_idx_pat,char(13)));
primer = convertCharsToStrings(extractBetween(txt_file,primer_pat,char(13)));
% primer_phoneme_nb = str2double(extractBetween(txt_file,primer_phoneme_pat,char(13)));
% primer_lexical_freq = str2double(extractBetween(txt_file,primer_lexical_pat,char(13)));
% primer_visual_freq = str2double(extractBetween(txt_file,primer_visual_pat,char(13)));
target = convertCharsToStrings(extractBetween(txt_file,target_pat,char(13)));
% target_phoneme_nb = str2double(extractBetween(txt_file,target_phoneme_pat,char(13)));
% target_lexical_freq = str2double(extractBetween(txt_file,target_lexical_pat,char(13)));
% target_visual_freq = str2double(extractBetween(txt_file,target_visual_pat,char(13)));

% append all the vectors to trl

if ~isempty(cond) && ~all(isnan(cond)) && length(trl)==length(cond) && length(trl)==length(primer) && length(trl)==length(target)
    trl_table.begsample = trl(:,1);
    trl_table.endsample = trl(:,2);
    trl_table.offset = trl(:,3);
    trl_table.condition = cond;
    trl_table.primer = primer;
    trl_table.target = target;
    
    p_idx = [];
    t_idx = [];
    for i = 1:length(primer)
        p_i = find(strcmpi(T.(var_names{1}),primer{i}));
        t_i = find(strcmpi(T.(var_names{1}),target{i}));
        if ~isempty(p_i) && ~isempty(t_i)
            p_idx = [p_idx,p_i];
            t_idx = [t_idx,t_i];
        elseif isempty(p_i) && ~isempty(t_i)
            p_idx = [p_idx,height(T)+1];
            t_idx = [t_idx,t_i];
        elseif ~isempty(p_i) && isempty(t_i)
            p_idx = [p_idx,p_i];
            t_idx = [t_idx,height(T)+1];
        else
            p_idx = [p_idx,height(T)+1];
            t_idx = [t_idx,height(T)+1];
        end
    end
    
    last_row = array2table(NaN(1,width(T)),'VariableNames',var_names);
    T = [T;last_row];
    
    for var = {var_names{4:end}}
        var = var{1};
        tmp = genvarname(['primer_' var]);
        trl_table.(tmp) = T.(var)(p_idx);
        tmp = genvarname(['target_' var]);
        trl_table.(tmp) = T.(var)(t_idx);
    end
    trl = struct2table(trl_table);
end

end