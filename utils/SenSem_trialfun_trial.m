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
%     trl = table(begsample,endsample,offset,condition,trial_idx,primer,target,primer_phoneme_nb,target_phoneme_nb,primer_lexical_freq,target_lexical_freq,primer_visual_freq,target_visual_freq);
end

%%
% trigger = [];
% for j = 1:length(event)
%     trigger = [trigger; event(j).value];
% end
% trigger = unique(trigger);
% if ~ or(all(trigger-20), all(trigger-30)) %if there are training trials
%     for j = 1:length(event)
%         if event(j).value < 10
%             event = event(j:end);
%             break;
%         end
%     end
% end
% 
% 
% % Remove bad trials
% min_delay = 0.3; % seconds
% max_delay = 1; % seconds
% trigger_val = [];
% for i = 1:length(event)
%     trigger_val = [trigger_val,event(i).value];
% end
% 
% trigger = unique(trigger_val);
% clc;
% fprintf('Triggers: ')
% fprintf('%d ', trigger)
% fprintf('\n')
% 
% % % idx = or(a==5,or(a==4,or(a==3,or(a==1,a==2))));
% idx = or(trigger_val==40,trigger_val==50); % Target onset
% diff = [];
% to_del = [];
% for i = drange(find(idx))
%     if i==1
%         continue
%     end
%     diff_sample = (event(i).sample - event(i-1).sample)/fs;
%     if diff_sample > max_delay
%         to_del = [to_del, i-2, i-1, i];
%     elseif diff_sample < min_delay
%         to_del = [to_del, i-2, i-1, i];
%     else
%         diff = [diff, diff_sample];
%     end
% end
% 
% new_event = event;
% new_event(to_del) = [];
% 
% % Sort trials by type (answer=yes/no; link=taxo/thema/none/neutral/filler)
% trigger_val = [];
% for i = 1:length(new_event)
%     trigger_val = [trigger_val,new_event(i).value];
% end
% 
% idx = trigger_val==40;
% idx_yes_thema = [];
% idx_yes_taxo = [];
% idx_yes_none = [];
% idx_yes_neutral = [];
% idx_yes_filler = [];
% 
% idx_thema_fab = [];
% idx_thema_nat = [];
% idx_taxo_fab = [];
% idx_taxo_nat = [];
% idx_none_fab = [];
% idx_none_nat = [];
% idx_neutral_fab = [];
% idx_neutral_nat = [];
% for i = drange(find(idx))
%     switch trigger_val(i-2)
%         case 1
%             idx_yes_thema = [idx_yes_thema, i-2, i-1, i];
%             if trigger_val(i-1) == 10
%                 idx_thema_fab = [idx_thema_fab, i-2, i-1, i];
%             elseif trigger_val(i-1) == 11
%                 idx_thema_nat= [idx_thema_nat, i-2, i-1, i];
%             end
%         case 2
%             idx_yes_taxo = [idx_yes_taxo, i-2, i-1, i];
%             if trigger_val(i-1) == 10
%                 idx_taxo_fab = [idx_taxo_fab, i-2, i-1, i];
%             elseif trigger_val(i-1) == 11
%                 idx_taxo_nat= [idx_taxo_nat, i-2, i-1, i];
%             end
%         case 3
%             idx_yes_none = [idx_yes_none, i-2, i-1, i];
%             if trigger_val(i-1) == 10
%                 idx_none_fab = [idx_none_fab, i-2, i-1, i];
%             elseif trigger_val(i-1) == 11
%                 idx_none_nat= [idx_none_nat, i-2, i-1, i];
%             end
%         case 4
%             idx_yes_neutral = [idx_yes_neutral, i-2, i-1, i];
%             if trigger_val(i-1) == 10
%                 idx_neutral_fab = [idx_neutral_fab, i-2, i-1, i];
%             elseif trigger_val(i-1) == 11
%                 idx_neutral_nat= [idx_neutral_nat, i-2, i-1, i];
%             end
%         case 5
%             idx_yes_filler = [idx_yes_filler, i-2, i-1, i];
%     end
% end
% idx = trigger_val==50;
% idx_no_filler = [];
% idx_no_any = [];
% for i = drange(find(idx))
%     switch trigger_val(i-2)
%         case 5
%             idx_no_filler = [idx_no_filler, i-2, i-1, i];
%         otherwise
%             idx_no_any= [idx_no_any, i-2, i-1, i];
%     end
% end
% 
% event_yes_thema = new_event(idx_yes_thema);
% event_yes_taxo = new_event(idx_yes_taxo);
% event_yes_none = new_event(idx_yes_none);
% event_yes_neutral = new_event(idx_yes_neutral);
% event_yes_filler = new_event(idx_yes_filler);
% event_no_filler = new_event(idx_no_filler);
% event_no_any = new_event(idx_no_any);
% 
% event_thema_fab = new_event(idx_thema_fab);
% event_thema_nat= new_event(idx_thema_nat);
% event_taxo_fab = new_event(idx_taxo_fab);
% event_taxo_nat= new_event(idx_taxo_nat);
% event_none_fab = new_event(idx_none_fab);
% event_none_nat= new_event(idx_none_nat);
% event_neutral_fab = new_event(idx_neutral_fab);
% event_neutral_nat= new_event(idx_neutral_nat);


end