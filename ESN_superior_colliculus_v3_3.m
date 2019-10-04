%% CLEAR 
clc;
clear;
close all;

%% Load data 
[file_name,file_path] = uigetfile([pwd filesep '*.mat'], 'Select MAT file');
fprintf(['Loading ', file_name, ' ... ']);
load([file_path filesep file_name]);
fprintf(' --> Completed. \n')

%% Extract variables from CH_DATA_ALL 
clearvars -except CH_DATA_ALL SMR_FILE DATA;
clearvars DATA;
fprintf(['Extracting variables ', SMR_FILE.file_name, ' ... ']);
CH_DATA_ALL_cell = struct2cell(CH_DATA_ALL);
SMR_FILE.field_names = CH_DATA_ALL_cell(4,:);
clearvars CH_DATA_ALL_cell;
DATA.NEURAL.waveform_time = (0 : CH_DATA_ALL(ismember(SMR_FILE.field_names,'Unit')).interval : CH_DATA_ALL(ismember(SMR_FILE.field_names,'Unit')).max_time)'; % 'Unit'
DATA.NEURAL.waveform      = CH_DATA_ALL(ismember(SMR_FILE.field_names,'Unit')).values; % 'Unit'
DATA.NEURAL.spike_time    = CH_DATA_ALL(ismember(SMR_FILE.field_names,'Spike')).values_times; % 'Spike'
DATA.NEURAL.spike_firing  = CH_DATA_ALL(ismember(SMR_FILE.field_names,'Spike Den')).values; % 'Spike Den'
DATA.TRIAL.eye_time       = (0 : CH_DATA_ALL(ismember(SMR_FILE.field_names,'HE')).interval : CH_DATA_ALL(ismember(SMR_FILE.field_names,'HE')).max_time)'; % 'HE'
DATA.TRIAL.eye_px         = CH_DATA_ALL(ismember(SMR_FILE.field_names,'HE')).values; % 'HE'
DATA.TRIAL.eye_py         = CH_DATA_ALL(ismember(SMR_FILE.field_names,'VE')).values; % 'VE'
DATA.TRIAL.tgt_px         = CH_DATA_ALL(ismember(SMR_FILE.field_names,'HT')).values; % 'HT'
DATA.TRIAL.tgt_py         = CH_DATA_ALL(ismember(SMR_FILE.field_names,'VT')).values; % 'VT'

% filter params
sampling_freq = 1 / (CH_DATA_ALL(ismember(SMR_FILE.field_names,'HE')).interval); % 'HE'
cutoff_freq = 100.0;
[b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
% filter eye_r data
DATA.TRIAL.eye_px_filt = filtfilt(b_butter,a_butter,DATA.TRIAL.eye_px);
DATA.TRIAL.eye_py_filt = filtfilt(b_butter,a_butter,DATA.TRIAL.eye_py);
DATA.TRIAL.eye_vx_filt = diff(DATA.TRIAL.eye_px_filt)./diff(DATA.TRIAL.eye_time); DATA.TRIAL.eye_vx_filt=[DATA.TRIAL.eye_vx_filt(1); DATA.TRIAL.eye_vx_filt];
DATA.TRIAL.eye_vy_filt = diff(DATA.TRIAL.eye_py_filt)./diff(DATA.TRIAL.eye_time); DATA.TRIAL.eye_vy_filt=[DATA.TRIAL.eye_vy_filt(1); DATA.TRIAL.eye_vy_filt];
DATA.TRIAL.eye_vm_filt = sqrt(DATA.TRIAL.eye_vx_filt.^2 + DATA.TRIAL.eye_vy_filt.^2);

eye_time            = DATA.TRIAL.eye_time;
length_time         = length(eye_time);
spike_time          = DATA.NEURAL.spike_time;
sac_str_time        = CH_DATA_ALL(ismember(SMR_FILE.field_names,'Saccade O')).values_times; % 'Saccade O'
sac_end_time        = CH_DATA_ALL(ismember(SMR_FILE.field_names,'Saccade E')).values_times; % 'Saccade E'
tgt_jmp_time        = CH_DATA_ALL(ismember(SMR_FILE.field_names,'Target Ti')).values_times; % 'Target Ti'
spike_time(  end+1) = max([eye_time(end), spike_time(end)])+1;
sac_str_time(end+1) = max([eye_time(end), sac_str_time(end)])+1;
sac_end_time(end+1) = max([eye_time(end), sac_end_time(end)])+1;
tgt_jmp_time(end+1) = max([eye_time(end), tgt_jmp_time(end)])+1;
spike_event         = false(length_time, 1);
sac_str_event       = false(length_time, 1);
sac_end_event       = false(length_time, 1);
tgt_jmp_event       = false(length_time, 1);
counter_spike_event   = find(spike_time    > eye_time(1), 1, 'first');
counter_sac_str_event = find(sac_str_time  > eye_time(1), 1, 'first');
counter_sac_end_event = find(sac_end_time  > eye_time(1), 1, 'first');
counter_tgt_jmp_event = find(tgt_jmp_time  > eye_time(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = eye_time(counter_time_point);
    if time_ponit_ >= spike_time(  counter_spike_event)
        spike_event(  counter_time_point) = true;
        counter_spike_event   = counter_spike_event + 1;
    end
    if time_ponit_ >= sac_str_time(counter_sac_str_event)
        sac_str_event(counter_time_point) = true;
        counter_sac_str_event = counter_sac_str_event + 1;
    end
    if time_ponit_ >= sac_end_time(counter_sac_end_event)
        sac_end_event(counter_time_point) = true;
        counter_sac_end_event = counter_sac_end_event + 1;
    end
    if time_ponit_ >= tgt_jmp_time(counter_tgt_jmp_event)
        tgt_jmp_event(counter_time_point) = true;
        counter_tgt_jmp_event = counter_tgt_jmp_event + 1;
    end
end
DATA.TRIAL.spike_event   = spike_event;
DATA.TRIAL.sac_str_event = sac_str_event;
DATA.TRIAL.sac_end_event = sac_end_event;
DATA.TRIAL.tgt_jmp_event = tgt_jmp_event;
fprintf(' --> Completed. \n')

%% Calc primary & corrective saccade params from DATA.TRIAL 
clearvars -except CH_DATA_ALL SMR_FILE DATA;
fprintf(['Calculating Prim&Corr ', SMR_FILE.file_name, ' ... ']);

% Finding params based on the histogram of the diff(tgt)
% if the file name does NOT contain SquirrelUnit
if ~contains(SMR_FILE.file_name,'SquirrelUnit')
tgt_diff_px = diff([DATA.TRIAL.tgt_px; DATA.TRIAL.tgt_px(end)]);
tgt_diff_py = diff([DATA.TRIAL.tgt_py; DATA.TRIAL.tgt_py(end)]);
tgt_diff_pm = sqrt(tgt_diff_px.^2 + tgt_diff_py.^2);
tgt_diff_angle = atan2d(tgt_diff_py, tgt_diff_px);
tgt_diff_angle = ESN_Round(tgt_diff_angle, 0.5);
inds_fixation = tgt_diff_pm < 0.01;
tgt_diff_pm(inds_fixation) = nan;
inds_pursuit = tgt_diff_pm < 1.5;
tgt_diff_pm(inds_pursuit) = nan;
tgt_diff_pm = ESN_Round(tgt_diff_pm, 0.1);
end

% Finding params by trusting the "Target Ti" in the CH_DATA_ALL
% if the file name contain SquirrelUnit
if contains(SMR_FILE.file_name,'SquirrelUnit')
tgt_jmp_event = find(DATA.TRIAL.tgt_jmp_event);
window_length_ = 10;
tgt_px = DATA.TRIAL.tgt_px;
tgt_py = DATA.TRIAL.tgt_py;
tgt_jmp_event((tgt_jmp_event-window_length_)<1)              = 1; % make sure all inds are greater than 1
tgt_jmp_event((tgt_jmp_event+window_length_)>length(tgt_px)) = length(tgt_px); % make sure all inds are within the bound
tgt_jmp_px = tgt_px(tgt_jmp_event+window_length_) - tgt_px(tgt_jmp_event-window_length_);
tgt_jmp_py = tgt_py(tgt_jmp_event+window_length_) - tgt_py(tgt_jmp_event-window_length_);
tgt_jmp_pm = sqrt(tgt_jmp_px.^2 + tgt_jmp_py.^2);
tgt_diff_pm = nan(size(tgt_px));
tgt_diff_pm(tgt_jmp_event) = tgt_jmp_pm;

tgt_jmp_angle = atan2d(tgt_jmp_py, tgt_jmp_px);
tgt_jmp_angle = ESN_Round(tgt_jmp_angle, 0.5);
tgt_diff_angle = nan(size(tgt_px));
tgt_diff_angle(tgt_jmp_event) = tgt_jmp_angle;
end

edges_tgt_diff_pm = -0.5 : 1 : 20.5;
[N_,~] = histcounts(tgt_diff_pm, edges_tgt_diff_pm);
% histogram('BinEdges',edges_tgt_diff_pm,'BinCounts',N_) % Check the
% % histogram to make sure there is prominent peak for prim_sac_amp
[~, ind_max_] = max(N_);
tgt_prim_sac_ind = (tgt_diff_pm > edges_tgt_diff_pm(ind_max_)) & (tgt_diff_pm < edges_tgt_diff_pm(ind_max_+1));

tgt_diff_pm(tgt_prim_sac_ind) = nan;

edges_tgt_diff_angle = -182.5 : 5 : 182.5;
[N_,~] = histcounts(tgt_diff_angle(~isnan(tgt_diff_pm)), edges_tgt_diff_angle);
% histogram('BinEdges',edges_tgt_diff_angle,'BinCounts',N_) % Check the
% % histogram to make sure there is prominent peak for corr_sac_ang
[~, ind_max_] = max(N_);
tgt_corr_sac_ind = (tgt_diff_angle > edges_tgt_diff_angle(ind_max_)) & (tgt_diff_angle < edges_tgt_diff_angle(ind_max_+1)) & (~isnan(tgt_diff_pm));
tgt_corr_sac_ang = tgt_diff_angle(tgt_corr_sac_ind);
if (ESN_Round(mean(tgt_corr_sac_ang), 0.5) >= 0 )
    tgt_corr_sac_ind = (tgt_diff_angle >= 0) & (~isnan(tgt_diff_pm)) ;
else
    tgt_corr_sac_ind = (tgt_diff_angle <  0) & (~isnan(tgt_diff_pm)) ;
end

tgt_prim_sac_ind_ = find(tgt_prim_sac_ind);
tgt_corr_sac_ind_ = find(tgt_corr_sac_ind);
num_trials_       = length(tgt_prim_sac_ind_);
tgt_corr_sac_ind_match_ = nan(num_trials_, 1);
tgt_prob_sac_ind_match_ = nan(num_trials_, 1);
validity_prim_sac_      = true( num_trials_, 1);
validity_corr_sac_      = false(num_trials_, 1);
validity_prob_sac_      = false(num_trials_, 1);

tgt_prim_sac_ind_(end+1)   = max([tgt_prim_sac_ind_(end), tgt_corr_sac_ind_(end)])+1;
for counter_trial = 1 : 1 : num_trials_
    ind_current_prim_ = tgt_prim_sac_ind_(counter_trial);
    ind_next_prim_    = tgt_prim_sac_ind_(counter_trial+1);
    ind_current_corr_  = tgt_corr_sac_ind_( (tgt_corr_sac_ind_>ind_current_prim_) & (tgt_corr_sac_ind_<ind_next_prim_) );
    if isempty(ind_current_corr_)
        validity_corr_sac_(counter_trial, 1) = false;
        validity_prob_sac_(counter_trial, 1) = false;
        tgt_corr_sac_ind_match_(counter_trial, 1) = ind_current_prim_;
        tgt_prob_sac_ind_match_(counter_trial, 1) = ind_current_prim_;
    elseif length(ind_current_corr_) == 1
        validity_corr_sac_(counter_trial, 1) = true;
        validity_prob_sac_(counter_trial, 1) = false;
        tgt_corr_sac_ind_match_(counter_trial, 1) = ind_current_corr_(1);
        tgt_prob_sac_ind_match_(counter_trial, 1) = ind_current_corr_(1);
    elseif length(ind_current_corr_) > 1
        validity_corr_sac_(counter_trial, 1) = true;
        validity_prob_sac_(counter_trial, 1) = true;
        tgt_corr_sac_ind_match_(counter_trial, 1) = ind_current_corr_(1);
        tgt_prob_sac_ind_match_(counter_trial, 1) = ind_current_corr_(2);
    end
end

tgt_prim_sac_ind        = find(tgt_prim_sac_ind);
tgt_corr_sac_ind        = tgt_corr_sac_ind_match_;
tgt_prob_sac_ind        = tgt_prob_sac_ind_match_;

tgt_prim_sac_validity   = validity_prim_sac_;
tgt_corr_sac_validity   = validity_corr_sac_;
tgt_prob_sac_validity   = validity_prob_sac_;

DATA.PRIM_SAC.ind_tgt_jump   = tgt_prim_sac_ind;
DATA.PRIM_SAC.tgt_amp_m = tgt_diff_pm(   tgt_prim_sac_ind);
DATA.PRIM_SAC.tgt_ang   = tgt_diff_angle(tgt_prim_sac_ind);
DATA.PRIM_SAC.validity  = tgt_prim_sac_validity;

DATA.CORR_SAC.ind_tgt_jump   = tgt_corr_sac_ind;
DATA.CORR_SAC.tgt_amp_m = tgt_diff_pm(   tgt_corr_sac_ind);
DATA.CORR_SAC.tgt_ang   = tgt_diff_angle(tgt_corr_sac_ind);
DATA.CORR_SAC.validity  = tgt_corr_sac_validity;

DATA.PROB_SAC.ind_tgt_jump   = tgt_prob_sac_ind;
DATA.PROB_SAC.tgt_amp_m = tgt_diff_pm(   tgt_prob_sac_ind);
DATA.PROB_SAC.tgt_ang   = tgt_diff_angle(tgt_prob_sac_ind);
DATA.PROB_SAC.validity  = tgt_prob_sac_validity;

fprintf(' --> Completed. \n')

%% rotate the eye and tgt data 
clearvars -except CH_DATA_ALL SMR_FILE DATA;
fprintf(['Rotating data ', SMR_FILE.file_name, ' ... ']);
eye_time = DATA.TRIAL.eye_time;
eye_px   = DATA.TRIAL.eye_px;
eye_py   = DATA.TRIAL.eye_py;
tgt_px   = DATA.TRIAL.tgt_px;
tgt_py   = DATA.TRIAL.tgt_py;

tgt_corr_sac_ang = DATA.CORR_SAC.tgt_ang;
main_dir = ESN_Round(mean(tgt_corr_sac_ang), 0.5);

rotation_matrix = [cosd(main_dir) -sind(main_dir); sind(main_dir) cosd(main_dir)];
eye_rot = [eye_px eye_py]*rotation_matrix;
tgt_rot = [tgt_px tgt_py]*rotation_matrix;
eye_px   = eye_rot(:,1);
eye_py   = eye_rot(:,2);
tgt_px   = tgt_rot(:,1);
tgt_py   = tgt_rot(:,2);

% filter params
sampling_freq = 1 / (CH_DATA_ALL(ismember(SMR_FILE.field_names,'HE')).interval); % HE
cutoff_freq = 100.0;
[b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
% filter eye_r data
eye_px_filt = filtfilt(b_butter,a_butter,eye_px);
eye_py_filt = filtfilt(b_butter,a_butter,eye_py);
eye_vx_filt = diff(eye_px_filt)./diff(eye_time); eye_vx_filt=[eye_vx_filt(1); eye_vx_filt];
eye_vy_filt = diff(eye_py_filt)./diff(eye_time); eye_vy_filt=[eye_vy_filt(1); eye_vy_filt];
eye_vm_filt = sqrt(eye_vx_filt.^2 + eye_vy_filt.^2);

DATA.TRIAL_ROT.main_dir      = main_dir;
DATA.TRIAL_ROT.eye_px        = eye_px;
DATA.TRIAL_ROT.eye_py        = eye_py;
DATA.TRIAL_ROT.eye_px_filt   = eye_px_filt;
DATA.TRIAL_ROT.eye_py_filt   = eye_py_filt;
DATA.TRIAL_ROT.eye_vx_filt   = eye_vx_filt;
DATA.TRIAL_ROT.eye_vy_filt   = eye_vy_filt;
DATA.TRIAL_ROT.eye_vm_filt   = eye_vm_filt;
DATA.TRIAL_ROT.tgt_px        = tgt_px;
DATA.TRIAL_ROT.tgt_py        = tgt_py;
fprintf(' --> Completed. \n')

%% Calc PRIM_SAC, CORR_SAC, PROB_SAC
clearvars -except CH_DATA_ALL SMR_FILE DATA;
fprintf(['Calculating Prim Sac ', SMR_FILE.file_name, ' ... ']);
% get all saccade data and label them prim, corr, prob based on the time of
% their occurance
ind_sac_str_event = find(DATA.TRIAL.sac_str_event);
ind_sac_end_event = find(DATA.TRIAL.sac_end_event);
% take care of the error that number of ind_sac_str_event is not the same
% as ind_sac_end_event
if length(ind_sac_str_event) ~= length(ind_sac_end_event)
    if length(ind_sac_str_event) < length(ind_sac_end_event)
        num_trials_ = length(ind_sac_str_event);
        condition_ = 1;
    else
        num_trials_ = length(ind_sac_end_event);
        condition_ = 2;
    end
    ind_sac_str_event_ = nan(num_trials_, 1);
    ind_sac_end_event_ = nan(num_trials_, 1);
    
    for counter_trial_ = 1 : 1 : num_trials_
        if condition_ == 1
            ind_sac_str_event_(counter_trial_) = ind_sac_str_event(counter_trial_);
            ind_sac_end_event_(counter_trial_) = ind_sac_end_event(find(ind_sac_end_event > ind_sac_str_event(counter_trial_), 1, 'first'));
        end
        if condition_ == 2
            ind_sac_end_event_(counter_trial_) = ind_sac_end_event(counter_trial_);
            ind_sac_str_event_(counter_trial_) = ind_sac_str_event(find(ind_sac_str_event < ind_sac_end_event(counter_trial_), 1, 'last'));
        end
    end
    ind_sac_str_event = ind_sac_str_event_;
    ind_sac_end_event = ind_sac_end_event_;
end

ind_tgt_jump_prim = DATA.PRIM_SAC.ind_tgt_jump;
ind_tgt_jump_corr = DATA.CORR_SAC.ind_tgt_jump;
ind_tgt_jump_prob = DATA.PROB_SAC.ind_tgt_jump;
num_trials   = length(ind_tgt_jump_prim);
validity_prim          = DATA.PRIM_SAC.validity;
validity_corr          = DATA.CORR_SAC.validity;
validity_prob          = DATA.PROB_SAC.validity;
ind_str_prim           = nan(  num_trials, 1);
ind_str_corr           = nan(  num_trials, 1);
ind_str_prob           = nan(  num_trials, 1);
ind_end_prim           = nan(  num_trials, 1);
ind_end_corr           = nan(  num_trials, 1);
ind_end_prob           = nan(  num_trials, 1);
ind_tgt_jump_prim(end+1)   = max([ind_tgt_jump_prim(end), ind_sac_end_event(end)])+1;
ind_tgt_jump_corr(end+1)   = max([ind_tgt_jump_corr(end), ind_sac_end_event(end)])+1;
ind_tgt_jump_prob(end+1)   = max([ind_tgt_jump_prob(end), ind_sac_end_event(end)])+1;
for counter_trial = 1 : 1 : num_trials
    flag_overlap = false;
    ind_current_prim_ = ind_tgt_jump_prim(counter_trial);
    ind_next_prim_    = ind_tgt_jump_prim(counter_trial+1);
    ind_current_corr_ = ind_tgt_jump_corr(counter_trial);
    ind_next_corr_    = ind_tgt_jump_corr(counter_trial+1);
    ind_current_prob_ = ind_tgt_jump_prob(counter_trial);
    ind_next_prob_    = ind_tgt_jump_prob(counter_trial+1);
    
    ind_sac_str_event_  = ind_sac_str_event( (ind_sac_str_event>ind_current_prim_) & (ind_sac_str_event<ind_next_prim_) );
    ind_sac_end_event_  = ind_sac_end_event( (ind_sac_str_event>ind_current_prim_) & (ind_sac_str_event<ind_next_prim_) );
    
    if isempty(ind_sac_str_event_)
        % there is no sac in current trial. prim, corr, prob are invalid
        validity_prim(counter_trial, 1) = validity_prim(counter_trial, 1) & false;
        ind_str_prim(counter_trial, 1)  = ind_current_prim_;
        ind_end_prim(counter_trial, 1)  = ind_current_prim_;
        validity_corr(counter_trial, 1) = validity_corr(counter_trial, 1) & false;
        ind_str_corr(counter_trial, 1)  = ind_current_corr_;
        ind_end_corr(counter_trial, 1)  = ind_current_corr_;
        validity_prob(counter_trial, 1) = validity_prob(counter_trial, 1) & false;
        ind_str_prob(counter_trial, 1)  = ind_current_prob_;
        ind_end_prob(counter_trial, 1)  = ind_current_prob_;
    else
        % there is at least 1 sac. The first sac is prim
        validity_prim(counter_trial, 1) = validity_prim(counter_trial, 1) & true;
        ind_str_prim(counter_trial, 1)  = ind_sac_str_event_(1);
        ind_end_prim(counter_trial, 1)  = ind_sac_end_event_(1);
        if length(ind_sac_str_event_) == 1
            % there exsits 1 and only 1 sac. corr, prob are invalid
            validity_corr(counter_trial, 1) = validity_corr(counter_trial, 1) & false;
            ind_str_corr(counter_trial, 1)  = ind_current_corr_;
            ind_end_corr(counter_trial, 1)  = ind_current_corr_;
            validity_prob(counter_trial, 1) = validity_prob(counter_trial, 1) & false;
            ind_str_prob(counter_trial, 1)  = ind_current_prob_;
            ind_end_prob(counter_trial, 1)  = ind_current_prob_;
        else
            % there are more than 1 sac. The 1st is already chosen to be
            % prim. let's check the conditions for 2nd sac and if it
            % matched the conditions then the 2nd sac is corr
            if ind_sac_str_event_(2) > ind_current_corr_
                % the start of potential corr sac happens after corr target
                % jump
                if validity_prob(counter_trial, 1)
                    % there exists a prob trial. let's check the
                    % possibility that the potential corr sac is actual
                    % corr sac or a prob sac
                    if ind_sac_str_event_(2) < ind_current_prob_
                        validity_corr(counter_trial, 1) = validity_corr(counter_trial, 1) & true;
                        ind_str_corr(counter_trial, 1)  = ind_sac_str_event_(2);
                        ind_end_corr(counter_trial, 1)  = ind_sac_end_event_(2);
                    else
                        % the potential corr sac has happened after prob
                        % target jump. it is a prob sac and not a corr sac
                        validity_corr(counter_trial, 1) = validity_corr(counter_trial, 1) & false;
                        ind_str_corr(counter_trial, 1)  = ind_current_corr_;
                        ind_end_corr(counter_trial, 1)  = ind_current_corr_;
                        validity_prob(counter_trial, 1) = validity_prob(counter_trial, 1) & true;
                        ind_str_prob(counter_trial, 1)  = ind_sac_str_event_(2);
                        ind_end_prob(counter_trial, 1)  = ind_sac_end_event_(2);
                        flag_overlap = true;
                    end
                else
                    validity_corr(counter_trial, 1) = validity_corr(counter_trial, 1) & true;
                    ind_str_corr(counter_trial, 1)  = ind_sac_str_event_(2);
                    ind_end_corr(counter_trial, 1)  = ind_sac_end_event_(2);
                end
            else
                validity_corr(counter_trial, 1) = validity_corr(counter_trial, 1) & false;
                ind_str_corr(counter_trial, 1)  = ind_current_corr_;
                ind_end_corr(counter_trial, 1)  = ind_current_corr_;
            end
            
            if length(ind_sac_str_event_) == 2
                % there exsits 2 and only 2 sac. prob is invalid
                if(~flag_overlap)
                    validity_prob(counter_trial, 1) = validity_prob(counter_trial, 1) & false;
                    ind_str_prob(counter_trial, 1)  = ind_current_prob_;
                    ind_end_prob(counter_trial, 1)  = ind_current_prob_;
                end
            else
                % there are more than 2 sac. The 1st is already chosen to be
                % prim and 2nd to be corr. let's check the conditions for 3rd sac and if it
                % matched the conditions then the 3rd sac is prob
                if ind_sac_str_event_(3) > ind_current_prob_
                    validity_prob(counter_trial, 1) = validity_prob(counter_trial, 1) & true;
                    ind_str_prob(counter_trial, 1)  = ind_sac_str_event_(3);
                    ind_end_prob(counter_trial, 1)  = ind_sac_end_event_(3);
                else
                    validity_prob(counter_trial, 1) = validity_prob(counter_trial, 1) & false;
                    ind_str_prob(counter_trial, 1)  = ind_current_prob_;
                    ind_end_prob(counter_trial, 1)  = ind_current_prob_;
                end
            end
        end
    end
   
end

ind_tgt_jump_prim      = DATA.PRIM_SAC.ind_tgt_jump;
ind_tgt_jump_str_prim  = ind_tgt_jump_prim;
ind_tgt_jump_end_prim  = ind_tgt_jump_prim + 1;
if contains(SMR_FILE.file_name,'SquirrelUnit')
    ind_tgt_jump_str_prim   = ind_tgt_jump_prim - 10;
    ind_tgt_jump_end_prim   = ind_tgt_jump_prim + 10;
end
DATA.PRIM_SAC = ESN_Build_SAC_DATA(DATA, validity_prim, ind_str_prim, ind_end_prim, ind_tgt_jump_prim, ind_tgt_jump_str_prim, ind_tgt_jump_end_prim);
DATA.PRIM_SAC.reaction       = DATA.PRIM_SAC.ind_str - DATA.PRIM_SAC.ind_tgt_jump;

ind_tgt_jump_corr      = DATA.CORR_SAC.ind_tgt_jump;
ind_tgt_jump_str_corr  = ind_tgt_jump_corr;
ind_tgt_jump_end_corr  = ind_tgt_jump_corr + 1;
if contains(SMR_FILE.file_name,'SquirrelUnit')
    ind_tgt_jump_str_corr   = ind_tgt_jump_corr - 10;
    ind_tgt_jump_end_corr   = ind_tgt_jump_corr + 10;
end
DATA.CORR_SAC = ESN_Build_SAC_DATA(DATA, validity_corr, ind_str_corr, ind_end_corr, ind_tgt_jump_corr, ind_tgt_jump_str_corr, ind_tgt_jump_end_corr);
DATA.CORR_SAC.reaction       = DATA.CORR_SAC.ind_str - DATA.PRIM_SAC.ind_end; % override the reaction variable

ind_tgt_jump_prob      = DATA.PROB_SAC.ind_tgt_jump;
ind_tgt_jump_str_prob  = ind_tgt_jump_prob;
ind_tgt_jump_end_prob  = ind_tgt_jump_prob + 1;
if contains(SMR_FILE.file_name,'SquirrelUnit')
    ind_tgt_jump_str_prob   = ind_tgt_jump_prob - 10;
    ind_tgt_jump_end_prob   = ind_tgt_jump_prob + 10;
end
DATA.PROB_SAC = ESN_Build_SAC_DATA(DATA, validity_prob, ind_str_prob, ind_end_prob, ind_tgt_jump_prob, ind_tgt_jump_str_prob, ind_tgt_jump_end_prob);
DATA.PROB_SAC.reaction       = DATA.PROB_SAC.ind_str - DATA.PROB_SAC.ind_tgt_jump;
fprintf(' --> Completed. \n')

%% Build Raster Plot Data
clearvars -except CH_DATA_ALL SMR_FILE DATA;
fprintf(['Building Raster Plot Data ', SMR_FILE.file_name, ' ... ']);

% Build inds_primSac_str
eye_time     = DATA.TRIAL.eye_time;
length_time = length(eye_time);
inds_span_primSac_str = (-50+1) : 1 : (250);
ind_primSac_str  = DATA.PRIM_SAC.ind_str;
inds_primSac_str = repmat( ind_primSac_str(:), 1, length(inds_span_primSac_str)) + repmat(inds_span_primSac_str(:)', length(ind_primSac_str), 1);
inds_primSac_str( inds_primSac_str < 1 ) = 1;
inds_primSac_str( inds_primSac_str > length_time ) = length_time;
DATA.RASTER.ind_primSac_str       = ind_primSac_str;
DATA.RASTER.inds_primSac_str      = inds_primSac_str;
DATA.RASTER.inds_span_primSac_str = inds_span_primSac_str;

% Build inds_primSac_end
inds_span_primSac_end = (-50+1) : 1 : (250);
ind_primSac_end  = DATA.PRIM_SAC.ind_end;
inds_primSac_end = repmat( ind_primSac_end(:), 1, length(inds_span_primSac_end)) + repmat(inds_span_primSac_end(:)', length(ind_primSac_end), 1);
inds_primSac_end( inds_primSac_end < 1 ) = 1;
inds_primSac_end( inds_primSac_end > length_time ) = length_time;
DATA.RASTER.ind_primSac_end       = ind_primSac_end;
DATA.RASTER.inds_primSac_end      = inds_primSac_end;
DATA.RASTER.inds_span_primSac_end = inds_span_primSac_end;

% Build inds_corrSac_str
eye_time     = DATA.TRIAL.eye_time;
length_time = length(eye_time);
inds_span_corrSac_str = (-50+1) : 1 : (250);
ind_corrSac_str  = DATA.CORR_SAC.ind_str;
inds_corrSac_str = repmat( ind_corrSac_str(:), 1, length(inds_span_corrSac_str)) + repmat(inds_span_corrSac_str(:)', length(ind_corrSac_str), 1);
inds_corrSac_str( inds_corrSac_str < 1 ) = 1;
inds_corrSac_str( inds_corrSac_str > length_time ) = length_time;
DATA.RASTER.ind_corrSac_str       = ind_corrSac_str;
DATA.RASTER.inds_corrSac_str      = inds_corrSac_str;
DATA.RASTER.inds_span_corrSac_str = inds_span_corrSac_str;

% Build inds_corrSac_end
inds_span_corrSac_end = (-50+1) : 1 : (250);
ind_corrSac_end  = DATA.CORR_SAC.ind_end;
inds_corrSac_end = repmat( ind_corrSac_end(:), 1, length(inds_span_corrSac_end)) + repmat(inds_span_corrSac_end(:)', length(ind_corrSac_end), 1);
inds_corrSac_end( inds_corrSac_end < 1 ) = 1;
inds_corrSac_end( inds_corrSac_end > length_time ) = length_time;
DATA.RASTER.ind_corrSac_end       = ind_corrSac_end;
DATA.RASTER.inds_corrSac_end      = inds_corrSac_end;
DATA.RASTER.inds_span_corrSac_end = inds_span_corrSac_end;

fprintf(' --> Completed. \n')

%% plot-1 learning 
clearvars -except CH_DATA_ALL SMR_FILE DATA;
validity      = DATA.PRIM_SAC.validity;
gain_u_bool   = DATA.PRIM_SAC.gain_u_bool;
gain_d_bool   = DATA.PRIM_SAC.gain_d_bool;
trial_num     = DATA.PRIM_SAC.trial_num;
eye_amp_rot_x = DATA.PRIM_SAC.eye_amp_rot_x;
eye_amp_rot_x(gain_d_bool) = -1 * eye_amp_rot_x(gain_d_bool);
reaction = DATA.PRIM_SAC.reaction;
% find invalid trials for prim sac
validity = validity & (reaction  < 500); % reaction time should be less than specified value
validity = validity & (eye_amp_rot_x < (mean(DATA.PRIM_SAC.tgt_amp_m)+10)); % saccade amplitude should be less than specified value
validity = validity & (eye_amp_rot_x > (mean(DATA.PRIM_SAC.tgt_amp_m)-10)); % saccade amplitude should be greater than specified value
% figure
fig_num_ = 1;
fig_handle_(fig_num_) = figure(fig_num_);
clf(fig_handle_(fig_num_))
hold on
plot(trial_num(validity&gain_u_bool) , eye_amp_rot_x(validity&gain_u_bool), '.b')
plot(trial_num(validity&gain_d_bool) , eye_amp_rot_x(validity&gain_d_bool), '.r')
xlabel('Trial')
ylabel('Prim sac projected amplitude (deg)')
ESN_Beautify_Plot

%% plot-2 raster PRIM
clearvars -except CH_DATA_ALL SMR_FILE DATA;
% figure
fig_num_ = 2;
fig_handle_(fig_num_) = figure(fig_num_);
clf(fig_handle_(fig_num_))
hold on
validity      = DATA.PRIM_SAC.validity;
gain_u_bool   = DATA.PRIM_SAC.gain_u_bool;
gain_d_bool   = DATA.PRIM_SAC.gain_d_bool;
trial_num     = DATA.PRIM_SAC.trial_num;
eye_amp_rot_x = DATA.PRIM_SAC.eye_amp_rot_x;
eye_amp_rot_x(gain_d_bool) = -1 * eye_amp_rot_x(gain_d_bool);
reaction = DATA.PRIM_SAC.reaction;
% find invalid trials for prim sac
validity = validity & (reaction  < 500); % reaction time should be less than specified value
validity = validity & (eye_amp_rot_x < (mean(DATA.PRIM_SAC.tgt_amp_m)+10)); % saccade amplitude should be less than specified value
validity = validity & (eye_amp_rot_x > (mean(DATA.PRIM_SAC.tgt_amp_m)-10)); % saccade amplitude should be greater than specified value

% data primSac_str raster
inds_trials      = validity;
inds_span        = DATA.RASTER.inds_span_primSac_str;
trial_num_raster = trial_num(inds_trials);
inds_event       = DATA.RASTER.inds_primSac_str;
train_data_event = DATA.TRIAL.spike_event(inds_event(inds_trials,:));
[x_axis_raster, y_axis_raster] = ESN_raster_plot_axes(train_data_event, inds_span, trial_num_raster, 0.5);
train_data_firing = DATA.NEURAL.spike_firing(inds_event(inds_trials,:));
x_axis_surf  = inds_span;
y_axis_surf  = trial_num(inds_trials);
train_data_velocity = DATA.TRIAL_ROT.eye_vm_filt(inds_event(inds_trials,:));
% plot primSac_str raster
subplot(2,3,1)
plot(x_axis_raster(:), y_axis_raster(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
xlabel('Time from prim sac start (ms)')
ylabel('Trial')
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (trial_num_raster(end)+3)])
title('Raster')
subplot(2,3,2)
surf_h_ = surf(x_axis_surf(:)', y_axis_surf(:), train_data_firing);
surf_h_.EdgeColor = 'none';
view(2)
xlabel('Time from prim sac start (ms)')
ylabel('Trial')
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (y_axis_surf(end)+3)])
title('Firing Rate')
subplot(2,3,3)
surf_h_ = surf(x_axis_surf(:)', y_axis_surf(:), train_data_velocity);
surf_h_.EdgeColor = 'none';
view(2)
xlabel('Time from prim sac start (ms)')
ylabel('Trial')
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (y_axis_surf(end)+3)])
title('Eye Velocity')

% data primSac_end raster
inds_trials      = validity;
inds_span        = DATA.RASTER.inds_span_primSac_end;
trial_num_raster = trial_num(inds_trials);
inds_event       = DATA.RASTER.inds_primSac_end;
train_data_event = DATA.TRIAL.spike_event(inds_event(inds_trials,:));
[x_axis_raster, y_axis_raster] = ESN_raster_plot_axes(train_data_event, inds_span, trial_num_raster, 0.5);
train_data_firing = DATA.NEURAL.spike_firing(inds_event(inds_trials,:));
x_axis_surf  = inds_span;
y_axis_surf  = trial_num(inds_trials);
train_data_velocity = DATA.TRIAL_ROT.eye_vm_filt(inds_event(inds_trials,:));
% plot primSac_end raster
subplot(2,3,4)
plot(x_axis_raster(:), y_axis_raster(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
xlabel('Time from prim sac end (ms)')
ylabel('Trial')
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (trial_num_raster(end)+3)])
subplot(2,3,5)
surf_h_ = surf(x_axis_surf(:)', y_axis_surf(:), train_data_firing);
surf_h_.EdgeColor = 'none';
view(2)
xlabel('Time from prim sac end (ms)')
ylabel('Trial')
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (y_axis_surf(end)+3)])
subplot(2,3,6)
surf_h_ = surf(x_axis_surf(:)', y_axis_surf(:), train_data_velocity);
surf_h_.EdgeColor = 'none';
view(2)
xlabel('Time from prim sac end (ms)')
ylabel('Trial')
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (y_axis_surf(end)+3)])

hFig = fig_handle_(fig_num_);
figure_size  = [8.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'Color', [1 1 1]);
set(hFig, 'Units', 'inches');
set(hFig, 'PaperUnits', 'inches');
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

%% plot-3 raster CORR
clearvars -except CH_DATA_ALL SMR_FILE DATA;
% figure
fig_num_ = 3;
fig_handle_(fig_num_) = figure(fig_num_);
clf(fig_handle_(fig_num_))
hold on
validity      = DATA.CORR_SAC.validity;
gain_u_bool   = DATA.CORR_SAC.gain_u_bool;
gain_d_bool   = DATA.CORR_SAC.gain_d_bool;
trial_num     = DATA.CORR_SAC.trial_num;
eye_amp_rot_x = DATA.CORR_SAC.eye_amp_rot_x;
eye_amp_rot_x(gain_d_bool) = -1 * eye_amp_rot_x(gain_d_bool);
reaction = DATA.CORR_SAC.reaction;
% find invalid trials for corr sac
validity = validity & (reaction  < 500); % reaction time should be less than specified value
validity = validity & (eye_amp_rot_x < (mean(DATA.CORR_SAC.tgt_amp_m)+10)); % saccade amplitude should be less than specified value
validity = validity & (eye_amp_rot_x > (mean(DATA.CORR_SAC.tgt_amp_m)-10)); % saccade amplitude should be greater than specified value

% data corrSac_str raster
inds_trials      = validity;
inds_span        = DATA.RASTER.inds_span_corrSac_str;
trial_num_raster = trial_num(inds_trials);
inds_event       = DATA.RASTER.inds_corrSac_str;
train_data_event = DATA.TRIAL.spike_event(inds_event(inds_trials,:));
[x_axis_raster, y_axis_raster] = ESN_raster_plot_axes(train_data_event, inds_span, trial_num_raster, 0.5);
train_data_firing = DATA.NEURAL.spike_firing(inds_event(inds_trials,:));
x_axis_surf  = inds_span;
y_axis_surf  = trial_num(inds_trials);
train_data_velocity = DATA.TRIAL_ROT.eye_vm_filt(inds_event(inds_trials,:));
% plot corrSac_str raster
subplot(2,3,1)
plot(x_axis_raster(:), y_axis_raster(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
xlabel('Time from corr sac start (ms)')
ylabel('Trial')
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (trial_num_raster(end)+3)])
title('Raster')
subplot(2,3,2)
surf_h_ = surf(x_axis_surf(:)', y_axis_surf(:), train_data_firing);
surf_h_.EdgeColor = 'none';
view(2)
xlabel('Time from corr sac start (ms)')
ylabel('Trial')
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (y_axis_surf(end)+3)])
title('Firing Rate')
subplot(2,3,3)
surf_h_ = surf(x_axis_surf(:)', y_axis_surf(:), train_data_velocity);
surf_h_.EdgeColor = 'none';
view(2)
xlabel('Time from corr sac start (ms)')
ylabel('Trial')
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (y_axis_surf(end)+3)])
title('Eye Velocity')

% data corrSac_end raster
inds_trials      = validity;
inds_span        = DATA.RASTER.inds_span_corrSac_end;
trial_num_raster = trial_num(inds_trials);
inds_event       = DATA.RASTER.inds_corrSac_end;
train_data_event = DATA.TRIAL.spike_event(inds_event(inds_trials,:));
[x_axis_raster, y_axis_raster] = ESN_raster_plot_axes(train_data_event, inds_span, trial_num_raster, 0.5);
train_data_firing = DATA.NEURAL.spike_firing(inds_event(inds_trials,:));
x_axis_surf  = inds_span;
y_axis_surf  = trial_num(inds_trials);
train_data_velocity = DATA.TRIAL_ROT.eye_vm_filt(inds_event(inds_trials,:));
% plot corrSac_end raster
subplot(2,3,4)
plot(x_axis_raster(:), y_axis_raster(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
xlabel('Time from corr sac end (ms)')
ylabel('Trial')
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (trial_num_raster(end)+3)])
subplot(2,3,5)
surf_h_ = surf(x_axis_surf(:)', y_axis_surf(:), train_data_firing);
surf_h_.EdgeColor = 'none';
view(2)
xlabel('Time from corr sac end (ms)')
ylabel('Trial')
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (y_axis_surf(end)+3)])
subplot(2,3,6)
surf_h_ = surf(x_axis_surf(:)', y_axis_surf(:), train_data_velocity);
surf_h_.EdgeColor = 'none';
view(2)
xlabel('Time from corr sac end (ms)')
ylabel('Trial')
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (y_axis_surf(end)+3)])

hFig = fig_handle_(fig_num_);
figure_size  = [8.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'Color', [1 1 1]);
set(hFig, 'Units', 'inches');
set(hFig, 'PaperUnits', 'inches');
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

%% plot-4 raw data 
% figure
fig_num_ = 4;
fig_handle_(fig_num_) = figure(fig_num_);
clf(fig_handle_(fig_num_))
hold on
plot(DATA.TRIAL_ROT.tgt_px, '.-')
plot(DATA.TRIAL_ROT.eye_px_filt, '.-')
% plot(DATA.TRIAL.tgt_px, DATA.TRIAL.tgt_py, '.')

%% function ESN_raster_plot_axes 
function [x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, x_axis_values, y_axis_values, line_half_len)
if nargin < 2
    x_axis_values = 1 : 1 : size(train_data_logic, 2);
    line_half_len = 0.5;
end
if nargin < 3
    line_half_len = 0.5;
end

train_data_row_number = nan(size(train_data_logic));
for counter_row = 1 : size(train_data_logic, 1)
    train_data_row_number(counter_row, train_data_logic(counter_row,:)) = y_axis_values(counter_row);
end
train_data_col_number = repmat(x_axis_values(:)', size(train_data_logic,1), 1);
x_axis = [train_data_col_number(:)'; train_data_col_number(:)'; nan(length(train_data_col_number(:)), 1)'];
y_axis = [(train_data_row_number(:)-line_half_len)'; (train_data_row_number(:)+line_half_len)'; nan(length(train_data_row_number(:)), 1)'];
x_axis = x_axis(:);
y_axis = y_axis(:);
end
%% function ESN_Build_SAC_DATA
function SAC_DATA = ESN_Build_SAC_DATA(DATA, validity, ind_str, ind_end, ind_tgt_jump, ind_tgt_jump_str, ind_tgt_jump_end)
SAC_DATA.validity          = validity;
SAC_DATA.ind_str           = ind_str;
SAC_DATA.ind_end           = ind_end;
SAC_DATA.ind_tgt_jump      = ind_tgt_jump;
SAC_DATA.ind_tgt_jump_str  = ind_tgt_jump_str;
SAC_DATA.ind_tgt_jump_end  = ind_tgt_jump_end;
SAC_DATA.eye_px_str     = DATA.TRIAL.eye_px_filt(    ind_str);
SAC_DATA.eye_px_end     = DATA.TRIAL.eye_px_filt(    ind_end);
SAC_DATA.eye_px_rot_str = DATA.TRIAL_ROT.eye_px_filt(ind_str);
SAC_DATA.eye_px_rot_end = DATA.TRIAL_ROT.eye_px_filt(ind_end);
SAC_DATA.eye_py_str     = DATA.TRIAL.eye_py_filt(    ind_str);
SAC_DATA.eye_py_end     = DATA.TRIAL.eye_py_filt(    ind_end);
SAC_DATA.eye_py_rot_str = DATA.TRIAL_ROT.eye_py_filt(ind_str);
SAC_DATA.eye_py_rot_end = DATA.TRIAL_ROT.eye_py_filt(ind_end);
SAC_DATA.eye_amp_x      = SAC_DATA.eye_px_end - SAC_DATA.eye_px_str;
SAC_DATA.eye_amp_y      = SAC_DATA.eye_py_end - SAC_DATA.eye_py_str;
SAC_DATA.eye_amp_m      = sqrt(SAC_DATA.eye_amp_x.^2 + SAC_DATA.eye_amp_y.^2);
SAC_DATA.eye_amp_rot_x  = SAC_DATA.eye_px_rot_end - SAC_DATA.eye_px_rot_str;
SAC_DATA.eye_amp_rot_y  = SAC_DATA.eye_py_rot_end - SAC_DATA.eye_py_rot_str;
SAC_DATA.eye_amp_rot_m  = sqrt(SAC_DATA.eye_amp_rot_x.^2 + SAC_DATA.eye_amp_rot_y.^2);
SAC_DATA.tgt_px_str     = DATA.TRIAL.tgt_px(    ind_tgt_jump_str);
SAC_DATA.tgt_px_end     = DATA.TRIAL.tgt_px(    ind_tgt_jump_end);
SAC_DATA.tgt_px_rot_str = DATA.TRIAL_ROT.tgt_px(ind_tgt_jump_str);
SAC_DATA.tgt_px_rot_end = DATA.TRIAL_ROT.tgt_px(ind_tgt_jump_end);
SAC_DATA.tgt_py_str     = DATA.TRIAL.tgt_py(    ind_tgt_jump_str);
SAC_DATA.tgt_py_end     = DATA.TRIAL.tgt_py(    ind_tgt_jump_end);
SAC_DATA.tgt_py_rot_str = DATA.TRIAL_ROT.tgt_py(ind_tgt_jump_str);
SAC_DATA.tgt_py_rot_end = DATA.TRIAL_ROT.tgt_py(ind_tgt_jump_end);
SAC_DATA.tgt_amp_x      = SAC_DATA.tgt_px_end - SAC_DATA.tgt_px_str;
SAC_DATA.tgt_amp_y      = SAC_DATA.tgt_py_end - SAC_DATA.tgt_py_str;
SAC_DATA.tgt_amp_m      = sqrt(SAC_DATA.tgt_amp_x.^2 + SAC_DATA.tgt_amp_y.^2);
SAC_DATA.tgt_amp_rot_x  = SAC_DATA.tgt_px_rot_end - SAC_DATA.tgt_px_rot_str;
SAC_DATA.tgt_amp_rot_y  = SAC_DATA.tgt_py_rot_end - SAC_DATA.tgt_py_rot_str;
SAC_DATA.tgt_amp_rot_m  = sqrt(SAC_DATA.tgt_amp_rot_x.^2 + SAC_DATA.tgt_amp_rot_y.^2);
SAC_DATA.tgt_ang        = atan2d(SAC_DATA.tgt_amp_y, SAC_DATA.tgt_amp_x);
SAC_DATA.reaction       = ind_str - ind_tgt_jump;
SAC_DATA.trial_num      = ( 1 : 1 : length(ind_str) )';
SAC_DATA.eye_px_rot_end_centered = SAC_DATA.eye_px_rot_end - SAC_DATA.tgt_px_rot_str;
SAC_DATA.eye_py_rot_end_centered = SAC_DATA.eye_py_rot_end - SAC_DATA.tgt_py_rot_str;
SAC_DATA.gain_u_bool    = SAC_DATA.tgt_amp_rot_x >= 0;
SAC_DATA.gain_d_bool    = SAC_DATA.tgt_amp_rot_x <  0;
end