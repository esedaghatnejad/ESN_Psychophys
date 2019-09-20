function ESN_plot_trial_by_trial
%% load EPHYS EVENT DATA
[file_name,file_path] = uigetfile([pwd filesep '*_aligned.mat'], 'Select EVENT DATA file');
fprintf(['Loading ', file_name, ' ... ']);
EPHYS.CH_EVE = load([file_path filesep file_name]);
fprintf(' --> Completed. \n')

%% load BEHAVE DATA
[file_name,file_path] = uigetfile([file_path filesep '*_ANALYZED.mat'], 'Select _ANALYZED file');
fprintf(['Loading ', file_name, ' ... ']);
BEHAVE = load([file_path filesep file_name]);
fprintf(' --> Completed. \n')

%% load EPHYS sorted DATA
% The fiel should have ch_data, CS_data, SS_data format
[file_name,file_path] = uigetfile([file_path filesep '*_sorted.mat'], 'Select _sorted file');
fprintf(['Loading ', file_name, ' ... ']);
EPHYS.CH_sorted = load([file_path filesep file_name]);
EPHYS.CH_sorted_file_name = file_name;
EPHYS.CH_sorted_file_path = file_path;
fprintf(' --> Completed. \n')

%% Extract Trial Data
clearvars -except EPHYS BEHAVE
num_trials = length(EPHYS.CH_EVE.BEHAVE_trial_nums_aligned);
inds_fixation_begin_EPHYS  = EPHYS.CH_EVE.EPHYS_ind_fixation_begin;
inds_fixation_begin_BEHAVE = EPHYS.CH_EVE.BEHAVE_ind_fixation_begin;
% inds_fixation_end_EPHYS    = EPHYS.CH_EVE.EPHYS_ind_fixation_end;
% inds_fixation_end_BEHAVE   = EPHYS.CH_EVE.BEHAVE_ind_fixation_end;

% for counter_trial = 1 : 1 : num_trials
    counter_trial = 33;
    time_fixation_begin_BEHAVE_ = EPHYS.CH_EVE.BEHAVE_time_aligned(inds_fixation_begin_BEHAVE(counter_trial) );
    time_cueonset_begin_BEHAVE_ = BEHAVE.TRIALS_DATA.time_state_sac_detect_on{1, EPHYS.CH_EVE.BEHAVE_trial_nums_aligned(counter_trial)}(end);
    time_trial_end_BEHAVE_      = BEHAVE.TRIALS_DATA.time_state_next_trial(1, EPHYS.CH_EVE.BEHAVE_trial_nums_aligned(counter_trial));
    ind_cueonset_begin_BEHAVE_  = find(BEHAVE.TRIALS_DATA.time_1K{1, EPHYS.CH_EVE.BEHAVE_trial_nums_aligned(counter_trial)} >= time_cueonset_begin_BEHAVE_, 1, 'first');%find(EPHYS.CH_EVE.BEHAVE_time >= time_cueonset_begin_BEHAVE_, 1, 'first');
    ind_trial_end_BEHAVE_       = find(BEHAVE.TRIALS_DATA.time_1K{1, EPHYS.CH_EVE.BEHAVE_trial_nums_aligned(counter_trial)} >= time_trial_end_BEHAVE_,      1, 'last');%find(EPHYS.CH_EVE.BEHAVE_time <= time_trial_end_BEHAVE_,      1, 'last');
    time_fixation_begin_EPHYS_  = EPHYS.CH_EVE.EPHYS_time_aligned(inds_fixation_begin_EPHYS(counter_trial) );
    time_cueonset_begin_EPHYS_  = time_fixation_begin_EPHYS_ - (time_fixation_begin_BEHAVE_ - time_cueonset_begin_BEHAVE_);
    time_trial_end_EPHYS_       = time_fixation_begin_EPHYS_ + (time_trial_end_BEHAVE_ - time_fixation_begin_BEHAVE_);
    ind_cueonset_begin_EPHYS_   = find(EPHYS.CH_EVE.EPHYS_time  >= time_cueonset_begin_EPHYS_,  1, 'first');
    ind_trial_end_EPHYS_        = find(EPHYS.CH_EVE.EPHYS_time  <= time_trial_end_EPHYS_,       1, 'last');
    
    time_shift_ = 0.000;
    inds_trial_EPHYS_  = (round(ind_cueonset_begin_EPHYS_-(time_shift_*30e3))  : 1 : ind_trial_end_EPHYS_)';
    inds_trial_BEHAVE_ = (round(ind_cueonset_begin_BEHAVE_-(time_shift_*1e3))  : 1 : ind_trial_end_BEHAVE_)';
    
    trial_time_EPHYS_ = EPHYS.CH_EVE.EPHYS_time(inds_trial_EPHYS_)';
    trial_time_EPHYS_ = trial_time_EPHYS_ - trial_time_EPHYS_(1)-time_shift_;
    
    trial_time_BEHAVE_ = EPHYS.CH_EVE.BEHAVE_time(inds_trial_BEHAVE_)';
    trial_time_BEHAVE_ = trial_time_BEHAVE_ - trial_time_BEHAVE_(1)-time_shift_;
    
    trial_eye_r_px = BEHAVE.TRIALS_DATA.eye_r_px_filt{1, EPHYS.CH_EVE.BEHAVE_trial_nums_aligned(counter_trial)}(inds_trial_BEHAVE_);
    trial_eye_r_py = BEHAVE.TRIALS_DATA.eye_r_py_filt{1, EPHYS.CH_EVE.BEHAVE_trial_nums_aligned(counter_trial)}(inds_trial_BEHAVE_);
    trial_tgt_px   = BEHAVE.TRIALS_DATA.tgt_px{       1, EPHYS.CH_EVE.BEHAVE_trial_nums_aligned(counter_trial)}(inds_trial_BEHAVE_);
    trial_tgt_py   = BEHAVE.TRIALS_DATA.tgt_py{       1, EPHYS.CH_EVE.BEHAVE_trial_nums_aligned(counter_trial)}(inds_trial_BEHAVE_);
    trial_state_fixation_ephys  = EPHYS.CH_EVE.EPHYS_state_fixation(inds_trial_EPHYS_);

clf(figure(1))
num_rows = 1;
num_cols = 1;
subplot_index = reshape(1:num_cols*num_rows, num_cols, num_rows).';
clearvars axis_handle;

counter_fig = 1;
axis_handle(counter_fig) = subplot(num_rows+2,num_cols,subplot_index(counter_fig));
plot(trial_time_EPHYS_, [ EPHYS.CH_sorted.ch_data.ch_data_hipass(inds_trial_EPHYS_) EPHYS.CH_sorted.ch_data.ch_data_lopass(inds_trial_EPHYS_) ]);
title( EPHYS.CH_sorted_file_name , 'interpret', 'none');
grid on

counter_fig = counter_fig + 1;
axis_handle(counter_fig) = subplot(num_rows+2,num_cols,counter_fig);
plot(trial_time_BEHAVE_, [ trial_eye_r_px trial_tgt_px ]);
ylabel( 'eye_r_px' , 'interpret', 'none');
grid on

counter_fig = counter_fig + 1;
axis_handle(counter_fig) = subplot(num_rows+2,num_cols,counter_fig);
hold on
plot(trial_time_BEHAVE_, [ trial_eye_r_py trial_tgt_py ]);
plot(trial_time_EPHYS_, trial_state_fixation_ephys);
ylabel( 'eye_r_py' , 'interpret', 'none');
xlabel(['trial # ' num2str(counter_trial)], 'interpret', 'none')
grid on

linkaxes(axis_handle,'x')
ESN_Beautify_Plot

% waitforbuttonpress
% end
