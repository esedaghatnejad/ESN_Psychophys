%% function ESN_population_coding_sac_sorter
function ESN_population_coding_sac_sorter
%% Global variables
global data_type_eye_list data_type_BEHAVE_list data_type_neuro_list data_type_EPHYS_list event_type_list ...
    length_trace inds_span ...
    ang_step ang_edges ang_values amp_edges vel_edges ...
    range_cell_with_4dir_behave
data_type_eye_list    = {'eye_vx', 'eye_vy'};
data_type_BEHAVE_list = {'eye_r_vx_filt', 'eye_r_vy_filt'};
data_type_neuro_list  = {'neuro_SS', 'neuro_CS'};
data_type_EPHYS_list  = {'EPHYS_SS_train_1K', 'EPHYS_CS_train_1K'};
event_type_list       = {'visual', 'onset', 'vmax', 'offset', 'auditory'};
length_trace = 500;
inds_span    = ((-(length_trace/2)+1) : 1 : (length_trace/2))';
ang_step     = 45;
ang_edges    = (0 - (ang_step/2)) : ang_step : (360 + (ang_step/2));
ang_values   = (0) : ang_step : (360 - ang_step);
amp_edges    = [-0.5 2 4 6 8 10 50];
vel_edges    = [0 200 300 400 500 600 10000];
range_cell_with_4dir_behave = [1 51]; % This is to correct the data for the first round of recordings from Mirza, pCell_list_Mirza_pre201906
% plese set the range_cell_with_4dir_behave to [-1 -1] if you do not have 4dir sessions

%% Steps to prepare the combined eye & neuro data
tic
% (1) re_run_ESN_sac_sorter(); % this func will re-analyze the _ANALYZED files and add the SACS_ALL_DATA with tags to them
% (2) re_run_add_ephys_sac_sorter(); % this func will fuse the eye & neuro data and save the results with _sac name
% (3) combine_sac_files(); % combine the _sac files and form one file per cell. save the data in the ALL_PCELL_### folder.
% (4) CS_on_analysis(); % load cell_data files (_combine_) and add the CS_on_analysis  to them.
% (5) 
build_population_data(); % load cell_data files (_combine_) and form population data.
toc
end

%% function RE-RUN ESN_Sac_Sorter
function re_run_ESN_sac_sorter()
clc; close all;
pCell_list = ESN_build_pCell_list();
path_data_monkey_sorted = uigetdir;

if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);
%% Loop over pCells
for counter_pCell = 1 : 1 : num_pCells
    fprintf(['### ' 'Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ###' '\n']);
    num_recording = sum(pCell_list_isstr(counter_pCell, :));
    for counter_recording = 1 : 1 : num_recording
        %% build plot_data address
        file_name_cell = pCell_list{counter_pCell, counter_recording}; % '190423_142023_01_sorted_ESN_plot_data';
        year_ = file_name_cell(1:2);
        month_ = file_name_cell(3:4);
        day_ = file_name_cell(5:6);
        hour_ = file_name_cell(8:9);
        minute_ = file_name_cell(10:11);
        second_ = file_name_cell(12:13);
        subFolder_month = ['20' year_ '-' month_ filesep];
        subFolder_day = ['20' year_ '-' month_ '-' day_ filesep];
        subFolder_recording = ['20' year_ '-' month_ '-' day_ '_' hour_ '-' minute_ '-' second_ filesep];
        subFolder_figs = ['analyzed_data' filesep];
        %% load BEHAVE data
        file_path = [path_data_monkey_sorted subFolder_month subFolder_day subFolder_recording subFolder_figs];
        file_name = [file_name_cell(1:13) '_ANALYZED.mat'];
        if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
        BEHAVE = load([file_path file_name], 'EXPERIMENT_PARAMS', 'TRIALS_DATA');
        %% Build SACS_ALL_DATA using ESN_Sac_Sorter
        flag_session_figure = true;
        [SACS_ALL_DATA, TRIALS_DATA, EXPERIMENT_PARAMS] = ESN_Sac_Sorter(BEHAVE.TRIALS_DATA, BEHAVE.EXPERIMENT_PARAMS, flag_session_figure);

        %% Save _ANALYZED.mat Data to disk
        save([file_path file_name_cell(1:13) '_ANALYZED.mat'], ...
            'EXPERIMENT_PARAMS', 'TRIALS_DATA', 'SACS_ALL_DATA', '-v7.3');

        %% Save _REDUCED.mat Data to disk
        rmfields_list = {'eye_l_vm_filt', 'eye_l_vy_filt', 'eye_l_vx_filt', 'eye_l_py_filt', 'eye_l_px_filt', ...
            'eye_r_vm_filt', 'eye_r_vy_filt', 'eye_r_vx_filt', 'eye_r_py_filt', 'eye_r_px_filt', ...
            'time', 'time_1K', 'target_visible', 'reward', 'tgt_py', 'tgt_px', 'time_tgt', ...
            'eye_l_vm', 'eye_r_vm', 'eye_l_vy', 'eye_l_vx', 'eye_r_vy', 'eye_r_vx', ...
            'eye_l_py', 'eye_l_px', 'eye_r_py', 'eye_r_px', 'time_eyelink', 'inds_invalid', 'inds_trial'};
        TRIALS_DATA = rmfield(TRIALS_DATA,rmfields_list);
        save([file_path file_name_cell(1:13) '_REDUCED.mat'], ...
            'EXPERIMENT_PARAMS', 'TRIALS_DATA', 'SACS_ALL_DATA', '-v7.3');

        %% Save Fig
        hFig_ = gcf;
        file_name_plot_ = file_name_cell(1:13);
        file_path_plot_ = [file_path '..' filesep 'analyzed_figs' filesep];
        saveas(hFig_,[file_path_plot_ file_name_plot_ '_sac_sorter'], 'pdf');
        saveas(hFig_,[file_path_plot_ file_name_plot_ '_sac_sorter'], 'png');
        close(hFig_);
    end
end
fprintf('### ALL DONE. ###\n')
end

%% function RE-RUN add_ephys_sac_sorter
function re_run_add_ephys_sac_sorter()
clc; close all;
pCell_list = ESN_build_pCell_list();
path_data_monkey_sorted = uigetdir;

if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);
%% Loop over pCells
for counter_pCell = 1 : 1 : num_pCells
    fprintf(['### ' 'Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ###' '\n']);
    num_recording = sum(pCell_list_isstr(counter_pCell, :));
    for counter_recording = 1 : 1 : num_recording
        %% build plot_data address
        file_name_cell = pCell_list{counter_pCell, counter_recording}; % '190423_142023_01_sorted_ESN_plot_data';
        year_ = file_name_cell(1:2);
        month_ = file_name_cell(3:4);
        day_ = file_name_cell(5:6);
        hour_ = file_name_cell(8:9);
        minute_ = file_name_cell(10:11);
        second_ = file_name_cell(12:13);
        subFolder_month = ['20' year_ '-' month_ filesep];
        subFolder_day = ['20' year_ '-' month_ '-' day_ filesep];
        subFolder_recording = ['20' year_ '-' month_ '-' day_ '_' hour_ '-' minute_ '-' second_ filesep];
        subFolder_figs = ['analyzed_data' filesep];
        %% RE-RUN add_ephys_sac_sorter
        file_path = [path_data_monkey_sorted subFolder_month subFolder_day subFolder_recording subFolder_figs];
        if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
        file_name = file_name_cell;
        [SACS_ALL_DATA, Neural_Properties, EXPERIMENT_PARAMS] = add_ephys_sac_sorter(file_path, file_name);
        %% Save _sac Data to disk
        save([file_path file_name_cell(1:end-4) '_sac' file_name_cell(end-3:end) '.mat'], ...
            'SACS_ALL_DATA', 'Neural_Properties', 'EXPERIMENT_PARAMS', '-v7.3');
        
    end
end
fprintf('### ALL DONE. ###\n')
end

%% function add_ephys_sac_sorter
function [SACS_ALL_DATA, Neural_Properties, EXPERIMENT_PARAMS] = add_ephys_sac_sorter(file_path, file_name)
%% load EPHYS sorted DATA
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
file_name = [file_name '.psort'];
fprintf(['Loading ', file_name, ' ... ']);
DATA_PSORT = Psort_read_psort([file_path file_name]);
EPHYS.CH_sorted_file_name = file_name;
EPHYS.CH_sorted_file_path = file_path;
fprintf(' --> Completed. \n')

%% load EPHYS EVENT DATA
file_name = [EPHYS.CH_sorted_file_name(1:13) '_EVE1_aligned.mat'];
EPHYS.CH_EVE = load([file_path file_name]);
if isfield(EPHYS.CH_EVE, 'EPHYS_time_15K')
    EPHYS.CH_EVE.EPHYS_time_30K = EPHYS.CH_EVE.EPHYS_time_15K(:);
else
    EPHYS.CH_EVE.EPHYS_time_30K = EPHYS.CH_EVE.EPHYS_time_30K(:);
end
EPHYS.CH_EVE.EPHYS_time_1K  = reshape(EPHYS.CH_EVE.EPHYS_time_1K ,[], 1);
EPHYS.CH_EVE.BEHAVE_time_1K = reshape(EPHYS.CH_EVE.BEHAVE_time_1K,[], 1);

%% load BEHAVE DATA
file_name = [EPHYS.CH_sorted_file_name(1:13) '_ANALYZED.mat'];
BEHAVE = load([file_path file_name]);
BEHAVE.EXPERIMENT_PARAMS.EPHYS_file_name = EPHYS.CH_sorted_file_name;
BEHAVE.EXPERIMENT_PARAMS.EPHYS_file_path = EPHYS.CH_sorted_file_path;

%% build EPHYS.CH_sorted from DATA_PSORT
ch_data = double(DATA_PSORT.topLevel_data.ch_data);
ch_time = double(DATA_PSORT.topLevel_data.ch_time);
SS_index = find(logical(double(DATA_PSORT.topLevel_data.ss_index)));
CS_index = find(logical(double(DATA_PSORT.topLevel_data.cs_index)));
SS_time = ch_time(SS_index);
CS_time = ch_time(CS_index);

waveform_inds_span = ((-60+1) : 1 : (120));
SS_inds = repmat(waveform_inds_span(:)', length(SS_index), 1) + repmat(SS_index(:), 1, length(waveform_inds_span));
SS_inds(SS_inds < 1) = 1;
SS_inds(SS_inds > length(ch_data)) = length(ch_data);
CS_inds = repmat(waveform_inds_span(:)', length(CS_index), 1) + repmat(CS_index(:), 1, length(waveform_inds_span));
CS_inds(CS_inds < 1) = 1;
CS_inds(CS_inds > length(ch_data)) = length(ch_data);
SS_waveform = ch_data(SS_inds);
CS_waveform = ch_data(CS_inds);

SS_waveform  = reshape(SS_waveform,[], length(waveform_inds_span));
CS_waveform  = reshape(CS_waveform,[], length(waveform_inds_span));

EPHYS.CH_sorted.SS_data.SS_ind = SS_index;
EPHYS.CH_sorted.CS_data.CS_ind = CS_index;
EPHYS.CH_sorted.SS_data.SS_time = SS_time;
EPHYS.CH_sorted.CS_data.CS_time = CS_time;
EPHYS.CH_sorted.SS_data.SS_waveform = SS_waveform;
EPHYS.CH_sorted.CS_data.CS_waveform = CS_waveform;
EPHYS.CH_sorted.duration = ch_time(end) - ch_time(1);

%% build SSxSS AUTO PROBABILITY
clearvars -except EPHYS BEHAVE
SS_time   = EPHYS.CH_sorted.SS_data.SS_time;
CS_time   = EPHYS.CH_sorted.CS_data.CS_time;
Corr_data = ESN_correlogram(SS_time, CS_time);
EPHYS.CH_sorted.Corr_data = Corr_data;

%% SS & CS train_aligned
clearvars -except EPHYS BEHAVE
EPHYS_time_1K     = EPHYS.CH_EVE.EPHYS_time_1K;
length_time_ = length(EPHYS_time_1K);
CS_time = EPHYS.CH_sorted.CS_data.CS_time;
if isempty(CS_time)
    CS_time = EPHYS_time_1K(1);
end
CS_time(end+1) = max([EPHYS_time_1K(end), CS_time(end)])+1;
SS_time = EPHYS.CH_sorted.SS_data.SS_time;
if isempty(SS_time)
    SS_time = EPHYS_time_1K(1);
end
SS_time(end+1) = max([EPHYS_time_1K(end), SS_time(end)])+1;
EPHYS_CS_train_1K = false(size(EPHYS_time_1K));
EPHYS_SS_train_1K = false(size(EPHYS_time_1K));
counter_CS = find(CS_time >= EPHYS_time_1K(1), 1, 'first');
counter_SS = find(SS_time >= EPHYS_time_1K(1), 1, 'first');
for counter_time_point = 1 : length_time_
    time_ponit_ = EPHYS_time_1K(counter_time_point);
    if time_ponit_>=CS_time(counter_CS)
        EPHYS_CS_train_1K(counter_time_point) = true;
        counter_CS = counter_CS + 1;
    end
    if time_ponit_>=SS_time(counter_SS)
        EPHYS_SS_train_1K(counter_time_point) = true;
        counter_SS = counter_SS + 1;
    end
end
EPHYS.CH_EVE.EPHYS_CS_train_1K = EPHYS_CS_train_1K;
EPHYS.CH_EVE.EPHYS_SS_train_1K = EPHYS_SS_train_1K;

%% extract BEHAVE_ind from BEHAVE_EB_xcorr_time_1K
clearvars -except EPHYS BEHAVE
global data_type_eye_list data_type_BEHAVE_list data_type_neuro_list data_type_EPHYS_list event_type_list inds_span
fprintf(['Analyzing ', EPHYS.CH_sorted_file_name, ' ... ']);
REFRENCE_TIME = EPHYS.CH_EVE.BEHAVE_time_1K;
length_time_  = length(REFRENCE_TIME);
num_sacs      = length(BEHAVE.SACS_ALL_DATA.tag);

% Init variables
% use nan for eye data and logical false for neuro data
for counter_event_type = 1 : length(event_type_list)
    event_type_name = event_type_list{counter_event_type};
    for counter_data_type_eye = 1 : length(data_type_eye_list)
        data_type_eye_name = data_type_eye_list{counter_data_type_eye};
        BEHAVE.SACS_ALL_DATA.([data_type_eye_name '_' event_type_name]) = nan(length_trace, num_sacs);
    end
    for counter_data_type_neuro = 1 : length(data_type_neuro_list)
        data_type_neuro_name = data_type_neuro_list{counter_data_type_neuro};
        BEHAVE.SACS_ALL_DATA.([data_type_neuro_name '_' event_type_name]) = false(length_trace, num_sacs);
    end
end

% Build stream dataset. 
% concatenate the data using the cell2mat and then interpolate the potential missing data
BEHAVE.stream = struct;
time_1K_cell2mat = cell2mat(BEHAVE.TRIALS_DATA.('time_1K')(:))';
BEHAVE.stream.('time_1K') = REFRENCE_TIME;
for counter_variable = 1 : 1 : length(data_type_BEHAVE_list)
    variable_name = data_type_BEHAVE_list{counter_variable};
    variable_cell2mat = cell2mat(BEHAVE.TRIALS_DATA.(variable_name)(:))';
    BEHAVE.stream.(variable_name) = interp1(time_1K_cell2mat, variable_cell2mat, REFRENCE_TIME, 'nearest', 'extrap');
end

% Loop over saccades and add the eye traces to the SACS_ALL_DATA
for counter_sac = 1 : num_sacs
    for counter_event_type = 1 : length(event_type_list)
        event_type_name = event_type_list{counter_event_type};
        event_time = BEHAVE.SACS_ALL_DATA.(['time' '_' event_type_name])(1,counter_sac);
        if isnan(event_time)
            % if the event_time is nan, then skip the event.
            continue;
        end
        event_ind = find(REFRENCE_TIME >= event_time, 1, 'first');
        if isempty(event_ind)
            % if the event_ind is empty, then skip the event.
            BEHAVE.SACS_ALL_DATA.validity(1, counter_sac) = false;
            continue;
        end
        event_inds = repmat( event_ind, 1, length(inds_span)) + repmat(inds_span(:)', length(event_ind), 1);
        event_inds( event_inds < 1 ) = 1;
        event_inds( event_inds > length_time_ ) = length_time_;
        for counter_data_type_eye = 1 : length(data_type_eye_list)
            data_type_eye_name = data_type_eye_list{counter_data_type_eye};
            data_type_BEHAVE_name = data_type_BEHAVE_list{counter_data_type_eye};
            BEHAVE.SACS_ALL_DATA.([data_type_eye_name '_' event_type_name])(:,counter_sac) = ...
                reshape(BEHAVE.stream.(data_type_BEHAVE_name)(event_inds), length_trace, 1);
        end
    end
end

% Loop over saccades and add neuro_SS & neuro_CS to the SACS_ALL_DATA
REFRENCE_TIME = EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_time_1K(:); % the initial ind will be drawn from BEHAVE_EB_xcorr_time_1K
length_time_  = length(EPHYS.CH_EVE.EPHYS_time_1K); % the converted ind will be applied to EPHYS_time_1K
for counter_sac = 1 : num_sacs
    for counter_event_type = 1 : length(event_type_list)
        event_type_name = event_type_list{counter_event_type};
        event_time = BEHAVE.SACS_ALL_DATA.(['time' '_' event_type_name])(1,counter_sac);
        if isnan(event_time)
            % if the event_time is nan, then skip the event.
            continue;
        end
        event_ind = find(REFRENCE_TIME >= event_time, 1, 'first');
        if isempty(event_ind)
            % if the event_ind is empty, then skip the event.
            BEHAVE.SACS_ALL_DATA.validity(1, counter_sac) = false;
            continue;
        end
        % set the event_ind_converted based on align_states
        event_ind_converted = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(event_ind);
        % re-write the event_ind_converted to be based on align_photodiode if
        % the event is a cue presentation
        if strcmp(event_type_name,'visual')
            tag_ = BEHAVE.SACS_ALL_DATA.tag(1, counter_sac);
            if (tag_ == 1) || (tag_ == 2) || (tag_ == 3) || (tag_ == 6) || (tag_ == 7)
                % 'prim_success' tag 1 % 'prim_attempt' tag 2 % 'prim_fail' tag 3 % 'back_center_success' tag 6 % 'back_center_prim' tag 7
                event_ind_converted = EPHYS.CH_EVE.align_photodiode.EPHYS_PD_aligned_ind_1K(event_ind);
            end
        end
        if isempty(event_ind_converted)
            % if the event_ind_converted is empty, then skip the event.
            BEHAVE.SACS_ALL_DATA.validity(1, counter_sac) = false;
            continue;
        end
        event_inds_converted = repmat( event_ind_converted, 1, length(inds_span)) + repmat(inds_span(:)', length(event_ind_converted), 1);
        event_inds_converted( event_inds_converted < 1 ) = 1;
        event_inds_converted( event_inds_converted > length_time_ ) = length_time_;
        for counter_data_type_eye = 1 : length(data_type_eye_list)
            data_type_neuro_name = data_type_neuro_list{counter_data_type_eye};
            data_type_EPHYS_name = data_type_EPHYS_list{counter_data_type_eye};
            BEHAVE.SACS_ALL_DATA.([data_type_neuro_name '_' event_type_name])(:,counter_sac) = ...
                logical(reshape(EPHYS.CH_EVE.(data_type_EPHYS_name)(event_inds_converted), length_trace, 1));
        end
    end
end

% Add neuro_CS_count_visual & neuro_CS_count_auditory to BEHAVE.SACS_ALL_DATA
% count the number of the the CS after the visual/auditory event.
% from the event till 200ms or sac onset, whichever happen first
BEHAVE.SACS_ALL_DATA.neuro_CS_count_visual   = zeros(size(BEHAVE.SACS_ALL_DATA.tag));
BEHAVE.SACS_ALL_DATA.neuro_CS_count_auditory = zeros(size(BEHAVE.SACS_ALL_DATA.tag));
for counter_sac = 1 : num_sacs
    ind_start     = (length_trace/2);
    time_onset    = BEHAVE.SACS_ALL_DATA.time_onset(   1, counter_sac);
    time_visual   = BEHAVE.SACS_ALL_DATA.time_visual(  1, counter_sac);
    time_auditory = BEHAVE.SACS_ALL_DATA.time_auditory(1, counter_sac);
    reaction_visual   = round((time_onset - time_visual  )*1000.0);
    reaction_auditory = round((time_onset - time_auditory)*1000.0);
    
    if (reaction_visual > 0) && (reaction_visual < 200)
        ind_end_visual = ind_start + reaction_visual;
    elseif (reaction_visual >= 200)
        ind_end_visual = ind_start + 200;
    elseif (reaction_visual <= 0)
        ind_end_visual = ind_start;
    else
        % this condition covers the nan values
        ind_end_visual = ind_start;
    end
    
    if (reaction_auditory > 0) && (reaction_auditory < 200)
        ind_end_auditory = ind_start + reaction_auditory;
    elseif (reaction_auditory >= 200)
        ind_end_auditory = ind_start + 200;
    elseif (reaction_auditory <= 0)
        ind_end_auditory = ind_start;
    else
        % this condition covers the nan values
        ind_end_auditory = ind_start;
    end
    
    BEHAVE.SACS_ALL_DATA.neuro_CS_count_visual(  1, counter_sac) = sum(BEHAVE.SACS_ALL_DATA.neuro_CS_visual(  ind_start:ind_end_visual,   counter_sac));
    BEHAVE.SACS_ALL_DATA.neuro_CS_count_auditory(1, counter_sac) = sum(BEHAVE.SACS_ALL_DATA.neuro_CS_auditory(ind_start:ind_end_auditory, counter_sac));
end

fprintf(' --> Completed. \n')

%% Add Neural_Properties
clearvars -except EPHYS BEHAVE
EPHYS.Neural_Properties = struct();
if length(EPHYS.CH_sorted.SS_data.SS_time)>1
    EPHYS.Neural_Properties.SS_num = length(EPHYS.CH_sorted.SS_data.SS_time);
    EPHYS.Neural_Properties.SS_duration = EPHYS.CH_sorted.duration;
    EPHYS.Neural_Properties.SS_firing_rate = EPHYS.Neural_Properties.SS_num / EPHYS.CH_sorted.duration;
    EPHYS.Neural_Properties.SS_time = EPHYS.CH_sorted.SS_data.SS_time;
    EPHYS.Neural_Properties.SS_waveform = nanmean(EPHYS.CH_sorted.SS_data.SS_waveform);
else
    EPHYS.Neural_Properties.SS_num = 0;
    EPHYS.Neural_Properties.SS_duration = EPHYS.CH_sorted.duration;
    EPHYS.Neural_Properties.SS_firing_rate = 0;
    EPHYS.Neural_Properties.SS_time = [];
    EPHYS.Neural_Properties.SS_waveform = zeros(1, size(EPHYS.CH_sorted.CS_data.CS_waveform, 2));
end

if length(EPHYS.CH_sorted.CS_data.CS_time)>1
    EPHYS.Neural_Properties.CS_num = length(EPHYS.CH_sorted.CS_data.CS_time);
    EPHYS.Neural_Properties.CS_firing_rate = EPHYS.Neural_Properties.CS_num / EPHYS.CH_sorted.duration;
    EPHYS.Neural_Properties.CS_time = EPHYS.CH_sorted.CS_data.CS_time;
    EPHYS.Neural_Properties.CS_waveform = nanmean(EPHYS.CH_sorted.CS_data.CS_waveform);
else
    EPHYS.Neural_Properties.CS_num = 0;
    EPHYS.Neural_Properties.CS_firing_rate = 0;
    EPHYS.Neural_Properties.CS_time = [];
    EPHYS.Neural_Properties.CS_waveform = zeros(1, size(EPHYS.CH_sorted.SS_data.SS_waveform, 2));
end
EPHYS.Neural_Properties.Corr_data_CS_inds_span     = nanmean(EPHYS.CH_sorted.Corr_data.CS_inds_span);
EPHYS.Neural_Properties.Corr_data_CS_bin_size_time = nanmean(EPHYS.CH_sorted.Corr_data.CS_bin_size_time);
EPHYS.Neural_Properties.Corr_data_SS_inds_span     = nanmean(EPHYS.CH_sorted.Corr_data.SS_inds_span);
EPHYS.Neural_Properties.Corr_data_SS_bin_size_time = nanmean(EPHYS.CH_sorted.Corr_data.SS_bin_size_time);
EPHYS.Neural_Properties.Corr_data_SS_SSxSS_AUTO    = nanmean(EPHYS.CH_sorted.Corr_data.SS_SSxSS_AUTO);
EPHYS.Neural_Properties.Corr_data_CS_CSxSS_AUTO    = nanmean(EPHYS.CH_sorted.Corr_data.CS_CSxSS_AUTO);

%% outputs
SACS_ALL_DATA = BEHAVE.SACS_ALL_DATA;
Neural_Properties = EPHYS.Neural_Properties;
EXPERIMENT_PARAMS = BEHAVE.EXPERIMENT_PARAMS;

end

%% function combine_sac_mat_files()
function combine_sac_files()
clc; close all;
pCell_list = ESN_build_pCell_list();
path_data_monkey_sorted = uigetdir;
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);
ALL_PCELL_name = ['ALL_PCELL_' num2str(num_pCells)];
if ~exist([path_data_monkey_sorted ALL_PCELL_name filesep 'cell_data'], 'dir')
    mkdir([path_data_monkey_sorted ALL_PCELL_name filesep 'cell_data']);
end
if ~exist([path_data_monkey_sorted ALL_PCELL_name filesep 'cell_figs_behave'], 'dir')
    mkdir([path_data_monkey_sorted ALL_PCELL_name filesep 'cell_figs_behave']);
end
%% Loop over pCells
for counter_pCell = 1 : 1 : num_pCells
    fprintf(['### ' 'Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ###' '\n']);
    num_recording = sum(pCell_list_isstr(counter_pCell, :));
    clearvars data_recordings
    %% Loop over recordings
    for counter_recording = 1 : 1 : num_recording
        %% build plot_data address
        file_name_cell = pCell_list{counter_pCell, counter_recording}; % '190423_142023_01_sorted_ESN_plot_data';
        year_ = file_name_cell(1:2);
        month_ = file_name_cell(3:4);
        day_ = file_name_cell(5:6);
        hour_ = file_name_cell(8:9);
        minute_ = file_name_cell(10:11);
        second_ = file_name_cell(12:13);
        subFolder_month = ['20' year_ '-' month_ filesep];
        subFolder_day = ['20' year_ '-' month_ '-' day_ filesep];
        subFolder_recording = ['20' year_ '-' month_ '-' day_ '_' hour_ '-' minute_ '-' second_ filesep];
        subFolder_figs = ['analyzed_data' filesep];
        %% load recording data
        file_path = [path_data_monkey_sorted subFolder_month subFolder_day subFolder_recording subFolder_figs];
        if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
        data_recording = load([file_path file_name_cell(1:end-4) '_sac' file_name_cell(end-3:end) '.mat'], ...
            'SACS_ALL_DATA', 'Neural_Properties', 'EXPERIMENT_PARAMS');
        if file_name_cell(18) == 's'
            data_recording.id          = file_name_cell(1:16);
        elseif file_name_cell(18) == '2'
            data_recording.id          = file_name_cell(1:18);
        else
            error('Build plot_data_compress: cell id does not follow the standards')
        end
        data_recordings(counter_recording) = data_recording;
    end
    
    %% Init data_cell
    data_cell = struct;
    
    %% id
    data_cell.id = cell(num_recording, 1);
    for counter_recording = 1 : 1 : num_recording
        data_cell.id{counter_recording, 1} = data_recordings(counter_recording).id;
    end
    
    %% EXPERIMENT_PARAMS
    field_names_EXPERIMENT_PARAMS = fieldnames(data_recordings(1).EXPERIMENT_PARAMS);
    for counter_field = 1 : 1 : length(field_names_EXPERIMENT_PARAMS)
        field_name_EXPERIMENT_PARAMS = field_names_EXPERIMENT_PARAMS{counter_field};
        data_cell.EXPERIMENT_PARAMS.(field_name_EXPERIMENT_PARAMS) = cell(num_recording, 1);
        for counter_recording = 1 : 1 : num_recording
            data_cell.EXPERIMENT_PARAMS.(field_name_EXPERIMENT_PARAMS){counter_recording, 1} = data_recordings(counter_recording).EXPERIMENT_PARAMS.(field_name_EXPERIMENT_PARAMS);
        end
    end
    
    %% Neural_Properties
    % Init variables
    data_cell.Neural_Properties.SS_num = 0;
    data_cell.Neural_Properties.SS_duration = 0;
    data_cell.Neural_Properties.SS_firing_rate = 0;
    data_cell.Neural_Properties.CS_num = 0;
    data_cell.Neural_Properties.CS_firing_rate = 0;
    data_cell.Neural_Properties.SS_time = [];
    data_cell.Neural_Properties.CS_time = [];
    variable_names = {'SS_waveform', 'CS_waveform', ...
        'Corr_data_CS_inds_span', 'Corr_data_CS_bin_size_time', 'Corr_data_SS_inds_span', 'Corr_data_SS_bin_size_time', ...
        'Corr_data_SS_SSxSS_AUTO','Corr_data_CS_CSxSS_AUTO'};
    for counter_variable = 1 : length(variable_names)
        variable_name = variable_names{counter_variable};
        data_cell.Neural_Properties.(variable_name) = zeros(size(data_recordings(1).Neural_Properties.(variable_name)));
    end
    % Loop over recordings
    for counter_recording = 1 : 1 : num_recording
        data_cell.Neural_Properties.SS_num = data_cell.Neural_Properties.SS_num + ...
            data_recordings(counter_recording).Neural_Properties.SS_num;
        data_cell.Neural_Properties.SS_duration = data_cell.Neural_Properties.SS_duration + ...
            data_recordings(counter_recording).Neural_Properties.SS_duration;
        data_cell.Neural_Properties.CS_num = data_cell.Neural_Properties.CS_num + ...
            data_recordings(counter_recording).Neural_Properties.CS_num;
        
        data_cell.Neural_Properties.SS_time = vertcat(data_cell.Neural_Properties.SS_time ,...
            data_recordings(counter_recording).Neural_Properties.SS_time);
        data_cell.Neural_Properties.CS_time = vertcat(data_cell.Neural_Properties.CS_time, ...
            data_recordings(counter_recording).Neural_Properties.CS_time);
        % compute weighted average for these variables
        for counter_variable = 1 : length(variable_names)
            variable_name = variable_names{counter_variable};
            if contains(variable_name, 'Corr_data_SS') || contains(variable_name, 'SS_waveform')
                num_spike = data_recordings(counter_recording).Neural_Properties.SS_num;
            elseif contains(variable_name, 'Corr_data_CS') || contains(variable_name, 'CS_waveform')
                num_spike = data_recordings(counter_recording).Neural_Properties.CS_num;
            end
            data_cell.Neural_Properties.(variable_name) = data_cell.Neural_Properties.(variable_name) + ...
                ( data_recordings(counter_recording).Neural_Properties.(variable_name) .* num_spike);
        end
    end
    data_cell.Neural_Properties.SS_firing_rate = ...
            data_cell.Neural_Properties.SS_num ./ data_cell.Neural_Properties.SS_duration;
    data_cell.Neural_Properties.CS_firing_rate = ...
            data_cell.Neural_Properties.CS_num ./ data_cell.Neural_Properties.SS_duration;
    % devide by number of events to compute the weigted average
    for counter_variable = 1 : length(variable_names)
        variable_name = variable_names{counter_variable};
        if contains(variable_name, 'Corr_data_SS') || contains(variable_name, 'SS_waveform')
            num_spike = data_cell.Neural_Properties.SS_num;
        elseif contains(variable_name, 'Corr_data_CS') || contains(variable_name, 'CS_waveform')
            num_spike = data_cell.Neural_Properties.CS_num;
        end
        data_cell.Neural_Properties.(variable_name) = data_cell.Neural_Properties.(variable_name) ./ num_spike;
    end
    
    %% SACS_ALL_DATA
    field_names_SACS_ALL_DATA = fieldnames(data_recordings(1).SACS_ALL_DATA);
    for counter_field = 1 : 1 : length(field_names_SACS_ALL_DATA)
        field_name_SACS_ALL_DATA = field_names_SACS_ALL_DATA{counter_field};
        data_cell.SACS_ALL_DATA.(field_name_SACS_ALL_DATA) = [];
        for counter_recording = 1 : 1 : num_recording
            data_cell.SACS_ALL_DATA.(field_name_SACS_ALL_DATA) = horzcat(...
                data_cell.SACS_ALL_DATA.(field_name_SACS_ALL_DATA), ...
                data_recordings(counter_recording).SACS_ALL_DATA.(field_name_SACS_ALL_DATA));
        end
    end
    
    %% Save data_cell
    id = data_cell.id;
    EXPERIMENT_PARAMS = data_cell.EXPERIMENT_PARAMS;
    Neural_Properties = data_cell.Neural_Properties;
    SACS_ALL_DATA = data_cell.SACS_ALL_DATA;
    cell_name = [data_cell.id{1} '_' 'combine' '_' num2str(length(data_cell.id))];
    fprintf([cell_name ': Saving .mat file ...'])
    save([path_data_monkey_sorted ALL_PCELL_name filesep 'cell_data' filesep cell_name '.mat'],...
        'SACS_ALL_DATA', 'Neural_Properties', 'EXPERIMENT_PARAMS', 'id','-v7.3')
    fprintf(' --> Completed. \n')
    
    %% plot_sac_sorter
    params.cell_name    = cell_name;
    params.duration     = Neural_Properties.SS_duration;
    params.sac_tag_list = EXPERIMENT_PARAMS.sac_tag_list{1};
    params.num_trials   = sum(cell2mat(EXPERIMENT_PARAMS.num_trials));
    fprintf([cell_name ': Saving .png plot ...'])
    plot_sac_sorter(SACS_ALL_DATA, params)
    hFig_ = gcf;
%     saveas(hFig_,[path_data_monkey_sorted ALL_PCELL_name filesep 'cell_figs_behave' filesep cell_name], 'pdf');
    saveas(hFig_,[path_data_monkey_sorted ALL_PCELL_name filesep 'cell_figs_behave' filesep cell_name], 'png');
    close(hFig_)
    fprintf(' --> Completed. \n')
    
end
fprintf('### ALL DONE. ###\n')
end

%% plot_sac_sorter
function plot_sac_sorter(SACS_ALL_DATA, params)
%% Set parameters
amp_edges = -.25 : 0.5 : 15.25;
ang_edges = (-pi-(pi/16)) : (pi/8) : (pi-(pi/16));
react_edges = -12.5: 25 : 512.5;
num_row = 9;
num_col = 9;

%% Init plot
hFig = figure(1);
clf(hFig)
hold on

%% Plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for counter_tag = 1 : 9
if counter_tag == 8
    idx_tag = (SACS_ALL_DATA.tag==8) | (SACS_ALL_DATA.tag==9);
    title_ = [params.sac_tag_list{8} ' & ' params.sac_tag_list{9}];
elseif counter_tag == 9
    idx_tag = (SACS_ALL_DATA.tag==10);
    title_ = params.sac_tag_list{10};
else
    idx_tag = (SACS_ALL_DATA.tag==counter_tag);
    title_ = params.sac_tag_list{counter_tag};
end

axes_minor_nums = reshape(1:num_row*num_col, num_row, num_col)';
axes_main_row = floor((counter_tag - 1) / 3)+1;
axes_main_col = mod(counter_tag, 3); if (axes_main_col==0); axes_main_col=3; end
row1_ = ((axes_main_row-1)*3)+1;
row2_ = ((axes_main_row-1)*3)+2;
row3_ = ((axes_main_row-1)*3)+3;
col1_ = ((axes_main_col-1)*3)+1;
col2_ = ((axes_main_col-1)*3)+2;
col3_ = ((axes_main_col-1)*3)+3;

axes_trace = [axes_minor_nums(row1_,col1_), axes_minor_nums(row1_,col2_), axes_minor_nums(row1_,col3_)...
              axes_minor_nums(row2_,col1_), axes_minor_nums(row2_,col2_), axes_minor_nums(row2_,col3_) ];
axes_amp_dis = axes_minor_nums(row3_,col1_);
axes_ang_dis = axes_minor_nums(row3_,col2_);
axes_react_dis = axes_minor_nums(row3_,col3_);

subplot(num_row,num_col,axes_trace)
hold on
plot(SACS_ALL_DATA.eye_r_px(:,idx_tag), ...
     SACS_ALL_DATA.eye_r_py(:,idx_tag), 'k')
plot(SACS_ALL_DATA.eye_r_px_offset(:,idx_tag), ...
     SACS_ALL_DATA.eye_r_py_offset(:,idx_tag), 'om')
title([title_ ': ' num2str(sum(idx_tag)) ' sac'], 'interpret', 'none');
xlim([-17, 17])
ylim([-15, 15])
% axis equal;

subplot(num_row,num_col,axes_amp_dis)
hold on
histogram(SACS_ALL_DATA.eye_r_amp_m(:,idx_tag), amp_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
histogram(SACS_ALL_DATA.eye_r_amp_m(:,idx_tag), amp_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
xlim([0, 15])
set(gca, 'XTick', 0:3:15)
ylabel('Amplitude')

subplot(num_row,num_col,axes_ang_dis)
polarhistogram(deg2rad(SACS_ALL_DATA.eye_r_ang(:,idx_tag)), (ang_edges), 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
hold on
polarhistogram(deg2rad(SACS_ALL_DATA.eye_r_ang(:,idx_tag)), (ang_edges), 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
set(gca, 'ThetaTick', [])
set(gca, 'RTick', [])
set(gca, 'Title', [])

subplot(num_row,num_col,axes_react_dis)
hold on
histogram(SACS_ALL_DATA.reaction(:,idx_tag)/1000, react_edges/1000, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
histogram(SACS_ALL_DATA.reaction(:,idx_tag)/1000, react_edges/1000, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
xlim([0, 500]/1000)
set(gca, 'XTick', (0:200:500)/1000)
ylabel('Reaction')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add the title info
sgtitle([params.cell_name ', ' ...
    'trial: ' num2str(params.num_trials) ', ' ...
    'sac: ' num2str(length(SACS_ALL_DATA.validity)) ', ' ...
    'dur: ' num2str(params.duration/60,3) 'min' ...
    ], ...
    'interpret', 'none');
ESN_Beautify_Plot(hFig, [13 13], 8)

end

%% function CS_on_analysis()
function CS_on_analysis()
clc; close all;
global ang_edges ang_values range_cell_with_4dir_behave
path_cell_data = uigetdir;
if ~strcmp(path_cell_data(end), filesep);path_cell_data = [path_cell_data filesep];end
cell_file_names = dir([path_cell_data '*_combine_*.mat']);
num_pCells = length(cell_file_names);
%% Loop over pCells
for counter_pCell = 1 : num_pCells
    fprintf(['### ' 'Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ###' '\n']);
    %% load SACS_ALL_DATA
    cell_file_name = cell_file_names(counter_pCell).name;
    load([path_cell_data cell_file_name], 'SACS_ALL_DATA');
    %% compute CS-on
    visual_px_offset = SACS_ALL_DATA.visual_px_offset;
    visual_py_offset = SACS_ALL_DATA.visual_py_offset;
    eye_r_px_onset   = SACS_ALL_DATA.eye_r_px_onset;
    eye_r_py_onset   = SACS_ALL_DATA.eye_r_py_onset;
    delta_x = visual_px_offset - eye_r_px_onset;
    delta_y = visual_py_offset - eye_r_py_onset;
    visual_ang = wrapTo360(atan2d(delta_y, delta_x));
    visual_ang_bin = discretize(visual_ang, ang_edges);
    last_bin_id = length(ang_edges) - 1;
    visual_ang_bin(visual_ang_bin == last_bin_id) = 1; % wrap the circle around
    % 1: 0deg % 2: 45deg % 3: 90deg % 4: 135deg % 5: 180deg % 6: 225deg % 7: 270deg % 8: 315deg
    num_ang_bin = length(ang_edges) - 2;
    tag_bin = unique(SACS_ALL_DATA.tag);
    num_tag_bin = length(tag_bin);
    CS_count  = zeros(num_tag_bin, num_ang_bin);
    sac_count = zeros(num_tag_bin, num_ang_bin);
    CS_prob   = zeros(num_tag_bin, num_ang_bin);
    for counter_tag = 1 : num_tag_bin
        for counter_ang = 1 : num_ang_bin
            idx_tag = (SACS_ALL_DATA.tag == tag_bin(counter_tag));
            idx_ang = (visual_ang_bin == counter_ang);
            CS_count(counter_tag, counter_ang) = sum(SACS_ALL_DATA.neuro_CS_count_visual(1, idx_tag&idx_ang));
            sac_count(counter_tag, counter_ang) = sum(idx_tag&idx_ang);
            CS_prob(counter_tag, counter_ang) = CS_count(counter_tag, counter_ang) ./ sac_count(counter_tag, counter_ang);
        end
    end
    
    r = nansum(CS_prob.* repmat(exp(1i*deg2rad(ang_values)), num_tag_bin, 1) , 2); % compute weighted sum of cos and sin of angles
    CS_ang =  wrapTo360(rad2deg(angle(r))); % Computes the mean direction for circular data.
    CS_prob_sum = nansum(CS_prob,2); % sum of weights
    CS_rho = abs(r) ./ CS_prob_sum; % Computes mean resultant vector length for circular data.
    
    CS_count_avg  = CS_count( 1, :) + CS_count( 4, :);
    sac_count_avg = sac_count(1, :) + sac_count(4, :);
    % 'prim_success' tag 1 % 'corr_success' tag 4
    CS_prob_avg = CS_count_avg ./ sac_count_avg;
    r_avg = nansum(CS_prob_avg.* exp(1i*deg2rad(ang_values)) , 2); % compute weighted sum of cos and sin of angles
    CS_ang_avg =  wrapTo360(rad2deg(angle(r_avg)));
    CS_rho_avg = abs(r_avg) ./ nansum(CS_prob_avg,2);
    
    if (counter_pCell >= range_cell_with_4dir_behave(1)) && (counter_pCell <= range_cell_with_4dir_behave(2))
        idx_CS_on_dir = discretize( ESN_Round(CS_ang_avg, 90.0, 'round'), ang_edges);
    else
        idx_CS_on_dir = discretize(CS_ang_avg, ang_edges);
    end
    
    if idx_CS_on_dir == last_bin_id; idx_CS_on_dir = 1; end
    
    idx_ = idx_CS_on_dir - 1; % make it 0-index format
    if (idx_ == 8); idx_ = 0; end
    idx_CS_tuning = mod((idx_ : 1 : idx_+7), 8) + 1;
    
    %% Build CS_on_data
    CS_on_data.CS_prob = CS_prob;
    CS_on_data.CS_ang  = CS_ang;
    CS_on_data.CS_rho  = CS_rho;
    CS_on_data.CS_prob_avg = CS_prob_avg;
    CS_on_data.CS_ang_avg  = CS_ang_avg;
    CS_on_data.CS_rho_avg  = CS_rho_avg;
    CS_on_data.idx_CS_on_dir  = idx_CS_on_dir;
    CS_on_data.idx_CS_tuning  = idx_CS_tuning;
    CS_on_data.visual_ang_bin       = visual_ang_bin;
    CS_on_data.visual_ang_bin_tuned = idx_CS_tuning(visual_ang_bin);
    
    %% Append CS_on_data to cell_data
    save([path_cell_data cell_file_name], 'CS_on_data', '-append');
end
fprintf('### ALL DONE. ###\n')
end

%% function build_population_data()
function build_population_data()
clc; close all;
global event_type_list inds_span amp_edges vel_edges
path_cell_data = uigetdir;
if ~strcmp(path_cell_data(end), filesep);path_cell_data = [path_cell_data filesep];end
cell_file_names = dir([path_cell_data '*_combine_*.mat']);
cell_file_name = cell_file_names(1).name;
load([path_cell_data cell_file_name], 'SACS_ALL_DATA', 'CS_on_data');
tag_bin = unique(SACS_ALL_DATA.tag);
length_trace = length(inds_span);
num_pCells = length(cell_file_names);
num_tag_bin = length(tag_bin);
num_ang_bin = length(CS_on_data.idx_CS_tuning);
num_amp_bin = length(amp_edges) - 1;
num_vel_bin = length(vel_edges) - 1;
%% Init variables
% STRUCT (SS / CS / VM) -> 10x1 tag struct (amp / vel) -> variable name (onset / vmax / offset / visual) -> 6x8 cell (ampXang / velXang) -> 138x500 double (pCellXtrace)
for counter_tag = 1 : num_tag_bin
    for counter_event_type = 1 : length(event_type_list)
        event_type_name = event_type_list{counter_event_type};
        SS_population_absol.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        SS_population_tuned.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        SS_population_absol.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        SS_population_tuned.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        
        CS_population_absol.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        CS_population_tuned.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        CS_population_absol.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        CS_population_tuned.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        
        VM_population_absol.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        VM_population_tuned.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        VM_population_absol.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        VM_population_tuned.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        
        num_sac_absol.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        num_sac_tuned.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        num_sac_absol.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        num_sac_tuned.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        
        for counter_ang = 1 : num_ang_bin
            for counter_amp = 1 : num_amp_bin
                SS_population_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                SS_population_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                
                CS_population_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                CS_population_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                
                VM_population_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                VM_population_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                
                num_sac_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = zeros(num_pCells, 1);
                num_sac_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = zeros(num_pCells, 1);
            end
            for counter_vel = 1 : num_vel_bin
                SS_population_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                SS_population_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                
                CS_population_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                CS_population_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                
                VM_population_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                VM_population_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                
                num_sac_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = zeros(num_pCells, 1);
                num_sac_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = zeros(num_pCells, 1);
            end
        end
    end
end

%% Loop over pCells
for counter_pCell = 1 : num_pCells
    fprintf(['### ' 'Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ###' '\n']);
    %% load SACS_ALL_DATA
    cell_file_name = cell_file_names(counter_pCell).name;
    load([path_cell_data cell_file_name], 'SACS_ALL_DATA', 'CS_on_data');
    
    %% Compute data
    SACS_amp_bin = discretize(SACS_ALL_DATA.eye_r_amp_m,  amp_edges);
    SACS_vel_bin = discretize(SACS_ALL_DATA.eye_r_vm_max, amp_edges);
    for counter_event_type = 1 : length(event_type_list)
        event_type_name = event_type_list{counter_event_type};
        SACS_ALL_DATA.(['eye_vm' '_' event_type_name]) = sqrt(...
            SACS_ALL_DATA.(['eye_vx' '_' event_type_name]).^2 + SACS_ALL_DATA.(['eye_vy' '_' event_type_name]).^2 );
        for counter_tag = 1 : num_tag_bin
            idx_tag = (SACS_ALL_DATA.tag == tag_bin(counter_tag));
            for counter_ang = 1 : num_ang_bin
                idx_ang_absol = (CS_on_data.visual_ang_bin == counter_ang);
                idx_ang_tuned = (CS_on_data.visual_ang_bin_tuned == counter_ang);
                for counter_amp = 1 : num_amp_bin
                    idx_amp = (SACS_amp_bin == counter_amp);
                    idx_tuned = idx_tag & idx_amp & idx_ang_tuned;
                    idx_absol = idx_tag & idx_amp & idx_ang_absol;
                    
                    SS_population_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['neuro_SS' '_' event_type_name])(:,idx_absol), 2), 1, length_trace);
                    SS_population_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['neuro_SS' '_' event_type_name])(:,idx_tuned), 2), 1, length_trace);
                    CS_population_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['neuro_CS' '_' event_type_name])(:,idx_absol), 2), 1, length_trace);
                    CS_population_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['neuro_CS' '_' event_type_name])(:,idx_tuned), 2), 1, length_trace);
                    VM_population_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['eye_vm' '_' event_type_name])(:,idx_absol), 2), 1, length_trace);
                    VM_population_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['eye_vm' '_' event_type_name])(:,idx_tuned), 2), 1, length_trace);
                    num_sac_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell, :) = ...
                        nansum(idx_absol);
                    num_sac_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell, :) = ...
                        nansum(idx_tuned);
                end
                for counter_vel = 1 : num_vel_bin
                    idx_vel = (SACS_vel_bin == counter_vel);
                    idx_tuned = idx_tag & idx_vel & idx_ang_tuned;
                    idx_absol = idx_tag & idx_vel & idx_ang_absol;
                    
                    SS_population_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['neuro_SS' '_' event_type_name])(:,idx_absol), 2), 1, length_trace);
                    SS_population_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['neuro_SS' '_' event_type_name])(:,idx_tuned), 2), 1, length_trace);
                    CS_population_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['neuro_CS' '_' event_type_name])(:,idx_absol), 2), 1, length_trace);
                    CS_population_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['neuro_CS' '_' event_type_name])(:,idx_tuned), 2), 1, length_trace);
                    VM_population_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['eye_vm' '_' event_type_name])(:,idx_absol), 2), 1, length_trace);
                    VM_population_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['eye_vm' '_' event_type_name])(:,idx_tuned), 2), 1, length_trace);
                    num_sac_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell, :) = ...
                        nansum(idx_absol);
                    num_sac_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell, :) = ...
                        nansum(idx_tuned);
                end
            end
        end
    end

end
fprintf('### ALL DONE. ###\n')
%% Save data
fprintf(['Saving .mat files' ' ...'])
save([path_cell_data '..' filesep 'SS_population_absol' '.mat'], 'SS_population_absol', '-v7.3');
save([path_cell_data '..' filesep 'SS_population_tuned' '.mat'], 'SS_population_tuned', '-v7.3');
save([path_cell_data '..' filesep 'CS_population_absol' '.mat'], 'CS_population_absol', '-v7.3');
save([path_cell_data '..' filesep 'CS_population_tuned' '.mat'], 'CS_population_tuned', '-v7.3');
save([path_cell_data '..' filesep 'VM_population_absol' '.mat'], 'VM_population_absol', '-v7.3');
save([path_cell_data '..' filesep 'VM_population_tuned' '.mat'], 'VM_population_tuned', '-v7.3');
save([path_cell_data '..' filesep 'num_sac_absol' '.mat'], 'num_sac_absol', '-v7.3');
save([path_cell_data '..' filesep 'num_sac_tuned' '.mat'], 'num_sac_tuned', '-v7.3');
fprintf(' --> Completed. \n')
end
