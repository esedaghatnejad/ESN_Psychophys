%% READ ME
% Author: Ehsan Sedaghat-Nejad (esedaghatnejad@gmail.com)
% Analyzing v9_15
% this task is looking at the effct of illuminancy on the learning during a random gain-up gain-down
% task. the type of error and type of image content switch in a stable fashion.
% There are three types of sequences: (1)Perturb+ErrorClamp
% (2)Perturb+Perturb+ErrorClamp (3)Perturb+ErrorClamp+ErrorClamp
% ErrorClamps have been implemented via no feedback

%% Clear workspace
clc;
clear;
close all;
warning('off', 'signal:findpeaks:largeMinPeakHeight');

%% EXPERIMENT_PARAMS
EXPERIMENT_PARAMS.State_ENUM = {'PRE', 'INIT', 'START', 'ITI',...
    'SACCADE', 'END', 'FIXATE', 'ERRCLAMP', 'NEXT_TRIAL', 'BLOCK_FINISHED', 'FIXATE_ERRCLMP'};

%% Get file path
[EXPERIMENT_PARAMS.mat_FileName, EXPERIMENT_PARAMS.mat_PathName] = uigetfile('*.mat', 'Select mat source file');
[~,EXPERIMENT_PARAMS.file_name,~] = fileparts(EXPERIMENT_PARAMS.mat_FileName);

%% Check the file name
if ~isempty(strfind(EXPERIMENT_PARAMS.mat_FileName, 'ANALYZED'))
    error('sac_adapt_8_3 :: ERROR, file name contains ''ANALYZED'' word.');
end

%% Load source file
clearvars -except EXPERIMENT_PARAMS
load([EXPERIMENT_PARAMS.mat_PathName EXPERIMENT_PARAMS.mat_FileName]);

%% Init vars
clearvars('TRIALS');
field_names = fieldnames(DATA_ALL);
trial_number = 0;

%% ## MAIN FOR LOOP ##
for counter_main = 1 : length(field_names)
    fprintf([EXPERIMENT_PARAMS.mat_FileName ': ' field_names{counter_main} ' ...'])
    %% Extract raw data
    if ~exist('counter_main','var'); counter_main = 1; end;
    clearvars -except DATA_ALL TRIALS field_names counter_main trial_number EXPERIMENT_PARAMS
    dataRaw = DATA_ALL.(['block' num2str(counter_main)]);
    
    %% Extract data, interpolate (sample-up), filter
    data.EL.length_data = min([length(dataRaw.right_horizontal_eye), ...
        length(dataRaw.right_vertical_eye) ]);
    data.PC.length_data = min([length(dataRaw.PC_t), ...
        length(dataRaw.PC_target_x), length(dataRaw.PC_target_y) ]);
    
    data.PC.time = dataRaw.PC_t(1:data.PC.length_data); data.PC.time = data.PC.time(:);
    % %% find backward (irregular) time points
    data.PC.temp = data.PC.time;
    data.PC.length_temp = length(find(diff(data.PC.temp)<eps));
    if data.PC.length_temp > 0
        for i = 1 : data.PC.length_temp
            inds = find( diff(data.PC.temp) < eps, 1, 'first');
            data.PC.temp( data.PC.temp(1:inds) >= data.PC.temp(inds+1) ) = NaN;
        end
    end
    data.PC.time = data.PC.temp;
    
    ind_invalid = isnan(data.PC.time);
    data.PC.time(ind_invalid) = []; dataRaw.t(ind_invalid) = [];
    dataRaw.PC_t(ind_invalid) = []; dataRaw.PC_state(ind_invalid) = [];
    dataRaw.PC_target_x(ind_invalid) = []; dataRaw.PC_target_y(ind_invalid) = [];
    
    data.PC.length_data = min([length(dataRaw.PC_t), ...
        length(dataRaw.PC_target_x), length(dataRaw.PC_target_y) ]);
    
    data.EL.time = dataRaw.eyelink_time(1:data.EL.length_data); data.EL.time = data.EL.time(:);
    % %% find backward (irregular) time points
    data.EL.temp = data.EL.time;
    data.EL.length_temp = length(find(diff(data.EL.temp)<eps));
    if data.EL.length_temp > 0
        for i = 1 : data.EL.length_temp
            inds = find( diff(data.EL.temp) < eps, 1, 'first');
            data.EL.temp( data.EL.temp(1:inds) >= data.EL.temp(inds+1) ) = NaN;
        end
    end
    data.EL.time = data.EL.temp;
    
    ind_invalid = isnan(data.EL.time);
    data.EL.time(ind_invalid) = [];
    dataRaw.eyelink_time(ind_invalid) = [];
    dataRaw.eyelink_sample(ind_invalid) = [];
    dataRaw.right_horizontal_eye(ind_invalid) = [];
    dataRaw.right_vertical_eye(ind_invalid) = [];
    dataRaw.right_pupil_area(ind_invalid) = [];
    
    data.EL.length_data = min([length(dataRaw.right_horizontal_eye), ...
        length(dataRaw.right_vertical_eye) ]);
    
    % %% sync the end point of EL and PC wrt PC
    data.EL.time = data.EL.time - (dataRaw.t(1) - dataRaw.PC_t(1));
    % sync the last point of EL and PC wrt PC
    %data.EL.time = data.EL.time .* ( data.PC.time(end) ./ data.EL.time(end) );
    % %% sample up time vector
    data.EL.interp_time_step = 1e-3;
    data.EL.time_1K = (ESN_Round(data.EL.time(2),data.EL.interp_time_step) : data.EL.interp_time_step : data.EL.time(end-1))';
    
    sampling_freq = 1000.0;
    cutoff_freq = 100.0;
    [b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
    
    % %% Load EL data and calculate velocity from position
    data.EL.raw_pupil    = dataRaw.right_pupil_area(1:data.EL.length_data);     data.EL.raw_pupil  = data.EL.raw_pupil(:);
    data.EL.raw_px       = dataRaw.right_horizontal_eye(1:data.EL.length_data); data.EL.raw_px     = data.EL.raw_px(:);
    data.EL.raw_py       = dataRaw.right_vertical_eye(1:data.EL.length_data);   data.EL.raw_py     = data.EL.raw_py(:);
    
    data.EL.interp_pupil = interp1(data.EL.time(~isnan(data.EL.time)), data.EL.raw_pupil(~isnan(data.EL.time)), data.EL.time_1K, 'linear', 'extrap');
    data.EL.interp_px    = interp1(data.EL.time(~isnan(data.EL.time)), data.EL.raw_px(   ~isnan(data.EL.time)), data.EL.time_1K, 'linear', 'extrap');
    data.EL.interp_py    = interp1(data.EL.time(~isnan(data.EL.time)), data.EL.raw_py(   ~isnan(data.EL.time)), data.EL.time_1K, 'linear', 'extrap');
    
    data.EL.px_ = dataRaw.right_horizontal_eye(1:data.EL.length_data); data.EL.px_ = data.EL.px_(:);
    data.EL.py_ = dataRaw.right_vertical_eye(1:data.EL.length_data);   data.EL.py_ = data.EL.py_(:);
    data.EL.px  = interp1(data.EL.time(~isnan(data.EL.time)), data.EL.px_(~isnan(data.EL.time)), data.EL.time_1K, 'linear', 'extrap');
    data.EL.py  = interp1(data.EL.time(~isnan(data.EL.time)), data.EL.py_(~isnan(data.EL.time)), data.EL.time_1K, 'linear', 'extrap');
    data.EL.px_ = data.EL.px; % filtfilt(b_butter,a_butter,data.EL.px); % sgolayfilt(data.EL.px, 3, 21);
    data.EL.py_ = data.EL.py; % filtfilt(b_butter,a_butter,data.EL.py); % sgolayfilt(data.EL.py, 3, 21);
    data.EL.vx_ = diff(data.EL.px_)./diff(data.EL.time_1K); data.EL.vx_=[data.EL.vx_(1) ;data.EL.vx_];
    data.EL.vy_ = diff(data.EL.py_)./diff(data.EL.time_1K); data.EL.vy_=[data.EL.vy_(1) ;data.EL.vy_];
    data.EL.pm_ = sqrt(data.EL.px_.^2 + data.EL.py_.^2);
    data.EL.vm_ = sqrt(data.EL.vx_.^2 + data.EL.vy_.^2);
    
    ind_invalid = isnan(data.EL.pm_) | isnan(data.EL.vm_) | (data.EL.vm_ > 1000.0);
    num_data_points_around_nan_ = 11;
    ind__ = transpose( repmat(find(ind_invalid), 1, num_data_points_around_nan_+1) +...
        repmat(-round(.5*num_data_points_around_nan_):round(.5*num_data_points_around_nan_)-1,...
        length(find(ind_invalid)), 1) );
    ind__ = sort(ind__(:));
    ind__(diff(ind__)==0) = [];
    ind__(ind__<1) = [];
    ind__(ind__>length(ind_invalid)) = [];
    ind_invalid(ind__) = true;
    ind_invalid = ind_invalid | isnan(data.EL.px) | isnan(data.EL.py) | isnan(data.EL.time_1K);
    data.EL.px(ind_invalid) = [];    data.EL.py(ind_invalid) = [];
    time_temp = data.EL.time_1K; time_temp(ind_invalid) = [];
    
    data.EL.px_1K   = interp1(time_temp, data.EL.px, data.EL.time_1K, 'linear', 'extrap');
    data.EL.py_1K   = interp1(time_temp, data.EL.py, data.EL.time_1K, 'linear', 'extrap');
    data.EL.px_filt = filtfilt(b_butter,a_butter,data.EL.px_1K); % sgolayfilt(data.EL.px_1K, 3, 21);
    data.EL.py_filt = filtfilt(b_butter,a_butter,data.EL.py_1K); % sgolayfilt(data.EL.py_1K, 3, 21);
    data.EL.vx_filt = diff(data.EL.px_filt)./diff(data.EL.time_1K); data.EL.vx_filt=[data.EL.vx_filt(1) ;data.EL.vx_filt];
    data.EL.vy_filt = diff(data.EL.py_filt)./diff(data.EL.time_1K); data.EL.vy_filt=[data.EL.vy_filt(1) ;data.EL.vy_filt];
    data.EL.pm_filt = sqrt(data.EL.px_filt.^2 + data.EL.py_filt.^2);
    data.EL.vm_filt = sqrt(data.EL.vx_filt.^2 + data.EL.vy_filt.^2);
    
    % %% Load PC data
    data.PC.px       = dataRaw.PC_target_x(1:data.PC.length_data); data.PC.px = data.PC.px(:);
    data.PC.py       = dataRaw.PC_target_y(1:data.PC.length_data); data.PC.py = data.PC.py(:);
    data.PC.state    = double(dataRaw.PC_state(1:data.PC.length_data) + 1); data.PC.state = data.PC.state(:);
    ind_invalid      = isnan(data.PC.time) | isnan(data.PC.px) | isnan(data.PC.py) | isnan(data.PC.state);
    data.PC.time(ind_invalid) = []; data.PC.px(ind_invalid) = []; data.PC.py(ind_invalid) = []; data.PC.state(ind_invalid) = [];
    data.PC.px_1K    = interp1(data.PC.time(~isnan(data.PC.px)),    data.PC.px(~isnan(data.PC.px)),       data.EL.time_1K, 'nearest', 'extrap');
    data.PC.py_1K    = interp1(data.PC.time(~isnan(data.PC.py)),    data.PC.py(~isnan(data.PC.py)),       data.EL.time_1K, 'nearest', 'extrap');
    data.PC.state_1K = interp1(data.PC.time(~isnan(data.PC.state)), data.PC.state(~isnan(data.PC.state)), data.EL.time_1K, 'nearest', 'extrap');
    
    data.PC.time_START   = dataRaw.PC_time_start; data.PC.time_START = data.PC.time_START(:);
    data.PC.rand_num     = dataRaw.PC_rand_num; data.PC.rand_num = data.PC.rand_num(:);
    
    data.PC.start_x      = dataRaw.PC_start_x; data.PC.start_x = data.PC.start_x(:);
    data.PC.start_y      = dataRaw.PC_start_y; data.PC.start_y = data.PC.start_y(:);
    data.PC.cue_x        = dataRaw.PC_cue_x;   data.PC.cue_x = data.PC.cue_x(:);
    data.PC.cue_y        = dataRaw.PC_cue_y;   data.PC.cue_y = data.PC.cue_y(:);
    data.PC.end_x        = dataRaw.PC_end_x;   data.PC.end_x = data.PC.end_x(:);
    data.PC.end_y        = dataRaw.PC_end_y;
    data.PC.amp_x        = data.PC.end_x - data.PC.start_x;
    data.PC.amp_y        = data.PC.end_y - data.PC.start_y;
    data.PC.err_clmp     = dataRaw.PC_err_clmp;
    
    data.PC.img_start    = dataRaw.PC_img_start; data.PC.img_start = data.PC.img_start(:);
    data.PC.img_cue      = dataRaw.PC_img_cue;   data.PC.img_cue = data.PC.img_cue(:);
    data.PC.img_end      = dataRaw.PC_img_end;   data.PC.img_end = data.PC.img_end(:);
    
    data.PC.FileName     = EXPERIMENT_PARAMS.mat_FileName;
    data.PC.State_ENUM   = EXPERIMENT_PARAMS.State_ENUM;
    data.inds.ind_state_saccade    = find(ismember(EXPERIMENT_PARAMS.State_ENUM, 'SACCADE'));
    data.inds.ind_state_fixate     = find(ismember(EXPERIMENT_PARAMS.State_ENUM, 'FIXATE'));
    data.inds.ind_state_fixate_errclmp = find(ismember(EXPERIMENT_PARAMS.State_ENUM, 'FIXATE_ERRCLMP'));
    
    data.inds.time_START = knnsearch(data.EL.time_1K, data.PC.time_START);
    data.inds.time_START(end+1) = length(data.EL.time_1K);
    
    %% Trial Analyzing
    for counter_trials = 1 : length(data.PC.time_START)
        trial_number = trial_number + 1;
        TRIAL = ESN_Sac_Analyser(data, counter_trials, counter_main, trial_number);
        TRIALS(trial_number) = TRIAL;
    end
    fprintf('--> Completed.\n')
    
end

%% Arrange 'EXPERIMENT'
clearvars -except DATA_ALL TRIALS EXPERIMENT_PARAMS
fprintf([EXPERIMENT_PARAMS.mat_FileName ': Arranging EXPERIMENT ...'])
clearvars('EXPERIMENT'); EXPERIMENT = struct;
field_names_TRIALS = fieldnames(TRIALS);
TRIALS_cell = struct2cell(TRIALS);
for counter_fields = 1 : 1 : length(field_names_TRIALS)
    TRIALS_field_cell = TRIALS_cell(counter_fields,:,:);
    TRIALS_field_cell = reshape(TRIALS_field_cell, 1, []);
    nz = max(cellfun(@numel,TRIALS_field_cell));
    TRIALS_field_mat = cell2mat(cellfun(@(x) vertcat(double(x(:)),NaN(nz-numel(x), 1)),TRIALS_field_cell,'uni',false));
    EXPERIMENT.all.(field_names_TRIALS{counter_fields}) = TRIALS_field_mat;
end
EXPERIMENT.all.sac_validity  = logical(EXPERIMENT.all.sac_validity);
EXPERIMENT.all.sac2_validity = logical(EXPERIMENT.all.sac2_validity);
fprintf('--> Completed.\n')

%% Reviewing Invalid Trials
clearvars -except DATA_ALL TRIALS EXPERIMENT EXPERIMENT_PARAMS
invalid_trials = find(~EXPERIMENT.all.sac_validity);
want_to_review = input(['Num of invalid trials: ', num2str(length(invalid_trials)) ', Do you want to review invalid trials? (0/1) ']);
if want_to_review
    for counter_trial = 1 : length(invalid_trials)
        invalid_trial_ind = invalid_trials(1, counter_trial);
        EXPERIMENT = ESN_Sac_Review(EXPERIMENT, invalid_trial_ind);
    end
end

%% Learning
sac_px_sacFinish_s    = double(EXPERIMENT.all.sac_px_sacFinish_s);
inds_err_clp          = find(EXPERIMENT.all.tgt_err_clp == 1); % extract error clamp indeces
inds_after_pert       = inds_err_clp(2:end); % error clamp after perturbation trial, this method falsely include double error clamp trials which we are going to exclude shortly
inds_before_pert      = inds_err_clp(1:end-1); % error clamp trials before perturbation
sac_learning_signed   = nan(1,length(sac_px_sacFinish_s)); % initialize with nan
sac_learning_signed(:,inds_after_pert) = sac_px_sacFinish_s(:,inds_after_pert) - sac_px_sacFinish_s(:,inds_before_pert); % learning is after minus before
inds_err_clp_double = inds_err_clp(find(diff(inds_err_clp)==1)+1); % find the double error clamp trials and extract the second error clamp as the ones to be eliminated
sac_learning_signed(:,inds_err_clp_double) = nan; % eliminate the double error clamp trials
sac_learning_signed = [sac_learning_signed(:, 2:end) nan]; % shift the index so learning will be for the perturbation trials and not error clamp trials, for double perturbation the learning is considered for the second perturbation
inds_nan = isnan(sac_learning_signed); % trials without learning value
tgt_end_x             = double(EXPERIMENT.all.tgt_end_x);
tgt_cue_x             = double(EXPERIMENT.all.tgt_cue_x);
sac_error             = (sac_px_sacFinish_s) - (tgt_end_x); sac_error(:,inds_nan) = nan; % sensory prediction error during learning trials
sac_error_sign        = sign( (tgt_end_x) - (tgt_cue_x) ); sac_error_sign(:,inds_nan) = nan;
sac_learning_unsigned = sac_learning_signed .* sac_error_sign;

sac_validity    = logical(EXPERIMENT.all.sac_validity);
sac_learning_validity = false(1, length(sac_px_sacFinish_s));
sac_learning_validity(1,inds_after_pert) = sac_validity(1,inds_after_pert) & sac_validity(1,inds_before_pert);
sac_learning_validity(1,inds_err_clp_double) = false;
sac_learning_validity = [sac_learning_validity(:, 2:end) false];

EXPERIMENT.all.sac_learning_signed   = sac_learning_signed;
EXPERIMENT.all.sac_learning_unsigned = sac_learning_unsigned;
EXPERIMENT.all.sac_error             = sac_error;
EXPERIMENT.all.sac_error_sign        = sac_error_sign;
EXPERIMENT.all.sac_learning_validity = sac_learning_validity;

%% Save 'EXPERIMENT'
clearvars -except DATA_ALL TRIALS EXPERIMENT EXPERIMENT_PARAMS
fprintf([EXPERIMENT_PARAMS.mat_FileName ': Saving MAT file ...\n'])
save([EXPERIMENT_PARAMS.mat_PathName EXPERIMENT_PARAMS.file_name '_ANALYZED.mat'], 'EXPERIMENT', 'EXPERIMENT_PARAMS', '-v7.3');
fprintf([EXPERIMENT_PARAMS.mat_FileName ': MAT file saved.\n'])

%% Plot 1 {'sac_reaction', 'sac_vm_max', 'sac2_reaction', 'sac2_vm_max'}
clearvars -except EXPERIMENT EXPERIMENT_PARAMS
variable_list = {'sac_reaction', 'sac_vm_max', 'sac2_reaction', 'sac2_vm_max'};
ylabel_list = {'Primary Reaction (ms)',...
    'Primary V_{max} (deg/s)',...
    'Corrective Reaction (ms)',...
    'Corrective V_{max} (deg/s)'};

clf(figure(1))
hold on

for counter_variable = 1 : length(variable_list)
    %         inds_not_1st_block = true(1, length(EXPERIMENT.all.trial_num));
    %         inds_1st_block(1, [1:sum(EXPERIMENT.all.block_num==1)]) = true; % all
    %         sac_vm_max_mean = mean(EXPERIMENT.all.sac_vm_max((EXPERIMENT.all.sac_validity)));
    %         EXPERIMENT.all.sac_vm_max = EXPERIMENT.all.sac_vm_max ./ sac_vm_max_mean;
    %         EXPERIMENT.all.sac2_vm_max = EXPERIMENT.all.sac2_vm_max ./ sac_vm_max_mean;
    
    inds_not_1st_block = true(1, length(EXPERIMENT.all.trial_num));
    inds_not_1st_block(1, [1:sum(EXPERIMENT.all.block_num==1)]) = false; % all
    inds_err_clp = EXPERIMENT.all.tgt_err_clp == 1;
    
    inds_F_Prim = EXPERIMENT.all.tgt_img_cue == 1;
    inds_N_Prim = EXPERIMENT.all.tgt_img_cue == 2;
    inds_F_Corr = EXPERIMENT.all.tgt_img_end == 1;
    inds_N_Corr = EXPERIMENT.all.tgt_img_end == 2;
    
    
    variable_ = variable_list{counter_variable};
    if (strfind(variable_,'sac_learning_'))
        sac_case_ = 3;
    elseif (strfind(variable_,'sac2_'))
        sac_case_ = 2;
    elseif (strfind(variable_,'sac_'))
        sac_case_ = 1;
    else
        error('cannot find sac_ or sac2_');
    end
    
    if sac_case_ == 1
        validity_ = EXPERIMENT.all.sac_validity;
    elseif sac_case_ == 2
        validity_ = EXPERIMENT.all.sac2_validity;
    elseif sac_case_ == 3
        validity_ = EXPERIMENT.all.sac_learning_validity;
    end
    variable__ = EXPERIMENT.all.(variable_);
    
    validity_ = validity_ & inds_not_1st_block & (~inds_err_clp);
    
    mean_F_Prim = mean(variable__(validity_ & inds_F_Prim));
    mean_N_Prim = mean(variable__(validity_ & inds_N_Prim));
    mean_F_Corr = mean(variable__(validity_ & inds_F_Corr));
    mean_N_Corr = mean(variable__(validity_ & inds_N_Corr));
    mean_FF     = mean(variable__(validity_ & inds_F_Prim & inds_F_Corr));
    mean_NF     = mean(variable__(validity_ & inds_N_Prim & inds_F_Corr));
    mean_FN     = mean(variable__(validity_ & inds_F_Prim & inds_N_Corr));
    mean_NN     = mean(variable__(validity_ & inds_N_Prim & inds_N_Corr));
    stdv_F_Prim = std( variable__(validity_ & inds_F_Prim));
    stdv_N_Prim = std( variable__(validity_ & inds_N_Prim));
    stdv_F_Corr = std( variable__(validity_ & inds_F_Corr));
    stdv_N_Corr = std( variable__(validity_ & inds_N_Corr));
    stdv_FF     = std( variable__(validity_ & inds_F_Prim & inds_F_Corr));
    stdv_NF     = std( variable__(validity_ & inds_N_Prim & inds_F_Corr));
    stdv_FN     = std( variable__(validity_ & inds_F_Prim & inds_N_Corr));
    stdv_NN     = std( variable__(validity_ & inds_N_Prim & inds_N_Corr));
    
    mean_max = max([mean_FF mean_NF mean_FN mean_NN mean_F_Prim mean_N_Prim mean_F_Corr mean_N_Corr]);
    mean_min = min([mean_FF mean_NF mean_FN mean_NN mean_F_Prim mean_N_Prim mean_F_Corr mean_N_Corr]);
    stdv_max = max([stdv_FF stdv_NF stdv_FN stdv_NN stdv_F_Prim stdv_N_Prim stdv_F_Corr stdv_N_Corr]);
    
    y_lim_max_ = mean_max + 1.5 * stdv_max;
    y_lim_min_ = mean_min - 1.5 * stdv_max;
    
    subplot(4,2,2*counter_variable-1)
    hold on
    bar([mean_F_Prim mean_N_Prim NaN mean_F_Corr mean_N_Corr],'FaceColor',[0 .5 .5])
    errorbar(1:5, ...
        [mean_F_Prim mean_N_Prim NaN mean_F_Corr mean_N_Corr],...
        [stdv_F_Prim stdv_N_Prim NaN stdv_F_Corr stdv_N_Corr],...
        '.k', 'LineWidth', 2)
    set(gca, 'XTick', 1:5, 'XTickLabel', {'H-Prim', 'L-Prim', '', 'H-Corr', 'L-Corr'})
    set(gca,'XTickLabelRotation',30);
    ylim([y_lim_min_ y_lim_max_]);
    ylabel(ylabel_list{counter_variable})
    
    subplot(4,2,2*counter_variable)
    hold on
    bar([mean_FF mean_NF mean_FN mean_NN],'FaceColor',[0 .5 .5])
    errorbar(1:4, ...
        [mean_FF mean_NF mean_FN mean_NN],...
        [stdv_FF stdv_NF stdv_FN stdv_NN],...
        '.k', 'LineWidth', 2)
    set(gca, 'XTick', 1:4, 'XTickLabel', {'H-Prim,H-Corr', 'L-Prim,H-Corr', 'H-Prim,L-Corr', 'L-Prim,L-Corr'})
    set(gca,'XTickLabelRotation',30);
    ylim([y_lim_min_ y_lim_max_]);
end



ESN_Beautify_Plot
hFig = gcf;
figure_size  = [8.5 11.0];
paper_margin = [0.2 0.2];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'Clipping', 'off');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'Renderer', 'painters');
set(hFig, 'PaperOrientation', 'portrait');

%% Plot 2 BETWEEN sac_learning_unsigned
clearvars -except EXPERIMENT EXPERIMENT_PARAMS
variable_list = {'sac_learning_unsigned', 'sac_learning_unsigned', 'sac_learning_unsigned'};
ylabel_ = 'Change in Primary Amp (deg)';
title_list = {'Overall',...
    'Gain Down',...
    'Gain Up'};
clf(figure(2))
hold on

for counter_variable = 1 : length(variable_list)
    
    %         inds_1st_block = false(1, length(EXPERIMENT.all.trial_num));
    %         inds_1st_block(1, [1:sum(EXPERIMENT.all.block_num==1)]) = true; % all
    %         sac_vm_max_mean = mean(EXPERIMENT.all.sac_vm_max((EXPERIMENT.all.sac_validity)));
    %         EXPERIMENT.all.sac_vm_max = EXPERIMENT.all.sac_vm_max ./ sac_vm_max_mean;
    %         EXPERIMENT.all.sac2_vm_max = EXPERIMENT.all.sac2_vm_max ./ sac_vm_max_mean;
    
    inds_not_1st_block = true(1, length(EXPERIMENT.all.trial_num));
    inds_not_1st_block(1, [1:sum(EXPERIMENT.all.block_num==1)]) = false; % all
    inds_err_clp = EXPERIMENT.all.tgt_err_clp == 1;
    
    inds_Up     = EXPERIMENT.all.sac_error_sign == 1;
    inds_Down   = EXPERIMENT.all.sac_error_sign == -1;
    inds_F_Prim = EXPERIMENT.all.tgt_img_cue == 1;
    inds_N_Prim = EXPERIMENT.all.tgt_img_cue == 2;
    inds_F_Corr = EXPERIMENT.all.tgt_img_end == 1;
    inds_N_Corr = EXPERIMENT.all.tgt_img_end == 2;
    
    if (strfind(title_list{counter_variable},'Overall'))
        %inds_not_1st_block = inds_not_1st_block;
    elseif (strfind(title_list{counter_variable},'Down'))
        inds_not_1st_block = inds_not_1st_block&inds_Down;
    elseif (strfind(title_list{counter_variable},'Up'))
        inds_not_1st_block = inds_not_1st_block&inds_Up;
    else
        error('cannot find sac_ or sac2_');
    end
    
    variable_ = variable_list{counter_variable};
    if (strfind(variable_,'sac_learning_'))
        sac_case_ = 3;
    elseif (strfind(variable_,'sac2_'))
        sac_case_ = 2;
    elseif (strfind(variable_,'sac_'))
        sac_case_ = 1;
    else
        error('cannot find sac_ or sac2_');
    end
    
    if sac_case_ == 1
        validity_ = EXPERIMENT.all.sac_validity;
    elseif sac_case_ == 2
        validity_ = EXPERIMENT.all.sac2_validity;
    elseif sac_case_ == 3
        validity_ = EXPERIMENT.all.sac_learning_validity;
    end
    variable__ = EXPERIMENT.all.(variable_);
    
    validity_ = validity_ & inds_not_1st_block & (~inds_err_clp);
%     validity_ = EXPERIMENT.all.sac_validity & EXPERIMENT.all.sac_learning_validity & inds_not_1st_block & (~inds_err_clp);
    
    mean_F_Prim = mean(variable__(validity_ & inds_F_Prim));
    mean_N_Prim = mean(variable__(validity_ & inds_N_Prim));
    mean_F_Corr = mean(variable__(validity_ & inds_F_Corr));
    mean_N_Corr = mean(variable__(validity_ & inds_N_Corr));
    mean_FF     = mean(variable__(validity_ & inds_F_Prim & inds_F_Corr));
    mean_NF     = mean(variable__(validity_ & inds_N_Prim & inds_F_Corr));
    mean_FN     = mean(variable__(validity_ & inds_F_Prim & inds_N_Corr));
    mean_NN     = mean(variable__(validity_ & inds_N_Prim & inds_N_Corr));
    stdv_F_Prim = std( variable__(validity_ & inds_F_Prim));
    stdv_N_Prim = std( variable__(validity_ & inds_N_Prim));
    stdv_F_Corr = std( variable__(validity_ & inds_F_Corr));
    stdv_N_Corr = std( variable__(validity_ & inds_N_Corr));
    stdv_FF     = std( variable__(validity_ & inds_F_Prim & inds_F_Corr));
    stdv_NF     = std( variable__(validity_ & inds_N_Prim & inds_F_Corr));
    stdv_FN     = std( variable__(validity_ & inds_F_Prim & inds_N_Corr));
    stdv_NN     = std( variable__(validity_ & inds_N_Prim & inds_N_Corr));
    
    mean_max = max([mean_FF mean_NF mean_FN mean_NN mean_F_Prim mean_N_Prim mean_F_Corr mean_N_Corr]);
    mean_min = min([mean_FF mean_NF mean_FN mean_NN mean_F_Prim mean_N_Prim mean_F_Corr mean_N_Corr]);
    stdv_max = max([stdv_FF stdv_NF stdv_FN stdv_NN stdv_F_Prim stdv_N_Prim stdv_F_Corr stdv_N_Corr]);
    
    y_lim_max_ = mean_max + 1.5 * stdv_max;
    y_lim_min_ = mean_min - 1.5 * stdv_max;
    
    subplot(3,2,2*counter_variable-1)
    hold on
    bar([mean_F_Prim mean_N_Prim NaN mean_F_Corr mean_N_Corr],'FaceColor',[0 .5 .5])
    errorbar(1:5, ...
        [mean_F_Prim mean_N_Prim NaN mean_F_Corr mean_N_Corr],...
        [stdv_F_Prim stdv_N_Prim NaN stdv_F_Corr stdv_N_Corr],...
        '.k', 'LineWidth', 2)
    set(gca, 'XTick', 1:5, 'XTickLabel', {'H-Prim', 'L-Prim', '', 'H-Corr', 'L-Corr'})
    set(gca,'XTickLabelRotation',30);
    ylim([y_lim_min_ y_lim_max_]);
    title(title_list{counter_variable})
    ylabel(ylabel_)
    
    subplot(3,2,2*counter_variable)
    hold on
    bar([mean_FF mean_NF mean_FN mean_NN],'FaceColor',[0 .5 .5])
    errorbar(1:4, ...
        [mean_FF mean_NF mean_FN mean_NN],...
        [stdv_FF stdv_NF stdv_FN stdv_NN],...
        '.k', 'LineWidth', 2)
    set(gca, 'XTick', 1:4, 'XTickLabel', {'H-Prim,H-Corr', 'L-Prim,H-Corr', 'H-Prim,L-Corr', 'L-Prim,L-Corr'})
    set(gca,'XTickLabelRotation',30);
    ylim([y_lim_min_ y_lim_max_]);
    title(title_list{counter_variable})
end



ESN_Beautify_Plot
hFig = gcf;
figure_size  = [8.5 11.0];
paper_margin = [0.2 0.2];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'Clipping', 'off');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'Renderer', 'painters');
set(hFig, 'PaperOrientation', 'portrait');

%% Saveas figs and Close
fprintf([EXPERIMENT_PARAMS.mat_FileName ': Saving plots ...\n']);
figHandles = findobj('Type','figure');
for counter_fig = 1 : length(figHandles)
    saveas(figHandles(counter_fig),[EXPERIMENT_PARAMS.mat_PathName 'Fig' num2str(figHandles(counter_fig).Number) ' ' EXPERIMENT_PARAMS.file_name], 'pdf');
    saveas(figHandles(counter_fig),[EXPERIMENT_PARAMS.mat_PathName 'Fig' num2str(figHandles(counter_fig).Number) ' ' EXPERIMENT_PARAMS.file_name], 'png');
    close(figHandles(counter_fig));
end
close all
fprintf([EXPERIMENT_PARAMS.mat_FileName ': Plots saved.\n']);

