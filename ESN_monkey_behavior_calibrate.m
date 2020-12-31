function ESN_monkey_behavior_calibrate(file_path_list, file_name_list)
%% Clear
%{
clc;clear;close all;
file_path_list = cell(0,1);
file_name_list = cell(0,1);
file_path_list{1,1} = ...
'/Users/ehsan/Downloads/data_184d_sorted/2019-04/2019-04-03/2019-04-03_14-13-13/analyzed_data';
file_name_list{1,1} = ...
'190403_141313_ANALYZED.mat';
file_path_list{2,1} = ...
'/Users/ehsan/Downloads/data_184d_sorted/2019-04/2019-04-03/2019-04-03_14-38-42/analyzed_data';
file_name_list{2,1} = ...
'190403_143842_ANALYZED.mat';
%}
%% Handle inputs
if nargin < 1
    [file_name_,file_path_] = uigetfile([pwd filesep '*_ANALYZED.mat'], 'Select behavior file');
    file_path_list = {file_path_};
    file_name_list = {file_name_};
end
%% Load data
fprintf(['Loading files ' ' ... ']);
num_files = length(file_name_list);
clearvars data_ANALYZED
for counter_file = 1 : 1 : num_files
    file_path_ = file_path_list{counter_file};
    file_name_ = file_name_list{counter_file};
    if ~strcmp(file_path_(end), filesep);file_path_ = [file_path_ filesep];end
    file_full_path_ = [file_path_ file_name_];
    data_ANALYZED(counter_file) = load(file_full_path_);
end
fprintf(' --> Completed. \n');

%% RE-CALIBRATE
fprintf(['Calibrating files ' ' ... ']);
% perform the re-calibration 5 times to make sure we have used all of the data
num_recalibration_iteration = 5;
for counter_iteration =  1 : 1 : num_recalibration_iteration
    %% Get variables
    sac_prim_validity = false(1,0);
    sac_corr_validity = false(1,0);
    sac_prim_reaction = [];
    sac_corr_reaction = [];
    sac_prim_px_start = [];
    sac_prim_py_start = [];
    sac_prim_px_finish = [];
    sac_prim_py_finish = [];
    sac_corr_px_finish = [];
    sac_corr_py_finish = [];
    tgt_start_x = [];
    tgt_start_y = [];
    tgt_cue_x = [];
    tgt_cue_y = [];
    tgt_end_x = [];
    tgt_end_y = [];
    for counter_file = 1 : 1 : num_files
        sac_prim_validity  = horzcat(sac_prim_validity,  data_ANALYZED(counter_file).SACS_PRIM_DATA.validity);
        sac_corr_validity  = horzcat(sac_corr_validity,  data_ANALYZED(counter_file).SACS_CORR_DATA.validity);
        sac_prim_reaction  = horzcat(sac_prim_reaction,  data_ANALYZED(counter_file).SACS_PRIM_DATA.reaction);
        sac_corr_reaction  = horzcat(sac_corr_reaction,  data_ANALYZED(counter_file).SACS_CORR_DATA.reaction);
        sac_prim_px_start  = horzcat(sac_prim_px_start,  data_ANALYZED(counter_file).SACS_PRIM_DATA.eye_r_px_start);
        sac_prim_py_start  = horzcat(sac_prim_py_start,  data_ANALYZED(counter_file).SACS_PRIM_DATA.eye_r_py_start);
        sac_prim_px_finish = horzcat(sac_prim_px_finish, data_ANALYZED(counter_file).SACS_PRIM_DATA.eye_r_px_finish);
        sac_prim_py_finish = horzcat(sac_prim_py_finish, data_ANALYZED(counter_file).SACS_PRIM_DATA.eye_r_py_finish);
        sac_corr_px_finish = horzcat(sac_corr_px_finish, data_ANALYZED(counter_file).SACS_CORR_DATA.eye_r_px_finish);
        sac_corr_py_finish = horzcat(sac_corr_py_finish, data_ANALYZED(counter_file).SACS_CORR_DATA.eye_r_py_finish);
        tgt_start_x = horzcat(tgt_start_x, data_ANALYZED(counter_file).TRIALS_DATA.start_x);
        tgt_start_y = horzcat(tgt_start_y, data_ANALYZED(counter_file).TRIALS_DATA.start_y);
        tgt_cue_x   = horzcat(tgt_cue_x,   data_ANALYZED(counter_file).TRIALS_DATA.cue_x  );
        tgt_cue_y   = horzcat(tgt_cue_y,   data_ANALYZED(counter_file).TRIALS_DATA.cue_y  );
        tgt_end_x   = horzcat(tgt_end_x,   data_ANALYZED(counter_file).TRIALS_DATA.end_x  );
        tgt_end_y   = horzcat(tgt_end_y,   data_ANALYZED(counter_file).TRIALS_DATA.end_y  );
    end
    
    %% Finding Invali Trials
    start_error = sqrt(((sac_prim_px_start -tgt_start_x).^2)+((sac_prim_py_start -tgt_start_y).^2));
    prim_error  = sqrt(((sac_prim_px_finish-tgt_cue_x  ).^2)+((sac_prim_py_finish-tgt_cue_y  ).^2));
    corr_error  = sqrt(((sac_corr_px_finish-tgt_end_x  ).^2)+((sac_corr_py_finish-tgt_end_y  ).^2));
    sac_validity = sac_prim_validity & sac_corr_validity & (start_error<1.5) & (prim_error<1.5) & (corr_error<1.5) & (sac_prim_reaction < 500) & (sac_corr_reaction<500);
    
    %% Finding the Calibration Matrix
    sac_prim_start = [sac_prim_px_start(sac_validity)', sac_prim_py_start(sac_validity)', ones(size(sac_prim_px_start(sac_validity)'))];
    tgt_start      = [tgt_start_x(sac_validity)',       tgt_start_y(sac_validity)',       ones(size(tgt_start_x(sac_validity)'))];
    num_data_points = round(sum(sac_validity) / 8) + 1;
    
    [sac_prim_start,idx] = datasample(sac_prim_start,num_data_points,1,'Replace',false);
    tgt_start = tgt_start(idx, :);
    
    sac_corr_finish = [sac_corr_px_finish(sac_validity)', sac_corr_py_finish(sac_validity)', ones(size(sac_corr_px_finish(sac_validity)'))];
    tgt_end         = [tgt_end_x(sac_validity)',          tgt_end_y(sac_validity)',          ones(size(tgt_end_x(sac_validity)'))];
    
    eye_position = [sac_corr_finish; sac_prim_start];
    tgt_position = [tgt_end;         tgt_start];
    calib_matrix = ((eye_position' * eye_position)^(-1)) * (eye_position') * tgt_position;
    
    %% Recalibrate the data
    main_field_list = {'TRIALS_DATA', 'SACS_PRIM_DATA', 'SACS_CORR_DATA'};
    for counter_file = 1 : 1 : num_files
        num_trials = length(data_ANALYZED(counter_file).TRIALS_DATA.start_x);
        for counter_main_field = 1 : length(main_field_list)
            main_field = main_field_list{counter_main_field};
            if counter_main_field == 1
                variable_list = { ...
                    'eye_r_px', 'eye_r_py'; ...
                    'eye_r_px_filt', 'eye_r_py_filt'; ...
                    };
            else
                variable_list = { ...
                    'eye_r_px', 'eye_r_py'; ...
                    'eye_r_px_centered', 'eye_r_py_centered'; ...
                    'eye_r_px_start', 'eye_r_py_start'; ...
                    'eye_r_px_finish', 'eye_r_py_finish'; ...
                    'eye_r_px_start_centered', 'eye_r_py_start_centered'; ...
                    'eye_r_px_finish_centered', 'eye_r_py_finish_centered'; ...
                    };
            end
            for counter_variable = 1 : size(variable_list, 1)
                variable_x = data_ANALYZED(counter_file).(main_field).(variable_list{counter_variable, 1});
                variable_y = data_ANALYZED(counter_file).(main_field).(variable_list{counter_variable, 2});
                for counter_trials = 1 : num_trials
                    if counter_main_field == 1
                        variable_ = [variable_x{:,counter_trials}, variable_y{:,counter_trials}, ones(size(variable_x{:,counter_trials}))];
                        variable_calib = variable_ * calib_matrix;
                        variable_x{:,counter_trials} = variable_calib(:,1);
                        variable_y{:,counter_trials} = variable_calib(:,2);
                    else
                        variable_ = [variable_x(:,counter_trials), variable_y(:,counter_trials), ones(size(variable_x(:,counter_trials)))];
                        variable_calib = variable_ * calib_matrix;
                        variable_x(:,counter_trials) = variable_calib(:,1);
                        variable_y(:,counter_trials) = variable_calib(:,2);
                    end
                end
                data_ANALYZED(counter_file).(main_field).(variable_list{counter_variable, 1}) = variable_x;
                data_ANALYZED(counter_file).(main_field).(variable_list{counter_variable, 2}) = variable_y;
            end
        end
    end
end

fprintf(' --> Completed. \n');

%% Chech validity again
for counter_file = 1 : 1 : num_files
    %% Get variables
    sac_prim_reaction  = data_ANALYZED(counter_file).SACS_PRIM_DATA.reaction;
    sac_corr_reaction  = data_ANALYZED(counter_file).SACS_CORR_DATA.reaction;
    sac_prim_vm_max    = data_ANALYZED(counter_file).SACS_PRIM_DATA.eye_r_vm_max;
    sac_corr_vm_max    = data_ANALYZED(counter_file).SACS_CORR_DATA.eye_r_vm_max;
    sac_prim_px_start  = data_ANALYZED(counter_file).SACS_PRIM_DATA.eye_r_px_start;
    sac_prim_py_start  = data_ANALYZED(counter_file).SACS_PRIM_DATA.eye_r_py_start;
    sac_prim_px_finish = data_ANALYZED(counter_file).SACS_PRIM_DATA.eye_r_px_finish;
    sac_prim_py_finish = data_ANALYZED(counter_file).SACS_PRIM_DATA.eye_r_py_finish;
    sac_corr_px_finish = data_ANALYZED(counter_file).SACS_CORR_DATA.eye_r_px_finish;
    sac_corr_py_finish = data_ANALYZED(counter_file).SACS_CORR_DATA.eye_r_py_finish;
    tgt_start_x = data_ANALYZED(counter_file).TRIALS_DATA.start_x;
    tgt_start_y = data_ANALYZED(counter_file).TRIALS_DATA.start_y;
    tgt_cue_x   = data_ANALYZED(counter_file).TRIALS_DATA.cue_x;
    tgt_cue_y   = data_ANALYZED(counter_file).TRIALS_DATA.cue_y;
    tgt_end_x   = data_ANALYZED(counter_file).TRIALS_DATA.end_x;
    tgt_end_y   = data_ANALYZED(counter_file).TRIALS_DATA.end_y;
    
    %% re-compute validity
    start_error = sqrt(((sac_prim_px_start -tgt_start_x).^2)+((sac_prim_py_start -tgt_start_y).^2));
    prim_error  = sqrt(((sac_prim_px_finish-tgt_cue_x  ).^2)+((sac_prim_py_finish-tgt_cue_y  ).^2));
    corr_error  = sqrt(((sac_corr_px_finish-tgt_end_x  ).^2)+((sac_corr_py_finish-tgt_end_y  ).^2));
    sac_prim_validity = ...
        (start_error<3.0) & ...
        (prim_error <3.0) & ...
        (sac_prim_reaction < 500) & (sac_prim_reaction > 50) & ...
        (sac_prim_vm_max < 1000)  & (sac_prim_vm_max > 50);
    sac_corr_validity = ...
        (corr_error<3.0) & ...
        (sac_corr_reaction < 500) & (sac_corr_reaction > 50) & ...
        (sac_corr_vm_max < 1000)  & (sac_corr_vm_max > 50);
    data_ANALYZED(counter_file).SACS_PRIM_DATA.validity = sac_prim_validity;
    data_ANALYZED(counter_file).SACS_CORR_DATA.validity = sac_corr_validity;
    
    %% re-compute amp
    data_ANALYZED(counter_file).SACS_PRIM_DATA.eye_r_amp_x = (data_ANALYZED(counter_file).SACS_PRIM_DATA.eye_r_px_finish - data_ANALYZED(counter_file).SACS_PRIM_DATA.eye_r_px_start);
    data_ANALYZED(counter_file).SACS_PRIM_DATA.eye_r_amp_y = (data_ANALYZED(counter_file).SACS_PRIM_DATA.eye_r_py_finish - data_ANALYZED(counter_file).SACS_PRIM_DATA.eye_r_py_start);
    data_ANALYZED(counter_file).SACS_PRIM_DATA.eye_r_amp_m = (sqrt(data_ANALYZED(counter_file).SACS_PRIM_DATA.eye_r_amp_x.^2+data_ANALYZED(counter_file).SACS_PRIM_DATA.eye_r_amp_y.^2));
    data_ANALYZED(counter_file).SACS_CORR_DATA.eye_r_amp_x = (data_ANALYZED(counter_file).SACS_CORR_DATA.eye_r_px_finish - data_ANALYZED(counter_file).SACS_CORR_DATA.eye_r_px_start);
    data_ANALYZED(counter_file).SACS_CORR_DATA.eye_r_amp_y = (data_ANALYZED(counter_file).SACS_CORR_DATA.eye_r_py_finish - data_ANALYZED(counter_file).SACS_CORR_DATA.eye_r_py_start);
    data_ANALYZED(counter_file).SACS_CORR_DATA.eye_r_amp_m = (sqrt(data_ANALYZED(counter_file).SACS_CORR_DATA.eye_r_amp_x.^2+data_ANALYZED(counter_file).SACS_CORR_DATA.eye_r_amp_y.^2));
    
    %% plot outcome
    hFig = figure(1);
    clf(hFig);
    hold('on')
    tgt_start_unique = unique([tgt_start_x(:) tgt_start_y(:)], 'rows');
    tgt_cue_unique = unique([tgt_cue_x(:) tgt_cue_y(:)], 'rows');
    tgt_end_unique = unique([tgt_end_x(:) tgt_end_y(:)], 'rows');
    plot(sac_prim_px_start(sac_prim_validity), sac_prim_py_start(sac_prim_validity), 'ok','MarkerFaceColor', 'g')
    plot(sac_prim_px_finish(sac_prim_validity), sac_prim_py_finish(sac_prim_validity), 'ok','MarkerFaceColor', 'b')
    plot(sac_corr_px_finish(sac_corr_validity), sac_corr_py_finish(sac_corr_validity), 'ok','MarkerFaceColor', 'r')
    plot(tgt_start_unique(:,1), tgt_start_unique(:,2), 'om','MarkerFaceColor', 'k')
    plot(tgt_cue_unique(:,1), tgt_cue_unique(:,2), 'om','MarkerFaceColor', 'k')
    plot(tgt_end_unique(:,1), tgt_end_unique(:,2), 'om','MarkerFaceColor', 'k')
    title(file_name_list{counter_file}(1:end-4), 'interpret', 'none');
    xlabel('Horz position (deg)')
    ylabel('Vert position (deg)')
    axis('equal')
    ESN_Beautify_Plot(hFig, [5 4])

    %% Save _ANALYZED.mat Data to disk
    EXPERIMENT_PARAMS = data_ANALYZED(counter_file).EXPERIMENT_PARAMS;
    TRIALS_DATA = data_ANALYZED(counter_file).TRIALS_DATA;
    SACS_PRIM_DATA = data_ANALYZED(counter_file).SACS_PRIM_DATA;
    SACS_CORR_DATA = data_ANALYZED(counter_file).SACS_CORR_DATA;
    file_path_ = file_path_list{counter_file};
    if ~strcmp(file_path_(end), filesep);file_path_ = [file_path_ filesep];end
    file_name_ = file_name_list{counter_file};
    fprintf([file_name_ ': Saving ...'])
    save([file_path_ file_name_], ...
        'EXPERIMENT_PARAMS', 'TRIALS_DATA', 'SACS_PRIM_DATA', 'SACS_CORR_DATA', '-v7.3');
    fprintf(' --> Completed. \n')
    
    %% Save _REDUCED.mat Data to disk
    rmfields_list = {'eye_l_vm_filt', 'eye_l_vy_filt', 'eye_l_vx_filt', 'eye_l_py_filt', 'eye_l_px_filt', ...
        'eye_r_vm_filt', 'eye_r_vy_filt', 'eye_r_vx_filt', 'eye_r_py_filt', 'eye_r_px_filt', ...
        'time', 'time_1K', 'target_visible', 'reward', 'tgt_py', 'tgt_px', 'time_tgt', ...
        'eye_l_vm', 'eye_r_vm', 'eye_l_vy', 'eye_l_vx', 'eye_r_vy', 'eye_r_vx', ...
        'eye_l_py', 'eye_l_px', 'eye_r_py', 'eye_r_px', 'time_eyelink', 'inds_invalid', 'inds_trial'};
    TRIALS_DATA = rmfield(TRIALS_DATA,rmfields_list);
    file_name_REDUCED_ = [file_name_(1:13) '_REDUCED.mat'];
    fprintf([file_name_REDUCED_ ': Saving ...'])
    save([file_path_ file_name_REDUCED_], ...
        'EXPERIMENT_PARAMS', 'TRIALS_DATA', 'SACS_PRIM_DATA', 'SACS_CORR_DATA', '-v7.3');
    fprintf(' --> Completed. \n')
    
    %% Save Fig
    file_name_plot_ = file_name_(1:end-4);
    file_path_plot_ = [file_path_ '..' filesep 'analyzed_figs' filesep];
    saveas(hFig,[file_path_plot_ file_name_plot_ '_plot'], 'pdf');
    saveas(hFig,[file_path_plot_ file_name_plot_ '_plot'], 'png');
    close(hFig)
end


