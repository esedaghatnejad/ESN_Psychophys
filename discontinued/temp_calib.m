%% load _ANALYZED
[file_name_,file_path_] = uigetfile([pwd filesep '*_ANALYZED.mat'], 'Select behavior file');
file_path_list = {file_path_};
file_name_list = {file_name_};
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

data_ANALYZED_orig = data_ANALYZED;

%% Extract RAW data Recalibrate the data
affine_tf = zeros(3,3);
affine_tf(1,1) = data.eyelink_right_affine_transformation.p0;
affine_tf(1,2) = data.eyelink_right_affine_transformation.p1;
affine_tf(1,3) = data.eyelink_right_affine_transformation.p2;
affine_tf(2,1) = data.eyelink_right_affine_transformation.p3;
affine_tf(2,2) = data.eyelink_right_affine_transformation.p4;
affine_tf(2,3) = data.eyelink_right_affine_transformation.p5;
affine_tf(3,1) = 0;
affine_tf(3,2) = 0;
affine_tf(3,3) = 1;

calib_matrix = (affine_tf^-1)';

% Extract RAW data Recalibrate the data
data_ANALYZED.SACS_PRIM_DATA.start_x = data_ANALYZED.TRIALS_DATA.start_x;
data_ANALYZED.SACS_PRIM_DATA.start_y = data_ANALYZED.TRIALS_DATA.start_y;
data_ANALYZED.SACS_PRIM_DATA.cue_x   = data_ANALYZED.TRIALS_DATA.cue_x;
data_ANALYZED.SACS_PRIM_DATA.cue_y   = data_ANALYZED.TRIALS_DATA.cue_y;
data_ANALYZED.SACS_PRIM_DATA.end_x   = data_ANALYZED.TRIALS_DATA.end_x;
data_ANALYZED.SACS_PRIM_DATA.end_y   = data_ANALYZED.TRIALS_DATA.end_y;

data_ANALYZED.SACS_CORR_DATA.start_x = data_ANALYZED.TRIALS_DATA.start_x;
data_ANALYZED.SACS_CORR_DATA.start_y = data_ANALYZED.TRIALS_DATA.start_y;
data_ANALYZED.SACS_CORR_DATA.cue_x   = data_ANALYZED.TRIALS_DATA.cue_x;
data_ANALYZED.SACS_CORR_DATA.cue_y   = data_ANALYZED.TRIALS_DATA.cue_y;
data_ANALYZED.SACS_CORR_DATA.end_x   = data_ANALYZED.TRIALS_DATA.end_x;
data_ANALYZED.SACS_CORR_DATA.end_y   = data_ANALYZED.TRIALS_DATA.end_y;

% Extract RAW data Recalibrate the data
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
                'start_x', 'start_y'; ...
                'cue_x', 'cue_y'; ...
                'end_x', 'end_y'; ...
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

%
data_ANALYZED.TRIALS_DATA.start_x = data_ANALYZED.SACS_PRIM_DATA.start_x;
data_ANALYZED.TRIALS_DATA.start_y = data_ANALYZED.SACS_PRIM_DATA.start_y;
data_ANALYZED.TRIALS_DATA.cue_x   = data_ANALYZED.SACS_PRIM_DATA.cue_x;
data_ANALYZED.TRIALS_DATA.cue_y   = data_ANALYZED.SACS_PRIM_DATA.cue_y;
data_ANALYZED.TRIALS_DATA.end_x   = data_ANALYZED.SACS_PRIM_DATA.end_x;
data_ANALYZED.TRIALS_DATA.end_y   = data_ANALYZED.SACS_PRIM_DATA.end_y;

data_ANALYZED_raw = data_ANALYZED;

%% Calculate Calibration Matrix
% trial_list_start = [...
% 2,3,4,5,6,7,8,9,10,12,13,14,15,18,19,21,22,23,25,26,28,29,30,33,34,35,36,37,40,41,43,44,45,46,47,48   ...
% ];
% 
% trial_list_prim = [...
% 2,4,6,9,10,12,13,14,18,19,21,22,23,25,26,28,29,30,33,34,35,36,37,40,41,43,44,45,46,47,48,49   ...
% ];
% 
% trial_list_corr = [...
% 2,3,4,5,6,7,8,9,10,12,13,14,15,18,19,21,22,23,25,28,29,31,33,35,36,39,40,43,47,48   ...
% ];

sac_validity = data_ANALYZED_raw.SACS_PRIM_DATA.validity & data_ANALYZED_raw.SACS_CORR_DATA.validity;
sac_validity = find(sac_validity);
trial_list_start = sac_validity;
trial_list_prim = sac_validity;
trial_list_corr = sac_validity;

eye_position_s = [data_ANALYZED_raw.SACS_PRIM_DATA.eye_r_px_start(:,trial_list_start)', ...
                  data_ANALYZED_raw.SACS_PRIM_DATA.eye_r_py_start(:,trial_list_start)', ...
                  ones(length(trial_list_start), 1)];
tgt_position_s = [data_ANALYZED_orig.TRIALS_DATA.start_x(:,trial_list_start)', ...
                  data_ANALYZED_orig.TRIALS_DATA.start_y(:,trial_list_start)', ...
                  ones(length(trial_list_start), 1)];

eye_position_p = [data_ANALYZED_raw.SACS_PRIM_DATA.eye_r_px_finish(:,trial_list_prim)', ...
                  data_ANALYZED_raw.SACS_PRIM_DATA.eye_r_py_finish(:,trial_list_prim)', ...
                  ones(length(trial_list_prim), 1)];
tgt_position_p = [data_ANALYZED_orig.TRIALS_DATA.cue_x(:,trial_list_prim)', ...
                  data_ANALYZED_orig.TRIALS_DATA.cue_y(:,trial_list_prim)', ...
                  ones(length(trial_list_prim), 1)];

eye_position_c = [data_ANALYZED_raw.SACS_CORR_DATA.eye_r_px_finish(:,trial_list_corr)', ...
                  data_ANALYZED_raw.SACS_CORR_DATA.eye_r_py_finish(:,trial_list_corr)', ...
                  ones(length(trial_list_corr), 1)];
tgt_position_c = [data_ANALYZED_orig.TRIALS_DATA.cue_x(:,trial_list_corr)', ...
                  data_ANALYZED_orig.TRIALS_DATA.cue_y(:,trial_list_corr)', ...
                  ones(length(trial_list_corr), 1)];

eye_position = [eye_position_s; eye_position_p; eye_position_c];
tgt_position = [tgt_position_s; tgt_position_p; tgt_position_c];
calib_matrix = ((eye_position' * eye_position)^(-1)) * (eye_position') * tgt_position;


fprintf([ ...
     '\n' ...
     '   \"right_calibration\": {\n' ...
     '     \"p0\": ' num2str(calib_matrix(1,1),16) ',\n' ...
     '     \"p1\": ' num2str(calib_matrix(2,1),16) ',\n' ...
     '     \"p2\": ' num2str(calib_matrix(3,1),16) ',\n' ...
     '     \"p3\": ' num2str(calib_matrix(1,2),16) ',\n' ...
     '     \"p4\": ' num2str(calib_matrix(2,2),16) ',\n' ...
     '     \"p5\": ' num2str(calib_matrix(3,2),16) ',\n' ...
])

%% Recalibrate the data
data_ANALYZED = data_ANALYZED_raw;
data_ANALYZED.TRIALS_DATA.start_x = data_ANALYZED_orig.TRIALS_DATA.start_x;
data_ANALYZED.TRIALS_DATA.start_y = data_ANALYZED_orig.TRIALS_DATA.start_y;
data_ANALYZED.TRIALS_DATA.cue_x   = data_ANALYZED_orig.TRIALS_DATA.cue_x;
data_ANALYZED.TRIALS_DATA.cue_y   = data_ANALYZED_orig.TRIALS_DATA.cue_y;
data_ANALYZED.TRIALS_DATA.end_x   = data_ANALYZED_orig.TRIALS_DATA.end_x;
data_ANALYZED.TRIALS_DATA.end_y   = data_ANALYZED_orig.TRIALS_DATA.end_y;

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
                'start_x', 'start_y'; ...
                'cue_x', 'cue_y'; ...
                'end_x', 'end_y'; ...
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

% data_ANALYZED.TRIALS_DATA.start_x = data_ANALYZED.SACS_PRIM_DATA.start_x;
% data_ANALYZED.TRIALS_DATA.start_y = data_ANALYZED.SACS_PRIM_DATA.start_y;
% data_ANALYZED.TRIALS_DATA.cue_x   = data_ANALYZED.SACS_PRIM_DATA.cue_x;
% data_ANALYZED.TRIALS_DATA.cue_y   = data_ANALYZED.SACS_PRIM_DATA.cue_y;
% data_ANALYZED.TRIALS_DATA.end_x   = data_ANALYZED.SACS_PRIM_DATA.end_x;
% data_ANALYZED.TRIALS_DATA.end_y   = data_ANALYZED.SACS_PRIM_DATA.end_y;

data_ANALYZED.TRIALS_DATA.start_x = data_ANALYZED_orig.TRIALS_DATA.start_x;
data_ANALYZED.TRIALS_DATA.start_y = data_ANALYZED_orig.TRIALS_DATA.start_y;
data_ANALYZED.TRIALS_DATA.cue_x   = data_ANALYZED_orig.TRIALS_DATA.cue_x;
data_ANALYZED.TRIALS_DATA.cue_y   = data_ANALYZED_orig.TRIALS_DATA.cue_y;
data_ANALYZED.TRIALS_DATA.end_x   = data_ANALYZED_orig.TRIALS_DATA.end_x;
data_ANALYZED.TRIALS_DATA.end_y   = data_ANALYZED_orig.TRIALS_DATA.end_y;

%% Plot & Review trials
dataset_ = data_ANALYZED;
% dataset_ = data_ANALYZED_raw;
% dataset_ = data_ANALYZED_orig;

% for counter_trial = [51 67 70 112 116 124]
for counter_trial = 1:50
    clf(figure(4))
    hold on
    
    
    plot(dataset_.SACS_PRIM_DATA.eye_r_px(:,counter_trial), ...
         dataset_.SACS_PRIM_DATA.eye_r_py(:,counter_trial), '-b')
    plot(dataset_.SACS_PRIM_DATA.eye_r_px_start(:,counter_trial), ...
         dataset_.SACS_PRIM_DATA.eye_r_py_start(:,counter_trial), 'og','MarkerFaceColor', 'b', 'linewidth', 2)
    plot(dataset_.SACS_PRIM_DATA.eye_r_px_finish(:,counter_trial), ...
         dataset_.SACS_PRIM_DATA.eye_r_py_finish(:,counter_trial), 'ob','MarkerFaceColor', 'b', 'linewidth', 2)
    
    plot(dataset_.SACS_CORR_DATA.eye_r_px(:,counter_trial), ...
         dataset_.SACS_CORR_DATA.eye_r_py(:,counter_trial), '-r')
    plot(dataset_.SACS_CORR_DATA.eye_r_px_start(:,counter_trial), ...
         dataset_.SACS_CORR_DATA.eye_r_py_start(:,counter_trial), 'og','MarkerFaceColor', 'r', 'linewidth', 2)
    plot(dataset_.SACS_CORR_DATA.eye_r_px_finish(:,counter_trial), ...
         dataset_.SACS_CORR_DATA.eye_r_py_finish(:,counter_trial), 'or','MarkerFaceColor', 'r', 'linewidth', 2)
    
    plot(dataset_.TRIALS_DATA.start_x(:,counter_trial), ...
         dataset_.TRIALS_DATA.start_y(:,counter_trial), 'og','MarkerFaceColor', 'c', 'linewidth', 2)
    plot(dataset_.TRIALS_DATA.cue_x(:,counter_trial), ...
         dataset_.TRIALS_DATA.cue_y(:,counter_trial), 'ob','MarkerFaceColor', 'c', 'linewidth', 2)
    plot(dataset_.TRIALS_DATA.end_x(:,counter_trial), ...
         dataset_.TRIALS_DATA.end_y(:,counter_trial), 'or','MarkerFaceColor', 'c', 'linewidth', 2)
    title(['Trial no.: ' num2str(counter_trial)])
%     xlim([0 10000]); ylim([4000 12000])
    xlim([-8 12]); ylim([-8 8])
    
    waitforbuttonpress
end

%% Plot error field
clf(figure(1))
hold on

dataset_ = data_ANALYZED;
% dataset_ = data_ANALYZED_raw;
% dataset_ = data_ANALYZED_orig;

plot(dataset_.SACS_PRIM_DATA.eye_r_px_start(:,trial_list_start)', ...
     dataset_.SACS_PRIM_DATA.eye_r_py_start(:,trial_list_start)', 'og','MarkerFaceColor', 'g', 'linewidth', 2)

plot(dataset_.TRIALS_DATA.start_x(:,trial_list_start)', ...
     dataset_.TRIALS_DATA.start_y(:,trial_list_start)', 'dg', 'MarkerSize', 10, 'MarkerFaceColor', 'c','linewidth', 2)
 
plot(dataset_.SACS_PRIM_DATA.eye_r_px_finish(:,trial_list_prim)', ...
     dataset_.SACS_PRIM_DATA.eye_r_py_finish(:,trial_list_prim)', 'ob','MarkerFaceColor', 'b', 'linewidth', 2)

plot(dataset_.TRIALS_DATA.cue_x(:,trial_list_prim)', ...
     dataset_.TRIALS_DATA.cue_y(:,trial_list_prim)', 'db', 'MarkerSize', 10, 'MarkerFaceColor', 'c','linewidth', 2)
 
plot(dataset_.SACS_CORR_DATA.eye_r_px_finish(:,trial_list_corr)', ...
     dataset_.SACS_CORR_DATA.eye_r_py_finish(:,trial_list_corr)', 'or','MarkerFaceColor', 'r', 'linewidth', 2)
 
plot(dataset_.TRIALS_DATA.end_x(:,trial_list_corr)', ...
     dataset_.TRIALS_DATA.end_y(:,trial_list_corr)', 'dr', 'MarkerSize', 10 , 'MarkerFaceColor', 'c', 'linewidth', 2)

 plot([dataset_.SACS_PRIM_DATA.eye_r_px_start(:,trial_list_start);dataset_.TRIALS_DATA.start_x(:,trial_list_start)], ...
      [dataset_.SACS_PRIM_DATA.eye_r_py_start(:,trial_list_start);dataset_.TRIALS_DATA.start_y(:,trial_list_start)], 'k')
 
 plot([dataset_.SACS_PRIM_DATA.eye_r_px_finish(:,trial_list_prim);dataset_.TRIALS_DATA.cue_x(:,trial_list_prim)], ...
      [dataset_.SACS_PRIM_DATA.eye_r_py_finish(:,trial_list_prim);dataset_.TRIALS_DATA.cue_y(:,trial_list_prim)],'k' )
 
 plot([dataset_.SACS_CORR_DATA.eye_r_px_finish(:,trial_list_corr);dataset_.TRIALS_DATA.end_x(:,trial_list_corr)], ...
      [dataset_.SACS_CORR_DATA.eye_r_py_finish(:,trial_list_corr);dataset_.TRIALS_DATA.end_y(:,trial_list_corr)],'k' )
 
xlim([-8 12]); ylim([-8 8])

%% Use data_ANALYZED_calib to find calibration matrix
sac_validity = find(sac_validity);
trial_list_start = sac_validity;
trial_list_prim = sac_validity;
trial_list_corr = sac_validity;

eye_position_s = [data_ANALYZED_raw.SACS_PRIM_DATA.eye_r_px_start(:,trial_list_start)', ...
                  data_ANALYZED_raw.SACS_PRIM_DATA.eye_r_py_start(:,trial_list_start)', ...
                  ones(length(trial_list_start), 1)];
tgt_position_s = [data_ANALYZED_calib.SACS_PRIM_DATA.eye_r_px_start(:,trial_list_start)', ...
                  data_ANALYZED_calib.SACS_PRIM_DATA.eye_r_py_start(:,trial_list_start)', ...
                  ones(length(trial_list_start), 1)];

eye_position_p = [data_ANALYZED_raw.SACS_PRIM_DATA.eye_r_px_finish(:,trial_list_prim)', ...
                  data_ANALYZED_raw.SACS_PRIM_DATA.eye_r_py_finish(:,trial_list_prim)', ...
                  ones(length(trial_list_prim), 1)];
tgt_position_p = [data_ANALYZED_calib.SACS_PRIM_DATA.eye_r_px_finish(:,trial_list_prim)', ...
                  data_ANALYZED_calib.SACS_PRIM_DATA.eye_r_py_finish(:,trial_list_prim)', ...
                  ones(length(trial_list_prim), 1)];

eye_position_c = [data_ANALYZED_raw.SACS_CORR_DATA.eye_r_px_finish(:,trial_list_corr)', ...
                  data_ANALYZED_raw.SACS_CORR_DATA.eye_r_py_finish(:,trial_list_corr)', ...
                  ones(length(trial_list_corr), 1)];
tgt_position_c = [data_ANALYZED_calib.SACS_CORR_DATA.eye_r_px_finish(:,trial_list_corr)', ...
                  data_ANALYZED_calib.SACS_CORR_DATA.eye_r_py_finish(:,trial_list_corr)', ...
                  ones(length(trial_list_corr), 1)];

eye_position = [eye_position_s; eye_position_p; eye_position_c];
tgt_position = [tgt_position_s; tgt_position_p; tgt_position_c];
calib_matrix = ((eye_position' * eye_position)^(-1)) * (eye_position') * tgt_position;

fprintf([ ...
     '\n' ...
     '   \"right_calibration\": {\n' ...
     '     \"p0\": ' num2str(calib_matrix(1,1),16) ',\n' ...
     '     \"p1\": ' num2str(calib_matrix(2,1),16) ',\n' ...
     '     \"p2\": ' num2str(calib_matrix(3,1),16) ',\n' ...
     '     \"p3\": ' num2str(calib_matrix(1,2),16) ',\n' ...
     '     \"p4\": ' num2str(calib_matrix(2,2),16) ',\n' ...
     '     \"p5\": ' num2str(calib_matrix(3,2),16) ',\n' ...
])



