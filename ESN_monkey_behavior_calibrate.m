function ESN_monkey_behavior_calibrate

%% Clear and Load data
clc;clear;close all;
% subject_name = 'Mirza';
% load('random_corrective_saccades_094041_REDUCED.mat');
% subject_name = 'Ramon';
% load('random_corrective_saccades_170202_REDUCED.mat');
% subject_name = 'Ramon2';
% load('random_corrective_saccades_090024_REDUCED.mat');



%% Get variables
clearvars -except TRIALS_DATA SACS_PRIM_DATA SACS_CORR_DATA EXPERIMENT_PARAMS subject_name
sac1_r_px = SACS_PRIM_DATA.eye_r_px;
sac1_r_py = SACS_PRIM_DATA.eye_r_py;
start_x = TRIALS_DATA.start_x;
start_y = TRIALS_DATA.start_y;
sac1_r_px_start = SACS_PRIM_DATA.eye_r_px_start;
sac1_r_py_start = SACS_PRIM_DATA.eye_r_py_start;
sac1_r_px_finish = SACS_PRIM_DATA.eye_r_px_finish;
sac1_r_py_finish = SACS_PRIM_DATA.eye_r_py_finish;
sac1_valid = SACS_PRIM_DATA.validity;
cue_x = TRIALS_DATA.cue_x;
cue_y = TRIALS_DATA.cue_y;

sac2_r_px = SACS_CORR_DATA.eye_r_px;
sac2_r_py = SACS_CORR_DATA.eye_r_py;
sac2_r_px_finish = SACS_CORR_DATA.eye_r_px_finish;
sac2_r_py_finish = SACS_CORR_DATA.eye_r_py_finish;
sac2_valid = SACS_CORR_DATA.validity;
end_x = TRIALS_DATA.end_x;
end_y = TRIALS_DATA.end_y;

inds_valid = SACS_PRIM_DATA.validity & SACS_CORR_DATA.validity;
error_prim = sqrt( (SACS_PRIM_DATA.eye_r_px_finish - TRIALS_DATA.cue_x).^2 + (SACS_PRIM_DATA.eye_r_py_finish - TRIALS_DATA.cue_y).^2 );
error_corr = sqrt( (SACS_CORR_DATA.eye_r_px_finish - TRIALS_DATA.end_x).^2 + (SACS_CORR_DATA.eye_r_py_finish - TRIALS_DATA.end_y).^2 );
inds_valid = inds_valid & (error_prim<3) & (error_corr<3);

%% Finding Invali Trials
sac1_error = sqrt(((sac1_r_px_finish-cue_x).^2)+((sac1_r_py_finish-cue_y).^2));
sac2_error = sqrt(((sac2_r_px_finish-end_x).^2)+((sac2_r_py_finish-end_y).^2));
sac1_reaction = SACS_PRIM_DATA.reaction;
sac2_reaction = SACS_CORR_DATA.reaction;
sac12_valid = sac1_valid & sac2_valid & (sac1_error<1.5) & (sac2_error<1.5) & (sac1_reaction < 500) & (sac2_reaction<500);

%% Plot-1 eye trajectories
h_fig(1) = figure(1);
clf(h_fig(1))
subplot(1,2,1)
hold on
inds_valid_plot = inds_valid;
plot(sac1_r_px(:,inds_valid_plot), sac1_r_py(:,inds_valid_plot), 'color', [0.7 0.7 0.7])
plot(sac2_r_px(:,inds_valid_plot), sac2_r_py(:,inds_valid_plot), 'color', [0.5 0.5 0.5])
plot(sac1_r_px_start(:,inds_valid_plot), sac1_r_py_start(:,inds_valid_plot), 'og','MarkerFaceColor', 'g')
plot(sac1_r_px_finish(:,inds_valid_plot), sac1_r_py_finish(:,inds_valid_plot), 'ob','MarkerFaceColor', 'b')
plot(sac2_r_px_finish(:,inds_valid_plot), sac2_r_py_finish(:,inds_valid_plot), 'or','MarkerFaceColor', 'r')
plot(start_x, start_y, 'dk','MarkerFaceColor', 'k')
plot(end_x, end_y, 'dk','MarkerFaceColor', 'k')
plot(cue_x, cue_y, 'dk','MarkerFaceColor', 'k')
xlabel('Horizontal Eye Posiotion (deg)')
ylabel('Vertical Eye Posiotion (deg)')
title('Before Calibration')
% ESN_Beautify_Plot
axis equal
xlim([-12 12])
ylim([-12 12])

%% Finding the Calibration Matrix
eye_position = [sac2_r_px_finish(sac12_valid)', sac2_r_py_finish(sac12_valid)', ones(size(sac2_r_px_finish(sac12_valid)'))];
tgt_position = [end_x(sac12_valid)',            end_y(sac12_valid)',            ones(size(end_x(sac12_valid)'))];
calib_matrix = ((eye_position' * eye_position)^(-1)) * (eye_position') * tgt_position;

%% Recalibrate the data
sac2_r_finish = [sac2_r_px_finish(:), sac2_r_py_finish(:), ones(size(sac2_r_px_finish(:)))];
sac2_r_finish = sac2_r_finish * calib_matrix;
sac2_r_px_finish = sac2_r_finish(:,1)';
sac2_r_py_finish = sac2_r_finish(:,2)';
sac1_r_finish = [sac1_r_px_finish(:), sac1_r_py_finish(:), ones(size(sac1_r_px_finish(:)))];
sac1_r_finish = sac1_r_finish * calib_matrix;
sac1_r_px_finish = sac1_r_finish(:,1)';
sac1_r_py_finish = sac1_r_finish(:,2)';
sac1_r_start = [sac1_r_px_start(:), sac1_r_py_start(:), ones(size(sac1_r_px_start(:)))];
sac1_r_start = sac1_r_start * calib_matrix;
sac1_r_px_start = sac1_r_start(:,1)';
sac1_r_py_start = sac1_r_start(:,2)';
num_trials = length(cue_x);

for counter_trials = 1 : num_trials
    sac2_r = [sac2_r_px(:,counter_trials), sac2_r_py(:,counter_trials), ones(size(sac2_r_px(:,counter_trials)))];
    sac2_r = sac2_r * calib_matrix;
    sac2_r_px(:,counter_trials) = sac2_r(:,1);
    sac2_r_py(:,counter_trials) = sac2_r(:,2);
    sac1_r = [sac1_r_px(:,counter_trials), sac1_r_py(:,counter_trials), ones(size(sac1_r_px(:,counter_trials)))];
    sac1_r = sac1_r * calib_matrix;
    sac1_r_px(:,counter_trials) = sac1_r(:,1);
    sac1_r_py(:,counter_trials) = sac1_r(:,2);
end

%% Plot-1 eye trajectories
% h_fig(2) = figure(2);
% clf(h_fig(2))
subplot(1,2,2)
hold on
inds_valid_plot = inds_valid;
plot(sac1_r_px(:,inds_valid_plot), sac1_r_py(:,inds_valid_plot), 'color', [0.7 0.7 0.7])
plot(sac2_r_px(:,inds_valid_plot), sac2_r_py(:,inds_valid_plot), 'color', [0.5 0.5 0.5])
plot(sac1_r_px_start(:,inds_valid_plot), sac1_r_py_start(:,inds_valid_plot), 'og','MarkerFaceColor', 'g')
plot(sac1_r_px_finish(:,inds_valid_plot), sac1_r_py_finish(:,inds_valid_plot), 'ob','MarkerFaceColor', 'b')
plot(sac2_r_px_finish(:,inds_valid_plot), sac2_r_py_finish(:,inds_valid_plot), 'or','MarkerFaceColor', 'r')
plot(start_x, start_y, 'dk','MarkerFaceColor', 'k')
plot(end_x, end_y, 'dk','MarkerFaceColor', 'k')
plot(cue_x, cue_y, 'dk','MarkerFaceColor', 'k')
xlabel('Horizontal Eye Posiotion (deg)')
ylabel('Vertical Eye Posiotion (deg)')
title('After Calibration')
axis equal
xlim([-12 12])
ylim([-12 12])

ESN_Beautify_Plot(h_fig, [20 10])


%% Plot-2 reaction time histogram
h_fig(2) = figure(2);
clf(h_fig(2))
hold on
reaction_edges = 80 : 10 : 400;
histogram(SACS_PRIM_DATA.reaction(sac12_valid), reaction_edges, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2)
histogram(SACS_CORR_DATA.reaction(sac12_valid), reaction_edges, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2)
xlabel('Reaction Time (ms)')
ylabel('Probability')
ESN_Beautify_Plot

disp( [ 'mean sac1_reaction ' num2str(mean(sac1_reaction(sac12_valid))) ] )
disp( [ 'stdv sac1_reaction ' num2str(std(sac1_reaction(sac12_valid))) ] )
disp( [ 'mean sac2_reaction ' num2str(mean(sac2_reaction(sac12_valid))) ] )
disp( [ 'stdv sac2_reaction ' num2str(std(sac2_reaction(sac12_valid))) ] )

%% Num trial performance
num_trial_performance = ...
length(TRIALS_DATA.time_state_next_trial(:)) ./ ...
length(cell2mat(TRIALS_DATA.time_state_sac_detect_off(:))) .* ...
100;
disp( [ 'Num trial performance: ' num2str(num_trial_performance) ...
    '%, ' num2str(length(TRIALS_DATA.time_state_next_trial(:))) ...
    ' out of ' num2str(length(cell2mat(TRIALS_DATA.time_state_sac_detect_off(:)))) ...
    '.' ] );

%% Plot-3 error histogram
sac1_error = sqrt(((sac1_r_px_finish-cue_x).^2)+((sac1_r_py_finish-cue_y).^2));
sac2_error = sqrt(((sac2_r_px_finish-end_x).^2)+((sac2_r_py_finish-end_y).^2));
sac12_valid = sac1_valid & sac2_valid & (sac1_error<2.5) & (sac2_error<2.5);
h_fig(3) = figure(3);
clf(h_fig(3))
hold on
sac_error_edges = 0 : 0.05 : 2.0;
histogram(sac1_error(sac12_valid), sac_error_edges, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2)
histogram(sac2_error(sac12_valid), sac_error_edges, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2)
xlabel('Distance to Target (deg)')
ylabel('Probability')
ESN_Beautify_Plot

disp( [ 'mean sac1_error ' num2str(mean(sac1_error(sac12_valid))) ] )
disp( [ 'stdv sac1_error ' num2str(std(sac1_error(sac12_valid))) ] )
disp( [ 'mean sac2_error ' num2str(mean(sac2_error(sac12_valid))) ] )
disp( [ 'stdv sac2_error ' num2str(std(sac2_error(sac12_valid))) ] )

%% Save Plot
saveas(h_fig(1),[subject_name '_eye_trajectory' ], 'pdf');
saveas(h_fig(2),[subject_name '_reaction_time'], 'pdf');
saveas(h_fig(3),[subject_name '_distance'], 'pdf');