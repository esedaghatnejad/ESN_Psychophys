%% function ESN_plot_neural_modulation
function ESN_plot_neural_modulation(num_data_set)
if nargin < 1
    num_data_set = 1;
end
%% Build EPHYS_ and BEHAVE_ for each single dataset
clearvars EPHYS BEHAVE
for counter_dataset = 1 : 1 : num_data_set
    [EPHYS_(counter_dataset), BEHAVE_(counter_dataset)] = build_EPHYS_BEHAVE_single_dataset;
end

%% Plot-1 SS & CS train cue_present 
clearvars raster_data_
for counter_dataset = 1 : 1 : num_data_set
    inds_event           = EPHYS_(counter_dataset).CH_EVE.BEHAVE_inds_cue_present;
    raster_data_(counter_dataset).data = single_dataset_raster(EPHYS_(counter_dataset), ...
        BEHAVE_(counter_dataset), inds_event, 'prim');
end

raster_data = concatenate_dataset(raster_data_, 'data', @vertcat);

plot_data.fig_num_               = 1;
plot_data.xlabel_text_raster_    = {'Time relative to cue presentation (ms)', 'Directions based on primary sac'};
plot_data.xlabel_text_CS_probab_ = {'CS probability based on [0 +200]ms'};
plot_data.inds_span = EPHYS_(1).CH_EVE.inds_span_cue_present;
fig_handle_(plot_data.fig_num_)  = plot_rasters_data(raster_data, plot_data);

%% Plot-2 SS & CS train primSac_onset 
clearvars raster_data_
for counter_dataset = 1 : 1 : num_data_set
    inds_event           = EPHYS_(counter_dataset).CH_EVE.BEHAVE_inds_primSac_onset;
    raster_data_(counter_dataset).data = single_dataset_raster(EPHYS_(counter_dataset), ...
        BEHAVE_(counter_dataset), inds_event, 'prim');
end

raster_data = concatenate_dataset(raster_data_, 'data', @vertcat);

plot_data.fig_num_               = 2;
plot_data.xlabel_text_raster_    = {'Time relative to prim sac onset (ms)', 'Directions based on primary sac'};
plot_data.xlabel_text_CS_probab_ = {'CS probability based on [-200 0]ms'};
plot_data.inds_span = EPHYS_(1).CH_EVE.inds_span_primSac_onset;
fig_handle_(plot_data.fig_num_)  = plot_rasters_data(raster_data, plot_data);

%% Plot-3 SS & CS train primSac_offset 
clearvars raster_data_
for counter_dataset = 1 : 1 : num_data_set
    inds_event           = EPHYS_(counter_dataset).CH_EVE.BEHAVE_inds_primSac_offset;
    raster_data_(counter_dataset).data = single_dataset_raster(EPHYS_(counter_dataset), ...
        BEHAVE_(counter_dataset), inds_event, 'corr');
end

raster_data = concatenate_dataset(raster_data_, 'data', @vertcat);

plot_data.fig_num_               = 3;
plot_data.xlabel_text_raster_    = {'Time relative to prim sac offset (ms)', 'Directions based on secondary sac'};
plot_data.xlabel_text_CS_probab_ = {'CS probability based on [0 +200]ms'};
plot_data.inds_span = EPHYS_(1).CH_EVE.inds_span_primSac_offset;
fig_handle_(plot_data.fig_num_)  = plot_rasters_data(raster_data, plot_data);

%% Plot-4 SS & CS train corrSac_onset 
clearvars raster_data_
for counter_dataset = 1 : 1 : num_data_set
    inds_event           = EPHYS_(counter_dataset).CH_EVE.BEHAVE_inds_corrSac_onset;
    raster_data_(counter_dataset).data = single_dataset_raster(EPHYS_(counter_dataset), ...
        BEHAVE_(counter_dataset), inds_event, 'corr');
end

raster_data = concatenate_dataset(raster_data_, 'data', @vertcat);

plot_data.fig_num_               = 4;
plot_data.xlabel_text_raster_    = {'Time relative to corr sac onset (ms)', 'Directions based on secondary sac'};
plot_data.xlabel_text_CS_probab_ = {'CS probability based on [-200 0]ms'};
plot_data.inds_span = EPHYS_(1).CH_EVE.inds_span_corrSac_onset;
fig_handle_(plot_data.fig_num_)  = plot_rasters_data(raster_data, plot_data);

%%
%{
%% Plot-1 SS & CS train cue_present 
clearvars -except EPHYS BEHAVE fig_handle_
% inds of interest
inds_event           = EPHYS.CH_EVE.BEHAVE_inds_cue_present;
% CS and SS data
SS_train_aligned     = EPHYS.CH_EVE.EPHYS_SS_train_aligned;
CS_train_aligned     = EPHYS.CH_EVE.EPHYS_CS_train_aligned;
inds_span            = EPHYS.CH_EVE.inds_span_cue_present(1,:);
% eye velocity
BEHAVE_eye_r_vm_filt = EPHYS.CH_EVE.BEHAVE_eye_r_vm_filt;
inds_valid = BEHAVE.SACS_PRIM_DATA.validity & BEHAVE.SACS_CORR_DATA.validity;
error_prim = sqrt( (BEHAVE.SACS_PRIM_DATA.eye_r_px_finish - BEHAVE.TRIALS_DATA.cue_x).^2 + (BEHAVE.SACS_PRIM_DATA.eye_r_py_finish - BEHAVE.TRIALS_DATA.cue_y).^2 );
error_corr = sqrt( (BEHAVE.SACS_CORR_DATA.eye_r_px_finish - BEHAVE.TRIALS_DATA.end_x).^2 + (BEHAVE.SACS_CORR_DATA.eye_r_py_finish - BEHAVE.TRIALS_DATA.end_y).^2 );
inds_valid = inds_valid & (error_prim<3) & (error_corr<3);
% inds directions
inds_top   = ((BEHAVE.TRIALS_DATA.cue_x-BEHAVE.TRIALS_DATA.start_x) == 0) & ...
             ((BEHAVE.TRIALS_DATA.cue_y-BEHAVE.TRIALS_DATA.start_y) >  0) & ...
             inds_valid;
inds_left  = ((BEHAVE.TRIALS_DATA.cue_x-BEHAVE.TRIALS_DATA.start_x) <  0) & ...
             ((BEHAVE.TRIALS_DATA.cue_y-BEHAVE.TRIALS_DATA.start_y) == 0) & ...
             inds_valid;
inds_right = ((BEHAVE.TRIALS_DATA.cue_x-BEHAVE.TRIALS_DATA.start_x) >  0) & ...
             ((BEHAVE.TRIALS_DATA.cue_y-BEHAVE.TRIALS_DATA.start_y) == 0) & ...
             inds_valid;
inds_down  = ((BEHAVE.TRIALS_DATA.cue_x-BEHAVE.TRIALS_DATA.start_x) == 0) & ...
             ((BEHAVE.TRIALS_DATA.cue_y-BEHAVE.TRIALS_DATA.start_y) <  0) & ...
             inds_valid;
% CS Probab
range_inds_probability = 101:300;
train_data_logic = CS_train_aligned(inds_event(inds_top,:));
prob_top   = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_top);
train_data_logic = CS_train_aligned(inds_event(inds_left,:));
prob_left  = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_left);
train_data_logic = CS_train_aligned(inds_event(inds_right,:));
prob_right = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_right);
train_data_logic = CS_train_aligned(inds_event(inds_down,:));
prob_down  = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_down);
amplitude = [prob_right prob_top prob_left prob_down prob_right];
% plot xlim and ylim
range_SS_Firing = [0 200];
range_vm        = [-50 600];
% xlabel text
xlabel_text_raster_ = {'Time relative to cue presentation (ms)', 'Directions based on primary sac'};
xlabel_text_CS_probab_ = {'CS probability based on [0 +200]ms'};
% figure
fig_num_ = 1;
fig_handle_(fig_num_) = figure(fig_num_);
clf(fig_handle_(fig_num_))

% Top Plot Raster
subplot(3,3,2)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_top,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis(:), y_axis(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic = CS_train_aligned(inds_event(inds_top,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 3);
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Trials')
title('Top')

% Top Plot Velocity
subplot(3,3,1)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_top,:));
SS_Firing = mean(train_data_logic) * 1000;
plot(inds_span, ESN_smooth(SS_Firing), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis right;
velocity_data = BEHAVE_eye_r_vm_filt(inds_event(inds_top,:));
velocity_data_mean = nanmean(velocity_data);
velocity_data_stdv = nanstd( velocity_data);
plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Top')

% Left Plot Raster
subplot(3,3,4)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_left,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis(:), y_axis(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic = CS_train_aligned(inds_event(inds_left,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 3);
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Trials')
title('Left')

% Left Plot Velocity
subplot(3,3,7)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_left,:));
SS_Firing = mean(train_data_logic) * 1000;
plot(inds_span, ESN_smooth(SS_Firing), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis right;
velocity_data = BEHAVE_eye_r_vm_filt(inds_event(inds_left,:));
velocity_data_mean = nanmean(velocity_data);
velocity_data_stdv = nanstd( velocity_data);
plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Left')

% Right Plot Raster
subplot(3,3,6)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_right,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis(:), y_axis(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic = CS_train_aligned(inds_event(inds_right,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 3);
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Trials')
title('Right')

% Right Plot Velocity
subplot(3,3,3)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_right,:));
SS_Firing = mean(train_data_logic) * 1000;
plot(inds_span, ESN_smooth(SS_Firing), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis right;
velocity_data = BEHAVE_eye_r_vm_filt(inds_event(inds_right,:));
velocity_data_mean = nanmean(velocity_data);
velocity_data_stdv = nanstd( velocity_data);
plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Right')

% Down Plot Raster
subplot(3,3,8)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_down,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis(:), y_axis(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic = CS_train_aligned(inds_event(inds_down,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 3);
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Trials')
title('Down')
xlabel(xlabel_text_raster_);

% Down Plot Velocity
subplot(3,3,9)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_down,:));
SS_Firing = mean(train_data_logic) * 1000;
plot(inds_span, ESN_smooth(SS_Firing), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis right;
velocity_data = BEHAVE_eye_r_vm_filt(inds_event(inds_down,:));
velocity_data_mean = nanmean(velocity_data);
velocity_data_stdv = nanstd( velocity_data);
plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Down')

% Probability
subplot(3,3,5)
hold on

plot(cosd(0:5:360)*0.10, sind(0:5:360)*0.10, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.20, sind(0:5:360)*0.20, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.30, sind(0:5:360)*0.30, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])

x_axis = [cosd(0) cosd(90) cosd(180) cosd(270) cosd(0)] .* amplitude;
y_axis = [sind(0) sind(90) sind(180) sind(270) sind(0)] .* amplitude;
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
axis equal
xlabel(xlabel_text_CS_probab_);

ESN_Beautify_Plot
hFig = fig_handle_(fig_num_);
figure_size  = [8.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

%% Plot-2 SS & CS train primSac_onset 
clearvars -except EPHYS BEHAVE fig_handle_
% inds of interest
inds_event           = EPHYS.CH_EVE.BEHAVE_inds_primSac_onset;
% CS and SS data
SS_train_aligned     = EPHYS.CH_EVE.EPHYS_SS_train_aligned;
CS_train_aligned     = EPHYS.CH_EVE.EPHYS_CS_train_aligned;
inds_span            = EPHYS.CH_EVE.inds_span_primSac_onset(1,:);
% eye velocity
BEHAVE_eye_r_vm_filt = EPHYS.CH_EVE.BEHAVE_eye_r_vm_filt;
inds_valid = BEHAVE.SACS_PRIM_DATA.validity & BEHAVE.SACS_CORR_DATA.validity;
error_prim = sqrt( (BEHAVE.SACS_PRIM_DATA.eye_r_px_finish - BEHAVE.TRIALS_DATA.cue_x).^2 + (BEHAVE.SACS_PRIM_DATA.eye_r_py_finish - BEHAVE.TRIALS_DATA.cue_y).^2 );
error_corr = sqrt( (BEHAVE.SACS_CORR_DATA.eye_r_px_finish - BEHAVE.TRIALS_DATA.end_x).^2 + (BEHAVE.SACS_CORR_DATA.eye_r_py_finish - BEHAVE.TRIALS_DATA.end_y).^2 );
inds_valid = inds_valid & (error_prim<3) & (error_corr<3);
% inds directions
inds_top   = ((BEHAVE.TRIALS_DATA.cue_x-BEHAVE.TRIALS_DATA.start_x) == 0) & ...
             ((BEHAVE.TRIALS_DATA.cue_y-BEHAVE.TRIALS_DATA.start_y) >  0) & ...
             inds_valid;
inds_left  = ((BEHAVE.TRIALS_DATA.cue_x-BEHAVE.TRIALS_DATA.start_x) <  0) & ...
             ((BEHAVE.TRIALS_DATA.cue_y-BEHAVE.TRIALS_DATA.start_y) == 0) & ...
             inds_valid;
inds_right = ((BEHAVE.TRIALS_DATA.cue_x-BEHAVE.TRIALS_DATA.start_x) >  0) & ...
             ((BEHAVE.TRIALS_DATA.cue_y-BEHAVE.TRIALS_DATA.start_y) == 0) & ...
             inds_valid;
inds_down  = ((BEHAVE.TRIALS_DATA.cue_x-BEHAVE.TRIALS_DATA.start_x) == 0) & ...
             ((BEHAVE.TRIALS_DATA.cue_y-BEHAVE.TRIALS_DATA.start_y) <  0) & ...
             inds_valid;
% CS Probab
range_inds_probability = 101:300;
train_data_logic = CS_train_aligned(inds_event(inds_top,:));
prob_top = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_top);
train_data_logic = CS_train_aligned(inds_event(inds_left,:));
prob_left = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_left);
train_data_logic = CS_train_aligned(inds_event(inds_right,:));
prob_right = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_right);
train_data_logic = CS_train_aligned(inds_event(inds_down,:));
prob_down = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_down);
amplitude = [prob_right prob_top prob_left prob_down prob_right];
% plot xlim and ylim
range_SS_Firing = [0  200];
range_vm        = [-50 600];
% xlabel text
xlabel_text_raster_ = {'Time relative to prim sac onset (ms)', 'Directions based on primary sac'};
xlabel_text_CS_probab_ = {'CS probability based on [-200 0]ms'};
% figure
fig_num_ = 2;
fig_handle_(fig_num_) = figure(fig_num_);
clf(fig_handle_(fig_num_))

% Top Plot Raster
subplot(3,3,2)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_top,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis(:), y_axis(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic = CS_train_aligned(inds_event(inds_top,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 3);
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Trials')
title('Top')

% Top Plot Velocity
subplot(3,3,1)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_top,:));
SS_Firing = mean(train_data_logic) * 1000;
plot(inds_span, ESN_smooth(SS_Firing), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis right;
velocity_data = BEHAVE_eye_r_vm_filt(inds_event(inds_top,:));
velocity_data_mean = nanmean(velocity_data);
velocity_data_stdv = nanstd( velocity_data);
plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Top')

% Left Plot Raster
subplot(3,3,4)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_left,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis(:), y_axis(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic = CS_train_aligned(inds_event(inds_left,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 3);
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Trials')
title('Left')

% Left Plot Velocity
subplot(3,3,7)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_left,:));
SS_Firing = mean(train_data_logic) * 1000;
plot(inds_span, ESN_smooth(SS_Firing), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis right;
velocity_data = BEHAVE_eye_r_vm_filt(inds_event(inds_left,:));
velocity_data_mean = nanmean(velocity_data);
velocity_data_stdv = nanstd( velocity_data);
plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Left')

% Right Plot Raster
subplot(3,3,6)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_right,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis(:), y_axis(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic = CS_train_aligned(inds_event(inds_right,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 3);
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Trials')
title('Right')

% Right Plot Velocity
subplot(3,3,3)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_right,:));
SS_Firing = mean(train_data_logic) * 1000;
plot(inds_span, ESN_smooth(SS_Firing), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis right;
velocity_data = BEHAVE_eye_r_vm_filt(inds_event(inds_right,:));
velocity_data_mean = nanmean(velocity_data);
velocity_data_stdv = nanstd( velocity_data);
plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Right')

% Down Plot Raster
subplot(3,3,8)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_down,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis(:), y_axis(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic = CS_train_aligned(inds_event(inds_down,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 3);
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Trials')
title('Down')
xlabel(xlabel_text_raster_);

% Down Plot Velocity
subplot(3,3,9)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_down,:));
SS_Firing = mean(train_data_logic) * 1000;
plot(inds_span, ESN_smooth(SS_Firing), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis right;
velocity_data = BEHAVE_eye_r_vm_filt(inds_event(inds_down,:));
velocity_data_mean = nanmean(velocity_data);
velocity_data_stdv = nanstd( velocity_data);
plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Down')

% Probability
subplot(3,3,5)
hold on

plot(cosd(0:5:360)*0.10, sind(0:5:360)*0.10, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.20, sind(0:5:360)*0.20, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.30, sind(0:5:360)*0.30, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])

x_axis = [cosd(0) cosd(90) cosd(180) cosd(270) cosd(0)] .* amplitude;
y_axis = [sind(0) sind(90) sind(180) sind(270) sind(0)] .* amplitude;
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
axis equal
xlabel(xlabel_text_CS_probab_);

ESN_Beautify_Plot
hFig = fig_handle_(fig_num_);
figure_size  = [8.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

%% Plot-3 SS & CS train primSac_offset 
clearvars -except EPHYS BEHAVE fig_handle_
% inds of interest
inds_event           = EPHYS.CH_EVE.BEHAVE_inds_primSac_offset;
% CS and SS data
SS_train_aligned     = EPHYS.CH_EVE.EPHYS_SS_train_aligned;
CS_train_aligned     = EPHYS.CH_EVE.EPHYS_CS_train_aligned;
inds_span            = EPHYS.CH_EVE.inds_span_primSac_offset(1,:);
% eye velocity
BEHAVE_eye_r_vm_filt = EPHYS.CH_EVE.BEHAVE_eye_r_vm_filt;
inds_valid = BEHAVE.SACS_PRIM_DATA.validity & BEHAVE.SACS_CORR_DATA.validity;
error_prim = sqrt( (BEHAVE.SACS_PRIM_DATA.eye_r_px_finish - BEHAVE.TRIALS_DATA.cue_x).^2 + (BEHAVE.SACS_PRIM_DATA.eye_r_py_finish - BEHAVE.TRIALS_DATA.cue_y).^2 );
error_corr = sqrt( (BEHAVE.SACS_CORR_DATA.eye_r_px_finish - BEHAVE.TRIALS_DATA.end_x).^2 + (BEHAVE.SACS_CORR_DATA.eye_r_py_finish - BEHAVE.TRIALS_DATA.end_y).^2 );
inds_valid = inds_valid & (error_prim<3) & (error_corr<3);
% inds directions
inds_top   = ((BEHAVE.TRIALS_DATA.end_x-BEHAVE.TRIALS_DATA.cue_x) == 0) & ...
             ((BEHAVE.TRIALS_DATA.end_y-BEHAVE.TRIALS_DATA.cue_y) >  0) & ...
             inds_valid;
inds_left  = ((BEHAVE.TRIALS_DATA.end_x-BEHAVE.TRIALS_DATA.cue_x) <  0) & ...
             ((BEHAVE.TRIALS_DATA.end_y-BEHAVE.TRIALS_DATA.cue_y) == 0) & ...
             inds_valid;
inds_right = ((BEHAVE.TRIALS_DATA.end_x-BEHAVE.TRIALS_DATA.cue_x) >  0) & ...
             ((BEHAVE.TRIALS_DATA.end_y-BEHAVE.TRIALS_DATA.cue_y) == 0) & ...
             inds_valid;
inds_down  = ((BEHAVE.TRIALS_DATA.end_x-BEHAVE.TRIALS_DATA.cue_x) == 0) & ...
             ((BEHAVE.TRIALS_DATA.end_y-BEHAVE.TRIALS_DATA.cue_y) <  0) & ...
             inds_valid;
% CS Probab
range_inds_probability = 101:300;
train_data_logic = CS_train_aligned(inds_event(inds_top,:));
prob_top = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_top);
train_data_logic = CS_train_aligned(inds_event(inds_left,:));
prob_left = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_left);
train_data_logic = CS_train_aligned(inds_event(inds_right,:));
prob_right = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_right);
train_data_logic = CS_train_aligned(inds_event(inds_down,:));
prob_down = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_down);
amplitude = [prob_right prob_top prob_left prob_down prob_right];
% plot xlim and ylim
range_SS_Firing = [0 200];
range_vm        = [-50 600];
% xlabel text
xlabel_text_raster_ = {'Time relative to prim sac offset (ms)', 'Directions based on secondary sac'};
xlabel_text_CS_probab_ = {'CS probability based on [0 +200]ms'};
% figure
fig_num_ = 3;
fig_handle_(fig_num_) = figure(fig_num_);
clf(fig_handle_(fig_num_))

% Top Plot Raster
subplot(3,3,2)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_top,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis(:), y_axis(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic = CS_train_aligned(inds_event(inds_top,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 3);
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Trials')
title('Top')

% Top Plot Velocity
subplot(3,3,1)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_top,:));
SS_Firing = mean(train_data_logic) * 1000;
plot(inds_span, ESN_smooth(SS_Firing), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis right;
velocity_data = BEHAVE_eye_r_vm_filt(inds_event(inds_top,:));
velocity_data_mean = nanmean(velocity_data);
velocity_data_stdv = nanstd( velocity_data);
plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Top')

% Left Plot Raster
subplot(3,3,4)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_left,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis(:), y_axis(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic = CS_train_aligned(inds_event(inds_left,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 3);
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Trials')
title('Left')

% Left Plot Velocity
subplot(3,3,7)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_left,:));
SS_Firing = mean(train_data_logic) * 1000;
plot(inds_span, ESN_smooth(SS_Firing), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis right;
velocity_data = BEHAVE_eye_r_vm_filt(inds_event(inds_left,:));
velocity_data_mean = nanmean(velocity_data);
velocity_data_stdv = nanstd( velocity_data);
plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Left')

% Right Plot Raster
subplot(3,3,6)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_right,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis(:), y_axis(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic = CS_train_aligned(inds_event(inds_right,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 3);
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Trials')
title('Right')

% Right Plot Velocity
subplot(3,3,3)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_right,:));
SS_Firing = mean(train_data_logic) * 1000;
plot(inds_span, ESN_smooth(SS_Firing), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis right;
velocity_data = BEHAVE_eye_r_vm_filt(inds_event(inds_right,:));
velocity_data_mean = nanmean(velocity_data);
velocity_data_stdv = nanstd( velocity_data);
plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Right')

% Down Plot Raster
subplot(3,3,8)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_down,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis(:), y_axis(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic = CS_train_aligned(inds_event(inds_down,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 3);
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Trials')
title('Down')
xlabel(xlabel_text_raster_);

% Down Plot Velocity
subplot(3,3,9)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_down,:));
SS_Firing = mean(train_data_logic) * 1000;
plot(inds_span, ESN_smooth(SS_Firing), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis right;
velocity_data = BEHAVE_eye_r_vm_filt(inds_event(inds_down,:));
velocity_data_mean = nanmean(velocity_data);
velocity_data_stdv = nanstd( velocity_data);
plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Down')

% Probability
subplot(3,3,5)
hold on

plot(cosd(0:5:360)*0.10, sind(0:5:360)*0.10, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.20, sind(0:5:360)*0.20, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.30, sind(0:5:360)*0.30, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])

x_axis = [cosd(0) cosd(90) cosd(180) cosd(270) cosd(0)] .* amplitude;
y_axis = [sind(0) sind(90) sind(180) sind(270) sind(0)] .* amplitude;
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
axis equal
xlabel(xlabel_text_CS_probab_);

ESN_Beautify_Plot
hFig = fig_handle_(fig_num_);
figure_size  = [8.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

%% Plot-4 SS & CS train corrSac_onset 
clearvars -except EPHYS BEHAVE fig_handle_
% inds of interest
inds_event           = EPHYS.CH_EVE.BEHAVE_inds_corrSac_onset;
% CS and SS data
SS_train_aligned     = EPHYS.CH_EVE.EPHYS_SS_train_aligned;
CS_train_aligned     = EPHYS.CH_EVE.EPHYS_CS_train_aligned;
inds_span            = EPHYS.CH_EVE.inds_span_corrSac_onset(1,:);
% eye velocity
BEHAVE_eye_r_vm_filt = EPHYS.CH_EVE.BEHAVE_eye_r_vm_filt;
inds_valid = BEHAVE.SACS_PRIM_DATA.validity & BEHAVE.SACS_CORR_DATA.validity;
error_prim = sqrt( (BEHAVE.SACS_PRIM_DATA.eye_r_px_finish - BEHAVE.TRIALS_DATA.cue_x).^2 + (BEHAVE.SACS_PRIM_DATA.eye_r_py_finish - BEHAVE.TRIALS_DATA.cue_y).^2 );
error_corr = sqrt( (BEHAVE.SACS_CORR_DATA.eye_r_px_finish - BEHAVE.TRIALS_DATA.end_x).^2 + (BEHAVE.SACS_CORR_DATA.eye_r_py_finish - BEHAVE.TRIALS_DATA.end_y).^2 );
inds_valid = inds_valid & (error_prim<3) & (error_corr<3);
% inds directions
inds_top   = ((BEHAVE.TRIALS_DATA.end_x-BEHAVE.TRIALS_DATA.cue_x) == 0) & ...
             ((BEHAVE.TRIALS_DATA.end_y-BEHAVE.TRIALS_DATA.cue_y) >  0) & ...
             inds_valid;
inds_left  = ((BEHAVE.TRIALS_DATA.end_x-BEHAVE.TRIALS_DATA.cue_x) <  0) & ...
             ((BEHAVE.TRIALS_DATA.end_y-BEHAVE.TRIALS_DATA.cue_y) == 0) & ...
             inds_valid;
inds_right = ((BEHAVE.TRIALS_DATA.end_x-BEHAVE.TRIALS_DATA.cue_x) >  0) & ...
             ((BEHAVE.TRIALS_DATA.end_y-BEHAVE.TRIALS_DATA.cue_y) == 0) & ...
             inds_valid;
inds_down  = ((BEHAVE.TRIALS_DATA.end_x-BEHAVE.TRIALS_DATA.cue_x) == 0) & ...
             ((BEHAVE.TRIALS_DATA.end_y-BEHAVE.TRIALS_DATA.cue_y) <  0) & ...
             inds_valid;
% CS Probab
range_inds_probability = 101:300;
train_data_logic = CS_train_aligned(inds_event(inds_top,:));
prob_top = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_top);
train_data_logic = CS_train_aligned(inds_event(inds_left,:));
prob_left = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_left);
train_data_logic = CS_train_aligned(inds_event(inds_right,:));
prob_right = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_right);
train_data_logic = CS_train_aligned(inds_event(inds_down,:));
prob_down = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_down);
amplitude = [prob_right prob_top prob_left prob_down prob_right];
% plot xlim and ylim
range_SS_Firing = [0  200];
range_vm        = [-50 600];
% xlabel text
xlabel_text_raster_ = {'Time relative to corr sac onset (ms)', 'Directions based on secondary sac'};
xlabel_text_CS_probab_ = {'CS probability based on [-200 0]ms'};
% figure
fig_num_ = 4;
fig_handle_(fig_num_) = figure(fig_num_);
clf(fig_handle_(fig_num_))

% Top Plot Raster
subplot(3,3,2)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_top,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis(:), y_axis(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic = CS_train_aligned(inds_event(inds_top,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 3);
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Trials')
title('Top')

% Top Plot Velocity
subplot(3,3,1)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_top,:));
SS_Firing = mean(train_data_logic) * 1000;
plot(inds_span, ESN_smooth(SS_Firing), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis right;
velocity_data = BEHAVE_eye_r_vm_filt(inds_event(inds_top,:));
velocity_data_mean = nanmean(velocity_data);
velocity_data_stdv = nanstd( velocity_data);
plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Top')

% Left Plot Raster
subplot(3,3,4)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_left,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis(:), y_axis(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic = CS_train_aligned(inds_event(inds_left,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 3);
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Trials')
title('Left')

% Left Plot Velocity
subplot(3,3,7)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_left,:));
SS_Firing = mean(train_data_logic) * 1000;
plot(inds_span, ESN_smooth(SS_Firing), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis right;
velocity_data = BEHAVE_eye_r_vm_filt(inds_event(inds_left,:));
velocity_data_mean = nanmean(velocity_data);
velocity_data_stdv = nanstd( velocity_data);
plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Left')

% Right Plot Raster
subplot(3,3,6)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_right,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis(:), y_axis(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic = CS_train_aligned(inds_event(inds_right,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 3);
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Trials')
title('Right')

% Right Plot Velocity
subplot(3,3,3)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_right,:));
SS_Firing = mean(train_data_logic) * 1000;
plot(inds_span, ESN_smooth(SS_Firing), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis right;
velocity_data = BEHAVE_eye_r_vm_filt(inds_event(inds_right,:));
velocity_data_mean = nanmean(velocity_data);
velocity_data_stdv = nanstd( velocity_data);
plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Right')

% Down Plot Raster
subplot(3,3,8)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_down,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 0.5);
plot(x_axis(:), y_axis(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic = CS_train_aligned(inds_event(inds_down,:));
[x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, inds_span, 3);
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic,1)+3)])
ylabel('Trials')
title('Down')
xlabel(xlabel_text_raster_);

% Down Plot Velocity
subplot(3,3,9)
hold on
train_data_logic = SS_train_aligned(inds_event(inds_down,:));
SS_Firing = mean(train_data_logic) * 1000;
plot(inds_span, ESN_smooth(SS_Firing), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis right;
velocity_data = BEHAVE_eye_r_vm_filt(inds_event(inds_down,:));
velocity_data_mean = nanmean(velocity_data);
velocity_data_stdv = nanstd( velocity_data);
plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Down')

% Probability
subplot(3,3,5)
hold on

plot(cosd(0:5:360)*0.10, sind(0:5:360)*0.10, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.20, sind(0:5:360)*0.20, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.30, sind(0:5:360)*0.30, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])

x_axis = [cosd(0) cosd(90) cosd(180) cosd(270) cosd(0)] .* amplitude;
y_axis = [sind(0) sind(90) sind(180) sind(270) sind(0)] .* amplitude;
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
axis equal
xlabel(xlabel_text_CS_probab_);

ESN_Beautify_Plot
hFig = fig_handle_(fig_num_);
figure_size  = [8.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');
%%
%}

%% Plot-5 Neural Properties
CH_sorted_ = concatenate_dataset(EPHYS_, 'CH_sorted', @horzcat);
EPHYS.CH_sorted.SS_data = concatenate_dataset(CH_sorted_.SS_data, [], @vertcat);
EPHYS.CH_sorted.CS_data = concatenate_dataset(CH_sorted_.CS_data, [], @vertcat);
EPHYS.CH_sorted.bin_size_time = CH_sorted_.bin_size_time(1);
EPHYS.CH_sorted.inds_span = CH_sorted_.inds_span(:,1);

if isfield(EPHYS.CH_sorted.SS_data, 'SS_waveform_hipass')
    EPHYS.CH_sorted.SS_data.SS_waveform = EPHYS.CH_sorted.SS_data.SS_waveform_hipass;
    EPHYS.CH_sorted.CS_data.CS_waveform = EPHYS.CH_sorted.CS_data.CS_waveform_hipass;
end

% figure
fig_num_ = 5;
fig_handle_(fig_num_) = figure(fig_num_);
clf(fig_handle_(fig_num_))

% subplot(2,2,1) Waveform
plot_handle_(1) = subplot(2,2,1);
hold on

plot((1:180)/30-2, mean(EPHYS.CH_sorted.SS_data.SS_waveform)+std(EPHYS.CH_sorted.SS_data.SS_waveform), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot((1:180)/30-2, mean(EPHYS.CH_sorted.SS_data.SS_waveform)-std(EPHYS.CH_sorted.SS_data.SS_waveform), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot((1:180)/30-2, mean(EPHYS.CH_sorted.SS_data.SS_waveform), ...
    'LineWidth', 2, 'Color', [0.1 0.1 0.9])

plot((1:180)/30-2, mean(EPHYS.CH_sorted.CS_data.CS_waveform)+std(EPHYS.CH_sorted.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot((1:180)/30-2, mean(EPHYS.CH_sorted.CS_data.CS_waveform)-std(EPHYS.CH_sorted.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot((1:180)/30-2, mean(EPHYS.CH_sorted.CS_data.CS_waveform), ...
    'LineWidth', 2, 'Color', [0.9 0.1 0.1])

xlabel('Time (ms)')
ylabel('Voltage (uv)')

% subplot(2,2,2) Probablities
plot_handle_(2) = subplot(2,2,2);
hold on

x_axis_ = linspace((-50+(EPHYS.CH_sorted.bin_size_time*1000)), 50, length(EPHYS.CH_sorted.inds_span))';

prob_value_ = mean(EPHYS.CH_sorted.SS_data.SSxSS_AUTO);
prob_value_(round(length(EPHYS.CH_sorted.inds_span)/2)) = 0;
y_axis_mean_ = prob_value_;
y_axis_stdv_ = sqrt(size(EPHYS.CH_sorted.SS_data.SSxSS_AUTO, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(EPHYS.CH_sorted.SS_data.SSxSS_AUTO, 1);
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', [0.1 0.1 0.9])

prob_value_ = mean(EPHYS.CH_sorted.SS_data.CSxSS_AUTO);
y_axis_mean_ = prob_value_;
y_axis_stdv_ = sqrt(size(EPHYS.CH_sorted.SS_data.CSxSS_AUTO, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(EPHYS.CH_sorted.SS_data.CSxSS_AUTO, 1);
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', [0.9 0.1 0.1])

ylim([0 inf])
xlabel('Time (ms)')
ylabel('Probability')

% subplot(2,2,3) ISI SS
plot_handle_(3) = subplot(2,2,3);
hold on

edges_SS = (0 : 0.002 : 0.050) *1000;
ISI_SS = diff(EPHYS.CH_sorted.SS_data.SS_time) * 1000;
histogram(ISI_SS,edges_SS, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', [0.1 0.1 0.9]);
set(plot_handle_(3), 'XTick', [0 0.025 0.050]*1000)
xlabel('Time (ms)')
ylabel('Probability')

% subplot(2,2,4) ISI CS
plot_handle_(4) = subplot(2,2,4);
hold on
edges_CS = 0 : 0.200 : 5.000;
ISI_CS = diff(EPHYS.CH_sorted.CS_data.CS_time);
histogram(ISI_CS,edges_CS, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', [0.9 0.1 0.1]);
set(plot_handle_(4), 'XTick', [0 2.5 5.0])
xlabel('Time (s)')
ylabel('Probability')

ESN_Beautify_Plot
hFig = fig_handle_(fig_num_);
figure_size  = [8.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

%% Save Fig
EPHYS.CH_sorted_file_name = EPHYS_(1).CH_sorted_file_name;
EPHYS.CH_sorted_file_path = EPHYS_(1).CH_sorted_file_path;
[~, file_name, ~]  = fileparts(EPHYS.CH_sorted_file_name);
if num_data_set > 1
    file_name = [file_name '_combine_' num2str(num_data_set)];
end

fprintf(['Saving plots', ' ...'])
saveas(fig_handle_(1),[EPHYS.CH_sorted_file_path file_name '_modulation_cue_present'], 'pdf');
saveas(fig_handle_(1),[EPHYS.CH_sorted_file_path file_name '_modulation_cue_present'], 'png');
saveas(fig_handle_(2),[EPHYS.CH_sorted_file_path file_name '_modulation_primSac_onset'], 'pdf');
saveas(fig_handle_(2),[EPHYS.CH_sorted_file_path file_name '_modulation_primSac_onset'], 'png');
saveas(fig_handle_(3),[EPHYS.CH_sorted_file_path file_name '_modulation_primSac_offset'], 'pdf');
saveas(fig_handle_(3),[EPHYS.CH_sorted_file_path file_name '_modulation_primSac_offset'], 'png');
saveas(fig_handle_(4),[EPHYS.CH_sorted_file_path file_name '_modulation_corrSac_onset'], 'pdf');
saveas(fig_handle_(4),[EPHYS.CH_sorted_file_path file_name '_modulation_corrSac_onset'], 'png');
saveas(fig_handle_(5),[EPHYS.CH_sorted_file_path file_name '_properties'], 'pdf');
saveas(fig_handle_(5),[EPHYS.CH_sorted_file_path file_name '_properties'], 'png');
close(fig_handle_(1))
close(fig_handle_(2))
close(fig_handle_(3))
close(fig_handle_(4))
close(fig_handle_(5))
fprintf(' --> Completed. \n')

end
%% function ESN_raster_plot_axes
function [x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, x_axis_values, line_half_len)
if nargin < 2
    x_axis_values = 1 : 1 : size(train_data_logic, 2);
    line_half_len = 0.5;
end
if nargin < 3
    line_half_len = 0.5;
end
train_data_logic = train_data_logic > 0.1;
train_data_row_number = nan(size(train_data_logic));
for counter_row = 1 : size(train_data_logic, 1)
    train_data_row_number(counter_row, train_data_logic(counter_row,:)) = counter_row;
end
train_data_col_number = repmat(x_axis_values(:)', size(train_data_logic,1), 1);
x_axis = [train_data_col_number(:)'; train_data_col_number(:)'; nan(length(train_data_col_number(:)), 1)'];
y_axis = [(train_data_row_number(:)-line_half_len)'; (train_data_row_number(:)+line_half_len)'; nan(length(train_data_row_number(:)), 1)'];
x_axis = x_axis(:);
y_axis = y_axis(:);
end
%% function ESN_smooth
function smooth_data_ = ESN_smooth(data_)
% method = 'moving';  % Moving average. A lowpass filter with filter coefficients equal to the reciprocal of the span.
% method = 'lowess';  % Local regression using weighted linear least squares and a 1st degree polynomial model.
% method = 'loess';   % Local regression using weighted linear least squares and a 2nd degree polynomial model.
% method = 'sgolay';  % Savitzky-Golay filter. A generalized moving average with filter coefficients determined by an unweighted linear least-squares regression and a polynomial model of specified degree (default is 2). The method can accept nonuniform predictor data.
% method = 'rlowess'; % A robust version of 'lowess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
% method = 'rloess';  % A robust version of 'loess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
% smooth_data_ = smooth(data_, method);
smooth_data_ = smooth(data_, 21, 'sgolay', 2);
end
%% function build_EPHYS_BEHAVE_single_dataset
function [EPHYS, BEHAVE] = build_EPHYS_BEHAVE_single_dataset
%% load EPHYS EVENT DATA
[file_name,file_path] = uigetfile([pwd filesep '*_aligned.mat'], 'Select EVENT DATA file');
fprintf(['Loading ', file_name, ' ... ']);
EPHYS.CH_EVE = load([file_path filesep file_name]);
EPHYS.CH_EVE.EPHYS_time = EPHYS.CH_EVE.EPHYS_time(:);
EPHYS.CH_EVE.BEHAVE_time = EPHYS.CH_EVE.BEHAVE_time(:);
fprintf(' --> Completed. \n')

%% load BEHAVE DATA
% [file_name,file_path] = uigetfile([file_path filesep '*_REDUCED.mat'], 'Select _REDUCED file');
[file_name,file_path] = uigetfile([file_path filesep '*_ANALYZED.mat'], 'Select _ANALYZED file');
fprintf(['Loading ', file_name, ' ... ']);
BEHAVE = load([file_path filesep file_name]);
fprintf(' --> Completed. \n')

%% load EPHYS sorted DATA
[file_name,file_path] = uigetfile([file_path filesep '*_sorted*.mat'], 'Select _sorted* file');
fprintf(['Loading ', file_name, ' ... ']);
EPHYS.CH_sorted = load([file_path filesep file_name], 'CS_data', 'SS_data');
EPHYS.CH_sorted_file_name = file_name;
EPHYS.CH_sorted_file_path = file_path;
fprintf(' --> Completed. \n')

%% build SSxSS AUTO PROBABILITY
fprintf(['Building SSxSS_AUTO & CSxSS_AUTO PROBABILITY ' ' ...'])

bin_size_time = 1e-3; % seconds
span_window_size = (1 / bin_size_time) * (100 / 1000);
span_window_size_half = round(span_window_size / 2);
inds_span = ((-span_window_size_half+1) : 1 : (span_window_size_half))';

ch_time_min = min([EPHYS.CH_sorted.SS_data.SS_time(1) EPHYS.CH_sorted.CS_data.CS_time(1)]);
ch_time_min = max([(ch_time_min-2.0) 0]);
ch_time_max = max([EPHYS.CH_sorted.SS_data.SS_time(end) EPHYS.CH_sorted.CS_data.CS_time(end)]) + 2.0;

CH__.SS_data.SS_time =  EPHYS.CH_sorted.SS_data.SS_time - ch_time_min;
CH__.CS_data.CS_time =  EPHYS.CH_sorted.CS_data.CS_time - ch_time_min;
CH__.SS_data.SS_time = ESN_Round(CH__.SS_data.SS_time, bin_size_time);
CH__.CS_data.CS_time = ESN_Round(CH__.CS_data.CS_time, bin_size_time);
CH__.SS_data.SS_ind  = round(CH__.SS_data.SS_time .* (1/bin_size_time));
CH__.SS_data.SS_ind( CH__.SS_data.SS_ind < 1 ) = 1;
CH__.CS_data.CS_ind  = round(CH__.CS_data.CS_time .* (1/bin_size_time));
CH__.CS_data.CS_ind( CH__.CS_data.CS_ind < 1 ) = 1;
CH__.CH_data.ch_time_reconstruct = ch_time_min : bin_size_time : ch_time_max;

% SSxSS_AUTO
CH__.SS_data.SS_inds_reconstruct = repmat( CH__.SS_data.SS_ind(:), 1, length(inds_span)) + repmat(inds_span(:)', length(CH__.SS_data.SS_ind), 1);
CH__.SS_data.SS_inds_reconstruct( CH__.SS_data.SS_inds_reconstruct < 1 ) = 1;
CH__.SS_data.SS_inds_reconstruct( CH__.SS_data.SS_inds_reconstruct > length( CH__.CH_data.ch_time_reconstruct ) ) = length( CH__.CH_data.ch_time_reconstruct );

CH__.SS_data.SS_event_trace = false( size(CH__.CH_data.ch_time_reconstruct) );
CH__.SS_data.SS_event_trace( CH__.SS_data.SS_ind ) = true ;
CH__.SS_data.SS_event_trace( 1   ) = false;
CH__.SS_data.SS_event_trace( end ) = false;

CH__.SS_data.SS_event_reconstruct = CH__.SS_data.SS_event_trace( CH__.SS_data.SS_inds_reconstruct );
EPHYS.CH_sorted.SS_data.SSxSS_AUTO = CH__.SS_data.SS_event_reconstruct;

% CSxSS_WITHIN
CH__.CS_data.CS_inds_reconstruct = repmat( CH__.CS_data.CS_ind(:), 1, length(inds_span)) + repmat(inds_span(:)', length(CH__.CS_data.CS_ind), 1);
CH__.CS_data.CS_inds_reconstruct( CH__.CS_data.CS_inds_reconstruct < 1 ) = 1;
CH__.CS_data.CS_inds_reconstruct( CH__.CS_data.CS_inds_reconstruct > length( CH__.CH_data.ch_time_reconstruct ) ) = length( CH__.CH_data.ch_time_reconstruct );

CH__.SS_data.SS_event_trace = false( size(CH__.CH_data.ch_time_reconstruct) );
CH__.SS_data.SS_event_trace( CH__.SS_data.SS_ind ) = true ;
CH__.SS_data.SS_event_trace( 1   ) = false;
CH__.SS_data.SS_event_trace( end ) = false;

CH__.SS_data.SS_event_reconstruct = CH__.SS_data.SS_event_trace( CH__.CS_data.CS_inds_reconstruct );
EPHYS.CH_sorted.SS_data.CSxSS_AUTO = CH__.SS_data.SS_event_reconstruct;

EPHYS.CH_sorted.inds_span = inds_span;
EPHYS.CH_sorted.bin_size_time = bin_size_time;
fprintf(' --> Completed. \n')

%% SS & CS train_aligned
clearvars -except EPHYS BEHAVE
fprintf(['Building CS & SS train_aligned', ' ... ']);
EPHYS_time_aligned     = EPHYS.CH_EVE.EPHYS_time_aligned;
length_time_ = length(EPHYS_time_aligned);
CS_time = EPHYS.CH_sorted.CS_data.CS_time;
CS_time(end+1) = max([EPHYS_time_aligned(end), CS_time(end)])+1;
SS_time = EPHYS.CH_sorted.SS_data.SS_time;
SS_time(end+1) = max([EPHYS_time_aligned(end), SS_time(end)])+1;
EPHYS_CS_train_aligned = false(size(EPHYS_time_aligned));
EPHYS_SS_train_aligned = false(size(EPHYS_time_aligned));
counter_CS = find(CS_time > EPHYS_time_aligned(1), 1, 'first');
counter_SS = find(SS_time > EPHYS_time_aligned(1), 1, 'first');
for counter_time_point = 1 : length_time_
    time_ponit_ = EPHYS_time_aligned(counter_time_point);
    if time_ponit_>CS_time(counter_CS)
        EPHYS_CS_train_aligned(counter_time_point) = true;
        counter_CS = counter_CS + 1;
    end
    if time_ponit_>SS_time(counter_SS)
        EPHYS_SS_train_aligned(counter_time_point) = true;
        counter_SS = counter_SS + 1;
    end
end
EPHYS.CH_EVE.EPHYS_CS_train_aligned = EPHYS_CS_train_aligned;
EPHYS.CH_EVE.EPHYS_SS_train_aligned = EPHYS_SS_train_aligned;

fprintf(' --> Completed. \n')

%% Build Raster Plot Data
clearvars -except EPHYS BEHAVE
fprintf(['Building Raster Plot Data', ' ... ']);
num_trials = length(BEHAVE.TRIALS_DATA.time_end);
BEHAVE_time_aligned     = EPHYS.CH_EVE.BEHAVE_time_aligned;
length_time_ = length(BEHAVE_time_aligned);
BEHAVE_time_cue_present    = nan(num_trials, 1);
BEHAVE_time_primSac_onset  = nan(num_trials, 1);
BEHAVE_time_primSac_offset = nan(num_trials, 1);
BEHAVE_time_corrSac_onset  = nan(num_trials, 1);
for counter_trial = 1 : 1 : num_trials
    BEHAVE_time_cue_present(counter_trial) = BEHAVE.TRIALS_DATA.time_state_cue_present{1,counter_trial}(end);
    ind_primSac_onset_  = (BEHAVE.SACS_PRIM_DATA.ind_start(1,counter_trial)  == BEHAVE.SACS_PRIM_DATA.inds(:,counter_trial));
    if sum(ind_primSac_onset_) ~= 1
        ind_primSac_onset_ = 1;
    end
    BEHAVE_time_primSac_onset(counter_trial)  = BEHAVE.SACS_PRIM_DATA.time(ind_primSac_onset_, counter_trial);
    ind_primSac_offset_ = (BEHAVE.SACS_PRIM_DATA.ind_finish(1,counter_trial) == BEHAVE.SACS_PRIM_DATA.inds(:,counter_trial));
    if sum(ind_primSac_offset_) ~= 1
        ind_primSac_offset_ = 150;
    end
    BEHAVE_time_primSac_offset(counter_trial) = BEHAVE.SACS_PRIM_DATA.time(ind_primSac_offset_, counter_trial);
    ind_corrSac_onset_  = (BEHAVE.SACS_CORR_DATA.ind_start(1,counter_trial)  == BEHAVE.SACS_CORR_DATA.inds(:,counter_trial));
    if sum(ind_corrSac_onset_) ~= 1
        ind_corrSac_onset_ = 1;
    end
    BEHAVE_time_corrSac_onset(counter_trial)  = BEHAVE.SACS_CORR_DATA.time(ind_corrSac_onset_, counter_trial);
end
BEHAVE_time_cue_present(end+1)    = max([BEHAVE_time_aligned(end), BEHAVE_time_cue_present(end)])+1;
BEHAVE_time_primSac_onset(end+1)  = max([BEHAVE_time_aligned(end), BEHAVE_time_primSac_onset(end)])+1;
BEHAVE_time_primSac_offset(end+1) = max([BEHAVE_time_aligned(end), BEHAVE_time_primSac_offset(end)])+1;
BEHAVE_time_corrSac_onset(end+1)  = max([BEHAVE_time_aligned(end), BEHAVE_time_corrSac_onset(end)])+1;
BEHAVE_ind_cue_present    = nan(num_trials, 1);
BEHAVE_ind_primSac_onset  = nan(num_trials, 1);
BEHAVE_ind_primSac_offset = nan(num_trials, 1);
BEHAVE_ind_corrSac_onset  = nan(num_trials, 1);
counter_cue_present    = find(BEHAVE_time_cue_present    > BEHAVE_time_aligned(1), 1, 'first');
counter_primSac_onset  = find(BEHAVE_time_primSac_onset  > BEHAVE_time_aligned(1), 1, 'first');
counter_primSac_offset = find(BEHAVE_time_primSac_offset > BEHAVE_time_aligned(1), 1, 'first');
counter_corrSac_onset  = find(BEHAVE_time_corrSac_onset  > BEHAVE_time_aligned(1), 1, 'first');
counter_trial_cue_present    = 1;
counter_trial_primSac_onset  = 1;
counter_trial_primSac_offset = 1;
counter_trial_corrSac_onset  = 1;
for counter_time_point = 1 : length_time_
    time_ponit_ = BEHAVE_time_aligned(counter_time_point);
    if time_ponit_>BEHAVE_time_cue_present(counter_cue_present)
        BEHAVE_ind_cue_present(counter_trial_cue_present) = counter_time_point;
        counter_cue_present = counter_cue_present + 1;
        counter_trial_cue_present = counter_trial_cue_present + 1;
    end
    if time_ponit_>BEHAVE_time_primSac_onset(counter_primSac_onset)
        BEHAVE_ind_primSac_onset(counter_trial_primSac_onset) = counter_time_point;
        counter_primSac_onset = counter_primSac_onset + 1;
        counter_trial_primSac_onset = counter_trial_primSac_onset + 1;
    end
    if time_ponit_>BEHAVE_time_primSac_offset(counter_primSac_offset)
        BEHAVE_ind_primSac_offset(counter_trial_primSac_offset) = counter_time_point;
        counter_primSac_offset = counter_primSac_offset + 1;
        counter_trial_primSac_offset = counter_trial_primSac_offset + 1;
    end
    if time_ponit_>BEHAVE_time_corrSac_onset(counter_corrSac_onset)
        BEHAVE_ind_corrSac_onset(counter_trial_corrSac_onset) = counter_time_point;
        counter_corrSac_onset = counter_corrSac_onset + 1;
        counter_trial_corrSac_onset = counter_trial_corrSac_onset + 1;
    end
end

% Build inds
inds_span_cue_present = ((-100+1) : 1 : (300))';
BEHAVE_time_aligned     = EPHYS.CH_EVE.BEHAVE_time_aligned;
length_time_ = length(BEHAVE_time_aligned);
BEHAVE_inds_cue_present = repmat( BEHAVE_ind_cue_present(:), 1, length(inds_span_cue_present)) + repmat(inds_span_cue_present(:)', length(BEHAVE_ind_cue_present), 1);
BEHAVE_inds_cue_present( BEHAVE_inds_cue_present < 1 ) = 1;
BEHAVE_inds_cue_present( BEHAVE_inds_cue_present > length_time_ ) = length_time_;

inds_span_primSac_onset = ((-300+1) : 1 : (100))';
BEHAVE_inds_primSac_onset = repmat( BEHAVE_ind_primSac_onset(:), 1, length(inds_span_primSac_onset)) + repmat(inds_span_primSac_onset(:)', length(BEHAVE_ind_primSac_onset), 1);
BEHAVE_inds_primSac_onset( BEHAVE_inds_primSac_onset < 1 ) = 1;
BEHAVE_inds_primSac_onset( BEHAVE_inds_primSac_onset > length_time_ ) = length_time_;

inds_span_primSac_offset = ((-100+1) : 1 : (300))';
BEHAVE_inds_primSac_offset = repmat( BEHAVE_ind_primSac_offset(:), 1, length(inds_span_primSac_offset)) + repmat(inds_span_primSac_offset(:)', length(BEHAVE_ind_primSac_offset), 1);
BEHAVE_inds_primSac_offset( BEHAVE_inds_primSac_offset < 1 ) = 1;
BEHAVE_inds_primSac_offset( BEHAVE_inds_primSac_offset > length_time_ ) = length_time_;

inds_span_corrSac_onset = ((-300+1) : 1 : (100))';
BEHAVE_inds_corrSac_onset = repmat( BEHAVE_ind_corrSac_onset(:), 1, length(inds_span_corrSac_onset)) + repmat(inds_span_corrSac_onset(:)', length(BEHAVE_ind_corrSac_onset), 1);
BEHAVE_inds_corrSac_onset( BEHAVE_inds_corrSac_onset < 1 ) = 1;
BEHAVE_inds_corrSac_onset( BEHAVE_inds_corrSac_onset > length_time_ ) = length_time_;

EPHYS.CH_EVE.BEHAVE_ind_cue_present     = BEHAVE_ind_cue_present;
EPHYS.CH_EVE.BEHAVE_ind_primSac_onset   = BEHAVE_ind_primSac_onset;
EPHYS.CH_EVE.BEHAVE_ind_primSac_offset  = BEHAVE_ind_primSac_offset;
EPHYS.CH_EVE.BEHAVE_ind_corrSac_onset   = BEHAVE_ind_corrSac_onset;
EPHYS.CH_EVE.BEHAVE_inds_cue_present    = BEHAVE_inds_cue_present;
EPHYS.CH_EVE.BEHAVE_inds_primSac_onset  = BEHAVE_inds_primSac_onset;
EPHYS.CH_EVE.BEHAVE_inds_primSac_offset = BEHAVE_inds_primSac_offset;
EPHYS.CH_EVE.BEHAVE_inds_corrSac_onset  = BEHAVE_inds_corrSac_onset;
EPHYS.CH_EVE.inds_span_cue_present      = inds_span_cue_present(:)';
EPHYS.CH_EVE.inds_span_primSac_onset    = inds_span_primSac_onset(:)';
EPHYS.CH_EVE.inds_span_primSac_offset   = inds_span_primSac_offset(:)';
EPHYS.CH_EVE.inds_span_corrSac_onset    = inds_span_corrSac_onset(:)';
fprintf(' --> Completed. \n')

%% Build BEHAVE_eye_r_vm_filt
clearvars -except EPHYS BEHAVE
fprintf(['Building BEHAVE_eye_r_vm_filt', ' ... ']);
eye_r_vm_filt = cell2mat(BEHAVE.TRIALS_DATA.eye_r_vm_filt(:));
time_1K = cell2mat(BEHAVE.TRIALS_DATA.time_1K(:));
BEHAVE_time_aligned     = EPHYS.CH_EVE.BEHAVE_time_aligned;
length_time_ = length(BEHAVE_time_aligned);
BEHAVE_eye_r_vm_filt = nan(size(BEHAVE_time_aligned));
time_1K(end+1)    = max([BEHAVE_time_aligned(end), time_1K(end)])+1;
counter_time_1K = find(time_1K >= BEHAVE_time_aligned(1), 1, 'first');
for counter_time_point = 1 : length_time_
    time_ponit_ = BEHAVE_time_aligned(counter_time_point);
    if time_ponit_>=time_1K(counter_time_1K)
        BEHAVE_eye_r_vm_filt(counter_time_point) = eye_r_vm_filt(counter_time_1K);
        counter_time_1K = counter_time_1K + 1;
    end
end
EPHYS.CH_EVE.BEHAVE_eye_r_vm_filt = BEHAVE_eye_r_vm_filt;
fprintf(' --> Completed. \n')
end
%% function concatenate_dataset
function upper_field_struct = concatenate_dataset(dataset_, upper_field_name, horz_OR_vert)
if ~isempty(upper_field_name)
dataset = struct(upper_field_name,struct());
field_names_ = fieldnames(dataset_(1).(upper_field_name));
for counter_fields = 1 : 1 : length(field_names_)
    for counter_dataset = 1 : 1 : length(dataset_)
        variable_TRIALS_DATA_ALL_ = dataset_(counter_dataset).(upper_field_name).(field_names_{counter_fields});
        % the field does not exist in TRIALS_DATA
        if ~isfield(dataset.(upper_field_name), field_names_{counter_fields})
            dataset.(upper_field_name).(field_names_{counter_fields}) = [];
        end
        variable_TRIALS_DATA_ = dataset.(upper_field_name).(field_names_{counter_fields});
        variable_TRIALS_DATA_ = horz_OR_vert(variable_TRIALS_DATA_, variable_TRIALS_DATA_ALL_);
        dataset.(upper_field_name).(field_names_{counter_fields}) = variable_TRIALS_DATA_;
    end
end
upper_field_struct = dataset.(upper_field_name);
elseif isempty(upper_field_name)
dataset = struct();
field_names_ = fieldnames(dataset_);
for counter_fields = 1 : 1 : length(field_names_)
    for counter_dataset = 1 : 1 : length(dataset_)
        variable_TRIALS_DATA_ALL_ = dataset_(counter_dataset).(field_names_{counter_fields});
        % the field does not exist in TRIALS_DATA
        if ~isfield(dataset, field_names_{counter_fields})
            dataset.(field_names_{counter_fields}) = [];
        end
        variable_TRIALS_DATA_ = dataset.(field_names_{counter_fields});
        variable_TRIALS_DATA_ = horz_OR_vert(variable_TRIALS_DATA_, variable_TRIALS_DATA_ALL_);
        dataset.(field_names_{counter_fields}) = variable_TRIALS_DATA_;
    end
end
upper_field_struct = dataset;
end
end
%% function plot_rasters_data
function fig_handle_ = plot_rasters_data(raster_data, plot_data)
fig_num_               = plot_data.fig_num_;
xlabel_text_raster_    = plot_data.xlabel_text_raster_;
xlabel_text_CS_probab_ = plot_data.xlabel_text_CS_probab_;
inds_span              = plot_data.inds_span;

train_data_logic_SS_top = raster_data.train_data_logic_SS_top;
train_data_logic_CS_top = raster_data.train_data_logic_CS_top;
velocity_data_top = raster_data.velocity_data_top;
% velocity_data_mean_top = raster_data.velocity_data_mean_top;
% velocity_data_stdv_top = raster_data.velocity_data_stdv_top;

train_data_logic_SS_left = raster_data.train_data_logic_SS_left;
train_data_logic_CS_left = raster_data.train_data_logic_CS_left;
velocity_data_left = raster_data.velocity_data_left;
% velocity_data_mean_left = raster_data.velocity_data_mean_left;
% velocity_data_stdv_left = raster_data.velocity_data_stdv_left;

train_data_logic_SS_right = raster_data.train_data_logic_SS_right;
train_data_logic_CS_right = raster_data.train_data_logic_CS_right;
velocity_data_right = raster_data.velocity_data_right;
% velocity_data_mean_right = raster_data.velocity_data_mean_right;
% velocity_data_stdv_right = raster_data.velocity_data_stdv_right;

train_data_logic_SS_down = raster_data.train_data_logic_SS_down;
train_data_logic_CS_down = raster_data.train_data_logic_CS_down;
velocity_data_down = raster_data.velocity_data_down;
% velocity_data_mean_down = raster_data.velocity_data_mean_down;
% velocity_data_stdv_down = raster_data.velocity_data_stdv_down;


% % CS Probab
range_inds_probability = 101:300;
prob_top   = sum( sum(train_data_logic_CS_top(  :,range_inds_probability) ,2) > 0 ) / size(train_data_logic_CS_top, 1);
prob_left  = sum( sum(train_data_logic_CS_left( :,range_inds_probability) ,2) > 0 ) / size(train_data_logic_CS_left, 1);
prob_right = sum( sum(train_data_logic_CS_right(:,range_inds_probability) ,2) > 0 ) / size(train_data_logic_CS_right, 1);
prob_down  = sum( sum(train_data_logic_CS_down( :,range_inds_probability) ,2) > 0 ) / size(train_data_logic_CS_down, 1);
prob_amplitude = [prob_right prob_top prob_left prob_down prob_right];
% % plot xlim and ylim
range_SS_Firing = [0 200];
range_vm        = [-50 600];
% % xlabel text
% xlabel_text_raster_ = {'Time relative to cue presentation (ms)', 'Directions based on primary sac'};
% xlabel_text_CS_probab_ = {'CS probability based on [0 +200]ms'};
% % figure
% fig_num_ = 1;
fig_handle_ = figure(fig_num_);
clf(fig_handle_)

% % Top Plot Raster
subplot(3,3,2)
hold on
% train_data_logic_SS_top = SS_train_aligned(inds_event(inds_top,:));
[x_axis_SS_top, y_axis_SS_top] = ESN_raster_plot_axes(train_data_logic_SS_top, inds_span, 0.5);
plot(x_axis_SS_top(:), y_axis_SS_top(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
% train_data_logic_CS_top = CS_train_aligned(inds_event(inds_top,:));
[x_axis_CS_top, y_axis_CS_top] = ESN_raster_plot_axes(train_data_logic_CS_top, inds_span, 3);
plot(x_axis_CS_top(:), y_axis_CS_top(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic_CS_top,1)+3)])
ylabel('Trials')
title('Top')

% % Top Plot Velocity
subplot(3,3,1)
hold on
yyaxis left;
% velocity_data_top = BEHAVE_eye_r_vm_filt(inds_event(inds_top,:));
velocity_data_mean_top = nanmean(velocity_data_top);
velocity_data_stdv_top = nanstd( velocity_data_top);
plot(inds_span, velocity_data_mean_top+velocity_data_stdv_top, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean_top-velocity_data_stdv_top, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean_top,                        '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis right;
% train_data_logic_ = SS_train_aligned(inds_event(inds_top,:));
firing_SS_top = mean(train_data_logic_SS_top) * 1000;
plot(inds_span, ESN_smooth(firing_SS_top), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Top')

% % Left Plot Raster
subplot(3,3,4)
hold on
% train_data_logic_SS_left = SS_train_aligned(inds_event(inds_left,:));
[x_axis_SS_left, y_axis_SS_left] = ESN_raster_plot_axes(train_data_logic_SS_left, inds_span, 0.5);
plot(x_axis_SS_left(:), y_axis_SS_left(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
% train_data_logic_CS_left = CS_train_aligned(inds_event(inds_left,:));
[x_axis_CS_left, y_axis_CS_left] = ESN_raster_plot_axes(train_data_logic_CS_left, inds_span, 3);
plot(x_axis_CS_left(:), y_axis_CS_left(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic_CS_left,1)+3)])
ylabel('Trials')
title('Left')

% % Left Plot Velocity
subplot(3,3,7)
hold on
yyaxis left;
% velocity_data_left = BEHAVE_eye_r_vm_filt(inds_event(inds_left,:));
velocity_data_mean_left = nanmean(velocity_data_left);
velocity_data_stdv_left = nanstd( velocity_data_left);
plot(inds_span, velocity_data_mean_left+velocity_data_stdv_left, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean_left-velocity_data_stdv_left, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean_left,                         '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis right;
% train_data_logic_ = SS_train_aligned(inds_event(inds_left,:));
firing_SS_left = mean(train_data_logic_SS_left) * 1000;
plot(inds_span, ESN_smooth(firing_SS_left), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Left')

% % Right Plot Raster
subplot(3,3,6)
hold on
% train_data_logic_SS_right = SS_train_aligned(inds_event(inds_right,:));
[x_axis_SS_right, y_axis_SS_right] = ESN_raster_plot_axes(train_data_logic_SS_right, inds_span, 0.5);
plot(x_axis_SS_right(:), y_axis_SS_right(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
% train_data_logic_CS_right = CS_train_aligned(inds_event(inds_right,:));
[x_axis_CS_right, y_axis_CS_right] = ESN_raster_plot_axes(train_data_logic_CS_right, inds_span, 3);
plot(x_axis_CS_right(:), y_axis_CS_right(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic_CS_right,1)+3)])
ylabel('Trials')
title('Right')

% % Right Plot Velocity
subplot(3,3,3)
hold on
yyaxis left;
% velocity_data_right = BEHAVE_eye_r_vm_filt(inds_event(inds_right,:));
velocity_data_mean_right = nanmean(velocity_data_right);
velocity_data_stdv_right = nanstd( velocity_data_right);
plot(inds_span, velocity_data_mean_right+velocity_data_stdv_right, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean_right-velocity_data_stdv_right, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean_right,                          '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis right;
% train_data_logic_ = SS_train_aligned(inds_event(inds_right,:));
firing_SS_right = mean(train_data_logic_SS_right) * 1000;
plot(inds_span, ESN_smooth(firing_SS_right), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Right')

% % Down Plot Raster
subplot(3,3,8)
hold on
% train_data_logic_SS_down = SS_train_aligned(inds_event(inds_down,:));
[x_axis_SS_down, y_axis_SS_down] = ESN_raster_plot_axes(train_data_logic_SS_down, inds_span, 0.5);
plot(x_axis_SS_down(:), y_axis_SS_down(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
% train_data_logic_CS_down = CS_train_aligned(inds_event(inds_down,:));
[x_axis_CS_down, y_axis_CS_down] = ESN_raster_plot_axes(train_data_logic_CS_down, inds_span, 3);
plot(x_axis_CS_down(:), y_axis_CS_down(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic_CS_down,1)+3)])
ylabel('Trials')
title('Down')
xlabel(xlabel_text_raster_);

% % Down Plot Velocity
subplot(3,3,9)
hold on
yyaxis left;
% velocity_data_down = BEHAVE_eye_r_vm_filt(inds_event(inds_down,:));
velocity_data_mean_down = nanmean(velocity_data_down);
velocity_data_stdv_down = nanstd( velocity_data_down);
plot(inds_span, velocity_data_mean_down+velocity_data_stdv_down, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean_down-velocity_data_stdv_down, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
plot(inds_span, velocity_data_mean_down,                         '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', [0.1 0.9 0.1])
yyaxis right;
% train_data_logic_ = SS_train_aligned(inds_event(inds_down,:));
firing_SS_down = mean(train_data_logic_SS_down) * 1000;
plot(inds_span, ESN_smooth(firing_SS_down), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Down')

% % Probability
subplot(3,3,5)
hold on

plot(cosd(0:5:360)*0.10, sind(0:5:360)*0.10, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.20, sind(0:5:360)*0.20, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.30, sind(0:5:360)*0.30, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])

x_axis = [cosd(0) cosd(90) cosd(180) cosd(270) cosd(0)] .* prob_amplitude;
y_axis = [sind(0) sind(90) sind(180) sind(270) sind(0)] .* prob_amplitude;
plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
axis equal
xlabel(xlabel_text_CS_probab_);

ESN_Beautify_Plot
hFig = fig_handle_;
figure_size  = [8.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');
end
%% function single_dataset_raster
function raster_data = single_dataset_raster(EPHYS, BEHAVE, inds_event, prim_OR_corr)
% inds of interest
% inds_event           = EPHYS.CH_EVE.BEHAVE_inds_cue_present;
% % CS and SS data
SS_train_aligned     = EPHYS.CH_EVE.EPHYS_SS_train_aligned;
CS_train_aligned     = EPHYS.CH_EVE.EPHYS_CS_train_aligned;
% inds_span            = EPHYS.CH_EVE.inds_span_cue_present(1,:);
% % eye velocity
BEHAVE_eye_r_vm_filt = EPHYS.CH_EVE.BEHAVE_eye_r_vm_filt;
inds_valid = BEHAVE.SACS_PRIM_DATA.validity & BEHAVE.SACS_CORR_DATA.validity;
error_prim = sqrt( (BEHAVE.SACS_PRIM_DATA.eye_r_px_finish - BEHAVE.TRIALS_DATA.cue_x).^2 + (BEHAVE.SACS_PRIM_DATA.eye_r_py_finish - BEHAVE.TRIALS_DATA.cue_y).^2 );
error_corr = sqrt( (BEHAVE.SACS_CORR_DATA.eye_r_px_finish - BEHAVE.TRIALS_DATA.end_x).^2 + (BEHAVE.SACS_CORR_DATA.eye_r_py_finish - BEHAVE.TRIALS_DATA.end_y).^2 );
inds_valid = inds_valid & (error_prim<3) & (error_corr<3);
if contains(prim_OR_corr, 'prim')
% % inds directions
inds_top   = ((BEHAVE.TRIALS_DATA.cue_x-BEHAVE.TRIALS_DATA.start_x) == 0) & ...
             ((BEHAVE.TRIALS_DATA.cue_y-BEHAVE.TRIALS_DATA.start_y) >  0) & ...
             inds_valid;
inds_left  = ((BEHAVE.TRIALS_DATA.cue_x-BEHAVE.TRIALS_DATA.start_x) <  0) & ...
             ((BEHAVE.TRIALS_DATA.cue_y-BEHAVE.TRIALS_DATA.start_y) == 0) & ...
             inds_valid;
inds_right = ((BEHAVE.TRIALS_DATA.cue_x-BEHAVE.TRIALS_DATA.start_x) >  0) & ...
             ((BEHAVE.TRIALS_DATA.cue_y-BEHAVE.TRIALS_DATA.start_y) == 0) & ...
             inds_valid;
inds_down  = ((BEHAVE.TRIALS_DATA.cue_x-BEHAVE.TRIALS_DATA.start_x) == 0) & ...
             ((BEHAVE.TRIALS_DATA.cue_y-BEHAVE.TRIALS_DATA.start_y) <  0) & ...
             inds_valid;
elseif contains(prim_OR_corr, 'corr')
% inds directions
inds_top   = ((BEHAVE.TRIALS_DATA.end_x-BEHAVE.TRIALS_DATA.cue_x) == 0) & ...
             ((BEHAVE.TRIALS_DATA.end_y-BEHAVE.TRIALS_DATA.cue_y) >  0) & ...
             inds_valid;
inds_left  = ((BEHAVE.TRIALS_DATA.end_x-BEHAVE.TRIALS_DATA.cue_x) <  0) & ...
             ((BEHAVE.TRIALS_DATA.end_y-BEHAVE.TRIALS_DATA.cue_y) == 0) & ...
             inds_valid;
inds_right = ((BEHAVE.TRIALS_DATA.end_x-BEHAVE.TRIALS_DATA.cue_x) >  0) & ...
             ((BEHAVE.TRIALS_DATA.end_y-BEHAVE.TRIALS_DATA.cue_y) == 0) & ...
             inds_valid;
inds_down  = ((BEHAVE.TRIALS_DATA.end_x-BEHAVE.TRIALS_DATA.cue_x) == 0) & ...
             ((BEHAVE.TRIALS_DATA.end_y-BEHAVE.TRIALS_DATA.cue_y) <  0) & ...
             inds_valid;
end
% % CS Probab
% range_inds_probability = 101:300;
% train_data_logic = CS_train_aligned(inds_event(inds_top,:));
% prob_top   = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_top);
% train_data_logic = CS_train_aligned(inds_event(inds_left,:));
% prob_left  = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_left);
% train_data_logic = CS_train_aligned(inds_event(inds_right,:));
% prob_right = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_right);
% train_data_logic = CS_train_aligned(inds_event(inds_down,:));
% prob_down  = sum( sum(train_data_logic(:,range_inds_probability),2) > 0 ) / sum(inds_down);
% prob_amplitude = [prob_right prob_top prob_left prob_down prob_right];
% % plot xlim and ylim
% range_SS_Firing = [0 200];
% range_vm        = [-50 600];
% xlabel text
% xlabel_text_raster_ = {'Time relative to cue presentation (ms)', 'Directions based on primary sac'};
% xlabel_text_CS_probab_ = {'CS probability based on [0 +200]ms'};
% figure
% fig_num_ = 1;
% fig_handle_(fig_num_) = figure(fig_num_);
% clf(fig_handle_(fig_num_))

% Top Plot Raster
% subplot(3,3,2)
% hold on
train_data_logic_SS_top = SS_train_aligned(inds_event(inds_top,:));
% [x_axis_SS_top, y_axis_SS_top] = ESN_raster_plot_axes(train_data_logic_SS_top, inds_span, 0.5);
% plot(x_axis_SS_top(:), y_axis_SS_top(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic_CS_top = CS_train_aligned(inds_event(inds_top,:));
% [x_axis_CS_top, y_axis_CS_top] = ESN_raster_plot_axes(train_data_logic_CS_top, inds_span, 3);
% plot(x_axis_CS_top(:), y_axis_CS_top(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
% xlim([min(inds_span)-1 max(inds_span)+1])
% ylim([(1-3) (size(train_data_logic,1)+3)])
% ylabel('Trials')
% title('Top')

% % Top Plot Velocity
% subplot(3,3,1)
% hold on
% train_data_logic_ = SS_train_aligned(inds_event(inds_top,:));
% firing_SS_top = mean(train_data_logic_) * 1000;
% plot(inds_span, ESN_smooth(firing_SS_top), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
% ylabel('SS Firing (spk/s)')
% ylim(range_SS_Firing)
% yyaxis right;
velocity_data_top = BEHAVE_eye_r_vm_filt(inds_event(inds_top,:));
% velocity_data_mean_top = nanmean(velocity_data_top);
% velocity_data_stdv_top = nanstd( velocity_data_top);
% plot(inds_span, velocity_data_mean_top+velocity_data_stdv_top, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
% plot(inds_span, velocity_data_mean_top-velocity_data_stdv_top, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
% plot(inds_span, velocity_data_mean_top,                        '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
% ylim(range_vm)
% set(gca, 'YColor', [0.1 0.9 0.1])
% yyaxis left;
% xlim([min(inds_span)-1 max(inds_span)+1])
% title('Top')

% % Left Plot Raster
% subplot(3,3,4)
% hold on
train_data_logic_SS_left = SS_train_aligned(inds_event(inds_left,:));
% [x_axis_SS_left, y_axis_SS_left] = ESN_raster_plot_axes(train_data_logic_SS_left, inds_span, 0.5);
% plot(x_axis_SS_left(:), y_axis_SS_left(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic_CS_left = CS_train_aligned(inds_event(inds_left,:));
% [x_axis_CS_left, y_axis_CS_left] = ESN_raster_plot_axes(train_data_logic_CS_left, inds_span, 3);
% plot(x_axis_CS_left(:), y_axis_CS_left(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
% xlim([min(inds_span)-1 max(inds_span)+1])
% ylim([(1-3) (size(train_data_logic,1)+3)])
% ylabel('Trials')
% title('Left')

% % Left Plot Velocity
% subplot(3,3,7)
% hold on
% train_data_logic_ = SS_train_aligned(inds_event(inds_left,:));
% firing_SS_left = mean(train_data_logic_) * 1000;
% plot(inds_span, ESN_smooth(firing_SS_left), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
% ylabel('SS Firing (spk/s)')
% ylim(range_SS_Firing)
% yyaxis right;
velocity_data_left = BEHAVE_eye_r_vm_filt(inds_event(inds_left,:));
% velocity_data_mean_left = nanmean(velocity_data_left);
% velocity_data_stdv_left = nanstd( velocity_data_left);
% plot(inds_span, velocity_data_mean_left+velocity_data_stdv_left, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
% plot(inds_span, velocity_data_mean_left-velocity_data_stdv_left, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
% plot(inds_span, velocity_data_mean_left,                         '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
% ylim(range_vm)
% set(gca, 'YColor', [0.1 0.9 0.1])
% yyaxis left;
% xlim([min(inds_span)-1 max(inds_span)+1])
% title('Left')

% % Right Plot Raster
% subplot(3,3,6)
% hold on
train_data_logic_SS_right = SS_train_aligned(inds_event(inds_right,:));
% [x_axis_SS_right, y_axis_SS_right] = ESN_raster_plot_axes(train_data_logic_SS_right, inds_span, 0.5);
% plot(x_axis_SS_right(:), y_axis_SS_right(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic_CS_right = CS_train_aligned(inds_event(inds_right,:));
% [x_axis_CS_right, y_axis_CS_right] = ESN_raster_plot_axes(train_data_logic_CS_right, inds_span, 3);
% plot(x_axis_CS_right(:), y_axis_CS_right(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
% xlim([min(inds_span)-1 max(inds_span)+1])
% ylim([(1-3) (size(train_data_logic,1)+3)])
% ylabel('Trials')
% title('Right')

% % Right Plot Velocity
% subplot(3,3,3)
% hold on
% train_data_logic_ = SS_train_aligned(inds_event(inds_right,:));
% firing_SS_right = mean(train_data_logic_) * 1000;
% plot(inds_span, ESN_smooth(firing_SS_right), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
% ylabel('SS Firing (spk/s)')
% ylim(range_SS_Firing)
% yyaxis right;
velocity_data_right = BEHAVE_eye_r_vm_filt(inds_event(inds_right,:));
% velocity_data_mean_right = nanmean(velocity_data_right);
% velocity_data_stdv_right = nanstd( velocity_data_right);
% plot(inds_span, velocity_data_mean_right+velocity_data_stdv_right, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
% plot(inds_span, velocity_data_mean_right-velocity_data_stdv_right, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
% plot(inds_span, velocity_data_mean_right,                          '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
% ylim(range_vm)
% set(gca, 'YColor', [0.1 0.9 0.1])
% yyaxis left;
% xlim([min(inds_span)-1 max(inds_span)+1])
% title('Right')

% % Down Plot Raster
% subplot(3,3,8)
% hold on
train_data_logic_SS_down = SS_train_aligned(inds_event(inds_down,:));
% [x_axis_SS_down, y_axis_SS_down] = ESN_raster_plot_axes(train_data_logic_SS_down, inds_span, 0.5);
% plot(x_axis_SS_down(:), y_axis_SS_down(:), 'LineWidth', 2, 'Color', [0.1 0.1 0.9])
train_data_logic_CS_down = CS_train_aligned(inds_event(inds_down,:));
% [x_axis_CS_down, y_axis_CS_down] = ESN_raster_plot_axes(train_data_logic_CS_down, inds_span, 3);
% plot(x_axis_CS_down(:), y_axis_CS_down(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
% xlim([min(inds_span)-1 max(inds_span)+1])
% ylim([(1-3) (size(train_data_logic,1)+3)])
% ylabel('Trials')
% title('Down')
% xlabel(xlabel_text_raster_);

% % Down Plot Velocity
% subplot(3,3,9)
% hold on
% train_data_logic_ = SS_train_aligned(inds_event(inds_down,:));
% firing_SS_down = mean(train_data_logic_) * 1000;
% plot(inds_span, ESN_smooth(firing_SS_down), 'LineWidth', 1, 'Color', [0.1 0.1 0.9])
% ylabel('SS Firing (spk/s)')
% ylim(range_SS_Firing)
% yyaxis right;
velocity_data_down = BEHAVE_eye_r_vm_filt(inds_event(inds_down,:));
% velocity_data_mean_down = nanmean(velocity_data_down);
% velocity_data_stdv_down = nanstd( velocity_data_down);
% plot(inds_span, velocity_data_mean_down+velocity_data_stdv_down, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
% plot(inds_span, velocity_data_mean_down-velocity_data_stdv_down, '-', 'LineWidth', 1, 'Color', [0.5 0.9 0.5])
% plot(inds_span, velocity_data_mean_down,                         '-', 'LineWidth', 2, 'Color', [0.1 0.9 0.1])
% ylim(range_vm)
% set(gca, 'YColor', [0.1 0.9 0.1])
% yyaxis left;
% xlim([min(inds_span)-1 max(inds_span)+1])
% title('Down')

% Probability
% subplot(3,3,5)
% hold on

% plot(cosd(0:5:360)*0.10, sind(0:5:360)*0.10, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
% plot(cosd(0:5:360)*0.20, sind(0:5:360)*0.20, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
% plot(cosd(0:5:360)*0.30, sind(0:5:360)*0.30, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])

% x_axis = [cosd(0) cosd(90) cosd(180) cosd(270) cosd(0)] .* prob_amplitude;
% y_axis = [sind(0) sind(90) sind(180) sind(270) sind(0)] .* prob_amplitude;
% plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', [0.9 0.1 0.1])
% axis equal
% xlabel(xlabel_text_CS_probab_);

% ESN_Beautify_Plot
% hFig = fig_handle_(fig_num_);
% figure_size  = [8.0 8.0];
% paper_margin = [0.1 0.1];
% paper_size = figure_size + 2 * paper_margin;
% set(hFig, 'PaperSize', paper_size);
% set(hFig, 'PaperPositionMode', 'manual');
% set(hFig, 'PaperPosition', [paper_margin figure_size]);
% set(hFig, 'Position', [[1 1] figure_size]);
% set(hFig, 'PaperOrientation', 'portrait');
if size(inds_event, 2) == size(train_data_logic_SS_top, 2)
raster_data.train_data_logic_SS_top = train_data_logic_SS_top;
raster_data.train_data_logic_CS_top = train_data_logic_CS_top;
raster_data.velocity_data_top       = velocity_data_top;
else
raster_data.train_data_logic_SS_top = train_data_logic_SS_top';
raster_data.train_data_logic_CS_top = train_data_logic_CS_top';
raster_data.velocity_data_top       = velocity_data_top';
end

if size(inds_event, 2) == size(train_data_logic_SS_left, 2)
raster_data.train_data_logic_SS_left = train_data_logic_SS_left;
raster_data.train_data_logic_CS_left = train_data_logic_CS_left;
raster_data.velocity_data_left       = velocity_data_left;
else
raster_data.train_data_logic_SS_left = train_data_logic_SS_left';
raster_data.train_data_logic_CS_left = train_data_logic_CS_left';
raster_data.velocity_data_left       = velocity_data_left';
end

if size(inds_event, 2) == size(train_data_logic_SS_right, 2)
raster_data.train_data_logic_SS_right = train_data_logic_SS_right;
raster_data.train_data_logic_CS_right = train_data_logic_CS_right;
raster_data.velocity_data_right       = velocity_data_right;
else
raster_data.train_data_logic_SS_right = train_data_logic_SS_right';
raster_data.train_data_logic_CS_right = train_data_logic_CS_right';
raster_data.velocity_data_right       = velocity_data_right';
end

if size(inds_event, 2) == size(train_data_logic_SS_down, 2)
raster_data.train_data_logic_SS_down = train_data_logic_SS_down;
raster_data.train_data_logic_CS_down = train_data_logic_CS_down;
raster_data.velocity_data_down       = velocity_data_down;
else
raster_data.train_data_logic_SS_down = train_data_logic_SS_down';
raster_data.train_data_logic_CS_down = train_data_logic_CS_down';
raster_data.velocity_data_down       = velocity_data_down';
end


end


