%% function temp_ESN_plot_data_compress_all_pcell
function temp_ESN_plot_data_compress_all_pcell
%% load ALL_PCELL_COMPRESSED_DATA
clearvars ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL RASTER_DATA_ALL_PCELL_TUNED bundle_inds
ALL_PCELL_COMPRESSED_DATA_file_name = 'ALL_PCELL_COMPRESSED_DATA.mat';
ALL_PCELL_COMPRESSED_DATA_file_path = './compress_pCells';
load([ALL_PCELL_COMPRESSED_DATA_file_path filesep ALL_PCELL_COMPRESSED_DATA_file_name]);

%% Build RASTER_DATA_ALL_PCELL
clearvars -except ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL RASTER_DATA_ALL_PCELL_TUNED bundle_inds
clearvars RASTER_DATA_ALL_PCELL
for counter_ALL_PCELL = 1 : length(ALL_PCELL_COMPRESSED_DATA)
    field_names_ALL_PCELL = fieldnames(ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL));
    field_names_ALL_PCELL = field_names_ALL_PCELL(contains(field_names_ALL_PCELL, 'raster_data_'));
    for counter_field_names_ALL_PCELL = 1 : 1 : length(field_names_ALL_PCELL)
        field_name_ALL_PCELL = field_names_ALL_PCELL{counter_field_names_ALL_PCELL};
        field_names_raster_data = fieldnames(ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).(field_name_ALL_PCELL));
        for counter_field_names_raster_data = 1 : 1 : length(field_names_raster_data)
            field_name_raster_data = field_names_raster_data{counter_field_names_raster_data};
            if contains(field_name_raster_data, '_sparse')
                continue;
            end
            RASTER_DATA_ALL_PCELL.(field_name_ALL_PCELL).(field_name_raster_data)(counter_ALL_PCELL, :) = ...
                ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).(field_name_ALL_PCELL).(field_name_raster_data);
        end
    end
end

train_data_names = {'train_data_logic_SS', 'train_data_logic_CS', 'velocity_data'};
field_names_RASTER_DATA_ALL_PCELL = fieldnames(RASTER_DATA_ALL_PCELL);
for counter_field_names_RASTER_DATA_ALL_PCELL = 1 : 1 : length(field_names_RASTER_DATA_ALL_PCELL)
    field_name_RASTER_DATA_ALL_PCELL = field_names_RASTER_DATA_ALL_PCELL{counter_field_names_RASTER_DATA_ALL_PCELL};
    for counter_train_data_names = 1 : 1 : length(train_data_names)
        train_data_name = train_data_names{counter_train_data_names};
        RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).([    train_data_name '_all']) = [ ...
            RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_right']) ; ...
            RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_top']) ; ...
            RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_left']) ; ...
            RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_down']) ; ...
            ];
    end
    train_data_name = 'train_data_logic_CS';
    RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).([    train_data_name '_all_numTrial']) = [ ...
        RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_right_numTrial']) ; ...
        RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_top_numTrial']) ; ...
        RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_left_numTrial']) ; ...
        RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_down_numTrial']) ; ...
        ];
    
    train_data_logic_SS_all = RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_SS_all;
    train_data_logic_CS_all = RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_CS_all;
    velocity_data_all       = RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).velocity_data_all;
    numTrial_all            = RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_CS_all_numTrial;
    
    RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_SS_all_mean = (numTrial_all' * train_data_logic_SS_all)./sum(numTrial_all);
    RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_CS_all_mean = (numTrial_all' * train_data_logic_CS_all)./sum(numTrial_all);
    RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).velocity_data_all_mean       = (numTrial_all' * velocity_data_all)      ./sum(numTrial_all);
    RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_SS_all_stdv = std(train_data_logic_SS_all, numTrial_all, 1);
    RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_CS_all_stdv = std(train_data_logic_CS_all, numTrial_all, 1);
    RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).velocity_data_all_stdv       = std(velocity_data_all,       numTrial_all, 1);
    
end

%% Plot-1 RASTER_DATA_ALL_PCELL
clearvars -except ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL RASTER_DATA_ALL_PCELL_TUNED bundle_inds
range_SS_Firing = [40 90];
range_vm        = [0 500];
Color_Vec = lines(7);
blue_color = Color_Vec(1, :);
red_color  = Color_Vec(7, :);
hFig = figure(1);
clf(hFig)
subplot(1, 4, 1)
hold on
yyaxis left;
% velocity_data_top = BEHAVE_eye_r_vm_filt(inds_event(inds_top,:));
inds_span          = nanmean(RASTER_DATA_ALL_PCELL.raster_data_cue_present.inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL.raster_data_cue_present.velocity_data_all_mean;
velocity_data_stdv = RASTER_DATA_ALL_PCELL.raster_data_cue_present.velocity_data_all_stdv;
firing_SS_mean     = RASTER_DATA_ALL_PCELL.raster_data_cue_present.train_data_logic_SS_all_mean * 1000;
firing_SS_stdv     = RASTER_DATA_ALL_PCELL.raster_data_cue_present.train_data_logic_SS_all_stdv * 1000 ./ sqrt(66);
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Cue Pres')
set(gca, 'YColor', red_color)

subplot(1, 4, 2)
hold on
yyaxis left;
% velocity_data_top = BEHAVE_eye_r_vm_filt(inds_event(inds_top,:));
inds_span          = nanmean(RASTER_DATA_ALL_PCELL.raster_data_primSac_onset.inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL.raster_data_primSac_onset.velocity_data_all_mean;
velocity_data_stdv = RASTER_DATA_ALL_PCELL.raster_data_primSac_onset.velocity_data_all_stdv;
firing_SS_mean     = RASTER_DATA_ALL_PCELL.raster_data_primSac_onset.train_data_logic_SS_all_mean * 1000;
firing_SS_stdv     = RASTER_DATA_ALL_PCELL.raster_data_primSac_onset.train_data_logic_SS_all_stdv * 1000 ./ sqrt(66);
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Prim Sac On')
set(gca, 'YColor', red_color)

subplot(1, 4, 3)
hold on
yyaxis left;
% velocity_data_top = BEHAVE_eye_r_vm_filt(inds_event(inds_top,:));
inds_span          = nanmean(RASTER_DATA_ALL_PCELL.raster_data_primSac_offset.inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL.raster_data_primSac_offset.velocity_data_all_mean;
velocity_data_stdv = RASTER_DATA_ALL_PCELL.raster_data_primSac_offset.velocity_data_all_stdv;
firing_SS_mean     = RASTER_DATA_ALL_PCELL.raster_data_primSac_offset.train_data_logic_SS_all_mean * 1000;
firing_SS_stdv     = RASTER_DATA_ALL_PCELL.raster_data_primSac_offset.train_data_logic_SS_all_stdv * 1000 ./ sqrt(66);
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Prim Sac Off')
set(gca, 'YColor', red_color)

subplot(1, 4, 4)
hold on
yyaxis left;
% velocity_data_top = BEHAVE_eye_r_vm_filt(inds_event(inds_top,:));
inds_span          = nanmean(RASTER_DATA_ALL_PCELL.raster_data_corrSac_onset.inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL.raster_data_corrSac_onset.velocity_data_all_mean;
velocity_data_stdv = RASTER_DATA_ALL_PCELL.raster_data_corrSac_onset.velocity_data_all_stdv;
firing_SS_mean     = RASTER_DATA_ALL_PCELL.raster_data_corrSac_onset.train_data_logic_SS_all_mean * 1000;
firing_SS_stdv     = RASTER_DATA_ALL_PCELL.raster_data_corrSac_onset.train_data_logic_SS_all_stdv * 1000 ./ sqrt(66);
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
title('Corr Sac On')
set(gca, 'YColor', red_color)

ESN_Beautify_Plot
figure_size  = [16.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

%% Build RASTER_DATA_ALL_PCELL_TUNED
clearvars -except ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL RASTER_DATA_ALL_PCELL_TUNED bundle_inds
clearvars RASTER_DATA_ALL_PCELL_TUNED
train_data_names = {'train_data_logic_SS', 'train_data_logic_CS', 'velocity_data'};
train_data_dirs  = {'_right', '_top', '_left', '_down'};
for counter_ALL_PCELL = 1 : length(ALL_PCELL_COMPRESSED_DATA)
    field_names_ALL_PCELL = fieldnames(ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL));
    field_names_ALL_PCELL = field_names_ALL_PCELL(contains(field_names_ALL_PCELL, 'raster_data_'));
    overall_prob_tuning_bundle = ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).CS_Tuning.overall_prob_tuning_bundle;
    for counter_field_names_ALL_PCELL = 1 : 1 : length(field_names_ALL_PCELL)
        field_name_ALL_PCELL = field_names_ALL_PCELL{counter_field_names_ALL_PCELL};
        for counter_train_data_dirs = 1 : 1 : length(train_data_dirs)
            train_data_dir = train_data_dirs{counter_train_data_dirs};
            train_data_dir_tuned = train_data_dirs{overall_prob_tuning_bundle(counter_train_data_dirs)};
            for counter_train_data_names = 1 : 1 : length(train_data_names)
                train_data_name = train_data_names{counter_train_data_names};
                
                RASTER_DATA_ALL_PCELL_TUNED.(field_name_ALL_PCELL).([train_data_name train_data_dir])(counter_ALL_PCELL, :) = ...
                    ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).(field_name_ALL_PCELL).([train_data_name train_data_dir_tuned]);
            end
            RASTER_DATA_ALL_PCELL_TUNED.(field_name_ALL_PCELL).(['train_data_logic_CS' train_data_dir '_numTrial'])(counter_ALL_PCELL, :) = ...
                ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).(field_name_ALL_PCELL).(['train_data_logic_CS' train_data_dir_tuned '_numTrial']);
        end
        RASTER_DATA_ALL_PCELL_TUNED.(field_name_ALL_PCELL).('inds_span')(counter_ALL_PCELL, :) = ...
            ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).(field_name_ALL_PCELL).('inds_span');
    end
end

field_names_RASTER_DATA_ALL_PCELL = fieldnames(RASTER_DATA_ALL_PCELL_TUNED);
for counter_field_names_RASTER_DATA_ALL_PCELL = 1 : 1 : length(field_names_RASTER_DATA_ALL_PCELL)
    field_name_RASTER_DATA_ALL_PCELL = field_names_RASTER_DATA_ALL_PCELL{counter_field_names_RASTER_DATA_ALL_PCELL};
    for counter_train_data_names = 1 : 1 : length(train_data_names)
        train_data_name = train_data_names{counter_train_data_names};
        RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([    train_data_name '_all']) = [ ...
            RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_right']) ; ...
            RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_top']) ; ...
            RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_left']) ; ...
            RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_down']) ; ...
            ];
    end
    train_data_name = 'train_data_logic_CS';
    RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([    train_data_name '_all_numTrial']) = [ ...
        RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_right_numTrial']) ; ...
        RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_top_numTrial']) ; ...
        RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_left_numTrial']) ; ...
        RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_down_numTrial']) ; ...
        ];
    
    train_data_logic_SS_all = RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_SS_all;
    train_data_logic_CS_all = RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_CS_all;
    velocity_data_all       = RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).velocity_data_all;
    numTrial_all            = RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_CS_all_numTrial;
    
    RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_SS_all_mean = (numTrial_all' * train_data_logic_SS_all)./sum(numTrial_all);
    RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_CS_all_mean = (numTrial_all' * train_data_logic_CS_all)./sum(numTrial_all);
    RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).velocity_data_all_mean       = (numTrial_all' * velocity_data_all)      ./sum(numTrial_all);
    RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_SS_all_stdv = std(train_data_logic_SS_all, numTrial_all, 1);
    RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_CS_all_stdv = std(train_data_logic_CS_all, numTrial_all, 1);
    RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).velocity_data_all_stdv       = std(velocity_data_all,       numTrial_all, 1);
    
end

%% Plot-01 RASTER_DATA_ALL_PCELL_TUNED SS Firing
clearvars -except ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL RASTER_DATA_ALL_PCELL_TUNED bundle_inds
range_SS_Firing = [40 100];
range_vm        = [0 500];
Color_Vec = lines(7);
blue_color = Color_Vec(1, :);
red_color  = Color_Vec(7, :);
hFig = figure(1);
clf(hFig)
subplot(1, 4, 1)
hold on
yyaxis left;
% velocity_data_top = BEHAVE_eye_r_vm_filt(inds_event(inds_top,:));
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.raster_data_cue_present.inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.raster_data_cue_present.velocity_data_all_mean;
velocity_data_stdv = RASTER_DATA_ALL_PCELL_TUNED.raster_data_cue_present.velocity_data_all_stdv;
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.raster_data_cue_present.train_data_logic_SS_all_mean * 1000;
firing_SS_stdv     = RASTER_DATA_ALL_PCELL_TUNED.raster_data_cue_present.train_data_logic_SS_all_stdv * 1000 ./ sqrt(66);
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Cue Pres')
set(gca, 'YColor', red_color)

subplot(1, 4, 2)
hold on
yyaxis left;
% velocity_data_top = BEHAVE_eye_r_vm_filt(inds_event(inds_top,:));
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_onset.inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_onset.velocity_data_all_mean;
velocity_data_stdv = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_onset.velocity_data_all_stdv;
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_onset.train_data_logic_SS_all_mean * 1000;
firing_SS_stdv     = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_onset.train_data_logic_SS_all_stdv * 1000 ./ sqrt(66);
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac On')
set(gca, 'YColor', red_color)

subplot(1, 4, 3)
hold on
yyaxis left;
% velocity_data_top = BEHAVE_eye_r_vm_filt(inds_event(inds_top,:));
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_offset.inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_offset.velocity_data_all_mean;
velocity_data_stdv = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_offset.velocity_data_all_stdv;
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_offset.train_data_logic_SS_all_mean * 1000;
firing_SS_stdv     = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_offset.train_data_logic_SS_all_stdv * 1000 ./ sqrt(66);
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac Off')
set(gca, 'YColor', red_color)

subplot(1, 4, 4)
hold on
yyaxis left;
% velocity_data_top = BEHAVE_eye_r_vm_filt(inds_event(inds_top,:));
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.raster_data_corrSac_onset.inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.raster_data_corrSac_onset.velocity_data_all_mean;
velocity_data_stdv = RASTER_DATA_ALL_PCELL_TUNED.raster_data_corrSac_onset.velocity_data_all_stdv;
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.raster_data_corrSac_onset.train_data_logic_SS_all_mean * 1000;
firing_SS_stdv     = RASTER_DATA_ALL_PCELL_TUNED.raster_data_corrSac_onset.train_data_logic_SS_all_stdv * 1000 ./ sqrt(66);
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Corr Sac On')
set(gca, 'YColor', red_color)

ESN_Beautify_Plot
figure_size  = [16.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

sgtitle(hFig, 'All directions', 'Interpreter', 'none');

%% Plot-02 RASTER_DATA_ALL_PCELL_TUNED CS Firing
clearvars -except ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL RASTER_DATA_ALL_PCELL_TUNED bundle_inds
range_CS_Firing = [0 5];
range_vm        = [0 500];
Color_Vec = lines(7);
blue_color = Color_Vec(1, :);
red_color  = Color_Vec(7, :);
hFig = figure(2);
clf(hFig)
subplot(1, 4, 1)
hold on
yyaxis left;
% velocity_data_top = BEHAVE_eye_r_vm_filt(inds_event(inds_top,:));
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.raster_data_cue_present.inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.raster_data_cue_present.velocity_data_all_mean;
velocity_data_stdv = RASTER_DATA_ALL_PCELL_TUNED.raster_data_cue_present.velocity_data_all_stdv;
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.raster_data_cue_present.train_data_logic_CS_all_mean * 1000;
firing_CS_stdv     = RASTER_DATA_ALL_PCELL_TUNED.raster_data_cue_present.train_data_logic_CS_all_stdv * 1000 ./ sqrt(66);
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Cue Pres')
set(gca, 'YColor', red_color)

subplot(1, 4, 2)
hold on
yyaxis left;
% velocity_data_top = BEHAVE_eye_r_vm_filt(inds_event(inds_top,:));
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_onset.inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_onset.velocity_data_all_mean;
velocity_data_stdv = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_onset.velocity_data_all_stdv;
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_onset.train_data_logic_CS_all_mean * 1000;
firing_CS_stdv     = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_onset.train_data_logic_CS_all_stdv * 1000 ./ sqrt(66);
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac On')
set(gca, 'YColor', red_color)

subplot(1, 4, 3)
hold on
yyaxis left;
% velocity_data_top = BEHAVE_eye_r_vm_filt(inds_event(inds_top,:));
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_offset.inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_offset.velocity_data_all_mean;
velocity_data_stdv = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_offset.velocity_data_all_stdv;
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_offset.train_data_logic_CS_all_mean * 1000;
firing_CS_stdv     = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_offset.train_data_logic_CS_all_stdv * 1000 ./ sqrt(66);
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac Off')
set(gca, 'YColor', red_color)

subplot(1, 4, 4)
hold on
yyaxis left;
% velocity_data_top = BEHAVE_eye_r_vm_filt(inds_event(inds_top,:));
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.raster_data_corrSac_onset.inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.raster_data_corrSac_onset.velocity_data_all_mean;
velocity_data_stdv = RASTER_DATA_ALL_PCELL_TUNED.raster_data_corrSac_onset.velocity_data_all_stdv;
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.raster_data_corrSac_onset.train_data_logic_CS_all_mean * 1000;
firing_CS_stdv     = RASTER_DATA_ALL_PCELL_TUNED.raster_data_corrSac_onset.train_data_logic_CS_all_stdv * 1000 ./ sqrt(66);
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Corr Sac On')
set(gca, 'YColor', red_color)

ESN_Beautify_Plot
figure_size  = [16.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

sgtitle(hFig, 'All directions', 'Interpreter', 'none');

%% Plot-03 RASTER_DATA_ALL_PCELL_TUNED CS-000 SS Firing
clearvars -except ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL RASTER_DATA_ALL_PCELL_TUNED bundle_inds
range_SS_Firing = [40 100];
range_vm        = [0 500];
Color_Vec = lines(7);
blue_color = Color_Vec(1, :);
red_color  = Color_Vec(7, :);
hFig = figure(3);
clf(hFig)
subplot(1, 4, 1)
hold on
yyaxis left;
raster_data_name = 'raster_data_cue_present';
raster_data_dir  = '_right';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_SS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Cue Pres')
set(gca, 'YColor', red_color)

subplot(1, 4, 2)
hold on
yyaxis left;
raster_data_name = 'raster_data_primSac_onset';
raster_data_dir  = '_right';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_SS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac On')
set(gca, 'YColor', red_color)

subplot(1, 4, 3)
hold on
yyaxis left;
raster_data_name = 'raster_data_primSac_offset';
raster_data_dir  = '_right';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_SS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac Off')
set(gca, 'YColor', red_color)

subplot(1, 4, 4)
hold on
yyaxis left;
raster_data_name = 'raster_data_corrSac_onset';
raster_data_dir  = '_right';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_SS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Corr Sac On')
set(gca, 'YColor', red_color)

ESN_Beautify_Plot
figure_size  = [16.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

sgtitle(hFig, 'CS-000', 'Interpreter', 'none');

%% Plot-04 RASTER_DATA_ALL_PCELL_TUNED CS-090 SS Firing
clearvars -except ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL RASTER_DATA_ALL_PCELL_TUNED bundle_inds
range_SS_Firing = [40 100];
range_vm        = [0 500];
Color_Vec = lines(7);
blue_color = Color_Vec(1, :);
red_color  = Color_Vec(7, :);
hFig = figure(4);
clf(hFig)
subplot(1, 4, 1)
hold on
yyaxis left;
raster_data_name = 'raster_data_cue_present';
raster_data_dir  = '_top';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_SS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Cue Pres')
set(gca, 'YColor', red_color)

subplot(1, 4, 2)
hold on
yyaxis left;
raster_data_name = 'raster_data_primSac_onset';
raster_data_dir  = '_top';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_SS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac On')
set(gca, 'YColor', red_color)

subplot(1, 4, 3)
hold on
yyaxis left;
raster_data_name = 'raster_data_primSac_offset';
raster_data_dir  = '_top';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_SS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac Off')
set(gca, 'YColor', red_color)

subplot(1, 4, 4)
hold on
yyaxis left;
raster_data_name = 'raster_data_corrSac_onset';
raster_data_dir  = '_top';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_SS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Corr Sac On')
set(gca, 'YColor', red_color)

ESN_Beautify_Plot
figure_size  = [16.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

sgtitle(hFig, 'CS-090', 'Interpreter', 'none');

%% Plot-05 RASTER_DATA_ALL_PCELL_TUNED CS-180 SS Firing
clearvars -except ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL RASTER_DATA_ALL_PCELL_TUNED bundle_inds
range_SS_Firing = [40 100];
range_vm        = [0 500];
Color_Vec = lines(7);
blue_color = Color_Vec(1, :);
red_color  = Color_Vec(7, :);
hFig = figure(5);
clf(hFig)
subplot(1, 4, 1)
hold on
yyaxis left;
raster_data_name = 'raster_data_cue_present';
raster_data_dir  = '_left';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_SS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Cue Pres')
set(gca, 'YColor', red_color)

subplot(1, 4, 2)
hold on
yyaxis left;
raster_data_name = 'raster_data_primSac_onset';
raster_data_dir  = '_left';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_SS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac On')
set(gca, 'YColor', red_color)

subplot(1, 4, 3)
hold on
yyaxis left;
raster_data_name = 'raster_data_primSac_offset';
raster_data_dir  = '_left';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_SS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac Off')
set(gca, 'YColor', red_color)

subplot(1, 4, 4)
hold on
yyaxis left;
raster_data_name = 'raster_data_corrSac_onset';
raster_data_dir  = '_left';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_SS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Corr Sac On')
set(gca, 'YColor', red_color)

ESN_Beautify_Plot
figure_size  = [16.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

sgtitle(hFig, 'CS-180', 'Interpreter', 'none');

%% Plot-06 RASTER_DATA_ALL_PCELL_TUNED CS-270 SS Firing
clearvars -except ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL RASTER_DATA_ALL_PCELL_TUNED bundle_inds
range_SS_Firing = [40 100];
range_vm        = [0 500];
Color_Vec = lines(7);
blue_color = Color_Vec(1, :);
red_color  = Color_Vec(7, :);
hFig = figure(6);
clf(hFig)
subplot(1, 4, 1)
hold on
yyaxis left;
raster_data_name = 'raster_data_cue_present';
raster_data_dir  = '_down';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_SS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Cue Pres')
set(gca, 'YColor', red_color)

subplot(1, 4, 2)
hold on
yyaxis left;
raster_data_name = 'raster_data_primSac_onset';
raster_data_dir  = '_down';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_SS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac On')
set(gca, 'YColor', red_color)

subplot(1, 4, 3)
hold on
yyaxis left;
raster_data_name = 'raster_data_primSac_offset';
raster_data_dir  = '_down';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_SS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac Off')
set(gca, 'YColor', red_color)

subplot(1, 4, 4)
hold on
yyaxis left;
raster_data_name = 'raster_data_corrSac_onset';
raster_data_dir  = '_down';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_SS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_SS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_SS_mean+firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_SS_mean-firing_SS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_SS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
ylabel('SS Firing (spk/s)')
ylim(range_SS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Corr Sac On')
set(gca, 'YColor', red_color)

ESN_Beautify_Plot
figure_size  = [16.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

sgtitle(hFig, 'CS-270', 'Interpreter', 'none');

%% Plot-07 RASTER_DATA_ALL_PCELL_TUNED CS-000 CS-Firing
clearvars -except ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL RASTER_DATA_ALL_PCELL_TUNED bundle_inds
range_CS_Firing = [0 5.0];
range_vm        = [0 500];
Color_Vec = lines(7);
blue_color = Color_Vec(1, :);
red_color  = Color_Vec(7, :);
hFig = figure(7);
clf(hFig)
subplot(1, 4, 1)
hold on
yyaxis left;
raster_data_name = 'raster_data_cue_present';
raster_data_dir  = '_right';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Cue Pres')
set(gca, 'YColor', red_color)

subplot(1, 4, 2)
hold on
yyaxis left;
raster_data_name = 'raster_data_primSac_onset';
raster_data_dir  = '_right';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac On')
set(gca, 'YColor', red_color)

subplot(1, 4, 3)
hold on
yyaxis left;
raster_data_name = 'raster_data_primSac_offset';
raster_data_dir  = '_right';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac Off')
set(gca, 'YColor', red_color)

subplot(1, 4, 4)
hold on
yyaxis left;
raster_data_name = 'raster_data_corrSac_onset';
raster_data_dir  = '_right';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Corr Sac On')
set(gca, 'YColor', red_color)

ESN_Beautify_Plot
figure_size  = [16.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

sgtitle(hFig, 'CS-000', 'Interpreter', 'none');

%% Plot-08 RASTER_DATA_ALL_PCELL_TUNED CS-090 CS Firing
clearvars -except ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL RASTER_DATA_ALL_PCELL_TUNED bundle_inds
range_CS_Firing = [0 5];
range_vm        = [0 500];
Color_Vec = lines(7);
blue_color = Color_Vec(1, :);
red_color  = Color_Vec(7, :);
hFig = figure(8);
clf(hFig)
subplot(1, 4, 1)
hold on
yyaxis left;
raster_data_name = 'raster_data_cue_present';
raster_data_dir  = '_top';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Cue Pres')
set(gca, 'YColor', red_color)

subplot(1, 4, 2)
hold on
yyaxis left;
raster_data_name = 'raster_data_primSac_onset';
raster_data_dir  = '_top';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac On')
set(gca, 'YColor', red_color)

subplot(1, 4, 3)
hold on
yyaxis left;
raster_data_name = 'raster_data_primSac_offset';
raster_data_dir  = '_top';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac Off')
set(gca, 'YColor', red_color)

subplot(1, 4, 4)
hold on
yyaxis left;
raster_data_name = 'raster_data_corrSac_onset';
raster_data_dir  = '_top';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Corr Sac On')
set(gca, 'YColor', red_color)

ESN_Beautify_Plot
figure_size  = [16.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

sgtitle(hFig, 'CS-090', 'Interpreter', 'none');

%% Plot-09 RASTER_DATA_ALL_PCELL_TUNED CS-180 CS Firing
clearvars -except ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL RASTER_DATA_ALL_PCELL_TUNED bundle_inds
range_CS_Firing = [0 5];
range_vm        = [0 500];
Color_Vec = lines(7);
blue_color = Color_Vec(1, :);
red_color  = Color_Vec(7, :);
hFig = figure(9);
clf(hFig)
subplot(1, 4, 1)
hold on
yyaxis left;
raster_data_name = 'raster_data_cue_present';
raster_data_dir  = '_left';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Cue Pres')
set(gca, 'YColor', red_color)

subplot(1, 4, 2)
hold on
yyaxis left;
raster_data_name = 'raster_data_primSac_onset';
raster_data_dir  = '_left';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac On')
set(gca, 'YColor', red_color)

subplot(1, 4, 3)
hold on
yyaxis left;
raster_data_name = 'raster_data_primSac_offset';
raster_data_dir  = '_left';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac Off')
set(gca, 'YColor', red_color)

subplot(1, 4, 4)
hold on
yyaxis left;
raster_data_name = 'raster_data_corrSac_onset';
raster_data_dir  = '_left';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Corr Sac On')
set(gca, 'YColor', red_color)

ESN_Beautify_Plot
figure_size  = [16.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

sgtitle(hFig, 'CS-180', 'Interpreter', 'none');

%% Plot-10 RASTER_DATA_ALL_PCELL_TUNED CS-270 CS Firing
clearvars -except ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL RASTER_DATA_ALL_PCELL_TUNED bundle_inds
range_CS_Firing = [0 5];
range_vm        = [0 500];
Color_Vec = lines(7);
blue_color = Color_Vec(1, :);
red_color  = Color_Vec(7, :);
hFig = figure(10);
clf(hFig)
subplot(1, 4, 1)
hold on
yyaxis left;
raster_data_name = 'raster_data_cue_present';
raster_data_dir  = '_down';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Cue Pres')
set(gca, 'YColor', red_color)

subplot(1, 4, 2)
hold on
yyaxis left;
raster_data_name = 'raster_data_primSac_onset';
raster_data_dir  = '_down';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac On')
set(gca, 'YColor', red_color)

subplot(1, 4, 3)
hold on
yyaxis left;
raster_data_name = 'raster_data_primSac_offset';
raster_data_dir  = '_down';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
% ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Prim Sac Off')
set(gca, 'YColor', red_color)

subplot(1, 4, 4)
hold on
yyaxis left;
raster_data_name = 'raster_data_corrSac_onset';
raster_data_dir  = '_down';
inds_span          = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).inds_span);
velocity_data_mean = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['velocity_data' raster_data_dir]) ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
firing_CS_mean     = RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial'])' * RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir]) .* 1000 ./ sum(RASTER_DATA_ALL_PCELL_TUNED.(raster_data_name).(['train_data_logic_CS' raster_data_dir '_numTrial']));
% plot(inds_span, velocity_data_mean+velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
% plot(inds_span, velocity_data_mean-velocity_data_stdv, '-', 'LineWidth', 1, 'Color', red_color)
plot(inds_span, velocity_data_mean,                    '-', 'LineWidth', 2, 'Color', red_color)
% ylabel('Eye velocity (deg/s)')
ylim(range_vm)
set(gca, 'YColor', red_color)
yyaxis right;
% plot(inds_span, ESN_smooth(firing_CS_mean+firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
% plot(inds_span, ESN_smooth(firing_CS_mean-firing_CS_stdv), '-', 'LineWidth', 1, 'Color', blue_color)
plot(inds_span, ESN_smooth(firing_CS_mean),                '-', 'LineWidth', 2, 'Color', blue_color)
ylabel('CS Firing (spk/s)')
ylim(range_CS_Firing)
set(gca, 'YColor', blue_color)
yyaxis left;
xlim([min(inds_span)-1 max(inds_span)+1])
xlabel('Corr Sac On')
set(gca, 'YColor', red_color)

ESN_Beautify_Plot
figure_size  = [16.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

sgtitle(hFig, 'CS-270', 'Interpreter', 'none');

%% Plot-11 CS_TUNING
cue_present_prob_ALL_PCELL    = nan(length(ALL_PCELL_COMPRESSED_DATA), 4);
primSac_onset_prob_ALL_PCELL  = nan(length(ALL_PCELL_COMPRESSED_DATA), 4);
primSac_offset_prob_ALL_PCELL = nan(length(ALL_PCELL_COMPRESSED_DATA), 4);
corrSac_onset_prob_ALL_PCELL  = nan(length(ALL_PCELL_COMPRESSED_DATA), 4);
numTrials_ALL_PCELL = nan(length(ALL_PCELL_COMPRESSED_DATA), 1);
for counter_ALL_PCELL_COMPRESSED_DATA = 1 : 1 : length(ALL_PCELL_COMPRESSED_DATA)
    cue_present_prob_ALL_PCELL(counter_ALL_PCELL_COMPRESSED_DATA, :) = ...
        ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL_COMPRESSED_DATA).CS_Tuning.cue_present_prob( ...
        ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL_COMPRESSED_DATA).CS_Tuning.overall_prob_tuning_bundle ... 
        );
    primSac_onset_prob_ALL_PCELL(counter_ALL_PCELL_COMPRESSED_DATA, :) = ...
        ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL_COMPRESSED_DATA).CS_Tuning.primSac_onset_prob( ...
        ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL_COMPRESSED_DATA).CS_Tuning.overall_prob_tuning_bundle ... 
        );
    primSac_offset_prob_ALL_PCELL(counter_ALL_PCELL_COMPRESSED_DATA, :) = ...
        ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL_COMPRESSED_DATA).CS_Tuning.primSac_offset_prob( ...
        ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL_COMPRESSED_DATA).CS_Tuning.overall_prob_tuning_bundle ... 
        );
    corrSac_onset_prob_ALL_PCELL(counter_ALL_PCELL_COMPRESSED_DATA, :) = ...
        ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL_COMPRESSED_DATA).CS_Tuning.corrSac_onset_prob( ...
        ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL_COMPRESSED_DATA).CS_Tuning.overall_prob_tuning_bundle ... 
        );
    numTrials_ALL_PCELL(counter_ALL_PCELL_COMPRESSED_DATA, :) = ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL_COMPRESSED_DATA).CS_Tuning.numTrials_bundle;
end
cue_present_prob_ALL_PCELL_mean    = numTrials_ALL_PCELL' * cue_present_prob_ALL_PCELL    ./ sum(numTrials_ALL_PCELL);
primSac_onset_prob_ALL_PCELL_mean  = numTrials_ALL_PCELL' * primSac_onset_prob_ALL_PCELL  ./ sum(numTrials_ALL_PCELL);
primSac_offset_prob_ALL_PCELL_mean = numTrials_ALL_PCELL' * primSac_offset_prob_ALL_PCELL ./ sum(numTrials_ALL_PCELL);
corrSac_onset_prob_ALL_PCELL_mean  = numTrials_ALL_PCELL' * corrSac_onset_prob_ALL_PCELL  ./ sum(numTrials_ALL_PCELL);


overall_prob_ALL_PCELL = nan(length(bundle_inds), 4);
numTrials_ALL_PCELL = nan(length(bundle_inds), 1);
for counter_bundle_inds = 1 : 1 : length(bundle_inds)
    bundle_ind = bundle_inds{counter_bundle_inds}(1);
    overall_prob_ALL_PCELL(counter_bundle_inds, :) = ...
        ALL_PCELL_COMPRESSED_DATA(bundle_ind).CS_Tuning.overall_prob_bundle(ALL_PCELL_COMPRESSED_DATA(bundle_ind).CS_Tuning.overall_prob_tuning_bundle);
    numTrials_ALL_PCELL(counter_bundle_inds, :) = ALL_PCELL_COMPRESSED_DATA(bundle_ind).CS_Tuning.numTrials_bundle;
end
overall_prob_ALL_PCELL_mean = numTrials_ALL_PCELL' * overall_prob_ALL_PCELL ./ sum(numTrials_ALL_PCELL);


plot_data_amp = [overall_prob_ALL_PCELL overall_prob_ALL_PCELL(:,1) nan(length(bundle_inds), 1)];
plot_data_deg = repmat([0, 90, 180, 270, 0, nan], length(bundle_inds), 1);
plot_data_x_axis = plot_data_amp .* cosd(plot_data_deg);
plot_data_y_axis = plot_data_amp .* sind(plot_data_deg);
plot_data_x_axis = plot_data_x_axis';
plot_data_y_axis = plot_data_y_axis';
plot_data_x_axis = plot_data_x_axis(:);
plot_data_y_axis = plot_data_y_axis(:);

Color_Vec = lines(7);
% blue_color = Color_Vec(7, :);
red_color  = Color_Vec(7, :);
hFig = figure(11);
clf(hFig)

subplot(1,5,1)
hold on
plot([0.4, -0.4, nan, 0, -0, nan,]', [0, 0, nan, 0.4, -0.4, nan], 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.10, sind(0:5:360)*0.10, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.20, sind(0:5:360)*0.20, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.30, sind(0:5:360)*0.30, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.40, sind(0:5:360)*0.40, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
line_width = 0.5;
line_alpha = 0.30;
% plot(plot_data_x_axis, plot_data_y_axis, '-', 'Color', red_color);
patch(plot_data_x_axis, plot_data_y_axis,'w',...
        'linewidth',line_width,...
        'EdgeAlpha',line_alpha,...
        'EdgeColor',red_color);
line_width = 2.5;
line_alpha = 0.95;
plot_data_amp_mean = [overall_prob_ALL_PCELL_mean, overall_prob_ALL_PCELL_mean(1), nan]';
plot_data_deg_mean = [0, 90, 180, 270, 0, nan]';
plot_data_x_axis_mean = plot_data_amp_mean .* cosd(plot_data_deg_mean);
plot_data_y_axis_mean = plot_data_amp_mean .* sind(plot_data_deg_mean);
% plot(plot_data_x_axis_mean, plot_data_y_axis_mean, '-', 'LineWidth', 3, 'Color', blue_color)
patch(plot_data_x_axis_mean, plot_data_y_axis_mean,'w',...
        'linewidth',line_width,...
        'EdgeAlpha',line_alpha,...
        'EdgeColor',red_color);
axis equal
xlim([-0.4 0.4])
ylim([-0.4 0.4])
xlabel('Overall');
ylabel('CS probability');
set(gca, 'XTick', -0.4:0.2:0.4, 'YTick', -0.4:0.2:0.4)

subplot(1,5,2)
hold on
plot([0.4, -0.4, nan, 0, -0, nan,]', [0, 0, nan, 0.4, -0.4, nan], 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.10, sind(0:5:360)*0.10, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.20, sind(0:5:360)*0.20, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.30, sind(0:5:360)*0.30, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.40, sind(0:5:360)*0.40, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
% line_width = 0.5;
% line_alpha = 0.30;
% plot(plot_data_x_axis, plot_data_y_axis, '-', 'Color', red_color);
% patch(plot_data_x_axis, plot_data_y_axis,'w',...
%         'linewidth',line_width,...
%         'EdgeAlpha',line_alpha,...
%         'EdgeColor',red_color);
line_width = 2.5;
line_alpha = 0.95;
plot_data_amp_mean = [cue_present_prob_ALL_PCELL_mean cue_present_prob_ALL_PCELL_mean(1), nan]';
plot_data_deg_mean = [0, 90, 180, 270, 0, nan]';
plot_data_x_axis_mean = plot_data_amp_mean .* cosd(plot_data_deg_mean);
plot_data_y_axis_mean = plot_data_amp_mean .* sind(plot_data_deg_mean);
% plot(plot_data_x_axis_mean, plot_data_y_axis_mean, '-', 'LineWidth', 3, 'Color', blue_color)
patch(plot_data_x_axis_mean, plot_data_y_axis_mean,'w',...
        'linewidth',line_width,...
        'EdgeAlpha',line_alpha,...
        'EdgeColor',red_color);
axis equal
xlim([-0.4 0.4])
ylim([-0.4 0.4])
xlabel('Cue Pres');
% ylabel('CS probability');
set(gca, 'XTick', -0.4:0.2:0.4, 'YTick', -0.4:0.2:0.4)

subplot(1,5,3)
hold on
plot([0.4, -0.4, nan, 0, -0, nan,]', [0, 0, nan, 0.4, -0.4, nan], 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.10, sind(0:5:360)*0.10, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.20, sind(0:5:360)*0.20, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.30, sind(0:5:360)*0.30, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.40, sind(0:5:360)*0.40, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
% line_width = 0.5;
% line_alpha = 0.30;
% plot(plot_data_x_axis, plot_data_y_axis, '-', 'Color', red_color);
% patch(plot_data_x_axis, plot_data_y_axis,'w',...
%         'linewidth',line_width,...
%         'EdgeAlpha',line_alpha,...
%         'EdgeColor',red_color);
line_width = 2.5;
line_alpha = 0.95;
plot_data_amp_mean = [primSac_onset_prob_ALL_PCELL_mean primSac_onset_prob_ALL_PCELL_mean(1), nan]';
plot_data_deg_mean = [0, 90, 180, 270, 0, nan]';
plot_data_x_axis_mean = plot_data_amp_mean .* cosd(plot_data_deg_mean);
plot_data_y_axis_mean = plot_data_amp_mean .* sind(plot_data_deg_mean);
% plot(plot_data_x_axis_mean, plot_data_y_axis_mean, '-', 'LineWidth', 3, 'Color', blue_color)
patch(plot_data_x_axis_mean, plot_data_y_axis_mean,'w',...
        'linewidth',line_width,...
        'EdgeAlpha',line_alpha,...
        'EdgeColor',red_color);
axis equal
xlim([-0.4 0.4])
ylim([-0.4 0.4])
xlabel('Prim On');
% ylabel('CS probability');
set(gca, 'XTick', -0.4:0.2:0.4, 'YTick', -0.4:0.2:0.4)

subplot(1,5,4)
hold on
plot([0.4, -0.4, nan, 0, -0, nan,]', [0, 0, nan, 0.4, -0.4, nan], 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.10, sind(0:5:360)*0.10, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.20, sind(0:5:360)*0.20, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.30, sind(0:5:360)*0.30, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.40, sind(0:5:360)*0.40, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
% line_width = 0.5;
% line_alpha = 0.30;
% plot(plot_data_x_axis, plot_data_y_axis, '-', 'Color', red_color);
% patch(plot_data_x_axis, plot_data_y_axis,'w',...
%         'linewidth',line_width,...
%         'EdgeAlpha',line_alpha,...
%         'EdgeColor',red_color);
line_width = 2.5;
line_alpha = 0.95;
plot_data_amp_mean = [primSac_offset_prob_ALL_PCELL_mean primSac_offset_prob_ALL_PCELL_mean(1), nan]';
plot_data_deg_mean = [0, 90, 180, 270, 0, nan]';
plot_data_x_axis_mean = plot_data_amp_mean .* cosd(plot_data_deg_mean);
plot_data_y_axis_mean = plot_data_amp_mean .* sind(plot_data_deg_mean);
% plot(plot_data_x_axis_mean, plot_data_y_axis_mean, '-', 'LineWidth', 3, 'Color', blue_color)
patch(plot_data_x_axis_mean, plot_data_y_axis_mean,'w',...
        'linewidth',line_width,...
        'EdgeAlpha',line_alpha,...
        'EdgeColor',red_color);
axis equal
xlim([-0.4 0.4])
ylim([-0.4 0.4])
xlabel('Prim Off');
% ylabel('CS probability');
set(gca, 'XTick', -0.4:0.2:0.4, 'YTick', -0.4:0.2:0.4)

subplot(1,5,5)
hold on
plot([0.4, -0.4, nan, 0, -0, nan,]', [0, 0, nan, 0.4, -0.4, nan], 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.10, sind(0:5:360)*0.10, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.20, sind(0:5:360)*0.20, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.30, sind(0:5:360)*0.30, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.40, sind(0:5:360)*0.40, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
% line_width = 0.5;
% line_alpha = 0.30;
% plot(plot_data_x_axis, plot_data_y_axis, '-', 'Color', red_color);
% patch(plot_data_x_axis, plot_data_y_axis,'w',...
%         'linewidth',line_width,...
%         'EdgeAlpha',line_alpha,...
%         'EdgeColor',red_color);
line_width = 2.5;
line_alpha = 0.95;
plot_data_amp_mean = [corrSac_onset_prob_ALL_PCELL_mean corrSac_onset_prob_ALL_PCELL_mean(1), nan]';
plot_data_deg_mean = [0, 90, 180, 270, 0, nan]';
plot_data_x_axis_mean = plot_data_amp_mean .* cosd(plot_data_deg_mean);
plot_data_y_axis_mean = plot_data_amp_mean .* sind(plot_data_deg_mean);
% plot(plot_data_x_axis_mean, plot_data_y_axis_mean, '-', 'LineWidth', 3, 'Color', blue_color)
patch(plot_data_x_axis_mean, plot_data_y_axis_mean,'w',...
        'linewidth',line_width,...
        'EdgeAlpha',line_alpha,...
        'EdgeColor',red_color);
axis equal
xlim([-0.4 0.4])
ylim([-0.4 0.4])
xlabel('Corr On');
% ylabel('CS probability');
set(gca, 'XTick', -0.4:0.2:0.4, 'YTick', -0.4:0.2:0.4)

ESN_Beautify_Plot
figure_size  = [16.0 8.0];
paper_margin = [0.1 0.1];
paper_size = figure_size + 2 * paper_margin;
set(hFig, 'PaperSize', paper_size);
set(hFig, 'PaperPositionMode', 'manual');
set(hFig, 'PaperPosition', [paper_margin figure_size]);
set(hFig, 'Position', [[1 1] figure_size]);
set(hFig, 'PaperOrientation', 'portrait');

sgtitle(hFig, 'CS Tuning', 'Interpreter', 'none');

%% Save combined file

fprintf(['Saving plots ' ' ... ']);
save_file_path = './compress_pCells/ALL_PCELL_46_figs';
save_file_name = 'ALL_PCELL_46';
saveas(figure(1), [save_file_path filesep save_file_name '_SS_all'], 'pdf');
saveas(figure(1), [save_file_path filesep save_file_name '_SS_all'], 'png');
saveas(figure(2), [save_file_path filesep save_file_name '_CS_all'], 'pdf');
saveas(figure(2), [save_file_path filesep save_file_name '_CS_all'], 'png');
saveas(figure(3), [save_file_path filesep save_file_name '_SS_000'], 'pdf');
saveas(figure(3), [save_file_path filesep save_file_name '_SS_000'], 'png');
saveas(figure(4), [save_file_path filesep save_file_name '_SS_090'], 'pdf');
saveas(figure(4), [save_file_path filesep save_file_name '_SS_090'], 'png');
saveas(figure(5), [save_file_path filesep save_file_name '_SS_180'], 'pdf');
saveas(figure(5), [save_file_path filesep save_file_name '_SS_180'], 'png');
saveas(figure(6), [save_file_path filesep save_file_name '_SS_270'], 'pdf');
saveas(figure(6), [save_file_path filesep save_file_name '_SS_270'], 'png');
saveas(figure(7), [save_file_path filesep save_file_name '_CS_000'], 'pdf');
saveas(figure(7), [save_file_path filesep save_file_name '_CS_000'], 'png');
saveas(figure(8), [save_file_path filesep save_file_name '_CS_090'], 'pdf');
saveas(figure(8), [save_file_path filesep save_file_name '_CS_090'], 'png');
saveas(figure(9), [save_file_path filesep save_file_name '_CS_180'], 'pdf');
saveas(figure(9), [save_file_path filesep save_file_name '_CS_180'], 'png');
saveas(figure(10),[save_file_path filesep save_file_name '_CS_270'], 'pdf');
saveas(figure(10),[save_file_path filesep save_file_name '_CS_270'], 'png');
saveas(figure(11),[save_file_path filesep save_file_name '_CS_tun'], 'pdf');
saveas(figure(11),[save_file_path filesep save_file_name '_CS_tun'], 'png');

fprintf(' --> Completed. \n');

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
