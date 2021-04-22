%% function temp_ESN_plot_data_compress_all_pcell
function temp_ESN_plot_data_compress_all_pcell
%% load ALL_PCELL_COMPRESSED_DATA
clear;
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
            RASTER_DATA_ALL_PCELL_TUNED.(field_name_ALL_PCELL).(['numTrial' train_data_dir])(counter_ALL_PCELL, :) = ...
                ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).(field_name_ALL_PCELL).(['train_data_logic_CS' train_data_dir_tuned '_numTrial']);
            RASTER_DATA_ALL_PCELL_TUNED.(field_name_ALL_PCELL).(['inds_span' train_data_dir])(counter_ALL_PCELL, :) = ...
                ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).(field_name_ALL_PCELL).('inds_span');
            RASTER_DATA_ALL_PCELL_TUNED.(field_name_ALL_PCELL).(['SS_firing_rate' train_data_dir])(counter_ALL_PCELL, :) = ...
                ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).('Neural_Properties_data').('SS_firing_rate');
            RASTER_DATA_ALL_PCELL_TUNED.(field_name_ALL_PCELL).(['CS_firing_rate' train_data_dir])(counter_ALL_PCELL, :) = ...
                ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).('Neural_Properties_data').('CS_firing_rate');
        end
        RASTER_DATA_ALL_PCELL_TUNED.(field_name_ALL_PCELL).('inds_span')(counter_ALL_PCELL, :) = ...
            ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).(field_name_ALL_PCELL).('inds_span');
    end
end

train_data_names = {'train_data_logic_SS', 'train_data_logic_CS', 'velocity_data', 'inds_span', 'numTrial', 'SS_firing_rate', 'CS_firing_rate'};
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
        RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([    train_data_name '_ON']) = [ ...
            RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_right']) ; ...
            RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_top']) ; ...
            ];
        RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([    train_data_name '_OFF']) = [ ...
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
    RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([    train_data_name '_ON_numTrial']) = [ ...
        RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_right_numTrial']) ; ...
        RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_top_numTrial']) ; ...
        ];
    RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([    train_data_name '_OFF_numTrial']) = [ ...
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

%% Build RASTER_DATA_ALL_PCELL_TUNED_BUNDLE
clearvars -except ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL RASTER_DATA_ALL_PCELL_TUNED bundle_inds
clearvars RASTER_DATA_ALL_PCELL_TUNED_BUNDLE
length_bundle = length(bundle_inds);
train_data_names = {'train_data_logic_SS', 'train_data_logic_CS', 'velocity_data'};
train_data_dirs  = {'_right', '_top', '_left', '_down'};
field_names = fieldnames(RASTER_DATA_ALL_PCELL_TUNED);
length_span = size(RASTER_DATA_ALL_PCELL_TUNED.(field_names{1}).([train_data_names{1} train_data_dirs{1}]), 2);
for counter_field_names_RASTER_DATA_ALL_PCELL = 1 : 1 : length(field_names)
    field_name = field_names{counter_field_names_RASTER_DATA_ALL_PCELL};
    for counter_train_data_names = 1 : 1 : length(train_data_names)
        train_data_name = train_data_names{counter_train_data_names};
        variable_data_all_dir = zeros(length_bundle, length_span);
        for counter_train_data_dirs = 1 : 1 : length(train_data_dirs)
            train_data_dir = train_data_dirs{counter_train_data_dirs};
            
            variable_data = RASTER_DATA_ALL_PCELL_TUNED.(field_name).([train_data_name train_data_dir]);
            numTrial = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['train_data_logic_CS' train_data_dir '_numTrial']);
            variable_data_TUNED_BUNDLE = nan(length_bundle, length_span);
            for counter_bundle = 1 : length_bundle
                inds_ = bundle_inds{counter_bundle};
                variable_data_bundle = variable_data(inds_, :);
                numTrial_bundle = numTrial(inds_, :);
                variable_data_bundle_mean = (numTrial_bundle(:)' * variable_data_bundle) ./ sum(numTrial_bundle);
                variable_data_TUNED_BUNDLE(counter_bundle, :) = variable_data_bundle_mean;
            end % counter_bundle
            RASTER_DATA_ALL_PCELL_TUNED_BUNDLE.(field_name).([train_data_name train_data_dir]) = variable_data_TUNED_BUNDLE;
            variable_data_all_dir = variable_data_all_dir + variable_data_TUNED_BUNDLE;
        end % counter_train_data_dirs
        variable_data_all_dir = variable_data_all_dir ./ length(train_data_dirs);
        RASTER_DATA_ALL_PCELL_TUNED_BUNDLE.(field_name).([train_data_name '_all']) = variable_data_all_dir;
        
        variable_data_right = RASTER_DATA_ALL_PCELL_TUNED_BUNDLE.(field_name).([train_data_name '_right']);
        variable_data_top   = RASTER_DATA_ALL_PCELL_TUNED_BUNDLE.(field_name).([train_data_name '_top']);
        variable_data_left  = RASTER_DATA_ALL_PCELL_TUNED_BUNDLE.(field_name).([train_data_name '_left']);
        variable_data_down  = RASTER_DATA_ALL_PCELL_TUNED_BUNDLE.(field_name).([train_data_name '_down']);
        variable_data_ON    = (variable_data_right + variable_data_top) ./ 2;
        variable_data_OFF   = (variable_data_left + variable_data_down) ./ 2;
        RASTER_DATA_ALL_PCELL_TUNED_BUNDLE.(field_name).([train_data_name '_ON'])  = variable_data_ON;
        RASTER_DATA_ALL_PCELL_TUNED_BUNDLE.(field_name).([train_data_name '_OFF']) = variable_data_OFF;
        
    end % counter_train_data_names
    inds_span = RASTER_DATA_ALL_PCELL_TUNED.(field_name).inds_span;
    inds_span_TUNED_BUNDLE = nan(length_bundle, length_span);
    for counter_bundle = 1 : length_bundle
        inds_ = bundle_inds{counter_bundle};
        inds_span_bundle = inds_span(inds_, :);
        inds_span_bundle_mean = mean(inds_span_bundle, 1);
        inds_span_TUNED_BUNDLE(counter_bundle, :) = inds_span_bundle_mean;
    end % counter_bundle
    RASTER_DATA_ALL_PCELL_TUNED_BUNDLE.(field_name).inds_span = inds_span_TUNED_BUNDLE;
end % counter_field_names_RASTER_DATA_ALL_PCELL

%% Plot-01 to 07 RASTER_DATA_ALL_PCELL_TUNED_BUNDLE
clearvars -except ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL RASTER_DATA_ALL_PCELL_TUNED bundle_inds RASTER_DATA_ALL_PCELL_TUNED_BUNDLE
close all
range_SS = [30 100];
range_SS_population = [-30 50];
range_CS = [0 4];
range_vm = [0 500];
length_bundle = length(bundle_inds);
clearvars hFig hAx

train_data_dirs  = {'_right', '_top', '_left', '_down', '_all', '_ON', '_OFF'};
train_data_names = {'train_data_logic_SS', 'train_data_logic_SS', 'train_data_logic_CS', 'velocity_data'};
field_names = {'raster_data_cue_present', 'raster_data_primSac_onset', 'raster_data_primSac_offset', 'raster_data_corrSac_onset' };
yLabels = {'Change in Population (spk/s)', 'SS Firing (spk/s)', 'CS Firing (spk/s)', 'Eye velocity (deg/s)'};
xLabels = {'Cue Pres (ms)', 'Prim Sac On (ms)', 'Prim Sac Off (ms)', 'Corr Sac On (ms)'};
sgTitles = {'CS-000', 'CS-090', 'CS-180', 'CS-270', 'All Direction', 'CS-ON', 'CS-OFF'};
for counter_train_data_dirs = 1 : 1 : length(train_data_dirs)
    hFig(counter_train_data_dirs) = figure(counter_train_data_dirs);
    clf(hFig(counter_train_data_dirs));
    train_data_dir = train_data_dirs{counter_train_data_dirs};
    for counter_train_data_names = 1 : 1 : length(train_data_names)
        train_data_name = train_data_names{counter_train_data_names};
        for counter_field_names = 1 : 1 : length(field_names)
            field_name = field_names{counter_field_names};
            inds_span = mean(RASTER_DATA_ALL_PCELL_TUNED_BUNDLE.(field_name).inds_span, 1);
            if counter_train_data_names == 1
                numTrial = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['numTrial' train_data_dir]);
                SS_firing_rate = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['SS_firing_rate' train_data_dir]);
                train_data_logic_SS = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['train_data_logic_SS' train_data_dir]);
                train_data_logic_SS = train_data_logic_SS - (SS_firing_rate ./ 1000);
                
                num_perm = 5000;
                num_cells = size(train_data_logic_SS, 1);
                length_span = size(train_data_logic_SS, 2);
                train_data_logic_SS_perm = nan(num_perm, length_span);
                for counter_perm = 1 : num_perm
                    inds_perm = randi(num_cells,[num_cells 1]);
                    train_data_logic_SS_ = train_data_logic_SS(inds_perm, :);
                    numTrial_ = numTrial(inds_perm, :);
                    firing_SS_mean = numTrial_' * train_data_logic_SS_ ./ sum(numTrial_) .* 1000;
                    train_data_logic_SS_perm(counter_perm, :) = firing_SS_mean;
                end
                variable_data_mean = ESN_smooth( nanmean(train_data_logic_SS_perm, 1) );
                variable_data_stdv = ESN_smooth( nanstd(train_data_logic_SS_perm, 0, 1) );
            end
            if counter_train_data_names ~= 1
                variable_data = RASTER_DATA_ALL_PCELL_TUNED_BUNDLE.(field_name).([train_data_name train_data_dir]);
                variable_data_mean = nanmean(variable_data, 1);
                variable_data_stdv = nanstd(variable_data, 0, 1) ./ sqrt(length_bundle);
                if(counter_train_data_names ~= 4)
                    variable_data_mean = ESN_smooth(variable_data_mean * 1000);
                    variable_data_stdv = ESN_smooth(variable_data_stdv * 1000);

                end
            end
            
            subplot_num = 4*(counter_train_data_names-1)+counter_field_names;
            hAx{counter_train_data_dirs}(subplot_num) = subplot(4, 4, subplot_num);
            hAx_ = hAx{counter_train_data_dirs}(subplot_num);
            hold(hAx_, 'on')
            plot(inds_span, variable_data_mean, '-k', 'LineWidth', 1)
            plot(inds_span, variable_data_mean+variable_data_stdv, '-k', 'LineWidth', 0.5)
            plot(inds_span, variable_data_mean-variable_data_stdv, '-k', 'LineWidth', 0.5)
            if(counter_field_names == 1)
                ylabel(hAx_, yLabels{counter_train_data_names});
            end
            if(counter_train_data_names == 4)
                xlabel(hAx_, xLabels{counter_field_names});
            end
            if(counter_train_data_names == 1)
                ylim(range_SS_population);
            elseif(counter_train_data_names == 2)
                ylim(range_SS);
            elseif(counter_train_data_names == 3)
                ylim(range_CS);
            elseif(counter_train_data_names == 4)
                ylim(range_vm);
            end
        end
    end
    ESN_Beautify_Plot(hFig(counter_train_data_dirs), [11, 8.5])
    sgtitle(hFig(counter_train_data_dirs), sgTitles{counter_train_data_dirs}, 'Interpreter', 'none');
end

hFig(8) = figure(8);
clf(hFig(8));
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

length_bundle = length(bundle_inds);
cue_present_prob_TUNED_BUNDLE    = nan(length_bundle, 4);
primSac_onset_prob_TUNED_BUNDLE  = nan(length_bundle, 4);
primSac_offset_prob_TUNED_BUNDLE = nan(length_bundle, 4);
corrSac_onset_prob_TUNED_BUNDLE  = nan(length_bundle, 4);
for counter_bundle = 1 : length_bundle
    inds_ = bundle_inds{counter_bundle};
    cue_present_prob_TUNED_BUNDLE(counter_bundle, :) = ...
        numTrials_ALL_PCELL(inds_, :)' * cue_present_prob_ALL_PCELL(inds_, :) ./ sum(numTrials_ALL_PCELL(inds_, :));
    primSac_onset_prob_TUNED_BUNDLE(counter_bundle, :) = ...
        numTrials_ALL_PCELL(inds_, :)' * primSac_onset_prob_ALL_PCELL(inds_, :) ./ sum(numTrials_ALL_PCELL(inds_, :));
    primSac_offset_prob_TUNED_BUNDLE(counter_bundle, :) = ...
        numTrials_ALL_PCELL(inds_, :)' * primSac_offset_prob_ALL_PCELL(inds_, :) ./ sum(numTrials_ALL_PCELL(inds_, :));
    corrSac_onset_prob_TUNED_BUNDLE(counter_bundle, :) = ...
        numTrials_ALL_PCELL(inds_, :)' * corrSac_onset_prob_ALL_PCELL(inds_, :) ./ sum(numTrials_ALL_PCELL(inds_, :));
end % counter_bundle

overall_prob_TUNED_BUNDLE = (cue_present_prob_TUNED_BUNDLE + ...
    primSac_onset_prob_TUNED_BUNDLE + ...
    primSac_offset_prob_TUNED_BUNDLE + ...
    corrSac_onset_prob_TUNED_BUNDLE) ./ 4;

overall_prob_TUNED_BUNDLE_mean = nanmean(overall_prob_TUNED_BUNDLE, 1);
overall_prob_TUNED_BUNDLE_stdv = nanstd(overall_prob_TUNED_BUNDLE, 0, 1) ./ sqrt(length_bundle);
overall_prob_TUNED_BUNDLE_stdv_plus = overall_prob_TUNED_BUNDLE_mean + overall_prob_TUNED_BUNDLE_stdv;
overall_prob_TUNED_BUNDLE_stdv_minus = overall_prob_TUNED_BUNDLE_mean - overall_prob_TUNED_BUNDLE_stdv;

plot_data_amp_mean = [overall_prob_TUNED_BUNDLE_mean, overall_prob_TUNED_BUNDLE_mean(1), nan]';
plot_data_deg_mean = [0, 90, 180, 270, 0, nan]';
plot_data_x_axis_mean = plot_data_amp_mean .* cosd(plot_data_deg_mean);
plot_data_y_axis_mean = plot_data_amp_mean .* sind(plot_data_deg_mean);

plot_data_amp_stdv_p = [overall_prob_TUNED_BUNDLE_stdv_plus, overall_prob_TUNED_BUNDLE_stdv_plus(1), nan]';
plot_data_deg_stdv_p = [0, 90, 180, 270, 0, nan]';
plot_data_x_axis_stdv_p = plot_data_amp_stdv_p .* cosd(plot_data_deg_stdv_p);
plot_data_y_axis_stdv_p = plot_data_amp_stdv_p .* sind(plot_data_deg_stdv_p);

plot_data_amp_stdv_m = [overall_prob_TUNED_BUNDLE_stdv_minus, overall_prob_TUNED_BUNDLE_stdv_minus(1), nan]';
plot_data_deg_stdv_m = [0, 90, 180, 270, 0, nan]';
plot_data_x_axis_stdv_m = plot_data_amp_stdv_m .* cosd(plot_data_deg_stdv_m);
plot_data_y_axis_stdv_m = plot_data_amp_stdv_m .* sind(plot_data_deg_stdv_m);

hold on
plot([0.4, -0.4, nan, 0, -0, nan,]', [0, 0, nan, 0.4, -0.4, nan], 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.05, sind(0:5:360)*0.05, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
plot(cosd(0:5:360)*0.10, sind(0:5:360)*0.10, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
plot(cosd(0:5:360)*0.15, sind(0:5:360)*0.15, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
plot(cosd(0:5:360)*0.20, sind(0:5:360)*0.20, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
plot(cosd(0:5:360)*0.25, sind(0:5:360)*0.25, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
plot(cosd(0:5:360)*0.30, sind(0:5:360)*0.30, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])

plot(plot_data_x_axis_mean,   plot_data_y_axis_mean, '-k', 'LineWidth', 1)
plot(plot_data_x_axis_stdv_p, plot_data_y_axis_stdv_p, '-k', 'LineWidth', 0.5)
plot(plot_data_x_axis_stdv_m, plot_data_y_axis_stdv_m, '-k', 'LineWidth', 0.5)

axis equal
xlim([-0.3 0.3])
ylim([-0.3 0.3])
xlabel('CS probability');
ylabel('CS probability');
set(gca, 'XTick', -0.3:0.1:0.3, 'YTick', -0.3:0.1:0.3)


ESN_Beautify_Plot(hFig(8), [4, 3])
title('CS Tuning', 'Interpreter', 'none');


%{
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

%}

%% Save combined file

fprintf(['Saving plots ' ' ... ']);
save_file_path = './compress_pCells/ALL_PCELL_48_figs';
save_file_name = 'ALL_PCELL_48';
saveas(hFig(1), [save_file_path filesep save_file_name '_CS_000'], 'pdf');
saveas(hFig(1), [save_file_path filesep save_file_name '_CS_000'], 'png');
saveas(hFig(2), [save_file_path filesep save_file_name '_CS_090'], 'pdf');
saveas(hFig(2), [save_file_path filesep save_file_name '_CS_090'], 'png');
saveas(hFig(3), [save_file_path filesep save_file_name '_CS_180'], 'pdf');
saveas(hFig(3), [save_file_path filesep save_file_name '_CS_180'], 'png');
saveas(hFig(4), [save_file_path filesep save_file_name '_CS_270'], 'pdf');
saveas(hFig(4), [save_file_path filesep save_file_name '_CS_270'], 'png');
saveas(hFig(5), [save_file_path filesep save_file_name '_CS_all'], 'pdf');
saveas(hFig(5), [save_file_path filesep save_file_name '_CS_all'], 'png');
saveas(hFig(6), [save_file_path filesep save_file_name '_CS_ON' ], 'pdf');
saveas(hFig(6), [save_file_path filesep save_file_name '_CS_ON' ], 'png');
saveas(hFig(7), [save_file_path filesep save_file_name '_CS_OFF'], 'pdf');
saveas(hFig(7), [save_file_path filesep save_file_name '_CS_OFF'], 'png');
saveas(hFig(8), [save_file_path filesep save_file_name '_CS_tun'], 'pdf');
saveas(hFig(8), [save_file_path filesep save_file_name '_CS_tun'], 'png');

% close all
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
