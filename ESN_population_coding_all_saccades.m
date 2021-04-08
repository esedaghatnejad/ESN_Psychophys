%% function ESN_population_coding_all_saccades
function ESN_population_coding_all_saccades
%% build and save ALL_PCELL_COMPRESSED_DATA
path_data_monkey_sorted = uigetdir;

% build pCell_list, this is a hard coded cell with the id of all of the pCells and the bundles
pCell_list_1 = build_pCell_list_Mirza_pre201906();
pCell_list_2 = build_pCell_list_Mirza_post201906();
pCell_list_3 = build_pCell_list_Ramon();
pCell_list_4 = build_pCell_list_Mirza_post202011();
pCell_list = vertcat(pCell_list_1, pCell_list_2, pCell_list_3, pCell_list_4);
%%
% build ALL_PCELL_COMPRESSED_DATA, open plot data and put them together
ALL_PCELL_all_saccades = build_ALL_PCELL_all_saccades(pCell_list, path_data_monkey_sorted);

num_pCells = length(pCell_list);
ALL_PCELL_name = ['ALL_PCELL_' num2str(num_pCells)];
if ~exist([path_data_monkey_sorted filesep ALL_PCELL_name], 'dir')
    mkdir([path_data_monkey_sorted filesep ALL_PCELL_name]);
end

% save ALL_PCELL_COMPRESSED_DATA
save([path_data_monkey_sorted filesep ALL_PCELL_name filesep 'ALL_PCELL_all_saccades.mat'],'ALL_PCELL_all_saccades','-v7.3')
 
%% build _tuned
clearvars -except ALL_PCELL_all_saccades idx_burster idx_pauser PAIR_DATA_all_saccades
ALL_PCELL_all_saccades_tuned = build_ALL_PCELL_all_saccades_tuned(ALL_PCELL_all_saccades);
%
ALL_PCELL_all_saccades_tuned_amp        = avg_over_angle(     ALL_PCELL_all_saccades_tuned);
ALL_PCELL_all_saccades_tuned_ang        = avg_over_amplitude( ALL_PCELL_all_saccades_tuned);
ALL_PCELL_all_saccades_tuned_90_270_amp = combine_90_270(     ALL_PCELL_all_saccades_tuned);
ALL_PCELL_all_saccades_amp              = avg_over_angle(     ALL_PCELL_all_saccades);
ALL_PCELL_all_saccades_ang              = avg_over_amplitude( ALL_PCELL_all_saccades);
ALL_PCELL_all_saccades_90_270_amp       = combine_90_270(     ALL_PCELL_all_saccades);
%
ALL_PCELL_all_saccades_tuned_ang_amp        = avg_over_angle(     ALL_PCELL_all_saccades_tuned_ang);
ALL_PCELL_all_saccades_tuned_amp_ang        = avg_over_amplitude( ALL_PCELL_all_saccades_tuned_amp);
ALL_PCELL_all_saccades_tuned_90_270_amp_ang = avg_over_amplitude( ALL_PCELL_all_saccades_tuned_90_270_amp);
ALL_PCELL_all_saccades_ang_amp              = avg_over_angle(     ALL_PCELL_all_saccades_ang);
ALL_PCELL_all_saccades_amp_ang              = avg_over_amplitude( ALL_PCELL_all_saccades_amp);
ALL_PCELL_all_saccades_90_270_amp_ang       = avg_over_amplitude( ALL_PCELL_all_saccades_90_270_amp);
%

%%
%%
% pCell_list = build_pCell_list_Mirza_post202011();
% re_run_ESN_monkey_behavior_all_saccades(pCell_list, path_data_monkey_sorted);

%% Build and save PAIR_DATA_all_saccades
path_data_monkey_sorted = uigetdir;
pair_list_full = build_pair_list_full();

PAIR_DATA_all_saccades = build_PAIR_all_saccades(pair_list_full, path_data_monkey_sorted);

num_pair = size(pair_list_full,1)/2;
PAIR_DATA_name = ['PAIR_DATA_' num2str(num_pair)];
if ~exist([path_data_monkey_sorted filesep PAIR_DATA_name], 'dir')
    mkdir([path_data_monkey_sorted filesep PAIR_DATA_name]);
end

% save ALL_PCELL_COMPRESSED_DATA
save([path_data_monkey_sorted filesep PAIR_DATA_name filesep 'PAIR_DATA_all_saccades.mat'],'PAIR_DATA_all_saccades','-v7.3')

%% build pair_data_tuned
clearvars -except ALL_PCELL_all_saccades idx_burster idx_pauser PAIR_DATA_all_saccades
PAIR_DATA_all_saccades_tuned = build_PAIR_DATA_all_saccades_tuned(PAIR_DATA_all_saccades);
%
PAIR_DATA_all_saccades_tuned_amp        = pair_data_avg_over_angle(     PAIR_DATA_all_saccades_tuned);
PAIR_DATA_all_saccades_tuned_ang        = pair_data_avg_over_amplitude( PAIR_DATA_all_saccades_tuned);
PAIR_DATA_all_saccades_amp              = pair_data_avg_over_angle(     PAIR_DATA_all_saccades);
PAIR_DATA_all_saccades_ang              = pair_data_avg_over_amplitude( PAIR_DATA_all_saccades);
%
PAIR_DATA_all_saccades_tuned_ang_amp        = pair_data_avg_over_angle(     PAIR_DATA_all_saccades_tuned_ang);
PAIR_DATA_all_saccades_tuned_amp_ang        = pair_data_avg_over_amplitude( PAIR_DATA_all_saccades_tuned_amp);
PAIR_DATA_all_saccades_ang_amp              = pair_data_avg_over_angle(     PAIR_DATA_all_saccades_ang);
PAIR_DATA_all_saccades_amp_ang              = pair_data_avg_over_amplitude( PAIR_DATA_all_saccades_amp);

end

%% function build_ALL_PCELL_all_saccades
function ALL_PCELL_all_saccades = build_ALL_PCELL_all_saccades(pCell_list, path_data_monkey_sorted)
%% init vars
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);
variable_list = {'primSac', 'corrSac', 'toTgt', 'fromTgt', 'lowAmp', 'toTgtStr', 'fromCenter', 'aroundCenter', 'nonTask', 'all', 'primCorrStr', 'tgtCuePres', 'tgtPrimSac', 'tgtCorrSac'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};
clearvars ALL_PCELL_COMPRESSED_DATA
%% Loop over pCells
for counter_pCell = 1 : 1 : num_pCells
    fprintf(['Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ... ']);
    num_recording = sum(pCell_list_isstr(counter_pCell, :));
    clearvars TRAIN_DATA_recordings
    for counter_recording = 1 : 1 : num_recording
        %% build plot_data address
        file_name = pCell_list{counter_pCell, counter_recording};
        year_ = file_name(1:2);
        month_ = file_name(3:4);
        day_ = file_name(5:6);
        hour_ = file_name(8:9);
        minute_ = file_name(10:11);
        second_ = file_name(12:13);
        subFolder_month = ['20' year_ '-' month_ filesep];
        subFolder_day = ['20' year_ '-' month_ '-' day_ filesep];
        subFolder_recording = ['20' year_ '-' month_ '-' day_ '_' hour_ '-' minute_ '-' second_ filesep];
        subFolder_data = ['analyzed_data' filesep];
        file_path = [path_data_monkey_sorted subFolder_month subFolder_day subFolder_recording subFolder_data];
        %% load plot_data
        TRAIN_DATA_recording = build_TRAIN_DATA(file_path, file_name);
        TRAIN_DATA_recording.file_name = file_name;
        TRAIN_DATA_recording.file_path = file_path;
        if TRAIN_DATA_recording.file_name(18) == 's'
            TRAIN_DATA_recording.id          = TRAIN_DATA_recording.file_name(1:16);
        elseif TRAIN_DATA_recording.file_name(18) == '2'
            TRAIN_DATA_recording.id          = TRAIN_DATA_recording.file_name(1:18);
        else
            error('Build plot_data_compress: cell id does not follow the standards')
        end
        TRAIN_DATA_recordings(counter_recording) = TRAIN_DATA_recording;
    end
    %% 
    TRAIN_DATA_cell = struct;
    if num_recording > 1
        
        for counter_variable  = 1 : length(variable_list)
        for counter_indType   = 1 : length(indType_list)
        for counter_spikeType = 1 : length(spikeType_list)
            variable_name = variable_list{counter_variable};
            indType_name = indType_list{counter_indType};
            spikeType_name = spikeType_list{counter_spikeType};
            num_amp_bin = size(TRAIN_DATA_recordings(1).(variable_name).([spikeType_name '_train_' indType_name]),1);
            num_ang_bin = size(TRAIN_DATA_recordings(1).(variable_name).([spikeType_name '_train_' indType_name]),2);

            for counter_amp = 1 : num_amp_bin
            for counter_ang = 1 : num_ang_bin
                train_data_recordings = zeros(size(TRAIN_DATA_recordings(1).(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang}));
                velocity_data_recordings = zeros(size(TRAIN_DATA_recordings(1).(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang}));
                num_saccades_recordings = 0;
                for counter_recording = 1 : 1 : num_recording
                    train_data = TRAIN_DATA_recordings(counter_recording).(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang};
                    velocity_data = TRAIN_DATA_recordings(counter_recording).(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang};
                    num_saccades = TRAIN_DATA_recordings(counter_recording).(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang);
                    train_data_recordings = train_data_recordings + (train_data * num_saccades);
                    velocity_data_recordings = velocity_data_recordings + (velocity_data * num_saccades);
                    num_saccades_recordings = num_saccades_recordings + num_saccades;
                end
                if num_saccades_recordings == 0
                    TRAIN_DATA_cell.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data_recordings;
                    TRAIN_DATA_cell.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data_recordings;
                    TRAIN_DATA_cell.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades_recordings;
                else
                    TRAIN_DATA_cell.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data_recordings ./ num_saccades_recordings;
                    TRAIN_DATA_cell.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data_recordings ./ num_saccades_recordings;
                    TRAIN_DATA_cell.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades_recordings;
                end
                
            end
            end
        end
        end
        end
        % concatenate cell id
        TRAIN_DATA_cell.id = cell(num_recording, 1);
        TRAIN_DATA_cell.file_name = cell(num_recording, 1);
        TRAIN_DATA_cell.file_path = cell(num_recording, 1);
        for counter_recording = 1 : 1 : num_recording
            TRAIN_DATA_cell.id{counter_recording, 1} = TRAIN_DATA_recordings(counter_recording).id;
            TRAIN_DATA_cell.file_name{counter_recording, 1} = TRAIN_DATA_recordings(counter_recording).file_name;
            TRAIN_DATA_cell.file_path{counter_recording, 1} = TRAIN_DATA_recordings(counter_recording).file_path;
        end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        TRAIN_DATA_cell.Neural_Properties.SS_num = 0;
        TRAIN_DATA_cell.Neural_Properties.SS_duration = 0;
        TRAIN_DATA_cell.Neural_Properties.SS_firing_rate = 0;
        TRAIN_DATA_cell.Neural_Properties.SS_time = [];
        TRAIN_DATA_cell.Neural_Properties.SS_waveform = [];
        TRAIN_DATA_cell.Neural_Properties.CS_num = 0;
        TRAIN_DATA_cell.Neural_Properties.CS_firing_rate = 0;
        TRAIN_DATA_cell.Neural_Properties.CS_time = [];
        TRAIN_DATA_cell.Neural_Properties.CS_waveform = [];
        for counter_recording = 1 : 1 : num_recording
            TRAIN_DATA_cell.Neural_Properties.SS_num = TRAIN_DATA_cell.Neural_Properties.SS_num + ...
                TRAIN_DATA_recordings(counter_recording).Neural_Properties.SS_num;
            TRAIN_DATA_cell.Neural_Properties.SS_duration = TRAIN_DATA_cell.Neural_Properties.SS_duration + ...
                TRAIN_DATA_recordings(counter_recording).Neural_Properties.SS_duration;
            TRAIN_DATA_cell.Neural_Properties.SS_time = vertcat(TRAIN_DATA_cell.Neural_Properties.SS_time ,...
                TRAIN_DATA_recordings(counter_recording).Neural_Properties.SS_time);
            TRAIN_DATA_cell.Neural_Properties.SS_waveform = vertcat(TRAIN_DATA_cell.Neural_Properties.SS_waveform, ...
                TRAIN_DATA_recordings(counter_recording).Neural_Properties.SS_waveform);

            TRAIN_DATA_cell.Neural_Properties.CS_num = TRAIN_DATA_cell.Neural_Properties.CS_num + ...
                TRAIN_DATA_recordings(counter_recording).Neural_Properties.CS_num;
            TRAIN_DATA_cell.Neural_Properties.CS_time = vertcat(TRAIN_DATA_cell.Neural_Properties.CS_time, ...
                TRAIN_DATA_recordings(counter_recording).Neural_Properties.CS_time);
            TRAIN_DATA_cell.Neural_Properties.CS_waveform = vertcat(TRAIN_DATA_cell.Neural_Properties.CS_waveform, ...
                TRAIN_DATA_recordings(counter_recording).Neural_Properties.CS_waveform);
            % concatenate Neural_Properties.Corr_data
            if ~isfield(TRAIN_DATA_cell.Neural_Properties, 'Corr_data')
                TRAIN_DATA_cell.Neural_Properties.Corr_data = struct();
            end
            field_names_Corr_data = fieldnames(TRAIN_DATA_recordings(counter_recording).Neural_Properties.Corr_data);
            for counter_field_names_Corr_data = 1 : 1 : length(field_names_Corr_data)
                field_name_Corr_data = field_names_Corr_data{counter_field_names_Corr_data};
                data_field_name_Corr_data = TRAIN_DATA_recordings(counter_recording).Neural_Properties.Corr_data.(field_name_Corr_data);
                if ~isfield(TRAIN_DATA_cell.Neural_Properties.Corr_data, field_name_Corr_data)
                    TRAIN_DATA_cell.Neural_Properties.Corr_data.(field_name_Corr_data) = [];
                end
                data_field_name_cell = TRAIN_DATA_cell.Neural_Properties.Corr_data.(field_name_Corr_data);
                data_field_name_cell = vertcat(data_field_name_cell, data_field_name_Corr_data);
                TRAIN_DATA_cell.Neural_Properties.Corr_data.(field_name_Corr_data) = data_field_name_cell;
            end
        end
        TRAIN_DATA_cell.Neural_Properties.SS_firing_rate = ...
                TRAIN_DATA_cell.Neural_Properties.SS_num / TRAIN_DATA_cell.Neural_Properties.SS_duration;
        TRAIN_DATA_cell.Neural_Properties.CS_firing_rate = ...
                TRAIN_DATA_cell.Neural_Properties.CS_num / TRAIN_DATA_cell.Neural_Properties.SS_duration;
        TRAIN_DATA_cell.Neural_Properties.SS_waveform = nanmean(TRAIN_DATA_cell.Neural_Properties.SS_waveform, 1);
        TRAIN_DATA_cell.Neural_Properties.CS_waveform = nanmean(TRAIN_DATA_cell.Neural_Properties.CS_waveform, 1);
        field_names_Corr_data = fieldnames(TRAIN_DATA_cell.Neural_Properties.Corr_data);
        for counter_field_names_Corr_data = 1 : 1 : length(field_names_Corr_data)
            field_name_Corr_data = field_names_Corr_data{counter_field_names_Corr_data};
            data_field_name_Corr_data = TRAIN_DATA_cell.Neural_Properties.Corr_data.(field_name_Corr_data);
            TRAIN_DATA_cell.Neural_Properties.Corr_data.(field_name_Corr_data) = nanmean(data_field_name_Corr_data, 1);
        end
    else
        TRAIN_DATA_cell = TRAIN_DATA_recordings;
        TRAIN_DATA_cell.Neural_Properties.SS_waveform = nanmean(TRAIN_DATA_cell.Neural_Properties.SS_waveform, 1);
        TRAIN_DATA_cell.Neural_Properties.CS_waveform = nanmean(TRAIN_DATA_cell.Neural_Properties.CS_waveform, 1);
        field_names_Corr_data = fieldnames(TRAIN_DATA_cell.Neural_Properties.Corr_data);
        for counter_field_names_Corr_data = 1 : 1 : length(field_names_Corr_data)
            field_name_Corr_data = field_names_Corr_data{counter_field_names_Corr_data};
            data_field_name_Corr_data = TRAIN_DATA_cell.Neural_Properties.Corr_data.(field_name_Corr_data);
            TRAIN_DATA_cell.Neural_Properties.Corr_data.(field_name_Corr_data) = nanmean(data_field_name_Corr_data, 1);
        end
    end
    
    %% Save TRAIN_DATA_cell into ALL_PCELL_COMPRESSED_DATA
    ALL_PCELL_COMPRESSED_DATA(counter_pCell) = TRAIN_DATA_cell;
    fprintf(' --> Completed. \n')
end

%% Loop over pCells
fprintf(['Building ALL_PCELL_all_saccades', ' ... ']);
clearvars ALL_PCELL_all_saccades
for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    num_amp_bin = size(ALL_PCELL_COMPRESSED_DATA(1).(variable_name).([spikeType_name '_train_' indType_name]),1);
    num_ang_bin = size(ALL_PCELL_COMPRESSED_DATA(1).(variable_name).([spikeType_name '_train_' indType_name]),2);

    for counter_amp = 1 : num_amp_bin
    for counter_ang = 1 : num_ang_bin
        for counter_pCell = 1 : 1 : num_pCells
            train_data = ALL_PCELL_COMPRESSED_DATA(counter_pCell).(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang};
            velocity_data = ALL_PCELL_COMPRESSED_DATA(counter_pCell).(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang};
            num_saccades = ALL_PCELL_COMPRESSED_DATA(counter_pCell).(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang);
            ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang}(counter_pCell,:) = train_data;
            ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang}(counter_pCell,:) = velocity_data;
            ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, counter_ang}(counter_pCell,:) = num_saccades;
        end
    end
    end
end
end
end
% Neural_Properties
for counter_pCell = 1 : 1 : num_pCells
    ALL_PCELL_all_saccades.Neural_Properties.SS_num(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.SS_num;
    ALL_PCELL_all_saccades.Neural_Properties.SS_duration(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.SS_duration;
    ALL_PCELL_all_saccades.Neural_Properties.SS_firing_rate(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.SS_firing_rate;
%     ALL_PCELL_all_saccades.Neural_Properties.SS_time(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.SS_time;
    ALL_PCELL_all_saccades.Neural_Properties.SS_waveform(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.SS_waveform;
    ALL_PCELL_all_saccades.Neural_Properties.CS_num(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.CS_num;
    ALL_PCELL_all_saccades.Neural_Properties.CS_firing_rate(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.CS_firing_rate;
%     ALL_PCELL_all_saccades.Neural_Properties.CS_time(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.CS_time;
    ALL_PCELL_all_saccades.Neural_Properties.CS_waveform(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.CS_waveform;
    field_names_Corr_data = fieldnames(ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.Corr_data);
    for counter_field_names_Corr_data = 1 : 1 : length(field_names_Corr_data)
        field_name_Corr_data = field_names_Corr_data{counter_field_names_Corr_data};
        data_field_name_Corr_data = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.Corr_data.(field_name_Corr_data);
        ALL_PCELL_all_saccades.Neural_Properties.Corr_data.(field_name_Corr_data)(counter_pCell,:) = data_field_name_Corr_data;
    end
end
fprintf(' --> Completed. \n')
end

%% function avg_over_amplitude
function ALL_PCELL_all_saccades_ang = avg_over_amplitude(ALL_PCELL_all_saccades)
%% Loop over pCells
fprintf(['Building ALL_PCELL_all_saccades_ang', ' ... ']);
clearvars ALL_PCELL_all_saccades_ang
variable_list = {'primSac', 'corrSac', 'toTgt', 'fromTgt', 'lowAmp', 'toTgtStr', 'fromCenter', 'aroundCenter', 'nonTask', 'all', 'primCorrStr', 'tgtCuePres', 'tgtPrimSac', 'tgtCorrSac'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};
num_pCells = size(ALL_PCELL_all_saccades.primSac.SS_train_start{1,1}, 1);
span_width = size(ALL_PCELL_all_saccades.primSac.SS_train_start{1,1}, 2);
for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    num_amp_bin = size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]),1);
    num_ang_bin = size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]),2);
    
    for counter_ang = 1 : num_ang_bin
        train_data_all_ang = zeros(size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]){1, counter_ang}));
        velocity_data_all_ang = zeros(size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_velocity_' indType_name]){1, counter_ang}));
        num_saccades_all_ang = zeros(size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){1, counter_ang}));
        if size(num_saccades_all_ang, 2) == 1
            num_saccades_all_ang = repmat(num_saccades_all_ang, 1, span_width);
        end
        for counter_amp = 1 : num_amp_bin
            train_data = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang};
            velocity_data = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang};
            num_saccades = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, counter_ang};
            if size(num_saccades, 2) == 1
                num_saccades = repmat(num_saccades, 1, span_width);
            end
            train_data_all_ang = train_data_all_ang + (train_data .* num_saccades);
            velocity_data_all_ang = velocity_data_all_ang + (velocity_data .* num_saccades);
            num_saccades_all_ang = num_saccades_all_ang + num_saccades;
        end
        train_data_all_ang = train_data_all_ang ./ num_saccades_all_ang;
        train_data_all_ang(isnan(train_data_all_ang)) = 0;
        velocity_data_all_ang = velocity_data_all_ang ./ num_saccades_all_ang;
        velocity_data_all_ang(isnan(velocity_data_all_ang)) = 0;
        ALL_PCELL_all_saccades_ang.(variable_name).([spikeType_name '_train_' indType_name]){1, counter_ang} = train_data_all_ang;
        ALL_PCELL_all_saccades_ang.(variable_name).([spikeType_name '_velocity_' indType_name]){1, counter_ang} = velocity_data_all_ang;
        ALL_PCELL_all_saccades_ang.(variable_name).([spikeType_name '_num_sac_' indType_name]){1, counter_ang} = num_saccades_all_ang;
    end
end
end
end
fprintf(' --> Completed. \n')
end

%% function avg_over_angle
function ALL_PCELL_all_saccades_amp = avg_over_angle(ALL_PCELL_all_saccades)
%% Loop over pCells
fprintf(['Building ALL_PCELL_all_saccades_amp', ' ... ']);
clearvars ALL_PCELL_all_saccades_amp
variable_list = {'primSac', 'corrSac', 'toTgt', 'fromTgt', 'lowAmp', 'toTgtStr', 'fromCenter', 'aroundCenter', 'nonTask', 'all', 'primCorrStr', 'tgtCuePres', 'tgtPrimSac', 'tgtCorrSac'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};
num_pCells = size(ALL_PCELL_all_saccades.primSac.SS_train_start{1,1}, 1);
span_width = size(ALL_PCELL_all_saccades.primSac.SS_train_start{1,1}, 2);
for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    num_amp_bin = size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]),1);
    num_ang_bin = size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]),2);
    
    for counter_amp = 1 : num_amp_bin
        train_data_all_ang = zeros(size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, 1}));
        velocity_data_all_ang = zeros(size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, 1}));
        num_saccades_all_ang = zeros(size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, 1}));
        if size(num_saccades_all_ang, 2) == 1
            num_saccades_all_ang = repmat(num_saccades_all_ang, 1, span_width);
        end
        for counter_ang = 1 : num_ang_bin
            train_data = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang};
            velocity_data = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang};
            num_saccades = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, counter_ang};
            if size(num_saccades, 2) == 1
                num_saccades = repmat(num_saccades, 1, span_width);
            end
            train_data_all_ang = train_data_all_ang + (train_data .* num_saccades);
            velocity_data_all_ang = velocity_data_all_ang + (velocity_data .* num_saccades);
            num_saccades_all_ang = num_saccades_all_ang + num_saccades;
        end
        train_data_all_ang = train_data_all_ang ./ num_saccades_all_ang;
        train_data_all_ang(isnan(train_data_all_ang)) = 0;
        velocity_data_all_ang = velocity_data_all_ang ./ num_saccades_all_ang;
        velocity_data_all_ang(isnan(velocity_data_all_ang)) = 0;
        ALL_PCELL_all_saccades_amp.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, 1} = train_data_all_ang;
        ALL_PCELL_all_saccades_amp.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, 1} = velocity_data_all_ang;
        ALL_PCELL_all_saccades_amp.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, 1} = num_saccades_all_ang;
    end
end
end
end
fprintf(' --> Completed. \n')
end

%% function build_ALL_PCELL_all_saccades_tuned
function ALL_PCELL_all_saccades_tuned = build_ALL_PCELL_all_saccades_tuned(ALL_PCELL_all_saccades)
%% find CS-on for each cell
ALL_PCELL_all_saccades_ang = avg_over_amplitude(ALL_PCELL_all_saccades);
num_ang_bin = size(ALL_PCELL_all_saccades_ang.tgtCuePres.CS_train_start,2);
num_pCells = size(ALL_PCELL_all_saccades_ang.tgtCuePres.CS_train_start{1,1},1);
%%
% METHOD 1
%{
CS_prob_tgtCuePres = zeros(num_pCells, num_ang_bin);
CS_prob_tgtPrimSac = zeros(num_pCells, num_ang_bin);
CS_prob_tgtEndPres = zeros(num_pCells, num_ang_bin);
CS_prob_tgtCorrSac = zeros(num_pCells, num_ang_bin);
for counter_ang = 1 : num_ang_bin
    CS_prob_tgtCuePres(:,counter_ang) = nansum(ALL_PCELL_all_saccades_ang.tgtCuePres.CS_train_start{1,counter_ang}(:,300:500), 2);
    CS_prob_tgtPrimSac(:,counter_ang) = nansum(ALL_PCELL_all_saccades_ang.tgtPrimSac.CS_train_start{1,counter_ang}(:,100:300), 2);
    CS_prob_tgtEndPres(:,counter_ang) = nansum(ALL_PCELL_all_saccades_ang.tgtPrimSac.CS_train_finish{1,counter_ang}(:,300:500), 2);
    CS_prob_tgtCorrSac(:,counter_ang) = nansum(ALL_PCELL_all_saccades_ang.tgtCorrSac.CS_train_start{1,counter_ang}(:,100:300), 2);
end
CS_prob_tgtCuePres(1:51, [2 4 6 8]) = nan; % Mirza 4-dir cells
CS_prob_tgtPrimSac(1:51, [2 4 6 8]) = nan; % Mirza 4-dir cells
CS_prob_tgtEndPres(1:51, [2 4 6 8]) = nan; % Mirza 4-dir cells
CS_prob_tgtCorrSac(1:51, [2 4 6 8]) = nan; % Mirza 4-dir cells
CS_prob_avg = (CS_prob_tgtCuePres + CS_prob_tgtPrimSac + CS_prob_tgtEndPres + CS_prob_tgtCorrSac) ./ 4;
[~, idx_CS_max] = nanmax(CS_prob_avg, [], 2);
hFig = figure(2); clf(hFig); hold('on');
plot(nanmean(CS_prob_tgtCuePres))
plot(nanmean(CS_prob_tgtPrimSac))
plot(nanmean(CS_prob_tgtEndPres))
plot(nanmean(CS_prob_tgtCorrSac))
plot(nanmean(CS_prob_avg), 'k')
ylim([0.14 0.24])
%}
% METHOD 2
%{
CS_prob_primSac = zeros(num_pCells, num_ang_bin);
CS_prob_corrSac = zeros(num_pCells, num_ang_bin);
for counter_ang = 1 : num_ang_bin
    CS_prob_primSac(:,counter_ang) = nansum(ALL_PCELL_all_saccades_ang.primSac.CS_train_start{1,counter_ang}(:,100:300), 2);
    CS_prob_corrSac(:,counter_ang) = nansum(ALL_PCELL_all_saccades_ang.corrSac.CS_train_start{1,counter_ang}(:,100:300), 2);
end
CS_prob_primSac(1:51, [2 4 6 8]) = nan;
CS_prob_corrSac(1:51, [2 4 6 8]) = nan;
CS_prob_avg = (CS_prob_primSac + CS_prob_corrSac) ./ 2;
[~, idx_CS_max] = nanmax(CS_prob_avg, [], 2);
hFig = figure(1); clf(hFig); hold('on');
plot(nanmean(CS_prob_primSac))
plot(nanmean(CS_prob_corrSac))
plot(nanmean(CS_prob_avg), 'k')
ylim([0.14 0.24])
%}
% METHOD 3
%
CS_prob_primSac = zeros(num_pCells, num_ang_bin);
CS_prob_corrSac = zeros(num_pCells, num_ang_bin);
for counter_ang = 1 : num_ang_bin
    CS_prob_primSac(:,counter_ang) = nansum(ALL_PCELL_all_saccades_ang.primSac.CS_train_start{1,counter_ang}(:,100:300), 2);
    CS_prob_corrSac(:,counter_ang) = nansum(ALL_PCELL_all_saccades_ang.corrSac.CS_train_start{1,counter_ang}(:,100:300), 2);
end

CS_prob_primSac(1:51, [2 4 6 8]) = nan;
CS_prob_corrSac(1:51, [2 4 6 8]) = nan;
CS_prob_avg = (CS_prob_primSac + CS_prob_corrSac) ./ 2;
% [~, idx_CS_max] = nanmax(CS_prob_avg, [], 2);
dir_angles = deg2rad([180, 225, 270, 315, 0, 45, 90, 135]);
% compute weighted sum of cos and sin of angles
r = nansum(CS_prob_avg.* repmat(exp(1i*dir_angles), num_pCells, 1) , 2);
% Computes the mean direction for circular data.
CS_ang_avg = rad2deg(angle(r));
% sum of weights
CS_prob_sum = nansum(CS_prob_avg,2);
% Computes mean resultant vector length for circular data.
CS_rho_avg = abs(r) ./ CS_prob_sum;

ang_edges = -202.5 : 45 : +202.5;
ang_edges_last_bin_id = length(ang_edges) - 1;
[~, ~, idx_CS_max] = histcounts(CS_ang_avg, ang_edges);
idx_CS_max(idx_CS_max == ang_edges_last_bin_id) = 1;

hFig = figure(1); clf(hFig);
step_size_ = pi/8;
ang_edges = -pi-(step_size_/2):step_size_:pi-(step_size_/2);
polarhistogram(deg2rad(CS_ang_avg), ang_edges, 'DisplayStyle', 'bar','FaceColor','red', 'EdgeColor', 'red')
set(gca, 'ThetaTick', 0:45:315, 'RTick', 0:5:25)
rlim([0 20])
%}
%%
idx_CS_tuning = zeros(num_pCells, num_ang_bin);
for counter_pCell = 1 : 1 : num_pCells
    CS_on_index = idx_CS_max(counter_pCell);
    CS_on_index = CS_on_index - 1; % CS_on_index should be in 0-index format
    if (CS_on_index == 8); CS_on_index = 0; end
    idx_CS_tuning(counter_pCell, :) = mod((CS_on_index : 1 : CS_on_index+7), 8) + 1;
end

%% build ALL_PCELL_all_saccades_tuned
fprintf(['Building ALL_PCELL_all_saccades_tuned', ' ... ']);
variable_list = {'primSac', 'corrSac', 'toTgt', 'fromTgt', 'lowAmp', 'toTgtStr', 'fromCenter', 'aroundCenter', 'nonTask', 'all', 'primCorrStr', 'tgtCuePres', 'tgtPrimSac', 'tgtCorrSac'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};
clearvars ALL_PCELL_all_saccades_tuned
for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    num_amp_bin = size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]),1);
    num_ang_bin = size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]),2);

    for counter_amp = 1 : num_amp_bin
    for counter_ang = 1 : num_ang_bin
        for counter_pCell = 1 : 1 : num_pCells
            idx_ang_ = idx_CS_tuning(counter_pCell, counter_ang);
            
            train_data = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, idx_ang_}(counter_pCell,:);
            velocity_data = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, idx_ang_}(counter_pCell,:);
            num_saccades = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, idx_ang_}(counter_pCell,:);
            
            ALL_PCELL_all_saccades_tuned.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang}(counter_pCell,:) = train_data;
            ALL_PCELL_all_saccades_tuned.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang}(counter_pCell,:) = velocity_data;
            ALL_PCELL_all_saccades_tuned.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, counter_ang}(counter_pCell,:) = num_saccades;
        end
    end
    end
end
end
end
ALL_PCELL_all_saccades_tuned.idx_CS_tuning = idx_CS_tuning;
fprintf(' --> Completed. \n')

end

%% function build_TRAIN_DATA
function [TRAIN_DATA] = build_TRAIN_DATA(file_path, file_name)
%% Handle inputs
if nargin < 1
    [file_name_,file_path] = uigetfile([pwd filesep '*.psort'], 'Select psort file');
    [~,file_name,~] = fileparts(file_name_);
end

%% load EPHYS sorted DATA
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
file_name = [file_name '.psort'];
fprintf(['Loading ', file_name, ' ... ']);
% EPHYS.CH_sorted = load([file_path filesep file_name], 'CS_data', 'SS_data');
DATA_PSORT = Psort_read_psort([file_path file_name]);
EPHYS.CH_sorted_file_name = file_name;
EPHYS.CH_sorted_file_path = file_path;
fprintf(' --> Completed. \n')

%% load EPHYS EVENT DATA
file_name = [EPHYS.CH_sorted_file_name(1:13) '_EVE1_aligned.mat'];
% fprintf(['Loading ', file_name, ' ... ']);
EPHYS.CH_EVE = load([file_path file_name]);
if isfield(EPHYS.CH_EVE, 'EPHYS_time_15K')
    EPHYS.CH_EVE.EPHYS_time_30K = EPHYS.CH_EVE.EPHYS_time_15K(:);
else
    EPHYS.CH_EVE.EPHYS_time_30K = EPHYS.CH_EVE.EPHYS_time_30K(:);
end
EPHYS.CH_EVE.EPHYS_time_1K  = reshape(EPHYS.CH_EVE.EPHYS_time_1K ,[], 1);
EPHYS.CH_EVE.BEHAVE_time_1K = reshape(EPHYS.CH_EVE.BEHAVE_time_1K,[], 1);
% fprintf(' --> Completed. \n')

%% load BEHAVE DATA
file_name = [EPHYS.CH_sorted_file_name(1:13) '_ANALYZED.mat'];
% fprintf(['Loading ', file_name, ' ... ']);
BEHAVE = load([file_path file_name]);
% fprintf(' --> Completed. \n')

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

%% build SSxSS AUTO PROBABILITY
clearvars -except EPHYS BEHAVE TRAIN_DATA
% fprintf(['Building SSxSS_AUTO & CSxSS_AUTO PROBABILITY ' ' ...'])
SS_time   = EPHYS.CH_sorted.SS_data.SS_time;
CS_time   = EPHYS.CH_sorted.CS_data.CS_time;
Corr_data = ESN_correlogram(SS_time, CS_time);
EPHYS.CH_sorted.Corr_data.CS_inds_span     = Corr_data.CS_inds_span;
EPHYS.CH_sorted.Corr_data.CS_bin_size_time = Corr_data.CS_bin_size_time;
EPHYS.CH_sorted.Corr_data.SS_inds_span     = Corr_data.SS_inds_span;
EPHYS.CH_sorted.Corr_data.SS_bin_size_time = Corr_data.SS_bin_size_time;
EPHYS.CH_sorted.Corr_data.SS_SSxSS_AUTO    = Corr_data.SS_SSxSS_AUTO;
EPHYS.CH_sorted.Corr_data.CS_CSxSS_AUTO    = Corr_data.CS_CSxSS_AUTO;
% fprintf(' --> Completed. \n')

%% SS & CS train_aligned
clearvars -except EPHYS BEHAVE TRAIN_DATA
% fprintf(['Building CS & SS train_aligned', ' ... ']);
EPHYS_time_1K     = EPHYS.CH_EVE.EPHYS_time_1K; % EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_time_1K; % 
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

% fprintf(' --> Completed. \n')

%% Building TRAIN_DATA
clearvars -except EPHYS BEHAVE TRAIN_DATA
fprintf(['Building TRAIN_DATA', ' ... ']);
TRAIN_DATA = struct;
inds_span    = ((-300+1) : 1 : (300))';
amp_edges = [-0.5 2 4 6 8 10 50];
vel_edges = [0 200 300 400 500 600 10000];
ang_edges = -202.5 : 45 : +202.5;
TRAIN_DATA.inds_span = inds_span;
TRAIN_DATA.amp_edges = amp_edges;
TRAIN_DATA.vel_edges = vel_edges;
TRAIN_DATA.ang_edges = ang_edges;
TRAIN_DATA.use_vel_instead_of_amp = false;
length_train_data_ = length(EPHYS.CH_EVE.EPHYS_time_1K);
length_velocity_data_ = length(BEHAVE.aligned.time_1K);
velocity_trace_data = BEHAVE.aligned.eye_r_vm_filt;

ang_edges_last_bin_id = length(ang_edges) - 1;

[~, ~, amp_bin] = histcounts(BEHAVE.SACS_ALL_DATA.eye_r_amp_m, amp_edges);
[~, ~, vel_bin] = histcounts(BEHAVE.SACS_ALL_DATA.eye_r_vm_max, vel_edges);
[~, ~, ang_bin] = histcounts(BEHAVE.SACS_ALL_DATA.eye_r_ang, ang_edges);
ang_bin(ang_bin == ang_edges_last_bin_id) = 1;

num_amp_bin = length(amp_edges) - 1;
num_vel_bin = length(vel_edges) - 1;
num_ang_bin = length(ang_edges) - 2;

if TRAIN_DATA.use_vel_instead_of_amp
    amp_bin = vel_bin;
    num_amp_bin = num_vel_bin;
end

BEHAVE.SACS_ALL_DATA.is_all = BEHAVE.SACS_ALL_DATA.validity;
BEHAVE.SACS_ALL_DATA.is_primCorrStr = BEHAVE.SACS_ALL_DATA.is_primSac | BEHAVE.SACS_ALL_DATA.is_corrSac | BEHAVE.SACS_ALL_DATA.is_toTgtStr;

variable_list = {'primSac', 'corrSac', 'toTgt', 'fromTgt', 'lowAmp', 'toTgtStr', 'fromCenter', 'aroundCenter', 'nonTask', 'all', 'primCorrStr'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};


for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    
    TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]) = cell(num_amp_bin, num_ang_bin);
    TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]) = cell(num_amp_bin, num_ang_bin);
    TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name]) = zeros(num_amp_bin, num_ang_bin);
    is_variable_name = BEHAVE.SACS_ALL_DATA.(['is_' variable_name]);
    spike_train_data = EPHYS.CH_EVE.(['EPHYS_' spikeType_name '_train_1K']);
    ind_indType_name = BEHAVE.SACS_ALL_DATA.(['ind_' indType_name]);
    ind_converted = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(ind_indType_name);
    
    for counter_amp = 1 : num_amp_bin
    for counter_ang = 1 : num_ang_bin
        is_ang_bin   = (ang_bin == counter_ang);
        is_amp_bin   = (amp_bin == counter_amp);
        ind_train    = ind_converted(is_variable_name & is_ang_bin & is_amp_bin);
        ind_velocity = ind_indType_name(is_variable_name & is_ang_bin & is_amp_bin);
        if isempty(ind_train)
            train_data = zeros(1, length(inds_span));
            velocity_data = zeros(1, length(inds_span));
            num_saccades = 0;
            TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
            continue;
        end
        inds_train = repmat( ind_train(:), 1, length(inds_span)) + repmat(inds_span(:)', length(ind_train), 1);
        inds_train( inds_train < 1 ) = 1;
        inds_train( inds_train > length_train_data_ ) = length_train_data_;
        inds_velocity = repmat( ind_velocity(:), 1, length(inds_span)) + repmat(inds_span(:)', length(ind_velocity), 1);
        inds_velocity( inds_velocity < 1 ) = 1;
        inds_velocity( inds_velocity > length_velocity_data_ ) = length_velocity_data_;
        if length(ind_train) == 1
            train_data = reshape( spike_train_data(inds_train) , 1, []);
            velocity_data = reshape( velocity_trace_data(inds_velocity) , 1, []);
            num_saccades = length(ind_train);
            TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
            continue;
        end
        train_data = nanmean( spike_train_data(inds_train) );
        velocity_data = nanmean( velocity_trace_data(inds_velocity) );
        num_saccades = length(ind_train);
        TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
        TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
        TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
    end
    end
end
end
end

%% Add tgtCuePres
clearvars -except EPHYS BEHAVE TRAIN_DATA
inds_span    = TRAIN_DATA.inds_span;
amp_edges = TRAIN_DATA.amp_edges;
vel_edges = TRAIN_DATA.vel_edges;
ang_edges = TRAIN_DATA.ang_edges;
length_train_data_ = length(EPHYS.CH_EVE.EPHYS_time_1K);
length_velocity_data_ = length(BEHAVE.aligned.time_1K);
velocity_trace_data = BEHAVE.aligned.eye_r_vm_filt;
ang_edges_last_bin_id = length(ang_edges) - 1;

range_trials = BEHAVE.aligned.BEHAVE_range_trials;
tgt_amp_m = sqrt( ...
    (BEHAVE.TRIALS_DATA.cue_y(range_trials) - BEHAVE.TRIALS_DATA.start_y(range_trials)).^2 + ...
    (BEHAVE.TRIALS_DATA.cue_x(range_trials) - BEHAVE.TRIALS_DATA.start_x(range_trials)).^2 ...
    );
tgt_ang = atan2d( ...
    (BEHAVE.TRIALS_DATA.cue_y(range_trials) - BEHAVE.TRIALS_DATA.start_y(range_trials)) , ...
    (BEHAVE.TRIALS_DATA.cue_x(range_trials) - BEHAVE.TRIALS_DATA.start_x(range_trials)) ...
    );

[~, ~, amp_bin] = histcounts(tgt_amp_m, amp_edges);
[~, ~, vel_bin] = histcounts(BEHAVE.SACS_PRIM_DATA.eye_r_vm_max(range_trials), vel_edges);
[~, ~, ang_bin] = histcounts(tgt_ang  , ang_edges);
ang_bin(ang_bin == ang_edges_last_bin_id) = 1;

num_amp_bin = length(amp_edges) - 1;
num_vel_bin = length(vel_edges) - 1;
num_ang_bin = length(ang_edges) - 2;

if TRAIN_DATA.use_vel_instead_of_amp
    amp_bin = vel_bin;
    num_amp_bin = num_vel_bin;
end

variable_list = {'tgtCuePres'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};

% TRAIN_DATA = struct;
for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    
    TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]) = cell(num_amp_bin, num_ang_bin);
    TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]) = cell(num_amp_bin, num_ang_bin);
    TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name]) = zeros(num_amp_bin, num_ang_bin);
    is_variable_name = reshape(BEHAVE.aligned.BEHAVE_validity_prim, 1, []);
    ind_indType_name = BEHAVE.aligned.BEHAVE_ind_cue_present;
    spike_train_data = EPHYS.CH_EVE.(['EPHYS_' spikeType_name '_train_1K']);
%     ind_converted = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(ind_indType_name);
    ind_converted = EPHYS.CH_EVE.align_photodiode.EPHYS_PD_aligned_ind_1K(ind_indType_name);
    
    for counter_amp = 1 : num_amp_bin
    for counter_ang = 1 : num_ang_bin
        is_ang_bin   = (ang_bin == counter_ang);
        is_amp_bin   = (amp_bin == counter_amp);
        ind_train    = ind_converted(is_variable_name & is_ang_bin & is_amp_bin);
        ind_velocity = ind_indType_name(is_variable_name & is_ang_bin & is_amp_bin);
        if isempty(ind_train)
            train_data = zeros(1, length(inds_span));
            velocity_data = zeros(1, length(inds_span));
            num_saccades = 0;
            TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
            continue;
        end
        inds_train = repmat( ind_train(:), 1, length(inds_span)) + repmat(inds_span(:)', length(ind_train), 1);
        inds_train( inds_train < 1 ) = 1;
        inds_train( inds_train > length_train_data_ ) = length_train_data_;
        inds_velocity = repmat( ind_velocity(:), 1, length(inds_span)) + repmat(inds_span(:)', length(ind_velocity), 1);
        inds_velocity( inds_velocity < 1 ) = 1;
        inds_velocity( inds_velocity > length_velocity_data_ ) = length_velocity_data_;
        if length(ind_train) == 1
            train_data = reshape( spike_train_data(inds_train) , 1, []);
            velocity_data = reshape( velocity_trace_data(inds_velocity) , 1, []);
            num_saccades = length(ind_train);
            TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
            continue;
        end
        train_data = nanmean( spike_train_data(inds_train) );
        velocity_data = nanmean( velocity_trace_data(inds_velocity) );
        num_saccades = length(ind_train);
        TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
        TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
        TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
    end
    end
end
end
end

%% Add tgtPrimSac tgtCorrSac
clearvars -except EPHYS BEHAVE TRAIN_DATA
inds_span    = TRAIN_DATA.inds_span;
amp_edges = TRAIN_DATA.amp_edges;
vel_edges = TRAIN_DATA.vel_edges;
ang_edges = TRAIN_DATA.ang_edges;
length_train_data_ = length(EPHYS.CH_EVE.EPHYS_time_1K);
length_velocity_data_ = length(BEHAVE.aligned.time_1K);
velocity_trace_data = BEHAVE.aligned.eye_r_vm_filt;
ang_edges_last_bin_id = length(ang_edges) - 1;

range_trials = BEHAVE.aligned.BEHAVE_range_trials;
tgtPrimSac_amp_m = sqrt( ...
    (BEHAVE.TRIALS_DATA.cue_y(range_trials) - BEHAVE.TRIALS_DATA.start_y(range_trials)).^2 + ...
    (BEHAVE.TRIALS_DATA.cue_x(range_trials) - BEHAVE.TRIALS_DATA.start_x(range_trials)).^2 ...
    );
tgtPrimSac_ang = atan2d( ...
    (BEHAVE.TRIALS_DATA.cue_y(range_trials) - BEHAVE.TRIALS_DATA.start_y(range_trials)) , ...
    (BEHAVE.TRIALS_DATA.cue_x(range_trials) - BEHAVE.TRIALS_DATA.start_x(range_trials)) ...
    );

tgtCorrSac_amp_m = sqrt( ...
    (BEHAVE.TRIALS_DATA.end_y(range_trials) - BEHAVE.TRIALS_DATA.cue_y(range_trials)).^2 + ...
    (BEHAVE.TRIALS_DATA.end_x(range_trials) - BEHAVE.TRIALS_DATA.cue_x(range_trials)).^2 ...
    );
tgtCorrSac_ang = atan2d( ...
    (BEHAVE.TRIALS_DATA.end_y(range_trials) - BEHAVE.TRIALS_DATA.cue_y(range_trials)) , ...
    (BEHAVE.TRIALS_DATA.end_x(range_trials) - BEHAVE.TRIALS_DATA.cue_x(range_trials)) ...
    );

[~, ~, tgtPrimSac_amp_m_bin] = histcounts(tgtPrimSac_amp_m, amp_edges);
[~, ~, tgtPrimSac_vel_m_bin] = histcounts(BEHAVE.SACS_PRIM_DATA.eye_r_vm_max(range_trials), vel_edges);
[~, ~, tgtPrimSac_ang_bin] = histcounts(tgtPrimSac_ang  , ang_edges);
tgtPrimSac_ang_bin(tgtPrimSac_ang_bin == ang_edges_last_bin_id) = 1;

[~, ~, tgtCorrSac_amp_m_bin] = histcounts(tgtCorrSac_amp_m, amp_edges);
[~, ~, tgtCorrSac_vel_m_bin] = histcounts(BEHAVE.SACS_CORR_DATA.eye_r_vm_max(range_trials), vel_edges);
[~, ~, tgtCorrSac_ang_bin] = histcounts(tgtCorrSac_ang  , ang_edges);
tgtCorrSac_ang_bin(tgtCorrSac_ang_bin == ang_edges_last_bin_id) = 1;

num_amp_bin = length(amp_edges) - 1;
num_vel_bin = length(vel_edges) - 1;
num_ang_bin = length(ang_edges) - 2;

if TRAIN_DATA.use_vel_instead_of_amp
    tgtPrimSac_amp_m_bin = tgtPrimSac_vel_m_bin;
    tgtCorrSac_amp_m_bin = tgtCorrSac_vel_m_bin;
    num_amp_bin = num_vel_bin;
end

bin_data.tgtPrimSac_amp_m_bin = tgtPrimSac_amp_m_bin;
bin_data.tgtPrimSac_ang_bin   = tgtPrimSac_ang_bin;
bin_data.tgtCorrSac_amp_m_bin = tgtCorrSac_amp_m_bin;
bin_data.tgtCorrSac_ang_bin   = tgtCorrSac_ang_bin;

variable_list = {'tgtPrimSac', 'tgtCorrSac'};
variable_list2 = {'primSac', 'corrSac'};
indType_list = {'start', 'vmax', 'finish'};
indType_list2 = {'onset', 'vmax', 'offset'};
spikeType_list = {'SS', 'CS'};

BEHAVE.aligned.BEHAVE_validity_primSac = BEHAVE.aligned.BEHAVE_validity_prim;
BEHAVE.aligned.BEHAVE_validity_corrSac = BEHAVE.aligned.BEHAVE_validity_corr;% & BEHAVE.aligned.BEHAVE_validity_prim;
% TRAIN_DATA = struct;
for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    variable_name2 = variable_list2{counter_variable};
    indType_name = indType_list{counter_indType};
    indType_name2 = indType_list2{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    
    TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]) = cell(num_amp_bin, num_ang_bin);
    TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]) = cell(num_amp_bin, num_ang_bin);
    TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name]) = zeros(num_amp_bin, num_ang_bin);
    is_variable_name = reshape(BEHAVE.aligned.(['BEHAVE_validity_' variable_name2]), 1, []);
    ind_indType_name = BEHAVE.aligned.(['BEHAVE_ind_' variable_name2 '_' indType_name2]);
    spike_train_data = EPHYS.CH_EVE.(['EPHYS_' spikeType_name '_train_1K']);
    ind_converted = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(ind_indType_name);
    
    for counter_amp = 1 : num_amp_bin
    for counter_ang = 1 : num_ang_bin
        is_ang_bin   = (bin_data.([variable_name '_ang_bin'])  == counter_ang );
        is_amp_bin   = (bin_data.([variable_name '_amp_m_bin'])  == counter_amp );
        if contains(variable_name,'tgtPrimSac') && contains(indType_name,'finish')
            % use tgtCorrSac angle for primSac_finish
            is_ang_bin   = (bin_data.(['tgtCorrSac' '_ang_bin'])  == counter_ang );
            is_amp_bin   = (bin_data.(['tgtPrimSac' '_amp_m_bin'])  == counter_amp );
        end
        ind_train    = ind_converted(is_variable_name & is_ang_bin & is_amp_bin);
        ind_velocity = ind_indType_name(is_variable_name & is_ang_bin & is_amp_bin);
        if isempty(ind_train)
            train_data = zeros(1, length(inds_span));
            velocity_data = zeros(1, length(inds_span));
            num_saccades = 0;
            TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
            continue;
        end
        inds_train = repmat( ind_train(:), 1, length(inds_span)) + repmat(inds_span(:)', length(ind_train), 1);
        inds_train( inds_train < 1 ) = 1;
        inds_train( inds_train > length_train_data_ ) = length_train_data_;
        inds_velocity = repmat( ind_velocity(:), 1, length(inds_span)) + repmat(inds_span(:)', length(ind_velocity), 1);
        inds_velocity( inds_velocity < 1 ) = 1;
        inds_velocity( inds_velocity > length_velocity_data_ ) = length_velocity_data_;
        if length(ind_train) == 1
            train_data = reshape( spike_train_data(inds_train) , 1, []);
            velocity_data = reshape( velocity_trace_data(inds_velocity) , 1, []);
            num_saccades = length(ind_train);
            TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
            continue;
        end
        train_data = nanmean( spike_train_data(inds_train) );
        velocity_data = nanmean( velocity_trace_data(inds_velocity) );
        num_saccades = length(ind_train);
        TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
        TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
        TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
    end
    end
end
end
end

%% Add Neural_Properties
clearvars -except EPHYS BEHAVE TRAIN_DATA
TRAIN_DATA.Neural_Properties = struct();
if ~isempty(EPHYS.CH_sorted.SS_data.SS_time)
    TRAIN_DATA.Neural_Properties.SS_num = length(EPHYS.CH_sorted.SS_data.SS_time);
    TRAIN_DATA.Neural_Properties.SS_duration = (EPHYS.CH_sorted.SS_data.SS_time(end)) - (EPHYS.CH_sorted.SS_data.SS_time(1));
    TRAIN_DATA.Neural_Properties.SS_firing_rate = TRAIN_DATA.Neural_Properties.SS_num / TRAIN_DATA.Neural_Properties.SS_duration;
    TRAIN_DATA.Neural_Properties.SS_time = EPHYS.CH_sorted.SS_data.SS_time;
    TRAIN_DATA.Neural_Properties.SS_waveform = EPHYS.CH_sorted.SS_data.SS_waveform;
else
    TRAIN_DATA.Neural_Properties.SS_num = 0;
    TRAIN_DATA.Neural_Properties.SS_duration = (EPHYS.CH_sorted.CS_data.CS_time(end)) - (EPHYS.CH_sorted.CS_data.CS_time(1));
    TRAIN_DATA.Neural_Properties.SS_firing_rate = 0;
    TRAIN_DATA.Neural_Properties.SS_time = [];
    TRAIN_DATA.Neural_Properties.SS_waveform = [];
end

if ~isempty(EPHYS.CH_sorted.CS_data.CS_time)
    TRAIN_DATA.Neural_Properties.CS_num = length(EPHYS.CH_sorted.CS_data.CS_time);
    TRAIN_DATA.Neural_Properties.CS_firing_rate = TRAIN_DATA.Neural_Properties.CS_num / TRAIN_DATA.Neural_Properties.SS_duration;
    TRAIN_DATA.Neural_Properties.CS_time = EPHYS.CH_sorted.CS_data.CS_time;
    TRAIN_DATA.Neural_Properties.CS_waveform = EPHYS.CH_sorted.CS_data.CS_waveform;
else
    TRAIN_DATA.Neural_Properties.CS_num = 0;
    TRAIN_DATA.Neural_Properties.CS_firing_rate = 0;
    TRAIN_DATA.Neural_Properties.CS_time = [];
    TRAIN_DATA.Neural_Properties.CS_waveform = [];
end
TRAIN_DATA.Neural_Properties.Corr_data = EPHYS.CH_sorted.Corr_data;

%% Remove some fields
rmfields_list = {'inds_span', 'amp_edges', 'ang_edges', 'vel_edges', 'use_vel_instead_of_amp'};
TRAIN_DATA = rmfield(TRAIN_DATA,rmfields_list);
fprintf(' --> Completed. \n')

end

%% function ESN_smooth
% function smooth_data_ = ESN_smooth(data_, dim)
% % smooth data using 2nd order Savitzky-Golay filter with 21 points
% % if data_ is a matrix, the method will smooth each column by default or smooth along dim.
% % method = 'moving';  % Moving average. A lowpass filter with filter coefficients equal to the reciprocal of the span.
% % method = 'lowess';  % Local regression using weighted linear least squares and a 1st degree polynomial model.
% % method = 'loess';   % Local regression using weighted linear least squares and a 2nd degree polynomial model.
% % method = 'sgolay';  % Savitzky-Golay filter. A generalized moving average with filter coefficients determined by an unweighted linear least-squares regression and a polynomial model of specified degree (default is 2). The method can accept nonuniform predictor data.
% % method = 'rlowess'; % A robust version of 'lowess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
% % method = 'rloess';  % A robust version of 'loess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
% % smooth_data_ = smooth(data_, method);
% if nargin < 2
%     dim = 1;
% end
% if size(data_, 1) == 1
%     smooth_data_ = reshape(smooth(data_, 21, 'sgolay', 2), 1, []);
% elseif size(data_, 2) == 1
%     smooth_data_ = reshape(smooth(data_, 21, 'sgolay', 2), [], 1);
% else
%     smooth_data_ = nan(size(data_));
%     if dim == 1
%         % smooth columns
%         for counter_col = 1 : size(data_, 2)
%             smooth_data_(:, counter_col) = reshape(smooth(data_(:, counter_col), 21, 'sgolay', 2), [], 1);
%         end
%     elseif dim == 2
%         % smooth rows
%         for counter_row = 1 : size(data_, 1)
%             smooth_data_(counter_row, :) = reshape(smooth(data_(counter_row, :), 21, 'sgolay', 2), 1, []);
%         end
%     end
%     
% end
% end

%% function ESN_correlogram
% function Corr_data = ESN_correlogram(SS_time, CS_time)
% bin_size_time = 1e-3; % seconds
% span_window_size = (1 / bin_size_time) * (100 / 1000);
% span_window_size_half = round(span_window_size / 2);
% inds_span = ((-span_window_size_half+1) : 1 : (span_window_size_half))';
% 
% if (~isempty(CS_time)) && (~isempty(SS_time))
%     ch_time_min = min([SS_time(1) CS_time(1)]);
%     ch_time_min = max([(ch_time_min-2.0) 0]);
%     ch_time_max = max([SS_time(end) CS_time(end)]) + 2.0;
%     
%     CH__.SS_data.SS_time =  SS_time - ch_time_min;
%     CH__.CS_data.CS_time =  CS_time - ch_time_min;
%     CH__.SS_data.SS_time = ESN_Round(CH__.SS_data.SS_time, bin_size_time);
%     CH__.CS_data.CS_time = ESN_Round(CH__.CS_data.CS_time, bin_size_time);
%     CH__.SS_data.SS_ind  = round(CH__.SS_data.SS_time .* (1/bin_size_time));
%     CH__.SS_data.SS_ind( CH__.SS_data.SS_ind < 1 ) = 1;
%     CH__.CS_data.CS_ind  = round(CH__.CS_data.CS_time .* (1/bin_size_time));
%     CH__.CS_data.CS_ind( CH__.CS_data.CS_ind < 1 ) = 1;
%     CH__.CH_data.ch_time_reconstruct = ch_time_min : bin_size_time : ch_time_max;
% elseif (~isempty(CS_time))
%     ch_time_min = min(  CS_time(1) );
%     ch_time_min = max([(ch_time_min-2.0) 0]);
%     ch_time_max = max(  CS_time(end) ) + 2.0;
%     
%     CH__.CS_data.CS_time =  CS_time - ch_time_min;
%     CH__.CS_data.CS_time = ESN_Round(CH__.CS_data.CS_time, bin_size_time);
%     CH__.CS_data.CS_ind  = round(CH__.CS_data.CS_time .* (1/bin_size_time));
%     CH__.CS_data.CS_ind( CH__.CS_data.CS_ind < 1 ) = 1;
%     CH__.CH_data.ch_time_reconstruct = ch_time_min : bin_size_time : ch_time_max;
% elseif (~isempty(SS_time))
%     ch_time_min = min( SS_time(1)  );
%     ch_time_min = max([(ch_time_min-2.0) 0]);
%     ch_time_max = max( SS_time(end)  ) + 2.0;
%     
%     CH__.SS_data.SS_time =  SS_time - ch_time_min;
%     CH__.SS_data.SS_time = ESN_Round(CH__.SS_data.SS_time, bin_size_time);
%     CH__.SS_data.SS_ind  = round(CH__.SS_data.SS_time .* (1/bin_size_time));
%     CH__.SS_data.SS_ind( CH__.SS_data.SS_ind < 1 ) = 1;
%     CH__.CH_data.ch_time_reconstruct = ch_time_min : bin_size_time : ch_time_max;
% end
% 
% % SSxSS_AUTO
% if (~isempty(SS_time))
%     CH__.SS_data.SS_inds_reconstruct = repmat( CH__.SS_data.SS_ind(:), 1, length(inds_span)) + repmat(inds_span(:)', length(CH__.SS_data.SS_ind), 1);
%     CH__.SS_data.SS_inds_reconstruct( CH__.SS_data.SS_inds_reconstruct < 1 ) = 1;
%     CH__.SS_data.SS_inds_reconstruct( CH__.SS_data.SS_inds_reconstruct > length( CH__.CH_data.ch_time_reconstruct ) ) = length( CH__.CH_data.ch_time_reconstruct );
%     
%     CH__.SS_data.SS_event_trace = false( size(CH__.CH_data.ch_time_reconstruct) );
%     CH__.SS_data.SS_event_trace( CH__.SS_data.SS_ind ) = true ;
%     CH__.SS_data.SS_event_trace( 1   ) = false;
%     CH__.SS_data.SS_event_trace( end ) = false;
%     
%     CH__.SS_data.SS_event_reconstruct = CH__.SS_data.SS_event_trace( CH__.SS_data.SS_inds_reconstruct );
%     % SSxSS correlogram
%     SSxSS_AUTO       = CH__.SS_data.SS_event_reconstruct;
%     ss_inds_span     = repmat(inds_span(:)',     size(SS_time(:),1), 1);
%     ss_bin_size_time = repmat(bin_size_time(:)', size(SS_time(:),1), 1);
% else
%     SSxSS_AUTO       = false(0, length(inds_span(:)'));
%     ss_inds_span     = nan(0, length(inds_span(:)'));
%     ss_bin_size_time = nan(0, 1);
% end
% 
% % CSxSS_WITHIN
% if (~isempty(CS_time)) && (~isempty(SS_time))
%     CH__.CS_data.CS_inds_reconstruct = repmat( CH__.CS_data.CS_ind(:), 1, length(inds_span)) + repmat(inds_span(:)', length(CH__.CS_data.CS_ind), 1);
%     CH__.CS_data.CS_inds_reconstruct( CH__.CS_data.CS_inds_reconstruct < 1 ) = 1;
%     CH__.CS_data.CS_inds_reconstruct( CH__.CS_data.CS_inds_reconstruct > length( CH__.CH_data.ch_time_reconstruct ) ) = length( CH__.CH_data.ch_time_reconstruct );
%     
%     CH__.SS_data.SS_event_trace = false( size(CH__.CH_data.ch_time_reconstruct) );
%     CH__.SS_data.SS_event_trace( CH__.SS_data.SS_ind ) = true ;
%     CH__.SS_data.SS_event_trace( 1   ) = false;
%     CH__.SS_data.SS_event_trace( end ) = false;
%     
%     CH__.SS_data.SS_event_reconstruct = CH__.SS_data.SS_event_trace( CH__.CS_data.CS_inds_reconstruct );
%     % CSxSS correlogram
%     CSxSS_AUTO       = CH__.SS_data.SS_event_reconstruct;
%     cs_inds_span     = repmat(inds_span(:)',     size(CS_time(:),1), 1);
%     cs_bin_size_time = repmat(bin_size_time(:)', size(CS_time(:),1), 1);
% else
%     CSxSS_AUTO       = false(0, length(inds_span(:)'));
%     cs_inds_span     = nan(0, length(inds_span(:)'));
%     cs_bin_size_time = nan(0, 1);
% end
% 
% Corr_data = struct;
% Corr_data.CS_inds_span     = cs_inds_span;
% Corr_data.CS_bin_size_time = cs_bin_size_time;
% Corr_data.SS_inds_span     = ss_inds_span;
% Corr_data.SS_bin_size_time = ss_bin_size_time;
% Corr_data.SS_SSxSS_AUTO    = SSxSS_AUTO;
% Corr_data.CS_CSxSS_AUTO    = CSxSS_AUTO;
% end

%% function RE-RUN ESN_monkey_behavior_all_saccades
function re_run_ESN_monkey_behavior_all_saccades(pCell_list, path_data_monkey_sorted)
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);
%% Loop over pCells
for counter_pCell = 1 : 1 : num_pCells
    fprintf(['Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ... ']);
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
        %% RE-RUN ESN_monkey_behavior_all_saccades
        file_path = [path_data_monkey_sorted subFolder_month subFolder_day subFolder_recording subFolder_figs];
        file_name = [file_name_cell(1:13) '_ANALYZED.mat'];
        ESN_monkey_behavior_all_saccades(file_path, file_name);
    end
end

end

%% function scratch_plot1
function scratch_plot1
%%
hFig = figure(1);
clf(hFig)
subplot(1, 2, 1)
hold on
variable_name1 =  'tgtPrimSac';
variable_name2 =  'CS_train_start';
variable_name2V = 'CS_velocity_start';
variable_name3 =  'CS_num_sac_start';
ind_ang = 1;
range_ = 001:600;
num_pCells = 110;

ALL_PCELL_DATA_1 = ALL_PCELL_all_saccades_tuned_ang;

data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){1, ind_ang}*1000;% - ALL_PCELL_DATA_1.(variable_name1).(variable_name2){1, 5}*1000;
data_ = data_(:, range_);
plot(ESN_smooth(nanmean(data_)+(nanstd(data_)./sqrt(num_pCells))), 'k', 'LineWidth', 1)
plot(ESN_smooth(nanmean(data_)-(nanstd(data_)./sqrt(num_pCells))), 'k', 'LineWidth', 1)
plot(ESN_smooth(nanmean(data_)), 'k', 'LineWidth', 2)
% set(gca, 'XTick', 0:50:400)
% ylim([40 90])
ylim([0 3])

subplot(1, 2, 2)
hold on
data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){1, ind_ang};% - ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){1, 5};
data_ = data_(:, range_);
plot((nanmean(data_)+(nanstd(data_)./sqrt(num_pCells))), 'k', 'LineWidth', 1)
plot((nanmean(data_)-(nanstd(data_)./sqrt(num_pCells))), 'k', 'LineWidth', 1)
plot((nanmean(data_)), 'k', 'LineWidth', 2)
% set(gca, 'XTick', 0:50:400)
ylim([0 600])

ESN_Beautify_Plot
end

%% function scratch_plot2
function scratch_plot2
%% 
hFig = figure(2);
clf(hFig)
subplot(1, 2, 1)
hold on
variable_name1 =  'toTgtStr';
variable_name2 =  'SS_train_start';
variable_name2V = 'SS_velocity_start';
variable_name3 =  'SS_num_sac_start';
ind_ang = 5;
range_ = 201:400;

ALL_PCELL_DATA_1 = ALL_PCELL_all_saccades_tuned;
ALL_PCELL_DATA_2 = ALL_PCELL_all_saccades_tuned_ang;

data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){2, ind_ang}(:,range_)*1000;
data_((ALL_PCELL_DATA_1.(variable_name1).(variable_name3){2, ind_ang}(:,1)==0), :) = nan;
plot(ESN_smooth(nanmean(data_)))
data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){3, ind_ang}(:,range_)*1000;
data_((ALL_PCELL_DATA_1.(variable_name1).(variable_name3){3, ind_ang}(:,1)==0), :) = nan;
plot(ESN_smooth(nanmean(data_)))
data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){4, ind_ang}(:,range_)*1000;
data_((ALL_PCELL_DATA_1.(variable_name1).(variable_name3){4, ind_ang}(:,1)==0), :) = nan;
plot(ESN_smooth(nanmean(data_)))
% data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){5, ind_ang}(:,range_)*1000;
% data_((ALL_PCELL_DATA_1.(variable_name1).(variable_name3){5, ind_ang}==0), :) = nan;
% plot(ESN_smooth(nanmean(data_)))
data_ = ALL_PCELL_DATA_2.(variable_name1).(variable_name2){1, ind_ang}(:,range_)*1000;
plot(ESN_smooth(nanmean(data_)), 'k', 'LineWidth', 2)
ylim([50 80])

subplot(1, 2, 2)
hold on
data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){2, ind_ang}(:,range_);
data_((ALL_PCELL_DATA_1.(variable_name1).(variable_name3){2, ind_ang}(:,1)==0), :) = nan;
plot((nanmean(data_)))
data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){3, ind_ang}(:,range_);
data_((ALL_PCELL_DATA_1.(variable_name1).(variable_name3){3, ind_ang}(:,1)==0), :) = nan;
plot((nanmean(data_)))
data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){4, ind_ang}(:,range_);
data_((ALL_PCELL_DATA_1.(variable_name1).(variable_name3){4, ind_ang}(:,1)==0), :) = nan;
plot((nanmean(data_)))
% data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){5, ind_ang}(:,range_);
% data_((ALL_PCELL_DATA_1.(variable_name1).(variable_name3){5, ind_ang}==0), :) = nan;
% plot(ESN_smooth(nanmean(data_)))
data_ = ALL_PCELL_DATA_2.(variable_name1).(variable_name2V){1, ind_ang}(:,range_);
plot((nanmean(data_)), 'k', 'LineWidth', 2)
ylim([0 600])

ESN_Beautify_Plot
end

%% function scratch_plot3
function scratch_plot3
%% 
hFig = figure(3);
clf(hFig)
subplot(1, 2, 1)
hold on
variable_name1 =  'toTgtStr';
variable_name2 =  'SS_train_start';
variable_name2V = 'SS_velocity_start';
variable_name3 =  'SS_num_sac_start';
ind_ang_1 = 5;
ind_ang_2 = 1;
range_ = 201:400;

ALL_PCELL_DATA_1 = ALL_PCELL_all_saccades_tuned;
ALL_PCELL_DATA_2 = ALL_PCELL_all_saccades_tuned_ang;

data_1 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){2, ind_ang_1}(:,range_);
data_1(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){2, ind_ang_1}(:,1)==0), :) = nan;
data_2 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){2, ind_ang_2}(:,range_);
data_2(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){2, ind_ang_2}(:,1)==0), :) = nan;
plot(ESN_smooth(nanmean(data_1-data_2)))
data_1 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){3, ind_ang_1}(:,range_);
data_1(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){3, ind_ang_1}(:,1)==0), :) = nan;
data_2 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){3, ind_ang_2}(:,range_);
data_2(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){3, ind_ang_2}(:,1)==0), :) = nan;
plot(ESN_smooth(nanmean(data_1-data_2)))
data_1 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){4, ind_ang_1}(:,range_);
data_1(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){4, ind_ang_1}(:,1)==0), :) = nan;
data_2 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){4, ind_ang_2}(:,range_);
data_2(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){4, ind_ang_2}(:,1)==0), :) = nan;
plot(ESN_smooth(nanmean(data_1-data_2)))
% data_1 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){5, ind_ang_1}(:,range_);
% data_1(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){5, ind_ang_1}(:,1)==0), :) = nan;
% data_2 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){5, ind_ang_2}(:,range_);
% data_2(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){5, ind_ang_2}(:,1)==0), :) = nan;
% plot(ESN_smooth(nanmean(data_1-data_2)))
data_1 = ALL_PCELL_DATA_2.(variable_name1).(variable_name2){1, ind_ang_1}(:,range_);
data_2 = ALL_PCELL_DATA_2.(variable_name1).(variable_name2){1, ind_ang_2}(:,range_);
plot(ESN_smooth(nanmean(data_1-data_2)), 'k')

subplot(1, 2, 2)
hold on
data_1 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){2, ind_ang_1}(:,range_);
data_1(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){2, ind_ang_1}(:,1)==0), :) = nan;
data_2 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){2, ind_ang_2}(:,range_);
data_2(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){2, ind_ang_2}(:,1)==0), :) = nan;
plot((nanmean((data_1+data_2)./2)))
data_1 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){3, ind_ang_1}(:,range_);
data_1(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){3, ind_ang_1}(:,1)==0), :) = nan;
data_2 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){3, ind_ang_2}(:,range_);
data_2(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){3, ind_ang_2}(:,1)==0), :) = nan;
plot((nanmean((data_1+data_2)./2)))
data_1 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){4, ind_ang_1}(:,range_);
data_1(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){4, ind_ang_1}(:,1)==0), :) = nan;
data_2 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){4, ind_ang_2}(:,range_);
data_2(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){4, ind_ang_2}(:,1)==0), :) = nan;
plot((nanmean((data_1+data_2)./2)))
% data_1 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){5, ind_ang_1}(:,range_);
% data_1(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){5, ind_ang_1}(:,1)==0), :) = nan;
% data_2 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){5, ind_ang_2}(:,range_);
% data_2(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){5, ind_ang_2}(:,1)==0), :) = nan;
% plot((nanmean((data_1+data_2)./2)))
data_1 = ALL_PCELL_DATA_2.(variable_name1).(variable_name2V){1, ind_ang_1}(:,range_);
data_2 = ALL_PCELL_DATA_2.(variable_name1).(variable_name2V){1, ind_ang_2}(:,range_);
plot((nanmean((data_1+data_2)./2)), 'k')
end


%% function scratch_plot4, USED for MLMC, Permutation
function scratch_plot4
%%
hFig = figure(4);
clf(hFig)

variable_name_category = 'all'; % 'corrSac'; % 'primSac'; % 'nonTask'; % 'toTgtStr'; % 'fromCenter'; % 
event_type = 'start'; % 'vmax'; % 'finish'; % 
spike_type = 'SS';
flag_cs_180_minus_cs_on = false;
variable_name_spike    = [spike_type '_train_' event_type];
variable_name_velocity = [spike_type '_velocity_' event_type];
variable_name_num_sac  = [spike_type '_num_sac_' event_type];
variable_name_firing   = [spike_type '_firing_rate'];
ind_amp = 6;
ind_ang = 5;
ind_ang_cs_180 = 5;
if contains(spike_type,'CS')
    range_ = 1:600;
else
    range_ = 201:400;
end

idx_pCells =  1:size(ALL_PCELL_all_saccades.all.SS_train_start{1,1},1); % idx_pauser; % idx_burster; % 
num_pCells = length(idx_pCells);

ALL_PCELL_DATA_1 = ALL_PCELL_all_saccades_tuned;
% ALL_PCELL_DATA_1 = ALL_PCELL_all_saccades_tuned_90_270_amp_ang;


num_sac = ALL_PCELL_DATA_1.(variable_name_category).(variable_name_num_sac){ind_amp, ind_ang}(idx_pCells, :);
firing_rate = ALL_PCELL_all_saccades.Neural_Properties.(variable_name_firing)(idx_pCells, :);
train_data = ALL_PCELL_DATA_1.(variable_name_category).(variable_name_spike){ind_amp, ind_ang}(idx_pCells, :)*1000;
train_data = train_data - repmat(firing_rate(:,1), 1, size(train_data, 2));
if flag_cs_180_minus_cs_on
num_sac_2 = ALL_PCELL_DATA_1.(variable_name_category).(variable_name_num_sac){ind_amp, ind_ang_cs_180}(idx_pCells, :);
train_data_2 = ALL_PCELL_DATA_1.(variable_name_category).(variable_name_spike){ind_amp, ind_ang_cs_180}(idx_pCells, :)*1000;
train_data_2 = train_data_2 - repmat(firing_rate(:,1), 1, size(train_data_2, 2));
end

num_perm = 5000;
num_samples = num_pCells ; % 50; % 
length_span = size(train_data, 2);
train_data_perm = nan(num_perm, length_span);
for counter_perm = 1 : num_perm
    inds_perm = randi(num_pCells,[num_samples 1]);
    train_data_ = train_data(inds_perm, :);
    train_data_(isnan(train_data_)) = 0;
    num_sac_ = num_sac(inds_perm, 1);
    firing_mean_ = num_sac_' * train_data_ ./ nansum(num_sac_);
    train_data_perm(counter_perm, :) = firing_mean_;
    
    if flag_cs_180_minus_cs_on
    inds_perm_2 = inds_perm; % randi(num_pCells,[num_samples 1]); % 
    train_data_2_ = train_data_2(inds_perm_2, :);
    train_data_2_(isnan(train_data_2_)) = 0;
    num_sac_2_ = num_sac_2(inds_perm_2, 1);
    firing_mean_2_ = num_sac_2_' * train_data_2_ ./ nansum(num_sac_2_);
    train_data_perm(counter_perm, :) = firing_mean_ - firing_mean_2_;
    end
    
end
data_mean = ESN_smooth( nanmean(train_data_perm, 1) );
data_stdv = ESN_smooth( nanstd( train_data_perm, 0, 1) ) ./ sqrt(num_samples) * 3; % 95% CI
data_mean = data_mean(:, range_);
data_stdv = data_stdv(:, range_);

subplot(1, 2, 1)
hold on
y_axis_data_stdv = [(data_mean + data_stdv) flip(data_mean - data_stdv)];
x_axis_data_stdv = [(1:1:length(data_mean)) (length(data_mean):-1:1)];
plot(x_axis_data_stdv, (y_axis_data_stdv), 'k', 'LineWidth', 0.25)
% plot((data_mean - data_stdv), 'k', 'LineWidth', 0.25)
plot((data_mean), 'k', 'LineWidth', 1)
if contains(spike_type,'CS')
    set(gca, 'XTick', 0:200:600)
    ylim([-1 2.5])
else
    ylim([-30 20])
end
xline(length(range_)/2)

subplot(1, 2, 2)
hold on
velocity_data = ALL_PCELL_DATA_1.(variable_name_category).(variable_name_velocity){ind_amp, ind_ang};% - ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){1, 5};
velocity_data((num_sac(:,1)==0), :) = nan;
velocity_data = velocity_data(:, range_);
% plot((nanmean(data_)+(nanstd(data_)./sqrt(num_pCells))), 'k', 'LineWidth', 1)
% plot((nanmean(data_)-(nanstd(data_)./sqrt(num_pCells))), 'k', 'LineWidth', 1)
plot((nanmean(velocity_data)), 'k', 'LineWidth', 1)
if contains(spike_type,'CS')
    set(gca, 'XTick', 0:200:600)
end
xline(length(range_)/2)
if contains(event_type,'vmax')
    ylim([0 750])
else
    ylim([0 650])
end

ESN_Beautify_Plot(hFig, [2.5 1.25], 8)

end

%% function classify_pCells
function classify_pCells
%%
clearvars -except ALL_PCELL_*
%% Build pca and umap data
data_ = ALL_PCELL_all_saccades_tuned_ang_amp.primSac.SS_train_start{1,1} * 1000;
SS_baseline = mean(data_(:,151:250), 2);
data_norm = data_ - repmat(SS_baseline, 1, size(data_, 2));
data_smooth = ESN_smooth(data_norm, 2);

[~, pca_mat, ~] = pca(data_smooth(:,250:400));
[reduction, umap, clusterIdentifiers, extras]=run_umap(data_smooth(:,250:400));

%% Cluster Using UMAP
data_gmm_ = [reduction(:, 1), reduction(:, 2)];
n_component_gmm_ = 2;
gmm_model_ = fitgmdist(data_gmm_,n_component_gmm_);
idx = cluster(gmm_model_,data_gmm_);
%% Plot UMAP
hFig = figure(2);
clf(hFig);
subplot(2, 2, 1)
hold on
plot((-49:100)',data_smooth(idx==1,251:400)', 'b')
plot((-49:100)',nanmean(data_smooth(idx==1,251:400))', 'k', 'LineWidth', 2)
plot((-49:100)',data_smooth(9,251:400)', 'm')
plot((-49:100)',data_smooth(27,251:400)', 'c')
ylim([-50 200])

subplot(2, 2, 3)
hold on
plot((-49:100)',data_smooth(idx==2,251:400)', 'r')
plot((-49:100)',nanmean(data_smooth(idx==2,251:400))', 'k', 'LineWidth', 2)
plot((-49:100)',data_smooth(110,251:400)', 'm')
plot((-49:100)',data_smooth(23,251:400)', 'c')
ylim([-125 125])
xlabel('Saccade onset (ms)')
ylabel('Normalized SS firing rate')
% ylim([0 inf])
MarkerSize_ = 3;
subplot(2, 2, [2 4])
hold on
plot(reduction(idx==1, 1), reduction(idx==1, 2), 'ob',...
    'MarkerFaceColor', 'b', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
plot(reduction(9, 1), reduction(9, 2), 'ob',...
    'MarkerFaceColor', 'm', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
plot(reduction(27, 1), reduction(27, 2), 'ob',...
    'MarkerFaceColor', 'c', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
plot(reduction(idx==2, 1), reduction(idx==2, 2), 'or', ...
    'MarkerFaceColor', 'r', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
plot(reduction(110, 1), reduction(110, 2), 'or',...
    'MarkerFaceColor', 'm', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
plot(reduction(23, 1), reduction(23, 2), 'or',...
    'MarkerFaceColor', 'c', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
xlabel('umap 1')
ylabel('umap 2')

ESN_Beautify_Plot(hFig, [4,2], 8)

%% Cluster Using PCA
data_gmm_ = [pca_mat(:, 1), pca_mat(:, 2)];
n_component_gmm_ = 2;
gmm_model_ = fitgmdist(data_gmm_,n_component_gmm_);
idx = cluster(gmm_model_,data_gmm_);
%% Plot PCA
hFig = figure(2);
clf(hFig);
subplot(2, 2, 1)
hold on
plot((-49:100)',data_smooth(idx==1,251:400)', 'b')
plot((-49:100)',nanmean(data_smooth(idx==1,251:400))', 'k', 'LineWidth', 2)
plot((-49:100)',data_smooth(9,251:400)', 'm')
plot((-49:100)',data_smooth(27,251:400)', 'c')
subplot(2, 2, 3)
hold on
plot((-49:100)',data_smooth(idx==2,251:400)', 'r')
plot((-49:100)',nanmean(data_smooth(idx==2,251:400))', 'k', 'LineWidth', 2)
plot((-49:100)',data_smooth(110,251:400)', 'm')
plot((-49:100)',data_smooth(23,251:400)', 'c')
xlabel('Saccade onset (ms)')
ylabel('Normalized SS firing rate')
% ylim([0 inf])
MarkerSize_ = 3;
subplot(2, 2, [2 4])
hold on
plot(pca_mat(idx==1, 1), pca_mat(idx==1, 2), 'ob',...
    'MarkerFaceColor', 'b', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
plot(pca_mat(9, 1), pca_mat(9, 2), 'ob',...
    'MarkerFaceColor', 'm', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
plot(pca_mat(27, 1), pca_mat(27, 2), 'ob',...
    'MarkerFaceColor', 'c', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
plot(pca_mat(idx==2, 1), pca_mat(idx==2, 2), 'or', ...
    'MarkerFaceColor', 'r', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
plot(pca_mat(110, 1), pca_mat(110, 2), 'or',...
    'MarkerFaceColor', 'm', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
plot(pca_mat(23, 1), pca_mat(23, 2), 'or',...
    'MarkerFaceColor', 'c', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
xlabel('pca 1')
ylabel('pca 2')

ESN_Beautify_Plot(hFig, [4,2], 8)
%%
idx_burster = find(idx==1);
idx_pauser = find(idx==2);
save('umap_results.mat', 'pca_mat', 'reduction', 'umap', 'clusterIdentifiers', 'extras', 'data_smooth', 'idx', 'idx_burster', 'idx_pauser', '-v7.3')

end

%% function sample_pCells
function sample_pCells
%% 
idx_cell_1 = 9;
idx_cell_2 = 106;
idx_cell_3 = 27;
idx_cell_4 = 23;

data_SS = ALL_PCELL_all_saccades_tuned_amp.primSac.SS_train_start{3,1} * 1000;
SS_baseline = mean(data_SS(:,151:250), 2);
data_SS = data_SS - repmat(SS_baseline, 1, size(data_SS, 2));
data_vel = ALL_PCELL_all_saccades_tuned_amp.primSac.SS_velocity_start{3,1};
data_SS_smooth = ESN_smooth(data_SS, 2);
hFig = figure(2);
clf(hFig);
subplot(1, 4, 1)
hold on
plot((-49:100)',data_SS_smooth(idx_cell_1,251:400)', 'b')
plot((-49:100)',data_SS_smooth(idx_cell_2,251:400)', 'r')
ylim([-100 200])
xlabel('Saccade onset (ms)')
ylabel('SS firing rate (change, Hz)')

subplot(1, 4, 3)
hold on
plot((-49:100)',data_SS_smooth(idx_cell_3,251:400)', 'b')
plot((-49:100)',data_SS_smooth(idx_cell_4,251:400)', 'r')
ylim([-100 200])
xlabel('Saccade onset (ms)')
ylabel('SS firing rate (change, Hz)')

subplot(1, 4, 2)
hold on
plot((-49:100)',data_vel(idx_cell_1,251:400)', 'b')
plot((-49:100)',data_vel(idx_cell_2,251:400)', 'r')
ylim([0 600])
xlabel('Saccade onset (ms)')
ylabel('Eye velocity (deg/s)')

subplot(1, 4, 4)
hold on
plot((-49:100)',data_vel(idx_cell_3,251:400)', 'b')
plot((-49:100)',data_vel(idx_cell_4,251:400)', 'r')
ylim([0 600])
xlabel('Saccade onset (ms)')
ylabel('Eye velocity (deg/s)')

ESN_Beautify_Plot(hFig, [8,1.5])
end

%% function combine_90_270
function ALL_PCELL_all_saccades_90_270 = combine_90_270(ALL_PCELL_all_saccades)
%% Loop over pCells
fprintf(['Building ALL_PCELL_all_saccades_90_270', ' ... ']);
clearvars ALL_PCELL_all_saccades_90_270
variable_list = {'primSac', 'corrSac', 'toTgt', 'fromTgt', 'lowAmp', 'toTgtStr', 'fromCenter', 'aroundCenter', 'nonTask', 'all', 'primCorrStr', 'tgtCuePres', 'tgtPrimSac', 'tgtCorrSac'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};
num_pCells = size(ALL_PCELL_all_saccades.primSac.SS_train_start{1,1}, 1);
span_width = size(ALL_PCELL_all_saccades.primSac.SS_train_start{1,1}, 2);
for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    num_amp_bin = size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]),1);
    num_ang_bin = size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]),2);
    
    for counter_amp = 1 : num_amp_bin
        train_data_all_ang = zeros(size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, 1}));
        velocity_data_all_ang = zeros(size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, 1}));
        num_saccades_all_ang = zeros(size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, 1}));
        if size(num_saccades_all_ang, 2) == 1
            num_saccades_all_ang = repmat(num_saccades_all_ang, 1, span_width);
        end
        for counter_ang = [3, 7]
            train_data = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang};
            velocity_data = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang};
            num_saccades = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, counter_ang};
            if size(num_saccades, 2) == 1
                num_saccades = repmat(num_saccades, 1, span_width);
            end
            train_data_all_ang = train_data_all_ang + (train_data .* num_saccades);
            velocity_data_all_ang = velocity_data_all_ang + (velocity_data .* num_saccades);
            num_saccades_all_ang = num_saccades_all_ang + num_saccades;
        end
        train_data_all_ang = train_data_all_ang ./ num_saccades_all_ang;
        train_data_all_ang(isnan(train_data_all_ang)) = 0;
        velocity_data_all_ang = velocity_data_all_ang ./ num_saccades_all_ang;
        velocity_data_all_ang(isnan(velocity_data_all_ang)) = 0;
        ALL_PCELL_all_saccades_90_270.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, 1} = train_data_all_ang;
        ALL_PCELL_all_saccades_90_270.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, 1} = velocity_data_all_ang;
        ALL_PCELL_all_saccades_90_270.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, 1} = num_saccades_all_ang;
    end
end
end
end
fprintf(' --> Completed. \n')
end

%%
function cs_on_circular_gaussian
%% METHOD 3: mean direction
ALL_PCELL_all_saccades_ang = avg_over_amplitude(ALL_PCELL_all_saccades);
num_ang_bin = size(ALL_PCELL_all_saccades_ang.tgtCuePres.CS_train_start,2);
num_pCells = size(ALL_PCELL_all_saccades_ang.tgtCuePres.CS_train_start{1,1},1);
CS_prob_primSac = zeros(num_pCells, num_ang_bin);
CS_prob_corrSac = zeros(num_pCells, num_ang_bin);
for counter_ang = 1 : num_ang_bin
    CS_prob_primSac(:,counter_ang) = nansum(ALL_PCELL_all_saccades_ang.primSac.CS_train_start{1,counter_ang}(:,100:300), 2);
    CS_prob_corrSac(:,counter_ang) = nansum(ALL_PCELL_all_saccades_ang.corrSac.CS_train_start{1,counter_ang}(:,100:300), 2);
end

CS_prob_primSac(1:51, [2 4 6 8]) = nan;
CS_prob_corrSac(1:51, [2 4 6 8]) = nan;
CS_prob_avg = (CS_prob_primSac + CS_prob_corrSac) ./ 2;
% [~, idx_CS_max] = nanmax(CS_prob_avg, [], 2);
dir_angles = deg2rad([180, 225, 270, 315, 0, 45, 90, 135]);
% compute weighted sum of cos and sin of angles
r = nansum(CS_prob_avg.* repmat(exp(1i*dir_angles), num_pCells, 1) , 2);
% Computes the mean direction for circular data.
CS_ang_avg = rad2deg(angle(r));
% sum of weights
CS_prob_sum = nansum(CS_prob_avg,2);
% Computes mean resultant vector length for circular data.
CS_rho_avg = abs(r) ./ CS_prob_sum;

%% rotate coordinate to get CS_prob_avg_tuned
ang_edges = -202.5 : 45 : +202.5;
ang_edges_last_bin_id = length(ang_edges) - 1;
[~, ~, idx_CS_max] = histcounts(CS_ang_avg, ang_edges);
idx_CS_max(idx_CS_max == ang_edges_last_bin_id) = 1;

idx_CS_tuning = zeros(num_pCells, num_ang_bin);
CS_prob_avg_tuned = nan(num_pCells, num_ang_bin);
for counter_pCell = 1 : 1 : num_pCells
    CS_on_index = idx_CS_max(counter_pCell);
    CS_on_index = CS_on_index - 1; % CS_on_index should be in 0-index format
    if (CS_on_index == 8); CS_on_index = 0; end
    idx_CS_tuning_ = mod((CS_on_index : 1 : CS_on_index+7), 8) + 1;
    idx_CS_tuning(counter_pCell, :) = idx_CS_tuning_;
    CS_prob_avg_tuned(counter_pCell, :) = CS_prob_avg(counter_pCell, idx_CS_tuning_);
end

%% Circular stats 
% Rayleigh's test for nonuniformity
R = CS_prob_sum.*CS_rho_avg;
% compute Rayleigh's z (equ. 27.2)
z = R.^2 ./ CS_prob_sum;
% compute p value using approxation in Zar, p. 617
pval = exp(sqrt(1+4.*CS_prob_sum+4*(CS_prob_sum.^2-R.^2))-(1+2.*CS_prob_sum));


% CS_rho_avg = sqrt(nansum(CS_vec_avg.^2,2)) ./ nansum(CS_prob_avg,2);
% ang_edges = -202.5 : 45 : +202.5;
% ang_edges_last_bin_id = length(ang_edges) - 1;
% [~, ~, idx_CS_max] = histcounts(CS_ang_avg, ang_edges);
% idx_CS_max(idx_CS_max == ang_edges_last_bin_id) = 1;

%% von Mises
ang_ = deg2rad([180 225 270 315 0 45 90 135]);
ang_exp_ = exp(1i*ang_);
theta_all = nan(num_pCells, 1);
r_all = nan(num_pCells, 1);
kappa_all = nan(num_pCells, 1);
circ_sig_all = nan(num_pCells, 1);
vonMises_pdf_all = nan(size(CS_prob_avg));
for counter_pcell = 1 : num_pCells
    CS_prob_ = CS_prob_avg(counter_pcell,:);
    
    % compute weighted sum of cos and sin of angles
    r = nansum(CS_prob_.* ang_exp_);
    % Computes the mean direction for circular data.
    thetahat = angle(r);
    % Computes mean resultant vector length for circular data.
    r = abs(r)./nansum(CS_prob_);
    
    R = r;
    if R < 0.53
        kappa = 2*R + R^3 + 5*R^5/6;
    elseif R>=0.53 && R<0.85
        kappa = -.4 + 1.39*R + 0.43/(1-R);
    else
        kappa = 1/(R^3 - 4*R^2 + 3*R);
    end
    
    % evaluate pdf
    C = 1/(2*pi*besseli(0,kappa));
    vonMises_pdf = C * exp(kappa*cos(ang_-thetahat));
    circ_var = 1 - (besseli(1,kappa) / besseli(0,kappa));
    
    theta_all(counter_pcell, :) = rad2deg(thetahat);
    r_all(counter_pcell, :) = r;
    kappa_all(counter_pcell, :) = kappa;
    circ_sig_all(counter_pcell, :) = rad2deg(sqrt(circ_var));
    vonMises_pdf_all(counter_pcell, :) = vonMises_pdf;
end

%% Plot CS-on distribution
hFig = figure(1); clf(hFig);
len_1 = length(pCell_list_1);
len_2 = length(pCell_list_2);
len_3 = length(pCell_list_3);
% len_4 = length(pCell_list_4);
len_tot = length(pCell_list);
% idx_pCells = 1:len_tot; % All
% idx_pCells = (len_1+len_2+len_3+1):len_tot; % Mirza L
% idx_pCells = (len_1+len_2+1):(len_1+len_2+len_3); % Ramon R
% idx_pCells = 1:22; % Mirza C
idx_pCells = 23:(len_1+len_2); % Mirza R-2

step_size_ = pi/8;
ang_edges = -pi-(step_size_/2):step_size_:pi-(step_size_/2);
polarhistogram(deg2rad(CS_ang_avg(idx_pCells)), ang_edges, 'DisplayStyle', 'bar','FaceColor','red', 'EdgeColor', 'none')
hold on
polarhistogram(deg2rad(CS_ang_avg(idx_pCells)), ang_edges, 'DisplayStyle', 'stairs','FaceColor','none', 'EdgeColor', 'r', 'linewidth', 1)
% polarhistogram(deg2rad(CS_ang_avg), ang_edges, 'DisplayStyle', 'stairs','FaceColor','none', 'EdgeColor', 'b', 'linewidth', 1)
set(gca, 'ThetaTick', 0:45:315, 'RTick', 0:5:25)
rlim([0 20])
ESN_Beautify_Plot(hFig, [2 2], 8)

%% Plot std distribution
hFig = figure(2); clf(hFig);
hold on
circ_sig_edges = 39:2:61;
histogram(circ_sig_all, circ_sig_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'r')
histogram(circ_sig_all, circ_sig_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'r', 'FaceColor', 'none', 'linewidth', 1)
xline(mean(circ_sig_all),'Color', 'r', 'linewidth', 1)
ESN_Beautify_Plot(hFig, [2, 2], 8)

%% Plot CS tuning
overall_prob_TUNED_mean = nanmean(CS_prob_avg_tuned, 1);
overall_prob_TUNED_stdv = nanstd(CS_prob_avg_tuned, 0, 1) ./ sqrt(num_pCells);
overall_prob_TUNED_stdv_plus = overall_prob_TUNED_mean + overall_prob_TUNED_stdv;
overall_prob_TUNED_stdv_minus = overall_prob_TUNED_mean - overall_prob_TUNED_stdv;

plot_data_amp_mean = [overall_prob_TUNED_mean, overall_prob_TUNED_mean(1), nan]';
plot_data_deg_mean = [0, 45, 90, 135, 180, 225, 270, 315, 0, nan]';
plot_data_x_axis_mean = plot_data_amp_mean .* cosd(plot_data_deg_mean);
plot_data_y_axis_mean = plot_data_amp_mean .* sind(plot_data_deg_mean);

plot_data_amp_stdv_p = [overall_prob_TUNED_stdv_plus, overall_prob_TUNED_stdv_plus(1), nan]';
plot_data_deg_stdv_p = [0, 45, 90, 135, 180, 225, 270, 315, 0, nan]';
plot_data_x_axis_stdv_p = plot_data_amp_stdv_p .* cosd(plot_data_deg_stdv_p);
plot_data_y_axis_stdv_p = plot_data_amp_stdv_p .* sind(plot_data_deg_stdv_p);

plot_data_amp_stdv_m = [overall_prob_TUNED_stdv_minus, overall_prob_TUNED_stdv_minus(1), nan]';
plot_data_deg_stdv_m = [0, 45, 90, 135, 180, 225, 270, 315, 0, nan]';
plot_data_x_axis_stdv_m = plot_data_amp_stdv_m .* cosd(plot_data_deg_stdv_m);
plot_data_y_axis_stdv_m = plot_data_amp_stdv_m .* sind(plot_data_deg_stdv_m);

hFig = figure(3);
clf(hFig);

hold on
plot([0.4, -0.4, nan, 0, -0, nan,]', [0, 0, nan, 0.4, -0.4, nan], 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.05, sind(0:5:360)*0.05, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
plot(cosd(0:5:360)*0.10, sind(0:5:360)*0.10, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
plot(cosd(0:5:360)*0.15, sind(0:5:360)*0.15, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
plot(cosd(0:5:360)*0.20, sind(0:5:360)*0.20, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
plot(cosd(0:5:360)*0.25, sind(0:5:360)*0.25, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
plot(cosd(0:5:360)*0.30, sind(0:5:360)*0.30, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])

plot_data_x_axis_mean   = plot_data_x_axis_mean(  ~isnan(plot_data_x_axis_mean)); 
plot_data_y_axis_mean   = plot_data_y_axis_mean(  ~isnan(plot_data_y_axis_mean));
plot_data_x_axis_stdv_p = plot_data_x_axis_stdv_p(~isnan(plot_data_x_axis_stdv_p)); 
plot_data_y_axis_stdv_p = plot_data_y_axis_stdv_p(~isnan(plot_data_y_axis_stdv_p));
plot_data_x_axis_stdv_m = plot_data_x_axis_stdv_m(~isnan(plot_data_x_axis_stdv_m)); 
plot_data_y_axis_stdv_m = plot_data_y_axis_stdv_m(~isnan(plot_data_y_axis_stdv_m));

plot(plot_data_x_axis_mean,   plot_data_y_axis_mean, '-k', 'LineWidth', 1)
plot(plot_data_x_axis_stdv_p, plot_data_y_axis_stdv_p, '-k', 'LineWidth', 0.5)
plot(plot_data_x_axis_stdv_m, plot_data_y_axis_stdv_m, '-k', 'LineWidth', 0.5)

axis equal
xlim([-0.3 0.3])
ylim([-0.3 0.3])
% xlabel('CS probability');
% ylabel('CS probability');
set(gca, 'XTick', -0.3:0.1:0.3, 'YTick', -0.3:0.1:0.3)

ESN_Beautify_Plot(hFig, [1.5, 1.5], 8)
title('CS Tuning', 'Interpreter', 'none');

hFig = figure(4);
clf(hFig);
hold on
% plot_order_ = [7 8 1 2 3 4 5 6];
plot_order_ = [5 6 7 8 1 2 3 4 5];
plot(overall_prob_TUNED_stdv_plus(plot_order_), '-k', 'LineWidth', 0.5)
plot(overall_prob_TUNED_stdv_minus(plot_order_), '-k', 'LineWidth', 0.5)
plot(overall_prob_TUNED_mean(plot_order_), '-k', 'LineWidth', 1)
ylabel('CS probability');
xlabel('Direction')
% set(gca, 'XTick', 1:1:8, 'XTickLabel', {'-90','-45','ON','45','90','135','180','225'})
set(gca, 'XTick', 1:1:9, 'XTickLabel', {'-180', '-135', '-90','-45','ON','45','90','135','180'})

avg_prob_TUNED_mean = nanmean(nanmean(CS_prob_avg_tuned, 2));
avg_prob_TUNED_stdv = nanstd(nanmean(CS_prob_avg_tuned, 2), 0, 1) ./ sqrt(num_pCells);
avg_prob_TUNED_mean = repmat(avg_prob_TUNED_mean, 1, size(CS_prob_avg_tuned,2));
avg_prob_TUNED_stdv = repmat(avg_prob_TUNED_stdv, 1, size(CS_prob_avg_tuned,2));
avg_prob_TUNED_stdv_plus = avg_prob_TUNED_mean + avg_prob_TUNED_stdv;
avg_prob_TUNED_stdv_minus = avg_prob_TUNED_mean - avg_prob_TUNED_stdv;


plot(avg_prob_TUNED_stdv_plus(plot_order_), '-k', 'LineWidth', 0.5)
plot(avg_prob_TUNED_stdv_minus(plot_order_), '-k', 'LineWidth', 0.5)
plot(avg_prob_TUNED_mean(plot_order_), '-k', 'LineWidth', 1)
ylim([0.10 0.30])

ESN_Beautify_Plot(hFig, [2, 1.5], 8)
title('CS Tuning', 'Interpreter', 'none');

%% PAIR_DATA find CS-on for each pair
PAIR_DATA_all_saccades_ang = pair_data_avg_over_amplitude(PAIR_DATA_all_saccades);
num_ang_bin = size(PAIR_DATA_all_saccades_ang.primSac.CS_Pr_CH1_start,2);
num_pairs   = size(PAIR_DATA_all_saccades_ang.primSac.CS_Pr_CH1_start{1,1},1);
%
span_length_ = size(PAIR_DATA_all_saccades_ang.primSac.CS_Pr_CH1_start{1,1},2);
range_before_sac = round(span_length_/6) : 3*round(span_length_/6);
CS_prob_primSac_1 = zeros(num_pairs, num_ang_bin);
CS_prob_corrSac_1 = zeros(num_pairs, num_ang_bin);
CS_prob_primSac_2 = zeros(num_pairs, num_ang_bin);
CS_prob_corrSac_2 = zeros(num_pairs, num_ang_bin);
for counter_ang = 1 : num_ang_bin
    CS_prob_primSac_1(:,counter_ang) = nansum(PAIR_DATA_all_saccades_ang.primSac.CS_Pr_CH1_start{1,counter_ang}(:,range_before_sac), 2);
    CS_prob_corrSac_1(:,counter_ang) = nansum(PAIR_DATA_all_saccades_ang.corrSac.CS_Pr_CH1_start{1,counter_ang}(:,range_before_sac), 2);
    CS_prob_primSac_2(:,counter_ang) = nansum(PAIR_DATA_all_saccades_ang.primSac.CS_Pr_CH2_start{1,counter_ang}(:,range_before_sac), 2);
    CS_prob_corrSac_2(:,counter_ang) = nansum(PAIR_DATA_all_saccades_ang.corrSac.CS_Pr_CH2_start{1,counter_ang}(:,range_before_sac), 2);
end

CS_prob_primSac_1(1:16, [2 4 6 8]) = nan;
CS_prob_corrSac_1(1:16, [2 4 6 8]) = nan;
CS_prob_primSac_2(1:16, [2 4 6 8]) = nan;
CS_prob_corrSac_2(1:16, [2 4 6 8]) = nan;
dir_angles = deg2rad([180, 225, 270, 315, 0, 45, 90, 135]);

CS_prob_avg = (CS_prob_primSac_1 + CS_prob_corrSac_1 + CS_prob_primSac_2 + CS_prob_corrSac_2) ./ 4;
r = nansum(CS_prob_avg.* repmat(exp(1i*dir_angles), num_pairs, 1) , 2);
CS_ang_avg = rad2deg(angle(r));
CS_prob_sum = nansum(CS_prob_avg,2);
CS_rho_avg = abs(r) ./ CS_prob_sum;

CS_prob_avg_1 = (CS_prob_primSac_1 + CS_prob_corrSac_1) ./ 2;
r_1 = nansum(CS_prob_avg_1.* repmat(exp(1i*dir_angles), num_pairs, 1) , 2);
CS_ang_avg_1 = rad2deg(angle(r_1));
CS_prob_sum_1 = nansum(CS_prob_avg_1,2);
CS_rho_avg_1 = abs(r_1) ./ CS_prob_sum_1;

CS_prob_avg_2 = (CS_prob_primSac_2 + CS_prob_corrSac_2) ./ 2;
r_2 = nansum(CS_prob_avg_2.* repmat(exp(1i*dir_angles), num_pairs, 1) , 2);
CS_ang_avg_2 = rad2deg(angle(r_2));
CS_prob_sum_2 = nansum(CS_prob_avg_2,2);
CS_rho_avg_2 = abs(r_2) ./ CS_prob_sum_2;

%% Plot diff between pair cs-on
hFig = figure(3); clf(hFig);
CS_ang_pairs = [CS_ang_avg_1 CS_ang_avg_2];
x_values_ = cosd(CS_ang_pairs); 
y_values_ = sind(CS_ang_pairs);
diff_theta_pairs = acosd( (x_values_(:,1) .* x_values_(:,2)) + (y_values_(:,1) .* y_values_(:,2)) );
% diff_theta_pairs = abs(theta_pairs(:,1) - theta_pairs(:,2));
% histogram(diff_theta_pairs)
step_size_ = pi/8;
ang_edges = -pi-(step_size_/2):step_size_:pi-(step_size_/2);
polarhistogram(deg2rad(diff_theta_pairs), ang_edges, 'DisplayStyle', 'bar','FaceColor','red', 'EdgeColor', 'none')
hold on
polarhistogram(deg2rad(diff_theta_pairs), ang_edges, 'DisplayStyle', 'stairs','FaceColor','none', 'EdgeColor', 'r', 'linewidth', 1)
% polarhistogram(deg2rad(CS_ang_avg), ang_edges, 'DisplayStyle', 'stairs','FaceColor','none', 'EdgeColor', 'b', 'linewidth', 1)
set(gca, 'ThetaTick', 0:45:315, 'RTick', 0:5:25)
rlim([0 15])
ESN_Beautify_Plot(hFig, [2 2], 8)

%% Plot CS_rho_avg
hFig = figure(2); clf(hFig); 
hold('on');
plot([0.5, -0.5, nan, 0, -0, nan,]', [0, 0, nan, 0.5, -0.5, nan], 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.10, sind(0:5:360)*0.10, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
plot(cosd(0:5:360)*0.20, sind(0:5:360)*0.20, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
plot(cosd(0:5:360)*0.30, sind(0:5:360)*0.30, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
plot(cosd(0:5:360)*0.40, sind(0:5:360)*0.40, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
plot(cosd(0:5:360)*0.50, sind(0:5:360)*0.50, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
% plot([zeros(1,num_pCells); CS_vec_avg(:,1)'], [zeros(1,num_pCells); CS_vec_avg(:,2)'], 'LineWidth', 2)
x_axis_ = CS_rho_avg .* cosd( CS_ang_avg );
y_axis_ = CS_rho_avg .* sind( CS_ang_avg );
plot([zeros(1,num_pCells); x_axis_'], [zeros(1,num_pCells); y_axis_'], 'LineWidth', 2)
axis equal

%% Plot sample pCell CS-on tuning
idx = 72;
CS_prob_ = CS_prob_avg(idx,:);
x_ = deg2rad([180 225 270 315 0 45 90 135]);
y_ = CS_prob_;
x_pref = deg2rad(CS_ang_avg(idx));
A_ = 0.0554;
B_ = 3.6416;
C_ = 0.1022;
params_init = [x_pref, A_, B_, C_];
params = fit_circular_gaussian(CS_prob_, params_init);
disp(['prefered dir: ' num2str(rad2deg(params(1)))]);
disp(['params: ' mat2str(params)]);

x_pref = params(1);
A_ = params(2);
B_ = params(3);
C_ = params(4);
y_hat_ = A_ * exp(B_ * cos( x_ - x_pref ) ) + C_;

[x_sort,idx_sort] = sort(x_);
y_sort = y_(idx_sort);
y_hat_sort = y_hat_(idx_sort);

clf(figure(1))

subplot(1,2,1)
hold on
plot(rad2deg(x_sort), y_sort, 'o-')
hold on
plot(rad2deg(x_sort), y_hat_sort, '.-')
xline(rad2deg(x_pref))

subplot(1,2,2)
polarplot([x_sort x_sort(1)],[y_sort y_sort(1)], 'o-')
hold on
polarplot([x_sort x_sort(1)],[y_hat_sort y_hat_sort(1)], '.-')
polarplot([x_pref x_pref],rlim, 'k')
end

%% von mises
function params = fit_circular_gaussian(CS_prob_, params_init)
x_ = deg2rad([180 225 270 315 0 45 90 135]);
y_ = CS_prob_;

function [rms] = f_search(params)
    x_pref = params(1);
    A_ = params(2);
    B_ = params(3);
    C_ = params(4);
    y_hat_ = A_ * exp(B_ * cos( x_ - x_pref ) ) + C_;
    
    rms = sqrt( nansum((y_hat_ - y_).^2) );
end

% Optimization
% makes sure that fmincon does not print results to screen
options = optimset('Display','off');

% performs minimization
% x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
lb_ = [0 0 0 0];
params = fmincon(@f_search,params_init,[],[],[],[],lb_,[],[],options);
% params = fminsearch(@f_search,params_init,options);

end

%% function build_pair_id_num
function pair_id_num = build_pair_id_num()
%%
id_Mirza_pre201906 = ...
[1, 2; 3, 4; 6, 7; 8, 9; 17, 18; 19, 20; 23, 24; 25, 26; 31, 32; 35, 36; 37, 38; 39, 40; 41, 42; 41, 43; 42, 43; 48, 49; ];
id_Mirza_post201906 = ...
[13, 14];
id_Ramon = ...
[2, 3; 4, 5; 15, 16; 15, 17; 16, 17; 21, 22; 24, 25; 27, 28; 32, 33; 36, 37; ];
id_Mirza_post202011 = ...
[1, 2; 9, 10; 11, 12; 14, 15; 14, 16; 15, 16; 17, 18; 20, 21; ];
%%
pCell_list_1 = build_pCell_list_Mirza_pre201906();
pCell_list_2 = build_pCell_list_Mirza_post201906();
pCell_list_3 = build_pCell_list_Ramon();
% pCell_list_4 = build_pCell_list_Mirza_post202011();
id_Mirza_post201906 = id_Mirza_post201906 + size(pCell_list_1, 1);
id_Ramon = id_Ramon + size(pCell_list_1, 1) + size(pCell_list_2, 1);
id_Mirza_post202011 = id_Mirza_post202011 + size(pCell_list_1, 1) + size(pCell_list_2, 1) + size(pCell_list_3, 1);
pair_id_num = [id_Mirza_pre201906; id_Mirza_post201906; id_Ramon; id_Mirza_post202011];
%%

end


%% function synchrony
function synchrony()
% The data refers to pairs of simultaneously recorded neurons.
%
% Our null hypothesis (H0):
% Pr(CS1=1,CS2=1 | err=1) = Pr(CS1=1 | err=1) Pr(CS2=1 | err=1)
%
% Main hypothesis (H1):
% Pr(CS1=1,CS2=1 | err=1) = Pr(CS1=1 | CS2=1, err=1) Pr(CS2=1 | err=1)
%
% Thus we compute the ratio of H1 to H0:
% H1/H0 = Pr(CS1=1 | CS2=1, err=1) / Pr(CS1=1 | err=1) 
%
% To do this, you compute the probabilities for say 2 ms bins from -100 to +200 ms 
% after the error event and then simply plot H0, H1, and the ratio. 
% Given that these are binomial random variables, if Pr(A=1) is q, then its variance 
% is q-q^2. Thus, one can compute z-scores from the ratio of H1 to H0.

%% build pair list
path_data_monkey_sorted = uigetdir;
cd(path_data_monkey_sorted);
% build pCell_list, this is a hard coded cell with the id of all of the pCells and the bundles
pair_list_1 = build_pair_list_Mirza_pre201906();
pair_list_2 = build_pair_list_Mirza_post201906();
pair_list_3 = build_pair_list_Ramon();
pair_list_4 = build_pair_list_Mirza_post202011();
pair_list = vertcat(pair_list_1, pair_list_2, pair_list_3, pair_list_4);
%% load pair data
num_pair = size(pair_list,1);
pair_data = cell(num_pair,2);
for counter_pair = 1 : num_pair
    fprintf(['Loading pair ' num2str(counter_pair) '/' num2str(num_pair) ' ...'])
    pair_data{counter_pair, 1} = load(pair_list{counter_pair, 1});
    pair_data{counter_pair, 2} = load(pair_list{counter_pair, 2});
    fprintf(' --> Completed. \n')
end
%% remove fields
sac_type_list = {'cue_present', 'primSac_onset', 'primSac_offset', 'corrSac_onset'};
sac_type_list2 = {'cue_present', 'primSac_onset', 'primSac_vmax', 'primSac_offset', 'corrSac_onset', 'corrSac_vmax', 'corrSac_offset'};
dir_list = {'000', '045', '090', '135', '180', '225', '270', '315'};
for counter_pair = 1 : num_pair
    for counter_ch = 1 : 2
        data_ch_ = pair_data{counter_pair, counter_ch};
        for counter_sac_type = 1 : length(sac_type_list)
            data_ch_ = rmfield(data_ch_,['plot_data' '_' sac_type_list{counter_sac_type}]);
        end
        for counter_sac_type = 1 : length(sac_type_list)
            data_raster_ = data_ch_.(['raster_data' '_' sac_type_list{counter_sac_type}]);
            for counter_dir = 1 : length(dir_list)
                for counter_sac_type2 = 1 : length(sac_type_list2)    
                    data_raster_ = rmfield(data_raster_,['train_data_logic' '_' sac_type_list2{counter_sac_type2} '_' dir_list{counter_dir}]);
                end
            data_raster_ = rmfield(data_raster_,['velocity_data' '_' dir_list{counter_dir}]);
            end
            data_ch_.(['raster_data' '_' sac_type_list{counter_sac_type}]) = data_raster_;
        end
        pair_data{counter_pair, counter_ch} = data_ch_;
    end
end

%% save pair_data
fprintf(['Saving pair data' ' ...'])
save([path_data_monkey_sorted filesep 'pair_data.mat'],...
    'pair_data', 'pair_list',...
    '-v7.3')
fprintf(' --> Completed. \n')

%% load pair_data
fprintf(['Loading pair data' ' ...'])
load('pair_data.mat','pair_data', 'pair_list')
fprintf(' --> Completed. \n')

%% find CS-on for each pair
dir_list = {'000', '045', '090', '135', '180', '225', '270', '315'};
sac_type_list = {'cue_present', 'primSac_onset', 'primSac_offset', 'corrSac_onset'};
sac_type_range = {301:500, 101:300, 301:500, 101:300};

num_dir = length(dir_list);
num_pairs = size(pair_data,1);
CS_prob = cell(0,2);
for counter_pair = 1 : num_pairs
    for counter_pCell = 1 : 2
        for counter_sac_type = 1 : length(sac_type_list)
            for counter_dir = 1 : length(dir_list)
                data_raster_ = pair_data{counter_pair, counter_pCell}.(['raster_data' '_' sac_type_list{counter_sac_type}]).(['train_data_logic_CS_' dir_list{counter_dir}]);
                if isempty(data_raster_)
                    prob_ = nan;
                else
                    prob_ = nansum( nansum(data_raster_(:,sac_type_range{counter_sac_type}), 2) ) ./ size(data_raster_, 1);
                end
                CS_prob{counter_sac_type, counter_pCell}(counter_pair, counter_dir) = prob_;
            end
        end
    end
end

CS_prob_avg = zeros(num_pairs, num_dir);
for counter_sac_type = 1 : length(sac_type_list)
    for counter_pCell = 1 : 2
        CS_prob_avg = CS_prob_avg + CS_prob{counter_sac_type, counter_pCell};
    end
end
CS_prob_avg = CS_prob_avg ./ 8;

dir_angles = deg2rad(str2double(dir_list));
% compute weighted sum of cos and sin of angles
r = nansum(CS_prob_avg.* repmat(exp(1i*dir_angles), num_pairs, 1) , 2);
% Computes the mean direction for circular data.
CS_ang_avg = rad2deg(angle(r));
% sum of weights
CS_prob_sum = nansum(CS_prob_avg,2);
% Computes mean resultant vector length for circular data.
CS_rho_avg = abs(r) ./ CS_prob_sum;

ang_edges = -202.5 : 45 : +202.5;
ang_edges_last_bin_id = length(ang_edges) - 1;
[~, ~, idx_CS_max] = histcounts(CS_ang_avg, ang_edges);
idx_CS_max(idx_CS_max == ang_edges_last_bin_id) = 1;

CS_dir.on = cell(num_pairs, 1);
CS_dir.off = cell(num_pairs, 1);
dir_list = {'180', '225', '270', '315', '000', '045', '090', '135'};
for counter_pair = 1 : num_pairs
    idx_on_ = idx_CS_max(counter_pair);
    idx_off_ = mod(idx_on_ + 3, 8) + 1;
    CS_dir.on{counter_pair, 1} = dir_list{idx_on_};
    CS_dir.off{counter_pair, 1} = dir_list{idx_off_};
end
dir_list = {'000', '045', '090', '135', '180', '225', '270', '315'};
CS_dir.all = repmat(dir_list, num_pairs, 1);

%% Building raster data
clearvars -except pair_data pair_list raster_data prob_data CS_dir
fprintf(['Building raster data' ' ...'])
num_pair = size(pair_list,1);
raster_data = cell(num_pair,2);
spike_type = 'SS'; % 'CS'; % 
event_type = 'primSac_onset'; % 'cue_present'; % 'primSac_offset'; % 
dir_type   = 'all'; % 'off'; % 'on'; % 
if contains(dir_type, 'on')
    dir_list = CS_dir.on;
elseif contains(dir_type, 'off')
    dir_list = CS_dir.off;
else
    dir_list = CS_dir.all;
end
for counter_pair = 1 : num_pair
    raster_CH1 = [];
    raster_CH2 = [];
%     idx_1 = double(rand>0.5);
%     idx_2 = mod(idx_1+1,2);
    for counter_dir = 1 : size(dir_list, 2)
        dir_name = dir_list{counter_pair, counter_dir};
        raster_CH1 = vertcat(raster_CH1, pair_data{counter_pair,1}.(['raster_data_' event_type]).(['train_data_logic_' spike_type '_' dir_name]) );
        raster_CH2 = vertcat(raster_CH2, pair_data{counter_pair,2}.(['raster_data_' event_type]).(['train_data_logic_' spike_type '_' dir_name]) );
    end
    raster_data{counter_pair, 1} = raster_CH1;
    raster_data{counter_pair, 2} = raster_CH2;
end
fprintf(' --> Completed. \n')

%% Building Prob data
clearvars -except pair_data pair_list raster_data prob_data CS_dir spike_type event_type dir_type
fprintf(['Building Prob data' ' ...'])
num_pair = size(pair_list,1);
bin_size = 5;
idx_bin = 1 : bin_size : 601;
prob_data = struct;
for counter_pair = 1 : num_pair
    raster_CH1 = raster_data{counter_pair, 1};
    raster_CH2 = raster_data{counter_pair, 2};
    num_trials = min([size(raster_CH1, 1) size(raster_CH2, 1)]);
    raster_CH1 = logical(raster_CH1(1:num_trials, :));
    raster_CH2 = logical(raster_CH2(1:num_trials, :));
    raster_CH1_bin = false(num_trials, length(idx_bin)-1);
    raster_CH2_bin = false(num_trials, length(idx_bin)-1);
    for counter_idx = 1 : 1 : length(idx_bin)-1
        idx_i = idx_bin(counter_idx);
        idx_f = idx_bin(counter_idx+1)-1;
        raster_CH1_bin(:, counter_idx) = logical(sum(raster_CH1(:,idx_i:idx_f),2));
        raster_CH2_bin(:, counter_idx) = logical(sum(raster_CH2(:,idx_i:idx_f),2));
    end
    Pr_CH1 = sum(raster_CH1_bin) ./ num_trials;
    Pr_CH2 = sum(raster_CH2_bin) ./ num_trials;
    Pr_CH1_and_CH2 = sum(raster_CH1_bin & raster_CH2_bin) ./ num_trials;
    Pr_CH1_or_CH2  = sum(raster_CH1_bin | raster_CH2_bin) ./ num_trials;
    Pr_CH1_given_CH2 = zeros(size(Pr_CH1));
    Pr_CH2_given_CH1 = zeros(size(Pr_CH2));
    for counter_col = 1 : 1 : length(Pr_CH1)
        CS_CH1_bin_col = raster_CH1_bin(:,counter_col);
        CS_CH2_bin_col = raster_CH2_bin(:,counter_col);
        Pr_CH1_given_CH2(counter_col) = sum(CS_CH1_bin_col(CS_CH2_bin_col)) ./ sum(CS_CH2_bin_col);
        Pr_CH2_given_CH1(counter_col) = sum(CS_CH2_bin_col(CS_CH1_bin_col)) ./ sum(CS_CH1_bin_col);
    end
    ratio = Pr_CH1_given_CH2./Pr_CH1;
    ratio(isnan(ratio)) = 0;
    prob_data.Pr_CH1(counter_pair, :) = Pr_CH1;
    prob_data.Pr_CH2(counter_pair, :) = Pr_CH2;
    prob_data.Pr_CH1_and_CH2(counter_pair, :) = Pr_CH1_and_CH2;
    prob_data.Pr_CH1_or_CH2(counter_pair, :) = Pr_CH1_or_CH2;
    prob_data.Pr_CH1_given_CH2(counter_pair, :) = Pr_CH1_given_CH2;
    prob_data.Pr_CH2_given_CH1(counter_pair, :) = Pr_CH2_given_CH1;
    prob_data.H0(counter_pair, :) = Pr_CH1.*Pr_CH2;
    prob_data.H1(counter_pair, :) = Pr_CH1_given_CH2.*Pr_CH2;
    prob_data.ratio(counter_pair, :) = ratio;
end
fprintf(' --> Completed. \n')

%% Plot
hFig = figure(1);
clf(hFig);

n_row = 3;
n_col = 3;
x_axis_ = idx_bin(1:end-1) - 300;

clearvars ax_
for counter_subplot = 1 : n_row*n_col
    ax_(counter_subplot) = subplot(n_row,n_col,counter_subplot);
end

subplot(n_row,n_col,1)
hold on
% plot(idx_bin, Pr_CH1, 'k')
plot(x_axis_, nanmean(prob_data.Pr_CH1, 1), 'k');
ylabel('Pr(CH1)')

subplot(n_row,n_col,4)
hold on
% plot(idx_bin, Pr_CH2, 'k')
plot(x_axis_, nanmean(prob_data.Pr_CH2, 1), 'k')
ylabel('Pr(CH2)')

subplot(n_row,n_col,7)
hold on
% plot(idx_bin, Pr_CH1_or_CH2, 'k')
plot(x_axis_, nanmean(prob_data.Pr_CH1_or_CH2, 1), 'k')
ylabel('Pr(CH1 or CH2)')

subplot(n_row,n_col,2)
hold on
% plot(idx_bin, Pr_CH1_given_CH2, 'k')
plot(x_axis_, nanmean(prob_data.Pr_CH1_given_CH2, 1), 'k')
ylabel('Pr(CH1 | CH2)')

subplot(n_row,n_col,5)
hold on
% plot(idx_bin, Pr_CH2_given_CH1, 'k')
plot(x_axis_, nanmean(prob_data.Pr_CH2_given_CH1, 1), 'k')
ylabel('Pr(CH2 | CH1)')

subplot(n_row,n_col,8)
hold on
% plot(idx_bin, Pr_CH1_and_CH2, 'k')
plot(x_axis_, nanmean(prob_data.Pr_CH1_and_CH2, 1), 'k')
ylabel('Pr(CH1, CH2)')

subplot(n_row,n_col,3)
hold on
% plot(idx_bin, Pr_CH1.*Pr_CH2, 'k')
plot(x_axis_, nanmean(prob_data.H0, 1), 'k')
ylabel('Pr(CH1) * Pr(CH2)')
title('H0')

subplot(n_row,n_col,6)
hold on
% plot(idx_bin, Pr_CH1_given_CH2.*Pr_CH2, 'k')
plot(x_axis_, nanmean(prob_data.H1, 1), 'k')
ylabel('Pr(CH1 | CH2) * Pr(CH2)')
title('H1')

subplot(n_row,n_col,9)
hold on
% ratio = Pr_CH1_given_CH2./Pr_CH1;
% ratio(isnan(ratio)) = 0;
% plot(idx_bin, ratio, 'k')
plot(x_axis_, nanmean(prob_data.ratio, 1), 'k')
ylabel('Pr(CH1 | CH2) / Pr(CH1)')
title('H1/H0')

% linkaxes([ax_(1) ax_(4) ax_(2) ax_(5)],'y')

% ESN_Beautify_Plot(hFig, [8 5], 8)
ESN_Beautify_Plot(hFig, [18 12], 12)
sgtitle(hFig, [spike_type ', ' event_type ', bin: ' num2str(bin_size) 'ms' ', ' 'dir: ' dir_type], 'Interpreter', 'none')

end

%% function build_PAIR_all_saccades
function ALL_PCELL_all_saccades = build_PAIR_all_saccades(pair_list_full, path_data_monkey_sorted)
%% init vars
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
pair_list_isstr = arrayfun(@iscellstr,pair_list_full);
num_pairs = size(pair_list_full, 1) / 2;
variable_list = {'primSac', 'corrSac', 'toTgt', 'fromTgt', 'lowAmp', 'toTgtStr', 'fromCenter', 'aroundCenter', 'nonTask', 'all', 'primCorrStr'};
Pr_list = {'Pr_CH1', 'Pr_CH2', 'Pr_CH1_and_CH2', 'Pr_CH1_or_CH2', 'Pr_CH1_given_CH2', 'Pr_CH2_given_CH1'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};
clearvars ALL_PCELL_COMPRESSED_DATA
%% Loop over pCells
for counter_pair = 1 : 1 : num_pairs
    fprintf(['Analyzing pair no. ', num2str(counter_pair), ' / ' num2str(num_pairs) ' ... \n']);
    num_recording = sum(pair_list_isstr(2*counter_pair-1, :));
    clearvars PAIR_DATA_recordings
    for counter_recording = 1 : 1 : num_recording
        %% build file_name_1 plot_data address
        file_name = pair_list_full{2*counter_pair-1, counter_recording};
        year_ = file_name(1:2);
        month_ = file_name(3:4);
        day_ = file_name(5:6);
        hour_ = file_name(8:9);
        minute_ = file_name(10:11);
        second_ = file_name(12:13);
        subFolder_month = ['20' year_ '-' month_ filesep];
        subFolder_day = ['20' year_ '-' month_ '-' day_ filesep];
        subFolder_recording = ['20' year_ '-' month_ '-' day_ '_' hour_ '-' minute_ '-' second_ filesep];
        subFolder_data = ['analyzed_data' filesep];
        file_path = [path_data_monkey_sorted subFolder_month subFolder_day subFolder_recording subFolder_data];
        file_name_1 = pair_list_full{2*counter_pair-1, counter_recording};
        file_name_2 = pair_list_full{2*counter_pair  , counter_recording};
        %% load plot_data
        PAIR_DATA_recording = build_PAIR_DATA(file_path, file_name_1, file_name_2);
        PAIR_DATA_recording.file_name_1 = file_name_1;
        PAIR_DATA_recording.file_name_2 = file_name_2;
        PAIR_DATA_recording.file_path = file_path;
        if PAIR_DATA_recording.file_name_1(18) == 's'
            PAIR_DATA_recording.id_1 = PAIR_DATA_recording.file_name_1(1:16);
        elseif PAIR_DATA_recording.file_name_1(18) == '2'
            PAIR_DATA_recording.id_1 = PAIR_DATA_recording.file_name_1(1:18);
        else
            error('Build plot_data_compress: cell id does not follow the standards')
        end
        if PAIR_DATA_recording.file_name_2(18) == 's'
            PAIR_DATA_recording.id_2 = PAIR_DATA_recording.file_name_2(1:16);
        elseif PAIR_DATA_recording.file_name_2(18) == '2'
            PAIR_DATA_recording.id_2 = PAIR_DATA_recording.file_name_2(1:18);
        else
            error('Build plot_data_compress: cell id does not follow the standards')
        end
        PAIR_DATA_recordings(counter_recording) = PAIR_DATA_recording;
    end
    %% 
    PAIR_DATA_session = struct;
    if num_recording > 1
        
        for counter_variable  = 1 : length(variable_list)
        for counter_indType   = 1 : length(indType_list)
        for counter_spikeType = 1 : length(spikeType_list)
            variable_name = variable_list{counter_variable};
            indType_name = indType_list{counter_indType};
            spikeType_name = spikeType_list{counter_spikeType};
            num_amp_bin = size(PAIR_DATA_recordings(1).(variable_name).([spikeType_name '_Pr_CH1_' indType_name]),1);
            num_ang_bin = size(PAIR_DATA_recordings(1).(variable_name).([spikeType_name '_Pr_CH1_' indType_name]),2);

            for counter_amp = 1 : num_amp_bin
            for counter_ang = 1 : num_ang_bin
                %%%%%%%%%%%%%%%%%%%%%%%%
                for counter_Pr = 1 : length(Pr_list)
                    Pr_name = Pr_list{counter_Pr};
                    Pr_data_recordings = zeros(size(PAIR_DATA_recordings(1).(variable_name).([spikeType_name '_' Pr_name '_' indType_name]){counter_amp, counter_ang}));
                    num_saccades_recordings = 0;
                    for counter_recording = 1 : 1 : num_recording
                        Pr_data = PAIR_DATA_recordings(counter_recording).(variable_name).([spikeType_name '_' Pr_name '_' indType_name]){counter_amp, counter_ang};
                        num_saccades = PAIR_DATA_recordings(counter_recording).(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang);
                        Pr_data_recordings = Pr_data_recordings + (Pr_data * num_saccades);
                        num_saccades_recordings = num_saccades_recordings + num_saccades;
                    end
                    if num_saccades_recordings == 0
                        PAIR_DATA_session.(variable_name).([spikeType_name '_' Pr_name '_' indType_name]){counter_amp, counter_ang} = Pr_data_recordings;
                    else
                        PAIR_DATA_session.(variable_name).([spikeType_name '_' Pr_name '_' indType_name]){counter_amp, counter_ang} = Pr_data_recordings ./ num_saccades_recordings;
                    end
                    PAIR_DATA_session.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades_recordings;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%
                Pr_CH1 = PAIR_DATA_session.(variable_name).([spikeType_name '_Pr_CH1_' indType_name]){counter_amp, counter_ang};
                Pr_CH2 = PAIR_DATA_session.(variable_name).([spikeType_name '_Pr_CH2_' indType_name]){counter_amp, counter_ang};
                Pr_CH1_given_CH2 = PAIR_DATA_session.(variable_name).([spikeType_name '_Pr_CH1_given_CH2_' indType_name]){counter_amp, counter_ang};
                H0 = Pr_CH1.*Pr_CH2;
                H1 = Pr_CH1_given_CH2.*Pr_CH2;
                ratio = Pr_CH1_given_CH2./Pr_CH1;
                ratio(isnan(ratio)) = 0;
                PAIR_DATA_session.(variable_name).([spikeType_name '_H0_' indType_name]){counter_amp, counter_ang} = H0;
                PAIR_DATA_session.(variable_name).([spikeType_name '_H1_' indType_name]){counter_amp, counter_ang} = H1;
                PAIR_DATA_session.(variable_name).([spikeType_name '_ratio_' indType_name]){counter_amp, counter_ang} = ratio;
            end
            end
        end
        end
        end
        % concatenate cell id
        PAIR_DATA_session.id_1 = cell(num_recording, 1);
        PAIR_DATA_session.id_2 = cell(num_recording, 1);
        PAIR_DATA_session.file_name_1 = cell(num_recording, 1);
        PAIR_DATA_session.file_name_2 = cell(num_recording, 1);
        PAIR_DATA_session.file_path = cell(num_recording, 1);
        for counter_recording = 1 : 1 : num_recording
            PAIR_DATA_session.id_1{counter_recording, 1} = PAIR_DATA_recordings(counter_recording).id_1;
            PAIR_DATA_session.id_2{counter_recording, 1} = PAIR_DATA_recordings(counter_recording).id_2;
            PAIR_DATA_session.file_name_1{counter_recording, 1} = PAIR_DATA_recordings(counter_recording).file_name_1;
            PAIR_DATA_session.file_name_2{counter_recording, 1} = PAIR_DATA_recordings(counter_recording).file_name_2;
            PAIR_DATA_session.file_path{counter_recording, 1} = PAIR_DATA_recordings(counter_recording).file_path;
        end
    else
        PAIR_DATA_session = PAIR_DATA_recordings;
    end
    
    %% Save PAIR_DATA_cell into ALL_PCELL_COMPRESSED_DATA
    ALL_PCELL_COMPRESSED_DATA(counter_pair) = PAIR_DATA_session;
    fprintf(' --> Completed. \n')
end

%% Loop over pairs
fprintf(['Building ALL_PCELL_all_saccades', ' ... ']);
clearvars ALL_PCELL_all_saccades
Pr_list = {'Pr_CH1', 'Pr_CH2', 'Pr_CH1_and_CH2', 'Pr_CH1_or_CH2', 'Pr_CH1_given_CH2', 'Pr_CH2_given_CH1', 'H0', 'H1', 'ratio'};
for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    num_amp_bin = size(ALL_PCELL_COMPRESSED_DATA(1).(variable_name).([spikeType_name '_Pr_CH1_' indType_name]),1);
    num_ang_bin = size(ALL_PCELL_COMPRESSED_DATA(1).(variable_name).([spikeType_name '_Pr_CH1_' indType_name]),2);

    for counter_amp = 1 : num_amp_bin
    for counter_ang = 1 : num_ang_bin
        for counter_pair = 1 : 1 : num_pairs
            %%%%%%%%%%%%%%%%%%%%%%%%
            for counter_Pr = 1 : length(Pr_list)
                Pr_name = Pr_list{counter_Pr};
                Pr_data = ALL_PCELL_COMPRESSED_DATA(counter_pair).(variable_name).([spikeType_name '_' Pr_name '_' indType_name]){counter_amp, counter_ang};
                ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_' Pr_name '_' indType_name]){counter_amp, counter_ang}(counter_pair,:) = Pr_data;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%
            num_saccades = ALL_PCELL_COMPRESSED_DATA(counter_pair).(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang);
            ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, counter_ang}(counter_pair,:) = num_saccades;
        end
    end
    end
end
end
end

fprintf(' --> Completed. \n')
end

%% function build_PAIR_DATA
function [PAIR_DATA] = build_PAIR_DATA(file_path, file_name_1, file_name_2)
%% Handle inputs
if nargin < 1
    [file_name_,~] = uigetfile([pwd filesep '*.psort'], 'Select psort file');
    [~,file_name_1,~] = fileparts(file_name_);
    [file_name_,file_path] = uigetfile([pwd filesep '*.psort'], 'Select psort file');
    [~,file_name_2,~] = fileparts(file_name_);
end

%% load EPHYS sorted DATA
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
file_name_1 = [file_name_1 '.psort'];
fprintf(['Loading ', file_name_1, ' ... ']);
DATA_PSORT_1 = Psort_read_psort([file_path file_name_1]);
EPHYS.CH_sorted_file_name_1 = file_name_1;
fprintf(' --> Completed. \n')
file_name_2 = [file_name_2 '.psort'];
fprintf(['Loading ', file_name_2, ' ... ']);
DATA_PSORT_2 = Psort_read_psort([file_path file_name_2]);
EPHYS.CH_sorted_file_name_2 = file_name_2;
fprintf(' --> Completed. \n')
EPHYS.CH_sorted_file_path = file_path;

%% load EPHYS EVENT DATA
file_name = [EPHYS.CH_sorted_file_name_1(1:13) '_EVE1_aligned.mat'];
% fprintf(['Loading ', file_name, ' ... ']);
EPHYS.CH_EVE = load([file_path file_name]);
if isfield(EPHYS.CH_EVE, 'EPHYS_time_15K')
    EPHYS.CH_EVE.EPHYS_time_30K = EPHYS.CH_EVE.EPHYS_time_15K(:);
else
    EPHYS.CH_EVE.EPHYS_time_30K = EPHYS.CH_EVE.EPHYS_time_30K(:);
end
EPHYS.CH_EVE.EPHYS_time_1K  = reshape(EPHYS.CH_EVE.EPHYS_time_1K ,[], 1);
EPHYS.CH_EVE.BEHAVE_time_1K = reshape(EPHYS.CH_EVE.BEHAVE_time_1K,[], 1);
% fprintf(' --> Completed. \n')

%% load BEHAVE DATA
file_name = [EPHYS.CH_sorted_file_name_1(1:13) '_ANALYZED.mat'];
% fprintf(['Loading ', file_name, ' ... ']);
BEHAVE = load([file_path file_name]);
% fprintf(' --> Completed. \n')

%% build EPHYS.CH_sorted_1 from DATA_PSORT_1
ch_data = double(DATA_PSORT_1.topLevel_data.ch_data);
ch_time = double(DATA_PSORT_1.topLevel_data.ch_time);
SS_index = find(logical(double(DATA_PSORT_1.topLevel_data.ss_index)));
CS_index = find(logical(double(DATA_PSORT_1.topLevel_data.cs_index)));
SS_time_1 = ch_time(SS_index);
CS_time_1 = ch_time(CS_index);

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

EPHYS.CH_sorted_1.SS_data.SS_ind = SS_index;
EPHYS.CH_sorted_1.CS_data.CS_ind = CS_index;
EPHYS.CH_sorted_1.SS_data.SS_time = SS_time_1;
EPHYS.CH_sorted_1.CS_data.CS_time = CS_time_1;
EPHYS.CH_sorted_1.SS_data.SS_waveform = SS_waveform;
EPHYS.CH_sorted_1.CS_data.CS_waveform = CS_waveform;

%% build EPHYS.CH_sorted from DATA_PSORT
ch_data = double(DATA_PSORT_2.topLevel_data.ch_data);
ch_time = double(DATA_PSORT_2.topLevel_data.ch_time);
SS_index = find(logical(double(DATA_PSORT_2.topLevel_data.ss_index)));
CS_index = find(logical(double(DATA_PSORT_2.topLevel_data.cs_index)));
SS_time_1 = ch_time(SS_index);
CS_time_1 = ch_time(CS_index);

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

EPHYS.CH_sorted_2.SS_data.SS_ind = SS_index;
EPHYS.CH_sorted_2.CS_data.CS_ind = CS_index;
EPHYS.CH_sorted_2.SS_data.SS_time = SS_time_1;
EPHYS.CH_sorted_2.CS_data.CS_time = CS_time_1;
EPHYS.CH_sorted_2.SS_data.SS_waveform = SS_waveform;
EPHYS.CH_sorted_2.CS_data.CS_waveform = CS_waveform;

%% SS & CS train_aligned
clearvars -except EPHYS BEHAVE PAIR_DATA
% fprintf(['Building CS & SS train_aligned', ' ... ']);
EPHYS_time_1K     = EPHYS.CH_EVE.EPHYS_time_1K; % EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_time_1K; % 
length_time_ = length(EPHYS_time_1K);
CS_time_1 = EPHYS.CH_sorted_1.CS_data.CS_time;
if isempty(CS_time_1)
    CS_time_1 = EPHYS_time_1K(1);
end
CS_time_1(end+1) = max([EPHYS_time_1K(end), CS_time_1(end)])+1;
CS_time_2 = EPHYS.CH_sorted_2.CS_data.CS_time;
if isempty(CS_time_2)
    CS_time_2 = EPHYS_time_1K(1);
end
CS_time_2(end+1) = max([EPHYS_time_1K(end), CS_time_2(end)])+1;
SS_time_1 = EPHYS.CH_sorted_1.SS_data.SS_time;
if isempty(SS_time_1)
    SS_time_1 = EPHYS_time_1K(1);
end
SS_time_1(end+1) = max([EPHYS_time_1K(end), SS_time_1(end)])+1;
SS_time_2 = EPHYS.CH_sorted_2.SS_data.SS_time;
if isempty(SS_time_2)
    SS_time_2 = EPHYS_time_1K(1);
end
SS_time_2(end+1) = max([EPHYS_time_1K(end), SS_time_2(end)])+1;
EPHYS_CS_train_1K_1 = false(size(EPHYS_time_1K));
EPHYS_SS_train_1K_1 = false(size(EPHYS_time_1K));
EPHYS_CS_train_1K_2 = false(size(EPHYS_time_1K));
EPHYS_SS_train_1K_2 = false(size(EPHYS_time_1K));
counter_CS_1 = find(CS_time_1 >= EPHYS_time_1K(1), 1, 'first');
counter_SS_1 = find(SS_time_1 >= EPHYS_time_1K(1), 1, 'first');
counter_CS_2 = find(CS_time_2 >= EPHYS_time_1K(1), 1, 'first');
counter_SS_2 = find(SS_time_2 >= EPHYS_time_1K(1), 1, 'first');
for counter_time_point = 1 : length_time_
    time_ponit_ = EPHYS_time_1K(counter_time_point);
    if time_ponit_>=CS_time_1(counter_CS_1)
        EPHYS_CS_train_1K_1(counter_time_point) = true;
        counter_CS_1 = counter_CS_1 + 1;
    end
    if time_ponit_>=SS_time_1(counter_SS_1)
        EPHYS_SS_train_1K_1(counter_time_point) = true;
        counter_SS_1 = counter_SS_1 + 1;
    end
    if time_ponit_>=CS_time_2(counter_CS_2)
        EPHYS_CS_train_1K_2(counter_time_point) = true;
        counter_CS_2 = counter_CS_2 + 1;
    end
    if time_ponit_>=SS_time_2(counter_SS_2)
        EPHYS_SS_train_1K_2(counter_time_point) = true;
        counter_SS_2 = counter_SS_2 + 1;
    end
end
EPHYS.CH_EVE.EPHYS_CS_train_1K_1 = EPHYS_CS_train_1K_1;
EPHYS.CH_EVE.EPHYS_SS_train_1K_1 = EPHYS_SS_train_1K_1;
EPHYS.CH_EVE.EPHYS_CS_train_1K_2 = EPHYS_CS_train_1K_2;
EPHYS.CH_EVE.EPHYS_SS_train_1K_2 = EPHYS_SS_train_1K_2;

% fprintf(' --> Completed. \n')

%% Building PAIR_DATA
clearvars -except EPHYS BEHAVE PAIR_DATA
fprintf(['Building PAIR_DATA', ' ... ']);
PAIR_DATA = struct;
inds_span    = ((-300+1) : 1 : (300))';
amp_edges = [-0.5 2 4 6 8 10 50];
vel_edges = [0 200 300 400 500 600 10000];
ang_edges = -202.5 : 45 : +202.5;
PAIR_DATA.inds_span = inds_span;
PAIR_DATA.amp_edges = amp_edges;
PAIR_DATA.vel_edges = vel_edges;
PAIR_DATA.ang_edges = ang_edges;

PAIR_DATA.use_vel_instead_of_amp = false;

length_train_data_ = length(EPHYS.CH_EVE.EPHYS_time_1K);
length_velocity_data_ = length(BEHAVE.aligned.time_1K);
velocity_trace_data = BEHAVE.aligned.eye_r_vm_filt;

ang_edges_last_bin_id = length(ang_edges) - 1;

[~, ~, amp_bin] = histcounts(BEHAVE.SACS_ALL_DATA.eye_r_amp_m, amp_edges);
[~, ~, vel_bin] = histcounts(BEHAVE.SACS_ALL_DATA.eye_r_vm_max, vel_edges);
[~, ~, ang_bin] = histcounts(BEHAVE.SACS_ALL_DATA.eye_r_ang, ang_edges);
ang_bin(ang_bin == ang_edges_last_bin_id) = 1;

num_amp_bin = length(amp_edges) - 1;
num_vel_bin = length(vel_edges) - 1;
num_ang_bin = length(ang_edges) - 2;

if PAIR_DATA.use_vel_instead_of_amp
    amp_bin = vel_bin;
    num_amp_bin = num_vel_bin;
end

BEHAVE.SACS_ALL_DATA.is_all = BEHAVE.SACS_ALL_DATA.validity;
BEHAVE.SACS_ALL_DATA.is_primCorrStr = BEHAVE.SACS_ALL_DATA.is_primSac | BEHAVE.SACS_ALL_DATA.is_corrSac | BEHAVE.SACS_ALL_DATA.is_toTgtStr;

variable_list = {'primSac', 'corrSac', 'toTgt', 'fromTgt', 'lowAmp', 'toTgtStr', 'fromCenter', 'aroundCenter', 'nonTask', 'all', 'primCorrStr'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};
bin_size = 5;
idx_bin = 1 : bin_size : length(PAIR_DATA.inds_span)+1;
% Main for loop
for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    
    PAIR_DATA.(variable_name).([spikeType_name '_Pr_CH1_' indType_name]) = cell(num_amp_bin, num_ang_bin);
    PAIR_DATA.(variable_name).([spikeType_name '_Pr_CH2_' indType_name]) = cell(num_amp_bin, num_ang_bin);
    PAIR_DATA.(variable_name).([spikeType_name '_Pr_CH1_and_CH2_'   indType_name]) = cell(num_amp_bin, num_ang_bin);
    PAIR_DATA.(variable_name).([spikeType_name '_Pr_CH1_or_CH2_'    indType_name]) = cell(num_amp_bin, num_ang_bin);
    PAIR_DATA.(variable_name).([spikeType_name '_Pr_CH1_given_CH2_' indType_name]) = cell(num_amp_bin, num_ang_bin);
    PAIR_DATA.(variable_name).([spikeType_name '_Pr_CH2_given_CH1_' indType_name]) = cell(num_amp_bin, num_ang_bin);
    PAIR_DATA.(variable_name).([spikeType_name '_H0_'    indType_name]) = cell(num_amp_bin, num_ang_bin);
    PAIR_DATA.(variable_name).([spikeType_name '_H1_'    indType_name]) = cell(num_amp_bin, num_ang_bin);
    PAIR_DATA.(variable_name).([spikeType_name '_ratio_' indType_name]) = cell(num_amp_bin, num_ang_bin);
    PAIR_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name]) = zeros(num_amp_bin, num_ang_bin);
    is_variable_name = BEHAVE.SACS_ALL_DATA.(['is_' variable_name]);
    spike_train_data_1 = EPHYS.CH_EVE.(['EPHYS_' spikeType_name '_train_1K_1']);
    spike_train_data_2 = EPHYS.CH_EVE.(['EPHYS_' spikeType_name '_train_1K_2']);
    ind_indType_name = BEHAVE.SACS_ALL_DATA.(['ind_' indType_name]);
    ind_converted = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(ind_indType_name);
    
    for counter_amp = 1 : num_amp_bin
    for counter_ang = 1 : num_ang_bin
        is_ang_bin   = (ang_bin == counter_ang);
        is_amp_bin   = (amp_bin == counter_amp);
        ind_train    = ind_converted(is_variable_name & is_ang_bin & is_amp_bin);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(ind_train)
            raster_CH1 = zeros(1, length(inds_span));
            raster_CH2 = zeros(1, length(inds_span));
            num_saccades = 0;
        else
            inds_train = repmat( ind_train(:), 1, length(inds_span)) + repmat(inds_span(:)', length(ind_train), 1);
            inds_train( inds_train < 1 ) = 1;
            inds_train( inds_train > length_train_data_ ) = length_train_data_;
            raster_CH1 = spike_train_data_1(inds_train);
            raster_CH2 = spike_train_data_2(inds_train);
            if length(ind_train) == 1
                raster_CH1 = reshape( raster_CH1 , 1, []);
                raster_CH2 = reshape( raster_CH2 , 1, []);
            end
            num_saccades = length(ind_train);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        num_trials = min([size(raster_CH1, 1) size(raster_CH2, 1)]);
        raster_CH1 = logical(raster_CH1(1:num_trials, :));
        raster_CH2 = logical(raster_CH2(1:num_trials, :));
        raster_CH1_bin = false(num_trials, length(idx_bin)-1);
        raster_CH2_bin = false(num_trials, length(idx_bin)-1);
        for counter_idx = 1 : 1 : length(idx_bin)-1
            idx_i = idx_bin(counter_idx);
            idx_f = idx_bin(counter_idx+1)-1;
            raster_CH1_bin(:, counter_idx) = logical(sum(raster_CH1(:,idx_i:idx_f),2));
            raster_CH2_bin(:, counter_idx) = logical(sum(raster_CH2(:,idx_i:idx_f),2));
        end
        if num_trials == 1
            Pr_CH1 = double(raster_CH1_bin) ./ num_trials;
            Pr_CH2 = double(raster_CH2_bin) ./ num_trials;
            Pr_CH1_and_CH2 = double(raster_CH1_bin & raster_CH2_bin) ./ num_trials;
            Pr_CH1_or_CH2  = double(raster_CH1_bin | raster_CH2_bin) ./ num_trials;
        else
            Pr_CH1 = sum(raster_CH1_bin) ./ num_trials;
            Pr_CH2 = sum(raster_CH2_bin) ./ num_trials;
            Pr_CH1_and_CH2 = sum(raster_CH1_bin & raster_CH2_bin) ./ num_trials;
            Pr_CH1_or_CH2  = sum(raster_CH1_bin | raster_CH2_bin) ./ num_trials;
        end
        Pr_CH1_given_CH2 = zeros(size(Pr_CH1));
        Pr_CH2_given_CH1 = zeros(size(Pr_CH2));
        for counter_col = 1 : 1 : length(Pr_CH1)
            CS_CH1_bin_col = raster_CH1_bin(:,counter_col);
            CS_CH2_bin_col = raster_CH2_bin(:,counter_col);
            Pr_CH1_given_CH2(counter_col) = sum(CS_CH1_bin_col(CS_CH2_bin_col)) ./ sum(CS_CH2_bin_col);
            Pr_CH2_given_CH1(counter_col) = sum(CS_CH2_bin_col(CS_CH1_bin_col)) ./ sum(CS_CH1_bin_col);
        end
        Pr_CH1_given_CH2(isnan(Pr_CH1_given_CH2)) = 0;
        Pr_CH2_given_CH1(isnan(Pr_CH2_given_CH1)) = 0;
        H0 = Pr_CH1.*Pr_CH2;
        H1 = Pr_CH1_given_CH2.*Pr_CH2;
        ratio = Pr_CH1_given_CH2./Pr_CH1;
        ratio(isnan(ratio)) = 0;
        
        PAIR_DATA.(variable_name).([spikeType_name '_Pr_CH1_' indType_name]){counter_amp, counter_ang} = Pr_CH1;
        PAIR_DATA.(variable_name).([spikeType_name '_Pr_CH2_' indType_name]){counter_amp, counter_ang} = Pr_CH2;
        PAIR_DATA.(variable_name).([spikeType_name '_Pr_CH1_and_CH2_'   indType_name]){counter_amp, counter_ang} = Pr_CH1_and_CH2;
        PAIR_DATA.(variable_name).([spikeType_name '_Pr_CH1_or_CH2_'    indType_name]){counter_amp, counter_ang} = Pr_CH1_or_CH2;
        PAIR_DATA.(variable_name).([spikeType_name '_Pr_CH1_given_CH2_' indType_name]){counter_amp, counter_ang} = Pr_CH1_given_CH2;
        PAIR_DATA.(variable_name).([spikeType_name '_Pr_CH2_given_CH1_' indType_name]){counter_amp, counter_ang} = Pr_CH2_given_CH1;
        PAIR_DATA.(variable_name).([spikeType_name '_H0_'    indType_name]){counter_amp, counter_ang} = H0;
        PAIR_DATA.(variable_name).([spikeType_name '_H1_'    indType_name]){counter_amp, counter_ang} = H1;
        PAIR_DATA.(variable_name).([spikeType_name '_ratio_' indType_name]){counter_amp, counter_ang} = ratio;
        PAIR_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
    end
    end
end
end
end

% Remove some fields
rmfields_list = {'inds_span', 'amp_edges', 'ang_edges', 'vel_edges', 'use_vel_instead_of_amp'};
PAIR_DATA = rmfield(PAIR_DATA,rmfields_list);
fprintf(' --> Completed. \n')

end

%% function pair_data_avg_over_amplitude
function PAIR_DATA_all_saccades_ang = pair_data_avg_over_amplitude(PAIR_DATA_all_saccades)
%% Loop over pCells
fprintf(['Building PAIR_DATA_all_saccades_ang', ' ... ']);
clearvars PAIR_DATA_all_saccades_ang
variable_list = {'primSac', 'corrSac', 'toTgt', 'fromTgt', 'lowAmp', 'toTgtStr', 'fromCenter', 'aroundCenter', 'nonTask', 'all', 'primCorrStr'};
Pr_list = {'Pr_CH1', 'Pr_CH2', 'Pr_CH1_and_CH2', 'Pr_CH1_or_CH2', 'Pr_CH1_given_CH2', 'Pr_CH2_given_CH1'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};
num_pCells = size(PAIR_DATA_all_saccades.primSac.SS_Pr_CH1_start{1,1}, 1);
span_width = size(PAIR_DATA_all_saccades.primSac.SS_Pr_CH1_start{1,1}, 2);
for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    num_amp_bin = size(PAIR_DATA_all_saccades.(variable_name).([spikeType_name '_Pr_CH1_' indType_name]),1);
    num_ang_bin = size(PAIR_DATA_all_saccades.(variable_name).([spikeType_name '_Pr_CH1_' indType_name]),2);
    
    for counter_ang = 1 : num_ang_bin
        %%%%%%%%%%%%%%%%%%%%%%%%
        for counter_Pr = 1 : length(Pr_list)
            Pr_name = Pr_list{counter_Pr};
            Pr_data_all_ang = zeros(size(PAIR_DATA_all_saccades.(variable_name).([spikeType_name '_' Pr_name '_' indType_name]){1, counter_ang}));
            num_saccades_all_ang = zeros(size(PAIR_DATA_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){1, counter_ang}));
            if size(num_saccades_all_ang, 2) == 1
                num_saccades_all_ang = repmat(num_saccades_all_ang, 1, span_width);
            end
            for counter_amp = 1 : num_amp_bin
                Pr_data = PAIR_DATA_all_saccades.(variable_name).([spikeType_name '_' Pr_name '_' indType_name]){counter_amp, counter_ang};
                num_saccades = PAIR_DATA_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, counter_ang};
                if size(num_saccades, 2) == 1
                    num_saccades = repmat(num_saccades, 1, span_width);
                end
                Pr_data_all_ang = Pr_data_all_ang + (Pr_data .* num_saccades);
                num_saccades_all_ang = num_saccades_all_ang + num_saccades;
            end
            Pr_data_all_ang = Pr_data_all_ang ./ num_saccades_all_ang;
            Pr_data_all_ang(isnan(Pr_data_all_ang)) = 0;
            PAIR_DATA_all_saccades_ang.(variable_name).([spikeType_name '_' Pr_name '_' indType_name]){1, counter_ang} = Pr_data_all_ang;
            PAIR_DATA_all_saccades_ang.(variable_name).([spikeType_name '_num_sac_' indType_name]){1, counter_ang} = num_saccades_all_ang;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
        Pr_CH1 = PAIR_DATA_all_saccades_ang.(variable_name).([spikeType_name '_Pr_CH1_' indType_name]){1, counter_ang};
        Pr_CH2 = PAIR_DATA_all_saccades_ang.(variable_name).([spikeType_name '_Pr_CH2_' indType_name]){1, counter_ang};
        Pr_CH1_given_CH2 = PAIR_DATA_all_saccades_ang.(variable_name).([spikeType_name '_Pr_CH1_given_CH2_' indType_name]){1, counter_ang};
        H0 = Pr_CH1.*Pr_CH2;
        H1 = Pr_CH1_given_CH2.*Pr_CH2;
        ratio = Pr_CH1_given_CH2./Pr_CH1;
        ratio(isnan(ratio)) = 0;
        PAIR_DATA_all_saccades_ang.(variable_name).([spikeType_name '_H0_' indType_name]){1, counter_ang} = H0;
        PAIR_DATA_all_saccades_ang.(variable_name).([spikeType_name '_H1_' indType_name]){1, counter_ang} = H1;
        PAIR_DATA_all_saccades_ang.(variable_name).([spikeType_name '_ratio_' indType_name]){1, counter_ang} = ratio;
    end
end
end
end
fprintf(' --> Completed. \n')
end

%% function pair_data_avg_over_angle
function PAIR_DATA_all_saccades_amp = pair_data_avg_over_angle(PAIR_DATA_all_saccades)
%% Loop over pCells
fprintf(['Building PAIR_DATA_all_saccades_amp', ' ... ']);
clearvars PAIR_DATA_all_saccades_amp
variable_list = {'primSac', 'corrSac', 'toTgt', 'fromTgt', 'lowAmp', 'toTgtStr', 'fromCenter', 'aroundCenter', 'nonTask', 'all', 'primCorrStr'};
Pr_list = {'Pr_CH1', 'Pr_CH2', 'Pr_CH1_and_CH2', 'Pr_CH1_or_CH2', 'Pr_CH1_given_CH2', 'Pr_CH2_given_CH1'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};
num_pCells = size(PAIR_DATA_all_saccades.primSac.SS_Pr_CH1_start{1,1}, 1);
span_width = size(PAIR_DATA_all_saccades.primSac.SS_Pr_CH1_start{1,1}, 2);
for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    num_amp_bin = size(PAIR_DATA_all_saccades.(variable_name).([spikeType_name '_Pr_CH1_' indType_name]),1);
    num_ang_bin = size(PAIR_DATA_all_saccades.(variable_name).([spikeType_name '_Pr_CH1_' indType_name]),2);
    
    for counter_amp = 1 : num_amp_bin
        %%%%%%%%%%%%%%%%%%%%%%%%
        for counter_Pr = 1 : length(Pr_list)
            Pr_name = Pr_list{counter_Pr};
            Pr_data_all_ang = zeros(size(PAIR_DATA_all_saccades.(variable_name).([spikeType_name '_' Pr_name '_' indType_name]){counter_amp, 1}));
            num_saccades_all_ang = zeros(size(PAIR_DATA_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, 1}));
            if size(num_saccades_all_ang, 2) == 1
                num_saccades_all_ang = repmat(num_saccades_all_ang, 1, span_width);
            end
            for counter_ang = 1 : num_ang_bin
                Pr_data = PAIR_DATA_all_saccades.(variable_name).([spikeType_name '_' Pr_name '_' indType_name]){counter_amp, counter_ang};
                num_saccades = PAIR_DATA_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, counter_ang};
                if size(num_saccades, 2) == 1
                    num_saccades = repmat(num_saccades, 1, span_width);
                end
                Pr_data_all_ang = Pr_data_all_ang + (Pr_data .* num_saccades);
                num_saccades_all_ang = num_saccades_all_ang + num_saccades;
            end
            Pr_data_all_ang = Pr_data_all_ang ./ num_saccades_all_ang;
            Pr_data_all_ang(isnan(Pr_data_all_ang)) = 0;
            PAIR_DATA_all_saccades_amp.(variable_name).([spikeType_name '_' Pr_name '_' indType_name]){counter_amp, 1} = Pr_data_all_ang;
            PAIR_DATA_all_saccades_amp.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, 1} = num_saccades_all_ang;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
        Pr_CH1 = PAIR_DATA_all_saccades_amp.(variable_name).([spikeType_name '_Pr_CH1_' indType_name]){counter_amp, 1};
        Pr_CH2 = PAIR_DATA_all_saccades_amp.(variable_name).([spikeType_name '_Pr_CH2_' indType_name]){counter_amp, 1};
        Pr_CH1_given_CH2 = PAIR_DATA_all_saccades_amp.(variable_name).([spikeType_name '_Pr_CH1_given_CH2_' indType_name]){counter_amp, 1};
        H0 = Pr_CH1.*Pr_CH2;
        H1 = Pr_CH1_given_CH2.*Pr_CH2;
        ratio = Pr_CH1_given_CH2./Pr_CH1;
        ratio(isnan(ratio)) = 0;
        PAIR_DATA_all_saccades_amp.(variable_name).([spikeType_name '_H0_' indType_name]){counter_amp, 1} = H0;
        PAIR_DATA_all_saccades_amp.(variable_name).([spikeType_name '_H1_' indType_name]){counter_amp, 1} = H1;
        PAIR_DATA_all_saccades_amp.(variable_name).([spikeType_name '_ratio_' indType_name]){counter_amp, 1} = ratio;
    end
end
end
end
fprintf(' --> Completed. \n')
end

%% function build_PAIR_DATA_all_saccades_tuned
function PAIR_DATA_all_saccades_tuned = build_PAIR_DATA_all_saccades_tuned(PAIR_DATA_all_saccades)
%% find CS-on for each pair
PAIR_DATA_all_saccades_ang = pair_data_avg_over_amplitude(PAIR_DATA_all_saccades);
num_ang_bin = size(PAIR_DATA_all_saccades_ang.primSac.CS_Pr_CH1_start,2);
num_pairs   = size(PAIR_DATA_all_saccades_ang.primSac.CS_Pr_CH1_start{1,1},1);
%%
span_length_ = size(PAIR_DATA_all_saccades_ang.primSac.CS_Pr_CH1_start{1,1},2);
range_before_sac = round(span_length_/6) : 3*round(span_length_/6);
CS_prob_primSac_1 = zeros(num_pairs, num_ang_bin);
CS_prob_corrSac_1 = zeros(num_pairs, num_ang_bin);
CS_prob_primSac_2 = zeros(num_pairs, num_ang_bin);
CS_prob_corrSac_2 = zeros(num_pairs, num_ang_bin);
for counter_ang = 1 : num_ang_bin
    CS_prob_primSac_1(:,counter_ang) = nansum(PAIR_DATA_all_saccades_ang.primSac.CS_Pr_CH1_start{1,counter_ang}(:,range_before_sac), 2);
    CS_prob_corrSac_1(:,counter_ang) = nansum(PAIR_DATA_all_saccades_ang.corrSac.CS_Pr_CH1_start{1,counter_ang}(:,range_before_sac), 2);
    CS_prob_primSac_2(:,counter_ang) = nansum(PAIR_DATA_all_saccades_ang.primSac.CS_Pr_CH2_start{1,counter_ang}(:,range_before_sac), 2);
    CS_prob_corrSac_2(:,counter_ang) = nansum(PAIR_DATA_all_saccades_ang.corrSac.CS_Pr_CH2_start{1,counter_ang}(:,range_before_sac), 2);
end

CS_prob_primSac_1(1:16, [2 4 6 8]) = nan;
CS_prob_corrSac_1(1:16, [2 4 6 8]) = nan;
CS_prob_primSac_2(1:16, [2 4 6 8]) = nan;
CS_prob_corrSac_2(1:16, [2 4 6 8]) = nan;
CS_prob_avg = (CS_prob_primSac_1 + CS_prob_corrSac_1 + CS_prob_primSac_2 + CS_prob_corrSac_2) ./ 4;
% [~, idx_CS_max] = nanmax(CS_prob_avg, [], 2);
dir_angles = deg2rad([180, 225, 270, 315, 0, 45, 90, 135]);
% compute weighted sum of cos and sin of angles
r = nansum(CS_prob_avg.* repmat(exp(1i*dir_angles), num_pairs, 1) , 2);
% Computes the mean direction for circular data.
CS_ang_avg = rad2deg(angle(r));
% sum of weights
CS_prob_sum = nansum(CS_prob_avg,2);
% Computes mean resultant vector length for circular data.
CS_rho_avg = abs(r) ./ CS_prob_sum;

ang_edges = -202.5 : 45 : +202.5;
ang_edges_last_bin_id = length(ang_edges) - 1;
[~, ~, idx_CS_max] = histcounts(CS_ang_avg, ang_edges);
idx_CS_max(idx_CS_max == ang_edges_last_bin_id) = 1;

hFig = figure(1); clf(hFig);
step_size_ = pi/8;
ang_edges = -pi-(step_size_/2):step_size_:pi-(step_size_/2);
polarhistogram(deg2rad(CS_ang_avg), ang_edges, 'DisplayStyle', 'bar','FaceColor','red', 'EdgeColor', 'red')
set(gca, 'ThetaTick', 0:45:315, 'RTick', 0:5:25)
rlim([0 5])
%}
%%
idx_CS_tuning = zeros(num_pairs, num_ang_bin);
for counter_pair = 1 : 1 : num_pairs
    CS_on_index = idx_CS_max(counter_pair);
    CS_on_index = CS_on_index - 1; % CS_on_index should be in 0-index format
    if (CS_on_index == 8); CS_on_index = 0; end
    idx_CS_tuning(counter_pair, :) = mod((CS_on_index : 1 : CS_on_index+7), 8) + 1;
end

%% build PAIR_DATA_all_saccades_tuned
fprintf(['Building PAIR_DATA_all_saccades_tuned', ' ... ']);
variable_list = {'primSac', 'corrSac', 'toTgt', 'fromTgt', 'lowAmp', 'toTgtStr', 'fromCenter', 'aroundCenter', 'nonTask', 'all', 'primCorrStr'};
Pr_list = {'Pr_CH1', 'Pr_CH2', 'Pr_CH1_and_CH2', 'Pr_CH1_or_CH2', 'Pr_CH1_given_CH2', 'Pr_CH2_given_CH1', 'H0', 'H1', 'ratio'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};
clearvars PAIR_DATA_all_saccades_tuned
for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    num_amp_bin = size(PAIR_DATA_all_saccades.(variable_name).([spikeType_name '_Pr_CH1_' indType_name]),1);
    num_ang_bin = size(PAIR_DATA_all_saccades.(variable_name).([spikeType_name '_Pr_CH1_' indType_name]),2);

    for counter_amp = 1 : num_amp_bin
    for counter_ang = 1 : num_ang_bin
        for counter_pair = 1 : 1 : num_pairs
            idx_ang_ = idx_CS_tuning(counter_pair, counter_ang);
            %%%%%%%%%%%%%%%%%%%%%%%%
            for counter_Pr = 1 : length(Pr_list)
                Pr_name = Pr_list{counter_Pr};
                Pr_data = PAIR_DATA_all_saccades.(variable_name).([spikeType_name '_' Pr_name '_' indType_name]){counter_amp, idx_ang_}(counter_pair,:);
                PAIR_DATA_all_saccades_tuned.(variable_name).([spikeType_name '_' Pr_name '_' indType_name]){counter_amp, counter_ang}(counter_pair,:) = Pr_data;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%
            num_saccades = PAIR_DATA_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, idx_ang_}(counter_pair,:);
            PAIR_DATA_all_saccades_tuned.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, counter_ang}(counter_pair,:) = num_saccades;
        end
    end
    end
end
end
end
PAIR_DATA_all_saccades_tuned.idx_CS_tuning = idx_CS_tuning;
fprintf(' --> Completed. \n')

end

%% function scratch_plot_pair_data
function scratch_plot_pair_data
%%
hFig = figure(5);
clf(hFig)

variable_name_category = 'all'; % 'primSac'; % 'nonTask'; % 'toTgtStr'; % 'corrSac'; % 'fromCenter'; % 
event_type = 'start'; % 'vmax'; % 'finish'; % 
spike_type = 'CS';
variable_name_H1    = [spike_type '_H1_' event_type];
variable_name_H0    = [spike_type '_H0_' event_type];
variable_name_ratio = [spike_type '_ratio_' event_type];
variable_name_num_sac  = [spike_type '_num_sac_' event_type];
ind_amp = 1;
ind_ang = 1;
span_length_ = size(PAIR_DATA_all_saccades.primSac.CS_Pr_CH1_start{1,1},2);
span_length_6th = round(span_length_/6);
if contains(spike_type,'CS')
    range_ = 1 : span_length_;
    range_XTick = 0:(2*span_length_6th):(span_length_);
else
%     range_ = 2*span_length_6th : 4*span_length_6th;
%     range_XTick = 0:(1*span_length_6th):(2*span_length_6th);
%     range_ = 1 : span_length_;
%     range_XTick = 0:(2*span_length_6th):(span_length_);
    range_ = 1*span_length_6th : 5*span_length_6th;
    range_XTick = 0:(2*span_length_6th):(4*span_length_6th);
end
num_pairs = 35;
idx_pairs =  1:35;

ALL_PCELL_DATA_1 = PAIR_DATA_all_saccades_ang_amp; % PAIR_DATA_all_saccades_tuned_ang; % PAIR_DATA_all_saccades; % PAIR_DATA_all_saccades_tuned; % 

num_sac = ALL_PCELL_DATA_1.(variable_name_category).(variable_name_num_sac){ind_amp, ind_ang}(idx_pairs, :);
H1_data = ALL_PCELL_DATA_1.(variable_name_category).(variable_name_H1){ind_amp, ind_ang}(idx_pairs, :);
H0_data = ALL_PCELL_DATA_1.(variable_name_category).(variable_name_H0){ind_amp, ind_ang}(idx_pairs, :);
ratio_data = ALL_PCELL_DATA_1.(variable_name_category).(variable_name_ratio){ind_amp, ind_ang}(idx_pairs, :);

num_perm = 5000;
num_samples = num_pairs ;
length_span = size(H1_data, 2);
H1_data_perm = nan(num_perm, length_span);
H0_data_perm = nan(num_perm, length_span);
ratio_data_perm = nan(num_perm, length_span);
for counter_perm = 1 : num_perm
    inds_perm = randi(num_pairs,[num_samples 1]);
    H1_data_ = H1_data(inds_perm, :);
    H1_data_(isnan(H1_data_)) = 0;
    H0_data_ = H0_data(inds_perm, :);
    H0_data_(isnan(H0_data_)) = 0;
    ratio_data_ = ratio_data(inds_perm, :);
    ratio_data_(isnan(ratio_data_)) = 0;
    num_sac_ = num_sac(inds_perm, 1);
    H1_mean_ = num_sac_' * H1_data_ ./ nansum(num_sac_);
    H0_mean_ = num_sac_' * H0_data_ ./ nansum(num_sac_);
    ratio_mean_ = num_sac_' * ratio_data_ ./ nansum(num_sac_);
    H1_data_perm(counter_perm, :) = H1_mean_;
    H0_data_perm(counter_perm, :) = H0_mean_;
    ratio_data_perm(counter_perm, :) = ratio_mean_;
end
H1_mean = ESN_smooth( nanmean(H1_data_perm, 1) );
H1_stdv = ESN_smooth( nanstd( H1_data_perm, 0, 1) ) ./ sqrt(num_samples) * 3; % 95% CI
H1_mean = H1_mean(:, range_);
H1_stdv = H1_stdv(:, range_);

H0_mean = ESN_smooth( nanmean(H0_data_perm, 1) );
H0_stdv = ESN_smooth( nanstd( H0_data_perm, 0, 1) ) ./ sqrt(num_samples) * 3; % 95% CI
H0_mean = H0_mean(:, range_);
H0_stdv = H0_stdv(:, range_);

ratio_mean = ESN_smooth( nanmean(ratio_data_perm, 1) );
ratio_stdv = ESN_smooth( nanstd( ratio_data_perm, 0, 1) ) ./ sqrt(num_samples) * 3; % 95% CI
ratio_mean = ratio_mean(:, range_);
ratio_stdv = ratio_stdv(:, range_);

y_axis_H1_stdv = [(H1_mean + H1_stdv) flip(H1_mean - H1_stdv)];
x_axis_H1_stdv = [(1:1:length(H1_mean)) (length(H1_mean):-1:1)];
y_axis_H0_stdv = [(H0_mean + H0_stdv) flip(H0_mean - H0_stdv)];
x_axis_H0_stdv = [(1:1:length(H0_mean)) (length(H0_mean):-1:1)];
y_axis_ratio_stdv = [(ratio_mean + ratio_stdv) flip(ratio_mean - ratio_stdv)];
x_axis_ratio_stdv = [(1:1:length(ratio_mean)) (length(ratio_mean):-1:1)];

ax_(1) = subplot(1, 3, 1);
hold on
plot(x_axis_H1_stdv, (y_axis_H1_stdv), 'k', 'LineWidth', 0.25)
plot((H1_mean), 'k', 'LineWidth', 1)
xline(length(range_)/2)
set(gca, 'XTick', range_XTick)
title('joint')

ax_(2) = subplot(1, 3, 2);
hold on
plot(x_axis_H0_stdv, (y_axis_H0_stdv), 'k', 'LineWidth', 0.25)
plot((H0_mean), 'k', 'LineWidth', 1)
xline(length(range_)/2)
set(gca, 'XTick', range_XTick)
title('Pr1*Pr2')

ax_(3) = subplot(1, 3, 3);
hold on
plot(x_axis_ratio_stdv, (y_axis_ratio_stdv), 'k', 'LineWidth', 0.25)
plot((ratio_mean), 'k', 'LineWidth', 1)
xline(length(range_)/2)
set(gca, 'XTick', range_XTick)
title('ratio')

linkaxes([ax_(1) ax_(2)], 'y')
ESN_Beautify_Plot(hFig, [3 1.3], 8)

end

