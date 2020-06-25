%% INSTRUCTION
% set the Matlab pwd to inside 'data_125d_sorted' where the program can see
% 'compress_pCells', 'single_pCells', and 'double_pCells'
% Add data to the ALL_PCELL_COMPRESSED_DATA using examples in ESN_rebuild_ALL_PCELL_COMPRESSED_DATA
% after adding data, update the bundle_inds in ESN_change_overall_prob_tuning_for_bundles
% run ESN_change_overall_prob_tuning_for_bundles, and save the updates ALL_PCELL_COMPRESSED_DATA


%% function temp_ESN_plot_data_compress
function temp_ESN_plot_data_compress(file_name, file_path, description_pCell)
%% Load plot_data
if nargin < 1
    [file_name,file_path] = uigetfile([pwd filesep '*_plot_data*.mat'], 'Select plot_data file');
end
fprintf(['Loading ', file_name, ' ... ']);
plot_data_raw = load([file_path filesep file_name]);
plot_data_raw.file_name = file_name;
plot_data_raw.file_path = file_path;
fprintf(' --> Completed. \n')

%% pCell description
if nargin < 1
    description_pCell = input('Provide pCell description in the following format "single/double-pCell-00, num-01, pair-01" :','s');
    if isempty(description_pCell)
        description_pCell = 'single-pCell-00, num-00, pair-00';
    end
end

%% Build plot_data_compress
clearvars plot_data_compress;
fprintf(['Build plot_data_compress ', ' ... ']);
plot_data_compress.id          = plot_data_raw.file_name(1:16);
plot_data_compress.description = description_pCell;
plot_data_compress.file_name   = plot_data_raw.file_name;
plot_data_compress.file_path   = plot_data_raw.file_path;
field_names_plot_data_raw = fieldnames(plot_data_raw);
field_names_plot_data_raw = field_names_plot_data_raw(contains(field_names_plot_data_raw, 'raster_data_'));
for counter_field_names_plot_data_raw = 1 : 1 : length(field_names_plot_data_raw)
    field_name_plot_data_raw = field_names_plot_data_raw{counter_field_names_plot_data_raw};
    field_names_raster_data = fieldnames(plot_data_raw.(field_name_plot_data_raw));
    for counter_field_names_raster_data = 1 : 1 : length(field_names_raster_data)
        field_name_raster_data = field_names_raster_data{counter_field_names_raster_data};
        data_field_name_raster_data = plot_data_raw.(field_name_plot_data_raw).(field_name_raster_data);
        plot_data_compress.(field_name_plot_data_raw).(field_name_raster_data) = nanmean(data_field_name_raster_data);
        if contains(field_name_raster_data, '_CS_')
            plot_data_compress.(field_name_plot_data_raw).([field_name_raster_data '_sparse']) = sparse(data_field_name_raster_data);
            plot_data_compress.(field_name_plot_data_raw).([field_name_raster_data '_numTrial']) = size(data_field_name_raster_data, 1);
        end
    end
end

plot_data_compress.raster_data_cue_present.inds_span    = plot_data_raw.plot_data_cue_present.inds_span;
plot_data_compress.raster_data_primSac_onset.inds_span  = plot_data_raw.plot_data_primSac_onset.inds_span;
plot_data_compress.raster_data_primSac_offset.inds_span = plot_data_raw.plot_data_primSac_offset.inds_span;
plot_data_compress.raster_data_corrSac_onset.inds_span  = plot_data_raw.plot_data_corrSac_onset.inds_span;

plot_data_compress.Neural_Properties_data.SS_num = ...
    length(plot_data_raw.Neural_Properties_data.CH_sorted.SS_data.SS_time);
plot_data_compress.Neural_Properties_data.SS_duration = ...
    (plot_data_raw.Neural_Properties_data.CH_sorted.SS_data.SS_time(end)) - (plot_data_raw.Neural_Properties_data.CH_sorted.SS_data.SS_time(1));
plot_data_compress.Neural_Properties_data.SS_firing_rate = ...
    plot_data_compress.Neural_Properties_data.SS_num / plot_data_compress.Neural_Properties_data.SS_duration;
plot_data_compress.Neural_Properties_data.SS_time = ...
    plot_data_raw.Neural_Properties_data.CH_sorted.SS_data.SS_time;
plot_data_compress.Neural_Properties_data.SS_waveform = ...
    nanmean(plot_data_raw.Neural_Properties_data.CH_sorted.SS_data.SS_waveform);

plot_data_compress.Neural_Properties_data.CS_num = ...
    length(plot_data_raw.Neural_Properties_data.CH_sorted.CS_data.CS_time);
plot_data_compress.Neural_Properties_data.CS_firing_rate = ...
    plot_data_compress.Neural_Properties_data.CS_num / plot_data_compress.Neural_Properties_data.SS_duration;
plot_data_compress.Neural_Properties_data.CS_time = ...
    plot_data_raw.Neural_Properties_data.CH_sorted.CS_data.CS_time;
plot_data_compress.Neural_Properties_data.CS_waveform = ...
    nanmean(plot_data_raw.Neural_Properties_data.CH_sorted.CS_data.CS_waveform);

field_names_Corr_data = fieldnames(plot_data_raw.Neural_Properties_data.CH_sorted.Corr_data);
for counter_field_names_Corr_data = 1 : 1 : length(field_names_Corr_data)
    field_name_Corr_data = field_names_Corr_data{counter_field_names_Corr_data};
    data_field_name_Corr_data = plot_data_raw.Neural_Properties_data.CH_sorted.Corr_data.(field_name_Corr_data);
    plot_data_compress.Neural_Properties_data.Corr_data.(field_name_Corr_data) = nanmean(data_field_name_Corr_data);
end
fprintf(' --> Completed. \n')

%% CS_Tuning
fprintf(['Build CS_Tuning ', ' ... ']);
clearvars CS_Tuning;
range_inds_probability = 101:300;

CS_Tuning.cue_present_right = full(plot_data_compress.raster_data_cue_present.train_data_logic_CS_right_sparse);
CS_Tuning.cue_present_top   = full(plot_data_compress.raster_data_cue_present.train_data_logic_CS_top_sparse);
CS_Tuning.cue_present_left  = full(plot_data_compress.raster_data_cue_present.train_data_logic_CS_left_sparse);
CS_Tuning.cue_present_down  = full(plot_data_compress.raster_data_cue_present.train_data_logic_CS_down_sparse);

prob_right = sum( sum(CS_Tuning.cue_present_right(:,range_inds_probability) ,2) > 0 ) / size(CS_Tuning.cue_present_right, 1);
prob_top   = sum( sum(CS_Tuning.cue_present_top(  :,range_inds_probability) ,2) > 0 ) / size(CS_Tuning.cue_present_top, 1);
prob_left  = sum( sum(CS_Tuning.cue_present_left( :,range_inds_probability) ,2) > 0 ) / size(CS_Tuning.cue_present_left, 1);
prob_down  = sum( sum(CS_Tuning.cue_present_down( :,range_inds_probability) ,2) > 0 ) / size(CS_Tuning.cue_present_down, 1);
CS_Tuning.cue_present_prob = [prob_right prob_top prob_left prob_down];

CS_Tuning.primSac_onset_right = full(plot_data_compress.raster_data_primSac_onset.train_data_logic_CS_right_sparse);
CS_Tuning.primSac_onset_top   = full(plot_data_compress.raster_data_primSac_onset.train_data_logic_CS_top_sparse);
CS_Tuning.primSac_onset_left  = full(plot_data_compress.raster_data_primSac_onset.train_data_logic_CS_left_sparse);
CS_Tuning.primSac_onset_down  = full(plot_data_compress.raster_data_primSac_onset.train_data_logic_CS_down_sparse);

prob_right = sum( sum(CS_Tuning.primSac_onset_right(:,range_inds_probability) ,2) > 0 ) / size(CS_Tuning.primSac_onset_right, 1);
prob_top   = sum( sum(CS_Tuning.primSac_onset_top(  :,range_inds_probability) ,2) > 0 ) / size(CS_Tuning.primSac_onset_top, 1);
prob_left  = sum( sum(CS_Tuning.primSac_onset_left( :,range_inds_probability) ,2) > 0 ) / size(CS_Tuning.primSac_onset_left, 1);
prob_down  = sum( sum(CS_Tuning.primSac_onset_down( :,range_inds_probability) ,2) > 0 ) / size(CS_Tuning.primSac_onset_down, 1);
CS_Tuning.primSac_onset_prob = [prob_right prob_top prob_left prob_down];

CS_Tuning.primSac_offset_right = full(plot_data_compress.raster_data_primSac_offset.train_data_logic_CS_right_sparse);
CS_Tuning.primSac_offset_top   = full(plot_data_compress.raster_data_primSac_offset.train_data_logic_CS_top_sparse);
CS_Tuning.primSac_offset_left  = full(plot_data_compress.raster_data_primSac_offset.train_data_logic_CS_left_sparse);
CS_Tuning.primSac_offset_down  = full(plot_data_compress.raster_data_primSac_offset.train_data_logic_CS_down_sparse);

prob_right = sum( sum(CS_Tuning.primSac_offset_right(:,range_inds_probability) ,2) > 0 ) / size(CS_Tuning.primSac_offset_right, 1);
prob_top   = sum( sum(CS_Tuning.primSac_offset_top(  :,range_inds_probability) ,2) > 0 ) / size(CS_Tuning.primSac_offset_top, 1);
prob_left  = sum( sum(CS_Tuning.primSac_offset_left( :,range_inds_probability) ,2) > 0 ) / size(CS_Tuning.primSac_offset_left, 1);
prob_down  = sum( sum(CS_Tuning.primSac_offset_down( :,range_inds_probability) ,2) > 0 ) / size(CS_Tuning.primSac_offset_down, 1);
CS_Tuning.primSac_offset_prob = [prob_right prob_top prob_left prob_down];

CS_Tuning.corrSac_onset_right = full(plot_data_compress.raster_data_corrSac_onset.train_data_logic_CS_right_sparse);
CS_Tuning.corrSac_onset_top   = full(plot_data_compress.raster_data_corrSac_onset.train_data_logic_CS_top_sparse);
CS_Tuning.corrSac_onset_left  = full(plot_data_compress.raster_data_corrSac_onset.train_data_logic_CS_left_sparse);
CS_Tuning.corrSac_onset_down  = full(plot_data_compress.raster_data_corrSac_onset.train_data_logic_CS_down_sparse);

prob_right = sum( sum(CS_Tuning.corrSac_onset_right(:,range_inds_probability) ,2) > 0 ) / size(CS_Tuning.corrSac_onset_right, 1);
prob_top   = sum( sum(CS_Tuning.corrSac_onset_top(  :,range_inds_probability) ,2) > 0 ) / size(CS_Tuning.corrSac_onset_top, 1);
prob_left  = sum( sum(CS_Tuning.corrSac_onset_left( :,range_inds_probability) ,2) > 0 ) / size(CS_Tuning.corrSac_onset_left, 1);
prob_down  = sum( sum(CS_Tuning.corrSac_onset_down( :,range_inds_probability) ,2) > 0 ) / size(CS_Tuning.corrSac_onset_down, 1);
CS_Tuning.corrSac_onset_prob = [prob_right prob_top prob_left prob_down];

CS_Tuning.overall_prob = nanmean([CS_Tuning.cue_present_prob; CS_Tuning.primSac_onset_prob; CS_Tuning.primSac_offset_prob; CS_Tuning.corrSac_onset_prob]);
overall_prob = CS_Tuning.overall_prob;
[~, CS_000] = max(CS_Tuning.overall_prob);
if     CS_000 == 1
    if overall_prob(2) >= overall_prob(4)
        CS_090 = 2; CS_270 = 4; CS_180 = 3;
    else
        CS_090 = 4; CS_270 = 2; CS_180 = 3;
    end
elseif CS_000 == 2
    if overall_prob(1) >= overall_prob(3)
        CS_090 = 1; CS_270 = 3; CS_180 = 4;
    else
        CS_090 = 3; CS_270 = 1; CS_180 = 4;
    end
elseif CS_000 == 3
    if overall_prob(2) >= overall_prob(4)
        CS_090 = 2; CS_270 = 4; CS_180 = 1;
    else
        CS_090 = 4; CS_270 = 2; CS_180 = 1;
    end
elseif CS_000 == 4
    if overall_prob(1) >= overall_prob(3)
        CS_090 = 1; CS_270 = 3; CS_180 = 2;
    else
        CS_090 = 3; CS_270 = 1; CS_180 = 2;
    end
else
    error('CS_000 is not defined');
end

CS_Tuning.overall_prob_tuning = [CS_000, CS_090, CS_180, CS_270];

field_names_CS_Tuning = fieldnames(CS_Tuning);
field_names_CS_Tuning = field_names_CS_Tuning(contains(field_names_CS_Tuning, '_prob'));
for counter_field_names_CS_Tuning = 1 : 1 : length(field_names_CS_Tuning)
    field_name_CS_Tuning = field_names_CS_Tuning{counter_field_names_CS_Tuning};
    plot_data_compress.CS_Tuning.(field_name_CS_Tuning) = CS_Tuning.(field_name_CS_Tuning);
end
plot_data_compress.CS_Tuning.numTrials = ...
    plot_data_compress.raster_data_cue_present.train_data_logic_CS_right_numTrial + ...
    plot_data_compress.raster_data_cue_present.train_data_logic_CS_top_numTrial + ...
    plot_data_compress.raster_data_cue_present.train_data_logic_CS_left_numTrial + ...
    plot_data_compress.raster_data_cue_present.train_data_logic_CS_down_numTrial ;
fprintf(' --> Completed. \n')

%% Save plot_data_raw & plot_data_compress into original plot_data_raw file
fprintf(['Saving ', plot_data_raw.file_name, ' ... ']);
Neural_Properties_data     = plot_data_raw.Neural_Properties_data;
plot_data_cue_present      = plot_data_raw.plot_data_cue_present;
plot_data_primSac_onset    = plot_data_raw.plot_data_primSac_onset;
plot_data_primSac_offset   = plot_data_raw.plot_data_primSac_offset;
plot_data_corrSac_onset    = plot_data_raw.plot_data_corrSac_onset;
raster_data_cue_present    = plot_data_raw.raster_data_cue_present;
raster_data_primSac_onset  = plot_data_raw.raster_data_primSac_onset;
raster_data_primSac_offset = plot_data_raw.raster_data_primSac_offset;
raster_data_corrSac_onset  = plot_data_raw.raster_data_corrSac_onset;

save([plot_data_raw.file_path, filesep, plot_data_raw.file_name], 'Neural_Properties_data',...
      'plot_data_cue_present',  'plot_data_primSac_onset',  'plot_data_primSac_offset',  'plot_data_corrSac_onset',...
      'raster_data_cue_present','raster_data_primSac_onset','raster_data_primSac_offset','raster_data_corrSac_onset',...
      'plot_data_compress',...
      '-v7.3')
fprintf(' --> Completed. \n')

%% Save plot_data_compress in ALL_PCELL_COMPRESSED_DATA
clearvars -except plot_data_compress;
ALL_PCELL_COMPRESSED_DATA_file_name = 'ALL_PCELL_COMPRESSED_DATA.mat';
ALL_PCELL_COMPRESSED_DATA_file_path = './compress_pCells';
fprintf(['Saving ', ALL_PCELL_COMPRESSED_DATA_file_name, ' ... ']);
if isfile([ALL_PCELL_COMPRESSED_DATA_file_path filesep ALL_PCELL_COMPRESSED_DATA_file_name])
    load([ALL_PCELL_COMPRESSED_DATA_file_path filesep ALL_PCELL_COMPRESSED_DATA_file_name]);
    row_num_total = length(ALL_PCELL_COMPRESSED_DATA);
    row_num = row_num_total + 1;
    ALL_PCELL_COMPRESSED_DATA(row_num) = plot_data_compress;
else
    ALL_PCELL_COMPRESSED_DATA(1) = plot_data_compress;
end
save([ALL_PCELL_COMPRESSED_DATA_file_path, filesep, ALL_PCELL_COMPRESSED_DATA_file_name],...
    'ALL_PCELL_COMPRESSED_DATA',...
    '-v7.3')
fprintf(' --> Completed. \n')

end

%% function LOOP OVER ALL_PCELL_COMPRESSED_DATA_bkp and rebiuld ALL_PCELL_COMPRESSED_DATA
function ESN_rebuild_ALL_PCELL_COMPRESSED_DATA
%% rebuild_ALL_PCELL_COMPRESSED_DATA

% first rename the 'ALL_PCELL_COMPRESSED_DATA.mat' to 'ALL_PCELL_COMPRESSED_DATA_bkp.mat'
% use ALL_PCELL_COMPRESSED_DATA_bkp and iterate from 1 to end and rebuild
% the ALL_PCELL_COMPRESSED_DATA
for counter_row = 1 : 66
    file_name         = ALL_PCELL_COMPRESSED_DATA(counter_row).file_name;
    file_path         = ALL_PCELL_COMPRESSED_DATA(counter_row).file_path;
    description_pCell = ALL_PCELL_COMPRESSED_DATA(counter_row).description;
    temp_ESN_plot_data_compress(file_name, file_path, description_pCell);
end


%% Manually add one file to the ALL_PCELL_COMPRESSED_DATA
file_name = '190423_142023_01_sorted_ESN_plot_data';
file_path = './single_pCells/single_pCell_22';
description_pCell = 'single-pCell-22, num-03, pair-01';
% file_name = '190516_150909_07_sorted_RSh_plot_data';
% file_path = './double_pCells/double_pCell_12';
% description_pCell = 'double-pCell-12, num-02, pair-02';
disp('##################################################################')
disp(file_name)
disp(description_pCell)
temp_ESN_plot_data_compress(file_name, file_path, description_pCell);

end

%% function ESN_change_overall_prob_tuning_for_bundles
function ESN_change_overall_prob_tuning_for_bundles
% change_overall_prob_tuning_for_bundles
%% Load ALL_PCELL_COMPRESSED_DATA
ALL_PCELL_COMPRESSED_DATA_file_name = 'ALL_PCELL_COMPRESSED_DATA.mat';
ALL_PCELL_COMPRESSED_DATA_file_path = './compress_pCells';
load([ALL_PCELL_COMPRESSED_DATA_file_path filesep ALL_PCELL_COMPRESSED_DATA_file_name]);

%% bundle_inds
bundle_inds = { ...
    01, ... % 190326_174047_06_sorted_plot_data, d00
    02, ... % 190326_174047_07_sorted_plot_data, d00
    03, ... % 190326_182755_03_sorted_plot_data, d01
    04, ... % 190326_182755_04_sorted_plot_data, d01
    [05, 07, 09], ... % 190327_135803_05_sorted_plot_data, d02
    [06, 08, 10], ... % 190327_135803_06_sorted_plot_data, d02
    11, ... % 190329_133538_03_sorted_plot_data, d03
    12, ... % 190329_133538_05_sorted_plot_data, d03
    13, ... % 190329_150427_03_sorted_plot_data, d04
    14, ... % 190329_150427_07_sorted_plot_data, d04
    15, ... % 190409_150223_03_sorted_ESN_plot_data, d05
    16, ... % 190409_150223_07_sorted_ESN_plot_data, d05
    17, ... % 190411_130256_03_sorted_plot_data, d06
    18, ... % 190411_130256_07_sorted_plot_data, d06
    [19, 21, 23, 25, 27], ... % 190412_115230_01_sorted_plot_data, d07
    [20, 22, 24, 26, 28], ... % 190412_115230_05_sorted_plot_data, d07
    [29, 31], ... % 190426_120355_05_sorted_plot_data, d08
    [30, 32], ... % 190426_120355_06_sorted_plot_data, d08
    [33, 35, 37], ... % 190426_125234_01_sorted_plot_data, d09
    [34, 36, 38], ... % 190426_125234_07_sorted_plot_data, d09
    [39, 41], ... % 190429_135343_05_sorted_plot_data, d10
    [40, 42], ... % 190429_135343_06_sorted_plot_data, d10
    [43, 45, 47], ... % 190429_142903_03_sorted_ESN_plot_data, d11
    [44, 46, 48], ... % 190429_142903_07_sorted_ESN_plot_data, d11
    [49, 51], ... % 190516_144904_02_sorted_RSh_plot_data, d12
    [50, 52], ... % 190516_144904_07_sorted_RSh_plot_data, d12
    [53, 54], ... % 190327_152414_03_sorted_plot_data, s01
    55, ... % 190402_150756_05_sorted_plot_data, s02
    [56, 57], ... % 190403_141313_04_sorted_plot_data, s03
    [58, 59, 60], ... % 190403_151506_03_sorted_plot_data, s04
    61, ... % 190404_110707_07_sorted_plot_data, s05
    62, ... % 190409_143035_04_sorted_plot_data, s06
    [63, 64, 65], ... % 190410_140235_04_sorted_plot_data, s07
    [66, 67, 68, 69], ... % 190415_135953_04_sorted_plot_data, s08
    [70, 71], ... % 190415_153818_07_sorted_plot_data, s09
    [72, 73], ... % 190416_142432_02_sorted_plot_data, s10
    [74, 75, 76], ... % 190424_151111_04_sorted_plot_data, s11
    77, ... % 190411_124025_04_sorted_plot_data, s12
    [78, 79, 80, 81], ... % 190430_132054_05_sorted_ESN_plot_data, s13
    82, ... % 190501_141616_03_sorted_RSh_plot_data, s14
    83, ... % 190515_135233_04_sorted_ESN_plot_data, s15
    [84, 85], ... % 190517_121917_05_sorted_RSh_plot_data, s16
    [86, 87, 88], ... % 190408_130716_04_sorted_ESN_plot_data, s17
    [89 90], ... % 190409_151632_03_sorted_ESN_plot_data, s18
    [91, 92, 93], ... % 190409_152524_07_sorted_ESN_plot_data, s19
    94, ... % 190411_133808_03_sorted_ESN_plot_data, s20
    [95, 96], ... % 190422_143306_04_sorted_plot_data, s21
    [97, 98, 99], ... % 190423_132324_01_sorted_ESN_plot_data, s22
    ...
    };

%% loop over bundle_inds
for counter_bundle_inds = 1 : length(bundle_inds)
index_nums = bundle_inds{counter_bundle_inds};
overall_prob_index_nums = [0, 0, 0, 0];
numTrials_index_nums = 0;
for counter_index = 1 : length(index_nums)
    index_num = index_nums(counter_index);
    overall_prob_index_num = ALL_PCELL_COMPRESSED_DATA(index_num).CS_Tuning.overall_prob;
    numTrials_index_num = ALL_PCELL_COMPRESSED_DATA(index_num).CS_Tuning.numTrials;
    overall_prob_index_nums = overall_prob_index_nums + (overall_prob_index_num * numTrials_index_num);
    numTrials_index_nums = numTrials_index_nums + numTrials_index_num;
end

overall_prob_bundle = overall_prob_index_nums ./ numTrials_index_nums;
numTrials_bundle = numTrials_index_nums;
[~, CS_000] = max(overall_prob_bundle);
if     CS_000 == 1
    if overall_prob_bundle(2) >= overall_prob_bundle(4)
        CS_090 = 2; CS_270 = 4; CS_180 = 3;
    else
        CS_090 = 4; CS_270 = 2; CS_180 = 3;
    end
elseif CS_000 == 2
    if overall_prob_bundle(1) >= overall_prob_bundle(3)
        CS_090 = 1; CS_270 = 3; CS_180 = 4;
    else
        CS_090 = 3; CS_270 = 1; CS_180 = 4;
    end
elseif CS_000 == 3
    if overall_prob_bundle(2) >= overall_prob_bundle(4)
        CS_090 = 2; CS_270 = 4; CS_180 = 1;
    else
        CS_090 = 4; CS_270 = 2; CS_180 = 1;
    end
elseif CS_000 == 4
    if overall_prob_bundle(1) >= overall_prob_bundle(3)
        CS_090 = 1; CS_270 = 3; CS_180 = 2;
    else
        CS_090 = 3; CS_270 = 1; CS_180 = 2;
    end
else
    error('CS_000 is not defined');
end
overall_prob_tuning_bundle = [CS_000, CS_090, CS_180, CS_270];
for counter_index = 1 : length(index_nums)
    index_num = index_nums(counter_index);
    ALL_PCELL_COMPRESSED_DATA(index_num).CS_Tuning.overall_prob_bundle        = overall_prob_bundle;
    ALL_PCELL_COMPRESSED_DATA(index_num).CS_Tuning.overall_prob_tuning_bundle = overall_prob_tuning_bundle;
    ALL_PCELL_COMPRESSED_DATA(index_num).CS_Tuning.numTrials_bundle           = numTrials_bundle;
end

end % for counter_bundle_inds = 1 : length(bundle_inds)

%% save ALL_PCELL_COMPRESSED_DATA
fprintf(['Saving ', ALL_PCELL_COMPRESSED_DATA_file_name, ' ... ']);
save([ALL_PCELL_COMPRESSED_DATA_file_path, filesep, ALL_PCELL_COMPRESSED_DATA_file_name],...
    'ALL_PCELL_COMPRESSED_DATA', 'bundle_inds', ...
    '-v7.3')
fprintf(' --> Completed. \n')

end

