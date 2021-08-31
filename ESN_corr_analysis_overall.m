function ESN_corr_analysis_overall()
%% Steps to prepare the data
clc; clear; close all;
tic
% (1) build_pair_data_corr(); % combine the .psort files and form one corr data per pair. save the data in the ALL_PCELL_### folder.
% (2) build_population_corr_data(); % load pair_data_corr files (_combine_) and form population_corr_data.
toc

%% Plot functions
% (1) plot_corr_data(); % load corr_data_population_ and plot interactions

end

%% function build_pair_data_corr()
function build_pair_data_corr()
clc; close all;
flag_pair_list = true; % This should be true. DO NOT change it to false
pCell_list = ESN_build_pCell_list(flag_pair_list);
path_data_monkey_sorted = uigetdir;
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);
ALL_PCELL_name = ['ALL_PCELL_' num2str(num_pCells)];
if ~exist([path_data_monkey_sorted ALL_PCELL_name filesep 'pair_data_corr'], 'dir')
    mkdir([path_data_monkey_sorted ALL_PCELL_name filesep 'pair_data_corr']);
end

ESN_global_variables();
global corr_data_list

pair_names = build_pair_names();
%% Loop over pCells
for counter_pCell = 1 : 2 : num_pCells-1
    fprintf(['### ' 'Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ###' '\n']);
    num_recording = nansum(pCell_list_isstr(counter_pCell, :));
    %% Init EPHYS_session
    EPHYS_session.CH1.SS_index = [];
    EPHYS_session.CH1.SS_time  = [];
    EPHYS_session.CH1.CS_index = [];
    EPHYS_session.CH1.CS_time  = [];
    EPHYS_session.CH2.SS_index = [];
    EPHYS_session.CH2.SS_time  = [];
    EPHYS_session.CH2.CS_index = [];
    EPHYS_session.CH2.CS_time  = [];
    for counter_variable1 = 1 : length(corr_data_list)
        var_name_1 = corr_data_list{counter_variable1};
        for counter_variable2 = 1 : length(corr_data_list)
            var_name_2 = corr_data_list{counter_variable2};
            EPHYS_session.Corr_data.([var_name_2 'x' var_name_1 '_inds_span'])     = zeros(0,100);
            EPHYS_session.Corr_data.([var_name_2 'x' var_name_1 '_bin_size_time']) = [];
            EPHYS_session.Corr_data.([var_name_2 'x' var_name_1 '_prob'])          = false(0,100);
            EPHYS_session.Corr_data.([var_name_2 'x' var_name_1 '_num_spikes'])    = [];
        end
    end
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
        %% RE-RUN add_ephys_sac_sorter
        file_path = [path_data_monkey_sorted subFolder_month subFolder_day subFolder_recording subFolder_figs];
        if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
        file_name_1 = pCell_list{counter_pCell,   counter_recording};
        file_name_2 = pCell_list{counter_pCell+1, counter_recording};
        DATA_PSORT_CH1 = Psort_read_psort([file_path file_name_1 '.psort']);
        DATA_PSORT_CH2 = Psort_read_psort([file_path file_name_2 '.psort']);
        EPHYS_recording = combine_CH2xCH1(DATA_PSORT_CH1, DATA_PSORT_CH2);
        %% Add EPHYS_recording to EPHYS_session
        EPHYS_session.CH1.SS_index = vertcat(EPHYS_session.CH1.SS_index, EPHYS_recording.CH1.SS_index);
        EPHYS_session.CH1.SS_time  = vertcat(EPHYS_session.CH1.SS_time,  EPHYS_recording.CH1.SS_time);
        EPHYS_session.CH1.CS_index = vertcat(EPHYS_session.CH1.CS_index, EPHYS_recording.CH1.CS_index);
        EPHYS_session.CH1.CS_time  = vertcat(EPHYS_session.CH1.CS_time,  EPHYS_recording.CH1.CS_time);
        EPHYS_session.CH2.SS_index = vertcat(EPHYS_session.CH2.SS_index, EPHYS_recording.CH2.SS_index);
        EPHYS_session.CH2.SS_time  = vertcat(EPHYS_session.CH2.SS_time,  EPHYS_recording.CH2.SS_time);
        EPHYS_session.CH2.CS_index = vertcat(EPHYS_session.CH2.CS_index, EPHYS_recording.CH2.CS_index);
        EPHYS_session.CH2.CS_time  = vertcat(EPHYS_session.CH2.CS_time,  EPHYS_recording.CH2.CS_time);
        for counter_variable1 = 1 : length(corr_data_list)
            var_name_1 = corr_data_list{counter_variable1};
            for counter_variable2 = 1 : length(corr_data_list)
                var_name_2 = corr_data_list{counter_variable2};
                EPHYS_session.Corr_data.([var_name_2 'x' var_name_1 '_inds_span'])     = ...
                    vertcat(EPHYS_session.Corr_data.([var_name_2 'x' var_name_1 '_inds_span']),  EPHYS_recording.Corr_data.([var_name_2 'x' var_name_1 '_inds_span']));
                EPHYS_session.Corr_data.([var_name_2 'x' var_name_1 '_bin_size_time']) = ...
                    vertcat(EPHYS_session.Corr_data.([var_name_2 'x' var_name_1 '_bin_size_time']),  EPHYS_recording.Corr_data.([var_name_2 'x' var_name_1 '_bin_size_time']));
                EPHYS_session.Corr_data.([var_name_2 'x' var_name_1 '_prob'])          = ...
                    vertcat(EPHYS_session.Corr_data.([var_name_2 'x' var_name_1 '_prob']),  EPHYS_recording.Corr_data.([var_name_2 'x' var_name_1 '_prob']));
                EPHYS_session.Corr_data.([var_name_2 'x' var_name_1 '_num_spikes'])    = ...
                    vertcat(EPHYS_session.Corr_data.([var_name_2 'x' var_name_1 '_num_spikes']),  EPHYS_recording.Corr_data.([var_name_2 'x' var_name_1 '_num_spikes']));
            end
        end
        
    end
    %% Save data
    counter_pair = floor(counter_pCell/2) + 1;
    pair_name = pair_names{counter_pair, 1};
    save([path_data_monkey_sorted ALL_PCELL_name filesep 'pair_data_corr' filesep pair_name '.mat'],...
        'EPHYS_session','-v7.3')
end
fprintf('### ALL DONE. ###\n')
end

%% function combine_CH2xCH1
function EPHYS = combine_CH2xCH1(DATA_PSORT_CH1, DATA_PSORT_CH2)

%% build EPHYS.CH1 from DATA_PSORT_CH1
ch_time = double(DATA_PSORT_CH1.topLevel_data.ch_time);
SS_index = find(logical(double(DATA_PSORT_CH1.topLevel_data.ss_index)));
SS_time = ch_time(SS_index);
EPHYS.CH1.SS_index = SS_index;
EPHYS.CH1.SS_time  = SS_time;

ch_time = double(DATA_PSORT_CH1.topLevel_data.ch_time);
CS_index = find(logical(double(DATA_PSORT_CH1.topLevel_data.cs_index)));
CS_time = ch_time(CS_index);
EPHYS.CH1.CS_index = CS_index;
EPHYS.CH1.CS_time  = CS_time;

%% build EPHYS.CH2 from DATA_PSORT_CH2
ch_time = double(DATA_PSORT_CH2.topLevel_data.ch_time);
SS_index = find(logical(double(DATA_PSORT_CH2.topLevel_data.ss_index)));
SS_time = ch_time(SS_index);
EPHYS.CH2.SS_index = SS_index;
EPHYS.CH2.SS_time  = SS_time;

ch_time = double(DATA_PSORT_CH2.topLevel_data.ch_time);
CS_index = find(logical(double(DATA_PSORT_CH2.topLevel_data.cs_index)));
CS_time = ch_time(CS_index);
EPHYS.CH2.CS_index = CS_index;
EPHYS.CH2.CS_time  = CS_time;

%% build EPHYS.Corr_data
fprintf(['Building Corr_data' ' ...'])
SS1 = EPHYS.CH1.SS_time;
SS2 = EPHYS.CH2.SS_time;
CS1 = EPHYS.CH1.CS_time;
CS2 = EPHYS.CH2.CS_time;

ESN_global_variables();
global corr_data_list

for counter_variable1 = 1 : length(corr_data_list)
    var_name_1 = corr_data_list{counter_variable1};
    for counter_variable2 = 1 : length(corr_data_list)
        var_name_2 = corr_data_list{counter_variable2};
        
        eval(['Corr_data = ESN_correlogram(' var_name_1 ', ' var_name_2 ');']);
        EPHYS.Corr_data.([var_name_2 'x' var_name_1 '_inds_span'])     = mean(Corr_data.CS_inds_span);
        EPHYS.Corr_data.([var_name_2 'x' var_name_1 '_bin_size_time']) = mean(Corr_data.CS_bin_size_time);
        EPHYS.Corr_data.([var_name_2 'x' var_name_1 '_prob'])          = Corr_data.CS_CSxSS_AUTO;
        EPHYS.Corr_data.([var_name_2 'x' var_name_1 '_num_spikes'])    = size(Corr_data.CS_CSxSS_AUTO, 1);
    end
end
fprintf(' --> Completed. \n')
end

%% function build_pair_names()
function pair_names = build_pair_names()
%% build_pair_names
flag_pair_list = true;
pCell_list = ESN_build_pCell_list(flag_pair_list);
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);
pair_names = cell(num_pCells/2, 1);
for counter_pCell = 1 : 2 : num_pCells-1
    file_name_1 = pCell_list{counter_pCell,   1};
    file_name_2 = pCell_list{counter_pCell+1, 1};
    pair_name_data = file_name_1(1:13);
    if file_name_1(18) == 's'
        pair_name_ch1  = file_name_1(15:16);
    elseif isstrprop(file_name_1(18),'digit')
        pair_name_ch1  = file_name_1(15:18);
    else
        error('cell id does not follow the standards.')
    end
    if file_name_2(18) == 's'
        pair_name_ch2  = file_name_2(15:16);
    elseif isstrprop(file_name_2(18),'digit')
        pair_name_ch2  = file_name_2(15:18);
    else
        error('cell id does not follow the standards.')
    end
    num_recording = nansum(pCell_list_isstr(counter_pCell, :));
    pair_name_combine = ['combine' '_' num2str(num_recording)];
    pair_name = [pair_name_data '_' pair_name_ch1 '_' pair_name_ch2 '_' pair_name_combine];
    counter_pair = floor(counter_pCell/2) + 1;
    pair_names{counter_pair, 1} = pair_name;
end

end


%% function build_population_corr_data
function build_population_corr_data()
%% clc;clear
clc; close all;
ESN_global_variables();
global corr_data_list
path_corr_data = uigetdir;
if ~strcmp(path_corr_data(end), filesep);path_corr_data = [path_corr_data filesep];end
pair_names = build_pair_names();
num_pairs = size(pair_names, 1);

%% Init variables
corr_data_population = struct;
for counter_variable1 = 1 : length(corr_data_list)
    var_name_1 = corr_data_list{counter_variable1};
    for counter_variable2 = 1 : length(corr_data_list)
        var_name_2 = corr_data_list{counter_variable2};
        corr_data_population.([var_name_2 'x' var_name_1 '_inds_span'])     = zeros(num_pairs,100);
        corr_data_population.([var_name_2 'x' var_name_1 '_bin_size_time']) = zeros(num_pairs,1);
        corr_data_population.([var_name_2 'x' var_name_1 '_prob'])          = cell(num_pairs,1);
        corr_data_population.([var_name_2 'x' var_name_1 '_num_spikes'])    = zeros(num_pairs,1);
        corr_data_population.([var_name_2 'x' var_name_1 '_prob_base'])     = zeros(num_pairs,1);
    end
end

%% Loop over pairs
for counter_pair = 1 : num_pairs
    fprintf(['### ' 'Analyzing pair no. ', num2str(counter_pair), ' / ' num2str(num_pairs) ' ###' '\n']);
    %% Load pair data
    pair_name = pair_names{counter_pair};
    load([path_corr_data pair_name '.mat'], 'EPHYS_session');
    %% Loop over variables
    for counter_variable1 = 1 : length(corr_data_list)
        var_name_1 = corr_data_list{counter_variable1};
        var_name_1_time = EPHYS_session.(['CH' var_name_1(3)]).([var_name_1(1:2) '_time']);
        for counter_variable2 = 1 : length(corr_data_list)
            var_name_2 = corr_data_list{counter_variable2};
            
            inds_span_ = EPHYS_session.Corr_data.([var_name_2 'x' var_name_1 '_inds_span']);
            bin_size_time_ = EPHYS_session.Corr_data.([var_name_2 'x' var_name_1 '_bin_size_time']);
            prob_ = EPHYS_session.Corr_data.([var_name_2 'x' var_name_1 '_prob']);
            num_spikes_ = EPHYS_session.Corr_data.([var_name_2 'x' var_name_1 '_num_spikes']);
            if size(inds_span_, 1) > 1
                inds_span_ = mean(inds_span_);
            end
            bin_size_time_ = mean(bin_size_time_);
            num_spikes_ = sum(num_spikes_);
            
            diff_var_name_1_time = diff(var_name_1_time);
            diff_var_name_1_time(abs(diff_var_name_1_time)> 5 ) = 0;
            duration = sum(diff_var_name_1_time);
            prob_base_ = size(var_name_1_time,1) * bin_size_time_ / duration ;
            
            corr_data_population.([var_name_2 'x' var_name_1 '_inds_span'])(counter_pair, :)     = inds_span_;
            corr_data_population.([var_name_2 'x' var_name_1 '_bin_size_time'])(counter_pair, 1) = bin_size_time_;
            corr_data_population.([var_name_2 'x' var_name_1 '_prob']){counter_pair, 1}          = prob_;
            corr_data_population.([var_name_2 'x' var_name_1 '_num_spikes'])(counter_pair, :)    = num_spikes_;
            corr_data_population.([var_name_2 'x' var_name_1 '_prob_base'])(counter_pair, :)     = prob_base_;
        end
    end
end

%% Save data
fprintf(['Saving .mat files' ' ...'])
save([path_corr_data '..' filesep 'corr_data_population' '.mat'], 'corr_data_population', '-v7.3');
fprintf(' --> Completed. \n')
end

%% function build_overall_sync
function build_overall_sync()
%% load corr_data_population
load(['corr_data_population' '.mat'], 'corr_data_population');
global expand_index
ESN_global_variables();
expand_index = 0;

%% overall_sync_index
fprintf(['overall_sync_index' ' ...'])
num_pairs = size(corr_data_population.SS1xSS1_num_spikes, 1);

overall_sync_index = zeros(num_pairs, 1);
for counter_pair = 1 : num_pairs
    raster_SS2xSS1 = corr_data_population.(['SS2' 'x' 'SS1' '_prob']){counter_pair,1}(:,50);
    raster_SS1xSS2 = corr_data_population.(['SS1' 'x' 'SS2' '_prob']){counter_pair,1}(:,50);
    SS1_prob_base = corr_data_population.(['SS1' 'x' 'SS1' '_prob_base'])(counter_pair,1);
    SS2_prob_base = corr_data_population.(['SS2' 'x' 'SS2' '_prob_base'])(counter_pair,1);
    SS1_num_spikes = corr_data_population.(['SS1' 'x' 'SS1' '_num_spikes'])(counter_pair,1);
    SS2_num_spikes = corr_data_population.(['SS2' 'x' 'SS2' '_num_spikes'])(counter_pair,1);
    SS_joint_num_spikes = max([sum(raster_SS2xSS1) sum(raster_SS1xSS2)]);
    SS1_duration = SS1_num_spikes ./ SS1_prob_base;
    SS2_duration = SS2_num_spikes ./ SS2_prob_base;
    SS_joint_duration = max([SS1_duration SS2_duration]);
    SS_joint_prob_base = SS_joint_num_spikes ./ SS_joint_duration;

    sync_index = SS_joint_prob_base ./ SS1_prob_base ./ SS2_prob_base;
    overall_sync_index(counter_pair, 1) = sync_index;
end

fprintf(' --> Completed. \n')
end

%% function plot_corr_data
function plot_corr_data()
%% clc;clear
clc; close all;
%% load corr_data_population
path_corr_data = uigetdir;
if ~strcmp(path_corr_data(end), filesep);path_corr_data = [path_corr_data filesep];end
load([path_corr_data '..' filesep 'corr_data_population' '.mat'], 'corr_data_population');
global expand_index
expand_index = 0;
% corr_data_population = exclude_overlap_CS(corr_data_population);
corr_data_population = average_prob(corr_data_population);

%% Within Cell SSxSS
SS1xSS1_prob = corr_data_population.SS1xSS1_prob(:,21:80);
SS1xSS1_prob_base = corr_data_population.SS1xSS1_prob_base;
SS1xSS1_norm = SS1xSS1_prob ./ repmat(SS1xSS1_prob_base, 1, size(SS1xSS1_prob, 2));
SS2xSS2_prob = corr_data_population.SS2xSS2_prob(:,21:80);
SS2xSS2_prob_base = corr_data_population.SS2xSS2_prob_base;
SS2xSS2_norm = SS2xSS2_prob ./ repmat(SS2xSS2_prob_base, 1, size(SS2xSS2_prob, 2));
SSxSS_norm = [SS1xSS1_norm; SS2xSS2_norm];
SSxSS_norm(:,30) = nan;

data_mean_x_axis = nanmean(corr_data_population.CS1xSS2_inds_span(:,21:80));
data_mean_ = nanmean(SSxSS_norm);
data_sem_ = nanstd(SSxSS_norm)./sqrt(size(SSxSS_norm,1));
data_sem_p_ = data_mean_ + data_sem_;
data_sem_m_ = data_mean_ - data_sem_;
data_sem_y_axis_ = [(data_sem_p_) flip(data_sem_m_)];
data_sem_x_axis_ = [(data_mean_x_axis) flip(data_mean_x_axis)];

hFig = figure(1);
clf(hFig);
hold on;
xline(0)
yline(1)
plot(data_sem_x_axis_ , data_sem_y_axis_,'-k', 'linewidth', 0.25)
plot(data_mean_x_axis , data_mean_,'-k', 'linewidth', 1)
xlabel('Time from SS1 (ms)')
ylabel('Normalized SS1 Firing')
set(gca, 'XTick', [-30 -20 -10 0 10 20 30])
xlim([-30 30])
ESN_Beautify_Plot(hFig, [2 1], 8)
% ESN_Beautify_Plot(hFig, [8 4], 12)

%% Within Cell CSxSS
CS1xSS1_prob = corr_data_population.CS1xSS1_prob;
CS1xSS1_prob_base = corr_data_population.CS1xSS1_prob_base;
CS1xSS1_norm = CS1xSS1_prob ./ repmat(CS1xSS1_prob_base, 1, size(CS1xSS1_prob, 2));
CS2xSS2_prob = corr_data_population.CS2xSS2_prob;
CS2xSS2_prob_base = corr_data_population.CS2xSS2_prob_base;
CS2xSS2_norm = CS2xSS2_prob ./ repmat(CS2xSS2_prob_base, 1, size(CS2xSS2_prob, 2));
CSxSS_norm = [CS1xSS1_norm; CS2xSS2_norm];

data_mean_x_axis = nanmean(corr_data_population.CS1xSS2_inds_span);
data_mean_ = nanmean(CSxSS_norm);
data_sem_ = nanstd(CSxSS_norm)./sqrt(size(CSxSS_norm,1));
data_sem_p_ = data_mean_ + data_sem_;
data_sem_m_ = data_mean_ - data_sem_;
data_sem_y_axis_ = [(data_sem_p_) flip(data_sem_m_)];
data_sem_x_axis_ = [(data_mean_x_axis) flip(data_mean_x_axis)];

hFig = figure(1);
clf(hFig);
hold on;
xline(0)
yline(1)
plot(data_sem_x_axis_ , data_sem_y_axis_,'-k', 'linewidth', 0.25)
plot(data_mean_x_axis , data_mean_,'-k', 'linewidth', 1)
xlabel('Time from CS1 (ms)')
ylabel('Normalized SS1 Firing')
set(gca, 'XTick', [-50 -25 0 25 50])
ESN_Beautify_Plot(hFig, [2 1], 8)
% ESN_Beautify_Plot(hFig, [8 4], 12)

%% Between Cell CSxSS
% CS1xSS2_prob = corr_data_population.CS1xSS2_prob(:,35:65);
% CS1xSS2_prob_base = corr_data_population.CS1xSS2_prob_base;
% CS1xSS2_norm = CS1xSS2_prob ./ repmat(CS1xSS2_prob_base, 1, size(CS1xSS2_prob, 2));
CS1xSS2_norm = corr_data_population.CS1xSS2_prob_norm(:,35:65);
% CS2xSS1_prob = corr_data_population.CS2xSS1_prob(:,35:65);
% CS2xSS1_prob_base = corr_data_population.CS2xSS1_prob_base;
% CS2xSS1_norm = CS2xSS1_prob ./ repmat(CS2xSS1_prob_base, 1, size(CS2xSS1_prob, 2));
CS2xSS1_norm = corr_data_population.CS2xSS1_prob_norm(:,35:65);
CSxSS_norm = (CS1xSS2_norm + CS2xSS1_norm) ./ 2;
% CSxSS_norm = (CS1xSS2_prob + CS2xSS1_prob) ./ 2;

data_mean_x_axis = nanmean(corr_data_population.CS1xSS2_inds_span(:,35:65));
data_mean_ = nanmean(CSxSS_norm);
data_sem_ = nanstd(CSxSS_norm)./sqrt(size(CSxSS_norm,1));
data_sem_p_ = data_mean_ + data_sem_;
data_sem_m_ = data_mean_ - data_sem_;
data_sem_y_axis_ = [(data_sem_p_) flip(data_sem_m_)];
data_sem_x_axis_ = [(data_mean_x_axis) flip(data_mean_x_axis)];

hFig = figure(2);
clf(hFig);
hold on;
xline(0)
yline(1)
plot(data_sem_x_axis_ , data_sem_y_axis_,'-k', 'linewidth', 0.25)
plot(data_mean_x_axis , data_mean_,'-k', 'linewidth', 1)
xlabel('Time from CS1 (ms)')
ylabel('Normalized SS2 Firing')
xlim([-10 10])
set(gca, 'XTick', -10:5:10)
ylim([0.70 1.10])
set(gca, 'YTick', 0.7:0.1:1.1)
ESN_Beautify_Plot(hFig, [1.5 1], 8)
% ESN_Beautify_Plot(hFig, [8 4], 12)


%% Between Cell SSxSS
SS1xSS2_prob = corr_data_population.SS1xSS2_prob(:,35:65);
SS1xSS2_prob_base = corr_data_population.SS1xSS2_prob_base;
SS1xSS2_norm = SS1xSS2_prob ./ repmat(SS1xSS2_prob_base, 1, size(SS1xSS2_prob, 2));
SS2xSS1_prob = corr_data_population.SS2xSS1_prob(:,35:65);
SS2xSS1_prob_base = corr_data_population.SS2xSS1_prob_base;
SS2xSS1_norm = SS2xSS1_prob ./ repmat(SS2xSS1_prob_base, 1, size(SS2xSS1_prob, 2));
SSxSS_norm = (SS1xSS2_norm + SS2xSS1_norm) ./ 2;

data_mean_x_axis = nanmean(corr_data_population.CS1xSS2_inds_span(:,35:65));
data_mean_ = nanmean(SSxSS_norm);
data_sem_ = nanstd(SSxSS_norm)./sqrt(size(SSxSS_norm,1));
data_sem_p_ = data_mean_ + data_sem_;
data_sem_m_ = data_mean_ - data_sem_;
data_sem_y_axis_ = [(data_sem_p_) flip(data_sem_m_)];
data_sem_x_axis_ = [(data_mean_x_axis) flip(data_mean_x_axis)];

hFig = figure(3);
clf(hFig);
hold on;
xline(0)
% yline(1)
plot(data_sem_x_axis_ , data_sem_y_axis_,'-k', 'linewidth', 0.25)
plot(data_mean_x_axis , data_mean_,'-k', 'linewidth', 1)
xlabel('Time from SS1 (ms)')
ylabel('Normalized SS2 Firing')
xlim([-10 10])
set(gca, 'XTick', -10:5:10)
ylim([1.00 1.40])
set(gca, 'YTick', 1.0:0.1:1.4)
ESN_Beautify_Plot(hFig, [1.5 1], 8)
% ESN_Beautify_Plot(hFig, [8 4], 12)

%% Between Cell CSxCS
CS1xCS2_prob = corr_data_population.CS1xCS2_prob;
CS1xCS2_prob_base = corr_data_population.CS1xCS2_prob_base;
CS1xCS2_norm = CS1xCS2_prob ./ repmat(CS1xCS2_prob_base, 1, size(CS1xCS2_prob, 2));
CS2xCS1_prob = corr_data_population.CS2xCS1_prob;
CS2xCS1_prob_base = corr_data_population.CS2xCS1_prob_base;
CS2xCS1_norm = CS2xCS1_prob ./ repmat(CS2xCS1_prob_base, 1, size(CS2xCS1_prob, 2));
CSxCS_norm = (CS1xCS2_norm + CS2xCS1_norm) ./ 2;
CSxCS_norm(:,50) = nan;

data_mean_x_axis = nanmean(corr_data_population.CS1xSS2_inds_span);
data_mean_ = nanmean(CSxCS_norm);
data_sem_ = nanstd(CSxCS_norm)./sqrt(size(CSxCS_norm,1));
data_sem_p_ = data_mean_ + data_sem_;
data_sem_m_ = data_mean_ - data_sem_;
data_sem_y_axis_ = [(data_sem_p_) flip(data_sem_m_)];
data_sem_x_axis_ = [(data_mean_x_axis) flip(data_mean_x_axis)];

hFig = figure(4);
clf(hFig);
hold on;
xline(0)
% yline(1)
plot(data_sem_x_axis_ , data_sem_y_axis_,'-k', 'linewidth', 0.25)
plot(data_mean_x_axis , data_mean_,'-k', 'linewidth', 1)
% plot(data_mean_x_axis , CSxCS_norm, 'linewidth', 0.5)
xlabel('Time from CS1 (ms)')
ylabel('Normalized CS2 Firing')
xlim([-40 40])
set(gca, 'XTick', -40:10:40)
ylim([1 3])
set(gca, 'YTick', 0.5:0.5:3)
ESN_Beautify_Plot(hFig, [1.5 1], 8)
% ESN_Beautify_Plot(hFig, [8 4], 12)

%% SSxSS vs CSxCS
CS1xCS2_prob = corr_data_population.CS1xCS2_prob;
CS1xCS2_prob_base = corr_data_population.CS1xCS2_prob_base;
CS1xCS2_norm = CS1xCS2_prob ./ repmat(CS1xCS2_prob_base, 1, size(CS1xCS2_prob, 2));
CS2xCS1_prob = corr_data_population.CS2xCS1_prob;
CS2xCS1_prob_base = corr_data_population.CS2xCS1_prob_base;
CS2xCS1_norm = CS2xCS1_prob ./ repmat(CS2xCS1_prob_base, 1, size(CS2xCS1_prob, 2));
CSxCS_norm = (CS1xCS2_norm + CS2xCS1_norm) ./ 2;

SS1xSS2_prob = corr_data_population.SS1xSS2_prob;
SS1xSS2_prob_base = corr_data_population.SS1xSS2_prob_base;
SS1xSS2_norm = SS1xSS2_prob ./ repmat(SS1xSS2_prob_base, 1, size(SS1xSS2_prob, 2));
SS2xSS1_prob = corr_data_population.SS2xSS1_prob;
SS2xSS1_prob_base = corr_data_population.SS2xSS1_prob_base;
SS2xSS1_norm = SS2xSS1_prob ./ repmat(SS2xSS1_prob_base, 1, size(SS2xSS1_prob, 2));
SSxSS_norm = (SS1xSS2_norm + SS2xSS1_norm) ./ 2;

data_x_axis = nanmean((CSxCS_norm(:,45:49) + CSxCS_norm(:,51:55)) ./ 2, 2);
data_y_axis = nanmean((SSxSS_norm(:,50) + SSxSS_norm(:,50)) ./ 2, 2);

% win_len = '15ms';
% load(['corr_data_population_' win_len '.mat'], 'corr_data_population');
% data_x_axis = corr_data_population.CS1xCS2_prob(:,50);
% win_len = '1ms';
% load(['corr_data_population_' win_len '.mat'], 'corr_data_population');
% data_y_axis = corr_data_population.SS1xSS2_prob(:,50);

mdl = fitlm(data_x_axis,data_y_axis, 'linear');
data_x_axis_hat = sort(data_x_axis);
data_y_axis_hat = polyval(flip(mdl.Coefficients.Estimate)', data_x_axis_hat);
txt_mdl = {[num2str(mdl.Coefficients.Estimate(2)) '*X + ' num2str(mdl.Coefficients.Estimate(1))],...
    ['R-sqr= ' num2str(mdl.Rsquared.Ordinary)], ...
    ['p-value= ' num2str(coefTest(mdl))]};

hFig = figure(5);
clf(hFig);
hold on;
plot(data_x_axis, data_y_axis, 'ob')
plot(data_x_axis_hat, data_y_axis_hat, '-r')
text(mean(data_x_axis),mean(data_y_axis),txt_mdl)

xlabel('CSxCS')
ylabel('SSxSS')
ESN_Beautify_Plot(hFig, [3 3], 8)
% ESN_Beautify_Plot(hFig, [5 5], 12)




end

%% function exclude_overlap_CS(corr_data_population)
function corr_data_population = exclude_overlap_CS(corr_data_population)
fprintf(['exclude_overlap_CS' ' ...'])
num_pairs = size(corr_data_population.SS1xSS1_num_spikes, 1);
for counter_pair = 1 : num_pairs
    CS1xSS2_prob = corr_data_population.CS1xSS2_prob{counter_pair, 1};
    CS1xCS2_prob = corr_data_population.CS1xCS2_prob{counter_pair, 1};
    inds_overlap_CS1xCS2 = find(CS1xCS2_prob(:,50));
    if ~isempty(inds_overlap_CS1xCS2)
        CS1xSS2_prob(inds_overlap_CS1xCS2, :) = [];
%         CS1xCS2_prob(inds_overlap_CS1xCS2, :) = [];
    end
    corr_data_population.CS1xSS2_prob{counter_pair, 1} = CS1xSS2_prob;
%     corr_data_population.CS1xCS2_prob{counter_pair, 1} = CS1xCS2_prob;
    
    CS2xSS1_prob = corr_data_population.CS2xSS1_prob{counter_pair, 1};
    CS2xCS1_prob = corr_data_population.CS2xCS1_prob{counter_pair, 1};
    inds_overlap_CS2xCS1 = find(CS2xCS1_prob(:,50));
    if ~isempty(inds_overlap_CS2xCS1)
        CS2xSS1_prob(inds_overlap_CS2xCS1, :) = [];
%         CS2xCS1_prob(inds_overlap_CS2xCS1, :) = [];
    end
    corr_data_population.CS2xSS1_prob{counter_pair, 1} = CS2xSS1_prob;
%     corr_data_population.CS2xCS1_prob{counter_pair, 1} = CS2xCS1_prob;
end
fprintf(' --> Completed. \n')
end

%% function average_prob(corr_data_population)
function corr_data_population = average_prob(corr_data_population)
fprintf(['average_prob' ' ...'])
ESN_global_variables();
global corr_data_list expand_index
num_pairs = size(corr_data_population.SS1xSS1_num_spikes, 1);
for counter_variable1 = 1 : length(corr_data_list)
    var_name_1 = corr_data_list{counter_variable1};
    for counter_variable2 = 1 : length(corr_data_list)
        var_name_2 = corr_data_list{counter_variable2};
        prob_cell = corr_data_population.([var_name_2 'x' var_name_1 '_prob']);
        
        prob_mat = zeros(num_pairs, 100);
        for counter_pair = 1 : num_pairs
            prob_ = prob_cell{counter_pair, 1};
            prob_ = ESN_expand_index_event_data(prob_, 2, expand_index, 'centered'); % dim=2, expand along row. prob_ is a nx100 matrix
            if size(prob_, 1) > 1
                prob_ = nanmean(prob_);
            end
            prob_mat(counter_pair, :) = prob_;
        end
        corr_data_population.([var_name_2 'x' var_name_1 '_prob']) = prob_mat;
        
        prob_base = corr_data_population.([var_name_2 'x' var_name_1 '_prob_base']);
        prob_base = prob_base .* ((expand_index*2)+1);
        corr_data_population.([var_name_2 'x' var_name_1 '_prob_base']) = prob_base;
        
        prob_norm = prob_mat ./ repmat(prob_base, 1, size(prob_mat, 2));
        corr_data_population.([var_name_2 'x' var_name_1 '_prob_norm']) = prob_norm;
    end
end
fprintf(' --> Completed. \n')
end

%% function sync_index(corr_data_population)
function corr_data_population = sync_index(corr_data_population)
fprintf(['average_prob' ' ...'])
ESN_global_variables();
global corr_data_list expand_index
num_pairs = size(corr_data_population.SS1xSS1_num_spikes, 1);

SS1xSS1_prob = corr_data_population.(['SS1' 'x' 'SS1' '_prob']);
SS1xSS2_prob = corr_data_population.(['SS1' 'x' 'SS2' '_prob']);
SS2xSS1_prob = corr_data_population.(['SS2' 'x' 'SS1' '_prob']);
SS2xSS2_prob = corr_data_population.(['SS2' 'x' 'SS2' '_prob']);
SS1xSS1_prob_base = corr_data_population.(['SS1' 'x' 'SS1' '_prob_base']);
SS2xSS2_prob_base = corr_data_population.(['SS2' 'x' 'SS2' '_prob_base']);


for counter_pair = 1 : num_pairs
    %%
    prob_SS1_given_SS1 = ESN_expand_index_event_data(SS1xSS1_prob{counter_pair, 1}, 2, expand_index, 'centered'); % dim=2, expand along row. prob_ is a nx100 matrix
    prob_SS2_given_SS1 = ESN_expand_index_event_data(SS1xSS2_prob{counter_pair, 1}, 2, expand_index, 'centered'); % dim=2, expand along row. prob_ is a nx100 matrix
    prob_SS1_given_SS2 = ESN_expand_index_event_data(SS2xSS1_prob{counter_pair, 1}, 2, expand_index, 'centered'); % dim=2, expand along row. prob_ is a nx100 matrix
    prob_SS2_given_SS2 = ESN_expand_index_event_data(SS2xSS2_prob{counter_pair, 1}, 2, expand_index, 'centered'); % dim=2, expand along row. prob_ is a nx100 matrix
    prob_SS1 = SS1xSS1_prob_base(counter_pair, 1) .* ((expand_index*2)+1);
    prob_SS2 = SS2xSS2_prob_base(counter_pair, 1) .* ((expand_index*2)+1);
    %%
    % sync_index_SS1_SS2 = P(SS1,SS2) / (P(SS1) * P(SS2))
    % P(SS1,SS2) = ( P(SS1,SS2|SS1) * P(SS1) + P(SS1,SS2|SS2) * P(SS2) ) / (P(SS1) + P(SS2))
    % P(SS1,SS2|SS1) = P(SS1|SS1) AND P(SS2|SS1)
    % P(SS1,SS2|SS2) = P(SS1|SS2) AND P(SS2|SS2)
    % P(SS1) = firing_rate SS1, P(SS2) = firing_rate SS2
    % ALTERNETIVELY
    % sync_index_SS1_SS2 = P(SS1,SS2) / (P(SS1) * P(SS2))
    % sync_index_SS1_SS2_1 = (P(SS1|SS2) * P(SS2)) / (P(SS1) * P(SS2)) = P(SS1|SS2) / P(SS1)
    % sync_index_SS1_SS2_2 = (P(SS2|SS1) * P(SS1)) / (P(SS1) * P(SS2)) = P(SS2|SS1) / P(SS2)
    % sync_index_SS1_SS2 = (sync_index_SS1_SS2_1+sync_index_SS1_SS2_2) / 2
    
    prob_SS1_and_SS2_given_SS1 = nanmean( prob_SS1_given_SS1 & prob_SS2_given_SS1 );
    prob_SS1_and_SS2_given_SS2 = nanmean( prob_SS1_given_SS2 & prob_SS2_given_SS2 );
    prob_SS1_and_SS2 = ( (prob_SS1_and_SS2_given_SS1 .* prob_SS1) + (prob_SS1_and_SS2_given_SS2 .* prob_SS2) ) ./ (prob_SS1 + prob_SS2);
    sync_index_SS1_and_SS2_v1 = prob_SS1_and_SS2 ./ (prob_SS1 .* prob_SS2);
    
    
    sync_index_SS1_and_SS2_v2 = ( (nanmean(prob_SS1_given_SS2)./prob_SS1) + (nanmean(prob_SS2_given_SS1)./prob_SS2)) / 2;
    
end


for counter_variable1 = 1 : length(corr_data_list)
    var_name_1 = corr_data_list{counter_variable1};
    for counter_variable2 = 1 : length(corr_data_list)
        var_name_2 = corr_data_list{counter_variable2};
        prob_cell = corr_data_population.([var_name_2 'x' var_name_1 '_prob']);
        
        prob_mat = zeros(num_pairs, 100);
        for counter_pair = 1 : num_pairs
            prob_ = prob_cell{counter_pair, 1};
            prob_ = ESN_expand_index_event_data(prob_, 2, expand_index, 'centered'); % dim=2, expand along row. prob_ is a nx100 matrix
            if size(prob_, 1) > 1
                prob_ = nanmean(prob_);
            end
            prob_mat(counter_pair, :) = prob_;
        end
        corr_data_population.([var_name_2 'x' var_name_1 '_prob']) = prob_mat;
        
        prob_base = corr_data_population.([var_name_2 'x' var_name_1 '_prob_base']);
        prob_base = prob_base .* ((expand_index*2)+1);
        corr_data_population.([var_name_2 'x' var_name_1 '_prob_base']) = prob_base;
        
        prob_norm = prob_mat ./ repmat(prob_base, 1, size(prob_mat, 2));
        corr_data_population.([var_name_2 'x' var_name_1 '_prob_norm']) = prob_norm;
        
        
        
        
    end
end
fprintf(' --> Completed. \n')
end
