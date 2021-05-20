function ESN_realign_sorted_CS_SS
%% load EPHYS sorted DATA
[file_name,file_path] = uigetfile([pwd filesep '*_sorted*.mat'], 'Select _sorted file');
fprintf(['Loading ', file_name, ' ... ']);
load([file_path filesep file_name]);
EPHYS.CH_sorted_file_name = file_name;
EPHYS.CH_sorted_file_path = file_path;
[~, file_name, ~]  = fileparts(EPHYS.CH_sorted_file_name);
EPHYS.file_name = file_name;
fprintf(' --> Completed. \n')

%% Init
ch_data = struct;
ch_data.ch_data_time   = ch_data_time;
ch_data.ch_data_hipass = ch_data_hipass;
ch_data.ch_data_lopass = ch_data_lopass;
ch_data.ch_data_cmn    = ch_data_cmn;

CS_data = struct;
CS_data.CS_num_total         = length(CS_time);
if CS_data.CS_num_total > 1
CS_data.CS_time              = CS_time;
CS_data.CS_ind               = nan(CS_data.CS_num_total, 1);
CS_data.CS_inds              = nan(CS_data.CS_num_total, 180);
CS_data.CS_waveform_hipass   = nan(CS_data.CS_num_total, 180);
CS_data.CS_peak_value_hipass = nan(CS_data.CS_num_total, 1);
CS_data.CS_prominence_hipass = nan(CS_data.CS_num_total, 1);
CS_data.CS_width_hipass      = nan(CS_data.CS_num_total, 1);
CS_data.CS_peak_value_lopass = CS_peak_value;
CS_data.CS_prominence_lopass = CS_prominence;
CS_data.CS_width_lopass      = CS_width;
end

SS_data = struct;
SS_data.SS_num_total         = length(SS_time);
SS_data.SS_time              = SS_time;
SS_data.SS_ind               = nan(SS_data.SS_num_total, 1);
SS_data.SS_inds              = nan(SS_data.SS_num_total, 180);
SS_data.SS_waveform_hipass   = nan(SS_data.SS_num_total, 180);
SS_data.SS_peak_value_hipass = SS_peak_value;
SS_data.SS_prominence_hipass = SS_prominence;
SS_data.SS_width_hipass      = SS_width;

%% CS Align
fprintf([EPHYS.file_name ': Aligning CS', ' ... ']);
clearvars -except ch_data CS_data SS_data EPHYS
data_length          = length(ch_data.ch_data_time);
CS_num_total         = CS_data.CS_num_total;
if CS_data.CS_num_total > 1
CS_time_hipass       = nan(CS_num_total, 1);
CS_ind               = nan(CS_num_total, 1);
CS_inds              = nan(CS_num_total, 180);
CS_waveform_hipass   = nan(CS_num_total, 180);
CS_peak_value_hipass = nan(CS_num_total, 1);
CS_prominence_hipass = nan(CS_num_total, 1);
CS_width_hipass      = nan(CS_num_total, 1);

for counter_CS = 1 : 1 : CS_num_total
    if counter_CS == 1
        CS_interval_begin_ = find(ch_data.ch_data_time >= CS_data.CS_time(counter_CS, 1), 1, 'first');
        CS_interval_begin_ = max([(CS_interval_begin_ - 30e3) 1]);
    else
        CS_interval_begin_ = CS_ind(counter_CS-1);
    end
    CS_interval_end_ = min([(CS_interval_begin_ + 30e4) data_length]);
    CS_interval_inds = CS_interval_begin_ : 1 : CS_interval_end_;
    CS_interval_time = ch_data.ch_data_time(CS_interval_inds);
    CS_interval_hipass = ch_data.ch_data_hipass(CS_interval_inds);
    
    CS_ind_ = find(CS_interval_time >= CS_data.CS_time(counter_CS, 1), 1, 'first');
    
    if isempty(CS_ind_)
        CS_interval_begin_ = find(ch_data.ch_data_time >= CS_data.CS_time(counter_CS, 1), 1, 'first');
        CS_interval_begin_ = max([(CS_interval_begin_ - 30e3) 1]);
        CS_interval_end_ = min([(CS_interval_begin_ + 30e4) data_length]);
        CS_interval_inds = CS_interval_begin_ : 1 : CS_interval_end_;
        CS_interval_time = ch_data.ch_data_time(CS_interval_inds);
        CS_interval_hipass = ch_data.ch_data_hipass(CS_interval_inds);
    
        CS_ind_ = find(CS_interval_time >= CS_data.CS_time(counter_CS, 1), 1, 'first');
    end
    
    CS_inds_ = (CS_ind_ - 149) : 1 : (CS_ind_ + 150);
    CS_waveform_hipass_ = CS_interval_hipass(CS_inds_, 1)';
    CS_waveform_time_ = CS_interval_time(CS_inds_, 1)';
    [CS_peak_value_hipass_, CS_time_hipass_, CS_width_hipass_, CS_prominence_hipass_] = findpeaks(CS_waveform_hipass_ , CS_waveform_time_, 'SortStr', 'descend', 'NPeaks', 1);%'MinPeakHeight',threshold_CS,'MinPeakProminence',threshold_CS,'MinPeakDistance',threshold_CS_time);
    CS_ind_local_ = find(CS_waveform_time_ >= CS_time_hipass_, 1, 'first');
    ind_upcross = find(CS_waveform_hipass_(1:CS_ind_local_-1) <= 0 & CS_waveform_hipass_(2:CS_ind_local_) > 0, 1, 'last');
    CS_time_hipass_ = CS_waveform_time_(ind_upcross);
    CS_ind_  = CS_inds_(ind_upcross);
    CS_inds_ = (CS_ind_ - 59) : 1 : (CS_ind_ + 120);
    CS_waveform_hipass_ = CS_interval_hipass(CS_inds_, 1)';
    CS_ind_  = CS_interval_inds(CS_ind_);
    CS_inds_ = CS_interval_inds(CS_inds_);
    
    CS_ind(              counter_CS, :) = CS_ind_;
    CS_peak_value_hipass(counter_CS, :) = CS_peak_value_hipass_;
    CS_time_hipass(      counter_CS, :) = CS_time_hipass_;
    CS_width_hipass(     counter_CS, :) = CS_width_hipass_;
    CS_prominence_hipass(counter_CS, :) = CS_prominence_hipass_;
    CS_waveform_hipass(  counter_CS, :) = CS_waveform_hipass_;
    CS_inds(             counter_CS, :) = CS_inds_;
    
    if mod(counter_CS, 1000) == 0
        disp([num2str(counter_CS) ' / ' num2str(CS_num_total)])
    end
end

CS_data.CS_peak_value_hipass = CS_peak_value_hipass;
CS_data.CS_prominence_hipass = CS_prominence_hipass;
CS_data.CS_width_hipass      = CS_width_hipass;
CS_data.CS_time              = CS_time_hipass;
CS_data.CS_ind               = CS_ind;
CS_data.CS_inds              = CS_inds;
CS_data.CS_waveform_hipass   = CS_waveform_hipass;

end
fprintf(' --> Completed. \n')

%% SS Align
fprintf([EPHYS.file_name ': Aligning SS', ' ... ']);
clearvars -except ch_data CS_data SS_data EPHYS
data_length          = length(ch_data.ch_data_time);
SS_num_total         = SS_data.SS_num_total;
SS_ind               = nan(SS_num_total, 1);
SS_waveform_hipass   = nan(SS_num_total, 180);
SS_inds              = nan(SS_num_total, 180);

for counter_SS = 1 : 1 : SS_num_total
    if counter_SS == 1
        SS_interval_begin_ = find(ch_data.ch_data_time >= SS_data.SS_time(counter_SS, 1), 1, 'first');
        SS_interval_begin_ = max([(SS_interval_begin_ - 30e3) 1]);
    else
        SS_interval_begin_ = SS_ind(counter_SS-1);
    end
    SS_interval_end_ = min([(SS_interval_begin_ + 15e4) data_length]);
    SS_interval_inds = SS_interval_begin_ : 1 : SS_interval_end_;
    SS_interval_time = ch_data.ch_data_time(SS_interval_inds);
    SS_interval_hipass = ch_data.ch_data_hipass(SS_interval_inds);
    
    SS_ind_ = find(SS_interval_time(1:end-120) >= SS_data.SS_time(counter_SS, 1), 1, 'first');
    
    if isempty(SS_ind_)
        SS_interval_begin_ = find(ch_data.ch_data_time >= SS_data.SS_time(counter_SS, 1), 1, 'first');
        SS_interval_begin_ = max([(SS_interval_begin_ - 30e3) 1]);
        SS_interval_end_ = min([(SS_interval_begin_ + 15e4) data_length]);
        SS_interval_inds = SS_interval_begin_ : 1 : SS_interval_end_;
        SS_interval_time = ch_data.ch_data_time(SS_interval_inds);
        SS_interval_hipass = ch_data.ch_data_hipass(SS_interval_inds);
    
        SS_ind_ = find(SS_interval_time >= SS_data.SS_time(counter_SS, 1), 1, 'first');
    end
    
    SS_inds_ = (SS_ind_ - 59) : 1 : (SS_ind_ + 120);
    if counter_SS == 1
        SS_inds_(SS_inds_<1) = 1;
    end
    if counter_SS == SS_num_total
        SS_inds_(SS_inds_>length(SS_interval_inds)) = length(SS_interval_inds);
    end
    SS_waveform_hipass_ = SS_interval_hipass(SS_inds_, 1)';
    SS_ind_  = SS_interval_inds(SS_ind_);
    SS_inds_ = SS_interval_inds(SS_inds_);
    
    SS_ind(              counter_SS, :) = SS_ind_;
    SS_inds(             counter_SS, :) = SS_inds_;
    SS_waveform_hipass(  counter_SS, :) = SS_waveform_hipass_;
    
    if mod(counter_SS, 1000) == 0
        disp([num2str(counter_SS) ' / ' num2str(SS_num_total)])
    end
end

SS_data.SS_ind               = SS_ind;
SS_data.SS_inds              = SS_inds;
SS_data.SS_waveform_hipass   = SS_waveform_hipass;
fprintf(' --> Completed. \n')

%% Save ch_data CS_data SS_data
clearvars -except ch_data CS_data SS_data EPHYS
fprintf([EPHYS.file_name ': Saving ch_data CS_data SS_data ...'])
save([EPHYS.CH_sorted_file_path EPHYS.CH_sorted_file_name], 'ch_data', 'CS_data', 'SS_data', '-v7.3');
fprintf(' --> Completed. \n')

%% Plot WaveForm
fig_handle_(1) = figure(1);
clf(fig_handle_(1))
plot_handle_(1) = subplot(1,2,1);
hold on
if CS_data.CS_num_total > 1
plot((1:180)/30, mean(CS_data.CS_waveform_hipass), 'LineWidth', 2, 'Color', [0.0 0.0 0.0])
plot((1:180)/30, mean(CS_data.CS_waveform_hipass)+std(CS_data.CS_waveform_hipass), 'LineWidth', 2, 'Color', [0.5 0.5 0.5])
plot((1:180)/30, mean(CS_data.CS_waveform_hipass)-std(CS_data.CS_waveform_hipass), 'LineWidth', 2, 'Color', [0.5 0.5 0.5])
end

plot_handle_(2) = subplot(1,2,2);
hold on
plot((1:180)/30, mean(SS_data.SS_waveform_hipass), 'LineWidth', 2, 'Color', [0.0 0.0 0.0])
plot((1:180)/30, mean(SS_data.SS_waveform_hipass)+std(SS_data.SS_waveform_hipass), 'LineWidth', 2, 'Color', [0.5 0.5 0.5])
plot((1:180)/30, mean(SS_data.SS_waveform_hipass)-std(SS_data.SS_waveform_hipass), 'LineWidth', 2, 'Color', [0.5 0.5 0.5])

ylim_plot_handle_1 = ylim(plot_handle_(1));
ylim_plot_handle_2 = ylim(plot_handle_(2));
ylim_plots = [min([ylim_plot_handle_1(1) ylim_plot_handle_2(1)]) max([ylim_plot_handle_1(2) ylim_plot_handle_2(2)])];
ylim(plot_handle_(1), ylim_plots);
ylim(plot_handle_(2), ylim_plots);

ESN_Beautify_Plot

%% Plot ISI
fig_handle_(2) = figure(2);
clf(fig_handle_(2))
plot_handle_(1) = subplot(1,2,1);
hold on
if CS_data.CS_num_total > 1
edges_CS = 0 : 0.200 : 5.000;
ISI_CS = diff(CS_data.CS_time);
histogram(ISI_CS,edges_CS, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', 'k');
set(plot_handle_(1), 'XTick', [0 2.5 5.0])
end

plot_handle_(2) = subplot(1,2,2);
hold on
% edges_SS = 0 : 0.001 : 0.050;
edges_SS = (0 : 0.002 : 0.050) *1000;
ISI_SS = diff(SS_data.SS_time) * 1000;
histogram(ISI_SS,edges_SS, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', 'k');
set(plot_handle_(2), 'XTick', [0 0.025 0.050]*1000)

ESN_Beautify_Plot

%% Plot Cs x SS correlation
CS_time = round((CS_data.CS_time) * 1000);
SS_time = round((SS_data.SS_time) * 1000);
min_time =  min([CS_time(1) SS_time(1)]) - 1;
CS_time = CS_time - min_time;
SS_time = SS_time - min_time;
max_time =  max([CS_time(end) SS_time(end)]);
CS_vec = zeros(max_time, 1);
SS_vec = zeros(max_time, 1);
CS_vec(CS_time)   = 1;
SS_vec(SS_time)   = 1;
[xcorr_value_CS_SS,xcorr_lag_CS_SS] = xcorr(SS_vec, CS_vec, 50); % cross-correlate signals with each other
[xcorr_value_SS_SS,xcorr_lag_SS_SS] = xcorr(SS_vec, 50); % cross-correlate signals with each other
xcorr_value_SS_SS((xcorr_lag_SS_SS==0)) = 0;

% xcorr_value_CS_SS = round(smooth(xcorr_value_CS_SS));
% xcorr_value_CS_SS = xcorr_value_CS_SS - min(xcorr_value_CS_SS);

fig_handle_(3) = figure(3);
clf(fig_handle_(3))
plot_handle_(1) = subplot(1,2,1);
plot(xcorr_lag_CS_SS, (xcorr_value_CS_SS), 'LineWidth', 2, 'Color', [0.0 0.0 0.0]); 
ylim([0 max(xcorr_value_CS_SS)])
xlim([-50 50]);
set(plot_handle_(1), 'XTick', -50:25:50);
ylabel('Cross-Correlogram');
xlabel('Time (ms)');
plot_handle_(2) = subplot(1,2,2);
plot(xcorr_lag_SS_SS, (xcorr_value_SS_SS), 'LineWidth', 2, 'Color', [0.0 0.0 0.0]); 
xlim([-50 50]);
set(plot_handle_(2), 'XTick', -50:25:50);
ylabel('Cross-Correlogram');
xlabel('Time (ms)');
ESN_Beautify_Plot

%% Save Fig
fprintf([EPHYS.file_name ': Saving plot CS SS', ' ...'])
saveas(fig_handle_(1),[EPHYS.CH_sorted_file_path EPHYS.file_name '_CS_SS_Waveform'], 'pdf');
saveas(fig_handle_(1),[EPHYS.CH_sorted_file_path EPHYS.file_name '_CS_SS_Waveform'], 'png');
saveas(fig_handle_(2),[EPHYS.CH_sorted_file_path EPHYS.file_name '_CS_SS_ISI'], 'pdf');
saveas(fig_handle_(2),[EPHYS.CH_sorted_file_path EPHYS.file_name '_CS_SS_ISI'], 'png');
saveas(fig_handle_(3),[EPHYS.CH_sorted_file_path EPHYS.file_name '_CS_SS_Correlogram'], 'pdf');
saveas(fig_handle_(3),[EPHYS.CH_sorted_file_path EPHYS.file_name '_CS_SS_Correlogram'], 'png');
close(fig_handle_(1))
close(fig_handle_(2))
close(fig_handle_(3))
fprintf(' --> Completed. \n')



