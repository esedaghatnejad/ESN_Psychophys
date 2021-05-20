file_path = '/home/kkarbasi/cssorter/125d_data_sorted/2019-04/2019-04-10/2019-04-10_14-45-41/analyzed_data';
file_name = '190410_144541_04_sorted.mat';
fprintf(['Loading ', file_name, ' ... ']);
load([file_path filesep file_name]);
fprintf(' --> Completed. \n')

%%
fprintf(['Power Spectrum ', ' ... ']);
CS_waveform_cmn = CH_data.CH_data_cmn(CS_data.CS_inds);
SS_waveform_cmn = CH_data.CH_data_cmn(SS_data.SS_inds);

% CS_waveform_cmn = CH_data.CH_hipass(CS_data.CS_inds);
% SS_waveform_cmn = CH_data.CH_hipass(SS_data.SS_inds);

[pxx_CS_waveform_cmn, freq] = pwelch(CS_waveform_cmn', [], [], [], 30e3);
[pxx_SS_waveform_cmn, ~] = pwelch(SS_waveform_cmn', [], [], [], 30e3);
fprintf(' --> Completed. \n')

%%
clf(figure(1))

y_axis_SS = pow2db(pxx_SS_waveform_cmn);
hold on
plot(freq, mean(y_axis_SS, 2)+std(y_axis_SS, 0, 2), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot(freq, mean(y_axis_SS, 2)-std(y_axis_SS, 0, 2), ...
    'LineWidth', 1, 'Color', [0.5 0.5 0.9])
plot(freq, mean(y_axis_SS, 2), ...
    'LineWidth', 2, 'Color', [0.1 0.1 0.9])

y_axis_CS = pow2db(pxx_CS_waveform_cmn);
hold on
plot(freq, mean(y_axis_CS, 2)+std(y_axis_CS, 0, 2), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot(freq, mean(y_axis_CS, 2)-std(y_axis_CS, 0, 2), ...
    'LineWidth', 1, 'Color', [0.9 0.5 0.5])
plot(freq, mean(y_axis_CS, 2), ...
    'LineWidth', 2, 'Color', [0.9 0.1 0.1])
set(gca, 'XScale', 'log')
xlabel('Freq (Hz)')
ylabel('Power/Freq (dB/Hz)')
title('Welch Power Spectral Density Estimate')

%%
clf(figure(2))

y_axis_CS = pow2db(pxx_CS_waveform_cmn);
hold on
plot(freq, y_axis_CS)

set(gca, 'XScale', 'log')
xlabel('Freq (Hz)')
ylabel('Power/Freq (dB/Hz)')
title('Welch Power Spectral Density Estimate')

%%
fprintf(['Building SSxSS_CROSS PROBABILITY ' ' ...'])
SS_time   = SS_data.SS_time;
CS_time   = CS_data.CS_time;
Corr_data = ESN_correlogram(SS_time, CS_time); % ESN_correlogram(SS_time, CS_time) 
% Corr_data.SS_inds_span; Corr_data.SS_bin_size_time; Corr_data.SS_SSxSS_AUTO; Corr_data.CS_inds_span; Corr_data.CS_bin_size_time; Corr_data.CS_CSxSS_AUTO;


fprintf(' --> Completed. \n')

%% function ESN_correlogram
function Corr_data = ESN_correlogram(SS_time, CS_time)
bin_size_time = 1e-3; % seconds
span_window_size = (1 / bin_size_time) * (100 / 1000);
span_window_size_half = round(span_window_size / 2);
inds_span = ((-span_window_size_half+1) : 1 : (span_window_size_half))';

if (~isempty(CS_time)) && (~isempty(SS_time))
    ch_time_min = min([SS_time(1) CS_time(1)]);
    ch_time_min = max([(ch_time_min-2.0) 0]);
    ch_time_max = max([SS_time(end) CS_time(end)]) + 2.0;
    
    CH__.SS_data.SS_time =  SS_time - ch_time_min;
    CH__.CS_data.CS_time =  CS_time - ch_time_min;
    CH__.SS_data.SS_time = ESN_Round(CH__.SS_data.SS_time, bin_size_time);
    CH__.CS_data.CS_time = ESN_Round(CH__.CS_data.CS_time, bin_size_time);
    CH__.SS_data.SS_ind  = round(CH__.SS_data.SS_time .* (1/bin_size_time));
    CH__.SS_data.SS_ind( CH__.SS_data.SS_ind < 1 ) = 1;
    CH__.CS_data.CS_ind  = round(CH__.CS_data.CS_time .* (1/bin_size_time));
    CH__.CS_data.CS_ind( CH__.CS_data.CS_ind < 1 ) = 1;
    CH__.CH_data.ch_time_reconstruct = ch_time_min : bin_size_time : ch_time_max;
elseif (~isempty(CS_time))
    ch_time_min = min(  CS_time(1) );
    ch_time_min = max([(ch_time_min-2.0) 0]);
    ch_time_max = max(  CS_time(end) ) + 2.0;
    
    CH__.CS_data.CS_time =  CS_time - ch_time_min;
    CH__.CS_data.CS_time = ESN_Round(CH__.CS_data.CS_time, bin_size_time);
    CH__.CS_data.CS_ind  = round(CH__.CS_data.CS_time .* (1/bin_size_time));
    CH__.CS_data.CS_ind( CH__.CS_data.CS_ind < 1 ) = 1;
    CH__.CH_data.ch_time_reconstruct = ch_time_min : bin_size_time : ch_time_max;
elseif (~isempty(SS_time))
    ch_time_min = min( SS_time(1)  );
    ch_time_min = max([(ch_time_min-2.0) 0]);
    ch_time_max = max( SS_time(end)  ) + 2.0;
    
    CH__.SS_data.SS_time =  SS_time - ch_time_min;
    CH__.SS_data.SS_time = ESN_Round(CH__.SS_data.SS_time, bin_size_time);
    CH__.SS_data.SS_ind  = round(CH__.SS_data.SS_time .* (1/bin_size_time));
    CH__.SS_data.SS_ind( CH__.SS_data.SS_ind < 1 ) = 1;
    CH__.CH_data.ch_time_reconstruct = ch_time_min : bin_size_time : ch_time_max;
end

% SSxSS_AUTO
if (~isempty(SS_time))
    CH__.SS_data.SS_inds_reconstruct = repmat( CH__.SS_data.SS_ind(:), 1, length(inds_span)) + repmat(inds_span(:)', length(CH__.SS_data.SS_ind), 1);
    CH__.SS_data.SS_inds_reconstruct( CH__.SS_data.SS_inds_reconstruct < 1 ) = 1;
    CH__.SS_data.SS_inds_reconstruct( CH__.SS_data.SS_inds_reconstruct > length( CH__.CH_data.ch_time_reconstruct ) ) = length( CH__.CH_data.ch_time_reconstruct );
    
    CH__.SS_data.SS_event_trace = false( size(CH__.CH_data.ch_time_reconstruct) );
    CH__.SS_data.SS_event_trace( CH__.SS_data.SS_ind ) = true ;
    CH__.SS_data.SS_event_trace( 1   ) = false;
    CH__.SS_data.SS_event_trace( end ) = false;
    
    CH__.SS_data.SS_event_reconstruct = CH__.SS_data.SS_event_trace( CH__.SS_data.SS_inds_reconstruct );
    % SSxSS correlogram
    SSxSS_AUTO       = CH__.SS_data.SS_event_reconstruct;
    ss_inds_span     = repmat(inds_span(:)',     size(SS_time(:),1), 1);
    ss_bin_size_time = repmat(bin_size_time(:)', size(SS_time(:),1), 1);
else
    SSxSS_AUTO       = false(0, length(inds_span(:)'));
    ss_inds_span     = nan(0, length(inds_span(:)'));
    ss_bin_size_time = nan(0, 1);
end

% CSxSS_WITHIN
if (~isempty(CS_time)) && (~isempty(SS_time))
    CH__.CS_data.CS_inds_reconstruct = repmat( CH__.CS_data.CS_ind(:), 1, length(inds_span)) + repmat(inds_span(:)', length(CH__.CS_data.CS_ind), 1);
    CH__.CS_data.CS_inds_reconstruct( CH__.CS_data.CS_inds_reconstruct < 1 ) = 1;
    CH__.CS_data.CS_inds_reconstruct( CH__.CS_data.CS_inds_reconstruct > length( CH__.CH_data.ch_time_reconstruct ) ) = length( CH__.CH_data.ch_time_reconstruct );
    
    CH__.SS_data.SS_event_trace = false( size(CH__.CH_data.ch_time_reconstruct) );
    CH__.SS_data.SS_event_trace( CH__.SS_data.SS_ind ) = true ;
    CH__.SS_data.SS_event_trace( 1   ) = false;
    CH__.SS_data.SS_event_trace( end ) = false;
    
    CH__.SS_data.SS_event_reconstruct = CH__.SS_data.SS_event_trace( CH__.CS_data.CS_inds_reconstruct );
    % CSxSS correlogram
    CSxSS_AUTO       = CH__.SS_data.SS_event_reconstruct;
    cs_inds_span     = repmat(inds_span(:)',     size(CS_time(:),1), 1);
    cs_bin_size_time = repmat(bin_size_time(:)', size(CS_time(:),1), 1);
else
    CSxSS_AUTO       = false(0, length(inds_span(:)'));
    cs_inds_span     = nan(0, length(inds_span(:)'));
    cs_bin_size_time = nan(0, 1);
end

Corr_data = struct;
Corr_data.CS_inds_span     = cs_inds_span;
Corr_data.CS_bin_size_time = cs_bin_size_time;
Corr_data.SS_inds_span     = ss_inds_span;
Corr_data.SS_bin_size_time = ss_bin_size_time;
Corr_data.SS_SSxSS_AUTO    = SSxSS_AUTO;
Corr_data.CS_CSxSS_AUTO    = CSxSS_AUTO;
end
