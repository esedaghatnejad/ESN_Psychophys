% load EPHYS sorted DATA
file_path = pwd;
[file_name, file_path] = uigetfile([file_path filesep '*.psort'], 'Select psort file');
fprintf(['Loading ', file_name, ' ... ']);
% EPHYS.CH_sorted = load([file_path filesep file_name], 'CS_data', 'SS_data');
DATA_PSORT = Psort_read_psort([file_path filesep file_name]);
EPHYS.CH_sorted_file_name = file_name;
EPHYS.CH_sorted_file_path = file_path;
fprintf(' --> Completed. \n')

% build EPHYS.CH_sorted from DATA_PSORT
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

EPHYS.CH_sorted.SS_data.SS_ind = SS_index;
EPHYS.CH_sorted.CS_data.CS_ind = CS_index;
EPHYS.CH_sorted.SS_data.SS_time = SS_time;
EPHYS.CH_sorted.CS_data.CS_time = CS_time;
EPHYS.CH_sorted.SS_data.SS_waveform = SS_waveform;
EPHYS.CH_sorted.CS_data.CS_waveform = CS_waveform;

% Report Properties
[~, file_name, ~]  = fileparts(EPHYS.CH_sorted_file_name);
duration = (DATA_PSORT.topLevel_data.ch_time(end)-DATA_PSORT.topLevel_data.ch_time(1));
numCS = length(EPHYS.CH_sorted.CS_data.CS_ind);
freqCS = numCS/duration;
numSS = length(EPHYS.CH_sorted.SS_data.SS_ind);
freqSS = numSS/duration;
numTrial = 0;
fprintf(['*******************************************' '\n'])
fprintf([file_name '\n'])
fprintf([...
    'ID' '\t'...
    'dur' '\t'...
    'numCS' '\t'...
    'freqCS' '\t'...
    'numSS' '\t'...
    'freqSS' '\t'...
    'numTrial' '\n'...
    ])
fprintf([...
    file_name '\t'...
    num2str(duration/60,'%.1f') '\t'...
    num2str(numCS,'%.0f') '\t'...
    num2str(freqCS,'%.2f') '\t'...
    num2str(numSS,'%.0f') '\t'...
    num2str(freqSS,'%.2f') '\t'...
    num2str(numTrial,'%.0f') '\n'...
    ])



