function ESN_smr2mat(folder_path)
% Author: Ehsan Sedaghat-Nejad (esedaghatnejad@gmail.com)
% this function get the folder path and convert all the smr files in that folder to mat files
% This function works based on CED 64-bit Son library, as used by Spike2. This library can read and 
% create both the older format 32-bit Son files (smr file extension) and the newer 64-bit files 
% (smrx file extension).
% The library is based on ceds64int.dll and is exclusive for Windows OS
% You need to change CEDS64ML_path in the initialize_CEDS64ML() so the library can be loaded and
% used properly. The code will search for all available channels in the smr file and save the data
% is a structural format.
% Please read "s64mat.pdf" or "s64mat.chm" for further information.


%% Initialise the library
initialize_CEDS64ML();

%% Get the folder_path
% if there is no inputs, then set folder_path to pwd
if nargin < 1
    folder_path = pwd;
end
% add '\' to the end of folder_path if there is none
if ~strcmp(folder_path(end), '\')
    folder_path = [folder_path '\'];
end
% get the list of smr files to analyze
FILES_SMR = dir([folder_path '*.smr']);

%% Loop over files
for counter_smrFiles = 1 : length(FILES_SMR)
    clearvars('SMR_FILE','CH_DATA_ALL');
    file_name = FILES_SMR(counter_smrFiles).name;
    fprintf([num2str(counter_smrFiles) '. ' file_name ' ... '])
    SMR_FILE = open_smr_file(folder_path, file_name);
    CH_DATA_ALL = extract_channels(SMR_FILE);
    save_mat_file(SMR_FILE, CH_DATA_ALL);
    close_smr_file(SMR_FILE);
    fprintf(' --> Done.\n')
end

%% Close all and tidy up
close_all()
fprintf('ALL Done. \n')
end

function initialize_CEDS64ML()
%% Initialise the library
CEDS64ML_path = 'C:\CEDMATLAB\CEDS64ML'; % the location of CED64ML library on this computer
addpath(CEDS64ML_path); % add CED64ML to the Matlab search path
setenv('CEDS64ML', CEDS64ML_path); % set environmental variable CED64ML
CEDS64LoadLib( CEDS64ML_path ); % load ceds64int.dll
end

function SMR_FILE = open_smr_file(folder_path, file_name)
% [SMR_FILE.file_name,SMR_FILE.folder_path] = uigetfile('*.smr');
%% Open a smr file to read
SMR_FILE.file_name = file_name;
SMR_FILE.folder_path = folder_path;
SMR_FILE.full_path = [SMR_FILE.folder_path SMR_FILE.file_name];
SMR_FILE.smr_fhand = CEDS64Open( SMR_FILE.full_path );
if (SMR_FILE.smr_fhand <= 0)
    unloadlibrary ceds64int;
    error(['ESN_smr2mat: Cannot open the smr file: ' SMR_FILE.full_path]);
end
end

function CH_DATA_ALL = extract_channels(SMR_FILE)
%% Loop through channels
clearvars -except SMR_FILE CH_DATA_ALL ch_num ch_counter
SMR_FILE.max_num_channels = CEDS64MaxChan( SMR_FILE.smr_fhand );
SMR_FILE.ch_type_description = {'unused', 'WaveForm', 'EventFall', 'EventRise', 'EventBoth', ...
    'Marker', 'WaveMark', 'RealMark', 'TextMark', 'RealWave'};
ch_counter = 1;
clearvars('CH_DATA_ALL');
for ch_num = 1 : 1 : SMR_FILE.max_num_channels
    %% Clear
    clearvars -except SMR_FILE CH_DATA_ALL ch_num ch_counter
    smr_fhand = SMR_FILE.smr_fhand;
    %% Extract channel data
    ch_type = CEDS64ChanType( smr_fhand, ch_num ); % an integer code indicating the type of a channel in a file
    % Channel types are:
    % 0=channel unused, 1=Adc,      2=EventFall, 3=EventRise, 4=EventBoth,
    % 5=Marker,         6=WaveMark, 7=RealMark,  8=TextMark,  9=RealWave
    if (ch_type <= 0)
        continue; % error in channel or the channel is unused
    end
    [ ~, ch_title ]               = CEDS64ChanTitle(   smr_fhand, ch_num ); % the channel title.
    [ ~, ch_unit ]                = CEDS64ChanUnits(   smr_fhand, ch_num ); % the channel units
    [ ~, ch_comment ]             = CEDS64ChanComment( smr_fhand, ch_num ); % the comment associated with a channel
    [ ~, ch_Y_lower, ch_Y_upper ] = CEDS64ChanYRange(  smr_fhand, ch_num ); % two floating point values for each channel that indicate suggested low and high values for displaying the channel
    [ ch_rate ]                   = CEDS64IdealRate(   smr_fhand, ch_num ); % the ideal rate, in Hz, for a channel. For a Waveform-based channel, this is the rate that the user who created the channel would have liked to sample at (the actual rate available may have been constrained by the data acquisition hardware), and is for information only.
    [ ch_resolution_ticks ]       = 1;
    [ ch_resolution_secs ]        = CEDS64TicksToSecs( smr_fhand, ch_resolution_ticks );
    [ ch_interval_ticks ]         = CEDS64ChanDiv(     smr_fhand, ch_num ); % channel divide (sample interval in ticks) for waveform and WaveMark channels. If used with a channel that does not use a channel divide, the result will be 0.
    [ ch_interval_secs ]          = CEDS64TicksToSecs( smr_fhand, ch_interval_ticks );
    [ ch_max_time_ticks ]         = CEDS64ChanMaxTime( smr_fhand, ch_num ) + 1; % the time in ticks of the last item in a channel, +1 so the read gets the last point
    [ ch_max_time_secs ]          = CEDS64TicksToSecs( smr_fhand, ch_max_time_ticks );
    [ ~, ch_offset ]              = CEDS64ChanOffset(  smr_fhand, ch_num ); % the offset used to convert between 16-bit values and user units for a channel
    [ ~, ch_scale ]               = CEDS64ChanScale(   smr_fhand, ch_num ); % the scale used to convert between 16-bit values and user units for a channel. This is only relevant for Waveform, WaveMark and RealWave channels.
    [ ch_scale ]                  = ch_scale / (2^16) * 10;
    [ ch_max_data_points ] = min([(ceil(1 / ch_interval_secs * ch_max_time_secs)+10) (ceil(1 / 2e-5 * ch_max_time_secs)+10)]);
    %% Waveform
    if (ch_type == 1) || (ch_type == 9)
        [ ch_length, ch_values, ch_start_ticks ] = CEDS64ReadWaveF( smr_fhand, ch_num, ch_max_data_points, 0, -1 ); % Read waveform data from waveform, realwave, wavemarker channels as 32-bit floating point values and returns them as a vector.
        if ch_length <= 0
%             disp([ 'ch_num: ' num2str(ch_num) ' is empty']);
            ch_start_secs   = [];
            ch_values_times = [];
            ch_values_codes = [];
            ch_values       = [];
        else
            ch_start_secs  = CEDS64TicksToSecs( smr_fhand, ch_start_ticks );
            ch_values_times = [];
            ch_values_codes = [];
        end
        
    end
    %% Event
    if ((ch_type >= 2) && (ch_type <= 4))
        [ ch_length, ch_values_ticks ] = CEDS64ReadEvents( smr_fhand, ch_num, ch_max_data_points, 0, -1 );
        if ch_length <= 0
%             disp([ 'ch_num: ' num2str(ch_num) ' is empty']);
            ch_start_secs   = [];
            ch_values_times = [];
            ch_values_codes = [];
            ch_values       = [];
        else
            ch_values_times  = CEDS64TicksToSecs( smr_fhand, ch_values_ticks );
            ch_values = [];
            ch_values_codes = [];
            ch_start_secs   = [];
        end
        
    end
    %% WaveMark
    if ((ch_type >= 6) && (ch_type <= 8))
        [ ch_length, ch_values_ticks ] = CEDS64ReadEvents( smr_fhand, ch_num, ch_max_data_points, 0, -1 );
        if ch_length <= 0
%             disp([ 'ch_num: ' num2str(ch_num) ' is empty']);
            ch_start_secs   = [];
            ch_values_times = [];
            ch_values_codes = [];
            ch_values       = [];
        else
            ch_values_times  = CEDS64TicksToSecs( smr_fhand, ch_values_ticks );
            [ ~, ch_values_CEDWaveMark ] = CEDS64ReadExtMarks( smr_fhand, ch_num, ch_length, 0, -1 );
            ch_values_codes = uint8(zeros(ch_length, 4));
            ch_values       = (zeros(ch_length, length(ch_values_CEDWaveMark(1, 1).m_Data)));
            for event_num = 1 : 1 : ch_length
                ch_values(event_num, :)       = double(ch_values_CEDWaveMark(event_num, 1).m_Data)/(2^15)*10;
                ch_values_codes(event_num, 1) = ch_values_CEDWaveMark(event_num, 1).m_Code1;
                ch_values_codes(event_num, 2) = ch_values_CEDWaveMark(event_num, 1).m_Code2;
                ch_values_codes(event_num, 3) = ch_values_CEDWaveMark(event_num, 1).m_Code3;
                ch_values_codes(event_num, 4) = ch_values_CEDWaveMark(event_num, 1).m_Code4;
            end
            ch_start_secs   = [];
        end
        
    end
    %% Marker
    if (ch_type == 5)
        [ ch_length, ch_values_ticks ] = CEDS64ReadEvents( smr_fhand, ch_num, ch_max_data_points, 0, -1 );
        if ch_length <= 0
%             disp([ 'ch_num: ' num2str(ch_num) ' is empty']);
            ch_start_secs   = [];
            ch_values_times = [];
            ch_values_codes = [];
            ch_values       = [];
        else
            ch_values_times  = CEDS64TicksToSecs( smr_fhand, ch_values_ticks );
            [ ~, ch_values_CEDWaveMark ] = CEDS64ReadMarkers( smr_fhand, ch_num, ch_length, 0, -1 );
            ch_values_codes = uint8(zeros(ch_length, 4));
            for event_num = 1 : 1 : ch_length
                ch_values_codes(event_num, 1) = ch_values_CEDWaveMark(event_num, 1).m_Code1;
                ch_values_codes(event_num, 2) = ch_values_CEDWaveMark(event_num, 1).m_Code2;
                ch_values_codes(event_num, 3) = ch_values_CEDWaveMark(event_num, 1).m_Code3;
                ch_values_codes(event_num, 4) = ch_values_CEDWaveMark(event_num, 1).m_Code4;
            end
            ch_start_secs   = [];
            ch_values       = [];
        end
        
    end
    %% save channel data
    ch_data_struct.ch_number        = ch_num;
    ch_data_struct.type             = ch_type;
    ch_data_struct.type_description = SMR_FILE.ch_type_description{ch_type+1};
    ch_data_struct.title            = ch_title;
    ch_data_struct.unit             = ch_unit;
    ch_data_struct.comment          = ch_comment;
    ch_data_struct.display_range    = [ch_Y_lower, ch_Y_upper];
    ch_data_struct.rate             = ch_rate;
    ch_data_struct.resolution       = ch_resolution_secs;
    ch_data_struct.interval         = ch_interval_secs;
    ch_data_struct.max_time         = ch_max_time_secs;
    ch_data_struct.offset           = ch_offset;
    ch_data_struct.scale            = ch_scale;
    ch_data_struct.start            = ch_start_secs;
    ch_data_struct.values_codes     = ch_values_codes;
    ch_data_struct.values_times     = double(ch_values_times);
    ch_data_struct.values           = double(ch_values);
    
    CH_DATA_ALL(ch_counter, 1) = ch_data_struct;
    ch_counter = ch_counter + 1;
end % end channel_num for loop
end

function save_mat_file(SMR_FILE, CH_DATA_ALL)
%% save .mat file
[~, file_name, ~] = fileparts(SMR_FILE.file_name);
folder_path = SMR_FILE.folder_path;
save([folder_path file_name '.mat'], 'CH_DATA_ALL','SMR_FILE','-v7.3');
end

function close_smr_file(SMR_FILE)
%% close the smr file
[ iOk ] = CEDS64Close( SMR_FILE.smr_fhand );
if (iOk < 0)
    unloadlibrary ceds64int;
    error(['ESN_smr2mat: Cannot close the smr file: ' SMR_FILE.full_path]);
end
end

function close_all()
%% Close all and tidy up
% Finish off by closing both files and unload the DLL.
CEDS64CloseAll(); % close all the files
unloadlibrary ceds64int; % unload ceds64int.dll
end
