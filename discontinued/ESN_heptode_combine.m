function ESN_heptode_combine(folder_path)
% Combine 4 channels together and delete unneccessary files

% if there is no inputs, then set folder_path to pwd
if nargin < 1
    folder_path = pwd;
end
% add filesep ('/' or '\') to the end of folder_path
if ~strcmp(folder_path(end), filesep)
    folder_path = [folder_path filesep];
end

FILE_NAMES = dir([folder_path '*.mat']);
port_num_string = FILE_NAMES(1).name(1:3);
% get the list of 100_CH*.mat files to analyze
FILE_NAMES = dir([folder_path port_num_string '_CH*.mat']);
% extract the file names
mat_file_list = extractfield(FILE_NAMES, 'name');
% extract the numbers
filenum = cellfun(@(x)sscanf(x,[port_num_string '_CH%d.mat']), mat_file_list);
% sort the numbers, and get the sorting order
[~,Sidx] = sort(filenum);
% use to this sorting order to sort the filenames
mat_file_list = mat_file_list(Sidx);
num_channels = length(mat_file_list);
counter_CH_CMN = 0;
for counter_heptode_channel = 1 : 1 : 7
    fprintf([num2str(counter_heptode_channel) '/' '9' '. Combining : ' port_num_string '_CH' num2str(counter_heptode_channel) ' ... '])
    
    CH_DATA_1 = load([folder_path mat_file_list{((counter_heptode_channel-1)*4)+1}]);
    CH_DATA_2 = load([folder_path mat_file_list{((counter_heptode_channel-1)*4)+2}]);
    CH_DATA_3 = load([folder_path mat_file_list{((counter_heptode_channel-1)*4)+3}]);
    CH_DATA_4 = load([folder_path mat_file_list{((counter_heptode_channel-1)*4)+4}]);
    clearvars('ch_data', 'ch_time', 'ch_info');
    ch_data = mean([CH_DATA_1.CH_DATA.ch_data(:), CH_DATA_2.CH_DATA.ch_data(:), CH_DATA_3.CH_DATA.ch_data(:), CH_DATA_4.CH_DATA.ch_data(:)], 2)';
    ch_time = mean([CH_DATA_1.CH_DATA.ch_time(:), CH_DATA_2.CH_DATA.ch_time(:), CH_DATA_3.CH_DATA.ch_time(:), CH_DATA_4.CH_DATA.ch_time(:)], 2)';
    ch_info = CH_DATA_1.CH_DATA.ch_info;
    delete([folder_path mat_file_list{((counter_heptode_channel-1)*4)+1}]);
    delete([folder_path mat_file_list{((counter_heptode_channel-1)*4)+2}]);
    delete([folder_path mat_file_list{((counter_heptode_channel-1)*4)+3}]);
    delete([folder_path mat_file_list{((counter_heptode_channel-1)*4)+4}]);
    save(  [folder_path port_num_string '_CH' num2str(counter_heptode_channel) '.mat'], 'ch_data', 'ch_time', 'ch_info', '-v7.3');
    
    if ~exist('CH_CMN_data','var')
        CH_CMN_data = zeros(size(ch_data));
        CH_CMN_time = ch_time;
        CH_CMN_info = ch_info;
    end
    
    CH_CMN_data = CH_CMN_data + ch_data;
    counter_CH_CMN = counter_CH_CMN + 1;
    
    fprintf(' --> Completed.\n')
end

fprintf(['8' '/' '9' '. Combining : ' port_num_string '_CH' '_CMN1' ' ... '])
CH_CMN_data = CH_CMN_data ./ counter_CH_CMN;
ch_data = CH_CMN_data;
ch_time = CH_CMN_time;
ch_info = CH_CMN_info;
save([folder_path port_num_string '_CMN1.mat'], 'ch_data', 'ch_time', 'ch_info', '-v7.3');
fprintf(' --> Completed.\n')

fprintf(['8' '/' '9' '. Combining : ' port_num_string '_CH' '_CMN1' ' ... '])
CH_EVE_data = load([folder_path port_num_string '_EVE1.mat']);
ch_time_EPHYS = CH_CMN_time;
ch_data = CH_EVE_data.CH_DATA.ch_data;
ch_time = CH_EVE_data.CH_DATA.ch_time;
ch_info = CH_EVE_data.CH_DATA.ch_info;
save([folder_path port_num_string '_EVE1.mat'], 'ch_data', 'ch_time', 'ch_info', 'ch_time_EPHYS', '-v7.3');
fprintf(' --> Completed.\n')

fprintf('Deleting unneccessary files ...')
% delete CH29 -> CH32
delete([folder_path mat_file_list{((8-1)*4)+1}]);
delete([folder_path mat_file_list{((8-1)*4)+2}]);
delete([folder_path mat_file_list{((8-1)*4)+3}]);
delete([folder_path mat_file_list{((8-1)*4)+4}]);

% delete AUX
FILE_NAMES = dir([folder_path port_num_string '_AUX*.mat']);
for counter_file = 1 : length(FILE_NAMES)
    delete([folder_path FILE_NAMES(counter_file).name]);
end
fprintf(' --> Completed.\n')

end