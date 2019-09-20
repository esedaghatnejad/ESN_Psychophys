function ESN_session_to_mat(session_folder_path)
% get the session folder path and convert all the *.continuous files for that session to mat files

% if there is no inputs, then set folder_path to pwd
if nargin < 1
    session_folder_path = pwd;
end
% add filesep ('/' or '\') to the end of folder_path
if ~strcmp(session_folder_path(end), filesep)
    session_folder_path = [session_folder_path filesep];
end

% session_folder_list = genpath(session_folder_path);
session_folder_list = dir(session_folder_path);

for counter_list = 3 : length(session_folder_list)
    if ~session_folder_list(counter_list).isdir
        continue;
    end
    data_folder_path = [session_folder_path session_folder_list(counter_list).name];
    continuous_file_list = dir([data_folder_path filesep '*.continuous']);
%     include_mat_file_folder = dir([data_folder_path filesep '*mat_files*']);
%     if (length(continuous_file_list) > 4) && (length(include_mat_file_folder)<1)
    if (length(continuous_file_list) > 4)
        try
            disp(data_folder_path);
            ESN_continuous_to_mat(data_folder_path);
            ESN_heptode_combine([data_folder_path filesep 'mat_files']);
        catch err_msg
            %open file
            fid = fopen([data_folder_path filesep 'errorLog.txt'],'a+');
            % write the error to file
            fprintf(fid,'%s\n',err_msg.message);
            % close file
            fclose(fid);
        end
    end
end




