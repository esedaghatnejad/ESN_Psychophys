function sac_analyze_monkey_x_axis_all_sessions(SESSION_DATA_ALL)
% Author: Ehsan Sedaghat-Nejad (esedaghatnejad@gmail.com)
% This function is meant for analyzing all the session data of x-axis task in monkeys


%% Build Data
clearvars -except SESSION_DATA_ALL
field_names_DATA_ALL = fieldnames(SESSION_DATA_ALL);
num_sessions = length(field_names_DATA_ALL);
sac_prim_py_r_dir_all = nan(num_sessions, 1100);
sac_prim_py_l_dir_all = nan(num_sessions, 1100);
sac_prim_py_u_dir_all = nan(num_sessions, 1100);
sac_prim_py_d_dir_all = nan(num_sessions, 1100);

trials_num_r_dir_all  = nan(num_sessions, 1100);
trials_num_l_dir_all  = nan(num_sessions, 1100);
trials_num_u_dir_all  = nan(num_sessions, 1100);
trials_num_d_dir_all  = nan(num_sessions, 1100);

sac_prim_py_r_dir_binned_all = nan(num_sessions, 110);
sac_prim_py_l_dir_binned_all = nan(num_sessions, 110);
sac_prim_py_u_dir_binned_all = nan(num_sessions, 110);
sac_prim_py_d_dir_binned_all = nan(num_sessions, 110);

trials_num_r_dir_binned_all  = nan(num_sessions, 110);
trials_num_l_dir_binned_all  = nan(num_sessions, 110);
trials_num_u_dir_binned_all  = nan(num_sessions, 110);
trials_num_d_dir_binned_all  = nan(num_sessions, 110);

for counter_fields = 1 : 1 : num_sessions
    sac_prim_py_r_dir_all(counter_fields, :) = ...
        SESSION_DATA_ALL.(field_names_DATA_ALL{counter_fields}).LEARNING.sac_prim_py_r_dir_;
    sac_prim_py_l_dir_all(counter_fields, :) = ...
        SESSION_DATA_ALL.(field_names_DATA_ALL{counter_fields}).LEARNING.sac_prim_py_l_dir_;
    sac_prim_py_u_dir_all(counter_fields, :) = ...
        SESSION_DATA_ALL.(field_names_DATA_ALL{counter_fields}).LEARNING.sac_prim_py_u_dir_;
    sac_prim_py_d_dir_all(counter_fields, :) = ...
        SESSION_DATA_ALL.(field_names_DATA_ALL{counter_fields}).LEARNING.sac_prim_py_d_dir_;
    
    sac_prim_py_r_dir_binned_all(counter_fields, :) = ...
        SESSION_DATA_ALL.(field_names_DATA_ALL{counter_fields}).LEARNING.sac_prim_py_r_dir_binned_;
    sac_prim_py_l_dir_binned_all(counter_fields, :) = ...
        SESSION_DATA_ALL.(field_names_DATA_ALL{counter_fields}).LEARNING.sac_prim_py_l_dir_binned_;
    sac_prim_py_u_dir_binned_all(counter_fields, :) = ...
        SESSION_DATA_ALL.(field_names_DATA_ALL{counter_fields}).LEARNING.sac_prim_py_u_dir_binned_;
    sac_prim_py_d_dir_binned_all(counter_fields, :) = ...
        SESSION_DATA_ALL.(field_names_DATA_ALL{counter_fields}).LEARNING.sac_prim_py_d_dir_binned_;
    
    trials_num_r_dir_all(counter_fields, :) = ...
        SESSION_DATA_ALL.(field_names_DATA_ALL{counter_fields}).LEARNING.trials_num_r_dir_;
    trials_num_l_dir_all(counter_fields, :) = ...
        SESSION_DATA_ALL.(field_names_DATA_ALL{counter_fields}).LEARNING.trials_num_l_dir_;
    trials_num_u_dir_all(counter_fields, :) = ...
        SESSION_DATA_ALL.(field_names_DATA_ALL{counter_fields}).LEARNING.trials_num_u_dir_;
    trials_num_d_dir_all(counter_fields, :) = ...
        SESSION_DATA_ALL.(field_names_DATA_ALL{counter_fields}).LEARNING.trials_num_d_dir_;
    
    trials_num_r_dir_binned_all(counter_fields, :) = ...
        SESSION_DATA_ALL.(field_names_DATA_ALL{counter_fields}).LEARNING.trials_num_r_dir_binned_;
    trials_num_l_dir_binned_all(counter_fields, :) = ...
        SESSION_DATA_ALL.(field_names_DATA_ALL{counter_fields}).LEARNING.trials_num_l_dir_binned_;
    trials_num_u_dir_binned_all(counter_fields, :) = ...
        SESSION_DATA_ALL.(field_names_DATA_ALL{counter_fields}).LEARNING.trials_num_u_dir_binned_;
    trials_num_d_dir_binned_all(counter_fields, :) = ...
        SESSION_DATA_ALL.(field_names_DATA_ALL{counter_fields}).LEARNING.trials_num_d_dir_binned_;
    
end

%% Plot
clf(figure(1))
hold on
plot(nanmean(trials_num_u_dir_binned_all), nanmean(sac_prim_py_u_dir_binned_all), 'o-');
plot(nanmean(trials_num_d_dir_binned_all), nanmean(sac_prim_py_d_dir_binned_all), 'o-');

end


