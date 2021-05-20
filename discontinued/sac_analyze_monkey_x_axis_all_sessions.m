function sac_analyze_monkey_x_axis_all_sessions(SESSION_DATA_ALL)
% Author: Ehsan Sedaghat-Nejad (esedaghatnejad@gmail.com)
% This function is meant for analyzing all the session data of x-axis task in monkeys

%% Remove sessions
rmfields_list = {...
    'SESSION_DATA_2018_08_10',...
    'SESSION_DATA_2018_08_27',...
    'SESSION_DATA_2018_09_06',...
    };
SESSION_DATA_ALL = rmfield(SESSION_DATA_ALL,rmfields_list);

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

%% Plot LEARNING
clf(figure(1))
hold on
plot_mean_std_color_vec = lines(7);
plot_mean_std_params.marker_on_off = 0;
plot_mean_std_params.std_boundry_on_off = 0;
plot_mean_std_params.linewidth = 2;
plot_mean_std_params.MarkerSize = 4;
plot_mean_std_params.Std_Alpha = 0.35;
plot_mean_std_params.Line_Alpha = 0.8;

plot([0 1100],  [0 0], '-k')
plot([100 100], [-0.31 0.31], '-k')
plot([800 800], [-0.31 0.31], '-k')

% plot(trials_num_u_dir_binned_all', sac_prim_py_u_dir_binned_all')
% plot(trials_num_d_dir_binned_all', sac_prim_py_d_dir_binned_all')

X_Values = nanmean(trials_num_u_dir_binned_all);
Y_Mean = nanmean(sac_prim_py_u_dir_binned_all);
Y_Std = nanstd(sac_prim_py_u_dir_binned_all) ./ sqrt(sum(~isnan(sac_prim_py_u_dir_binned_all)));
ax_1_ = ESN_Plot_MeanStd(X_Values, Y_Mean, Y_Std, plot_mean_std_color_vec(1,:), plot_mean_std_params);



X_Values = nanmean(trials_num_d_dir_binned_all);
Y_Mean   = nanmean(sac_prim_py_d_dir_binned_all);
Y_Std    = nanstd(sac_prim_py_d_dir_binned_all) ./ sqrt(sum(~isnan(sac_prim_py_d_dir_binned_all)));
ax_2_ = ESN_Plot_MeanStd(X_Values, Y_Mean, Y_Std, plot_mean_std_color_vec(2,:), plot_mean_std_params);

legend([ax_1_ ax_2_], {'X-Up', 'X-Down'})
ylim([-0.31 0.31])
xlim([0 1100])
xlabel('Trial number (#)')
ylabel('Sac Vert End Position (deg)')
title(['Num Session = ' num2str(num_sessions)])
ESN_Beautify_Plot

%% Plot Session-by-Session laerning
clf(figure(2))
hold on
num_sessions = length(field_names_DATA_ALL);
for counter_sessions = 1 : num_sessions
    if ( sum( ~isnan(sac_prim_py_r_dir_all(counter_sessions, :)) & ~isnan(sac_prim_py_u_dir_all(counter_sessions, :)) ) > ...
         sum( ~isnan(sac_prim_py_r_dir_all(counter_sessions, :)) & ~isnan(sac_prim_py_d_dir_all(counter_sessions, :)) ) )
        Label = ['090, ' (field_names_DATA_ALL{counter_sessions})];
        
    else
        Label = ['270, ' (field_names_DATA_ALL{counter_sessions})];
    end
    
    subplot(6,5, counter_sessions)
    hold on
    plot([0 1100], [0 0], '-k')
    plot(trials_num_u_dir_binned_all(counter_sessions, :), sac_prim_py_u_dir_binned_all(counter_sessions, :))
    plot(trials_num_d_dir_binned_all(counter_sessions, :), sac_prim_py_d_dir_binned_all(counter_sessions, :))
    ylim([-0.7 0.7])
    xlim([0 1100])
    title(Label, 'interpret', 'none')
end

end


