%% Global variables
global tag_name_list data_type_eye_list data_type_BEHAVE_list data_type_neuro_list data_type_EPHYS_list event_type_list ...
    waveform_inds_span length_trace inds_span ...
    ang_step ang_edges ang_values amp_edges vel_edges ...
    range_cell_with_4dir_behave
tag_name_list = { ...
    'prim_success', ... % tag 1
    'prim_attempt', ... % tag 2
    'prim_fail', ... % tag 3
    'corr_success', ... % tag 4
    'corr_fail', ... % tag 5
    'back_center_success', ... % tag 6
    'back_center_prim', ... % tag 7
    'back_center_irrelev', ... % tag 8
    'target_irrelev', ... % tag 9
    'other_irrelev', ... % tag 10
    'prim_no_corr',... % tag 11; prim. sac. that is not followed by corr. sac.
    'db_corr_success',... % tag 12; sac. that follows first corr. sac., back to 2nd jumped cue
    'corr_no_db_corr',... % tag 13; corr. sac. that is not followed by another corr. sac.
    };
data_type_eye_list    = {'eye_vx', 'eye_vy'};
data_type_BEHAVE_list = {'eye_r_vx_filt', 'eye_r_vy_filt'};
data_type_neuro_list  = {'neuro_SS', 'neuro_CS'};
data_type_EPHYS_list  = {'EPHYS_SS_train_1K', 'EPHYS_CS_train_1K'};
event_type_list       = {'visual', 'onset', 'vmax', 'offset', 'auditory'};
waveform_inds_span = ((-60+1) : 1 : (120));
length_trace = 500;
inds_span    = ((-(length_trace/2)+1) : 1 : (length_trace/2))';
ang_step     = 45;
ang_edges    = (0 - (ang_step/2)) : ang_step : (360 + (ang_step/2));
ang_values   = (0) : ang_step : (360 - ang_step);
amp_edges    = [0 1.5 2.5 3.5 4.5 5.5 7.5 100];
vel_edges    = [0 150 250 350 450 550 650 750 10000];
range_cell_with_4dir_behave = [-1 -1]; % This is to correct the data for the first round of recordings from Mirza, pCell_list_Mirza_pre201906
% plese set the range_cell_with_4dir_behave to [-1 -1] if you do not have 4dir sessions
%%
flag_pair_list = false; % This should be false. DO NOT change it to true
pCell_ids = build_pCell_ids();
cell_to_analyze = [1:length(pCell_list]; % array of cells to analyze; e.g., [1,5] analyzes 1st and 5th cells in your list

% Figure parameters
Line_Color = lines(7);
color_SS = Line_Color(1,:);
color_offset = Line_Color(2,:);
color_visual = Line_Color(3,:);
color_onset = Line_Color(4,:);
color_vm = Line_Color(5,:);
color_SS_firing = [0    0.3    0.5];
color_CS = Line_Color(7,:);
subplot_list = [6,3,2,1,4,7,8,9]; % in 3x3 plot, to plot from 0 deg counter-clockwise
range_SS_Firing = [0 200];
range_CS_Firing = [0,2];
range_vm        = [0 600];
data_type_list = {'SS', 'CS', 'VT', 'VM'};
CSYS_type_list = {'tuned', 'absol'};
% Which alignment to plot, depending on the tag
which_event_for_tag = cell(length(tag_name_list),length(event_type_list));
which_event_for_tag(1,1:2) = {'visual','onset'}; % tag 1
which_event_for_tag(2,1) = {'onset'}; % tag 2
which_event_for_tag(3,1) = {'onset'}; % tag 3
which_event_for_tag(4,1:2) = {'visual','onset'}; % tag 4
which_event_for_tag(5,1) = {'onset'}; % tag 5
which_event_for_tag(6,1) = {'onset'}; % tag 6
which_event_for_tag(7,1) = {'onset'}; % tag 7
which_event_for_tag(8,1) = {'onset'}; % tag 8
which_event_for_tag(9,1) = {'onset'}; % tag 9
which_event_for_tag(10,1) = {'onset'}; % tag 10
which_event_for_tag(11,1) = {'visual'}; % tag 11
which_event_for_tag(12,1:2) = {'visual','onset'}; % tag 12
which_event_for_tag(13,1) = {'visual'}; % tag 13
% Which event for plot for combination of tags 4+12 and tags 11+13
which_event_for_tag(14,1:2) = {'visual','onset'};
which_event_for_tag(15,1) = {'visual'};
which_event_for_tag_isstr = arrayfun(@iscellstr,which_event_for_tag);


path_cell_data = uigetdir;
if ~strcmp(path_cell_data(end), filesep)
    path_cell_data = [path_cell_data filesep];
end
% Whether to plot fig. in CS-tuned dir. or not
flag_build_absol = true;
flag_build_tuned = false;
if ~xor(flag_build_absol, flag_build_tuned)
    fprintf('ERROR, build_population_data, Please select either flag_build_absol or flag_build_tuned. Not both. The code will run much faster this way.\n')
    return;
end
if flag_build_absol
    CSYS_type = CSYS_type_list{2};
else
    CSYS_type = CSYS_type_list{1};
end
% Loop over pCells
clearvars SACS_ALL_DATA CS_on_data Neural_Properties EXPERIMENT_PARAMS id
for counter_pCell = cell_to_analyze
    fprintf(['### ' 'Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ###' '\n']);
    % load SACS_ALL_DATA
    cell_file_name = pCell_ids{counter_pCell, 1};
    load(fullfile(path_cell_data,all_pcell_folder_name,'cell_data',[cell_file_name,'.mat']));
    % build plot_data address
    year_ = cell_file_name(1:2);
    month_ = cell_file_name(3:4);
    day_ = cell_file_name(5:6);
    hour_ = cell_file_name(8:9);
    minute_ = cell_file_name(10:11);
    second_ = cell_file_name(12:13);
    subFolder_month = ['20' year_ '-' month_ filesep];
    subFolder_day = ['20' year_ '-' month_ '-' day_ filesep];
    subFolder_recording = ['20' year_ '-' month_ '-' day_ '_' hour_ '-' minute_ '-' second_ filesep];
    tag_bin = unique(SACS_ALL_DATA.tag);
    num_tag_bin = length(tag_bin);
    % Combine some tags if available
    % Tag 4 + 12; condition where corr. sac follows the previous sac.
    % See if tag 12 present
    if sum(tag_bin == 12) > 0
        num_tag_bin = num_tag_bin + 1;
    end
    % Tag 11 + 13; condition where corr. sac. could follow but doesn't
    % See if tag 13 present
    if sum(tag_bin == 13) > 0
        num_tag_bin = num_tag_bin + 1;
    end
    num_ang_bin = length(CS_on_data.idx_CS_tuned);
    % Loop over diff. tag and combinations of tags
    for counter_tag = 1:num_tag_bin
        % Combine tags if needed
        if counter_tag == length(tag_name_list) + 1
            idx_tag = (SACS_ALL_DATA.tag == 4) | (SACS_ALL_DATA.tag == 12);
        elseif counter_tag == length(tag_name_list) + 2
            idx_tag = (SACS_ALL_DATA.tag == 11) | (SACS_ALL_DATA.tag == 13);
        else
            idx_tag = (SACS_ALL_DATA.tag == tag_bin(counter_tag));
        end
        % Loop over diff. sac. event
        % Pick which alignment to use to plot
        num_event = sum(which_event_for_tag_isstr(counter_tag,:));
        for counter_event_type = 1 : num_event 
            event_type_name = which_event_for_tag{counter_tag,counter_event_type};
            fig_.('handle_1') = figure(1); clf; % SS firing
%             fig_.('handle_2') = figure(2); clf; % CS firing
%             fig_.('handle_3') = figure(3); clf; % VM
%             fig_.('handle_4') = figure(4); clf; %VT
            % Compute eye speed and tangential speed
            eye_vx = SACS_ALL_DATA.(['eye_vx' '_' event_type_name]);
            eye_vy = SACS_ALL_DATA.(['eye_vy' '_' event_type_name]);
            eye_vm = sqrt(eye_vx.^2 + eye_vy.^2 );
            eye_amp_x = repmat(SACS_ALL_DATA.eye_r_amp_x, length_trace, 1);
            eye_amp_y = repmat(SACS_ALL_DATA.eye_r_amp_y, length_trace, 1);
            eye_amp_m = repmat(SACS_ALL_DATA.eye_r_amp_m, length_trace, 1);
            eye_vt = ( (eye_vx.*eye_amp_x) + (eye_vy.*eye_amp_y) ) ./ eye_amp_m;
            % Loop over diff. angles
            for counter_ang = 1 : num_ang_bin
                if flag_build_absol
                    idx_ang = (CS_on_data.visual_ang_bin == counter_ang);
                elseif flag_build_tuned
                    idx_ang = (CS_on_data.visual_ang_bin == CS_on_data.idx_CS_tuned(counter_ang));
                end
                sac_idx = idx_tag & idx_ang;
                num_sac = sum(sac_idx);
                % Loop over figures; all plots have SS, CS, visual, sac. onset and offset rastor plots underlayed. 
                % Eye speed plotted by dividing it by 3
                train_data_logic_SS_ = SACS_ALL_DATA.(['neuro_SS' '_' event_type_name])(:,sac_idx)';
                train_data_logic_CS_ = SACS_ALL_DATA.(['neuro_CS' '_' event_type_name])(:,sac_idx)';
                for counter_fig = 1:1 % temporarily only plot SS
                    figure(fig_.(sprintf('handle_%d',counter_fig))); 
                    subplot(3,3,subplot_list(counter_ang)); hold on;
                    % SS rastor
                    [x_axis_SS_, y_axis_SS_] = ESN_raster_plot_axes(train_data_logic_SS_,inds_span,0.5);
                    plot(x_axis_SS_(:), y_axis_SS_(:), 'LineWidth', 1, 'Color', color_SS);
                    % CS rastor
                    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_CS_, inds_span, 1);
                    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 3, 'Color', color_CS);
                    % Visual rastor
                    time_visual_shifted = (SACS_ALL_DATA.time_visual(sac_idx) - SACS_ALL_DATA.(sprintf('time_%s',event_type_name))(sac_idx)).*1000;
                    time_visual_shifted = ESN_Round(time_visual_shifted);
                    idx_visual_shifted = time_visual_shifted + length_trace/2;
                    train_data_logic_visual = repmat(1:length_trace,num_sac,1) == idx_visual_shifted';
                    [x_axis_visual_, y_axis_visual_] = ESN_raster_plot_axes(train_data_logic_visual, inds_span, 1);
                    plot(x_axis_visual_(:), y_axis_visual_(:), 'LineWidth', 2, 'Color', color_visual);
                    % Onset rastor
                    time_onset_shifted = (SACS_ALL_DATA.time_onset(sac_idx) - SACS_ALL_DATA.(sprintf('time_%s',event_type_name))(sac_idx)).*1000;
                    time_onset_shifted = ESN_Round(time_onset_shifted);
                    idx_onset_shifted = time_onset_shifted + length_trace/2;
                    train_data_logic_onset = repmat(1:length_trace,num_sac,1) == idx_onset_shifted';
                    [x_axis_onset_, y_axis_onset_] = ESN_raster_plot_axes(train_data_logic_onset, inds_span, 1);
                    plot(x_axis_onset_(:), y_axis_onset_(:), 'LineWidth', 2, 'Color', color_onset);
                    % Offset rastor
                    time_offset_shifted = (SACS_ALL_DATA.time_offset(sac_idx) - SACS_ALL_DATA.(sprintf('time_%s',event_type_name))(sac_idx)).*1000;
                    time_offset_shifted = ESN_Round(time_offset_shifted);
                    idx_offset_shifted = time_offset_shifted + length_trace/2;
                    train_data_logic_offset = repmat(1:length_trace,num_sac,1) == idx_offset_shifted';
                    [x_axis_offset_, y_axis_offset_] = ESN_raster_plot_axes(train_data_logic_offset, inds_span, 1);
                    plot(x_axis_offset_(:), y_axis_offset_(:), 'LineWidth', 2, 'Color', color_offset);
                    
                    
                    
                    xlim([min(inds_span)-1 max(inds_span)+1]);
                    ylim([(1-3) (size(train_data_logic_CS_,1)+3)]);
                    if ismember(counter_ang, [4,5,6])
                        ylabel('Trials');
                    end
                    if ismember(counter_ang, 7)
                        xlabel(sprintf('Time relative to %s (ms)',event_type_name));
                    end
                    yyaxis right;
                    xlim([min(inds_span)-1 max(inds_span)+1]);
                    % Eye speed
                    plot(inds_span, mean(eye_vm(:,sac_idx),2,'omitnan')./3,'-','LineWidth',1.5,'Color',color_vm);
                    % SS
                    if counter_fig == 1
                        firing_SS_ = mean(train_data_logic_SS_,'omitnan') * 1000;
                        plot(inds_span, ESN_smooth(firing_SS_),'-', 'LineWidth', 2, 'Color', color_SS_firing);
                        ylim(range_SS_Firing);
                        set(gca, 'YColor', color_SS_firing);
                        if ismember(counter_ang,[1,2,8])
                            ylabel('SS Firing (spk/s)');
                        end
                    % CS
                    elseif counter_fig == 2
                        firing_CS_ = mean(train_data_logic_CS_, 'omitnan') * 1000;
                        plot(inds_span, ESN_smooth(firing_CS_), 'LineWidth', 2, 'Color', color_CS);
                        ylim(range_CS_Firing);
                        set(gca, 'YColor', color_CS);
                        if ismember(counter_ang, [1,2,8])
                            ylabel('CS Firing (spk/s)');
                        end
                    % Vm
                    elseif counter_fig == 3
                        plot(inds_span, mean(eye_vm(:,sac_idx)','omitnan'),'LineWidth',2,'Color', color_vm);
                        ylim(range_vm);
                        set(gca, 'YColor', color_vm);
                        if ismember(counter_ang, [1,2,8])
                            ylabel('Speed (deg/s)');
                        end
                    % Vt
                    elseif counter_fig == 4
                        plot(inds_span, mean(eye_vt(:,sac_idx)','omitnan'),'LineWidth',2,'Color', color_vm);
                        ylim(range_vm);
                        set(gca, 'YColor', color_vm);
                        if ismember(counter_ang, [1,2,8])
                            ylabel({'Tangential','speed (deg/s)'});
                        end
                    end
%                     % Eye speed
%                     plot(inds_span, mean(eye_vm(:,sac_idx),2,'omitnan')./3,'-','LineWidth',1.5,'Color',color_vm);
                    
                    if counter_ang == num_ang_bin
                        % CS probability
                        subplot(3,3,5); hold on;
                        plot(cosd(0:5:360)*0.10, sind(0:5:360)*0.10, 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);
                        plot(cosd(0:5:360)*0.20, sind(0:5:360)*0.20, 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);
                        plot(cosd(0:5:360)*0.30, sind(0:5:360)*0.30, 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);
                        % Tag 4 + 12
                        if counter_tag == length(tag_name_list) + 1
                            prob_amplitude = (CS_on_data.CS_count(4,:) + CS_on_data.CS_count(12,:))./ ...
                                (CS_on_data.sac_count(4,:) + CS_on_data.sac_count(12,:));
                        % Tag 11 + 13    
                        elseif counter_tag == length(tag_name_list) + 2
                            prob_amplitude = (CS_on_data.CS_count(11,:) + CS_on_data.CS_count(13,:))./ ...
                                (CS_on_data.sac_count(11,:) + CS_on_data.sac_count(13,:));
                        else
                            if flag_build_absol
                                prob_amplitude = CS_on_data.CS_prob(counter_tag,:);
                            elseif flag_build_tuned
                                prob_amplitude = CS_on_data.CS_prob(counter_tag,CS_on_data.idx_CS_tuned);
                            end
                        end
                        prob_amplitude(end+1) = prob_amplitude(1);
                        x_axis = [cosd(0) cosd(45) cosd(90) cosd(135) cosd(180) cosd(225) cosd(270) cosd(315) cosd(0)] .* prob_amplitude;
                        y_axis = [sind(0) sind(45) sind(90) sind(135) sind(180) sind(225) sind(270) sind(315) sind(0)] .* prob_amplitude;
                        x_axis = x_axis(~isnan(x_axis)); y_axis = y_axis(~isnan(y_axis));
                        plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', color_CS);

                        CS_prob_avg = prob_amplitude(1:end-1);
                        dir_angles = deg2rad([0, 45, 90, 135, 180, 225, 270, 315]);
                        r = sum(CS_prob_avg.* exp(1i*dir_angles) , 2,'omitnan'); % compute weighted sum of cos and sin of angles
                        CS_ang_avg = rad2deg(angle(r)); % Computes the mean direction for circular data.
                        CS_prob_sum = sum(CS_prob_avg,2,'omitnan'); % sum of weights
                        CS_rho_avg = abs(r) ./ CS_prob_sum; % Computes mean resultant vector length for circular data.
                        plot([0 CS_rho_avg*cosd(CS_ang_avg)], [0 CS_rho_avg*sind(CS_ang_avg)], 'LineWidth', 3, 'Color', color_CS);

                        axis equal
                        set(gca, 'YColor', color_CS)
                        set(gca, 'XColor', color_CS)
                        xlabel({'CS probability based on [0,200] ms','from visual event'});
                        % Title
                        figure_title(1) = {cell_file_name};
                        if counter_tag == length(tag_name_list) + 1
                            figure_title(2) = {['corr_follow, tag 4 + 12, ' data_type_list{counter_fig}, ', ' CSYS_type, ', ' event_type_name]};
                        elseif counter_tag == length(tag_name_list) + 2
                            figure_title(2) = {['no_corr_follow, tag 11 + 13, ' data_type_list{counter_fig}, ', ' CSYS_type, ', ' event_type_name]};
                        else
                            figure_title(2) = {[tag_name_list{counter_tag} ', ' ...
                                                data_type_list{counter_fig}, ', ' ...
                                                CSYS_type, ', ' ...
                                                event_type_name]};
                        end
                        sgtitle(figure_title,'Interpreter','none');                          

                        ESN_Beautify_Plot(fig_.(sprintf('handle_%d',counter_fig)), [8.0 8.0])
                    end
                    
                end
            end
            % Save figure             
%             path_fig_ = fullfile(path_cell_data,subFolder_month,subFolder_day,subFolder_recording,['bundle_figs_' cell_file_name],...
%                 CSYS_type,data_type_list{counter_fig},[num2str(counter_tag) '_' tag_name_list{counter_tag}]);
            path_fig_ = fullfile(path_cell_data,subFolder_month,subFolder_day,subFolder_recording,['bundle_figs_' cell_file_name]);
            if ~exist(path_fig_,'dir')
                %multiple_file_dirs = dir(fullfile(path_cell_data, 'ALL_PCELL_33\cell_data\',[cell_file_name(1:6),'*']));
                mkdir(path_fig_);
            end
            if counter_tag == length(tag_name_list) + 1
                file_name_fig_ =  [CSYS_type '_' data_type_list{counter_fig} '_4_12_corr_follow_' event_type_name];
            elseif counter_tag == length(tag_name_list) + 2
                file_name_fig_ =  [CSYS_type '_' data_type_list{counter_fig} '_11_13_no_corr_follow_' event_type_name];
            else
                file_name_fig_ =  [CSYS_type '_' data_type_list{counter_fig} '_' num2str(counter_tag) '_' tag_name_list{counter_tag} '_' event_type_name];
            end
            saveas(figure(fig_.(sprintf('handle_%d',counter_fig))), fullfile(path_fig_,file_name_fig_),'pdf');
        end
    end

end
%% function build_pCell_ids()
function pCell_ids = build_pCell_ids()
global flag_pair_list
pCell_list = ESN_build_pCell_list(flag_pair_list);
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);
pCell_ids = cell(num_pCells, 1);
for counter_pCell = 1 : num_pCells
    file_name_cell = pCell_list{counter_pCell, 1};
    if file_name_cell(18) == 's'
        id_          = file_name_cell(1:16);
    elseif file_name_cell(18) == '2'
        id_          = file_name_cell(1:18);
    else
        error('Build plot_data_compress: cell id does not follow the standards')
    end
    num_recording = nansum(pCell_list_isstr(counter_pCell, :));
    cell_file_name = [id_ '_' 'combine' '_' num2str(num_recording)];
    pCell_ids{counter_pCell, 1} = cell_file_name;
end
end