function EXPERIMENT = ESN_Sac_Review(EXPERIMENT, trial_ind)


counter_main        = EXPERIMENT.all.block_num(:,trial_ind);
trial_number        = EXPERIMENT.all.trial_num(:,trial_ind);
inds_trial          = EXPERIMENT.all.inds_trial(:,trial_ind);
ind_stateSacBegin   = EXPERIMENT.all.ind_stateSacBegin_trl(:,trial_ind);
ind_stateSacEnd     = EXPERIMENT.all.ind_stateSacEnd_trl(:,trial_ind);
ind_stateFixateBegin= EXPERIMENT.all.ind_stateFixateBegin_trl(:,trial_ind);
ind_stateFixateEnd  = EXPERIMENT.all.ind_stateFixateEnd_trl(:,trial_ind);
ind_sac_start       = EXPERIMENT.all.sac_ind_start_trl(:,trial_ind);
ind_sac_vmax        = EXPERIMENT.all.sac_ind_vmax_trl(:,trial_ind);
ind_sac_finish      = EXPERIMENT.all.sac_ind_finish_trl(:,trial_ind);
ind_sac2_start      = EXPERIMENT.all.sac2_ind_start_trl(:,trial_ind);
ind_sac2_vmax       = EXPERIMENT.all.sac2_ind_vmax_trl(:,trial_ind);
ind_sac2_finish     = EXPERIMENT.all.sac2_ind_finish_trl(:,trial_ind);
sac_validity        = EXPERIMENT.all.sac_validity(:,trial_ind);
sac2_validity       = EXPERIMENT.all.sac2_validity(:,trial_ind);
eye_vm              = EXPERIMENT.all.eye_vm(:,trial_ind);
eye_px_s            = EXPERIMENT.all.eye_px_s(:,trial_ind);
eye_py_s            = EXPERIMENT.all.eye_py_s(:,trial_ind);
tgt_str_x           = EXPERIMENT.all.tgt_str_x(:,trial_ind);
tgt_cue_x           = EXPERIMENT.all.tgt_cue_x(:,trial_ind);
tgt_end_x           = EXPERIMENT.all.tgt_end_x(:,trial_ind);
perturbation_jump   = 3;


% %% Plot datacursors

inds_sac_1_2          = min([ind_stateSacBegin, ind_sac_start]-5):1:max([(ind_stateFixateEnd) (ind_sac2_finish) (ind_sac_finish)]);
sac_1_2_vm          = double(eye_vm(inds_sac_1_2));
sac_1_2_vm_filt = sac_1_2_vm; % filtfilt(b_butter,a_butter,sac_1_2_vm); % sgolayfilt(sac_vm, 3, 25);

h_fig_ = figure(1);
clf(figure(1));
dcm_obj = datacursormode(h_fig_);
set(dcm_obj,'UpdateFcn',@ESN_DcmUpdateFcn);

subplot(2,1,1)
hold('on');
h_plot_ = plot(inds_sac_1_2, sac_1_2_vm_filt, 'linewidth', 2);
hAxes(1) = gca;
ESN_Add_DataCursor(h_plot_, ind_sac_start   - inds_sac_1_2(1));
ESN_Add_DataCursor(h_plot_, ind_sac_vmax    - inds_sac_1_2(1));
ESN_Add_DataCursor(h_plot_, ind_sac_finish  - inds_sac_1_2(1));
ESN_Add_DataCursor(h_plot_, ind_sac2_start  - inds_sac_1_2(1));
ESN_Add_DataCursor(h_plot_, ind_sac2_vmax   - inds_sac_1_2(1));
ESN_Add_DataCursor(h_plot_, ind_sac2_finish - inds_sac_1_2(1));

area(((ind_sac_start:ind_sac_finish)), sac_1_2_vm_filt(((ind_sac_start:ind_sac_finish)-inds_sac_1_2(1)+1)))
area(((ind_sac2_start:ind_sac2_finish)), sac_1_2_vm_filt(((ind_sac2_start:ind_sac2_finish)-inds_sac_1_2(1)+1)))
plot([inds_sac_1_2(1); inds_sac_1_2(end)], [20; 20],   'r')
plot([inds_sac_1_2(1); inds_sac_1_2(end)], [150; 150], '--r')
plot([inds_sac_1_2(1); inds_sac_1_2(end)], [50; 50],   'r')
y_lims_ = get(gca, 'YLim');
plot([ind_stateSacBegin;    ind_stateSacBegin],   y_lims_(:), 'r')
plot([ind_stateSacEnd;      ind_stateSacEnd],     y_lims_(:), 'r')
plot([ind_stateFixateBegin; ind_stateFixateBegin],y_lims_(:), 'r')
plot([ind_stateFixateEnd;   ind_stateFixateEnd],  y_lims_(:), 'r')
title(['validity (sac1,sac2): ' num2str(sac_validity) num2str(sac2_validity) ...
    ', block: ' num2str(counter_main) ', trial: ' num2str(trial_number)]);

subplot(2,1,2)
hold('on');
sac_1_2_px =  double(eye_px_s(inds_sac_1_2));
sac_1_2_py =  double(eye_py_s(inds_sac_1_2));
plot(inds_sac_1_2, sac_1_2_px, 'linewidth', 2);
area(((ind_sac_start:ind_sac_finish)), sac_1_2_px(((ind_sac_start:ind_sac_finish)-inds_sac_1_2(1)+1)))
area(((ind_sac2_start:ind_sac2_finish)), sac_1_2_px(((ind_sac2_start:ind_sac2_finish)-inds_sac_1_2(1)+1)))

plot(inds_sac_1_2, sac_1_2_py, 'linewidth', 2);
area(((ind_sac_start:ind_sac_finish)), sac_1_2_py(((ind_sac_start:ind_sac_finish)-inds_sac_1_2(1)+1)))
area(((ind_sac2_start:ind_sac2_finish)), sac_1_2_py(((ind_sac2_start:ind_sac2_finish)-inds_sac_1_2(1)+1)))

plot([inds_sac_1_2(1); inds_sac_1_2(end)], [0; 0], 'r')
plot([inds_sac_1_2(1); inds_sac_1_2(end)], repmat((tgt_cue_x - tgt_str_x),2,1), 'r')
plot([inds_sac_1_2(1); inds_sac_1_2(end)], repmat((tgt_end_x - tgt_str_x),2,1), 'r')
plot([inds_sac_1_2(1); inds_sac_1_2(end)], repmat((tgt_cue_x - tgt_str_x),2,1)+perturbation_jump, '--r')
plot([inds_sac_1_2(1); inds_sac_1_2(end)], repmat((tgt_cue_x - tgt_str_x),2,1)-perturbation_jump, '--r')

y_lims_ = get(gca, 'YLim');
plot([ind_stateSacBegin;    ind_stateSacBegin],   y_lims_(:), 'r')
plot([ind_stateSacEnd;      ind_stateSacEnd],     y_lims_(:), 'r')
plot([ind_stateFixateBegin; ind_stateFixateBegin],y_lims_(:), 'r')
plot([ind_stateFixateEnd;   ind_stateFixateEnd],  y_lims_(:), 'r')

title(['validity (sac1,sac2): ' num2str(sac_validity) num2str(sac2_validity) ...
    ', block: ' num2str(counter_main) ', trial: ' num2str(trial_number)]);

hAxes(2) = gca;
linkaxes(hAxes, 'x');

% %% User Input
input_validity = input(['TRIAL: ' num2str(trial_number) ', validity (sac1,sac2): ' num2str(sac_validity) num2str(sac2_validity) ', :']);

c_info = getCursorInfo(dcm_obj);
field_names_ = fieldnames(c_info);
temp_cell_ = struct2cell(c_info);
for counter_fields = 2 : 1 : length(field_names_)
    field_cell_ = temp_cell_(counter_fields,:,:);
    field_cell_ = reshape(field_cell_, 1, []);
    max_numel_ = max(cellfun(@numel,field_cell_));
    field_mat_ = cell2mat(cellfun(@(x) vertcat(double(x(:)),NaN(max_numel_-numel(x), 1)),field_cell_,'uni',false));
    c_info_.(field_names_{counter_fields}) = field_mat_;
end

c_info = c_info_;
[DataIndex, ind_sort_] = sort(c_info.DataIndex);
datacursor_pos_ = c_info.Position(:, ind_sort_);
datacursor_ind_ = c_info.DataIndex(:, ind_sort_);

ind_sac_start   = datacursor_pos_(1,1);
ind_sac_vmax    = datacursor_pos_(1,2);
ind_sac_finish  = datacursor_pos_(1,3);
ind_sac2_start  = datacursor_pos_(1,4);
ind_sac2_vmax   = datacursor_pos_(1,5);
ind_sac2_finish = datacursor_pos_(1,6);
inds_sac        = (ind_sac_vmax -61) : 1 : (ind_sac_vmax -61+149);
inds_sac2       = (ind_sac2_vmax-61) : 1 : (ind_sac2_vmax-61+149);
sac_vmax        = datacursor_pos_(2,2);
sac2_vmax       = datacursor_pos_(2,5);
if(isempty(input_validity));  input_validity = -1; end;
switch input_validity
    case 00
        sac_validity  = false;
        sac2_validity = false;
    case 01
        sac_validity  = false;
        sac2_validity = true;
    case 10
        sac_validity  = true;
        sac2_validity = false;
    case 11
        sac_validity  = true;
        sac2_validity = true;
    otherwise
end


% Make Sure that inds_sac2 does not exceeds the dimension of inds_trial
length_inds_trial = sum(~isnan(inds_trial));
inds_sac2(inds_sac2>length_inds_trial) = length_inds_trial;

% %% Save sac & sac2
EXPERIMENT.all.sac_validity(:,trial_ind)        = sac_validity;
EXPERIMENT.all.sac_inds(:,trial_ind)            = double(inds_sac);
EXPERIMENT.all.sac_ind_start_trl(:,trial_ind)   = double(ind_sac_start);
EXPERIMENT.all.sac_ind_vmax_trl(:,trial_ind)    = double(ind_sac_vmax);
EXPERIMENT.all.sac_ind_finish_trl(:,trial_ind)  = double(ind_sac_finish);

EXPERIMENT.all.sac_px(:,trial_ind)              = double(EXPERIMENT.all.eye_px(inds_sac,trial_ind));
EXPERIMENT.all.sac_py(:,trial_ind)              = double(EXPERIMENT.all.eye_py(inds_sac,trial_ind));
EXPERIMENT.all.sac_vx(:,trial_ind)              = double(EXPERIMENT.all.eye_vx(inds_sac,trial_ind));
EXPERIMENT.all.sac_vy(:,trial_ind)              = double(EXPERIMENT.all.eye_vy(inds_sac,trial_ind));
EXPERIMENT.all.sac_vm(:,trial_ind)              = double(EXPERIMENT.all.eye_vm(inds_sac,trial_ind));
EXPERIMENT.all.sac_vm_max(:,trial_ind)          = double(EXPERIMENT.all.eye_vm(ind_sac_vmax,trial_ind));
EXPERIMENT.all.sac_px_s(:,trial_ind)            = EXPERIMENT.all.sac_px(:,trial_ind) - EXPERIMENT.all.tgt_str_x(:,trial_ind);
EXPERIMENT.all.sac_py_s(:,trial_ind)            = EXPERIMENT.all.sac_py(:,trial_ind) - EXPERIMENT.all.tgt_str_y(:,trial_ind);

EXPERIMENT.all.sac_px_sacStart(:,trial_ind)     = double(EXPERIMENT.all.eye_px(ind_sac_start,trial_ind));
EXPERIMENT.all.sac_px_sacFinish(:,trial_ind)    = double(EXPERIMENT.all.eye_px(ind_sac_finish,trial_ind));
EXPERIMENT.all.sac_py_sacStart(:,trial_ind)     = double(EXPERIMENT.all.eye_py(ind_sac_start,trial_ind));
EXPERIMENT.all.sac_py_sacFinish(:,trial_ind)    = double(EXPERIMENT.all.eye_py(ind_sac_finish,trial_ind));
EXPERIMENT.all.sac_px_sacStart_s(:,trial_ind)   = double(EXPERIMENT.all.eye_px_s(ind_sac_start,trial_ind));
EXPERIMENT.all.sac_px_sacFinish_s(:,trial_ind)  = double(EXPERIMENT.all.eye_px_s(ind_sac_finish,trial_ind));
EXPERIMENT.all.sac_py_sacStart_s(:,trial_ind)   = double(EXPERIMENT.all.eye_py_s(ind_sac_start,trial_ind));
EXPERIMENT.all.sac_py_sacFinish_s(:,trial_ind)  = double(EXPERIMENT.all.eye_py_s(ind_sac_finish,trial_ind));
EXPERIMENT.all.sac_amp_x(:,trial_ind)           = double(EXPERIMENT.all.sac_px_sacFinish(:,trial_ind) - EXPERIMENT.all.sac_px_sacStart(:,trial_ind));
EXPERIMENT.all.sac_amp_y(:,trial_ind)           = double(EXPERIMENT.all.sac_py_sacFinish(:,trial_ind) - EXPERIMENT.all.sac_py_sacStart(:,trial_ind));
EXPERIMENT.all.sac_amp_m(:,trial_ind)           = double(sqrt(EXPERIMENT.all.sac_amp_x(:,trial_ind)^2+EXPERIMENT.all.sac_amp_y(:,trial_ind)^2));
EXPERIMENT.all.sac_reaction(:,trial_ind)        = EXPERIMENT.all.sac_ind_start_trl(:,trial_ind) - EXPERIMENT.all.ind_stateSacBegin_trl(:,trial_ind);

EXPERIMENT.all.sac2_validity(:,trial_ind)       = sac2_validity;
EXPERIMENT.all.sac2_inds(:,trial_ind)           = double(inds_sac2);
EXPERIMENT.all.sac2_ind_start_trl(:,trial_ind)  = double(ind_sac2_start);
EXPERIMENT.all.sac2_ind_vmax_trl(:,trial_ind)   = double(ind_sac2_vmax);
EXPERIMENT.all.sac2_ind_finish_trl(:,trial_ind) = double(ind_sac2_finish);

EXPERIMENT.all.sac2_px(:,trial_ind)             = double(EXPERIMENT.all.eye_px(inds_sac2,trial_ind));
EXPERIMENT.all.sac2_py(:,trial_ind)             = double(EXPERIMENT.all.eye_py(inds_sac2,trial_ind));
EXPERIMENT.all.sac2_vx(:,trial_ind)             = double(EXPERIMENT.all.eye_vx(inds_sac2,trial_ind));
EXPERIMENT.all.sac2_vy(:,trial_ind)             = double(EXPERIMENT.all.eye_vy(inds_sac2,trial_ind));
EXPERIMENT.all.sac2_vm(:,trial_ind)             = double(EXPERIMENT.all.eye_vm(inds_sac2,trial_ind));
EXPERIMENT.all.sac2_vm_max(:,trial_ind)         = double(EXPERIMENT.all.eye_vm(ind_sac2_vmax,trial_ind));
EXPERIMENT.all.sac2_px_s(:,trial_ind)           = EXPERIMENT.all.sac2_px(:,trial_ind) - EXPERIMENT.all.tgt_str_x(:,trial_ind);
EXPERIMENT.all.sac2_py_s(:,trial_ind)           = EXPERIMENT.all.sac2_py(:,trial_ind) - EXPERIMENT.all.tgt_str_y(:,trial_ind);

EXPERIMENT.all.sac2_px_sacStart(:,trial_ind)    = double(EXPERIMENT.all.eye_px(ind_sac2_start,trial_ind));
EXPERIMENT.all.sac2_px_sacFinish(:,trial_ind)   = double(EXPERIMENT.all.eye_px(ind_sac2_finish,trial_ind));
EXPERIMENT.all.sac2_py_sacStart(:,trial_ind)    = double(EXPERIMENT.all.eye_py(ind_sac2_start,trial_ind));
EXPERIMENT.all.sac2_py_sacFinish(:,trial_ind)   = double(EXPERIMENT.all.eye_py(ind_sac2_finish,trial_ind));
EXPERIMENT.all.sac2_px_sacStart_s(:,trial_ind)  = double(EXPERIMENT.all.eye_px_s(ind_sac2_start,trial_ind));
EXPERIMENT.all.sac2_px_sacFinish_s(:,trial_ind) = double(EXPERIMENT.all.eye_px_s(ind_sac2_finish,trial_ind));
EXPERIMENT.all.sac2_py_sacStart_s(:,trial_ind)  = double(EXPERIMENT.all.eye_py_s(ind_sac2_start,trial_ind));
EXPERIMENT.all.sac2_py_sacFinish_s(:,trial_ind) = double(EXPERIMENT.all.eye_py_s(ind_sac2_finish,trial_ind));
EXPERIMENT.all.sac2_amp_x(:,trial_ind)          = double(abs(EXPERIMENT.all.sac2_px_sacFinish(:,trial_ind) - EXPERIMENT.all.sac2_px_sacStart(:,trial_ind)));
EXPERIMENT.all.sac2_amp_y(:,trial_ind)          = double(abs(EXPERIMENT.all.sac2_py_sacFinish(:,trial_ind) - EXPERIMENT.all.sac2_py_sacStart(:,trial_ind)));
EXPERIMENT.all.sac2_amp_m(:,trial_ind)          = double(sqrt(EXPERIMENT.all.sac2_amp_x(:,trial_ind)^2+EXPERIMENT.all.sac2_amp_y(:,trial_ind)^2));
EXPERIMENT.all.sac2_reaction(:,trial_ind)       = EXPERIMENT.all.sac2_ind_start_trl(:,trial_ind) - EXPERIMENT.all.sac_ind_finish_trl(:,trial_ind);
