function TRIAL = ESN_Sac_Analyser(data, trial_num_in_block, block_num, trial_num_in_experiment, params)
DEBUG = false;
if DEBUG
    trial_num_in_block = counter_trials;
    block_num = counter_main;
    trial_num_in_experiment = trial_number;
end


if nargin < 5
    params.Treshold_sac2_pm = 5.0; %deg maximum amplitude of corrective saccade
    params.Treshold_sac2_vmax = 600; %deg/s maximum velocity of corrective saccade
    params.Treshold_eye_px_start_max = 2.0; %deg start location of sac
    params.Treshold_eye_px_finish_min = 4.0; %deg minimum end location of sac
    params.Treshold_eye_px_finish_max = 14.0; %deg maximum end location of sac
end

%% Extract trial data
ind_START_current = data.inds.time_START(trial_num_in_block);
ind_START_next = data.inds.time_START(trial_num_in_block+1);
inds_trial             = ind_START_current : ind_START_next-1;
time_blk               = data.EL.time_1K(inds_trial);
states                 = data.PC.state_1K(inds_trial); states(end) = 2;
states_diff            = diff(states(~isnan(states)));
states_change          = find(states_diff);
ind_stateSacBegin      = states_change(find(states(states_change+1)==data.inds.ind_state_saccade, 1, 'last'))+1;
ind_stateSacEnd        = states_change(find(states(states_change)  ==data.inds.ind_state_saccade, 1, 'last'));
ind_stateFixateBegin   = states_change(find(states(states_change+1)==data.inds.ind_state_fixate,  1, 'last'))+1;
ind_stateFixateEnd     = states_change(find(states(states_change)  ==data.inds.ind_state_fixate,  1, 'last'));
if isempty(ind_stateFixateBegin)
    ind_stateFixateBegin   = states_change(find(states(states_change+1)==data.inds.ind_state_fixate_errclmp,  1, 'last'))+1;
    ind_stateFixateEnd     = states_change(find(states(states_change)  ==data.inds.ind_state_fixate_errclmp,  1, 'last'));
end

TRIAL.block_num    = uint32(block_num);
TRIAL.trial_num    = uint32(trial_num_in_experiment);

TRIAL.tgt_str_x                = double(data.PC.start_x(  trial_num_in_block));
TRIAL.tgt_str_y                = double(data.PC.start_y(  trial_num_in_block));
TRIAL.tgt_cue_x                = double(data.PC.cue_x(    trial_num_in_block));
TRIAL.tgt_cue_y                = double(data.PC.cue_y(    trial_num_in_block));
TRIAL.tgt_end_x                = double(data.PC.end_x(    trial_num_in_block));
TRIAL.tgt_end_y                = double(data.PC.end_y(    trial_num_in_block));
TRIAL.tgt_amp_x                = double(data.PC.amp_x(    trial_num_in_block));
TRIAL.tgt_amp_y                = double(data.PC.amp_y(    trial_num_in_block));
TRIAL.tgt_err_clp              = uint32(data.PC.err_clmp( trial_num_in_block));
TRIAL.tgt_img_str              = uint32(data.PC.img_start(trial_num_in_block));
TRIAL.tgt_img_cue              = uint32(data.PC.img_cue(  trial_num_in_block));
TRIAL.tgt_img_end              = uint32(data.PC.img_end(  trial_num_in_block));
TRIAL.aim_px                   = double(data.PC.px_1K(inds_trial));
TRIAL.aim_py                   = double(data.PC.py_1K(inds_trial));

TRIAL.inds_trial               = uint32(inds_trial);
TRIAL.ind_stateSacBegin_trl    = uint32(ind_stateSacBegin);
TRIAL.ind_stateSacEnd_trl      = uint32(ind_stateSacEnd);
TRIAL.ind_stateFixateBegin_trl = uint32(ind_stateFixateBegin);
TRIAL.ind_stateFixateEnd_trl   = uint32(ind_stateFixateEnd);
TRIAL.ind_START_blk            = uint32(inds_trial(1));
TRIAL.ind_FIXATE_blk           = uint32(inds_trial(ind_stateFixateBegin));
TRIAL.time_START_blk           = double(time_blk(1));
TRIAL.time_FIXATE_blk          = double(time_blk(ind_stateFixateBegin));
TRIAL.time_blk                 = double(time_blk);
TRIAL.states                   = uint32(states);

TRIAL.eye_interp_pupil         = double(data.EL.interp_pupil(inds_trial));
TRIAL.eye_interp_px            = double(data.EL.interp_px(inds_trial));
TRIAL.eye_interp_py            = double(data.EL.interp_py(inds_trial));

TRIAL.eye_px     = double(data.EL.px_filt(inds_trial));
TRIAL.eye_py     = double(data.EL.py_filt(inds_trial));
TRIAL.eye_vx     = double(data.EL.vx_filt(inds_trial));
TRIAL.eye_vy     = double(data.EL.vy_filt(inds_trial));
TRIAL.eye_vm     = double(data.EL.vm_filt(inds_trial));
TRIAL.eye_px_s   = TRIAL.eye_px - TRIAL.tgt_str_x;
TRIAL.eye_py_s   = TRIAL.eye_py - TRIAL.tgt_str_y;

%% Extract sac
sampling_freq = 1000.0; cutoff_freq = 25.0;
[b_butter,a_butter] = butter(2,(cutoff_freq/(sampling_freq/2)), 'low');

sac_validity         = true;
% sac_validity_err     = 0;
inds_sac        = (ind_stateSacEnd-100):1:min([(ind_stateFixateBegin+100), length(inds_trial)]);
sac_vm         = double(TRIAL.eye_vm(inds_sac));
sac_vm_filt    = filtfilt(b_butter,a_butter,sac_vm); % sgolayfilt(sac_vm, 3, 25);

[~, ind_sac_vmax] = findpeaks(sac_vm_filt, 'MinPeakProminence',100, 'MinPeakHeight', 150);
if(sum(diff(ind_sac_vmax)<80))
    sac_validity         = false;
%     sac_validity_err     = 1;
end

[sac_vmax, ind_sac_vmax] = findpeaks(sac_vm_filt, 'MinPeakProminence',100,'SortStr','descend', 'NPeaks', 1, 'MinPeakHeight', 150);

if(isempty(sac_vmax))
    sac_validity         = false;
%     sac_validity_err     = 2;
    ind_sac_vmax         = round((ind_stateSacEnd+ind_stateFixateBegin)/2);
    inds_sac        = (ind_sac_vmax -61) : 1 : min([(ind_sac_vmax -61+149), length(inds_trial)]);
    ind_sac_start        = ind_sac_vmax - 50;
    ind_sac_finish       = ind_sac_vmax + 50;
    sac_vmax             = double(TRIAL.eye_vm(ind_sac_vmax));
else
    if (ind_sac_vmax+10) > (length(sac_vm))
        ind_sac_vmax_ = 9; % handling an error in which the max happened in the end
    else
        [~, ind_sac_vmax_] = max(sac_vm(ind_sac_vmax-10:ind_sac_vmax+10));
    end
    
    ind_sac_vmax_   = ind_sac_vmax_ - 1 + (ind_sac_vmax-10);
    ind_sac_vmax    = ind_sac_vmax_ + (ind_stateSacEnd-100);
    inds_sac        = (ind_sac_vmax - 100) : 1 : min([(ind_sac_vmax - 100 + 299), length(inds_trial)]);
    sac_vm         = double(TRIAL.eye_vm(inds_sac));
    sac_vm_filt    = filtfilt(b_butter,a_butter,sac_vm);
    
    ind_sac_start_50        = find(sac_vm_filt(1:100)<50, 1, 'last');
    if(isempty(ind_sac_start_50))
        ind_sac_start     = 1;
        sac_validity      = false;
%         sac_validity_err  = 3;
    else
        ind_sac_start_20    = find(sac_vm_filt(1:100)<20, 1, 'last');
        if(isempty(ind_sac_start_20))
            ind_sac_start_       = find(sac_vm(max([(ind_sac_start_50-10), 1]) : min([(ind_sac_start_50+10), length(sac_vm)]) )<50, 1, 'last') - 1 + (max([(ind_sac_start_50-10), 1]) );
            if(isempty(ind_sac_start_))
                ind_sac_start    = ind_sac_start_50;
            else
                ind_sac_start    = ind_sac_start_;
            end
%             sac_validity_err     = 4;
        else
            ind_sac_start_       = find(sac_vm(max([(ind_sac_start_20-10), 1]) : min([(ind_sac_start_20+10), length(sac_vm)]) )<20, 1, 'last') - 1 + (max([(ind_sac_start_20-10), 1]) );
            if(isempty(ind_sac_start_))
                ind_sac_start    = ind_sac_start_20;
            else
                ind_sac_start    = ind_sac_start_;
            end
        end
    end
    
    ind_sac_finish_50       = find(sac_vm_filt(100:end)<50, 1, 'first') - 1 + 100;
    if(isempty(ind_sac_finish_50))
        ind_sac_finish      = length(sac_vm_filt);
        sac_validity        = false;
%         sac_validity_err    = 5;
    else
        ind_sac_finish_20   = find(sac_vm_filt( 100 : min([(ind_sac_finish_50+50), length(sac_vm_filt)]) )<20, 1, 'first') - 1 + 100;
        if(isempty(ind_sac_finish_20))
            ind_sac_finish_       = find(sac_vm(max([(ind_sac_finish_50-10), 1]) : min([(ind_sac_finish_50+10), length(sac_vm)]) )<50, 1, 'first') - 1 + (max([(ind_sac_finish_50-10), 1]) );
            if(isempty(ind_sac_finish_))
                ind_sac_finish    = ind_sac_finish_50;
            else
                ind_sac_finish    = ind_sac_finish_;
            end
        else
            ind_sac_finish_       = find(sac_vm(max([(ind_sac_finish_20-10), 1]) : min([(ind_sac_finish_20+10), length(sac_vm)]) )<20, 1, 'first') - 1 + (max([(ind_sac_finish_20-10), 1]));
            if(isempty(ind_sac_finish_))
                ind_sac_finish    = ind_sac_finish_20;
            else
                ind_sac_finish    = ind_sac_finish_;
            end
        end
    end
    
    ind_sac_start        = ind_sac_vmax - 100 + ind_sac_start;
    ind_sac_finish       = ind_sac_vmax - 100 + ind_sac_finish;
    inds_sac        = (ind_sac_vmax -61) : 1 : min([(ind_sac_vmax -61+149), length(inds_trial)]);
end

if(abs(double(TRIAL.eye_px_s(ind_sac_start))) > params.Treshold_eye_px_start_max)
    sac_validity     = false;
%     sac_validity_err = 9;
elseif(abs(double(TRIAL.eye_px_s(ind_sac_finish))) < params.Treshold_eye_px_finish_min)
    sac_validity     = false;
%     sac_validity_err = 10;
elseif(abs(double(TRIAL.eye_px_s(ind_sac_finish))) > params.Treshold_eye_px_finish_max)
    sac_validity     = false;
%     sac_validity_err = 11;
end

%% Extract sac2
sac2_validity        = true;
inds_sac2            = ind_sac_finish:1:ind_stateFixateEnd;
if isempty(inds_sac2)
    inds_sac2 = length(inds_trial) - 200 : 1 : length(inds_trial);
    sac2_vmax = [];
else
    sac2_vm        = double(TRIAL.eye_vm(inds_sac2));
    sac2_vm_filt   = filtfilt(b_butter,a_butter,sac2_vm);
    [sac2_vmax,ind_sac2_vmax]=findpeaks(sac2_vm_filt, 'MinPeakProminence',40,'SortStr','descend', 'NPeaks', 1, 'MinPeakHeight', 80);

end

if(isempty(sac2_vmax))
    sac2_validity        = false;
    ind_sac2_vmax        = round((inds_sac2(1,1)+inds_sac2(1,end))/2);
    inds_sac2       = (ind_sac2_vmax-61) : 1 : min([(ind_sac2_vmax-61+149), length(inds_trial)] );
    ind_sac2_start       = ind_sac2_vmax - 25;
    ind_sac2_finish      = min([(ind_sac2_vmax + 25) length(inds_trial)]);
    sac2_vmax            = double(TRIAL.eye_vm(ind_sac2_vmax));
else
    if (ind_sac2_vmax+10) > (length(sac2_vm))
        ind_sac2_vmax_ = 9; % handling an error in which the max happened in the end
    else
        [~, ind_sac2_vmax_]  = max(sac2_vm(ind_sac2_vmax-10:ind_sac2_vmax+10));
    end
    ind_sac2_vmax_       = ind_sac2_vmax_ - 1 + (ind_sac2_vmax-10);
    ind_sac2_vmax        = ind_sac2_vmax_ + ind_sac_finish;
    inds_sac2            = (ind_sac2_vmax - 100) : 1 : min([(ind_sac2_vmax - 100 + 299), length(inds_trial)] );
    sac2_vm        = double(TRIAL.eye_vm(inds_sac2));
    sac2_vm_filt   = filtfilt(b_butter,a_butter,sac2_vm);
    
    ind_sac2_start_30        = find(sac2_vm_filt(1:100)<30, 1, 'last');
    if(isempty(ind_sac2_start_30))
        ind_sac2_start=1;
        sac2_validity=false;
    else
        ind_sac2_start_20    = find(sac2_vm_filt(1:100)<20, 1, 'last');
        if(isempty(ind_sac2_start_20))
            ind_sac2_start_       = find(sac2_vm(max([(ind_sac2_start_30-10), 1]) : min([(ind_sac2_start_30+10), length(sac2_vm)]) )<30, 1, 'last') - 1 + (max([(ind_sac2_start_30-10), 1]));
            if(isempty(ind_sac2_start_))
                ind_sac2_start    = ind_sac2_start_30;
            else
                ind_sac2_start    = ind_sac2_start_;
            end
        else
            ind_sac2_start_       = find(sac2_vm(max([(ind_sac2_start_20-10), 1]) : min([(ind_sac2_start_20+10), length(sac2_vm)]) )<20, 1, 'last') - 1 + (max([(ind_sac2_start_20-10), 1]) );
            if(isempty(ind_sac2_start_))
                ind_sac2_start    = ind_sac2_start_20;
            else
                ind_sac2_start    = ind_sac2_start_;
            end
        end
    end
    
    ind_sac2_finish_30       = find(sac2_vm_filt(100:end)<30, 1, 'first') - 1 + 100;
    if(isempty(ind_sac2_finish_30))
        ind_sac2_finish = length(sac2_vm_filt);
    else
        ind_sac2_finish_20   = find(sac2_vm_filt( 100 : min([(ind_sac2_finish_30+50), length(sac2_vm_filt)]) )<20, 1, 'first') + 100 - 1;
        if(isempty(ind_sac2_finish_20))
            ind_sac2_finish_       = find(sac2_vm(max([(ind_sac2_finish_30-10), 1]) : min([(ind_sac2_finish_30+10), length(sac2_vm)]) )<30, 1, 'first') - 1 + (max([(ind_sac2_finish_30-10), 1]) );
            if(isempty(ind_sac2_finish_))
                ind_sac2_finish    = ind_sac2_finish_30;
            else
                ind_sac2_finish    = ind_sac2_finish_;
            end
        else
            ind_sac2_finish_       = find(sac2_vm(max([(ind_sac2_finish_20-10), 1]) : min([(ind_sac2_finish_20+10), length(sac2_vm)]) )<20, 1, 'first') - 1 + (max([(ind_sac2_finish_20-10), 1]) );
            if(isempty(ind_sac2_finish_))
                ind_sac2_finish    = ind_sac2_finish_20;
            else
                ind_sac2_finish    = ind_sac2_finish_;
            end
        end
    end
    ind_sac2_start       = ind_sac2_vmax - 100 + ind_sac2_start;
    ind_sac2_finish      = min([(ind_sac2_vmax - 100 + ind_sac2_finish) length(inds_trial)]);
    inds_sac2       = (ind_sac2_vmax-61) : 1 : min([(ind_sac2_vmax-61+149), length(inds_trial)] );
end

sac2_pm = abs( double(TRIAL.eye_px_s(ind_sac2_finish)) - double(TRIAL.eye_px_s(ind_sac2_start)) );
if( sac2_pm < 0.5 )
    sac2_validity=false;
elseif( sac2_pm > params.Treshold_sac2_pm )
    sac2_validity=false;
elseif( sac2_vmax > params.Treshold_sac2_vmax )
    sac2_validity=false;
end


% Make Sure that inds_sac2 does not exceeds the dimension of inds_trial
inds_sac2(inds_sac2>length(inds_trial)) = length(inds_trial);

%% Save sac
TRIAL.sac_validity        = sac_validity;
TRIAL.sac_inds            = uint32(inds_sac);
TRIAL.sac_ind_start_trl   = uint32(ind_sac_start);
TRIAL.sac_ind_vmax_trl    = uint32(ind_sac_vmax);
TRIAL.sac_ind_finish_trl  = uint32(ind_sac_finish);

TRIAL.sac_px              = double(TRIAL.eye_px(inds_sac));
TRIAL.sac_py              = double(TRIAL.eye_py(inds_sac));
TRIAL.sac_vx              = double(TRIAL.eye_vx(inds_sac));
TRIAL.sac_vy              = double(TRIAL.eye_vy(inds_sac));
TRIAL.sac_vm              = double(TRIAL.eye_vm(inds_sac));
TRIAL.sac_vm_max          = double(TRIAL.eye_vm(ind_sac_vmax));
TRIAL.sac_px_s            = TRIAL.sac_px - TRIAL.tgt_str_x;
TRIAL.sac_py_s            = TRIAL.sac_py - TRIAL.tgt_str_y;

TRIAL.sac_px_sacStart     = double(TRIAL.eye_px(ind_sac_start));
TRIAL.sac_px_sacFinish    = double(TRIAL.eye_px(ind_sac_finish));
TRIAL.sac_py_sacStart     = double(TRIAL.eye_py(ind_sac_start));
TRIAL.sac_py_sacFinish    = double(TRIAL.eye_py(ind_sac_finish));
TRIAL.sac_px_sacStart_s   = double(TRIAL.eye_px_s(ind_sac_start));
TRIAL.sac_px_sacFinish_s  = double(TRIAL.eye_px_s(ind_sac_finish));
TRIAL.sac_py_sacStart_s   = double(TRIAL.eye_py_s(ind_sac_start));
TRIAL.sac_py_sacFinish_s  = double(TRIAL.eye_py_s(ind_sac_finish));
TRIAL.sac_amp_x           = double(TRIAL.sac_px_sacFinish - TRIAL.sac_px_sacStart);
TRIAL.sac_amp_y           = double(TRIAL.sac_py_sacFinish - TRIAL.sac_py_sacStart);
TRIAL.sac_amp_m           = double(sqrt(TRIAL.sac_amp_x^2+TRIAL.sac_amp_y^2));
TRIAL.sac_reaction        = TRIAL.sac_ind_start_trl - TRIAL.ind_stateSacBegin_trl;

%% Save sac2
TRIAL.sac2_validity       = sac2_validity;
TRIAL.sac2_inds           = uint32(inds_sac2);
TRIAL.sac2_ind_start_trl  = uint32(ind_sac2_start);
TRIAL.sac2_ind_vmax_trl   = uint32(ind_sac2_vmax);
TRIAL.sac2_ind_finish_trl = uint32(ind_sac2_finish);

TRIAL.sac2_px             = double(TRIAL.eye_px(inds_sac2));
TRIAL.sac2_py             = double(TRIAL.eye_py(inds_sac2));
TRIAL.sac2_vx             = double(TRIAL.eye_vx(inds_sac2));
TRIAL.sac2_vy             = double(TRIAL.eye_vy(inds_sac2));
TRIAL.sac2_vm             = double(TRIAL.eye_vm(inds_sac2));
TRIAL.sac2_vm_max         = double(TRIAL.eye_vm(ind_sac2_vmax));
TRIAL.sac2_px_s           = TRIAL.sac2_px - TRIAL.tgt_str_x;
TRIAL.sac2_py_s           = TRIAL.sac2_py - TRIAL.tgt_str_y;

TRIAL.sac2_px_sacStart    = double(TRIAL.eye_px(ind_sac2_start));
TRIAL.sac2_px_sacFinish   = double(TRIAL.eye_px(ind_sac2_finish));
TRIAL.sac2_py_sacStart    = double(TRIAL.eye_py(ind_sac2_start));
TRIAL.sac2_py_sacFinish   = double(TRIAL.eye_py(ind_sac2_finish));
TRIAL.sac2_px_sacStart_s  = double(TRIAL.eye_px_s(ind_sac2_start));
TRIAL.sac2_px_sacFinish_s = double(TRIAL.eye_px_s(ind_sac2_finish));
TRIAL.sac2_py_sacStart_s  = double(TRIAL.eye_py_s(ind_sac2_start));
TRIAL.sac2_py_sacFinish_s = double(TRIAL.eye_py_s(ind_sac2_finish));
TRIAL.sac2_amp_x          = double(abs(TRIAL.sac2_px_sacFinish - TRIAL.sac2_px_sacStart));
TRIAL.sac2_amp_y          = double(abs(TRIAL.sac2_py_sacFinish - TRIAL.sac2_py_sacStart));
TRIAL.sac2_amp_m          = double(sqrt(TRIAL.sac2_amp_x^2+TRIAL.sac2_amp_y^2));
TRIAL.sac2_reaction       = TRIAL.sac2_ind_start_trl - TRIAL.sac_ind_finish_trl;

end