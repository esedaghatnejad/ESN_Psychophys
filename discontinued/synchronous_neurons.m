function [spiketrains] = synchronous_neurons()
% What is the timestep in seconds ? 
dt = 1/1000; % 1 ms

% How long should we run the simulation in seconds ? 
T = 2500; 

% What are the firing rates of our two neurons (in spks/second) ? 
firing_rates = [50, 50];

% What is the mandatory refractory period for each neuron (in S) ? 
refractory_period = [5/1000, 5/1000];

% What axes should we use for the cross-correlogram (in seconds) ? 
cross_correlogram_time_axis = [-50 / 1000, 50 / 1000];
    
    % [spiketrains] = gen_synchronous_spiketrains
    % Generates a set of two spike trains with the passed firing rates
    % The output is an NxNt binary matrix containing the times when each
    % of the neurons fired. The sychronous factor should be between 0 (no
    % synchrony, neurons fire independently) or 1 (complete synchrony). 
    % The spike trains are generated using a Poisson random process.
    % Each neuron obeys its respective refractory period.
    function [spiketrains] = gen_synchronous_spiketrains(synchrony_fraction)
        % Convert our inputs into units of dt
        num_timesteps = ceil(T / dt);
        refractory_period_timesteps = ceil(refractory_period / dt);
        
        % Create our empty binary spiketrain
        spiketrains = false(2, num_timesteps);
        
        % Adjust our firing rates to take into account the refractory period
        adjusted_firing_rates = firing_rates ./ (1 - firing_rates .* refractory_period);
        
        % Generate our first neuron using a standard Poisson process
        last_spike_index = -Inf;
        for i = 1:num_timesteps
            if (i - last_spike_index) < refractory_period_timesteps(1)
                continue; % In the refractory period
            end
            spiketrains(1, i) = rand() < adjusted_firing_rates(1) * dt;
            if spiketrains(1, i)
                last_spike_index = i;
            end
        end
        
        % Generate our second neuron with 0ms synchrony from the first neuron
        last_spike_index = -Inf;
        for i = 1:num_timesteps
            if (i - last_spike_index < refractory_period_timesteps(2))
                continue; 
            end
            if spiketrains(1, i) && rand() < synchrony_fraction
                spiketrains(2, i) = true;
            else
                spiketrains(2, i) = rand() < adjusted_firing_rates(2) * dt;
            end
            if spiketrains(2, i)
                last_spike_index = i;
            end
        end
    end

    % plot_cross_correlogram
    % Generates a cross correlogram given a set of spiketrains an an input
    % The spiketrain should be a binary matrix with size NxNt (N neurons
    % by Nt time steps)
    function [] = plot_cross_correlogram(spiketrains)
        cross_correlogram_time_axis_samples = ceil(cross_correlogram_time_axis / dt);
        cross_correlogram_time_axis_samples = cross_correlogram_time_axis_samples(1)+1:1:cross_correlogram_time_axis_samples(2)-1;
        
        % Create our output variable
        y = zeros(length(cross_correlogram_time_axis_samples), 1);
        
        % Create each of our cross correlograms by shifting and binary
        % and'ing
        for i = 1:length(cross_correlogram_time_axis_samples)
            if cross_correlogram_time_axis_samples(i) > 0
                right_shift_index = cross_correlogram_time_axis_samples(i) + 1;
                temp_spiketrain_1 = spiketrains(1, 1:end-right_shift_index+1);
                temp_spiketrain_2 = spiketrains(2, right_shift_index:end);
            elseif cross_correlogram_time_axis_samples(i) < 0
                left_shift_index = -1 * cross_correlogram_time_axis_samples(i) + 1;
                temp_spiketrain_1 = spiketrains(1, left_shift_index:end);
                temp_spiketrain_2 = spiketrains(2, 1:end-left_shift_index+1);
            else
                temp_spiketrain_1 = spiketrains(1, :);
                temp_spiketrain_2 = spiketrains(2, :);
                % if spiketrain_2 == spiketrain_1 (for auto-correlogram),
                % arbitrarily set to zero
                if (all(temp_spiketrain_1 == temp_spiketrain_2))
                    temp_spiketrain_1 = temp_spiketrain_1 * 0;
                end
            end
            temp = (temp_spiketrain_1 & temp_spiketrain_2);
            y(i) = sum(temp) / sum(temp_spiketrain_1) / dt;
        end
        
        figure()
        bar(cross_correlogram_time_axis_samples, y, 1);
        xlim([cross_correlogram_time_axis_samples(1), cross_correlogram_time_axis_samples(end)]);
        xlabel("Samples")
        ylabel("Firing rate (spks/s)")
    end

    spiketrains = gen_synchronous_spiketrains(0);
    
    % Plot the auto-correlogram for neuron 2
    plot_cross_correlogram([spiketrains(2, :); spiketrains(2, :)]);
    title("Auto-correlogram")
    
    % Plot various cross correlograms given fractions of synchronous spike
    % trains
    spiketrains = gen_synchronous_spiketrains(0.0);
    plot_cross_correlogram(spiketrains);
    title("0% synchrony")
    
    spiketrains = gen_synchronous_spiketrains(0.1);
    plot_cross_correlogram(spiketrains);
    title("10% synchrony")
    
    spiketrains = gen_synchronous_spiketrains(0.25);
    plot_cross_correlogram(spiketrains);
    title("25% synchrony")
    
    spiketrains = gen_synchronous_spiketrains(0.50);
    plot_cross_correlogram(spiketrains);
    title("50% synchrony")
end