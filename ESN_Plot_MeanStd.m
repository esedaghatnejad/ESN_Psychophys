function ax_handles = ESN_Plot_MeanStd(X_Values, Y_Mean, Y_Std, Color_Vec, params)
% Author: Ehsan Sedaghat-Nejad (esedaghatnejad@gmail.com)
% This function will draw mean and shade the std around the mean with the same color.
% Take a look at input structure 'params' for ajusting the plot parameters
% ax_handles can be used for the legend box since this function will generate multiple axes
% input vectors should prefentially be column-wise and in case of having more than one vector, each
% column will be considered as an entree (similar to plot function)
% You can use this function in the following forms
% ESN_Plot_MeanStd(Y_Mean)
% ESN_Plot_MeanStd(X_Values, Y_Mean)
% ESN_Plot_MeanStd(X_Values, Y_Mean, Y_Std)
% ESN_Plot_MeanStd([], Y_Mean, Y_Std)
% ESN_Plot_MeanStd(X_Values, Y_Mean, Y_Std, Color_Vec)
% ESN_Plot_MeanStd(X_Values, Y_Mean, Y_Std, Color_Vec, params)

% handling the color when there is no input Color_Vec so the function will change the line colors
% each time it is called
persistent line_color_index
if isempty(line_color_index)
    line_color_index = 0;
end
% hold on
hold on
%% Check the number of inputs
% Use matlab line colormap
if nargin < 2
    Y_Mean = X_Values;
    X_Values = [];
end
if nargin < 3
    Y_Std = [];
end
if nargin < 4
    Color_Vec = lines(7);
end
if nargin < 5
    params.marker_on_off = 1;
    params.std_boundry_on_off = 1;
    params.linewidth = 2;
    params.MarkerSize = 4;
    params.Std_Alpha = 0.35;
    params.Line_Alpha = 0.8;
end
if nargin > 3
    line_color_index = 0;
else
    line_color_index = line_color_index + 1;
end
%% Check the params structure
% check if the params stucture have all the necessary fields
if ~isfield(params, 'marker_on_off')
    params.marker_on_off = 1;
end
if ~isfield(params, 'linewidth')
    params.linewidth = 2;
end
if ~isfield(params, 'MarkerSize')
    params.MarkerSize = 4;
end
if ~isfield(params, 'Std_Alpha')
    params.Std_Alpha = 0.35;
end
if ~isfield(params, 'Line_Alpha')
    params.Line_Alpha = 0.8;
end
%% Check vector dimensions
% use the index number as x-values 
if isempty(X_Values)
    X_Values = (1:1:size(Y_Mean,1))';
end
if isempty(Y_Std)
    Y_Std = zeros(size(Y_Mean));
end
if size(X_Values, 2) < size(Y_Mean, 2)
    X_Values = repmat(X_Values, 1, size(Y_Mean, 2));
end
% make the vector column-wise
if (size(X_Values, 1) == 1)
    X_Values = X_Values';
end
% make the vector column-wise
if (size(Y_Mean, 1) == 1)
    Y_Mean = Y_Mean';
end
% make the vector column-wise
if (size(Y_Std, 1) == 1)
    Y_Std = Y_Std';
end
%% Plotting
for counter = 1 : size(Y_Mean, 2)
    if line_color_index == 0
        color_index = mod(counter-1, size(Color_Vec, 1))+1;
    else
        color_index = mod(line_color_index-1, size(Color_Vec, 1))+1;
        line_color_index = line_color_index + 1;
    end
    % prepare the variables for plotting
    X_Values_ = X_Values(:, counter);
    Y_Mean_ = Y_Mean(:, counter);
    Y_Std_ = Y_Std(:, counter);
    invalid_ind = isnan(X_Values_) | isnan(Y_Mean_) | isnan(Y_Std_);
    X_Values_ = X_Values_(~invalid_ind); Y_Mean_ = Y_Mean_(~invalid_ind); Y_Std_ = Y_Std_(~invalid_ind);
    y_min = Y_Mean_ - Y_Std_;
    y_max = Y_Mean_ + Y_Std_;
    x_fill = [X_Values_; X_Values_(end:-1:1)];
    y_fill = [y_max; y_min(end:-1:1)];
    % plotting std shades using fill function
    fill(x_fill, y_fill, 'c',...
        'FaceColor', Color_Vec(color_index, :),...
        'FaceAlpha', params.Std_Alpha,...
        'LineStyle','none');
    % Plotting the line using patch function which accepts transparency
    patch([X_Values_;NaN],[Y_Mean_;NaN],'w',...
        'linewidth',params.linewidth,...
        'EdgeAlpha',params.Line_Alpha,...
        'EdgeColor',Color_Vec(color_index, :));
    if params.std_boundry_on_off
        patch([X_Values_;NaN],[y_min;NaN],'w',...
        'linewidth',1,...
        'EdgeAlpha',params.Line_Alpha,...
        'EdgeColor',Color_Vec(color_index, :));
        patch([X_Values_;NaN],[y_max;NaN],'w',...
        'linewidth',1,...
        'EdgeAlpha',params.Line_Alpha,...
        'EdgeColor',Color_Vec(color_index, :));
    end
    % plotting markers
    if params.marker_on_off
        % plotting markers
        plot(X_Values_, Y_Mean_,'o',...
            'MarkerSize',params.MarkerSize, ...
            'MarkerFaceColor', Color_Vec(color_index, :), ...
            'MarkerEdgeColor', 'k');
    end
    % This plot is just for the purpose of a better representation in the legend box and does not
    % add any visualization to the figure
    ax_handles(counter, 1) = plot(nan, nan,'-',...
        'Color', Color_Vec(color_index, :),...
        'linewidth',params.linewidth);
end
% hold off
hold off

end