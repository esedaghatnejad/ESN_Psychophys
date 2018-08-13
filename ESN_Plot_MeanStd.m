function ax_mean = ESN_Plot_MeanStd(X_Values, Y_Mean, Y_Std, Color_Vec, marker_on_off)
hold on
if nargin < 5
    marker_on_off = 1;
end

if (size(X_Values, 1) == 1)
    X_Values = X_Values';
end
if (size(Y_Mean, 1) == 1)
    Y_Mean = Y_Mean';
end
if (size(Y_Std, 1) == 1)
    Y_Std = Y_Std';
end

if isempty(X_Values)
    X_Values = (1:1:size(Y_Mean,1))';
end

if isempty(Y_Std)
    Y_Std = zeros(size(Y_Mean));
end

if size(X_Values, 2) < size(Y_Mean, 2)
    X_Values = repmat(X_Values, 1, size(Y_Mean, 2));
end

% if size(X_Values, 2) > size(X_Values, 1)
%     X_Values = X_Values';
% end
% if size(Y_Mean, 2) > size(Y_Mean, 1)
%     Y_Mean = Y_Mean';
% end
% if size(Y_Std, 2) > size(Y_Std, 1)
%     Y_Std = Y_Std';
% end

for counter = 1 : size(Y_Mean, 2)
    color_index = mod(counter-1, size(Color_Vec, 1))+1;
    
    X_Values_ = X_Values(:, counter);
    Y_Mean_ = Y_Mean(:, counter);
    Y_Std_ = Y_Std(:, counter);
    invalid_ind = isnan(X_Values_) | isnan(Y_Mean_) | isnan(Y_Std_);
    X_Values_ = X_Values_(~invalid_ind); Y_Mean_ = Y_Mean_(~invalid_ind); Y_Std_ = Y_Std_(~invalid_ind);
    y_min = Y_Mean_ - Y_Std_;
    y_max = Y_Mean_ + Y_Std_;
    x_fill = [X_Values_; X_Values_(end:-1:1)];
    y_fill = [y_max; y_min(end:-1:1)];
    fill(x_fill, y_fill, 'c','FaceColor', Color_Vec(color_index, :),'FaceAlpha', 0.35,'LineStyle','none');
    if marker_on_off
        ax_mean(counter, 1) = plot(X_Values_, Y_Mean_,'o','MarkerSize',8, 'MarkerFaceColor', Color_Vec(color_index, :), 'MarkerEdgeColor', 'k');
    end
    patch([X_Values_;NaN],[Y_Mean_;NaN],'w','linewidth',2,'EdgeAlpha',0.8,'EdgeColor',Color_Vec(color_index, :));
end
hold off

end