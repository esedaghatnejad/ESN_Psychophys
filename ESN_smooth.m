function smooth_data_ = ESN_smooth(data_, dim)
% smooth data using 2nd order Savitzky-Golay filter with 21 points
% if data_ is a matrix, the method will smooth each column by default or smooth along dim.
% method = 'moving';  % Moving average. A lowpass filter with filter coefficients equal to the reciprocal of the span.
% method = 'lowess';  % Local regression using weighted linear least squares and a 1st degree polynomial model.
% method = 'loess';   % Local regression using weighted linear least squares and a 2nd degree polynomial model.
% method = 'sgolay';  % Savitzky-Golay filter. A generalized moving average with filter coefficients determined by an unweighted linear least-squares regression and a polynomial model of specified degree (default is 2). The method can accept nonuniform predictor data.
% method = 'rlowess'; % A robust version of 'lowess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
% method = 'rloess';  % A robust version of 'loess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
% smooth_data_ = smooth(data_, method);
if nargin < 2
    dim = 1;
end
if size(data_, 1) == 1
    smooth_data_ = reshape(smooth(data_, 21, 'sgolay', 2), 1, []);
elseif size(data_, 2) == 1
    smooth_data_ = reshape(smooth(data_, 21, 'sgolay', 2), [], 1);
else
    smooth_data_ = nan(size(data_));
    if dim == 1
        % smooth columns
        for counter_col = 1 : size(data_, 2)
            smooth_data_(:, counter_col) = reshape(smooth(data_(:, counter_col), 21, 'sgolay', 2), [], 1);
        end
    elseif dim == 2
        % smooth rows
        for counter_row = 1 : size(data_, 1)
            smooth_data_(counter_row, :) = reshape(smooth(data_(counter_row, :), 21, 'sgolay', 2), 1, []);
        end
    end
    
end
end
