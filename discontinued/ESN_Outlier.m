function data = ESN_Outlier(data, movmedian_length)
% Take out outliers using recursive evaluation of data utill no outlier remains
% This function uses isoutlier which has been added since MATLAB 2017a
% The primary invalid values should be nan and the function will ignore those values.
% The output will be the same length as input but with nan values as outliers.
% The input should be a vector array and the window length for moving median.

if nargin < 2
    movmedian_length = ceil(length(data)/10);
end

flag_while = true;
while flag_while
    inds_outlier = isoutlier(data,'movmedian',movmedian_length);
    data(inds_outlier) = NaN;
    flag_while = ~(sum(inds_outlier)==0);
end
flag_while = true;
while flag_while
    inds_outlier = isoutlier(data);
    data(inds_outlier) = NaN;
    flag_while = ~(sum(inds_outlier)==0);
end
