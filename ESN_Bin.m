function [x_data_BIN_mean, y_data_BIN_mean, x_data_BIN_stdv, y_data_BIN_stdv] = ESN_Bin(x_data_, y_data_, edges_, func1_, func2_)

if nargin < 4
    func1_ = @nanmean;
    func2_ = @nanstd;
elseif nargin == 4
    func2_ = @std;
end

[x_data_Ns,~,x_axis_data_bins]    = histcounts(x_data_(:)', edges_(:)');
inds_out_of_bound = x_axis_data_bins < 1;
x_axis_data_bins(inds_out_of_bound) = []; x_data_(inds_out_of_bound) = []; y_data_(inds_out_of_bound) = [];
x_data_BIN_mean  = accumarray(double(x_axis_data_bins(:)), double(x_data_(:)), [], func1_, NaN);
y_data_BIN_mean  = accumarray(double(x_axis_data_bins(:)), double(y_data_(:)), [], func1_, NaN);
x_data_BIN_stdv  = accumarray(double(x_axis_data_bins(:)), double(x_data_(:)), [], func2_, NaN);
y_data_BIN_stdv  = accumarray(double(x_axis_data_bins(:)), double(y_data_(:)), [], func2_, NaN);

length_x_data = length(x_data_BIN_mean);
length_Ns = length(x_data_Ns);

if(length_x_data < length_Ns )
    x_data_BIN_mean(length_x_data+1:length_Ns) = nan;
    y_data_BIN_mean(length_x_data+1:length_Ns) = nan;
    x_data_BIN_stdv(length_x_data+1:length_Ns) = nan;
    y_data_BIN_stdv(length_x_data+1:length_Ns) = nan;
end

% x_data_BIN_mean( x_data_Ns<(sum(x_data_Ns)./(2*length(x_data_Ns)) ) ) = nan;
% y_data_BIN_mean( x_data_Ns<(sum(x_data_Ns)./(2*length(x_data_Ns)) ) ) = nan;
% x_data_BIN_stdv( x_data_Ns<(sum(x_data_Ns)./(2*length(x_data_Ns)) ) ) = nan;
% y_data_BIN_stdv( x_data_Ns<(sum(x_data_Ns)./(2*length(x_data_Ns)) ) ) = nan;


end