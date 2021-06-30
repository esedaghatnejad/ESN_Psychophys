function event_data_expand = ESN_expand_index_event_data(event_data, dim, expand_index, mode)
% get a logical matrix as input and expands the 1s along dim

%% handle nargin
if nargin < 2
    dim = 2; % dim=2; %expand along row; dim=1; %expand along col; 
end
if nargin < 3
    expand_index = 0;
end
if nargin < 4
    mode = 'centered'; % mode='centered'; % mode='forward'; % mode='backward';
end
%% convert event_data to logical
event_data = logical(event_data);
event_data_expand = logical(event_data);
%% expand the 1s based on the mode
if strcmp(mode,'centered')
    for counter = (-expand_index) : 1 : (expand_index)
        event_data_expand = event_data_expand | circshift(event_data, counter, dim);
    end
elseif strcmp(mode,'forward')
    for counter = 1 : 1 : (expand_index)
        event_data_expand = event_data_expand | circshift(event_data, counter, dim);
    end
elseif strcmp(mode,'backward')
    for counter = (-expand_index) : 1 : -1
        event_data_expand = event_data_expand | circshift(event_data, counter, dim);
    end
else
    error('mode should be: centered, forward, or backward')
end

end