function ESN_Moving_Dots
total_run_time = 0.5; % seconds
num_dots = 100; % number, digits
coherency = 25; % percent 0-100%
movement_direction = 0; % degrees 0-360
movement_speed = 7; % degrees/sec
size = 5; % degrees
x_offset = 1; % degrees
y_offset = 1; % degrees

% dependent variables
num_structured_moving = round(num_dots * coherency / 100);
num_random_moving = num_dots - num_structured_moving;

half_size = size ./ 2;

% Init positions
direction_structured_moving = ones(num_structured_moving, 1) .* movement_direction;
direction_random_moving = randi(360, num_random_moving, 1);
x_increment_unity_structured = cosd(direction_structured_moving);
y_increment_unity_structured = sind(direction_structured_moving);
x_increment_unity_random = cosd(direction_random_moving);
y_increment_unity_random = sind(direction_random_moving);

x_values_random_moving     = (rand(num_random_moving, 1).*size)-(half_size);
y_values_random_moving     = (rand(num_random_moving, 1).*size)-(half_size);
x_values_structured_moving = (rand(num_structured_moving, 1).*size)-(half_size);
y_values_structured_moving = (rand(num_structured_moving, 1).*size)-(half_size);

timerVal_overall_time = tic;
timerVal_plot = tic;
hfig_ = figure(1);
clf(hfig_)
while_flag = true;
while while_flag
    t_elapsted = toc(timerVal_plot);
    timerVal_plot = tic;
    increment_ = movement_speed * t_elapsted;
    x_increment_ = x_increment_unity_structured * increment_;
    y_increment_ = y_increment_unity_structured * increment_;
    x_values_structured_moving = x_values_structured_moving + x_increment_;
    y_values_structured_moving = y_values_structured_moving + y_increment_;
    x_values_structured_moving(x_values_structured_moving<-half_size) = x_values_structured_moving(x_values_structured_moving<-half_size) + size;
    x_values_structured_moving(x_values_structured_moving> half_size) = x_values_structured_moving(x_values_structured_moving> half_size) - size;
    y_values_structured_moving(y_values_structured_moving<-half_size) = y_values_structured_moving(y_values_structured_moving<-half_size) + size;
    y_values_structured_moving(y_values_structured_moving> half_size) = y_values_structured_moving(y_values_structured_moving> half_size) - size;
    x_increment_ = x_increment_unity_random * increment_;
    y_increment_ = y_increment_unity_random * increment_;
    x_values_random_moving = x_values_random_moving + x_increment_;
    y_values_random_moving = y_values_random_moving + y_increment_;
    x_values_random_moving(x_values_random_moving<-half_size) = x_values_random_moving(x_values_random_moving<-half_size) + size;
    x_values_random_moving(x_values_random_moving> half_size) = x_values_random_moving(x_values_random_moving> half_size) - size;
    y_values_random_moving(y_values_random_moving<-half_size) = y_values_random_moving(y_values_random_moving<-half_size) + size;
    y_values_random_moving(y_values_random_moving> half_size) = y_values_random_moving(y_values_random_moving> half_size) - size;
    
    figure(hfig_)
    plot([x_values_structured_moving; x_values_random_moving]+x_offset,...
         [y_values_structured_moving; y_values_random_moving]+y_offset, '.w', 'MarkerSize', 20)
    set(gca, 'Color', [0 0 0]);
    xlim([-half_size half_size]+x_offset)
    ylim([-half_size half_size]+y_offset)
    
    while_flag = toc(timerVal_overall_time) < total_run_time;
end



end