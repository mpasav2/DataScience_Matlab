function p_emg = filtering(channel, torque, time, mode, plotting)
    
    % Subtract mean of initial portion
    filtered_channel = channel - mean(channel); % 2nd step

    % Apply absolute value
    filtered_channel = abs(filtered_channel); % 3rd step

    % Low-pass filter parameters
    F_low = 6; % cutoff frequency
    
    fs = length(channel)/time(end);
    % Apply low-pass filter
    [bemg, aemg] = butter(4, F_low / (fs/2), 'low'); % 4th step
    filtered_channel = filtfilt(bemg, aemg, filtered_channel'); 

    normalized_channel = filtered_channel./max(filtered_channel);
    
    [btorque, atorque] = butter(4, F_low / (100/2), 'low'); % 4th step
    filtered_torque = filtfilt(btorque, atorque, torque');
    time_torque = 0:1/100:10;
    time_torque = time_torque(3:end-3);

    if plotting == 1
        figure;
        subplot(3,1,1)
        plot(time,channel);
        xlim([0 10])
        subplot(3,1,2)
        plot(time,normalized_channel);
        xlim([0 10])
        subplot(3,1,3)
        plot(time_torque,filtered_torque);
        xlim([0 10])
    end
    if mode == 1
        t_min = min(filtered_torque);
        t_min_time = find(filtered_torque == t_min)/100;
        t_min_start = find(round(time,2) == round((t_min_time - 0.1),2));
        t_min_end = find(round(time,2) == t_min_time);
        p_emg = mean(normalized_channel(t_min_start(1):t_min_end(end)));
    elseif mode == 2
        t_max = max(filtered_torque);
        t_max_time = find(filtered_torque == t_max)/100;
        t_max_start = find(round(time,2) == round((t_max_time - 0.1),2));
        t_max_end = find(round(time,2) == t_max_time);
        p_emg = mean(normalized_channel(t_max_start(1):t_max_end(end)));
    end
end