function segmented_channel = filter_segment(channel, onset, pre, post, time, case_num, opt)
 
switch case_num
    case 1

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
        
        onset_indx = find(round(time,3) == onset/100);
        pre_indx = find(round(time,3) == pre/100);
        post_indx = find(round(time,3) == post/100);
    
        % Segment the filtered channel
        segmented_channel = filtered_channel(onset_indx(1) - pre_indx(1):onset_indx(1) + post_indx(1));
        %segmented_channel = channel(onset_indx(1) - pre_indx(1):onset_indx(1) + post_indx(1));
        segmented_time = time(onset_indx(1) - pre_indx(1):onset_indx(1) + post_indx(1));
        
    case 2
        if opt == 1
            filtered_channel = channel - mean(channel(1:150)); % 2nd step
    
            % Low-pass filter parameters
            F_low = 6; % cutoff frequency
            
            fs = length(channel)/time(end);
            % Apply low-pass filter
            [bemg, aemg] = butter(4, F_low / (fs/2), 'low'); % 4th step
            filtered_channel = filtfilt(bemg, aemg, filtered_channel'); 
            
            onset_indx = find(round(time,3) == onset/100);
            pre_indx = find(round(time,3) == pre/100);
            post_indx = find(round(time,3) == post/100);
        
            % Segment the filtered channel
            segmented_channel = filtered_channel(onset_indx(1) - pre_indx(1):onset_indx(1) + post_indx(1));
            %segmented_channel = channel(onset_indx(1) - pre_indx(1):onset_indx(1) + post_indx(1));
            segmented_time = time(onset_indx(1) - pre_indx(1):onset_indx(1) + post_indx(1));
        else
            % Low-pass filter parameters
            F_low = 6; % cutoff frequency
            
            fs = length(channel)/time(end);
            % Apply low-pass filter
            [bemg, aemg] = butter(4, F_low / (fs/2), 'low'); % 4th step
            filtered_channel = filtfilt(bemg, aemg, channel'); 
            
            onset_indx = find(round(time,3) == onset/100);
            pre_indx = find(round(time,3) == pre/100);
            post_indx = find(round(time,3) == post/100);
        
            % Segment the filtered channel
            segmented_channel = filtered_channel(onset_indx(1) - pre_indx(1):onset_indx(1) + post_indx(1));
            %segmented_channel = channel(onset_indx(1) - pre_indx(1):onset_indx(1) + post_indx(1));
            segmented_time = time(onset_indx(1) - pre_indx(1):onset_indx(1) + post_indx(1));
        end

end
end
