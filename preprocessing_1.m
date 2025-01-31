function output_data = preprocessing_1(param, vari, S0x, opt)

output_data = {};
if contains(vari, 'EMG')
    case_num = 1;
else
    case_num = 2;
end
trial = 0;
for ee = 1:length(param.exo)
    for ll = 1:length(param.leg)
        for ss = 1:length(param.step)
            cond = length(param.leg) * length(param.step) * (ee-1) + length(param.step) * (ll-1) + ss + 1;
            files = param.filenum{cond};

            output_data(cond-1).Exo = param.exo{ee};
            output_data(cond-1).Leg = param.leg{ll};
            output_data(cond-1).Step = param.step{ss};
            trial = trial + length(param.trialnum{cond-1});

            data = [];
            for ff = 1:length(files)
                if case_num == 1
                    processed_data = filter_segment(eval(['S0x(trial+ff).' vari]), S0x(trial+ff).onset, 300, 300, S0x(trial+ff).time, case_num, opt);
                elseif case_num == 2
                    if contains(vari, 'LHipMoment') && strcmp(S0x(trial+ff).Exo, 'exo')

                        loadcell = (5.0008 * (S0x(trial+ff).LLoadCell) + 0.0028211) * 4.44822161533;

                        LB1 = S0x(trial+ff).LB1(1:3, :);
                        LB2 = S0x(trial+ff).LB2(1:3, :);
                        LFT = S0x(trial+ff).LFT(1:3, :);

                        % Calculate the displacement vector from RB1 to RB2
                        torque = cross(LFT, (LB1 - LB2) .* loadcell);
                        HipMoment = [1; 1; 1] .* (S0x(trial+ff).LHipMoment);
                        torque_truncated = torque(1:size(HipMoment, 1), :);
                        LTorque = [-1; 1; 1] .*(torque_truncated(:, 3:end-3) + HipMoment);
                        indx = vari(end-3);
                        LTorque = eval(['LTorque(' indx ',:)']);

                        clearvars loadcell LB1 LB2 LFT torque HipMoment torque_truncated Total_torque

                        processed_data = filter_segment((LTorque), S0x(trial+ff).onset, 300, 300, 0 + (2/100):1/100:10 - (2/100), case_num, opt);

                    elseif contains(vari, 'RHipMoment') && strcmp(S0x(trial+ff).Exo, 'exo')
                        loadcell = (4.9984 * (S0x(trial + ff).RLoadCell) + 0.0062595) * 4.44822161533;

                        RB1 = S0x(trial+ff).RB1(1:3, :);
                        RB2 = S0x(trial+ff).RB2(1:3, :);
                        RFT = S0x(trial+ff).RFT(1:3, :);

                        % Calculate the displacement vector from RB1 to RB2
                        torque = cross(RFT, (RB1 - RB2) .* loadcell);
                        HipMoment = [1; 1; 1] .* (S0x(trial+ff).RHipMoment);
                        torque_truncated = torque(1:size(HipMoment, 1), :);
                        RTorque = [-1; 1; 1] .*(torque_truncated(:, 3:end-3) + HipMoment);
                        indx = vari(end-3);
                        RTorque = eval(['RTorque(' indx ',:)']);

                        clearvars loadcell RB1 RB2 RFT torque HipMoment torque_truncated

                        processed_data = filter_segment(RTorque, S0x(trial+ff).onset, 300, 300, 0 + (2/100):1/100:10 - (2/100), case_num, opt);

                    elseif contains(vari, 'LHipMoment') || contains(vari, 'LHipAngle')
                        indx = vari(end-3);
                        if strcmp(S0x(trial+ff).Exo, 'noexo') || strcmp(S0x(trial+ff).Exo, 'exo')
                            if contains(vari, 'RHipAngle')
                                processed_data = filter_segment(eval(['-1 * S0x(trial+ff).' vari]), S0x(trial+ff).onset, 300, 300, 0 + (2/100):1/100:10 - (2/100), case_num, opt);
                            else
                                processed_data = filter_segment(eval(['-1*S0x(trial+ff).' vari]), S0x(trial+ff).onset, 300, 300, 0 + (2/100):1/100:10 - (2/100), case_num, opt);
                            end
                        else
                            processed_data = filter_segment(eval(['-1*S0x(trial+ff).' vari]), S0x(trial+ff).onset, 300, 300, 0:1/100:10, case_num, opt);
                        end
                    else
                        processed_data = filter_segment(eval(['-1*S0x(trial+ff).' vari]), S0x(trial+ff).onset, 300, 300, 0:1/100:10, case_num, opt);
                    end
                end

                if ff ~= 1
                    if length(processed_data) ~= size(data, 2)
                        if length(processed_data) < size(data, 2)
                            processed_data = [processed_data; zeros(size(data, 2) - length(processed_data), 1)];
                        else
                            processed_data = processed_data(1:size(data, 2));
                        end
                    end
                end

                % if files(ff) ~= 6 % Keep it only for S04
                   data = [data; processed_data'];
                %end

                output_data(cond-1).raw_data = data;
                output_data(cond-1).avg_data = mean(data, 1);
                output_data(cond-1).std_data = std(data, 1);
            end
        end
    end
end

if case_num == 1
    max_channel = 0;
    for i = 1:length(output_data)
        max_cond = max(max(output_data(i).raw_data));
        if max_cond > max_channel
            max_channel = max_cond;
        end
    end
    for i = 1:length(output_data)
        output_data(i).normalized_data = output_data(i).raw_data ./ max_channel;
        output_data(i).norm_avg_data = mean(output_data(i).normalized_data, 1);
        output_data(i).norm_std_data = std(output_data(i).normalized_data, 1);
    end
end

end
