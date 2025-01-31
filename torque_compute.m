function S0x = torque_compute(param, S0x)

    for ee=1:length(param.exo)
        for ll=1:length(param.leg)
            for ss=1:length(param.step)
                cond = length(param.leg)*length(param.step)*(ee-1) + length(param.step)*(ll-1) + ss + 1;
                files= param.trialnum{cond};
                for ff=1:length(files)
                    if strcmp(S0x(files(ff)).Exo, 'exo')
    
                        % Computing for left torques
                        loadcell = (5.0008*(S0x(files(ff)).LLoadCell) + 0.0028211)*4.44822161533;
        
                        LB1 = S0x(files(ff)).LB1(1:3, :);
                        LB2 = S0x(files(ff)).LB2(1:3, :);
                        LFT= S0x(files(ff)).LFT(1:3, :);
        
                        % Calculate the displacement vector from RB1 to RB2
                        torque = cross(LFT, (LB1 - LB2) .* loadcell);
                        HipMoment = [-1;1;1].*(S0x(files(ff)).LHipMoment);
                        torque_truncated = torque(1:size(HipMoment, 1), :);
                        LTorque = torque_truncated(:,3:end-3) + HipMoment;
        
                        clearvars loadcell LB1 LB2 LFT torque HipMoment torque_truncated Total_torque
        
                        % Computing for left torques
                        loadcell = (4.9984*(S0x(files(ff)).RLoadCell) + 0.0062595)*4.44822161533;
        
                        RB1 = S0x(files(ff)).RB1(1:3, :);
                        RB2 = S0x(files(ff)).RB2(1:3, :);
                        RFT= S0x(files(ff)).RFT(1:3, :);
        
                        % Calculate the displacement vector from RB1 to RB2
                        torque = cross(RFT, (RB1 - RB2) .* loadcell);
                        HipMoment = [-1;1;1].*(S0x(files(ff)).RHipMoment);
                        torque_truncated = torque(1:size(HipMoment, 1), :);
                        RTorque= torque_truncated(:,3:end-3) + HipMoment;
        
                        clearvars loadcell RB1 RB2 RFT torque HipMoment torque_truncated
        
                        S0x(files(ff)).LTorque = LTorque;
                        S0x(files(ff)).RTorque = RTorque;
                    end
                end
            end
        end
    end
end