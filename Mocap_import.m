function S01 = Mocap_import(param, S01)

    trial = 1;
    for j = 1:length(S01)
        if strcmp(S01(j).Exo, 'base')
            trial = trial+1;
        end
    end
    cd(param.Motion_datafolder)
    files = dir;
    filenames = {};
    indx = 0;
    for i = 1:length(files)
        if files(i).bytes ~= 0 && files(i).isdir == 0 && (contains(files(i).name, '.mat') == 1)
            indx = indx + 1;
            filenames{indx} = files(i).name;
        end
    
    end

    for ee=1:length(param.exo)
        for ll=1:length(param.leg)
            for ss=1:length(param.step)            
                cond = length(param.leg)*length(param.step)*(ee-1) + length(param.step)*(ll-1) + ss + 1;
                files= param.filenum{cond};
                for ff=1:length(files)
                    datafile = filenames{files(ff)};
                        load(datafile)
                        if contains(datafile, param.exo{ee}) && contains(datafile, param.step{ss})
                            S01(trial).Exo = param.exo{ee};
                            S01(trial).Leg = param.leg{ll};
                            S01(trial).Step = param.step{ss};
                            for vv=1:length(param.MotVars)
                                if strcmp(param.MotVars{vv},'RHipMoment') || strcmp(param.MotVars{vv},'LHipMoment') 
                                    eval(['S01(trial).' param.MotVars{vv} '= transpose(' param.MotVars{vv} '{1,1}(3:end-3,:));']);
                                else
                                    eval(['S01(trial).' param.MotVars{vv} '= transpose(' param.MotVars{vv} '{1,1});']);
                                end
                                
    
                                %eval(['clear ' param.MotVars{vv}]);
                            end
                            if strcmp(S01(trial).Exo, 'exo')
                                try(LB1)
                                    S01(trial).LB1 = transpose(LB1{1,1});
                                catch
                                    S01(trial).LB1 = transpose(L_ASIS{1,1});
                                end
    
                                try(LB2)
                                    S01(trial).LB2 = transpose(LB2{1,1});
                                catch
                                    S01(trial).LB2 = transpose(L_Knee_Med{1,1});
                                end
    
                                try(RB1)
                                    S01(trial).RB1 = transpose(RB1{1,1});
                                catch
                                    S01(trial).RB1 = transpose(R_ASIS{1,1});
                                end
    
                                try(LB2)
                                    S01(trial).RB2 = transpose(RB2{1,1});
                                catch
                                    S01(trial).RB2 = transpose(R_Knee_Med{1,1});
                                end
    
                                S01(trial).RFT = transpose(RFT{1,1});
                                S01(trial).LFT = transpose(LFT{1,1});
                            end
                            trial = trial + 1;
                        end
                end
            end
        end
    end
    
    cd .. % move back to parent folder

end