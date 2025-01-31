function S02 = onset_import(params,S02)

    cd(params.Onset_datafolder)

    files = dir;
    filenames = {};
    indx = 0;
    for i = 1:length(files)
        if files(i).bytes ~= 0 && files(i).isdir == 0 && (contains(files(i).name, '.mat') == 1)
            indx = indx + 1;
            filenames{indx} = files(i).name;
        end
    
    end


    for tt = 1:length(filenames) % We are only selecting files twice
        datafile = filenames{tt};
        load(datafile)
    
        if contains(datafile, '_exo')
            condition = 'exo';
        elseif contains(datafile, 'noexo')
            condition = 'noexo';
        end
    
        if contains(datafile, 'left')
            leg = 'left';
        elseif contains(datafile, 'right')
            leg = 'right';
        end
    
        S10(2*tt-1).condition = condition;
        S10(2*tt-1).leg = leg;
        S10(2*tt-1).Step = 'max';
        S10(2*tt-1).onset = onsets;
        
        S10(2*tt).condition = condition;
        S10(2*tt).leg = leg;
        S10(2*tt).Step = 'normal';
        S10(2*tt).onset = onsets1;
    end
    
    for i = 1:length(S02)
        S02(i).onset = 0;
    end
    for tt = 1:length(S02)
        for ttt = 1:length(S10)
            if (strcmp(S10(ttt).condition,S02(tt).Exo)) && (strcmp(S10(ttt).leg,S02(tt).Leg)) && (strcmp(S02(tt).Step, S10(ttt).Step)) && (S02(tt).onset == 0)
                loops = length(S10(ttt).onset);
                for j = 1:loops
                    S02(tt+j-1).onset = S10(ttt).onset(j);
                end
            end
        end
        
    end
    cd ..
end