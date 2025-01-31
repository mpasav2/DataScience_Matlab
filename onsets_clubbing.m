S02 = S01;
for i = 1:length(S02)
    S02(i).onset = 0;
end
for tt = 1:length(S02)
    for ttt = 1:length(S10)
        if (ismember(convertCharsToStrings(S10(ttt).condition),convertCharsToStrings(S02(tt).Exo))) && (ismember(convertCharsToStrings(S10(ttt).leg),convertCharsToStrings(S02(tt).Leg))) && (strcmp(S02(tt).Step, S10(ttt).Step)) && (S02(tt).onset == 0)
            loops = length(S10(ttt).onset);
            for j = 1:loops
                S02(tt+j-1).onset = S10(ttt).onset(j);
            end
        end
    end
    
end

clearvars -except S01 S10 S02