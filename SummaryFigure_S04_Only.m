%% Combine data
clear; clc;
% Works from the '.../All_subjects_ext_abd/' folder
cd('C:\Users\mpasav2\OneDrive - University of Illinois Chicago\Desktop\All_particpants_GI_study\All Subjects_Data')

sublist = {'S01','S02','S03','S04','S05','C01','C02'};
%sublist={'S01'};
planelist = {'frontal','sagittal'};
exolist = {'noexo','exo'};
leglist = {'right','left','right','right','right','right','right'};
%leglist = {'right'};

for ss = 1:length(sublist)
    cd(sublist{ss})
    for pp = 1:length(planelist)
        for ee = 1:length(exolist)
            filename = [];
            filename = [sublist{ss}, 'baseline_', planelist{pp}, '_max_normal_', exolist{ee}, '_', leglist{ss},'.mat'];
            load(filename)
            Data(2*ee-1).Exo = exolist{ee};
            Data(2*ee).Exo = exolist{ee};
            Data(2*ee-1).Step = 'Normal';
            Data(2*ee).Step = 'Long';
            if pp==1 % frontal plane torque
                Data(2*ee-1).HABD = hippos;%  - repmat(mean(hiposs(:,1:100),2),1,size(hiposs,2)); % uncomment to baseline-correct
                Data(2*ee).HABD = hipos;% - repmat(mean(hipposs(:,1:100),2),1,size(hipposs,2));
            elseif pp==2    % sagittal plane torque
                Data(2*ee-1).HEF = hippos;% - repmat(mean(hiposs(:,1:100),2),1,size(hiposs,2));
                Data(2*ee).HEF = hipos;% - repmat(mean(hipposs(:,1:100),2),1,size(hipposs,2));
            end
            eval([sublist{ss}, '= Data;'])
        end
    end
    cd ..
end

% necessary variables saved as SummaryResults.mat
%% Plot
%clear
%load SummaryResults
clc; close all;
% Color scheme
Colors = [ 0, 84/255, 147/255; ...    % noexo normal
    83/255, 27/255, 147/255; ...    % noexo long
    255/255, 147/255, 0; ...    % exo normal
    148/255, 17/255, 0];       % exo long

maxHE = [92.57 74.74 77.88 79.45 28.82 166.03 84.55]; % Nm for {'S01','S02','S03','S04','C05','C06'};
fliporder = [4 3 2 1];
plotornot = 1; % Put 1 for plotting individual plots

timelength = 310;% S01=350,S02=350, S03=400, S04=350,S05=380, C01=400
for ss=6 % Choose subject to plot, from sublist = {'S01','S02','S03','S04','C05','C06'};
    figure;
    eval(['Data = ' sublist{ss} ';'])
    accumHABDpeak = [];
    accumHEFpeak = [];

    for cc=1:4
        HABDpeak_cond = [];
        HEFpeak_cond = [];
        for tt = 1:size(Data(cc).HABD)

            if plotornot == 1
                figure(10)
                subplot(2,1,1)
                plot(Data(cc).HABD(tt,:))
                hold on
                [~,p_indx]= min(Data(cc).HABD(tt,1:timelength));
                scatter(p_indx,min(Data(cc).HABD(tt,1:timelength)),'green','filled')
                hold off
                drawnow

                subplot(2,1,2)
                plot(Data(cc).HEF(tt,:))
                hold on
                [~,p_indx]= min(Data(cc).HEF(tt,1:timelength)); %Change sign
                scatter(p_indx,min(Data(cc).HEF(tt,1:timelength)),'red','filled') %Change sign
                hold off
                drawnow

                
                        end % Put breakpoint for looking at individual plots
    
            close(figure(10))

            accumHABDpeak = [accumHABDpeak, min((Data(cc).HABD(tt,1:timelength)'))];
            HABDpeak_cond = [HABDpeak_cond, min((Data(cc).HABD(tt,1:timelength)'))];
            
    
            accumHEFpeak = [accumHEFpeak, min((Data(cc).HEF(tt,1:timelength)'))]; %Change sign
            HEFpeak_cond = [HEFpeak_cond, min((Data(cc).HEF(tt,1:timelength)'))]; %Change sign
    
            subplot(2,2,2)
            scatter(abs(min(Data(cc).HEF(tt,1:timelength))'), abs(max(Data(cc).HABD(tt,1:timelength))'),250,Colors(cc,:),'.')   % peaks from individual trial %Change sign`
            hold on
           
        end
            subplot(2,2,1)
            bar(cc,mean(abs(HABDpeak_cond)),'EdgeColor','none','FaceColor',Colors(cc,:))
            hold on
            errorbar(cc, mean(abs(HABDpeak_cond)), std(abs(HABDpeak_cond)), ...
                'k','LineWidth',2,'CapSize',10)
            if cc==4
                box off
                ubhabd = abs(min(accumHABDpeak))*1.1;
                %xlim([0.5 4.5]); 
                ylim([0 ubhabd ])
                ylabel('HABD Torque (Nm)')
                xticks([1 2 3 4])
                xticklabels({'Norm','Long','Norm','Long'});
                fontsize(gca,16,"points")
            end

            subplot(2,2,4)
            barh(fliporder(cc),mean(abs(HEFpeak_cond)),'EdgeColor','none','FaceColor',Colors(cc,:))%abs hip ext 
            hold on
            errorbar(mean(abs(HEFpeak_cond)), fliporder(cc), std(abs(HEFpeak_cond)), "horizontal", ...%abs hip ext
                'k','LineWidth',2,'CapSize',10)
            if cc==4
                box off
                plot([maxHE(ss)*0.25 maxHE(ss)*0.25], ylim, ':', 'Color', [0.5 0.5 0.5],'LineWidth',2.5) % 25% MVT line
                ubhef = min(accumHEFpeak)*1.1;
                %ylim([0.5 4.5]); 
                xlim([0 abs(ubhef)])
                xlabel('HE Torque (Nm)')
                yticks([1 2 3 4])
                yticklabels({'Long','Norm','Long','Norm'});
                fontsize(gca,16,"points")
                ax.YAxis.FontName = 'Arial';
            end

             if cc==4 && tt == size(Data(cc).HABD,1)
                subplot(2,2,2)
                plot([maxHE(ss)*0.25 maxHE(ss)*0.25], [0 ubhabd], ':', 'Color', [0.5 0.5 0.5],'LineWidth',2.5)  % 25% MVT line
                text(maxHE(ss)*0.25, ubhabd*0.4, {'25%'; 'MVT'},'Color',[0.5 0.5 0.5],'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle','Rotation', -90)
                xlim([0 abs(ubhef)]); 
                ylim([0 abs(ubhabd)])
    
                % correlation analysis
                [r, pval] = corrcoef(accumHEFpeak,accumHABDpeak);
                CorrResults(ss,1) = r(1,2);
                CorrResults(ss,2) = pval(1,2);
    
                title(['r=', num2str(r(1,2),2), ' (p=', num2str(pval(1,2),2), ')'])
    
                lm = fitlm(abs(accumHEFpeak),abs(accumHABDpeak)); % linear regression to illustrate the linear relationship
                predictrange = min(abs(accumHEFpeak)):5:max(abs(accumHEFpeak));
                ypredict = predict(lm,predictrange');
                plot(predictrange,ypredict,'k')
                fontsize(gca,16,"points")
                hold off
            end
        end
    end
