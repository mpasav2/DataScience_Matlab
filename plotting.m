function plotting(param, data, vari, opt)

    if contains(vari, 'EMG')
        case_num = 1;
    else
        case_num = 2;
    end
    
    colors = {[0 84/255 147/255],[83/255 27/255 147/255],[255/255 147/255 0], [148/255 17/255 0]};
    
    switch opt
        case 'avg'
            figure;
            for ee = 1:length(param.exo)
                subplot(2, 2, ee)
                for ll = 1:length(param.leg)
                    for ss = 1:length(param.step)            
                        for tt = 1:length(data)
                            if strcmp(param.exo{ee}, data(tt).Exo) && strcmp(param.leg{ll}, data(tt).Leg) && strcmp(param.step{ss}, data(tt).Step)
                                if case_num == 1
                                    time = linspace(-3, 3.0, length(data(tt).norm_avg_data));
                                    plot(time, data(tt).norm_avg_data, 'Color', colors{tt}, 'LineWidth', 3.5);
                                    hold on;
                                    fill([time, fliplr(time)], [data(tt).norm_avg_data + data(tt).norm_std_data, fliplr(data(tt).norm_avg_data - data(tt).norm_std_data)], colors{tt}, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
                                else
                                    time = linspace(-3, 3.0, length(data(tt).avg_data));
                                    plot(time, data(tt).avg_data, 'Color', colors{tt}, 'LineWidth', 3.5);
                                    hold on;
                                    fill([time, fliplr(time)], [data(tt).avg_data + data(tt).std_data, fliplr(data(tt).avg_data - data(tt).std_data)], colors{tt}, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
                                end
                            end
                        end
                    end
                end
                
                xlabel('Time (s)');
               
                ylim([-80 80]);
                yticks(-80:40:80);
                %ylim([-60 40]);
                %yticks(-60:20:40);
                ax = gca;
                ax.XColor = 'none';
                ax.Box = 'off';
                ax.YAxis.LineWidth = 3;
                ax.YAxis.FontName = 'Times New Roman';
                ax.YAxis.FontSize = 14;
                hold off;
                
            end
    
        case 'ind'
            figure;
            for ee = 1:length(param.exo)
                subplot(1, 2, ee)
                for ll = 1:length(param.leg)
                    for ss = 1:length(param.step)            
                        for tt = 1:length(data)
                            if strcmp(param.exo{ee}, data(tt).Exo) && strcmp(param.leg{ll}, data(tt).Leg) && strcmp(param.step{ss}, data(tt).Step)
                                for pp = 1:size(data(tt).raw_data, 1)
                                    if case_num == 1
                                        time = linspace(-3, 3.0, length(data(tt).normalized_data(pp, :)));
                                        plot(time, data(tt).normalized_data(pp, :), 'Color', colors{tt}, 'LineWidth', 3.5);
                                        hold on;
                                    else
                                        time = linspace(-3, 3.0, length(data(tt).raw_data(pp, :)));
                                        plot(time, data(tt).raw_data(pp, :), 'Color', colors{tt}, 'LineWidth', 3.5);
                                        hold on;
                                    end
                                end
                            end
                        end
                    end
                end
                xlabel('Time (s)');
                ax = gca;
                ax.XColor = 'none';
                ax.Box = 'off';
                ax.YAxis.LineWidth = 3;
                ax.YAxis.FontName = 'Times New Roman';
                ax.YAxis.FontSize = 14;
                hold off;
            end
    end
end
