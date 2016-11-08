function [ ] = analyzeDuration( results, plotFigs, exportFigs )
%ANALYZEDURATION Performs analysis on correlates of duration of training in days

durations = [11, 10, 10, 10, 16, 13, 12, 19, 16, 14, 24, 22, 20, 18, 24, 10, 10, 17, 10, 14, 31, 9, 19, 10, 19, 31, 12, 15, 10, 25,...
    19, 16, 11, 10, 18, 14, 17, 10, 27, 12, 14, 11, 19, 35, 12, 15, 13, 13, 10, 16, 11, 16, 14, 12, 12, 23, 10, 13, 12, 17, 10, 12]; 

groups = {'group_1', 'group_2'};
dom = {'perception', 'memory'};
stim = {'trained', 'untrained'};
for g = 1:numel(groups)
    trainingDuration.(groups{g}) = [];
    for d = 1:numel(dom)
        for s = 1:numel(stim)
            M_ratio_logDiff.(groups{g}).(dom{d}).(stim{s}) = [];
            AUC_diff.(groups{g}).(dom{d}).(stim{s}) = [];
            meanConf_diff.(groups{g}).(dom{d}).(stim{s}) = [];
        end
    end
end

subjects = fieldnames(results);
for sub = 1:numel(subjects)
    group = sprintf('group_%d', results.(subjects{sub}).group);
    results.(subjects{sub}).duration = durations(sub);
    trainingDuration.(group) = vertcat(trainingDuration.(group), results.(subjects{sub}).duration);
    for d = 1:numel(dom)
        for s = 1:numel(stim)
            M_ratio_logDiff.(subjects{sub}).(dom{d}).(stim{s}) = log(results.(subjects{sub}).session_10.(dom{d}).(stim{s}).fit.M_ratio) - log(results.(subjects{sub}).session_01.(dom{d}).(stim{s}).fit.M_ratio);
            M_ratio_logDiff.(group).(dom{d}).(stim{s}) = vertcat(M_ratio_logDiff.(group).(dom{d}).(stim{s}), M_ratio_logDiff.(subjects{sub}).(dom{d}).(stim{s}));
            AUC_diff.(subjects{sub}).(dom{d}).(stim{s}) = results.(subjects{sub}).session_10.(dom{d}).(stim{s}).AUC - results.(subjects{sub}).session_01.(dom{d}).(stim{s}).AUC;
            AUC_diff.(group).(dom{d}).(stim{s}) = vertcat(AUC_diff.(group).(dom{d}).(stim{s}), AUC_diff.(subjects{sub}).(dom{d}).(stim{s}));
            meanConf_diff.(subjects{sub}).(dom{d}).(stim{s}) = results.(subjects{sub}).session_10.(dom{d}).(stim{s}).meanConf - results.(subjects{sub}).session_01.(dom{d}).(stim{s}).meanConf;
            meanConf_diff.(group).(dom{d}).(stim{s}) = vertcat(meanConf_diff.(group).(dom{d}).(stim{s}), meanConf_diff.(subjects{sub}).(dom{d}).(stim{s}));
        end
    end
end

if plotFigs
    MRatioCorrelationScatter = figure;
    plotCounter = 0;
    for d = 1:numel(dom)
        for s = 1:numel(stim)
            plotCounter = plotCounter+1;
            subplot(2,2,plotCounter);
            set(gca,'fontsize', 14);
            for sub = 1:numel(subjects)
                if results.(subjects{sub}).group == 1
                   hPoint(1) = plot(results.(subjects{sub}).duration, M_ratio_logDiff.(subjects{sub}).(dom{d}).(stim{s}), 'o', 'markeredgecolor', 'k', 'markersize', 10); hold on;
                elseif results.(subjects{sub}).group == 2
                   hPoint(2) = plot(results.(subjects{sub}).duration, M_ratio_logDiff.(subjects{sub}).(dom{d}).(stim{s}), 'd', 'markeredgecolor', 'k', 'markersize', 10); hold on;
                end
                if strcmp(dom{d},'perception') && strcmp(stim{s}, 'trained')
                    set(hPoint(results.(subjects{sub}).group), 'markerfacecolor', [.5 0 0]);
                    title('Perception, Trained');
                elseif strcmp(dom{d},'perception') && strcmp(stim{s}, 'untrained')
                    set(hPoint(results.(subjects{sub}).group), 'markerfacecolor', [1 0 0]);
                    title('Perception, Untrained');
                elseif strcmp(dom{d},'memory') && strcmp(stim{s}, 'trained')
                    set(hPoint(results.(subjects{sub}).group), 'markerfacecolor', [0 0 .5]);
                    title('Memory, Trained');
                elseif strcmp(dom{d},'memory') && strcmp(stim{s}, 'untrained')
                    set(hPoint(results.(subjects{sub}).group), 'markerfacecolor', [0 0 1]);
                    title('Memory, Untrained');
                end
            end
            plot([0 40], [0 0], '-k', 'linewidth', 1.5);
            % Linear Regression Lines
                x = trainingDuration.group_1;
                X = [ones(size(trainingDuration.group_1)), x];
                y = M_ratio_logDiff.group_1.(dom{d}).(stim{s});
                b = X\y;
                yCalc = X*b;
                plot(x, yCalc, '--k', 'linewidth', 2);
%                 [rho, p] = corr(x,y)
                x2 = trainingDuration.group_2;
                X2 = [ones(size(trainingDuration.group_2)), x2];
                y2 = M_ratio_logDiff.group_2.(dom{d}).(stim{s});
                b2 = X2\y2;
                yCalc2 = X2*b2;
                plot(x2, yCalc2, ':k', 'linewidth', 2);
%                 [rho2, p2] = corr(x2,real(y2))
            if d == 2
                xlabel('Training Duration (days)');
            end
            if s == 1
                ylabel('log(meta-d''/d'') diff');
            end
            xlim([0 40]); ylim([-4 4]);
            leg = legend(hPoint, 'Control Group', 'Experimental Group', 'location', 'se');
            set(leg, 'fontsize', 6);
            legend boxoff;
            box off;
        end
    end
    if exportFigs
        export_fig MRatioCorrelationScatter -png -transparent 'MRatioDurationCorrelation.png';
    end
    
    AUCCorrelationScatter = figure;
    plotCounter = 0;
    for d = 1:numel(dom)
        for s = 1:numel(stim)
            plotCounter = plotCounter+1;
            subplot(2,2,plotCounter);
            set(gca,'fontsize', 14);
            for sub = 1:numel(subjects)
                if results.(subjects{sub}).group == 1
                   hPoint(1) = plot(results.(subjects{sub}).duration, AUC_diff.(subjects{sub}).(dom{d}).(stim{s}), 'o', 'markeredgecolor', 'k', 'markersize', 10); hold on;
                elseif results.(subjects{sub}).group == 2
                   hPoint(2) = plot(results.(subjects{sub}).duration, AUC_diff.(subjects{sub}).(dom{d}).(stim{s}), 'd', 'markeredgecolor', 'k', 'markersize', 10); hold on;
                end
                if strcmp(dom{d},'perception') && strcmp(stim{s}, 'trained')
                    set(hPoint(results.(subjects{sub}).group), 'markerfacecolor', [.5 0 0]);
                    title('Perception, Trained');
                elseif strcmp(dom{d},'perception') && strcmp(stim{s}, 'untrained')
                    set(hPoint(results.(subjects{sub}).group), 'markerfacecolor', [1 0 0]);
                    title('Perception, Untrained');
                elseif strcmp(dom{d},'memory') && strcmp(stim{s}, 'trained')
                    set(hPoint(results.(subjects{sub}).group), 'markerfacecolor', [0 0 .5]);
                    title('Memory, Trained');
                elseif strcmp(dom{d},'memory') && strcmp(stim{s}, 'untrained')
                    set(hPoint(results.(subjects{sub}).group), 'markerfacecolor', [0 0 1]);
                    title('Memory, Untrained');
                end
            end
            plot([0 40], [0 0], '-k', 'linewidth', 1.5);
            % Linear Regression Lines
                x = trainingDuration.group_1;
                X = [ones(size(trainingDuration.group_1)), x];
                y = AUC_diff.group_1.(dom{d}).(stim{s});
                b = X\y;
                yCalc = X*b;
                plot(x, yCalc, '--k', 'linewidth', 2);
%                 [rho, p] = corr(x,y)
                x2 = trainingDuration.group_2;
                X2 = [ones(size(trainingDuration.group_2)), x2];
                y2 = AUC_diff.group_2.(dom{d}).(stim{s});
                b2 = X2\y2;
                yCalc2 = X2*b2;
                plot(x2, yCalc2, ':k', 'linewidth', 2);
%                 [rho2, p2] = corr(x2,real(y2))
            
            if d == 2
                xlabel('Training Duration (days)');
            end
            if s == 1
                ylabel('AUC diff');
            end
            xlim([0 40]); ylim([-.5 .5]);
            leg = legend(hPoint, 'Control Group', 'Experimental Group', 'location', 'se');
            set(leg, 'fontsize', 6);
            legend boxoff;
            box off;
        end
    end
    if exportFigs
        export_fig AUCCorrelationScatter -png -transparent 'AUCDurationCorrelation.png';
    end
    
    meanConfCorrelationScatter = figure;
    plotCounter = 0;
    for d = 1:numel(dom)
        for s = 1:numel(stim)
            plotCounter = plotCounter+1;
            subplot(2,2,plotCounter);
            set(gca,'fontsize', 14);
            for sub = 1:numel(subjects)
                if results.(subjects{sub}).group == 1
                   hPoint(1) = plot(results.(subjects{sub}).duration, meanConf_diff.(subjects{sub}).(dom{d}).(stim{s}), 'o', 'markeredgecolor', 'k', 'markersize', 10); hold on;
                elseif results.(subjects{sub}).group == 2
                   hPoint(2) = plot(results.(subjects{sub}).duration, meanConf_diff.(subjects{sub}).(dom{d}).(stim{s}), 'd', 'markeredgecolor', 'k', 'markersize', 10); hold on;
                end
                if strcmp(dom{d},'perception') && strcmp(stim{s}, 'trained')
                    set(hPoint(results.(subjects{sub}).group), 'markerfacecolor', [.5 0 0]);
                    title('Perception, Trained');
                elseif strcmp(dom{d},'perception') && strcmp(stim{s}, 'untrained')
                    set(hPoint(results.(subjects{sub}).group), 'markerfacecolor', [1 0 0]);
                    title('Perception, Untrained');
                elseif strcmp(dom{d},'memory') && strcmp(stim{s}, 'trained')
                    set(hPoint(results.(subjects{sub}).group), 'markerfacecolor', [0 0 .5]);
                    title('Memory, Trained');
                elseif strcmp(dom{d},'memory') && strcmp(stim{s}, 'untrained')
                    set(hPoint(results.(subjects{sub}).group), 'markerfacecolor', [0 0 1]);
                    title('Memory, Untrained');
                end
            end
            plot([0 40], [0 0], '-k', 'linewidth', 1.5);
            % Linear Regression Lines
                x = trainingDuration.group_1;
                X = [ones(size(trainingDuration.group_1)), x];
                y = meanConf_diff.group_1.(dom{d}).(stim{s});
                b = X\y;
                yCalc = X*b;
                plot(x, yCalc, '--k', 'linewidth', 2);
%                 [rho, p] = corr(x,y)
                x2 = trainingDuration.group_2;
                X2 = [ones(size(trainingDuration.group_2)), x2];
                y2 = meanConf_diff.group_2.(dom{d}).(stim{s});
                b2 = X2\y2;
                yCalc2 = X2*b2;
                plot(x2, yCalc2, ':k', 'linewidth', 2);
%                 [rho2, p2] = corr(x2,real(y2))
            if d == 2
                xlabel('Training Duration (days)');
            end
            if s == 1
                ylabel('Mean Conf diff');
            end
            xlim([0 40]); ylim([-2 2]);
            leg = legend(hPoint, 'Control Group', 'Experimental Group', 'location', 'se');
            set(leg, 'fontsize', 6);
            legend boxoff;
            box off;
        end
    end
    if exportFigs
        export_fig meanConfCorrelationScatter -png -transparent 'meanConfDurationCorrelation.png';
    end
end

