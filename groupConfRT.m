function [ confRT ] = groupConfRT( results, plotFigs, exportFigs )
%GROUPRT Runs group RT analysis and optionally plots and expoconfRTs the figures

dom = {'perception', 'memory'};
stim = {'trained', 'untrained'};
subjects = fieldnames(results);
% Initialize arrays
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    for sesh = 1:10
        session = sprintf('session_%.2d', sesh);
        if sesh == 1 || sesh == 10
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    confRT.(groups{g}).(session).(dom{d}).(stim{s}).raw = [];
                end
            end
        else  % sessions 2-9
            confRT.(groups{g}).(session).perception.trained.raw = [];
        end
    end
end

% Concatenate raw data
for sub = 1:numel(subjects)
    group = sprintf('group_%d', results.(subjects{sub}).group);
    sessions = fieldnames(results.(subjects{sub}));
    for sesh = 1:numel(sessions)
        if strncmp(sessions{sesh},'session',7)
            session = str2double(sessions{sesh}(end-1:end));
            if session == 1 || session == 10
                for d = 1:numel(dom)
                    for s = 1:numel(stim)
                        confRT.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw = vertcat(confRT.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).meanConfRT);
                    end
                end
            else % Sessions 2-9
                confRT.(group).(sessions{sesh}).perception.trained.raw = vertcat(confRT.(group).(sessions{sesh}).perception.trained.raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.meanConfRT);
            end
        end
    end    
end

% Take mean and standard error
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    [plots.confRT.(groups{g}).learningCurve.mean, plots.confRT.(groups{g}).learningCurve.sem] = deal([]);
    for sesh = 1:10
        session = sprintf('session_%.2d',sesh);
        if sesh == 1 || sesh == 10
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    confRT.(groups{g}).(session).(dom{d}).(stim{s}).mean = nanmean(confRT.(groups{g}).(session).(dom{d}).(stim{s}).raw);
                    confRT.(groups{g}).(session).(dom{d}).(stim{s}).sem = nanstd(confRT.(groups{g}).(session).(dom{d}).(stim{s}).raw)/sqrt(length(confRT.(groups{g}).(session).(dom{d}).(stim{s}).raw));
                end
            end
            plots.confRT.(groups{g}).(session).mean = [confRT.(groups{g}).(session).perception.trained.mean, confRT.(groups{g}).(session).perception.untrained.mean;...
                confRT.(groups{g}).(session).memory.trained.mean, confRT.(groups{g}).(session).memory.untrained.mean];
            plots.confRT.(groups{g}).(session).sem = [confRT.(groups{g}).(session).perception.trained.sem, confRT.(groups{g}).(session).perception.untrained.sem;...
                confRT.(groups{g}).(session).memory.trained.sem, confRT.(groups{g}).(session).memory.untrained.sem];
        else % sessions 2-9
            confRT.(groups{g}).(session).perception.trained.mean = nanmean(confRT.(groups{g}).(session).perception.trained.raw);
            confRT.(groups{g}).(session).perception.trained.sem = nanstd(confRT.(groups{g}).(session).perception.trained.raw)/sqrt(length(confRT.(groups{g}).(session).perception.trained.raw));
            plots.confRT.(groups{g}).learningCurve.mean = vertcat(plots.confRT.(groups{g}).learningCurve.mean, confRT.(groups{g}).(session).perception.trained.mean);
            plots.confRT.(groups{g}).learningCurve.sem = vertcat(plots.confRT.(groups{g}).learningCurve.sem, confRT.(groups{g}).(session).perception.trained.sem);
        end
    end
    
end

if plotFigs
    % Error Bar Comparison Plot
    confRTComparisonPlot = figure;
    set(gcf,'position', [200 200 450 300]);
    subplot(1,2,1); % Control Group
    [hBar hErrorbar] = barwitherr([[plots.confRT.group_1.session_01.sem(1,:), plots.confRT.group_1.session_01.sem(2,:)]',[plots.confRT.group_1.session_10.sem(1,:), plots.confRT.group_1.session_10.sem(2,:)]']',...
        [[plots.confRT.group_1.session_01.mean(1,:), plots.confRT.group_1.session_01.mean(2,:)]',[plots.confRT.group_1.session_10.mean(1,:), plots.confRT.group_1.session_10.mean(2,:)]']');
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    set(gca, 'fontsize', 14);
    ylabel('confRT', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 6); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    [hBar hErrorbar] = barwitherr([[plots.confRT.group_2.session_01.sem(1,:), plots.confRT.group_2.session_01.sem(2,:)]',[plots.confRT.group_2.session_10.sem(1,:), plots.confRT.group_2.session_10.sem(2,:)]']',...
        [[plots.confRT.group_2.session_01.mean(1,:), plots.confRT.group_2.session_01.mean(2,:)]',[plots.confRT.group_2.session_10.mean(1,:), plots.confRT.group_2.session_10.mean(2,:)]']');
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    set(gca, 'fontsize', 14);
    ylabel('confRT', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Experimental Group', 'fontsize', 14);
    box off;
    if exportFigs
        export_fig confRTComparisonPlot -png -transparent 'confRTComparison.png';
    end
    
    confRTCurvePlot = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    hBar(1) = errorbar(.6, confRT.group_1.session_01.perception.trained.mean, confRT.group_1.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, confRT.group_1.session_01.perception.untrained.mean, confRT.group_1.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, confRT.group_1.session_01.memory.trained.mean, confRT.group_1.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, confRT.group_1.session_01.memory.untrained.mean, confRT.group_1.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.confRT.group_1.learningCurve.mean, plots.confRT.group_1.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, confRT.group_1.session_10.perception.trained.mean, confRT.group_1.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, confRT.group_1.session_10.perception.untrained.mean, confRT.group_1.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, confRT.group_1.session_10.memory.trained.mean, confRT.group_1.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, confRT.group_1.session_10.memory.untrained.mean, confRT.group_1.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('confRT', 'fontsize', 14);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'ne');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    hBar(1) = errorbar(.6, confRT.group_2.session_01.perception.trained.mean, confRT.group_2.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, confRT.group_2.session_01.perception.untrained.mean, confRT.group_2.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, confRT.group_2.session_01.memory.trained.mean, confRT.group_2.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, confRT.group_2.session_01.memory.untrained.mean, confRT.group_2.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.confRT.group_2.learningCurve.mean, plots.confRT.group_2.learningCurve.sem,'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, confRT.group_2.session_10.perception.trained.mean, confRT.group_2.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, confRT.group_2.session_10.perception.untrained.mean, confRT.group_2.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, confRT.group_2.session_10.memory.trained.mean, confRT.group_2.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, confRT.group_2.session_10.memory.untrained.mean, confRT.group_2.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('confRT', 'fontsize', 14);
    title('Experimental Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'ne');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    if exportFigs
        export_fig confRTCurvePlot -png -transparent 'confRTCurve.png';
    end

end
end