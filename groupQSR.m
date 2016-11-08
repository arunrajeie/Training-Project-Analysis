function [ QSR ] = groupQSR( analysis, results, plotFigs, exportFigs )
%GROUPQSR Runs group QSR analysis and optionally plots and exports the figures

dom = {'perception', 'memory'};
stim = {'trained', 'untrained'};
subjects = fieldnames(results);
% Initialize arrays
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    for sesh = 1:10
        session = sprintf('session_%.2d', sesh);
        TxTQSR.(groups{g}).(session).raw = [];
        if sesh == 1 || sesh == 10
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    QSR.(groups{g}).(session).(dom{d}).(stim{s}).raw = [];
                end
            end
        else  % sessions 2-9
            QSR.(groups{g}).(session).perception.trained.raw = [];
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
                        TxTQSR.(group).(sessions{sesh}).raw = vertcat(TxTQSR.(group).(sessions{sesh}).raw, analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).QSR');
                        QSR.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw = vertcat(QSR.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).meanQSR);
                    end
                end
            else % Sessions 2-9
                TxTQSR.(group).(sessions{sesh}).raw = vertcat(TxTQSR.(group).(sessions{sesh}).raw, analysis.(subjects{sub}).(sessions{sesh}).perception.trained.QSR');
                QSR.(group).(sessions{sesh}).perception.trained.raw = vertcat(QSR.(group).(sessions{sesh}).perception.trained.raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.meanQSR);
            end
        end
    end    
end

% Take mean and standard error
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    [plots.QSR.(groups{g}).learningCurve.mean, plots.QSR.(groups{g}).learningCurve.sem] = deal([]);
    [plots.TxTQSR.(groups{g}).allSessions.mean, plots.TxTQSR.(groups{g}).allSessions.sem] = deal([]);
    for sesh = 1:10
        session = sprintf('session_%.2d',sesh);
        TxTQSR.(groups{g}).(session).mean = nanmean(TxTQSR.(groups{g}).(session).raw, 1);
        TxTQSR.(groups{g}).(session).sem = nanstd(TxTQSR.(groups{g}).(session).raw, [], 1)/sqrt(size(TxTQSR.(groups{g}).(session).raw,1));
        plots.TxTQSR.(groups{g}).allSessions.mean = horzcat(plots.TxTQSR.(groups{g}).allSessions.mean, TxTQSR.(groups{g}).(session).mean);
        plots.TxTQSR.(groups{g}).allSessions.sem = horzcat(plots.TxTQSR.(groups{g}).allSessions.sem, TxTQSR.(groups{g}).(session).sem);
        if sesh == 1 || sesh == 10
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    QSR.(groups{g}).(session).(dom{d}).(stim{s}).mean = nanmean(QSR.(groups{g}).(session).(dom{d}).(stim{s}).raw);
                    QSR.(groups{g}).(session).(dom{d}).(stim{s}).sem = nanstd(QSR.(groups{g}).(session).(dom{d}).(stim{s}).raw)/sqrt(length(QSR.(groups{g}).(session).(dom{d}).(stim{s}).raw));
                end
            end
            plots.QSR.(groups{g}).(session).mean = [QSR.(groups{g}).(session).perception.trained.mean, QSR.(groups{g}).(session).perception.untrained.mean;...
                QSR.(groups{g}).(session).memory.trained.mean, QSR.(groups{g}).(session).memory.untrained.mean];
            plots.QSR.(groups{g}).(session).sem = [QSR.(groups{g}).(session).perception.trained.sem, QSR.(groups{g}).(session).perception.untrained.sem;...
                QSR.(groups{g}).(session).memory.trained.sem, QSR.(groups{g}).(session).memory.untrained.sem];
        else % sessions 2-9
            QSR.(groups{g}).(session).perception.trained.mean = nanmean(QSR.(groups{g}).(session).perception.trained.raw);
            QSR.(groups{g}).(session).perception.trained.sem = nanstd(QSR.(groups{g}).(session).perception.trained.raw)/sqrt(length(QSR.(groups{g}).(session).perception.trained.raw));
            plots.QSR.(groups{g}).learningCurve.mean = vertcat(plots.QSR.(groups{g}).learningCurve.mean, QSR.(groups{g}).(session).perception.trained.mean);
            plots.QSR.(groups{g}).learningCurve.sem = vertcat(plots.QSR.(groups{g}).learningCurve.sem, QSR.(groups{g}).(session).perception.trained.sem);
        end
    end
    
end

if plotFigs
    TxTQSRPlot = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    lineProps = struct('color', 'm', 'width', 2);
    mseb(1:numel(plots.TxTQSR.group_1.allSessions.mean), plots.TxTQSR.group_1.allSessions.mean, plots.TxTQSR.group_1.allSessions.sem, lineProps); hold on;
    for t = 108.5:270:2268.5
        plot([t t], [0 1], '-k', 'linewidth', 1);
    end
    xlim([1 numel(plots.TxTQSR.group_1.allSessions.mean)]); ylim([.5 1]);
    set(gca,'fontsize',14);
    title('Control Group');
    xlabel('Trials');
    ylabel('QSR');
    subplot(1,2,2); % Experimental Group
    mseb(1:numel(plots.TxTQSR.group_2.allSessions.mean), plots.TxTQSR.group_2.allSessions.mean, plots.TxTQSR.group_2.allSessions.sem, lineProps); hold on;
    for t = 108.5:270:2268.5
        plot([t t], [0 1], '-k', 'linewidth', 1);
    end
    xlim([1 numel(plots.TxTQSR.group_2.allSessions.mean)]); ylim([.5 1]);
    set(gca,'fontsize',14);
    title('Experimental Group');
    xlabel('Trials');
    ylabel('QSR');
    if exportFigs
        export_fig TxTQSRPlot -png -transparent 'TxTQSR.png';
    end
    
    % Error Bar Comparison Plot
    QSRComparisonPlot = figure;
    set(gcf,'position', [200 200 450 300]);
    subplot(1,2,1); % Control Group
    [hBar hErrorbar] = barwitherr([[plots.QSR.group_1.session_01.sem(1,:), plots.QSR.group_1.session_01.sem(2,:)]',[plots.QSR.group_1.session_10.sem(1,:), plots.QSR.group_1.session_10.sem(2,:)]']',...
        [[plots.QSR.group_1.session_01.mean(1,:), plots.QSR.group_1.session_01.mean(2,:)]',[plots.QSR.group_1.session_10.mean(1,:), plots.QSR.group_1.session_10.mean(2,:)]']');
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    ylim([.65 .95]);
    set(gca, 'fontsize', 14);
    ylabel('QSR', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 6); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    [hBar hErrorbar] = barwitherr([[plots.QSR.group_2.session_01.sem(1,:), plots.QSR.group_2.session_01.sem(2,:)]',[plots.QSR.group_2.session_10.sem(1,:), plots.QSR.group_2.session_10.sem(2,:)]']',...
        [[plots.QSR.group_2.session_01.mean(1,:), plots.QSR.group_2.session_01.mean(2,:)]',[plots.QSR.group_2.session_10.mean(1,:), plots.QSR.group_2.session_10.mean(2,:)]']');
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    ylim([.65 .95]);
    set(gca, 'fontsize', 14);
    ylabel('QSR', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Experimental Group', 'fontsize', 14);
    box off;
    if exportFigs
        export_fig QSRComparisonPlot -png -transparent 'QSRComparison.png';
    end
    
    QSRCurvePlot = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    hBar(1) = errorbar(.6, QSR.group_1.session_01.perception.trained.mean, QSR.group_1.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, QSR.group_1.session_01.perception.untrained.mean, QSR.group_1.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, QSR.group_1.session_01.memory.trained.mean, QSR.group_1.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, QSR.group_1.session_01.memory.untrained.mean, QSR.group_1.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.QSR.group_1.learningCurve.mean, plots.QSR.group_1.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, QSR.group_1.session_10.perception.trained.mean, QSR.group_1.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, QSR.group_1.session_10.perception.untrained.mean, QSR.group_1.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, QSR.group_1.session_10.memory.trained.mean, QSR.group_1.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, QSR.group_1.session_10.memory.untrained.mean, QSR.group_1.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    ylim([.65 .95]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('QSR', 'fontsize', 14);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    hBar(1) = errorbar(.6, QSR.group_2.session_01.perception.trained.mean, QSR.group_2.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, QSR.group_2.session_01.perception.untrained.mean, QSR.group_2.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, QSR.group_2.session_01.memory.trained.mean, QSR.group_2.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, QSR.group_2.session_01.memory.untrained.mean, QSR.group_2.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.QSR.group_2.learningCurve.mean, plots.QSR.group_2.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, QSR.group_2.session_10.perception.trained.mean, QSR.group_2.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, QSR.group_2.session_10.perception.untrained.mean, QSR.group_2.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, QSR.group_2.session_10.memory.trained.mean, QSR.group_2.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, QSR.group_2.session_10.memory.untrained.mean, QSR.group_2.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    ylim([.65 .95]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('QSR', 'fontsize', 14);
    title('Experimental Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    if exportFigs
        export_fig QSRCurvePlot -png -transparent 'QSRCurve.png';
    end
    
end
end