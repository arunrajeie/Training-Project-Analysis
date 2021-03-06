function [ meanConf ] = groupMeanConf( analysis, results, plotFigs, exportFigs )
%GROUPMEANCONF Runs group meanConf analysis and optionally plots and exports the figures

dom = {'perception', 'memory'};
stim = {'trained', 'untrained'};
subjects = fieldnames(results);
% Initialize arrays
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    for sesh = 1:10
        session = sprintf('session_%.2d', sesh);
        confResp.(groups{g}).(session).raw = [];
        if sesh == 1 || sesh == 10
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    meanConf.(groups{g}).(session).(dom{d}).(stim{s}).raw = [];
                end
            end
        else  % sessions 2-9
            meanConf.(groups{g}).(session).perception.trained.raw = [];
        end
    end
end

% Concatenate raw data
for sub = 1:numel(subjects)
    group = sprintf('group_%d', results.(subjects{sub}).group);
    sessions = fieldnames(results.(subjects{sub}));
    meanConf.(group).(subjects{sub}).correct.raw = [];
    meanConf.(group).(subjects{sub}).incorrect.raw = [];
    for sesh = 1:numel(sessions)
        if strncmp(sessions{sesh},'session',7)
            session = str2double(sessions{sesh}(end-1:end));
            if session == 1 || session == 10
                for d = 1:numel(dom)
                    for s = 1:numel(stim)
                        confResp.(group).(sessions{sesh}).raw = vertcat(confResp.(group).(sessions{sesh}).raw, analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confResp');
                        meanConf.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw = vertcat(meanConf.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).meanConf);
                        meanConf.(group).(subjects{sub}).correct.raw = vertcat(meanConf.(group).(subjects{sub}).correct.raw, nanmean(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confResp(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).acc == 1)));
                        meanConf.(group).(subjects{sub}).incorrect.raw = vertcat(meanConf.(group).(subjects{sub}).incorrect.raw,  nanmean(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confResp(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).acc == 0)));
                    end
                end
            else % Sessions 2-9
                confResp.(group).(sessions{sesh}).raw = vertcat(confResp.(group).(sessions{sesh}).raw, analysis.(subjects{sub}).(sessions{sesh}).perception.trained.confResp');
                meanConf.(group).(sessions{sesh}).perception.trained.raw = vertcat(meanConf.(group).(sessions{sesh}).perception.trained.raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.meanConf);
                meanConf.(group).(subjects{sub}).correct.raw = vertcat(meanConf.(group).(subjects{sub}).correct.raw, nanmean(analysis.(subjects{sub}).(sessions{sesh}).perception.trained.confResp(analysis.(subjects{sub}).(sessions{sesh}).perception.trained.acc == 1)));
                meanConf.(group).(subjects{sub}).incorrect.raw = vertcat(meanConf.(group).(subjects{sub}).incorrect.raw,  nanmean(analysis.(subjects{sub}).(sessions{sesh}).perception.trained.confResp(analysis.(subjects{sub}).(sessions{sesh}).perception.trained.acc == 0)));
            end
        end
    end    
end

% Take mean and standard error
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    [plots.meanConf.(groups{g}).learningCurve.mean, plots.meanConf.(groups{g}).learningCurve.sem] = deal([]);
    [plots.confResp.(groups{g}).allSessions.mean, plots.confResp.(groups{g}).allSessions.sem] = deal([]);
    [meanConf.(groups{g}).correct.raw, meanConf.(groups{g}).incorrect.raw] = deal([]);
    subjects = fieldnames(meanConf.(groups{g}));
    for sub = 1:numel(subjects)
        if strncmp(subjects{sub},'subject',7)
            meanConf.(groups{g}).(subjects{sub}).correct.mean = nanmean(meanConf.(groups{g}).(subjects{sub}).correct.raw);
            meanConf.(groups{g}).(subjects{sub}).incorrect.mean = nanmean(meanConf.(groups{g}).(subjects{sub}).incorrect.raw);
            meanConf.(groups{g}).correct.raw = vertcat(meanConf.(groups{g}).correct.raw, meanConf.(groups{g}).(subjects{sub}).correct.mean);
            meanConf.(groups{g}).incorrect.raw = vertcat(meanConf.(groups{g}).incorrect.raw, meanConf.(groups{g}).(subjects{sub}).incorrect.mean);
        end
    end
    meanConf.(groups{g}).correct.mean = nanmean(meanConf.(groups{g}).correct.raw);
    meanConf.(groups{g}).correct.sem = nanstd(meanConf.(groups{g}).correct.raw)/sqrt(length(meanConf.(groups{g}).correct.raw));
    meanConf.(groups{g}).incorrect.mean = nanmean(meanConf.(groups{g}).incorrect.raw);
    meanConf.(groups{g}).incorrect.sem = nanstd(meanConf.(groups{g}).incorrect.raw)/sqrt(length(meanConf.(groups{g}).incorrect.raw));
    [meanConf.(groups{g}).h, meanConf.(groups{g}).p, meanConf.(groups{g}).ci, meanConf.(groups{g}).stats] = ttest(meanConf.(groups{g}).correct.raw, meanConf.(groups{g}).incorrect.raw);
    for sesh = 1:10
        session = sprintf('session_%.2d',sesh);
        confResp.(groups{g}).(session).mean = nanmean(confResp.(groups{g}).(session).raw,1);
        confResp.(groups{g}).(session).sem = nanstd(confResp.(groups{g}).(session).raw,[],1)/sqrt(size(confResp.(groups{g}).(session).raw,1));
        plots.confResp.(groups{g}).allSessions.mean = horzcat(plots.confResp.(groups{g}).allSessions.mean, confResp.(groups{g}).(session).mean);
        plots.confResp.(groups{g}).allSessions.sem = horzcat(plots.confResp.(groups{g}).allSessions.sem, confResp.(groups{g}).(session).sem);
        if sesh == 1 || sesh == 10
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    meanConf.(groups{g}).(session).(dom{d}).(stim{s}).mean = nanmean(meanConf.(groups{g}).(session).(dom{d}).(stim{s}).raw);
                    meanConf.(groups{g}).(session).(dom{d}).(stim{s}).sem = nanstd(meanConf.(groups{g}).(session).(dom{d}).(stim{s}).raw)/sqrt(length(meanConf.(groups{g}).(session).(dom{d}).(stim{s}).raw));
                    meanConf.(groups{g}).(session).(dom{d}).(stim{s}).std = nanstd(meanConf.(groups{g}).(session).(dom{d}).(stim{s}).raw);
                end
            end
            plots.meanConf.(groups{g}).(session).mean = [meanConf.(groups{g}).(session).perception.trained.mean, meanConf.(groups{g}).(session).perception.untrained.mean;...
                meanConf.(groups{g}).(session).memory.trained.mean, meanConf.(groups{g}).(session).memory.untrained.mean];
            plots.meanConf.(groups{g}).(session).sem = [meanConf.(groups{g}).(session).perception.trained.sem, meanConf.(groups{g}).(session).perception.untrained.sem;...
                meanConf.(groups{g}).(session).memory.trained.sem, meanConf.(groups{g}).(session).memory.untrained.sem];
        else % sessions 2-9
            meanConf.(groups{g}).(session).perception.trained.mean = nanmean(meanConf.(groups{g}).(session).perception.trained.raw);
            meanConf.(groups{g}).(session).perception.trained.sem = nanstd(meanConf.(groups{g}).(session).perception.trained.raw)/sqrt(length(meanConf.(groups{g}).(session).perception.trained.raw));
            plots.meanConf.(groups{g}).learningCurve.mean = vertcat(plots.meanConf.(groups{g}).learningCurve.mean, meanConf.(groups{g}).(session).perception.trained.mean);
            plots.meanConf.(groups{g}).learningCurve.sem = vertcat(plots.meanConf.(groups{g}).learningCurve.sem, meanConf.(groups{g}).(session).perception.trained.sem);
        end
    end
end

if plotFigs
    confRespPlot_twoSessions = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    lineProps = struct('color', 'c', 'width', 2);
    mseb(1:numel([confResp.group_1.session_01.mean, confResp.group_1.session_02.mean]), [confResp.group_1.session_01.mean, confResp.group_1.session_02.mean], [confResp.group_1.session_01.sem, confResp.group_1.session_02.sem], lineProps); hold on;
    plot([108.5 108.5], [1 4], '-k', 'linewidth', 1);
    xlim([1 numel([confResp.group_1.session_01.mean, confResp.group_1.session_02.mean])]); ylim([1 4]);
    set(gca,'fontsize',14);
    t(1) = text(3,1.225, 'Session 1');
    t(2) = text(197,1.225, 'Session 2');
    set(t, 'fontsize', 14);
    title('Control Group');
    xlabel('Trials');
    ylabel('Confidence');
    subplot(1,2,2); % Experimental Group
    mseb(1:numel([confResp.group_2.session_01.mean, confResp.group_2.session_02.mean]), [confResp.group_2.session_01.mean, confResp.group_2.session_02.mean], [confResp.group_2.session_01.sem, confResp.group_2.session_02.sem], lineProps); hold on;
    plot([108.5 108.5], [1 4], '-k', 'linewidth', 1);
    xlim([1 numel([confResp.group_2.session_01.mean, confResp.group_2.session_02.mean])]); ylim([1 4]);
    set(gca,'fontsize',14);
    t(1) = text(3,1.225, 'Session 1');
    t(2) = text(197,1.225, 'Session 2');
    set(t, 'fontsize', 14);
    title('Experimental Group');
    xlabel('Trials');
    ylabel('Confidence');
    if exportFigs
        export_fig confRespPlot_twoSessions -png -transparent 'confResp_twoSessions.png';
    end
    
    confRespPlot_allSessions = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    lineProps = struct('color', 'c', 'width', 2);
    mseb(1:numel(plots.confResp.group_1.allSessions.mean), plots.confResp.group_1.allSessions.mean, plots.confResp.group_1.allSessions.sem, lineProps); hold on;
    for t = 108.5:270:2268.5
        plot([t t], [1 4], '-k', 'linewidth', 1);
    end
    xlim([1 numel(plots.confResp.group_1.allSessions.mean)]); ylim([1 4]);
    set(gca,'fontsize',14);
    title('Control Group');
    xlabel('Trials');
    ylabel('Confidence');
    subplot(1,2,2); % Experimental Group
    mseb(1:numel(plots.confResp.group_2.allSessions.mean), plots.confResp.group_2.allSessions.mean, plots.confResp.group_2.allSessions.sem, lineProps); hold on;
    for t = 108.5:270:2268.5
        plot([t t], [1 4], '-k', 'linewidth', 1);
    end
    xlim([1 numel(plots.confResp.group_2.allSessions.mean)]); ylim([1 4]);
    set(gca,'fontsize',14);
    title('Experimental Group');
    xlabel('Trials');
    ylabel('Confidence');
    if exportFigs
        export_fig confRespPlot_allSessions -png -transparent 'confResp_allSessions.png';
    end
    
    
    % Error Bar Comparison Plot
    meanConfComparisonPlot = figure;
    set(gcf,'position', [200 200 450 300]);
    subplot(1,2,1); % Control Group
    [hBar hErrorbar] = barwitherr([[plots.meanConf.group_1.session_01.sem(1,:), plots.meanConf.group_1.session_01.sem(2,:)]',[plots.meanConf.group_1.session_10.sem(1,:), plots.meanConf.group_1.session_10.sem(2,:)]']',...
        [[plots.meanConf.group_1.session_01.mean(1,:), plots.meanConf.group_1.session_01.mean(2,:)]',[plots.meanConf.group_1.session_10.mean(1,:), plots.meanConf.group_1.session_10.mean(2,:)]']');
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    ylim([0 4]);
    set(gca, 'fontsize', 14);
    ylabel('Mean Confidence', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 6); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    [hBar hErrorbar] = barwitherr([[plots.meanConf.group_2.session_01.sem(1,:), plots.meanConf.group_2.session_01.sem(2,:)]',[plots.meanConf.group_2.session_10.sem(1,:), plots.meanConf.group_2.session_10.sem(2,:)]']',...
        [[plots.meanConf.group_2.session_01.mean(1,:), plots.meanConf.group_2.session_01.mean(2,:)]',[plots.meanConf.group_2.session_10.mean(1,:), plots.meanConf.group_2.session_10.mean(2,:)]']');
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    ylim([0 4]);
    set(gca, 'fontsize', 14);
    ylabel('Mean Confidence', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Experimental Group', 'fontsize', 14);
    box off;
    if exportFigs
        export_fig meanConfComparisonPlot -png -transparent 'meanConfComparison.png';
    end
    
    % Scatter Plot
    meanConfScatterPlot = figure;
    set(gcf, 'position', [200 200 600 450]);
    subplot(2,2,1); % Perception, Control Group
    hMarker(1) = scatter(meanConf.group_1.session_01.perception.trained.raw, meanConf.group_1.session_10.perception.trained.raw, 100, 'o', 'markerfacecolor', [.5 0 0], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(meanConf.group_1.session_01.perception.untrained.raw, meanConf.group_1.session_10.perception.untrained.raw, 100, 'o', 'markerfacecolor', [1 0 0], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([0 4]); ylim([0 4]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('CG: Perception confidence');
    leg = legend('Trained', 'Untrained', 'location', 'nw');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    subplot(2,2,2); % Perception, Experimental Group
    hMarker(1) = scatter(meanConf.group_2.session_01.perception.trained.raw, meanConf.group_2.session_10.perception.trained.raw, 100, 'd', 'markerfacecolor', [.5 0 0], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(meanConf.group_2.session_01.perception.untrained.raw, meanConf.group_2.session_10.perception.untrained.raw, 100, 'd', 'markerfacecolor', [1 0 0], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([0 4]); ylim([0 4]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('EG: Perception confidence');
    leg = legend('Trained', 'Untrained', 'location', 'se');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    subplot(2,2,3); % Memory, Control Group
    hMarker(1) = scatter(meanConf.group_1.session_01.memory.trained.raw, meanConf.group_1.session_10.memory.trained.raw, 100, 'o', 'markerfacecolor', [0 0 0.5], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(meanConf.group_1.session_01.memory.untrained.raw, meanConf.group_1.session_10.memory.untrained.raw, 100, 'o', 'markerfacecolor', [0 0 1], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([0 4]); ylim([0 4]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('CG: Memory confidence');
    leg = legend('Trained', 'Untrained', 'location', 'nw');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    subplot(2,2,4); % Memory, Experimental Group
    hMarker(1) = scatter(meanConf.group_2.session_01.memory.trained.raw, meanConf.group_2.session_10.memory.trained.raw, 100, 'd', 'markerfacecolor', [0 0 .5], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(meanConf.group_2.session_01.memory.untrained.raw, meanConf.group_2.session_10.memory.untrained.raw, 100, 'd', 'markerfacecolor', [0 0 1], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([0 4]); ylim([0 4]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('EG: Memory confidence');
    leg = legend('Trained', 'Untrained', 'location', 'se');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    if exportFigs
        export_fig meanConfScatterPlot -png -transparent 'meanConfScatter.png';
    end
    
    meanConfCurvePlot = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    hBar(1) = errorbar(.6, meanConf.group_1.session_01.perception.trained.mean, meanConf.group_1.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, meanConf.group_1.session_01.perception.untrained.mean, meanConf.group_1.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, meanConf.group_1.session_01.memory.trained.mean, meanConf.group_1.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, meanConf.group_1.session_01.memory.untrained.mean, meanConf.group_1.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.meanConf.group_1.learningCurve.mean, plots.meanConf.group_1.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, meanConf.group_1.session_10.perception.trained.mean, meanConf.group_1.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, meanConf.group_1.session_10.perception.untrained.mean, meanConf.group_1.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, meanConf.group_1.session_10.memory.trained.mean, meanConf.group_1.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, meanConf.group_1.session_10.memory.untrained.mean, meanConf.group_1.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    ylim([2 4]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('Mean Confidence', 'fontsize', 14);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    hBar(1) = errorbar(.6, meanConf.group_2.session_01.perception.trained.mean, meanConf.group_2.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, meanConf.group_2.session_01.perception.untrained.mean, meanConf.group_2.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, meanConf.group_2.session_01.memory.trained.mean, meanConf.group_2.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, meanConf.group_2.session_01.memory.untrained.mean, meanConf.group_2.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.meanConf.group_2.learningCurve.mean, plots.meanConf.group_2.learningCurve.sem,'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, meanConf.group_2.session_10.perception.trained.mean, meanConf.group_2.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, meanConf.group_2.session_10.perception.untrained.mean, meanConf.group_2.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, meanConf.group_2.session_10.memory.trained.mean, meanConf.group_2.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, meanConf.group_2.session_10.memory.untrained.mean, meanConf.group_2.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    ylim([2 4]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('Mean Confidence', 'fontsize', 14);
    title('Experimental Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    if exportFigs
        export_fig meanConfCurvePlot -png -transparent 'meanConfCurve.png';
    end
    
end
end