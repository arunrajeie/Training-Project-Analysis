function [ t2HR, t2FAR ] = groupT2HRFAR( results, plotFigs, exportFigs )
%GROUPT2HRFAR Runs group type 2 hit rate/false alarm rate analysis and optionally plots and exports the figures

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
                    t2HR.(groups{g}).(session).(dom{d}).(stim{s}).raw = [];
                    t2FAR.(groups{g}).(session).(dom{d}).(stim{s}).raw = [];
                end
            end
        else  % sessions 2-9
            t2HR.(groups{g}).(session).perception.trained.raw = [];
            t2FAR.(groups{g}).(session).perception.trained.raw = [];
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
                        t2HR.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw = vertcat(t2HR.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).t2HR);
                        t2FAR.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw = vertcat(t2FAR.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).t2FAR);
                    end
                end
            else % Sessions 2-9
                t2HR.(group).(sessions{sesh}).perception.trained.raw = vertcat(t2HR.(group).(sessions{sesh}).perception.trained.raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.t2HR);
                t2FAR.(group).(sessions{sesh}).perception.trained.raw = vertcat(t2FAR.(group).(sessions{sesh}).perception.trained.raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.t2FAR);
            end
        end
    end    
end

% Take mean and standard error
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    [plots.t2HR.(groups{g}).learningCurve.mean, plots.t2HR.(groups{g}).learningCurve.sem,...
        plots.t2FAR.(groups{g}).learningCurve.mean, plots.t2FAR.(groups{g}).learningCurve.sem] = deal([]);
    for sesh = 1:10
        session = sprintf('session_%.2d',sesh);
        if sesh == 1 || sesh == 10
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    t2HR.(groups{g}).(session).(dom{d}).(stim{s}).mean = nanmean(t2HR.(groups{g}).(session).(dom{d}).(stim{s}).raw);
                    t2HR.(groups{g}).(session).(dom{d}).(stim{s}).sem = nanstd(t2HR.(groups{g}).(session).(dom{d}).(stim{s}).raw)/sqrt(length(t2HR.(groups{g}).(session).(dom{d}).(stim{s}).raw));
                    t2FAR.(groups{g}).(session).(dom{d}).(stim{s}).mean = nanmean(t2FAR.(groups{g}).(session).(dom{d}).(stim{s}).raw);
                    t2FAR.(groups{g}).(session).(dom{d}).(stim{s}).sem = nanstd(t2FAR.(groups{g}).(session).(dom{d}).(stim{s}).raw)/sqrt(length(t2FAR.(groups{g}).(session).(dom{d}).(stim{s}).raw));
                end
            end
            plots.t2HR.(groups{g}).(session).mean = [t2HR.(groups{g}).(session).perception.trained.mean, t2HR.(groups{g}).(session).perception.untrained.mean;...
                t2HR.(groups{g}).(session).memory.trained.mean, t2HR.(groups{g}).(session).memory.untrained.mean];
            plots.t2HR.(groups{g}).(session).sem = [t2HR.(groups{g}).(session).perception.trained.sem, t2HR.(groups{g}).(session).perception.untrained.sem;...
                t2HR.(groups{g}).(session).memory.trained.sem, t2HR.(groups{g}).(session).memory.untrained.sem];
            plots.t2FAR.(groups{g}).(session).mean = [t2FAR.(groups{g}).(session).perception.trained.mean, t2FAR.(groups{g}).(session).perception.untrained.mean;...
                t2FAR.(groups{g}).(session).memory.trained.mean, t2FAR.(groups{g}).(session).memory.untrained.mean];
            plots.t2FAR.(groups{g}).(session).sem = [t2FAR.(groups{g}).(session).perception.trained.sem, t2FAR.(groups{g}).(session).perception.untrained.sem;...
                t2FAR.(groups{g}).(session).memory.trained.sem, t2FAR.(groups{g}).(session).memory.untrained.sem];
        else % sessions 2-9
            t2HR.(groups{g}).(session).perception.trained.mean = nanmean(t2HR.(groups{g}).(session).perception.trained.raw);
            t2HR.(groups{g}).(session).perception.trained.sem = nanstd(t2HR.(groups{g}).(session).perception.trained.raw)/sqrt(length(t2HR.(groups{g}).(session).perception.trained.raw));
            t2FAR.(groups{g}).(session).perception.trained.mean = nanmean(t2FAR.(groups{g}).(session).perception.trained.raw);
            t2FAR.(groups{g}).(session).perception.trained.sem = nanstd(t2HR.(groups{g}).(session).perception.trained.raw)/sqrt(length(t2FAR.(groups{g}).(session).perception.trained.raw));
            plots.t2HR.(groups{g}).learningCurve.mean = vertcat(plots.t2HR.(groups{g}).learningCurve.mean, t2HR.(groups{g}).(session).perception.trained.mean);
            plots.t2HR.(groups{g}).learningCurve.sem = vertcat(plots.t2HR.(groups{g}).learningCurve.sem, t2HR.(groups{g}).(session).perception.trained.sem);
            plots.t2FAR.(groups{g}).learningCurve.mean = vertcat(plots.t2FAR.(groups{g}).learningCurve.mean, t2FAR.(groups{g}).(session).perception.trained.mean);
            plots.t2FAR.(groups{g}).learningCurve.sem = vertcat(plots.t2FAR.(groups{g}).learningCurve.sem, t2FAR.(groups{g}).(session).perception.trained.sem);
        end
    end 
end

% Plot outcomes
if plotFigs
    % Error Bar Comparison Plot
    t2HRFARComparisonPlot = figure;
    set(gcf,'position', [200 200 450 300]);
    subplot(1,2,1); % Control Group
    hErrorBar(1) = errorbar([.7,6.7], [t2HR.group_1.session_01.perception.trained.mean, t2HR.group_1.session_10.perception.trained.mean], [t2HR.group_1.session_01.perception.trained.sem, t2HR.group_1.session_10.perception.trained.sem], '-o', 'linewidth',2); hold on;
    hErrorBar(2) = errorbar([.9,6.9], [t2HR.group_1.session_01.perception.untrained.mean, t2HR.group_1.session_10.perception.untrained.mean], [t2HR.group_1.session_01.perception.untrained.sem, t2HR.group_1.session_10.perception.untrained.sem], '-o', 'linewidth', 2);
    hErrorBar(3) = errorbar([1.1,7.1], [t2HR.group_1.session_01.memory.trained.mean, t2HR.group_1.session_10.memory.trained.mean], [t2HR.group_1.session_01.memory.trained.sem, t2HR.group_1.session_10.memory.trained.sem], '-o', 'linewidth', 2);
    hErrorBar(4) = errorbar([1.3,7.3], [t2HR.group_1.session_01.memory.untrained.mean, t2HR.group_1.session_10.memory.untrained.mean], [t2HR.group_1.session_01.memory.untrained.sem, t2HR.group_1.session_10.memory.untrained.sem], '-o', 'linewidth', 2);
    set(hErrorBar(1), 'color', [.5 0 0]);
    set(hErrorBar(2), 'color', [1 0 0]);
    set(hErrorBar(3), 'color', [0 0 .5]);
    set(hErrorBar(4), 'color', [0 0 1]);
    hErrorBar(1) = errorbar([2.7,8.7], [t2FAR.group_1.session_01.perception.trained.mean, t2FAR.group_1.session_10.perception.trained.mean], [t2FAR.group_1.session_01.perception.trained.sem, t2FAR.group_1.session_10.perception.trained.sem], '-o', 'linewidth',2); hold on;
    hErrorBar(2) = errorbar([2.9,8.9], [t2FAR.group_1.session_01.perception.untrained.mean, t2FAR.group_1.session_10.perception.untrained.mean], [t2FAR.group_1.session_01.perception.untrained.sem, t2FAR.group_1.session_10.perception.untrained.sem], '-o', 'linewidth', 2);
    hErrorBar(3) = errorbar([3.1,9.1], [t2FAR.group_1.session_01.memory.trained.mean, t2FAR.group_1.session_10.memory.trained.mean], [t2FAR.group_1.session_01.memory.trained.sem, t2FAR.group_1.session_10.memory.trained.sem], '-o', 'linewidth', 2);
    hErrorBar(4) = errorbar([3.3,9.3], [t2FAR.group_1.session_01.memory.untrained.mean, t2FAR.group_1.session_10.memory.untrained.mean], [t2FAR.group_1.session_01.memory.untrained.sem, t2FAR.group_1.session_10.memory.untrained.sem], '-o', 'linewidth', 2);
    set(hErrorBar(1), 'color', [.5 0 0]);
    set(hErrorBar(2), 'color', [1 0 0]);
    set(hErrorBar(3), 'color', [0 0 .5]);
    set(hErrorBar(4), 'color', [0 0 1]);
    ylim([0 1]);
    set(gca, 'fontsize', 14);
    ylabel('T2 HR/FAR', 'fontsize', 14);
    set(gca,'xtick', [1,3,7,9], 'xticklabel', {'HR', 'FAR', 'HR', 'FAR'}, 'fontsize', 11);
    title('Control Group', 'fontsize', 14);
    t(1) = text(.9,.05,'Pre');
    t(2) = text(6.5,.05,'Post');
    set(t, 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 6); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    hErrorBar(1) = errorbar([.7,6.7], [t2HR.group_2.session_01.perception.trained.mean, t2HR.group_2.session_10.perception.trained.mean], [t2HR.group_2.session_01.perception.trained.sem, t2HR.group_2.session_10.perception.trained.sem], '-o', 'linewidth',2); hold on;
    hErrorBar(2) = errorbar([.9,6.9], [t2HR.group_2.session_01.perception.untrained.mean, t2HR.group_2.session_10.perception.untrained.mean], [t2HR.group_2.session_01.perception.untrained.sem, t2HR.group_2.session_10.perception.untrained.sem], '-o', 'linewidth', 2);
    hErrorBar(3) = errorbar([1.1,7.1], [t2HR.group_2.session_01.memory.trained.mean, t2HR.group_2.session_10.memory.trained.mean], [t2HR.group_2.session_01.memory.trained.sem, t2HR.group_2.session_10.memory.trained.sem], '-o', 'linewidth', 2);
    hErrorBar(4) = errorbar([1.3,7.3], [t2HR.group_2.session_01.memory.untrained.mean, t2HR.group_2.session_10.memory.untrained.mean], [t2HR.group_2.session_01.memory.untrained.sem, t2HR.group_2.session_10.memory.untrained.sem], '-o', 'linewidth', 2);
    set(hErrorBar(1), 'color', [.5 0 0]);
    set(hErrorBar(2), 'color', [1 0 0]);
    set(hErrorBar(3), 'color', [0 0 .5]);
    set(hErrorBar(4), 'color', [0 0 1]);
    hErrorBar(1) = errorbar([2.7,8.7], [t2FAR.group_2.session_01.perception.trained.mean, t2FAR.group_2.session_10.perception.trained.mean], [t2FAR.group_2.session_01.perception.trained.sem, t2FAR.group_2.session_10.perception.trained.sem], '-o', 'linewidth',2); hold on;
    hErrorBar(2) = errorbar([2.9,8.9], [t2FAR.group_2.session_01.perception.untrained.mean, t2FAR.group_2.session_10.perception.untrained.mean], [t2FAR.group_2.session_01.perception.untrained.sem, t2FAR.group_2.session_10.perception.untrained.sem], '-o', 'linewidth', 2);
    hErrorBar(3) = errorbar([3.1,9.1], [t2FAR.group_2.session_01.memory.trained.mean, t2FAR.group_2.session_10.memory.trained.mean], [t2FAR.group_2.session_01.memory.trained.sem, t2FAR.group_2.session_10.memory.trained.sem], '-o', 'linewidth', 2);
    hErrorBar(4) = errorbar([3.3,9.3], [t2FAR.group_2.session_01.memory.untrained.mean, t2FAR.group_2.session_10.memory.untrained.mean], [t2FAR.group_2.session_01.memory.untrained.sem, t2FAR.group_2.session_10.memory.untrained.sem], '-o', 'linewidth', 2);
    set(hErrorBar(1), 'color', [.5 0 0]);
    set(hErrorBar(2), 'color', [1 0 0]);
    set(hErrorBar(3), 'color', [0 0 .5]);
    set(hErrorBar(4), 'color', [0 0 1]);
    ylim([0 1]);
    set(gca, 'fontsize', 14);
    ylabel('T2 HR/FAR', 'fontsize', 14);
    set(gca,'xtick', [1,3,7,9], 'xticklabel', {'HR', 'FAR', 'HR', 'FAR'}, 'fontsize', 11);
    title('Experimental Group', 'fontsize', 14);
    t(1) = text(.9,.05,'Pre');
    t(2) = text(6.5,.05,'Post');
    set(t, 'fontsize', 14);
    box off;
    if exportFigs
        export_fig t2HRFARComparisonPlot -png -transparent 't2HRFARComparison.png';
    end
    
    t2HRFARCurvePlot = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    % T2HR
    hBar(1) = errorbar(.6, t2HR.group_1.session_01.perception.trained.mean, t2HR.group_1.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, t2HR.group_1.session_01.perception.untrained.mean, t2HR.group_1.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, t2HR.group_1.session_01.memory.trained.mean, t2HR.group_1.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, t2HR.group_1.session_01.memory.untrained.mean, t2HR.group_1.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.t2HR.group_1.learningCurve.mean, plots.t2HR.group_1.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, t2HR.group_1.session_10.perception.trained.mean, t2HR.group_1.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, t2HR.group_1.session_10.perception.untrained.mean, t2HR.group_1.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, t2HR.group_1.session_10.memory.trained.mean, t2HR.group_1.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, t2HR.group_1.session_10.memory.untrained.mean, t2HR.group_1.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    % T2FAR
    hBar(1) = errorbar(.6, t2FAR.group_1.session_01.perception.trained.mean, t2FAR.group_1.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, t2FAR.group_1.session_01.perception.untrained.mean, t2FAR.group_1.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, t2FAR.group_1.session_01.memory.trained.mean, t2FAR.group_1.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, t2FAR.group_1.session_01.memory.untrained.mean, t2FAR.group_1.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.t2FAR.group_1.learningCurve.mean, plots.t2FAR.group_1.learningCurve.sem,'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, t2FAR.group_1.session_10.perception.trained.mean, t2FAR.group_1.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, t2FAR.group_1.session_10.perception.untrained.mean, t2FAR.group_1.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, t2FAR.group_1.session_10.memory.trained.mean, t2FAR.group_1.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, t2FAR.group_1.session_10.memory.untrained.mean, t2FAR.group_1.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    ylim([0 1]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('T2 HR/FAR', 'fontsize', 14);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 8); 
    t(1) = text(4.33,.7,'T2 HR');
    t(2) = text(4.33,.4,'T2 FAR');
    set(t, 'fontsize', 14);
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    % T2HR
    hBar(1) = errorbar(.6, t2HR.group_2.session_01.perception.trained.mean, t2HR.group_2.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, t2HR.group_2.session_01.perception.untrained.mean, t2HR.group_2.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, t2HR.group_2.session_01.memory.trained.mean, t2HR.group_2.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, t2HR.group_2.session_01.memory.untrained.mean, t2HR.group_2.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.t2HR.group_2.learningCurve.mean, plots.t2HR.group_2.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, t2HR.group_2.session_10.perception.trained.mean, t2HR.group_2.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, t2HR.group_2.session_10.perception.untrained.mean, t2HR.group_2.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, t2HR.group_2.session_10.memory.trained.mean, t2HR.group_2.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, t2HR.group_2.session_10.memory.untrained.mean, t2HR.group_2.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    % T2FAR
    hBar(1) = errorbar(.6, t2FAR.group_2.session_01.perception.trained.mean, t2FAR.group_2.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, t2FAR.group_2.session_01.perception.untrained.mean, t2FAR.group_2.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, t2FAR.group_2.session_01.memory.trained.mean, t2FAR.group_2.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, t2FAR.group_2.session_01.memory.untrained.mean, t2FAR.group_2.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.t2FAR.group_2.learningCurve.mean, plots.t2FAR.group_2.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, t2FAR.group_2.session_10.perception.trained.mean, t2FAR.group_2.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, t2FAR.group_2.session_10.perception.untrained.mean, t2FAR.group_2.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, t2FAR.group_2.session_10.memory.trained.mean, t2FAR.group_2.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, t2FAR.group_2.session_10.memory.untrained.mean, t2FAR.group_2.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    ylim([0 1]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('T2 HR/FAR', 'fontsize', 14);
    title('Experimental Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'se');
    set(leg, 'FontSize', 8); 
    t(1) = text(4.33,.85,'T2 HR');
    t(2) = text(4.33,.55,'T2 FAR');
    set(t, 'fontsize', 14);
    legend boxoff; box off;
    if exportFigs
        export_fig t2HRFARCurvePlot -png -transparent 't2HRFARCurve.png';
    end

end
end

