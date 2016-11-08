function [ logM_ratio ] = groupLogM_ratio( results, plotFigs, exportFigs )
%GROUPLOGM_RATIO Runs group log(M_ratio) analysis and optionally plots and exports
%the figures

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
                    logM_ratio.(groups{g}).(session).(dom{d}).(stim{s}).raw = [];
                end
            end
        else  % sessions 2-9
            logM_ratio.(groups{g}).(session).perception.trained.raw = [];
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
                        logM_ratio.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw = vertcat(logM_ratio.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw, real(log(results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).fit.M_ratio)));
                    end
                end
            else % Sessions 2-9
                logM_ratio.(group).(sessions{sesh}).perception.trained.raw = vertcat(logM_ratio.(group).(sessions{sesh}).perception.trained.raw, real(log(results.(subjects{sub}).(sessions{sesh}).perception.trained.fit.M_ratio)));
            end
        end
    end    
end

% Take mean and standard error
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    [plots.logM_ratio.(groups{g}).learningCurve.mean, plots.logM_ratio.(groups{g}).learningCurve.sem] = deal([]);
    for sesh = 1:10
        session = sprintf('session_%.2d',sesh);
        if sesh == 1 || sesh == 10
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    logM_ratio.(groups{g}).(session).(dom{d}).(stim{s}).mean = nanmean(logM_ratio.(groups{g}).(session).(dom{d}).(stim{s}).raw);
                    logM_ratio.(groups{g}).(session).(dom{d}).(stim{s}).sem = nanstd(logM_ratio.(groups{g}).(session).(dom{d}).(stim{s}).raw)/sqrt(length(logM_ratio.(groups{g}).(session).(dom{d}).(stim{s}).raw));
                end
            end
            plots.logM_ratio.(groups{g}).(session).mean = [logM_ratio.(groups{g}).(session).perception.trained.mean, logM_ratio.(groups{g}).(session).perception.untrained.mean;...
                logM_ratio.(groups{g}).(session).memory.trained.mean, logM_ratio.(groups{g}).(session).memory.untrained.mean];
            plots.logM_ratio.(groups{g}).(session).sem = [logM_ratio.(groups{g}).(session).perception.trained.sem, logM_ratio.(groups{g}).(session).perception.untrained.sem;...
                logM_ratio.(groups{g}).(session).memory.trained.sem, logM_ratio.(groups{g}).(session).memory.untrained.sem];
        else % sessions 2-9
            logM_ratio.(groups{g}).(session).perception.trained.mean = nanmean(logM_ratio.(groups{g}).(session).perception.trained.raw);
            logM_ratio.(groups{g}).(session).perception.trained.sem = nanstd(logM_ratio.(groups{g}).(session).perception.trained.raw)/sqrt(length(logM_ratio.(groups{g}).(session).perception.trained.raw));
            plots.logM_ratio.(groups{g}).learningCurve.mean = vertcat(plots.logM_ratio.(groups{g}).learningCurve.mean, logM_ratio.(groups{g}).(session).perception.trained.mean);
            plots.logM_ratio.(groups{g}).learningCurve.sem = vertcat(plots.logM_ratio.(groups{g}).learningCurve.sem, logM_ratio.(groups{g}).(session).perception.trained.sem);
        end
    end
    
end

% Plot outcomes
if plotFigs
    % Error Bar Comparison Plot
    LogMRatioComparisonPlot = figure;
    set(gcf,'position', [200 200 450 300]);
    subplot(1,2,1); % Control Group
    [hBar hErrorbar] = barwitherr([[plots.logM_ratio.group_1.session_01.sem(1,:), plots.logM_ratio.group_1.session_01.sem(2,:)]',[plots.logM_ratio.group_1.session_10.sem(1,:), plots.logM_ratio.group_1.session_10.sem(2,:)]']',...
        [[plots.logM_ratio.group_1.session_01.mean(1,:), plots.logM_ratio.group_1.session_01.mean(2,:)]',[plots.logM_ratio.group_1.session_10.mean(1,:), plots.logM_ratio.group_1.session_10.mean(2,:)]']');
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    ylim([-.5 .8]);
    set(gca, 'fontsize', 14);
    ylabel('log(meta-d''/d'')', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 6); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    [hBar hErrorbar] = barwitherr([[plots.logM_ratio.group_2.session_01.sem(1,:), plots.logM_ratio.group_2.session_01.sem(2,:)]',[plots.logM_ratio.group_2.session_10.sem(1,:), plots.logM_ratio.group_2.session_10.sem(2,:)]']',...
        [[plots.logM_ratio.group_2.session_01.mean(1,:), plots.logM_ratio.group_2.session_01.mean(2,:)]',[plots.logM_ratio.group_2.session_10.mean(1,:), plots.logM_ratio.group_2.session_10.mean(2,:)]']');
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    ylim([-.5 .8]);
    set(gca, 'fontsize', 14);
    ylabel('log(meta-d''/d'')', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Experimental Group', 'fontsize', 14);
    box off;
    if exportFigs
        export_fig LogMRatioComparisonPlot -png -transparent 'LogMRatioComparison.png';
    end
    
    % Scatter Plot
    LogMRatioScatterPlot = figure;
    set(gcf, 'position', [200 200 600 450]);
    subplot(2,2,1); % Perception, Control Group
    hMarker(1) = scatter(logM_ratio.group_1.session_01.perception.trained.raw, logM_ratio.group_1.session_10.perception.trained.raw, 100, 'o', 'markerfacecolor', [.5 0 0], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(logM_ratio.group_1.session_01.perception.untrained.raw, logM_ratio.group_1.session_10.perception.untrained.raw, 100, 'o', 'markerfacecolor', [1 0 0], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([-3 3]); ylim([-3 3]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('CG: Perception log(meta-d''/d'')');
    leg = legend('Trained', 'Untrained', 'location', 'nw');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    subplot(2,2,2); % Perception, Experimental Group
    hMarker(1) = scatter(logM_ratio.group_2.session_01.perception.trained.raw, logM_ratio.group_2.session_10.perception.trained.raw, 100, 'd', 'markerfacecolor', [.5 0 0], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(logM_ratio.group_2.session_01.perception.untrained.raw, logM_ratio.group_2.session_10.perception.untrained.raw, 100, 'd', 'markerfacecolor', [1 0 0], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([-3 3]); ylim([-3 3]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('EG: Perception log(meta-d''/d'')');
    leg = legend('Trained', 'Untrained', 'location', 'se');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    subplot(2,2,3); % Memory, Control Group
    hMarker(1) = scatter(logM_ratio.group_1.session_01.memory.trained.raw, logM_ratio.group_1.session_10.memory.trained.raw, 100, 'o', 'markerfacecolor', [0 0 0.5], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(logM_ratio.group_1.session_01.memory.untrained.raw, logM_ratio.group_1.session_10.memory.untrained.raw, 100, 'o', 'markerfacecolor', [0 0 1], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([-3 3]); ylim([-3 3]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('CG: Memory log(meta-d''/d'')');
    leg = legend('Trained', 'Untrained', 'location', 'nw');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    subplot(2,2,4); % Memory, Experimental Group
    hMarker(1) = scatter(logM_ratio.group_2.session_01.memory.trained.raw, logM_ratio.group_2.session_10.memory.trained.raw, 100, 'd', 'markerfacecolor', [0 0 .5], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(logM_ratio.group_2.session_01.memory.untrained.raw, logM_ratio.group_2.session_10.memory.untrained.raw, 100, 'd', 'markerfacecolor', [0 0 1], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([-3 3]); ylim([-3 3]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('EG: Memory log(meta-d''/d'')');
    leg = legend('Trained', 'Untrained', 'location', 'se');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    if exportFigs
        export_fig LogMRatioScatterPlot -png -transparent 'LogMRatioScatter.png';
    end
    
    LogMCurvePlot = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    hBar(1) = errorbar(.6, logM_ratio.group_1.session_01.perception.trained.mean, logM_ratio.group_1.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, logM_ratio.group_1.session_01.perception.untrained.mean, logM_ratio.group_1.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, logM_ratio.group_1.session_01.memory.trained.mean, logM_ratio.group_1.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, logM_ratio.group_1.session_01.memory.untrained.mean, logM_ratio.group_1.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.logM_ratio.group_1.learningCurve.mean, plots.logM_ratio.group_1.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, logM_ratio.group_1.session_10.perception.trained.mean, logM_ratio.group_1.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, logM_ratio.group_1.session_10.perception.untrained.mean, logM_ratio.group_1.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, logM_ratio.group_1.session_10.memory.trained.mean, logM_ratio.group_1.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, logM_ratio.group_1.session_10.memory.untrained.mean, logM_ratio.group_1.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    ylim([-.5, .8]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('log(meta-d''/d'')', 'fontsize', 14);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    hBar(1) = errorbar(.6, logM_ratio.group_2.session_01.perception.trained.mean, logM_ratio.group_2.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, logM_ratio.group_2.session_01.perception.untrained.mean, logM_ratio.group_2.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, logM_ratio.group_2.session_01.memory.trained.mean, logM_ratio.group_2.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, logM_ratio.group_2.session_01.memory.untrained.mean, logM_ratio.group_2.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.logM_ratio.group_2.learningCurve.mean, plots.logM_ratio.group_2.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, logM_ratio.group_2.session_10.perception.trained.mean, logM_ratio.group_2.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, logM_ratio.group_2.session_10.perception.untrained.mean, logM_ratio.group_2.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, logM_ratio.group_2.session_10.memory.trained.mean, logM_ratio.group_2.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, logM_ratio.group_2.session_10.memory.untrained.mean, logM_ratio.group_2.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    ylim([-.5, .8]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('log(meta-d''/d'')', 'fontsize', 14);
    title('Experimental Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    if exportFigs
        export_fig LogMCurvePlot -png -transparent 'LogMCurve.png';
    end
end
       
end

