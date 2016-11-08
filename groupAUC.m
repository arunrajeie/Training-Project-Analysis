function [ AUC ] = groupAUC( results, plotFigs, exportFigs )
%GROUPAUC Runs group AUC analysis and optionally plots and exports the figures

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
                    AUC.(groups{g}).(session).(dom{d}).(stim{s}).raw = [];
                end
            end
        else  % sessions 2-9
            AUC.(groups{g}).(session).perception.trained.raw = [];
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
                        AUC.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw = vertcat(AUC.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).AUC);
                    end
                end
            else % Sessions 2-9
                AUC.(group).(sessions{sesh}).perception.trained.raw = vertcat(AUC.(group).(sessions{sesh}).perception.trained.raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.AUC);
            end
        end
    end    
end

% Take mean and standard error
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    [plots.AUC.(groups{g}).learningCurve.mean, plots.AUC.(groups{g}).learningCurve.sem] = deal([]);
    for sesh = 1:10
        session = sprintf('session_%.2d',sesh);
        if sesh == 1 || sesh == 10
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    AUC.(groups{g}).(session).(dom{d}).(stim{s}).mean = nanmean(AUC.(groups{g}).(session).(dom{d}).(stim{s}).raw);
                    AUC.(groups{g}).(session).(dom{d}).(stim{s}).sem = nanstd(AUC.(groups{g}).(session).(dom{d}).(stim{s}).raw)/sqrt(length(AUC.(groups{g}).(session).(dom{d}).(stim{s}).raw));
                end
            end
            plots.AUC.(groups{g}).(session).mean = [AUC.(groups{g}).(session).perception.trained.mean, AUC.(groups{g}).(session).perception.untrained.mean;...
                AUC.(groups{g}).(session).memory.trained.mean, AUC.(groups{g}).(session).memory.untrained.mean];
            plots.AUC.(groups{g}).(session).sem = [AUC.(groups{g}).(session).perception.trained.sem, AUC.(groups{g}).(session).perception.untrained.sem;...
                AUC.(groups{g}).(session).memory.trained.sem, AUC.(groups{g}).(session).memory.untrained.sem];
        else % sessions 2-9
            AUC.(groups{g}).(session).perception.trained.mean = nanmean(AUC.(groups{g}).(session).perception.trained.raw);
            AUC.(groups{g}).(session).perception.trained.sem = nanstd(AUC.(groups{g}).(session).perception.trained.raw)/sqrt(length(AUC.(groups{g}).(session).perception.trained.raw));
            plots.AUC.(groups{g}).learningCurve.mean = vertcat(plots.AUC.(groups{g}).learningCurve.mean, AUC.(groups{g}).(session).perception.trained.mean);
            plots.AUC.(groups{g}).learningCurve.sem = vertcat(plots.AUC.(groups{g}).learningCurve.sem, AUC.(groups{g}).(session).perception.trained.sem);
        end
    end
    
end

if plotFigs
    % Error Bar Comparison Plot
    AUCComparisonPlot = figure;
    set(gcf,'position', [200 200 450 300]);
    subplot(1,2,1); % Control Group
    [hBar hErrorbar] = barwitherr([[plots.AUC.group_1.session_01.sem(1,:), plots.AUC.group_1.session_01.sem(2,:)]',[plots.AUC.group_1.session_10.sem(1,:), plots.AUC.group_1.session_10.sem(2,:)]']',...
        [[plots.AUC.group_1.session_01.mean(1,:), plots.AUC.group_1.session_01.mean(2,:)]',[plots.AUC.group_1.session_10.mean(1,:), plots.AUC.group_1.session_10.mean(2,:)]']');
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    ylim([.6 .825]);
    set(gca, 'fontsize', 14);
    ylabel('AUC', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 6); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    [hBar hErrorbar] = barwitherr([[plots.AUC.group_2.session_01.sem(1,:), plots.AUC.group_2.session_01.sem(2,:)]',[plots.AUC.group_2.session_10.sem(1,:), plots.AUC.group_2.session_10.sem(2,:)]']',...
        [[plots.AUC.group_2.session_01.mean(1,:), plots.AUC.group_2.session_01.mean(2,:)]',[plots.AUC.group_2.session_10.mean(1,:), plots.AUC.group_2.session_10.mean(2,:)]']');
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    ylim([.6 .825]);
    set(gca, 'fontsize', 14);
    ylabel('AUC', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Experimental Group', 'fontsize', 14);
    box off;
    if exportFigs
        export_fig AUCComparisonPlot -png -transparent 'AUCComparison.png';
    end
    
    % Scatter Plot
    AUCScatterPlot = figure;
    set(gcf, 'position', [200 200 600 450]);
    subplot(2,2,1); % Perception, Control Group
    hMarker(1) = scatter(AUC.group_1.session_01.perception.trained.raw, AUC.group_1.session_10.perception.trained.raw, 100, 'o', 'markerfacecolor', [.5 0 0], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(AUC.group_1.session_01.perception.untrained.raw, AUC.group_1.session_10.perception.untrained.raw, 100, 'o', 'markerfacecolor', [1 0 0], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([.5 1]); ylim([.5 1]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('CG: Perception AUC');
    leg = legend('Trained', 'Untrained', 'location', 'nw');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    subplot(2,2,2); % Perception, Experimental Group
    hMarker(1) = scatter(AUC.group_2.session_01.perception.trained.raw, AUC.group_2.session_10.perception.trained.raw, 100, 'd', 'markerfacecolor', [.5 0 0], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(AUC.group_2.session_01.perception.untrained.raw, AUC.group_2.session_10.perception.untrained.raw, 100, 'd', 'markerfacecolor', [1 0 0], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([.5 1]); ylim([.5 1]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('EG: Perception AUC');
    leg = legend('Trained', 'Untrained', 'location', 'se');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    subplot(2,2,3); % Memory, Control Group
    hMarker(1) = scatter(AUC.group_1.session_01.memory.trained.raw, AUC.group_1.session_10.memory.trained.raw, 100, 'o', 'markerfacecolor', [0 0 0.5], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(AUC.group_1.session_01.memory.untrained.raw, AUC.group_1.session_10.memory.untrained.raw, 100, 'o', 'markerfacecolor', [0 0 1], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([.5 1]); ylim([.5 1]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('CG: Memory AUC');
    leg = legend('Trained', 'Untrained', 'location', 'nw');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    subplot(2,2,4); % Memory, Experimental Group
    hMarker(1) = scatter(AUC.group_2.session_01.memory.trained.raw, AUC.group_2.session_10.memory.trained.raw, 100, 'd', 'markerfacecolor', [0 0 .5], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(AUC.group_2.session_01.memory.untrained.raw, AUC.group_2.session_10.memory.untrained.raw, 100, 'd', 'markerfacecolor', [0 0 1], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([.5 1]); ylim([.5 1]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('EG: Memory AUC');
    leg = legend('Trained', 'Untrained', 'location', 'se');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    if exportFigs
        export_fig AUCScatterPlot -png -transparent 'AUCScatter.png';
    end
    
    AUCCurvePlot = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    hBar(1) = errorbar(.6, AUC.group_1.session_01.perception.trained.mean, AUC.group_1.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, AUC.group_1.session_01.perception.untrained.mean, AUC.group_1.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, AUC.group_1.session_01.memory.trained.mean, AUC.group_1.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, AUC.group_1.session_01.memory.untrained.mean, AUC.group_1.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.AUC.group_1.learningCurve.mean, plots.AUC.group_1.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, AUC.group_1.session_10.perception.trained.mean, AUC.group_1.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, AUC.group_1.session_10.perception.untrained.mean, AUC.group_1.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, AUC.group_1.session_10.memory.trained.mean, AUC.group_1.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, AUC.group_1.session_10.memory.untrained.mean, AUC.group_1.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    ylim([.6, .85]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('AUC', 'fontsize', 14);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    hBar(1) = errorbar(.6, AUC.group_2.session_01.perception.trained.mean, AUC.group_2.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, AUC.group_2.session_01.perception.untrained.mean, AUC.group_2.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, AUC.group_2.session_01.memory.trained.mean, AUC.group_2.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, AUC.group_2.session_01.memory.untrained.mean, AUC.group_2.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.AUC.group_2.learningCurve.mean, plots.AUC.group_2.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, AUC.group_2.session_10.perception.trained.mean, AUC.group_2.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, AUC.group_2.session_10.perception.untrained.mean, AUC.group_2.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, AUC.group_2.session_10.memory.trained.mean, AUC.group_2.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, AUC.group_2.session_10.memory.untrained.mean, AUC.group_2.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    ylim([.6, .85]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('AUC', 'fontsize', 14);
    title('Experimental Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    if exportFigs
        export_fig AUCCurvePlot -png -transparent 'AUCCurve.png';
    end
    
end
end