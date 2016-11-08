function [ M_ratio ] = groupM_ratio( results, plotFigs, exportFigs )
%GROUPM_RATIO Runs group M_ratio analysis and optionally plots and exports
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
                    M_ratio.(groups{g}).(session).(dom{d}).(stim{s}).raw = [];
                end
            end
        else  % sessions 2-9
            M_ratio.(groups{g}).(session).perception.trained.raw = [];
        end
    end
end

% Concatenate raw data
for sub = 1:numel(subjects)
    group = sprintf('group_%d', results.(subjects{sub}).group);
    sessions = fieldnames(results.(subjects{sub}));
    M_ratio.(group).(subjects{sub}).raw = [];
    for sesh = 1:numel(sessions)
        if strncmp(sessions{sesh},'session',7)
            session = str2double(sessions{sesh}(end-1:end));
            if session == 1 || session == 10
                for d = 1:numel(dom)
                    for s = 1:numel(stim)
                        M_ratio.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw = vertcat(M_ratio.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).fit.M_ratio);
                        M_ratio.(group).(subjects{sub}).raw = vertcat(M_ratio.(group).(subjects{sub}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).fit.M_ratio);
                    end
                end
            else % Sessions 2-9
                M_ratio.(group).(sessions{sesh}).perception.trained.raw = vertcat(M_ratio.(group).(sessions{sesh}).perception.trained.raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.fit.M_ratio);
                M_ratio.(group).(subjects{sub}).raw = vertcat(M_ratio.(group).(subjects{sub}).raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.fit.M_ratio);
            end
        end
    end    
end

% Take mean and standard error
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    [plots.M_ratio.(groups{g}).learningCurve.mean, plots.M_ratio.(groups{g}).learningCurve.sem] = deal([]);
    M_ratio.(groups{g}).raw = [];
    subjects = fieldnames(M_ratio.(groups{g}));
    for sub = 1:numel(subjects)
        if strncmp(subjects{sub},'subject',7)
            M_ratio.(groups{g}).(subjects{sub}).mean = nanmean(M_ratio.(groups{g}).(subjects{sub}).raw);
            M_ratio.(groups{g}).raw = vertcat(M_ratio.(groups{g}).raw, M_ratio.(groups{g}).(subjects{sub}).mean);
        end
    end
    M_ratio.(groups{g}).mean = nanmean(M_ratio.(groups{g}).raw);
    M_ratio.(groups{g}).sem = nanstd(M_ratio.(groups{g}).raw)/sqrt(length(M_ratio.(groups{g}).raw));
    for sesh = 1:10
        session = sprintf('session_%.2d',sesh);
        if sesh == 1 || sesh == 10
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    M_ratio.(groups{g}).(session).(dom{d}).(stim{s}).mean = nanmean(M_ratio.(groups{g}).(session).(dom{d}).(stim{s}).raw);
                    M_ratio.(groups{g}).(session).(dom{d}).(stim{s}).sem = nanstd(M_ratio.(groups{g}).(session).(dom{d}).(stim{s}).raw)/sqrt(length(M_ratio.(groups{g}).(session).(dom{d}).(stim{s}).raw));
                    M_ratio.(groups{g}).(session).(dom{d}).(stim{s}).std = nanstd(M_ratio.(groups{g}).(session).(dom{d}).(stim{s}).raw);
                end
            end
            plots.M_ratio.(groups{g}).(session).mean = [M_ratio.(groups{g}).(session).perception.trained.mean, M_ratio.(groups{g}).(session).perception.untrained.mean;...
                M_ratio.(groups{g}).(session).memory.trained.mean, M_ratio.(groups{g}).(session).memory.untrained.mean];
            plots.M_ratio.(groups{g}).(session).sem = [M_ratio.(groups{g}).(session).perception.trained.sem, M_ratio.(groups{g}).(session).perception.untrained.sem;...
                M_ratio.(groups{g}).(session).memory.trained.sem, M_ratio.(groups{g}).(session).memory.untrained.sem];
        else % sessions 2-9
            M_ratio.(groups{g}).(session).perception.trained.mean = nanmean(M_ratio.(groups{g}).(session).perception.trained.raw);
            M_ratio.(groups{g}).(session).perception.trained.sem = nanstd(M_ratio.(groups{g}).(session).perception.trained.raw)/sqrt(length(M_ratio.(groups{g}).(session).perception.trained.raw));
            plots.M_ratio.(groups{g}).learningCurve.mean = vertcat(plots.M_ratio.(groups{g}).learningCurve.mean, M_ratio.(groups{g}).(session).perception.trained.mean);
            plots.M_ratio.(groups{g}).learningCurve.sem = vertcat(plots.M_ratio.(groups{g}).learningCurve.sem, M_ratio.(groups{g}).(session).perception.trained.sem);
        end
    end
    
end

% Plot outcomes
if plotFigs
    % Error Bar Comparison Plot
    MRatioComparisonPlot = figure;
    set(gcf,'position', [200 200 450 300]);
    subplot(1,2,1); % Control Group
    [hBar hErrorbar] = barwitherr([[plots.M_ratio.group_1.session_01.sem(1,:), plots.M_ratio.group_1.session_01.sem(2,:)]',[plots.M_ratio.group_1.session_10.sem(1,:), plots.M_ratio.group_1.session_10.sem(2,:)]']',...
        [[plots.M_ratio.group_1.session_01.mean(1,:), plots.M_ratio.group_1.session_01.mean(2,:)]',[plots.M_ratio.group_1.session_10.mean(1,:), plots.M_ratio.group_1.session_10.mean(2,:)]']');
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    ylim([0 2.3]);
    set(gca, 'fontsize', 14);
    ylabel('meta-d''/d''', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 6); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    [hBar hErrorbar] = barwitherr([[plots.M_ratio.group_2.session_01.sem(1,:), plots.M_ratio.group_2.session_01.sem(2,:)]',[plots.M_ratio.group_2.session_10.sem(1,:), plots.M_ratio.group_2.session_10.sem(2,:)]']',...
        [[plots.M_ratio.group_2.session_01.mean(1,:), plots.M_ratio.group_2.session_01.mean(2,:)]',[plots.M_ratio.group_2.session_10.mean(1,:), plots.M_ratio.group_2.session_10.mean(2,:)]']');
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    ylim([0 2.3]);
    set(gca, 'fontsize', 14);
    ylabel('meta-d''/d''', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Experimental Group', 'fontsize', 14);
    box off;
    if exportFigs
        export_fig MRatioComparisonPlot -png -transparent 'MRatioComparison.png';
    end
    
    % Scatter Plot
    MRatioScatterPlot = figure;
    set(gcf, 'position', [200 200 600 450]);
    subplot(2,2,1); % Perception, Control Group
    hMarker(1) = scatter(M_ratio.group_1.session_01.perception.trained.raw, M_ratio.group_1.session_10.perception.trained.raw, 100, 'o', 'markerfacecolor', [.5 0 0], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(M_ratio.group_1.session_01.perception.untrained.raw, M_ratio.group_1.session_10.perception.untrained.raw, 100, 'o', 'markerfacecolor', [1 0 0], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([0 4]); ylim([0 4]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('CG: Perception meta-d''/d''');
    leg = legend('Trained', 'Untrained', 'location', 'nw');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    subplot(2,2,2); % Perception, Experimental Group
    hMarker(1) = scatter(M_ratio.group_2.session_01.perception.trained.raw, M_ratio.group_2.session_10.perception.trained.raw, 100, 'd', 'markerfacecolor', [.5 0 0], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(M_ratio.group_2.session_01.perception.untrained.raw, M_ratio.group_2.session_10.perception.untrained.raw, 100, 'd', 'markerfacecolor', [1 0 0], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([0 4]); ylim([0 4]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('EG: Perception meta-d''/d''');
    leg = legend('Trained', 'Untrained', 'location', 'se');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    subplot(2,2,3); % Memory, Control Group
    hMarker(1) = scatter(M_ratio.group_1.session_01.memory.trained.raw, M_ratio.group_1.session_10.memory.trained.raw, 100, 'o', 'markerfacecolor', [0 0 0.5], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(M_ratio.group_1.session_01.memory.untrained.raw, M_ratio.group_1.session_10.memory.untrained.raw, 100, 'o', 'markerfacecolor', [0 0 1], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([0 4]); ylim([0 4]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('CG: Memory meta-d''/d''');
    leg = legend('Trained', 'Untrained', 'location', 'nw');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    subplot(2,2,4); % Memory, Experimental Group
    hMarker(1) = scatter(M_ratio.group_2.session_01.memory.trained.raw, M_ratio.group_2.session_10.memory.trained.raw, 100, 'd', 'markerfacecolor', [0 0 .5], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(M_ratio.group_2.session_01.memory.untrained.raw, M_ratio.group_2.session_10.memory.untrained.raw, 100, 'd', 'markerfacecolor', [0 0 1], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([0 4]); ylim([0 4]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('EG: Memory meta-d''/d''');
    leg = legend('Trained', 'Untrained', 'location', 'se');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    if exportFigs
        export_fig MRatioScatterPlot -png -transparent 'MRatioScatter.png';
    end
    
    MRatioTransferScatterPlot = figure;
    set(gcf, 'position', [200 200 600 450]);
    subplot(2,2,1); % Perception, Trained
    hMarker(1) = scatter(M_ratio.group_2.session_01.perception.trained.raw(M_ratio.group_2.session_10.perception.trained.raw > M_ratio.group_2.session_01.perception.trained.raw), M_ratio.group_2.session_10.perception.trained.raw(M_ratio.group_2.session_10.perception.trained.raw > M_ratio.group_2.session_01.perception.trained.raw), 100, 'd', 'markerfacecolor', [.5 0 0], 'markeredgecolor', 'k'); hold on;
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([0 3.25]); ylim([0 3.25]);
    set(gca, 'fontsize', 14);
    ylabel('Post meta-d''/d''');
    title('EG: Perception, trained');
    box off;
    subplot(2,2,2); % Perception, Untrained
    hMarker(2) = scatter(M_ratio.group_2.session_01.perception.untrained.raw(M_ratio.group_2.session_10.perception.trained.raw > M_ratio.group_2.session_01.perception.trained.raw), M_ratio.group_2.session_10.perception.untrained.raw(M_ratio.group_2.session_10.perception.trained.raw > M_ratio.group_2.session_01.perception.trained.raw), 100, 'd', 'markerfacecolor', [1 0 0], 'markeredgecolor', 'k'); hold on;
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([0 3.25]); ylim([0 3.25]);
    set(gca, 'fontsize', 14);
    title('EG: Perception, untrained');
    box off;
    subplot(2,2,3); % Memory, Trained
    hMarker(3) = scatter(M_ratio.group_2.session_01.memory.trained.raw(M_ratio.group_2.session_10.perception.trained.raw > M_ratio.group_2.session_01.perception.trained.raw), M_ratio.group_2.session_10.memory.trained.raw(M_ratio.group_2.session_10.perception.trained.raw > M_ratio.group_2.session_01.perception.trained.raw), 100, 'd', 'markerfacecolor', [0 0 0.5], 'markeredgecolor', 'k'); hold on;
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([0 3.25]); ylim([0 3.25]);
    set(gca, 'fontsize', 14);
    xlabel('Pre meta-d''/d'''); ylabel('Post meta-d''/d''');
    title('EG: Memory, trained');
    box off;
    subplot(2,2,4); % Memory, Untrained
    hMarker(4) = scatter(M_ratio.group_2.session_01.memory.untrained.raw(M_ratio.group_2.session_10.perception.trained.raw > M_ratio.group_2.session_01.perception.trained.raw), M_ratio.group_2.session_10.memory.untrained.raw(M_ratio.group_2.session_10.perception.trained.raw > M_ratio.group_2.session_01.perception.trained.raw), 100, 'd', 'markerfacecolor', [0 0 1], 'markeredgecolor', 'k'); hold on;
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([0 3.25]); ylim([0 3.25]);
    set(gca, 'fontsize', 14);
    xlabel('Pre meta-d''/d'''); 
    title('EG: Memory, untrained');
    box off;
    if exportFigs
        export_fig MRatioTransferScatterPlot -png -transparent MRatioTransferScatter.png;
    end

    
    MCurvePlot = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    hBar(1) = errorbar(.6, M_ratio.group_1.session_01.perception.trained.mean, M_ratio.group_1.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, M_ratio.group_1.session_01.perception.untrained.mean, M_ratio.group_1.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, M_ratio.group_1.session_01.memory.trained.mean, M_ratio.group_1.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, M_ratio.group_1.session_01.memory.untrained.mean, M_ratio.group_1.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.M_ratio.group_1.learningCurve.mean, plots.M_ratio.group_1.learningCurve.sem,'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, M_ratio.group_1.session_10.perception.trained.mean, M_ratio.group_1.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, M_ratio.group_1.session_10.perception.untrained.mean, M_ratio.group_1.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, M_ratio.group_1.session_10.memory.trained.mean, M_ratio.group_1.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, M_ratio.group_1.session_10.memory.untrained.mean, M_ratio.group_1.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    ylim([.8, 2.5]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('meta-d''/d''', 'fontsize', 14);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    hBar(1) = errorbar(.6, M_ratio.group_2.session_01.perception.trained.mean, M_ratio.group_2.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, M_ratio.group_2.session_01.perception.untrained.mean, M_ratio.group_2.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, M_ratio.group_2.session_01.memory.trained.mean, M_ratio.group_2.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, M_ratio.group_2.session_01.memory.untrained.mean, M_ratio.group_2.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.M_ratio.group_2.learningCurve.mean, plots.M_ratio.group_2.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, M_ratio.group_2.session_10.perception.trained.mean, M_ratio.group_2.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, M_ratio.group_2.session_10.perception.untrained.mean, M_ratio.group_2.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, M_ratio.group_2.session_10.memory.trained.mean, M_ratio.group_2.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, M_ratio.group_2.session_10.memory.untrained.mean, M_ratio.group_2.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    ylim([.8, 2.5]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('meta-d''/d''', 'fontsize', 14);
    title('Experimental Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    if exportFigs
        export_fig MCurvePlot -png -transparent 'MCurve.png';
    end
    
end
       
end

