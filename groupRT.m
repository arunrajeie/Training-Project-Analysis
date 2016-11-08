function [ rt ] = groupRT( results, plotFigs, exportFigs )
%GROUPRT Runs group RT analysis and optionally plots and exports the figures

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
                    rt.(groups{g}).(session).(dom{d}).(stim{s}).raw = [];
                end
            end
        else  % sessions 2-9
            rt.(groups{g}).(session).perception.trained.raw = [];
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
                        rt.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw = vertcat(rt.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).meanRT);
                    end
                end
            else % Sessions 2-9
                rt.(group).(sessions{sesh}).perception.trained.raw = vertcat(rt.(group).(sessions{sesh}).perception.trained.raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.meanRT);
            end
        end
    end    
end

% Take mean and standard error
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    [plots.rt.(groups{g}).learningCurve.mean, plots.rt.(groups{g}).learningCurve.sem] = deal([]);
    for sesh = 1:10
        session = sprintf('session_%.2d',sesh);
        if sesh == 1 || sesh == 10
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    rt.(groups{g}).(session).(dom{d}).(stim{s}).mean = nanmean(rt.(groups{g}).(session).(dom{d}).(stim{s}).raw);
                    rt.(groups{g}).(session).(dom{d}).(stim{s}).sem = nanstd(rt.(groups{g}).(session).(dom{d}).(stim{s}).raw)/sqrt(length(rt.(groups{g}).(session).(dom{d}).(stim{s}).raw));
                end
            end
            plots.rt.(groups{g}).(session).mean = [rt.(groups{g}).(session).perception.trained.mean, rt.(groups{g}).(session).perception.untrained.mean;...
                rt.(groups{g}).(session).memory.trained.mean, rt.(groups{g}).(session).memory.untrained.mean];
            plots.rt.(groups{g}).(session).sem = [rt.(groups{g}).(session).perception.trained.sem, rt.(groups{g}).(session).perception.untrained.sem;...
                rt.(groups{g}).(session).memory.trained.sem, rt.(groups{g}).(session).memory.untrained.sem];
        else % sessions 2-9
            rt.(groups{g}).(session).perception.trained.mean = nanmean(rt.(groups{g}).(session).perception.trained.raw);
            rt.(groups{g}).(session).perception.trained.sem = nanstd(rt.(groups{g}).(session).perception.trained.raw)/sqrt(length(rt.(groups{g}).(session).perception.trained.raw));
            plots.rt.(groups{g}).learningCurve.mean = vertcat(plots.rt.(groups{g}).learningCurve.mean, rt.(groups{g}).(session).perception.trained.mean);
            plots.rt.(groups{g}).learningCurve.sem = vertcat(plots.rt.(groups{g}).learningCurve.sem, rt.(groups{g}).(session).perception.trained.sem);
        end
    end
    
end

if plotFigs
    % Error Bar Comparison Plot
    rtComparisonPlot = figure;
    set(gcf,'position', [200 200 450 300]);
    subplot(1,2,1); % Control Group
    [hBar hErrorbar] = barwitherr([[plots.rt.group_1.session_01.sem(1,:), plots.rt.group_1.session_01.sem(2,:)]',[plots.rt.group_1.session_10.sem(1,:), plots.rt.group_1.session_10.sem(2,:)]']',...
        [[plots.rt.group_1.session_01.mean(1,:), plots.rt.group_1.session_01.mean(2,:)]',[plots.rt.group_1.session_10.mean(1,:), plots.rt.group_1.session_10.mean(2,:)]']');
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    set(gca, 'fontsize', 14);
    ylabel('rt', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 6); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    [hBar hErrorbar] = barwitherr([[plots.rt.group_2.session_01.sem(1,:), plots.rt.group_2.session_01.sem(2,:)]',[plots.rt.group_2.session_10.sem(1,:), plots.rt.group_2.session_10.sem(2,:)]']',...
        [[plots.rt.group_2.session_01.mean(1,:), plots.rt.group_2.session_01.mean(2,:)]',[plots.rt.group_2.session_10.mean(1,:), plots.rt.group_2.session_10.mean(2,:)]']');
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    set(gca, 'fontsize', 14);
    ylabel('rt', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Experimental Group', 'fontsize', 14);
    box off;
    if exportFigs
        export_fig rtComparisonPlot -png -transparent 'rtComparison.png';
    end
    
    rtCurvePlot = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    hBar(1) = errorbar(.6, rt.group_1.session_01.perception.trained.mean, rt.group_1.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, rt.group_1.session_01.perception.untrained.mean, rt.group_1.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, rt.group_1.session_01.memory.trained.mean, rt.group_1.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, rt.group_1.session_01.memory.untrained.mean, rt.group_1.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.rt.group_1.learningCurve.mean, plots.rt.group_1.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, rt.group_1.session_10.perception.trained.mean, rt.group_1.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, rt.group_1.session_10.perception.untrained.mean, rt.group_1.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, rt.group_1.session_10.memory.trained.mean, rt.group_1.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, rt.group_1.session_10.memory.untrained.mean, rt.group_1.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('rt', 'fontsize', 14);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'north');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    hBar(1) = errorbar(.6, rt.group_2.session_01.perception.trained.mean, rt.group_2.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, rt.group_2.session_01.perception.untrained.mean, rt.group_2.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, rt.group_2.session_01.memory.trained.mean, rt.group_2.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, rt.group_2.session_01.memory.untrained.mean, rt.group_2.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.rt.group_2.learningCurve.mean, plots.rt.group_2.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, rt.group_2.session_10.perception.trained.mean, rt.group_2.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, rt.group_2.session_10.perception.untrained.mean, rt.group_2.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, rt.group_2.session_10.memory.trained.mean, rt.group_2.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, rt.group_2.session_10.memory.untrained.mean, rt.group_2.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('rt', 'fontsize', 14);
    title('Experimental Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'north');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    if exportFigs
        export_fig rtCurvePlot -png -transparent 'rtCurve.png';
    end

end
end