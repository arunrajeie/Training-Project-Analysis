function [ diff ] = groupDifficulty( results, plotFigs, exportFigs )
%GROUPDIFFICULTY Runs group difficulty analysis and optionally plots and exports the figures

stim = {'trained', 'untrained'};
subjects = fieldnames(results);
% Initialize arrays
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    for sesh = 1:10
        session = sprintf('session_%.2d', sesh);
        if sesh == 1 || sesh == 10
            for s = 1:numel(stim)
                diff.(groups{g}).(session).perception.(stim{s}).raw = [];
            end
        else  % sessions 2-9
            diff.(groups{g}).(session).perception.trained.raw = [];
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
                for s = 1:numel(stim)
                    diff.(group).(sessions{sesh}).perception.(stim{s}).raw = vertcat(diff.(group).(sessions{sesh}).perception.(stim{s}).raw, results.(subjects{sub}).(sessions{sesh}).perception.(stim{s}).meanDifficulty);
                end
            else % Sessions 2-9
                diff.(group).(sessions{sesh}).perception.trained.raw = vertcat(diff.(group).(sessions{sesh}).perception.trained.raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.meanDifficulty);
            end
        end
    end    
end

% Take mean and standard error
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    [plots.diff.(groups{g}).learningCurve.mean, plots.diff.(groups{g}).learningCurve.sem] = deal([]);
    for sesh = 1:10
        session = sprintf('session_%.2d',sesh);
        if sesh == 1 || sesh == 10
            for s = 1:numel(stim)
                diff.(groups{g}).(session).perception.(stim{s}).mean = nanmean(diff.(groups{g}).(session).perception.(stim{s}).raw);
                diff.(groups{g}).(session).perception.(stim{s}).sem = nanstd(diff.(groups{g}).(session).perception.(stim{s}).raw)/sqrt(length(diff.(groups{g}).(session).perception.(stim{s}).raw));
            end
            plots.diff.(groups{g}).(session).mean = [diff.(groups{g}).(session).perception.trained.mean, diff.(groups{g}).(session).perception.untrained.mean];
            plots.diff.(groups{g}).(session).sem = [diff.(groups{g}).(session).perception.trained.sem, diff.(groups{g}).(session).perception.untrained.sem];
        else % sessions 2-9
            diff.(groups{g}).(session).perception.trained.mean = nanmean(diff.(groups{g}).(session).perception.trained.raw);
            diff.(groups{g}).(session).perception.trained.sem = nanstd(diff.(groups{g}).(session).perception.trained.raw)/sqrt(length(diff.(groups{g}).(session).perception.trained.raw));
            plots.diff.(groups{g}).learningCurve.mean = vertcat(plots.diff.(groups{g}).learningCurve.mean, diff.(groups{g}).(session).perception.trained.mean);
            plots.diff.(groups{g}).learningCurve.sem = vertcat(plots.diff.(groups{g}).learningCurve.sem, diff.(groups{g}).(session).perception.trained.sem);
        end
    end
    
end

if plotFigs
    % Learning curve plot
    diffCurvePlot = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    hBar(1) = errorbar(.6, diff.group_1.session_01.perception.trained.mean, diff.group_1.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, diff.group_1.session_01.perception.untrained.mean, diff.group_1.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    hCurve = errorbar(2:9, plots.diff.group_1.learningCurve.mean, plots.diff.group_1.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, diff.group_1.session_10.perception.trained.mean, diff.group_1.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, diff.group_1.session_10.perception.untrained.mean, diff.group_1.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    ylim([114 124]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('Difficulty Level', 'fontsize', 14);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained','location', 'nw');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    hBar(1) = errorbar(.6, diff.group_2.session_01.perception.trained.mean, diff.group_2.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, diff.group_2.session_01.perception.untrained.mean, diff.group_2.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    hCurve = errorbar(2:9, plots.diff.group_2.learningCurve.mean, plots.diff.group_2.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, diff.group_2.session_10.perception.trained.mean, diff.group_2.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, diff.group_2.session_10.perception.untrained.mean, diff.group_2.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    ylim([114 124]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('Difficulty Level', 'fontsize', 14);
    title('Experimental Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'location', 'nw');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    if exportFigs
        export_fig diffCurvePlot -png -transparent 'difficultyCurve.png';
    end
    
end
end