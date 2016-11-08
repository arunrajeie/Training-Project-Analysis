function [ meta_da ] = groupMetaDa( results, plotFigs, exportFigs )
%GROUPDA Runs group meta-d' analysis and optionally plots and exports the figures

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
                    meta_da.(groups{g}).(session).(dom{d}).(stim{s}).raw = [];
                end
            end
        else  % sessions 2-9
            meta_da.(groups{g}).(session).perception.trained.raw = [];
        end
    end
end

% Concatenate raw data
for sub = 1:numel(subjects)
    group = sprintf('group_%d', results.(subjects{sub}).group);
    sessions = fieldnames(results.(subjects{sub}));
    meta_da.(group).(subjects{sub}).raw = [];
    for sesh = 1:numel(sessions)
        if strncmp(sessions{sesh},'session',7)
            session = str2double(sessions{sesh}(end-1:end));
            if session == 1 || session == 10
                for d = 1:numel(dom)
                    for s = 1:numel(stim)
                        meta_da.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw = vertcat(meta_da.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).fit.meta_da);
                        meta_da.(group).(subjects{sub}).raw = vertcat(meta_da.(group).(subjects{sub}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).fit.meta_da);
                    end
                end
            else % Sessions 2-9
                meta_da.(group).(sessions{sesh}).perception.trained.raw = vertcat(meta_da.(group).(sessions{sesh}).perception.trained.raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.fit.meta_da);
                meta_da.(group).(subjects{sub}).raw = vertcat(meta_da.(group).(subjects{sub}).raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.fit.meta_da);
            end
        end
    end    
end

% Take mean and standard error
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    [plots.meta_da.(groups{g}).learningCurve.mean, plots.meta_da.(groups{g}).learningCurve.sem] = deal([]);
    meta_da.(groups{g}).raw = [];
    subjects = fieldnames(meta_da.(groups{g}));
    for sub = 1:numel(subjects)
        if strncmp(subjects{sub},'subject',7)
            meta_da.(groups{g}).(subjects{sub}).mean = nanmean(meta_da.(groups{g}).(subjects{sub}).raw);
            meta_da.(groups{g}).raw = vertcat(meta_da.(groups{g}).raw, meta_da.(groups{g}).(subjects{sub}).mean);
        end
    end
    meta_da.(groups{g}).mean = nanmean(meta_da.(groups{g}).raw);
    meta_da.(groups{g}).sem = nanstd(meta_da.(groups{g}).raw)/sqrt(length(meta_da.(groups{g}).raw));
    for sesh = 1:10
        session = sprintf('session_%.2d',sesh);
        if sesh == 1 || sesh == 10
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    meta_da.(groups{g}).(session).(dom{d}).(stim{s}).mean = nanmean(meta_da.(groups{g}).(session).(dom{d}).(stim{s}).raw);
                    meta_da.(groups{g}).(session).(dom{d}).(stim{s}).sem = nanstd(meta_da.(groups{g}).(session).(dom{d}).(stim{s}).raw)/sqrt(length(meta_da.(groups{g}).(session).(dom{d}).(stim{s}).raw));
                end
            end
            plots.meta_da.(groups{g}).(session).mean = [meta_da.(groups{g}).(session).perception.trained.mean, meta_da.(groups{g}).(session).perception.untrained.mean;...
                meta_da.(groups{g}).(session).memory.trained.mean, meta_da.(groups{g}).(session).memory.untrained.mean];
            plots.meta_da.(groups{g}).(session).sem = [meta_da.(groups{g}).(session).perception.trained.sem, meta_da.(groups{g}).(session).perception.untrained.sem;...
                meta_da.(groups{g}).(session).memory.trained.sem, meta_da.(groups{g}).(session).memory.untrained.sem];
        else % sessions 2-9
            meta_da.(groups{g}).(session).perception.trained.mean = nanmean(meta_da.(groups{g}).(session).perception.trained.raw);
            meta_da.(groups{g}).(session).perception.trained.sem = nanstd(meta_da.(groups{g}).(session).perception.trained.raw)/sqrt(length(meta_da.(groups{g}).(session).perception.trained.raw));
            plots.meta_da.(groups{g}).learningCurve.mean = vertcat(plots.meta_da.(groups{g}).learningCurve.mean, meta_da.(groups{g}).(session).perception.trained.mean);
            plots.meta_da.(groups{g}).learningCurve.sem = vertcat(plots.meta_da.(groups{g}).learningCurve.sem, meta_da.(groups{g}).(session).perception.trained.sem);
        end
    end
end

[meta_da.h, meta_da.p, meta_da.ci, meta_da.stats] = ttest2(meta_da.group_1.raw, meta_da.group_2.raw);

if plotFigs
    % Error Bar Comparison Plot
    metaDaComparisonPlot = figure;
    set(gcf,'position', [200 200 450 300]);
    subplot(1,2,1); % Control Group
    [hBar hErrorbar] = barwitherr([[plots.meta_da.group_1.session_01.sem(1,:), plots.meta_da.group_1.session_01.sem(2,:)]',[plots.meta_da.group_1.session_10.sem(1,:), plots.meta_da.group_1.session_10.sem(2,:)]']',...
        [[plots.meta_da.group_1.session_01.mean(1,:), plots.meta_da.group_1.session_01.mean(2,:)]',[plots.meta_da.group_1.session_10.mean(1,:), plots.meta_da.group_1.session_10.mean(2,:)]']');
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    ylim([0 3]);
    set(gca, 'fontsize', 14);
    ylabel('meta-d''', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 6); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    [hBar hErrorbar] = barwitherr([[plots.meta_da.group_2.session_01.sem(1,:), plots.meta_da.group_2.session_01.sem(2,:)]',[plots.meta_da.group_2.session_10.sem(1,:), plots.meta_da.group_2.session_10.sem(2,:)]']',...
        [[plots.meta_da.group_2.session_01.mean(1,:), plots.meta_da.group_2.session_01.mean(2,:)]',[plots.meta_da.group_2.session_10.mean(1,:), plots.meta_da.group_2.session_10.mean(2,:)]']');
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    ylim([0 3]);
    set(gca, 'fontsize', 14);
    ylabel('meta-d''', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Experimental Group', 'fontsize', 14);
    box off;
    if exportFigs
        export_fig metaDaComparisonPlot -png -transparent 'metaDaComparison.png';
    end
    
    metaDaCurvePlot = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    hBar(1) = errorbar(.6, meta_da.group_1.session_01.perception.trained.mean, meta_da.group_1.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, meta_da.group_1.session_01.perception.untrained.mean, meta_da.group_1.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, meta_da.group_1.session_01.memory.trained.mean, meta_da.group_1.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, meta_da.group_1.session_01.memory.untrained.mean, meta_da.group_1.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.meta_da.group_1.learningCurve.mean, plots.meta_da.group_1.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, meta_da.group_1.session_10.perception.trained.mean, meta_da.group_1.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, meta_da.group_1.session_10.perception.untrained.mean, meta_da.group_1.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, meta_da.group_1.session_10.memory.trained.mean, meta_da.group_1.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, meta_da.group_1.session_10.memory.untrained.mean, meta_da.group_1.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    ylim([1, 3]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('meta-d''', 'fontsize', 14);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    hBar(1) = errorbar(.6, meta_da.group_2.session_01.perception.trained.mean, meta_da.group_2.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, meta_da.group_2.session_01.perception.untrained.mean, meta_da.group_2.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, meta_da.group_2.session_01.memory.trained.mean, meta_da.group_2.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, meta_da.group_2.session_01.memory.untrained.mean, meta_da.group_2.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.meta_da.group_2.learningCurve.mean, plots.meta_da.group_2.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, meta_da.group_2.session_10.perception.trained.mean, meta_da.group_2.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, meta_da.group_2.session_10.perception.untrained.mean, meta_da.group_2.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, meta_da.group_2.session_10.memory.trained.mean, meta_da.group_2.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, meta_da.group_2.session_10.memory.untrained.mean, meta_da.group_2.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    ylim([1, 3]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('meta-d''', 'fontsize', 14);
    title('Experimental Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    if exportFigs
        export_fig metaDaCurvePlot -png -transparent 'metaDaCurve.png';
    end
    
end
end