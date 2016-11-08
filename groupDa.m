function [ da ] = groupDa( results, plotFigs, exportFigs )
%GROUPDA Runs group d' analysis and optionally plots and exports the figures

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
                    da.(groups{g}).(session).(dom{d}).(stim{s}).raw = [];
                end
            end
        else  % sessions 2-9
            da.(groups{g}).(session).perception.trained.raw = [];
        end
    end
end

% Concatenate raw data
for sub = 1:numel(subjects)
    group = sprintf('group_%d', results.(subjects{sub}).group);
    sessions = fieldnames(results.(subjects{sub}));
    da.(group).(subjects{sub}).raw = [];
    for sesh = 1:numel(sessions)
        if strncmp(sessions{sesh},'session',7)
            session = str2double(sessions{sesh}(end-1:end));
            if session == 1 || session == 10
                for d = 1:numel(dom)
                    for s = 1:numel(stim)
                        da.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw = vertcat(da.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).fit.da);
                        da.(group).(subjects{sub}).raw = vertcat(da.(group).(subjects{sub}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).fit.da);
                    end
                end
            else % Sessions 2-9
                da.(group).(sessions{sesh}).perception.trained.raw = vertcat(da.(group).(sessions{sesh}).perception.trained.raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.fit.da);
                da.(group).(subjects{sub}).raw = vertcat(da.(group).(subjects{sub}).raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.fit.da);
            end
        end
    end    
end

% Take mean and standard error
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    [plots.da.(groups{g}).learningCurve.mean, plots.da.(groups{g}).learningCurve.sem] = deal([]);
    da.(groups{g}).raw = [];
    subjects = fieldnames(da.(groups{g}));
    for sub = 1:numel(subjects)
        if strncmp(subjects{sub},'subject',7)
            da.(groups{g}).(subjects{sub}).mean = nanmean(da.(groups{g}).(subjects{sub}).raw);
            da.(groups{g}).raw = vertcat(da.(groups{g}).raw, da.(groups{g}).(subjects{sub}).mean);
        end
    end
    da.(groups{g}).mean = nanmean(da.(groups{g}).raw);
    da.(groups{g}).sem = nanstd(da.(groups{g}).raw)/sqrt(length(da.(groups{g}).raw));
    for sesh = 1:10
        session = sprintf('session_%.2d',sesh);
        if sesh == 1 || sesh == 10
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    da.(groups{g}).(session).(dom{d}).(stim{s}).mean = nanmean(da.(groups{g}).(session).(dom{d}).(stim{s}).raw);
                    da.(groups{g}).(session).(dom{d}).(stim{s}).sem = nanstd(da.(groups{g}).(session).(dom{d}).(stim{s}).raw)/sqrt(length(da.(groups{g}).(session).(dom{d}).(stim{s}).raw));
                end
            end
            plots.da.(groups{g}).(session).mean = [da.(groups{g}).(session).perception.trained.mean, da.(groups{g}).(session).perception.untrained.mean;...
                da.(groups{g}).(session).memory.trained.mean, da.(groups{g}).(session).memory.untrained.mean];
            plots.da.(groups{g}).(session).sem = [da.(groups{g}).(session).perception.trained.sem, da.(groups{g}).(session).perception.untrained.sem;...
                da.(groups{g}).(session).memory.trained.sem, da.(groups{g}).(session).memory.untrained.sem];
        else % sessions 2-9
            da.(groups{g}).(session).perception.trained.mean = nanmean(da.(groups{g}).(session).perception.trained.raw);
            da.(groups{g}).(session).perception.trained.sem = nanstd(da.(groups{g}).(session).perception.trained.raw)/sqrt(length(da.(groups{g}).(session).perception.trained.raw));
            plots.da.(groups{g}).learningCurve.mean = vertcat(plots.da.(groups{g}).learningCurve.mean, da.(groups{g}).(session).perception.trained.mean);
            plots.da.(groups{g}).learningCurve.sem = vertcat(plots.da.(groups{g}).learningCurve.sem, da.(groups{g}).(session).perception.trained.sem);
        end
    end
end

[da.h, da.p, da.ci, da.stats] = ttest2(da.group_1.raw, da.group_2.raw);

if plotFigs
%     Error Bar Comparison Plot
    daComparisonPlot = figure;
    set(gcf,'position', [200 200 450 300]);
    subplot(1,2,1); % Control Group
    [hBar hErrorbar] = barwitherr([[plots.da.group_1.session_01.sem(1,:), plots.da.group_1.session_01.sem(2,:)]',[plots.da.group_1.session_10.sem(1,:), plots.da.group_1.session_10.sem(2,:)]']',...
        [[plots.da.group_1.session_01.mean(1,:), plots.da.group_1.session_01.mean(2,:)]',[plots.da.group_1.session_10.mean(1,:), plots.da.group_1.session_10.mean(2,:)]']');
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    ylim([0 2.3]);
    set(gca, 'fontsize', 14);
    ylabel('d''', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 6); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    [hBar hErrorbar] = barwitherr([[plots.da.group_2.session_01.sem(1,:), plots.da.group_2.session_01.sem(2,:)]',[plots.da.group_2.session_10.sem(1,:), plots.da.group_2.session_10.sem(2,:)]']',...
        [[plots.da.group_2.session_01.mean(1,:), plots.da.group_2.session_01.mean(2,:)]',[plots.da.group_2.session_10.mean(1,:), plots.da.group_2.session_10.mean(2,:)]']');
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    ylim([0 2.3]);
    set(gca, 'fontsize', 14);
    ylabel('d''', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Experimental Group', 'fontsize', 14);
    box off;
    if exportFigs
        export_fig daComparisonPlot -png -transparent 'daComparison.png';
    end
    
    daCurvePlot = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    hBar(1) = errorbar(.6, da.group_1.session_01.perception.trained.mean, da.group_1.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, da.group_1.session_01.perception.untrained.mean, da.group_1.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, da.group_1.session_01.memory.trained.mean, da.group_1.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, da.group_1.session_01.memory.untrained.mean, da.group_1.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.da.group_1.learningCurve.mean, plots.da.group_1.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, da.group_1.session_10.perception.trained.mean, da.group_1.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, da.group_1.session_10.perception.untrained.mean, da.group_1.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, da.group_1.session_10.memory.trained.mean, da.group_1.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, da.group_1.session_10.memory.untrained.mean, da.group_1.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    ylim([1.1, 1.8]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('d''', 'fontsize', 14);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    hBar(1) = errorbar(.6, da.group_2.session_01.perception.trained.mean, da.group_2.session_01.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, da.group_2.session_01.perception.untrained.mean, da.group_2.session_01.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, da.group_2.session_01.memory.trained.mean, da.group_2.session_01.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, da.group_2.session_01.memory.untrained.mean, da.group_2.session_01.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.da.group_2.learningCurve.mean, plots.da.group_2.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, da.group_2.session_10.perception.trained.mean, da.group_2.session_10.perception.trained.sem,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, da.group_2.session_10.perception.untrained.mean, da.group_2.session_10.perception.untrained.sem,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, da.group_2.session_10.memory.trained.mean, da.group_2.session_10.memory.trained.sem,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, da.group_2.session_10.memory.untrained.mean, da.group_2.session_10.memory.untrained.sem,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    ylim([1.1, 1.8]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('d''', 'fontsize', 14);
    title('Experimental Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    if exportFigs
        export_fig daCurvePlot -png -transparent 'daCurve.png';
    end
     
end
end