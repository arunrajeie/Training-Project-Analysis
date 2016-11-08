function [ pC ] = groupPercentCorrect( results, plotFigs, exportFigs )
%GROUPPERCENTCORRECT Runs group percent correct analysis and optionally plots and exports the figures

dom = {'perception', 'memory'};
stim = {'abstract', 'words'};
tStim = {'trained', 'untrained'};
for d = 1:numel(dom)
    for s = 1:numel(stim)
        pC.(dom{d}).(stim{s}).raw = [];
    end
end
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    for sesh = 1:10
        session = sprintf('session_%.2d',sesh);
        if sesh == 1 || sesh == 10
            for d = 1:numel(dom)
                for  s = 1:numel(tStim)
                    pC.(groups{g}).(session).(dom{d}).(tStim{s}).raw = [];
                end
            end
        else % Sessions 2-9
            pC.(groups{g}).(session).perception.trained.raw = [];
        end
    end
end

subjects = fieldnames(results);
% Concatenate raw data
for sub = 1:numel(subjects)
    group = sprintf('group_%d', results.(subjects{sub}).group);
    sessions = fieldnames(results.(subjects{sub}));
    pC.(group).(subjects{sub}).raw = [];
    for sesh = 1:numel(sessions)
        if strncmp(sessions{sesh},'session',7)
            session = str2double(sessions{sesh}(end-1:end));
            if session == 1 || session == 10
                for d = 1:numel(dom)
                    for s = 1:numel(stim)
                        pC.(group).(sessions{sesh}).(dom{d}).(tStim{s}).raw = vertcat(pC.(group).(sessions{sesh}).(dom{d}).(tStim{s}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(tStim{s}).percentCorrect);
                        pC.(group).(subjects{sub}).raw = vertcat(pC.(group).(subjects{sub}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(tStim{s}).percentCorrect);
                        pC.(dom{d}).(stim{s}).raw = vertcat(pC.(dom{d}).(stim{s}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).percentCorrect);
                    end
                end
            else % Sessions 2-9
                pC.(group).(subjects{sub}).raw = vertcat(pC.(group).(subjects{sub}).raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.percentCorrect);
                pC.(group).(sessions{sesh}).perception.trained.raw = vertcat(pC.(group).(sessions{sesh}).perception.trained.raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.percentCorrect);
                if strcmp(results.(subjects{sub}).stimCond, 'abstract')
                    pC.perception.abstract.raw = vertcat(pC.perception.abstract.raw, results.(subjects{sub}).(sessions{sesh}).perception.abstract.percentCorrect);
                elseif strcmp(results.(subjects{sub}).stimCond, 'words')
                    pC.perception.words.raw = vertcat(pC.perception.words.raw, results.(subjects{sub}).(sessions{sesh}).perception.words.percentCorrect);
                end
            end
        end
    end    
end

% Take mean and standard error
for d = 1:numel(dom)
    for s = 1:numel(stim)
        pC.(dom{d}).(stim{s}).mean = nanmean(pC.(dom{d}).(stim{s}).raw);
        pC.(dom{d}).(stim{s}).sem = nanstd(pC.(dom{d}).(stim{s}).raw)/sqrt(length(pC.(dom{d}).(stim{s}).raw));
    end
end

% Take mean and standard error
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    [plots.pC.(groups{g}).learningCurve.mean, plots.pC.(groups{g}).learningCurve.sem] = deal([]);
    pC.(groups{g}).raw = [];
    subjects = fieldnames(pC.(groups{g}));
    for sub = 1:numel(subjects)
        if strncmp(subjects{sub},'subject',7)
            pC.(groups{g}).(subjects{sub}).mean = nanmean(pC.(groups{g}).(subjects{sub}).raw);
            pC.(groups{g}).raw = vertcat(pC.(groups{g}).raw, pC.(groups{g}).(subjects{sub}).mean);
        end
    end
    pC.(groups{g}).mean = nanmean(pC.(groups{g}).raw);
    pC.(groups{g}).sem = nanstd(pC.(groups{g}).raw)/sqrt(length(pC.(groups{g}).raw));
    for sesh = 1:10
        session = sprintf('session_%.2d',sesh);
        if sesh == 1 || sesh == 10
            for d = 1:numel(dom)
                for s = 1:numel(tStim)
                    pC.(groups{g}).(session).(dom{d}).(tStim{s}).mean = nanmean(pC.(groups{g}).(session).(dom{d}).(tStim{s}).raw);
                    pC.(groups{g}).(session).(dom{d}).(tStim{s}).sem = nanstd(pC.(groups{g}).(session).(dom{d}).(tStim{s}).raw)/sqrt(length(pC.(groups{g}).(session).(dom{d}).(tStim{s}).raw));
                end
            end
            plots.pC.(groups{g}).(session).mean = [pC.(groups{g}).(session).perception.trained.mean, pC.(groups{g}).(session).perception.untrained.mean;...
                pC.(groups{g}).(session).memory.trained.mean, pC.(groups{g}).(session).memory.untrained.mean];
            plots.pC.(groups{g}).(session).sem = [pC.(groups{g}).(session).perception.trained.sem, pC.(groups{g}).(session).perception.untrained.sem;...
                pC.(groups{g}).(session).memory.trained.sem, pC.(groups{g}).(session).memory.untrained.sem];
        else % sessions 2-9
            pC.(groups{g}).(session).perception.trained.mean = nanmean(pC.(groups{g}).(session).perception.trained.raw);
            pC.(groups{g}).(session).perception.trained.sem = nanstd(pC.(groups{g}).(session).perception.trained.raw)/sqrt(length(pC.(groups{g}).(session).perception.trained.raw));
            plots.pC.(groups{g}).learningCurve.mean = vertcat(plots.pC.(groups{g}).learningCurve.mean, pC.(groups{g}).(session).perception.trained.mean);
            plots.pC.(groups{g}).learningCurve.sem = vertcat(plots.pC.(groups{g}).learningCurve.sem, pC.(groups{g}).(session).perception.trained.sem);
        end
    end    
end

[pC.h, pC.p, pC.ci, pC.stats] = ttest2(pC.group_1.raw, pC.group_2.raw);

if plotFigs
    % Error Bar Comparison Plot
    pCComparisonPlot = figure;
    set(gcf,'position', [200 200 450 300]);
    subplot(1,2,1); % Control Group
    [hBar hErrorbar] = barwitherr([[plots.pC.group_1.session_01.sem(1,:), plots.pC.group_1.session_01.sem(2,:)]',[plots.pC.group_1.session_10.sem(1,:), plots.pC.group_1.session_10.sem(2,:)]']' .* 100,...
        [[plots.pC.group_1.session_01.mean(1,:), plots.pC.group_1.session_01.mean(2,:)]',[plots.pC.group_1.session_10.mean(1,:), plots.pC.group_1.session_10.mean(2,:)]']' .* 100);
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    ylim([65 85]);
    set(gca, 'fontsize', 14);
    ylabel('Percent Correct', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 6); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    [hBar hErrorbar] = barwitherr([[plots.pC.group_2.session_01.sem(1,:), plots.pC.group_2.session_01.sem(2,:)]',[plots.pC.group_2.session_10.sem(1,:), plots.pC.group_2.session_10.sem(2,:)]']' .* 100,...
        [[plots.pC.group_2.session_01.mean(1,:), plots.pC.group_2.session_01.mean(2,:)]',[plots.pC.group_2.session_10.mean(1,:), plots.pC.group_2.session_10.mean(2,:)]']' .* 100);
    set(hBar(1),'FaceColor',[.5 0 0]);
    set(hBar(2),'FaceColor',[1 0 0]);
    set(hBar(3),'FaceColor',[0 0 .5]);
    set(hBar(4),'FaceColor',[0 0 1]);
    ylim([65 85]);
    set(gca, 'fontsize', 14);
    ylabel('Percent Correct', 'fontsize', 14);
    set(gca,'xticklabel', {'Pre', 'Post'}, 'fontsize', 11);
    title('Experimental Group', 'fontsize', 14);
    box off;
    if exportFigs
        export_fig pCComparisonPlot -png -transparent 'percentCorrectComparison.png';
    end
    
    pCCurvePlot = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    hBar(1) = errorbar(.6, pC.group_1.session_01.perception.trained.mean .* 100, pC.group_1.session_01.perception.trained.sem .* 100,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, pC.group_1.session_01.perception.untrained.mean .* 100, pC.group_1.session_01.perception.untrained.sem .* 100,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, pC.group_1.session_01.memory.trained.mean .* 100, pC.group_1.session_01.memory.trained.sem .* 100,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, pC.group_1.session_01.memory.untrained.mean .* 100, pC.group_1.session_01.memory.untrained.sem .* 100,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.pC.group_1.learningCurve.mean .* 100, plots.pC.group_1.learningCurve.sem .* 100, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, pC.group_1.session_10.perception.trained.mean .* 100, pC.group_1.session_10.perception.trained.sem .* 100,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, pC.group_1.session_10.perception.untrained.mean .* 100, pC.group_1.session_10.perception.untrained.sem .* 100,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, pC.group_1.session_10.memory.trained.mean .* 100, pC.group_1.session_10.memory.trained.sem .* 100,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, pC.group_1.session_10.memory.untrained.mean .* 100, pC.group_1.session_10.memory.untrained.sem .* 100,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    ylim([65 85]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('Percent Correct', 'fontsize', 14);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    hBar(1) = errorbar(.6, pC.group_2.session_01.perception.trained.mean .* 100, pC.group_2.session_01.perception.trained.sem .* 100,'o'); hold on; % P, TS
    hBar(2) = errorbar(.8, pC.group_2.session_01.perception.untrained.mean .* 100, pC.group_2.session_01.perception.untrained.sem .* 100,'o'); hold on; % P, UTS
    hBar(3) = errorbar(1, pC.group_2.session_01.memory.trained.mean .* 100, pC.group_2.session_01.memory.trained.sem .* 100,'o'); hold on; % M, TS
    hBar(4) = errorbar(1.2, pC.group_2.session_01.memory.untrained.mean .* 100, pC.group_2.session_01.memory.untrained.sem .* 100,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    hCurve = errorbar(2:9, plots.pC.group_2.learningCurve.mean .* 100, plots.pC.group_2.learningCurve.sem .* 100, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(9.8, pC.group_2.session_10.perception.trained.mean .* 100, pC.group_2.session_10.perception.trained.sem .* 100,'o'); hold on; % P, TS
    hBar(2) = errorbar(10, pC.group_2.session_10.perception.untrained.mean .* 100, pC.group_2.session_10.perception.untrained.sem .* 100,'o'); hold on; % P, UTS
    hBar(3) = errorbar(10.2, pC.group_2.session_10.memory.trained.mean .* 100, pC.group_2.session_10.memory.trained.sem .* 100,'o'); hold on; % M, TS
    hBar(4) = errorbar(10.4, pC.group_2.session_10.memory.untrained.mean .* 100, pC.group_2.session_10.memory.untrained.sem .* 100,'o'); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [1 0 0], 'markerfacecolor', [1 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 1], 'markerfacecolor', [0 0 1], 'markersize', 10);
    ylim([65 85]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', 1:10, 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlabel('Session', 'fontsize', 14);
    ylabel('Percent Correct', 'fontsize', 14);
    title('Experimental Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained', 'M, trained', 'M, untrained', 'location', 'nw');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    if exportFigs
        export_fig pCCurvePlot -png -transparent 'percentCorrectCurve.png';
    end
    
end

end