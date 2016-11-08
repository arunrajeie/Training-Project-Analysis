function [ criteria ] = groupCriteria( results, plotFigs, exportFigs )
%GROUPCRITERIA Runs group criteria analysis and optionally plots and exports the figures

% Initialize arrays
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    for sesh = 1:10
        session = sprintf('session_%.2d', sesh);
        criteria.(groups{g}).(session).meta_da.raw = [];
        criteria.(groups{g}).(session).meta_ca.raw = [];
        criteria.(groups{g}).(session).t2ca_rS1.raw = [];
        criteria.(groups{g}).(session).t2ca_rS2.raw = [];
    end
end

dom = {'perception', 'memory'};
stim = {'trained', 'untrained'};
subjects = fieldnames(results);
[numSubjects.group_1, numSubjects.group_2] = deal(0);
% Concatenate raw data
for sub = 1:numel(subjects)
    if strncmp(subjects{sub},'subject',7)
        group = sprintf('group_%d', results.(subjects{sub}).group);
        numSubjects.(group) = numSubjects.(group) + 1;
        for sesh = 1:10
            session = sprintf('session_%.2d', sesh);
            if sesh == 1 || sesh == 10
                for d = 1:numel(dom)
                    for s = 1:numel(stim)
                        criteria.(group).(session).meta_da.raw = vertcat(criteria.(group).(session).meta_da.raw, results.(subjects{sub}).(session).(dom{d}).(stim{s}).fit.meta_da);
                        criteria.(group).(session).meta_ca.raw = vertcat(criteria.(group).(session).meta_ca.raw, results.(subjects{sub}).(session).(dom{d}).(stim{s}).fit.meta_ca);
                        criteria.(group).(session).t2ca_rS1.raw = vertcat(criteria.(group).(session).t2ca_rS1.raw, results.(subjects{sub}).(session).(dom{d}).(stim{s}).fit.t2ca_rS1); 
                        criteria.(group).(session).t2ca_rS2.raw = vertcat(criteria.(group).(session).t2ca_rS2.raw, results.(subjects{sub}).(session).(dom{d}).(stim{s}).fit.t2ca_rS2); 
                    end
                end
            else % Training sessions 2-9
                criteria.(group).(session).meta_da.raw = vertcat(criteria.(group).(session).meta_da.raw, results.(subjects{sub}).(session).perception.trained.fit.meta_da);
                criteria.(group).(session).meta_ca.raw = vertcat(criteria.(group).(session).meta_ca.raw, results.(subjects{sub}).(session).perception.trained.fit.meta_ca);
                criteria.(group).(session).t2ca_rS1.raw = vertcat(criteria.(group).(session).t2ca_rS1.raw, results.(subjects{sub}).(session).perception.trained.fit.t2ca_rS1); 
                criteria.(group).(session).t2ca_rS2.raw = vertcat(criteria.(group).(session).t2ca_rS2.raw, results.(subjects{sub}).(session).perception.trained.fit.t2ca_rS2); 
            end
        end  
    end
end

% Take mean and standard error
groups = {'group_1', 'group_2'};
x = -3:.01:3;
for g = 1:numel(groups)
    [plotLearningCurve.(groups{g}).mean, plotLearningCurve.(groups{g}).sem] = deal([]);
    for sesh = 1:10
        session = sprintf('session_%.2d',sesh);
        criteria.(groups{g}).(session).meta_da.mean = nanmean(criteria.(groups{g}).(session).meta_da.raw,1);
        criteria.(groups{g}).(session).meta_da.sem = nanstd(criteria.(groups{g}).(session).meta_da.raw,0,1)/sqrt(numSubjects.(groups{g}));
        criteria.(groups{g}).(session).meta_ca.mean = nanmean(criteria.(groups{g}).(session).meta_ca.raw,1);
        criteria.(groups{g}).(session).meta_ca.sem = nanstd(criteria.(groups{g}).(session).meta_ca.raw,0,1)/sqrt(numSubjects.(groups{g}));
        criteria.(groups{g}).(session).t2ca_rS1.mean = nanmean(criteria.(groups{g}).(session).t2ca_rS1.raw,1);    
        criteria.(groups{g}).(session).t2ca_rS1.sem = nanstd(criteria.(groups{g}).(session).t2ca_rS1.raw,0,1)/sqrt(numSubjects.(groups{g}));
        criteria.(groups{g}).(session).t2ca_rS2.mean = nanmean(criteria.(groups{g}).(session).t2ca_rS2.raw,1);    
        criteria.(groups{g}).(session).t2ca_rS2.sem = nanstd(criteria.(groups{g}).(session).t2ca_rS2.raw,0,1)/sqrt(numSubjects.(groups{g}));
        plotCriteria.(groups{g}).(session).mean = [criteria.(groups{g}).(session).t2ca_rS1.mean, criteria.(groups{g}).(session).meta_ca.mean, criteria.(groups{g}).(session).t2ca_rS2.mean];
        plotCriteria.(groups{g}).(session).sem = [criteria.(groups{g}).(session).t2ca_rS1.sem, criteria.(groups{g}).(session).meta_ca.sem, criteria.(groups{g}).(session).t2ca_rS2.sem];
        plotLearningCurve.(groups{g}).mean = vertcat(plotLearningCurve.(groups{g}).mean, plotCriteria.(groups{g}).(session).mean);
        plotLearningCurve.(groups{g}).sem = vertcat(plotLearningCurve.(groups{g}).sem, plotCriteria.(groups{g}).(session).sem);
        plotS1.(groups{g}).(session) = normpdf(x, criteria.(groups{g}).(session).meta_ca.mean - criteria.(groups{g}).(session).meta_da.mean/2, 1);
        plotS2.(groups{g}).(session) = normpdf(x, criteria.(groups{g}).(session).meta_ca.mean + criteria.(groups{g}).(session).meta_da.mean/2, 1);
    end    
end

if plotFigs
    criteriaShiftPlot = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    H = errorbarxy(plotCriteria.group_1.session_01.mean, -ones(7,1), plotCriteria.group_1.session_01.sem, ones(7,1), {'ok', 'b', 'k'}); hold on;
    set(H.hMain, 'linewidth', 3);
    set(H.hErrorbar, 'linewidth', 3);
    H = errorbarxy(plotCriteria.group_1.session_10.mean, ones(7,1), plotCriteria.group_1.session_10.sem, ones(7,1), {'ok', 'b', 'k'});
    plot([-3 3], [0 0], '-k', 'linewidth', 5);
    H2 = plot(x, plotS1.group_1.session_01 - 1, '-b', x, plotS2.group_1.session_01 - 1, '-b', 'linewidth', 2);
    H2 = plot(x, plotS1.group_1.session_10 + 1, '-b', x, plotS2.group_1.session_10 + 1, '-b', 'linewidth', 2);
    set(H.hMain, 'linewidth', 3);
    set(H.hErrorbar, 'linewidth', 3);
    set(gca,'fontsize',14);
    title('Control Group'); 
    set(gca, 'ytick', [-1,1], 'yticklabel', {'Pre', 'Post'});
    xlabel('Criteria Locations'); ylabel('Session');
    xlim([-3, 3]); ylim([-2, 2]);
    subplot(1,2,2); % Experimental Group
    H = errorbarxy(plotCriteria.group_2.session_01.mean, -ones(7,1), plotCriteria.group_2.session_01.sem, ones(7,1), {'ok', 'b', 'k'}); hold on;
    set(H.hMain, 'linewidth', 3);
    set(H.hErrorbar, 'linewidth', 3);
    H = errorbarxy(plotCriteria.group_2.session_10.mean, ones(7,1), plotCriteria.group_2.session_10.sem, ones(7,1), {'ok', 'b', 'k'});
    plot([-3 3], [0 0], '-k', 'linewidth', 5);
    H2 = plot(x, plotS1.group_2.session_01 - 1, '-b', x, plotS2.group_2.session_01 - 1, '-b', 'linewidth', 2);
    H2 = plot(x, plotS1.group_2.session_10 + 1, '-b', x, plotS2.group_2.session_10 + 1, '-b', 'linewidth', 2);
    set(H.hMain, 'linewidth', 3);
    set(H.hErrorbar, 'linewidth', 3);
    xlim([-3, 3]); ylim([-2, 2]);
    set(gca,'fontsize',14);
    title('Experimental Group'); 
    set(gca, 'ytick', [-1,1], 'yticklabel', {'Pre', 'Post'});
    xlabel('Criteria Locations'); ylabel('Session');

    if exportFigs
        export_fig criteriaShiftPlot -png -transparent 'criteriaShift.png';
    end    
    
    criteriaCurvePlot = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    errorbar(1:10, plotLearningCurve.group_1.mean(:,4), plotLearningCurve.group_1.sem(:,4), '-k', 'linewidth', 2);hold on;
    errorbar(1:10, plotLearningCurve.group_1.mean(:,3), plotLearningCurve.group_1.sem(:,3), '--k', 'linewidth', 2);
    errorbar(1:10, plotLearningCurve.group_1.mean(:,5), plotLearningCurve.group_1.sem(:,5), '--k', 'linewidth', 2);
    errorbar(1:10, plotLearningCurve.group_1.mean(:,2), plotLearningCurve.group_1.sem(:,2), '-.k', 'linewidth', 2);
    errorbar(1:10, plotLearningCurve.group_1.mean(:,6), plotLearningCurve.group_1.sem(:,6), '-.k', 'linewidth', 2);
    errorbar(1:10, plotLearningCurve.group_1.mean(:,1), plotLearningCurve.group_1.sem(:,1), ':k', 'linewidth', 2); 
    errorbar(1:10, plotLearningCurve.group_1.mean(:,7), plotLearningCurve.group_1.sem(:,7), ':k', 'linewidth', 2);
    set(gca,'fontsize',14);
    title('Control Group');
    xlabel('Session'); ylabel('Criteria');
    set(gca,'xtick', 1:10, 'xticklabel', {'Pre', 2:9, 'Post'});
     xlim([0 13.5]); ylim([-2.5, 2.5]);
    legend('c_1', 'c_2', 'c_3', 'c_4', 'location', 'ne');
    legend boxoff; box off;
    subplot(1,2,2); % Experimental Group
    errorbar(1:10, plotLearningCurve.group_2.mean(:,4), plotLearningCurve.group_2.sem(:,4), '-k', 'linewidth', 2);hold on;
    errorbar(1:10, plotLearningCurve.group_2.mean(:,3), plotLearningCurve.group_2.sem(:,3), '--k', 'linewidth', 2);
    errorbar(1:10, plotLearningCurve.group_2.mean(:,5), plotLearningCurve.group_2.sem(:,5), '--k', 'linewidth', 2);
    errorbar(1:10, plotLearningCurve.group_2.mean(:,2), plotLearningCurve.group_2.sem(:,2), '-.k', 'linewidth', 2);
    errorbar(1:10, plotLearningCurve.group_2.mean(:,6), plotLearningCurve.group_2.sem(:,6), '-.k', 'linewidth', 2);
    errorbar(1:10, plotLearningCurve.group_2.mean(:,1), plotLearningCurve.group_2.sem(:,1), ':k', 'linewidth', 2); 
    errorbar(1:10, plotLearningCurve.group_2.mean(:,7), plotLearningCurve.group_2.sem(:,7), ':k', 'linewidth', 2);
    set(gca,'fontsize',14);
    title('Experimental Group');
    xlabel('Session'); ylabel('Criteria');
    set(gca,'xtick', 1:10, 'xticklabel', {'Pre', 2:9, 'Post'});
    xlim([0 13.5]); ylim([-2.5, 2.5]);
    legend('c_1', 'c_2', 'c_3', 'c_4', 'location', 'ne');
    legend boxoff; box off;
    
    if exportFigs
        export_fig criteriaCurvePlot -png -transparent 'criteriaCurve.png';
    end
end

end

