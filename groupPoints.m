function [ points ] = groupPoints( data, plotFigs, exportFigs )
%GROUPPOINTS Runs group points analysis and optionally plots and exports
%the figures

subjects = fieldnames(data);
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    numSubjects.(groups{g}) = 0;
    for sesh = 2:9
        session = sprintf('session_%.2d', sesh);
        points.(groups{g}).(session).raw = [];
    end
end


for sub = 1:numel(subjects)
    group = sprintf('group_%d', data.(subjects{sub}).session_02.feedbackCond(1));
    numSubjects.(group) = numSubjects.(group) + 1;
    for sesh = 2:9
        session = sprintf('session_%.2d', sesh);
        points.(group).(session).raw = vertcat(points.(group).(session).raw, data.(subjects{sub}).(session).pointsEarned(end));
    end
end

groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    [plots.pointsCurve.(groups{g}).mean, plots.pointsCurve.(groups{g}).sem] = deal([]);
    for sesh = 2:9
        session = sprintf('session_%.2d', sesh);
        points.(groups{g}).(session).mean = nanmean(points.(groups{g}).(session).raw, 1);
        points.(groups{g}).(session).sem = nanstd(points.(groups{g}).(session).raw, [], 1)/sqrt(numSubjects.(groups{g}));
        plots.pointsCurve.(groups{g}).mean = vertcat(plots.pointsCurve.(groups{g}).mean, points.(groups{g}).(session).mean);
        plots.pointsCurve.(groups{g}).sem = vertcat(plots.pointsCurve.(groups{g}).sem, points.(groups{g}).(session).sem);
    end
end

if plotFigs 
    pointsFeedbackCurve = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    hCurve = errorbar(2:9, plots.pointsCurve.group_1.mean, plots.pointsCurve.group_1.sem,  'linewidth', 3);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    ylim([100000 210000]);
    set(gca, 'xtick', 2:9);
    set(gca,'fontsize',14);
    title('Control Group');
    xlabel('Session');
    ylabel('Points Feedback');
    subplot(1,2,2); % Experimental Group
    lineProps = struct('color', 'k', 'width', 2);
    hCurve = errorbar(2:9, plots.pointsCurve.group_2.mean, plots.pointsCurve.group_2.sem, 'linewidth', 3);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    ylim([100000 210000]);
    set(gca, 'xtick', 2:9);
    set(gca,'fontsize',14);
    title('Experimental Group');
    xlabel('Session');
    ylabel('Points Feedback');
    if exportFigs
        export_fig pointsFeedbackCurve -png -transparent 'pointsFeedbackCurve.png';
    end
end
end