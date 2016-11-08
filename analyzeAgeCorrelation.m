function [ ] = analyzeAgeCorrelation( results, plotFigs, exportFigs )
%ANALYZEAGECORRELATION Performs analysis on correlatation between age
%and log(meta-d'/d')(efficiency)


groups = {'group_1', 'group_2'};
dom = {'perception', 'memory'};
stim = {'trained', 'untrained'};
for g = 1:numel(groups)
    age.(groups{g}) = [];
    for d = 1:numel(dom)
        for s = 1:numel(stim)
            M_ratio_logDiff.(groups{g}).(dom{d}).(stim{s}) = [];
        end
    end
end

subjects = fieldnames(results);
for sub = 1:numel(subjects)
    group = sprintf('group_%d', results.(subjects{sub}).group);
    age.(group) = vertcat(age.(group), results.(subjects{sub}).age);
    for d = 1:numel(dom)
        for s = 1:numel(stim)
            M_ratio_logDiff.(subjects{sub}).(dom{d}).(stim{s}) = log(results.(subjects{sub}).session_10.(dom{d}).(stim{s}).fit.M_ratio) - log(results.(subjects{sub}).session_01.(dom{d}).(stim{s}).fit.M_ratio);
            M_ratio_logDiff.(group).(dom{d}).(stim{s}) = vertcat(M_ratio_logDiff.(group).(dom{d}).(stim{s}), M_ratio_logDiff.(subjects{sub}).(dom{d}).(stim{s}));
        end
    end
end

if plotFigs
    ageCorrelationScatterPlot = figure;
    set(gcf, 'position', [200 200 600 450]);
    subplot(2,2,1); % Perception, Control Group
    hMarker(1) = scatter(age.group_1, M_ratio_logDiff.group_1.perception.trained, 100, 'o', 'markerfacecolor', [.5 0 0], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(age.group_1, M_ratio_logDiff.group_1.perception.untrained, 100, 'o', 'markerfacecolor', [1 0 0], 'markeredgecolor', 'k'); 
        x = [age.group_1; age.group_1];
        y = [M_ratio_logDiff.group_1.perception.trained; M_ratio_logDiff.group_1.perception.untrained];
        X = [ones(size(x,1)), x];
        b = X\y;
        yCalc = X*b;
        plot(x, yCalc, '--k', 'linewidth', 2);
        [rho, p] = corr(x,y);
        t = text(45, -2, sprintf('\\rho = %.3f\np = %.3f', rho, p));
    xlim([0 65]); ylim([-3 3]);
    set(gca, 'fontsize', 14);
    ylabel('log(meta-d''/d'') diff');
%     xlabel('Mean Conf Diff'); ylabel('log(meta-d''/d'') diff');
    title('CG: Perception');
    leg = legend('Trained', 'Untrained', 'location', 'sw');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    subplot(2,2,2); % Perception, Experimental Group
    hMarker(1) = scatter(age.group_2, M_ratio_logDiff.group_2.perception.trained, 100, 'd', 'markerfacecolor', [.5 0 0], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(age.group_2, M_ratio_logDiff.group_2.perception.untrained, 100, 'd', 'markerfacecolor', [1 0 0], 'markeredgecolor', 'k'); 
        x = [age.group_2; age.group_2];
        y = [M_ratio_logDiff.group_2.perception.trained; M_ratio_logDiff.group_2.perception.untrained];
        X = [ones(size(x,1)), x];
        b = X\y;
        yCalc = X*b;
        plot(x, yCalc, '--k', 'linewidth', 2);
        [rho, p] = corr(x,y);
        t = text(45, -2, sprintf('\\rho = %.3f\np = %.3f', rho, p));
    xlim([0 65]); ylim([-3 3]);
    set(gca, 'fontsize', 14);
%     xlabel('Mean Conf Diff'); ylabel('log(meta-d''/d'') diff');
    title('EG: Perception');
    leg = legend('Trained', 'Untrained', 'location', 'sw');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    subplot(2,2,3); % Memory, Control Group
    hMarker(1) = scatter(age.group_1, M_ratio_logDiff.group_1.memory.trained, 100, 'o', 'markerfacecolor', [0 0 .5], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(age.group_1, M_ratio_logDiff.group_1.memory.untrained, 100, 'o', 'markerfacecolor', [0 0 1], 'markeredgecolor', 'k'); 
        x = [age.group_1; age.group_1];
        y = [M_ratio_logDiff.group_1.memory.trained; M_ratio_logDiff.group_1.memory.untrained];
        X = [ones(size(x,1)), x];
        b = X\y;
        yCalc = X*b;
        plot(x, yCalc, '--k', 'linewidth', 2);
        [rho, p] = corr(x,y);
        t = text(45, -2, sprintf('\\rho = %.3f\np = %.3f', rho, p));
    xlim([0 65]); ylim([-3 3]);
    set(gca, 'fontsize', 14);
    xlabel('Age'); ylabel('log(meta-d''/d'') diff');
    title('CG: Memory');
    leg = legend('Trained', 'Untrained', 'location', 'sw');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    subplot(2,2,4); % Memory, Experimental Group
    hMarker(1) = scatter(age.group_2, M_ratio_logDiff.group_2.memory.trained, 100, 'd', 'markerfacecolor', [0 0 .5], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(age.group_2, M_ratio_logDiff.group_2.memory.untrained, 100, 'd', 'markerfacecolor', [0 0 1], 'markeredgecolor', 'k'); 
        x = [age.group_2; age.group_2];
        y = [M_ratio_logDiff.group_2.memory.trained; M_ratio_logDiff.group_2.memory.untrained];
        X = [ones(size(x,1)), x];
        b = X\y;
        yCalc = X*b;
        plot(x, yCalc, '--k', 'linewidth', 2);
        [rho, p] = corr(x,real(y));
        t = text(45, -2, sprintf('\\rho = %.3f\np = %.3f', rho, p));
    xlim([0 65]); ylim([-3 3]);
    set(gca, 'fontsize', 14);
    xlabel('Age'); % ylabel('log(meta-d''/d'') diff');
    title('EG: Memory');
    leg = legend('Trained', 'Untrained', 'location', 'sw');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    if exportFigs
        export_fig ageCorrelationScatterPlot -png -transparent 'ageCorrelationScatter.png';
    end

    ageCorrelationGroup2ScatterPlot = figure;
    set(gcf, 'position', [200 200 600 450]);
    subplot(2,2,1); % Perception, Trained
    hMarker(1) = scatter(age.group_2, M_ratio_logDiff.group_2.perception.trained, 100, 'd', 'markerfacecolor', [.5 0 0], 'markeredgecolor', 'k'); hold on;
        x = age.group_2;
        y = M_ratio_logDiff.group_2.perception.trained;
        X = [ones(size(x,1)), x];
        b = X\y;
        yCalc = X*b;
        plot(x, yCalc, '--k', 'linewidth', 2);
        [rho, p] = corr(x,y);
        t = text(45, -2, sprintf('\\rho = %.3f\np = %.3f', rho, p));
    xlim([0 65]); ylim([-3 3]);
    set(gca, 'fontsize', 14);
    ylabel('log(meta-d''/d'') diff');
%     xlabel('Mean Conf Diff'); ylabel('log(meta-d''/d'') diff');
    title('EG: Perception, trained');
    box off;
    subplot(2,2,2); % Perception, Experimental Group
    hMarker(2) = scatter(age.group_2, M_ratio_logDiff.group_2.perception.untrained, 100, 'd', 'markerfacecolor', [1 0 0], 'markeredgecolor', 'k'); hold on;
        x = age.group_2;
        y = M_ratio_logDiff.group_2.perception.untrained;
        X = [ones(size(x,1)), x];
        b = X\y;
        yCalc = X*b;
        plot(x, yCalc, '--k', 'linewidth', 2);
        [rho, p] = corr(x,y);
        t = text(45, -2, sprintf('\\rho = %.3f\np = %.3f', rho, p));
    xlim([0 65]); ylim([-3 3]);
    set(gca, 'fontsize', 14);
%     xlabel('Mean Conf Diff'); ylabel('log(meta-d''/d'') diff');
    title('EG: Perception, untrained');
    box off;
    subplot(2,2,3); % Memory, Control Group
    hMarker(3) = scatter(age.group_2, M_ratio_logDiff.group_2.memory.trained, 100, 'd', 'markerfacecolor', [0 0 .5], 'markeredgecolor', 'k'); hold on;
        x = age.group_2;
        y = M_ratio_logDiff.group_2.memory.trained;
        X = [ones(size(x,1)), x];
        b = X\y;
        yCalc = X*b;
        plot(x, yCalc, '--k', 'linewidth', 2);
        [rho, p] = corr(x,real(y));
        t = text(45, -2, sprintf('\\rho = %.3f\np = %.3f', rho, p));
    xlim([0 65]); ylim([-3 3]);
    set(gca, 'fontsize', 14);
    xlabel('Age'); ylabel('log(meta-d''/d'') diff');
    title('EG: Memory, trained');
    box off;
    subplot(2,2,4); % Memory, Experimental Group
    hMarker(4) = scatter(age.group_2, M_ratio_logDiff.group_2.memory.untrained, 100, 'd', 'markerfacecolor', [0 0 1], 'markeredgecolor', 'k'); hold on;
        x = age.group_2;
        y = M_ratio_logDiff.group_2.memory.untrained;
        X = [ones(size(x,1)), x];
        b = X\y;
        yCalc = X*b;
        plot(x, yCalc, '--k', 'linewidth', 2);
        [rho, p] = corr(x,real(y));
        t = text(45, -2, sprintf('\\rho = %.3f\np = %.3f', rho, p));
    xlim([0 65]); ylim([-3 3]);
    set(gca, 'fontsize', 14);
    xlabel('Age'); % ylabel('log(meta-d''/d'') diff');
    title('EG: Memory, untrained');
    box off;
    if exportFigs
        export_fig ageCorrelationGroup2ScatterPlot -png -transparent 'ageCorrelationGroup2Scatter.png';
    end

end
end