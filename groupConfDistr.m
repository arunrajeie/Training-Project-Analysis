function [ confDistr ] = groupConfDistr( analysis, plotFigs, exportFigs )
%GROUPCONFDISTR Runs group confidence distribution analysis and optionally plots and exports the figures

% Initialize arrays
groups = {'group_1', 'group_2'};
dom = {'perception', 'memory'};
for g = 1:numel(groups)
    numSubjects.(groups{g}) = 0;
    for sesh = [1,10]
        session = sprintf('session_%.2d', sesh);
        confDistr.(groups{g}).(session).correct.raw = [];
        confDistr.(groups{g}).(session).incorrect.raw = [];
        for d = 1:numel(dom)
            confDistr.(groups{g}).(session).(dom{d}).correct.raw = [];
            confDistr.(groups{g}).(session).(dom{d}).incorrect.raw = [];
        end
    end
end

dom = {'perception', 'memory'};
stim = {'abstract', 'words'};
subjects = fieldnames(analysis);
% Concatenate raw data
for sub = 1:numel(subjects)
    if strncmp(subjects{sub},'subject',7)
        group = sprintf('group_%d', analysis.(subjects{sub}).group);
        numSubjects.(group) = numSubjects.(group) + 1;
        for sesh = [1,10]
            session = sprintf('session_%.2d', sesh);
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    confDistr.(group).(session).correct.raw = vertcat(confDistr.(group).(session).correct.raw, analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned); % ./ (analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned + analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).numTrialsBinned));
                    confDistr.(group).(session).incorrect.raw = vertcat(confDistr.(group).(session).incorrect.raw, (analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).numTrialsBinned - analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned)); % ./ (analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned + analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).numTrialsBinned));
                    confDistr.(group).(session).(dom{d}).correct.raw = vertcat(confDistr.(group).(session).(dom{d}).correct.raw, analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned); % ./ (analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned + analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).numTrialsBinned));
                    confDistr.(group).(session).(dom{d}).incorrect.raw = vertcat(confDistr.(group).(session).(dom{d}).incorrect.raw, (analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).numTrialsBinned - analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned)); % ./ (analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned + analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).numTrialsBinned));
                end
            end
        end  
    end
end

% Take mean and standard error
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    for sesh = [1,10]
        session = sprintf('session_%.2d',sesh);
        confDistr.(groups{g}).(session).correct.mean = nanmean(confDistr.(groups{g}).(session).correct.raw,1);
        confDistr.(groups{g}).(session).correct.sem = nanstd(confDistr.(groups{g}).(session).correct.raw,0,1)/sqrt(numSubjects.(groups{g}));
        confDistr.(groups{g}).(session).incorrect.mean = nanmean(confDistr.(groups{g}).(session).incorrect.raw,1);    
        confDistr.(groups{g}).(session).incorrect.sem = nanstd(confDistr.(groups{g}).(session).incorrect.raw,0,1)/sqrt(numSubjects.(groups{g}));
        for d = 1:numel(dom)
            confDistr.(groups{g}).(session).(dom{d}).correct.mean = nanmean(confDistr.(groups{g}).(session).(dom{d}).correct.raw,1);
            confDistr.(groups{g}).(session).(dom{d}).correct.sem = nanstd(confDistr.(groups{g}).(session).(dom{d}).correct.raw,0,1)/sqrt(numSubjects.(groups{g}));
            confDistr.(groups{g}).(session).(dom{d}).incorrect.mean = nanmean(confDistr.(groups{g}).(session).(dom{d}).incorrect.raw,1);    
            confDistr.(groups{g}).(session).(dom{d}).incorrect.sem = nanstd(confDistr.(groups{g}).(session).(dom{d}).incorrect.raw,0,1)/sqrt(numSubjects.(groups{g}));
            
            plotConfDistr.(groups{g}).(session).(dom{d}) = [confDistr.(groups{g}).(session).(dom{d}).incorrect.mean/sum(confDistr.(groups{g}).(session).(dom{d}).incorrect.mean + confDistr.(groups{g}).(session).(dom{d}).correct.mean); confDistr.(groups{g}).(session).(dom{d}).correct.mean/sum(confDistr.(groups{g}).(session).(dom{d}).incorrect.mean + confDistr.(groups{g}).(session).(dom{d}).correct.mean)]';
        end
    end
end

if plotFigs
    confDistrPlot = figure;
    set(gcf,'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    hBar = bar(1:4,[confDistr.group_1.session_01.incorrect.mean/sum(confDistr.group_1.session_01.incorrect.mean + confDistr.group_1.session_01.correct.mean); confDistr.group_1.session_01.correct.mean/sum(confDistr.group_1.session_01.incorrect.mean + confDistr.group_1.session_01.correct.mean)]'); hold on;
    hBar2 = bar(7:10,[confDistr.group_1.session_10.incorrect.mean/sum(confDistr.group_1.session_10.incorrect.mean + confDistr.group_1.session_10.correct.mean); confDistr.group_1.session_10.correct.mean/sum(confDistr.group_1.session_10.incorrect.mean + confDistr.group_1.session_10.correct.mean)]');
    set(hBar(1),'facecolor', 'r');
    set(hBar(2),'facecolor', 'g');
    set(hBar2(1),'facecolor', 'r');
    set(hBar2(2),'facecolor', 'g');
    set(gca, 'xtick', [1:4,7:10],'xticklabel', [1:4,1:4]);
    xlim([0 11]); ylim([0 .5]);
    set(gca,'fontsize',14);
    ylabel('Proportion Responded');
    xlabel('Confidence Rating');
    title('Control Group');
    t(1) = text(1.75,.3,'Pre');
    t(2) = text(7.75,.3,'Post');
    set(t, 'fontsize', 14);
    legend(hBar(2:-1:1), 'Correct','Incorrect','location','nw');
    legend boxoff;
    box off;
    subplot(1,2,2); % Experimental Group
    hBar = bar(1:4,[confDistr.group_2.session_01.incorrect.mean/sum(confDistr.group_2.session_01.incorrect.mean + confDistr.group_2.session_01.correct.mean); confDistr.group_2.session_01.correct.mean/sum(confDistr.group_2.session_01.incorrect.mean + confDistr.group_2.session_01.correct.mean)]'); hold on;
    hBar2 = bar(7:10,[confDistr.group_2.session_10.incorrect.mean/sum(confDistr.group_2.session_10.incorrect.mean + confDistr.group_2.session_10.correct.mean); confDistr.group_2.session_10.correct.mean/sum(confDistr.group_2.session_10.incorrect.mean + confDistr.group_2.session_10.correct.mean)]');
    set(hBar(1),'facecolor', 'r');
    set(hBar(2),'facecolor', 'g');
    set(hBar2(1),'facecolor', 'r');
    set(hBar2(2),'facecolor', 'g');
    set(gca, 'xtick', [1:4,7:10],'xticklabel', [1:4,1:4]);
    xlim([0 11]); ylim([0 .5]);
    set(gca,'fontsize',14);
    ylabel('Proportion Responded');
    xlabel('Confidence Rating');
    title('Experimental Group');
    t(1) = text(1.75,.3,'Pre');
    t(2) = text(7.75,.3,'Post');
    set(t, 'fontsize', 14);
    box off;
    if exportFigs
        export_fig confDistrPlot -png -transparent 'confDistr.png';
    end    
    
    confDistrxDomainPlot = figure;
    set(gcf,'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    hBar = bar(1:4,plotConfDistr.group_1.session_01.memory); hold on;
    hBar2 = bar(6:9, plotConfDistr.group_1.session_01.perception);
    hBar3 = bar(12:15,plotConfDistr.group_1.session_10.memory);
    hBar4 = bar(17:20,plotConfDistr.group_1.session_10.perception);
    set(hBar(1),'facecolor', 'r');
    set(hBar(2),'facecolor', 'g');
    set(hBar2(1),'facecolor', 'r');
    set(hBar2(2),'facecolor', 'g');
    set(hBar3(1),'facecolor', 'r');
    set(hBar3(2),'facecolor', 'g');
    set(hBar4(1),'facecolor', 'r');
    set(hBar4(2),'facecolor', 'g');
    set(gca, 'xtick', [1:4,6:9,12:15,17:20],'xticklabel', [1:4,1:4,1:4,1:4]);
    xlim([0 21]); ylim([0 .6]);
    set(gca,'fontsize',14);
    ylabel('Proportion Responded');
    xlabel('Confidence Rating');
    title('Control Group');
    t(1) = text(1,.3,'Mem');
    t(2) = text(6,.3,'Per');
    t(3) = text(12,.3,'Mem');
    t(4) = text(17,.3,'Per');
    set(t, 'fontsize', 14);
    t2(1) = text(4, .375, 'Pre');
    t2(2) = text(14, .375, 'Post');
    set(t2, 'fontsize', 18);
    legend(hBar(2:-1:1), 'Correct','Incorrect','location','nw');
    legend boxoff;
    box off;
    subplot(1,2,2); % Experimental Group
    hBar = bar(1:4,plotConfDistr.group_2.session_01.memory); hold on;
    hBar2 = bar(6:9, plotConfDistr.group_2.session_01.perception);
    hBar3 = bar(12:15,plotConfDistr.group_2.session_10.memory);
    hBar4 = bar(17:20,plotConfDistr.group_2.session_10.perception);
    set(hBar(1),'facecolor', 'r');
    set(hBar(2),'facecolor', 'g');
    set(hBar2(1),'facecolor', 'r');
    set(hBar2(2),'facecolor', 'g');
    set(hBar3(1),'facecolor', 'r');
    set(hBar3(2),'facecolor', 'g');
    set(hBar4(1),'facecolor', 'r');
    set(hBar4(2),'facecolor', 'g');
    set(gca, 'xtick', [1:4,6:9,12:15,17:20],'xticklabel', [1:4,1:4,1:4,1:4]);
    xlim([0 21]); ylim([0 .6]);
    set(gca,'fontsize',14);
    ylabel('Proportion Responded');
    xlabel('Confidence Rating');
    title('Experimental Group');
    t(1) = text(1,.3,'Mem');
    t(2) = text(6,.3,'Per');
    t(3) = text(12,.3,'Mem');
    t(4) = text(17,.3,'Per');
    set(t, 'fontsize', 14);
    t2(1) = text(4, .375, 'Pre');
    t2(2) = text(14, .375, 'Post');
    set(t2, 'fontsize', 18);
    legend(hBar(2:-1:1), 'Correct','Incorrect','location','nw');
    legend boxoff;
    box off;
    if exportFigs
        export_fig confDistrxDomainPlot -png -transparent 'confDistrxDomain.png';
    end    
end
end