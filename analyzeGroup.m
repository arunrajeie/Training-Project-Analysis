function [ groupResults ] = analyzeGroup( data, analysis, results, analysesToRun, plotFigs, exportFigs )
%GROUPANALYSIS Allows the user to select which group analyses to run,
%whether to plot the outputs, and whether to export those figures
%   Input 1: data structure from readData function
%   Input 2: analysis structure from analyzeData function
%   Input 3: results structure from analyzeData function
%   Input 4: cell string of analyses to run
%       e.g. analysesToRun = {'meanConf', 'M_ratio', 'rt'};
%       e.g. analysesToRun = {'da', 'difficulty'};
%   Optional input 2: boolean determining whether or not to plot the output
%   of the analyses (default is false)
%   Optional input 3: boolean determining whether or not to export the
%   figures (default is false)
%
%   Output: groupResults structure including all selected analyses

if nargin < 5
    exportFigs = false;
    if nargin < 4
        plotFigs = false;
    end
end

groupResults = struct();

if strcmp(analysesToRun, 'all')
    groupResults.pC = groupPercentCorrect(results, plotFigs, exportFigs);
    groupResults.da = groupDa(results, plotFigs, exportFigs);
    groupResults.metaDa = groupMetaDa(results, plotFigs, exportFigs);
    groupResults.M_ratio = groupM_ratio(results, plotFigs, exportFigs);
    groupResults.logM_ratio = groupLogM_ratio(results, plotFigs, exportFigs);
    groupResults.meanConf = groupMeanConf(analysis, results, plotFigs, exportFigs);
    groupResults.criteria = groupCriteria(results, plotFigs, exportFigs);
    groupResults.AUC = groupAUC(results, plotFigs, exportFigs);
    [groupResults.t2HR, groupResults.t2FAR] = groupT2HRFAR(results, plotFigs, exportFigs);
    groupResults.confDistr = groupConfDistr(analysis, plotFigs, exportFigs);
    groupResults.rt = groupRT(results, plotFigs, exportFigs);
    groupResults.confRT = groupConfRT(results, plotFigs, exportFigs);
    groupResults.difficulty = groupDifficulty(results, plotFigs, exportFigs);
    groupResults.QSR = groupQSR(analysis, results, plotFigs, exportFigs);
    groupResults.points = groupPoints(data, plotFigs, exportFigs);
else
    if any(strcmp(analysesToRun,'pC'))
        groupResults.pC = groupPercentCorrect(results, plotFigs, exportFigs);
    end
    if any(strcmp(analysesToRun,'da'))
        groupResults.da = groupDa(results, plotFigs, exportFigs);
    end
    if any(strcmp(analysesToRun,'metaDa'))
        groupResults.metaDa = groupMetaDa(results, plotFigs, exportFigs);
    end
    if any(strcmp(analysesToRun,'M_ratio'))
        groupResults.M_ratio = groupM_ratio(results, plotFigs, exportFigs);
    end
    if any(strcmp(analysesToRun,'logM_ratio'))
        groupResults.logM_ratio = groupLogM_ratio(results, plotFigs, exportFigs);
    end
    if any(strcmp(analysesToRun,'meanConf'))
        groupResults.meanConf = groupMeanConf(analysis, results, plotFigs, exportFigs);
    end
    if any(strcmp(analysesToRun, 'criteria'))
        groupResults.criteria = groupCriteria(results, plotFigs, exportFigs);
    end
    if any(strcmp(analysesToRun,'AUC'))
        groupResults.AUC = groupAUC(results, plotFigs, exportFigs);
    end
    if any(strcmp(analysesToRun,'t2HRFAR'))
        [groupResults.t2HR, groupResults.t2FAR] = groupT2HRFAR(results, plotFigs, exportFigs);
    end
    if any(strcmp(analysesToRun,'confDistr'))
       groupResults.confDistr = groupConfDistr(analysis, plotFigs, exportFigs); 
    end
    if any(strcmp(analysesToRun,'rt'))
        groupResults.rt = groupRT(results, plotFigs, exportFigs);
    end
    if any(strcmp(analysesToRun, 'confRT'))
        groupResults.confRT = groupConfRT(results, plotFigs, exportFigs);
    end
    if any(strcmp(analysesToRun, 'difficulty'))
        groupResults.difficulty = groupDifficulty(results, plotFigs, exportFigs);
    end
    if any(strcmp(analysesToRun, 'QSR'))
        groupResults.QSR = groupQSR(analysis, results, plotFigs, exportFigs);
    end
    if any(strcmp(analysesToRun, 'points'))
        groupResults.points = groupPoints(data, plotFigs, exportFigs);
    end
end
    
end

