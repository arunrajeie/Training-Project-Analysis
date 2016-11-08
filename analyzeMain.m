load('results_62subs.mat')

% set paths for code we need
rootdir = fileparts(which('analyzeMain')); % Directory in which we put analyzeMain.m
addpath(rootdir); % Add rootdir path
addpath(fullfile(rootdir,'altmany-export_fig-4c015d5')); % Add export_fig path

% All complete, non-excluded subjects (62)
subjects = [1, 2, 4, 5, 8, 12, 13, 14, 15, 16, 17, 18, 19, 20, 23,...
    24, 25, 26, 27, 33, 34, 35, 40, 41, 42, 46, 47, 48, 51, 52, 53,...
    56, 57, 58, 60, 61, 62, 66, 67, 69, 70, 71, 72, 73, 76, 77, 80,...
    83, 84, 85, 86, 88, 89, 90, 91, 93, 94, 96, 98, 99, 100, 102];
% Exclude subjects 6, 9, 11, 29, 30, 74, 78 based on previous analysis
sessions = 1:10;

fitMetaD = true;
analysesToRun = {'all'}; % pC, da, metaDa, M_ratio, logM_ratio, meanConf, criteria, AUC, t2HRFAR, confDistr, rt, confRT, difficulty, QSR, GLM (or regression)|| all
plotFigs = true;
exportFigs = false;
saveResults = false;

%% Code is commented out because we cannot provide access to raw data on GitHub. Skip ahead to group-level and supplementary analyses with optional plotting
% data = readData(subjects, sessions);
% [analysis, results] = analyzeSubject(data, fitMetaD);
% [analysis, results] = reformatResults(analysis, results); % Reformat in terms of trained/untrained stimulus
groupResults = analyzeGroup(data, analysis, results, analysesToRun, plotFigs, exportFigs);
analyzeStrategyPT(results, plotFigs, exportFigs); 
[results, supplementary] = analyzeSupplementary(data, results);
analyzeAgeCorrelation(results, plotFigs, exportFigs);
analyzeDuration(results, plotFigs, exportFigs);

if saveResults
    fName = 'results_62subs.mat';
    save(fName, 'data', 'analysis', 'results', 'groupResults', 'supplementary', 'subjects', 'sessions');
    fprintf('DONE -- saved to %s', fName);
end
