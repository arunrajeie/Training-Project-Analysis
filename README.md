# Training-Project-Analysis
Analysis scripts for a psychophysics experiment which demonstrates that adaptive training of sensory metacognition leads to confidence gain and enhanced insight. 

The data included in this file is anonymized for privacy purposes.

------------------ Adaptive Training of Sensory Metacognition -------------------
---------------------------- Analysis Scripts README ----------------------------

analyzeMain.m -- Main script you need to run analysis. Makes calls to primary analysis functions: readData, analyzeSubject, reformatResults, and analyzeGroup. Also calls some secondary analysis functions: analyzeStrategy, analyzeSupplementary, analyzeAgeCorrelation, and analyzeDuration. NOTE: requires that subfolder 'data' contains all the .csv files with data for each subject/session selected
  Variables that need to be set in analyzeMain.m:
		subjects        -- the subjects you wish to analyze (be sure not to include subjects that were excluded by previous analysis)
		sessions        -- the sessions you wish to analyze (analyzeGroup and the secondary analysis functions all require that this is all of the sessions = 1:10)
		fitMetaD        -- true/false if you wish to use fit meta-d' to the behavioral data
		analysesToRun 	-- N x 1 cell of strings listing all the analyses you wish to run at the group level. If you wish to run all of them, for simplicity just set analysesToRun = {'all'}. Can also input any combination of the following: {'pC, 'da', 'metaDa', 'M_ratio', 'logM_ratio', 'meanConf', 'criteria', 'AUC', 't2HRFAR', 'confDistr', 'rt', 'confRT', 'difficulty', 'QSR', 'points'}
		plotFigs        -- true/false if you wish to show figures displaying the analyses that were run.
		exportFigs      -- true/false if you wish to use export_fig to export the figures to the folder. NOTE: requires that you set plotFigs = true AND that you have the subfolder 'altmany-export_fig-4c015d5' which contains all the scripts necessary to run export_fig
		saveResults 	  -- true/false if you wish to save the output variables from analyzeMain.m to file 'results_62subs.mat'
		
	Important output variables that will be in your workspace after running analyzeMain.m:
		data            -- Structure containing raw data for each subject/session defined in analyzeMain.m
		analysis        -- Structure containing intermediate analysis values for each subject/session, split by each domain/stimulus type condition
		results         -- Structure containing final summary results values for each subject/session, split by each domain/stimulus type condition
		groupResults 	  -- Structure containing group level analyses, with a substructure for each analysis defined in analysesToRun. Primarily contains group-level mean/sem per session, split by each domain/stimulus type condition.
		supplementary 	-- Structure containing supplementary group level analyses, with a substructure for each analysis. The numeric supplementary values contain t-tests for a difference between groups.
		
Primary analysis functions:
	readData.m 		     -- Function for reading in data for defined subject/session combination. Makes a call to readSubjectSession.m, which converts the .csv files into substructures of the data variable. If multiple subjects and/or sessions are defined, creates a substructure for each subject/session combination. If only one subject/session is defined, output variable data directly references the raw data.
	analyzeSubject.m 	  -- Function for producing intermediate analysis and summary statistics (results) for individual subject data. Produces analysis and results variables for all subjects/sessions in data structure. Can handle data either with substrucutres for each subject/session or data as defined for only a single subject/session from readData.m
	reformatResults.m 	-- Preparatory function for group-level analyses. Reformats analysis and results variables with respect to the trained and untrained stimulus during the training sessions 2-9. E.g. if the trained stimulus was words, then perception.trained = perception.words, perception.untrained = perception.abstract. 
	analyzeGroup.m 		  -- Messenger function which makes calls to each group-level analysis specified in analysesToRun and saves their outputs as substructures to groupResults variable. Optionally plots and exports figures.
		Individual group-level functions:
			groupPercentCorrect.m 	-- Percent correct analysis
			groupDa.m           	  -- d' analysis
      groupMetaDa.m       	  -- meta-d' analysis
			groupM_ratio.m 		      -- meta-d'/d' analysis
			grouplogM_ratio.m 	    -- log(meta-d'/d') analysis
			groupMeanConf.m 	      -- Mean confidence analysis
			groupCriteria.m 	      -- Criteria shift analysis
			groupAUC.m            	-- Area under the Curve for Type 2 ROC analysis
			groupT2HRFAR.m 		      -- Type 2 hit rate and false alarm rate analysis
			groupConfDistr.m 	      -- Confidence Distribution analysis
			groupRT.m           	  -- (Type 1) Reaction time analysis
			groupConfRT.m 		      -- (Type 2) Confidence reaction time analysis
			groupDifficulty.m 	    -- Difficulty level analysis
       groupQSR.m          	  -- QSR score analysis  
      groupPoints.m       	  -- Points feedback analysis
	
Secondary analysis functions:
	analyzeStrategyPT.m         -- Function for analyzing dissociation between mean confidence and (bias) and meta-d'/d' (efficiency)
	analyzeSupplementary.m 		  -- Function for analyzing supplementary factors such as age, gender, bonus earned, stimulus trained condition, work time (minutes), duration of training (days)
	analyzeAgeCorrelation.m 	  -- Function for analyzing correlation between age and log(meta-d'/d')
	analyzeDuration.m           -- Function for analyzing correlations between duration of training and measures of training such as log(meta-d'/d') difference, AUC difference, and mean confidence difference.
    analyzePercentCorrect.m   -- Function for analyzing percent correct across groups, sessions, and tasks
    analyzeNumMissedTrials.m  -- Function for reporting how many trials were not analyzed in the main analysis

NOTE: Once you have run analyzeMain.m once, in order to see the output for any particular group-level or secondary analysis, simply call that function with the appropriate parameters (usually analysis and/or results, plotFigs, and exportFigs). No need to run analyzeMain.m all over again.
