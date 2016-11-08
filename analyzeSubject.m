function [ analysis, results ] = analyzeSubject( data, fitMetaD )
%ANALYZESUBJECT Creates analysis and summary statistics for individual subject data
%   Inputs: The data structure from the readData file and (optional) boolean
%   determining whether or not to fit meta-d (default = false)
%   Output: Intermediate analysis, and summary statistics results for all 
%   subject(s) and session(s) in the data structure

if nargin == 1
    fitMetaD = false;
end

dataFields = fieldnames(data);
if any(strncmp(dataFields,'subject',7)) % Read in multiple subjects -> analyze all
    subjects = fieldnames(data);
    for sub = 1:numel(subjects)
        sessions = fieldnames(data.(subjects{sub}));        
        for sesh = 1:numel(sessions)
            analysis.(subjects{sub}).(sessions{sesh}) = struct();
            results.(subjects{sub}).(sessions{sesh}) = struct();
            
            dom = {'perception', 'memory'};
            type = {'response', 'trial'};
            stim = {'abstract', 'words'};
            
            session = str2double(sessions{sesh}(end-1:end));
            if session == 1 || session == 10
                numTrials = 108;              
            else % Sessions 2-9
                numTrials = 270;
                analysis.(subjects{sub}).group = data.(subjects{sub}).(sessions{sesh}).feedbackCond(1);
                analysis.(subjects{sub}).stimCond = data.(subjects{sub}).(sessions{sesh}).stimCond(1);
                results.(subjects{sub}).group = data.(subjects{sub}).(sessions{sesh}).feedbackCond(1);
                results.(subjects{sub}).stimCond = data.(subjects{sub}).(sessions{sesh}).stimCond(1);
            end
            
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d}) = find(strcmp(data.(subjects{sub}).(sessions{sesh}).domain,dom{d}) & strcmp(data.(subjects{sub}).(sessions{sesh}).stimType, stim{s}) & strcmp(data.(subjects{sub}).(sessions{sesh}).trialType, type{d}));
                    index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating = find(strcmp(data.(subjects{sub}).(sessions{sesh}).domain,dom{d}) & strcmp(data.(subjects{sub}).(sessions{sesh}).stimType, stim{s}) & strcmp(data.(subjects{sub}).(sessions{sesh}).trialType, 'confidence_rating'));
                    
                    if size(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) ~= 0
                        %% Analysis
                        % Initialize the variables to empty/NaNs
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).stimID = nan(numTrials,1); % Where the stimulus was presented on each trial (0 for left, 1 for right)
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).resp = nan(numTrials,1); % What the subject's response was for each trial (0 for left, 1 for right)
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).acc = nan(numTrials,1); % Array tracking trial-by-trial accuracy (0 for incorrect, 1 for correct)
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).rt = nan(numTrials,1); % Array tracking trial-by-trial rt
                        if strcmp(dom{d}, 'perception')
                            analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).brightness = nan(numTrials,1); % Array tracking trial-by-trial brightness
                            analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).difficulty = nan(numTrials,1); % Array tracking trial-by-trial difficulty
                        end
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).conf = nan(numTrials,1); % Array tracking trial-by-trial confidence high/lo (0 for low, 1 for high)
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confResp = nan(numTrials,1); % Array tracking trial-by-trial confidence responses (values = 1-4)
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).QSR = nan(numTrials,1); % Array tracking trial-by-trial QSR score
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRT = nan(numTrials,1); % Array tracking trial-by-trial confidence rt
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).correct = 0; % Count tracking how many total correct trials the subject got in that condition
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).accBinned = zeros(1,4); % Array tracking total accuracy by confidence rating (each bin contains the number of trials they got correct by pressing the corresponding button)
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).numTrials = 0; % Count tracking total number of trials the subject got in that condition
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).numTrialsBinned = zeros(1,4); % Array tracking total number of trials by confidence rating (each bin contains the number of trials they responded to by pressing the corresponding button)
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).numMissedTrials = 0; % Number of NaNs in resp array (corresponds to failure to respond on trial or missing data)

                        % Convert data to format we can use
                        %% Response(per)/trial(mem) trials
                        % StimID
                        stimID = nan(size(data.(subjects{sub}).(sessions{sesh}).trialNum(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})),1),1); 
                        stimID((data.(subjects{sub}).(sessions{sesh}).rt(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) > 200) & ... % If they responded in time, but not too fast
                               (data.(subjects{sub}).(sessions{sesh}).correct(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) == 1) & ... % And they were correct
                               (data.(subjects{sub}).(sessions{sesh}).key_press(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) == 79)) = 0; % And they pressed 'O', set stimID = 0
                        stimID((data.(subjects{sub}).(sessions{sesh}).rt(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) > 200) & ... % If they responded in time, but not too fast
                               (data.(subjects{sub}).(sessions{sesh}).correct(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) == 1) & ... % And they were correct
                               (data.(subjects{sub}).(sessions{sesh}).key_press(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) == 80)) = 1; % And they pressed 'P', set stimID = 1
                        stimID((data.(subjects{sub}).(sessions{sesh}).rt(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) > 200) & ... % If they responded in time, but not too fast
                               (data.(subjects{sub}).(sessions{sesh}).correct(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) == 0) & ... % And they were incorrect
                               (data.(subjects{sub}).(sessions{sesh}).key_press(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) == 79)) = 1; % And they pressed 'O', set stimID = 1
                        stimID((data.(subjects{sub}).(sessions{sesh}).rt(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) > 200) & ... % If they responded in time, but not too fast
                               (data.(subjects{sub}).(sessions{sesh}).correct(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) == 0) & ... % And they were incorrect
                               (data.(subjects{sub}).(sessions{sesh}).key_press(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) == 80)) = 0; % And they pressed 'P', set stimID = 0
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).stimID(data.(subjects{sub}).(sessions{sesh}).trialNum(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d}))) = stimID;

                        % Response
                        resp = nan(size(data.(subjects{sub}).(sessions{sesh}).trialNum(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})),1),1);
                        resp((data.(subjects{sub}).(sessions{sesh}).rt(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) > 200) & ... % If they responded in time, but not too fast
                             (data.(subjects{sub}).(sessions{sesh}).key_press(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) == 79)) = 0; % And they pressed 'O', set resp = 0
                        resp((data.(subjects{sub}).(sessions{sesh}).rt(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) > 200) & ... % If they responded in time, but not too fast
                             (data.(subjects{sub}).(sessions{sesh}).key_press(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) == 80)) = 1; % And they pressed 'P', set resp = 1
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).resp(data.(subjects{sub}).(sessions{sesh}).trialNum(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d}))) = resp;

                        % Accuracy
                        acc = nan(size(data.(subjects{sub}).(sessions{sesh}).trialNum(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})),1),1);
                        acc((data.(subjects{sub}).(sessions{sesh}).rt(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) > 200) & ... % If they responded in time, but not too fast
                            (data.(subjects{sub}).(sessions{sesh}).correct(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) == 1)) = 1; % And they were correct, set acc = 1
                        acc((data.(subjects{sub}).(sessions{sesh}).rt(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) > 200) & ... % If they responded in time, but not too fast
                            (data.(subjects{sub}).(sessions{sesh}).correct(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})) == 0)) = 0; % And they were incorrect, set acc = 0
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).acc(data.(subjects{sub}).(sessions{sesh}).trialNum(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d}))) = acc; 

                        % Reaction Time
                        rt = data.(subjects{sub}).(sessions{sesh}).rt(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d})); % Get RT for that condition
                        rt(rt<=200) = nan;
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).rt(data.(subjects{sub}).(sessions{sesh}).trialNum(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d}))) = rt;

                        % Correct count
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).correct = nansum(acc); % Num Correct = Sum of accuracy
                        % Num Trials
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).numTrials = size(acc(~isnan(acc)),1); % Num trials = size of acc

                        %% Confidence rating trials
                        % Confidence (high=1/low=0)
                        conf = nan(size(data.(subjects{sub}).(sessions{sesh}).trialNum(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating),1),1);
                        conf((data.(subjects{sub}).(sessions{sesh}).rt(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating) > 0) & ... % If they responded in time
                             (data.(subjects{sub}).(sessions{sesh}).key_press(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating) == 49)) = 0; % And they pressed CR = 1, set conf = 0
                        conf((data.(subjects{sub}).(sessions{sesh}).rt(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating) > 0) & ... % If they responded in time
                             (data.(subjects{sub}).(sessions{sesh}).key_press(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating) == 50)) = 0; % And they pressed CR = 2, set conf = 0
                        conf((data.(subjects{sub}).(sessions{sesh}).rt(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating) > 0) & ... % If they responded in time
                             (data.(subjects{sub}).(sessions{sesh}).key_press(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating) == 51)) = 1; % And they pressed CR = 3, set conf = 1
                        conf((data.(subjects{sub}).(sessions{sesh}).rt(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating) > 0) & ... % If they responded in time
                             (data.(subjects{sub}).(sessions{sesh}).key_press(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating) == 52)) = 1; % And they pressed CR = 4, set conf = 1
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).conf(data.(subjects{sub}).(sessions{sesh}).trialNum(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating)) = conf;
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).conf(isnan(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).resp)) = nan;
                        
                        % confResp (1,2,3,4)
                        confResp = nan(size(data.(subjects{sub}).(sessions{sesh}).trialNum(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating),1),1);
                        confResp((data.(subjects{sub}).(sessions{sesh}).rt(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating) > 0) & ... % If they responded in time
                             (data.(subjects{sub}).(sessions{sesh}).key_press(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating) == 49)) = 1; % And they pressed CR = 1, set confResp = 1
                        confResp((data.(subjects{sub}).(sessions{sesh}).rt(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating) > 0) & ... % If they responded in time
                             (data.(subjects{sub}).(sessions{sesh}).key_press(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating) == 50)) = 2; % And they pressed CR = 2, set confResp = 2
                        confResp((data.(subjects{sub}).(sessions{sesh}).rt(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating) > 0) & ... % If they responded in time
                             (data.(subjects{sub}).(sessions{sesh}).key_press(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating) == 51)) = 3; % And they pressed CR = 3, set confResp = 3
                        confResp((data.(subjects{sub}).(sessions{sesh}).rt(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating) > 0) & ... % If they responded in time
                             (data.(subjects{sub}).(sessions{sesh}).key_press(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating) == 52)) = 4; % And they pressed CR = 4, set confResp = 4
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confResp(data.(subjects{sub}).(sessions{sesh}).trialNum(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating)) = confResp;
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confResp(isnan(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).resp)) = nan;
                        
                        % Quadratic Scoring Rule -- Assuming confidence 1-4 <--> 0-1 
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).QSR = 1-((analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).acc-(-1/3+analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confResp./3)).^2);
                        
                        for bin = 1:4
                            % Accuracy binned by confResp
                            analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).accBinned(bin) = nansum(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).acc(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confResp == bin));
                            % Num Trials binned by confResp
                            analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).numTrialsBinned(bin) = nansum(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confResp == bin);
                        end
                        
                        % Confidence Reaction Time
                        confRT = data.(subjects{sub}).(sessions{sesh}).rt(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating); % Get RT for that confidence rating
                        confRT(confRT<0) = nan;
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRT(data.(subjects{sub}).(sessions{sesh}).trialNum(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRating)) = confRT;
                        
                        % Trials failed to respond in time to or failed to respond at all or missing data
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).numMissedTrials = sum(isnan(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).resp));
                       
                        % Brightness and difficulty
                        if strcmp(dom{d}, 'perception')
                            if strcmp(stim{s}, 'abstract')
                                brightness = data.(subjects{sub}).(sessions{sesh}).avgBrightness_abstract(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d}));
                            elseif strcmp(stim{s}, 'words')
                                brightness = data.(subjects{sub}).(sessions{sesh}).avgBrightness_words(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d}));
                            end
                            analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).brightness(data.(subjects{sub}).(sessions{sesh}).trialNum(index.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).(type{d}))) = brightness;
                            analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).difficulty = 128-(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).brightness-128);
                        end
                        
                        %% Results
                        
                        % Accuracy
                        results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).percentCorrect = analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).correct / analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).numTrials;
                        results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).percentBinned = analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).accBinned ./ analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).numTrialsBinned;
                        
                        % Mean RTs
                        results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).meanRT = nanmean(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).rt);
                        results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).meanConfRT = nanmean(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confRT);
                        
                        % Mean Difficulty
                        if strcmp(dom{d}, 'perception')
                            results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).meanDifficulty = nanmean(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).difficulty);
                        end
                        
                        % Mean QSR
                        results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).meanQSR = nanmean(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).QSR);
                        
                        % Mean Conf
                        results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).meanConf = nanmean(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confResp);
                        
                        % Area Under the Curve (type2roc)
                        results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).AUC = type2roc(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).acc, analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confResp, 4);
                        
                        % Type 2 HR/FAR
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).hits = nansum(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).acc(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).conf == 1));
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).misses = nansum(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).acc(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).conf == 0));
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).FAs = nansum(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).acc(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).conf == 1) == 0);
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).CRs = nansum(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).acc(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).conf == 0) == 0);
                        results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).t2HR = analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).hits / (analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).hits +  analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).misses);
                        results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).t2FAR = analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).FAs / (analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).FAs +  analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).CRs);

                        
                        % Prepare the data for fit_meta_d_MLE by using the trials2counts function
                        [analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).nR_S1, analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).nR_S2] = trials2counts(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).stimID, analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).resp, analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confResp, 4, 1);
                        results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).nR_S1 = analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).nR_S1;
                        results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).nR_S2 = analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).nR_S2;
                        if fitMetaD
                            % Get the metaD analysis by using the fit_meta_d_MLE function
                            results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).fit = fit_meta_d_MLE(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).nR_S1, analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).nR_S2);
                        end
                    end
                end
            end
        end
        results.(subjects{sub}) = orderfields(results.(subjects{sub})); % Make output a bit prettier but ordering alphabetically
    end    
else % Only read in one subject -> only analyze that subject
    analysis = struct();
    results = struct();
    
    dom = {'perception', 'memory'};
    type = {'response', 'trial'};
    stim = {'abstract','words'};
    
    if data.trialNum(end) > 90 && data.trialNum(end) < 110 % Sessions 1 & 10
        numTrials = 108;
    else % Sessions 2-9
        numTrials = 270;
        analysis.group = data.feedbackCond(1);
        analysis.stimCond = data.stimCond(1);
        results.group = data.feedbackCond(1);
        results.stimCond = data.stimCond(1);
    end

    for d = 1:numel(dom)
        for s = 1:numel(stim)
            index.(dom{d}).(stim{s}).(type{d}) = find(strcmp(data.domain,dom{d}) & strcmp(data.stimType, stim{s}) & strcmp(data.trialType, type{d}));
            index.(dom{d}).(stim{s}).confRating = find(strcmp(data.domain,dom{d}) & strcmp(data.stimType, stim{s}) & strcmp(data.trialType, 'confidence_rating'));

            if size(index.(dom{d}).(stim{s}).(type{d})) ~= 0
                %% Analysis
                % Initialize the variables to empty/NaNs
                analysis.(dom{d}).(stim{s}).stimID = nan(numTrials,1); % Where the stimulus was presented on each trial (0 for left, 1 for right)
                analysis.(dom{d}).(stim{s}).resp = nan(numTrials,1); % What the subject's response was for each trial (0 for left, 1 for right)
                analysis.(dom{d}).(stim{s}).acc = nan(numTrials,1); % Array tracking trial-by-trial accuracy (0 for incorrect, 1 for correct)
                analysis.(dom{d}).(stim{s}).rt = nan(numTrials,1); % Array tracking trial-by-trial rt
                if strcmp(dom{d}, 'perception')
                    analysis.(dom{d}).(stim{s}).brightness = nan(numTrials,1); % Array tracking trial-by-trial brightness
                    analysis.(dom{d}).(stim{s}).difficulty = nan(numTrials,1); % Array tracking trial-by-trial difficulty
                end
                analysis.(dom{d}).(stim{s}).conf = nan(numTrials,1); % Array tracking trial-by-trial confidence high/lo (0 for low, 1 for high)
                analysis.(dom{d}).(stim{s}).confResp = nan(numTrials,1); % Array tracking trial-by-trial confidence responses (values = 1-4)
                analysis.(dom{d}).(stim{s}).QSR = nan(numTrials,1); % Array tracking trial-by-trial QSR score
                analysis.(dom{d}).(stim{s}).confRT = nan(numTrials,1); % Array tracking trial-by-trial confidence rt
                analysis.(dom{d}).(stim{s}).correct = 0; % Count tracking how many total correct trials the subject got in that condition
                analysis.(dom{d}).(stim{s}).accBinned = zeros(1,4); % Array tracking total accuracy by confidence rating (each bin contains the number of trials they got correct by pressing the corresponding button)
                analysis.(dom{d}).(stim{s}).numTrials = 0; % Count tracking total number of trials the subject got in that condition
                analysis.(dom{d}).(stim{s}).numTrialsBinned = zeros(1,4); % Array tracking total number of trials by confidence rating (each bin contains the number of trials they responded to by pressing the corresponding button)
                analysis.(dom{d}).(stim{s}).numMissedTrials = 0; % Number of NaNs in resp array (corresponds to failure to respond on trial or missing data)

                % Convert data to format we can use
                %% Response(per)/trial(mem) trials
                % StimID
                stimID = nan(size(data.trialNum(index.(dom{d}).(stim{s}).(type{d})),1),1); 
                stimID((data.rt(index.(dom{d}).(stim{s}).(type{d})) > 200) & ... % If they responded in time, but not too fast
                       (data.correct(index.(dom{d}).(stim{s}).(type{d})) == 1) & ... % And they were correct
                       (data.key_press(index.(dom{d}).(stim{s}).(type{d})) == 79)) = 0; % And they pressed 'O', set stimID = 0
                stimID((data.rt(index.(dom{d}).(stim{s}).(type{d})) > 200) & ... % If they responded in time, but not too fast
                       (data.correct(index.(dom{d}).(stim{s}).(type{d})) == 1) & ... % And they were correct
                       (data.key_press(index.(dom{d}).(stim{s}).(type{d})) == 80)) = 1; % And they pressed 'P', set stimID = 1
                stimID((data.rt(index.(dom{d}).(stim{s}).(type{d})) > 200) & ... % If they responded in time, but not too fast
                       (data.correct(index.(dom{d}).(stim{s}).(type{d})) == 0) & ... % And they were incorrect
                       (data.key_press(index.(dom{d}).(stim{s}).(type{d})) == 79)) = 1; % And they pressed 'O', set stimID = 1
                stimID((data.rt(index.(dom{d}).(stim{s}).(type{d})) > 200) & ... % If they responded in time, but not too fast
                       (data.correct(index.(dom{d}).(stim{s}).(type{d})) == 0) & ... % And they were incorrect
                       (data.key_press(index.(dom{d}).(stim{s}).(type{d})) == 80)) = 0; % And they pressed 'P', set stimID = 0
                analysis.(dom{d}).(stim{s}).stimID(data.trialNum(index.(dom{d}).(stim{s}).(type{d}))) = stimID;

                % Response
                resp = nan(size(data.trialNum(index.(dom{d}).(stim{s}).(type{d})),1),1);
                resp((data.rt(index.(dom{d}).(stim{s}).(type{d})) > 200) & ... % If they responded in time, but not too fast
                     (data.key_press(index.(dom{d}).(stim{s}).(type{d})) == 79)) = 0; % And they pressed 'O', set resp = 0
                resp((data.rt(index.(dom{d}).(stim{s}).(type{d})) > 200) & ... % If they responded in time, but not too fast
                     (data.key_press(index.(dom{d}).(stim{s}).(type{d})) == 80)) = 1; % And they pressed 'P', set resp = 1
                analysis.(dom{d}).(stim{s}).resp(data.trialNum(index.(dom{d}).(stim{s}).(type{d}))) = resp;

                % Accuracy
                acc = nan(size(data.trialNum(index.(dom{d}).(stim{s}).(type{d})),1),1);
                acc((data.rt(index.(dom{d}).(stim{s}).(type{d})) > 200) & ... % If they responded in time, but not too fast
                    (data.correct(index.(dom{d}).(stim{s}).(type{d})) == 1)) = 1; % And they were correct, set acc = 1
                acc((data.rt(index.(dom{d}).(stim{s}).(type{d})) > 200) & ... % If they responded in time, but not too fast
                    (data.correct(index.(dom{d}).(stim{s}).(type{d})) == 0)) = 0; % And they were incorrect, set acc = 0
                analysis.(dom{d}).(stim{s}).acc(data.trialNum(index.(dom{d}).(stim{s}).(type{d}))) = acc; 

                % Reaction Time
                rt = data.rt(index.(dom{d}).(stim{s}).(type{d})); % Get RT for that condition
                rt(rt<0) = nan;
                analysis.(dom{d}).(stim{s}).rt(data.trialNum(index.(dom{d}).(stim{s}).(type{d}))) = rt;

                % Correct count
                analysis.(dom{d}).(stim{s}).correct = nansum(acc); % Num Correct = Sum of accuracy
                % Num Trials
                analysis.(dom{d}).(stim{s}).numTrials = size(acc(~isnan(acc)),1); % Num trials = size of acc

                %% Confidence rating trials
                % Confidence (high=1/low=0)
                conf = nan(size(data.trialNum(index.(dom{d}).(stim{s}).confRating),1),1);
                conf((data.rt(index.(dom{d}).(stim{s}).confRating) > 0) & ... % If they responded in time
                     (data.key_press(index.(dom{d}).(stim{s}).confRating) == 49)) = 0; % And they pressed CR = 1, set conf = 0
                conf((data.rt(index.(dom{d}).(stim{s}).confRating) > 0) & ... % If they responded in time
                     (data.key_press(index.(dom{d}).(stim{s}).confRating) == 50)) = 0; % And they pressed CR = 2, set conf = 0
                conf((data.rt(index.(dom{d}).(stim{s}).confRating) > 0) & ... % If they responded in time
                     (data.key_press(index.(dom{d}).(stim{s}).confRating) == 51)) = 1; % And they pressed CR = 3, set conf = 1
                conf((data.rt(index.(dom{d}).(stim{s}).confRating) > 0) & ... % If they responded in time
                     (data.key_press(index.(dom{d}).(stim{s}).confRating) == 52)) = 1; % And they pressed CR = 4, set conf = 1
                analysis.(dom{d}).(stim{s}).conf(data.trialNum(index.(dom{d}).(stim{s}).confRating)) = conf;
                analysis.(dom{d}).(stim{s}).conf(isnan(analysis.(dom{d}).(stim{s}).resp)) = nan;

                % confResp (1,2,3,4)
                confResp = nan(size(data.trialNum(index.(dom{d}).(stim{s}).confRating),1),1);
                confResp((data.rt(index.(dom{d}).(stim{s}).confRating) > 0) & ... % If they responded in time
                     (data.key_press(index.(dom{d}).(stim{s}).confRating) == 49)) = 1; % And they pressed CR = 1, set confResp = 1
                confResp((data.rt(index.(dom{d}).(stim{s}).confRating) > 0) & ... % If they responded in time
                     (data.key_press(index.(dom{d}).(stim{s}).confRating) == 50)) = 2; % And they pressed CR = 2, set confResp = 2
                confResp((data.rt(index.(dom{d}).(stim{s}).confRating) > 0) & ... % If they responded in time
                     (data.key_press(index.(dom{d}).(stim{s}).confRating) == 51)) = 3; % And they pressed CR = 3, set confResp = 3
                confResp((data.rt(index.(dom{d}).(stim{s}).confRating) > 0) & ... % If they responded in time
                     (data.key_press(index.(dom{d}).(stim{s}).confRating) == 52)) = 4; % And they pressed CR = 4, set confResp = 4
                analysis.(dom{d}).(stim{s}).confResp(data.trialNum(index.(dom{d}).(stim{s}).confRating)) = confResp;
                analysis.(dom{d}).(stim{s}).confResp(isnan(analysis.(dom{d}).(stim{s}).resp)) = nan;

                % Quadratic Scoring Rule -- Assuming confidence 1-4 <--> 0-1 
                analysis.(dom{d}).(stim{s}).QSR = 1-((analysis.(dom{d}).(stim{s}).acc-(-1/3+analysis.(dom{d}).(stim{s}).confResp./3)).^2);

                for bin = 1:4
                    % Accuracy binned by confResp
                    analysis.(dom{d}).(stim{s}).accBinned(bin) = nansum(analysis.(dom{d}).(stim{s}).acc(analysis.(dom{d}).(stim{s}).confResp == bin));
                    % Num Trials binned by confResp
                    analysis.(dom{d}).(stim{s}).numTrialsBinned(bin) = nansum(analysis.(dom{d}).(stim{s}).confResp == bin);
                end

                % Confidence Reaction Time
                confRT = data.rt(index.(dom{d}).(stim{s}).confRating); % Get RT for that confidence rating
                confRT(confRT<0) = nan;
                analysis.(dom{d}).(stim{s}).confRT(data.trialNum(index.(dom{d}).(stim{s}).confRating)) = confRT;

                % Trials failed to respond in time to or failed to respond at all or missing data
                analysis.(dom{d}).(stim{s}).numMissedTrials = sum(isnan(analysis.(dom{d}).(stim{s}).resp));

                % Brightness and difficulty
                if strcmp(dom{d}, 'perception')
                    if strcmp(stim{s}, 'abstract')
                        brightness = data.avgBrightness_abstract(index.(dom{d}).(stim{s}).(type{d}));
                    elseif strcmp(stim{s}, 'words')
                        brightness = data.avgBrightness_words(index.(dom{d}).(stim{s}).(type{d}));
                    end
                    analysis.(dom{d}).(stim{s}).brightness(data.trialNum(index.(dom{d}).(stim{s}).(type{d}))) = brightness;
                    analysis.(dom{d}).(stim{s}).difficulty = 128-(analysis.(dom{d}).(stim{s}).brightness-128);
                end

                %% Results

                % Accuracy
                results.(dom{d}).(stim{s}).percentCorrect = analysis.(dom{d}).(stim{s}).correct / analysis.(dom{d}).(stim{s}).numTrials;
                results.(dom{d}).(stim{s}).percentBinned = analysis.(dom{d}).(stim{s}).accBinned ./ analysis.(dom{d}).(stim{s}).numTrialsBinned;

                % Mean RTs
                results.(dom{d}).(stim{s}).meanRT = nanmean(analysis.(dom{d}).(stim{s}).rt);
                results.(dom{d}).(stim{s}).meanConfRT = nanmean(analysis.(dom{d}).(stim{s}).confRT);

                % Mean Difficulty
                if strcmp(dom{d}, 'perception')
                    results.(dom{d}).(stim{s}).meanDifficulty = nanmean(analysis.(dom{d}).(stim{s}).difficulty);
                end

                % Mean QSR
                results.(dom{d}).(stim{s}).meanQSR = nanmean(analysis.(dom{d}).(stim{s}).QSR);

                % Mean Conf
                results.(dom{d}).(stim{s}).meanConf = nanmean(analysis.(dom{d}).(stim{s}).confResp);

                % Area Under the Curve (type2roc)
                results.(dom{d}).(stim{s}).AUC = type2roc(analysis.(dom{d}).(stim{s}).acc, analysis.(dom{d}).(stim{s}).confResp, 4);

                % Prepare the data for fit_meta_d_MLE by using the trials2counts function
                [analysis.(dom{d}).(stim{s}).nR_S1, analysis.(dom{d}).(stim{s}).nR_S2] = trials2counts(analysis.(dom{d}).(stim{s}).stimID, analysis.(dom{d}).(stim{s}).resp, analysis.(dom{d}).(stim{s}).confResp, 4, 1);
                results.(dom{d}).(stim{s}).nR_S1 = analysis.(dom{d}).(stim{s}).nR_S1;
                results.(dom{d}).(stim{s}).nR_S2 = analysis.(dom{d}).(stim{s}).nR_S2;
                if fitMetaD
                    % Get the metaD analysis by using the fit_meta_d_MLE function
                    results.(dom{d}).(stim{s}).fit = fit_meta_d_MLE(analysis.(dom{d}).(stim{s}).nR_S1, analysis.(dom{d}).(stim{s}).nR_S2);
                end
            end
        end
    end
end

end

