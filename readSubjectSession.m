function [ data ] = readSubjectSession( dataPath )
%READSUBJECTSESSION Reads in data from CSV for one subject/session
    if isequal(exist(dataPath,'file'),2) % 2 means it's a file.
        display(['Found file: ' dataPath]);
        [num, txt, raw] = xlsread(dataPath);
        dataMatrix = raw;
    elseif isequal(exist(filename, 'dir'),7) % 7 = directory.
        display('Found the folder');
    else
        display(['Error! ' fileName ' not a file or folder!']);
    end
    
    if (size(dataMatrix) ~= 0)
        data.rt                     = cell2mat(dataMatrix(:,1)); % Reaction time
        data.key_press              = cell2mat(dataMatrix(:,2)); % Key press
        data.trial_type             = dataMatrix(:,3); % Trial Type (as defined by jsPsych)
        data.workerID               = dataMatrix(:,4); % WorkerID
        data.domain                 = dataMatrix(:,5); % Domain (Perception/Memory)
        data.stimType               = dataMatrix(:,6); % Stimulus Type (Words/Abstract)
        data.trialType              = dataMatrix(:,7); % Trial Type (as defined in code by Jason Carpenter)
        data.correct                = cell2mat(dataMatrix(:,8)); % Correctness of the trial 
        data.trialNum               = cell2mat(dataMatrix(:,9)); % Trial number of the current trial (each trial starts with type 1 response and includes type 2 response)
        data.avgBrightness_abstract = cell2mat(dataMatrix(:,10)); % Average perceptual brightness in the abstract stimulus type
        data.avgBrightness_words    = cell2mat(dataMatrix(:,11)); % Average perceptual brightness in the words stimulus type
        data.mem_blocksize_abstract = cell2mat(dataMatrix(:,12)); % Average memory blocksize in the abstract stimulus type
        data.mem_blocksize_words    = cell2mat(dataMatrix(:,13)); % Average memory blocksize in the words stimulus type
        data.time_elapsed           = cell2mat(dataMatrix(:,14)); % Time elapsed since the beginning of the experiment
        data.feedbackCond           = cell2mat(dataMatrix(:,15)); % Feedback condition (Type 1/Type 2)
        data.stimCond               = dataMatrix(:,16); % Stimulus training condition (Abstract/Words)
        data.pointsEarned           = cell2mat(dataMatrix(:,17)); % Points earned
    end
end

