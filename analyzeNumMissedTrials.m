function [ pctMissedTrials ] = analyzeNumMissedTrials( analysis )
%ANALYZENUMMISSEDTRIALS Reports percent of trials that were not analyzed (FFX)
dom = {'perception', 'memory'};
stim = {'trained', 'untrained'};
subjects = fieldnames(analysis);

numTrials = 0;
numMissedTrials = 0;

for sub = 1:numel(subjects)
    if strncmp(subjects{sub},'subject',7)
        sessions = fieldnames(analysis.(subjects{sub}));
        for sesh = 1:10
            session = sprintf('session_%.2d', sesh);
            if sesh == 1 || sesh == 10
                for d = 1:numel(dom)
                    for s = 1:numel(stim)
                        numTrials = numTrials + analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).numTrials;
                        numMissedTrials = numMissedTrials + analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).numMissedTrials;
                    end
                end
            else % training sessions 2-9
                numTrials = numTrials + analysis.(subjects{sub}).(session).perception.trained.numTrials;
                numMissedTrials = numMissedTrials + analysis.(subjects{sub}).(session).perception.trained.numMissedTrials;
            end
        end
    end
end

pctMissedTrials = numMissedTrials / numTrials;
