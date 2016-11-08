function [ analysis, results ] = reformatResults( analysis, results )
%REFORMATRESULTS Converts stimulus from abstract/words to trained/untrained
%for each subject individually

dom = {'perception', 'memory'};

subjects = fieldnames(results);
for sub = 1:numel(subjects)
    if strcmp(char(results.(subjects{sub}).stimCond),'abstract')
        sessions = fieldnames(results.(subjects{sub}));
        for sesh = 1:numel(sessions)
            if strncmp(sessions{sesh},'session',7)
                session = str2double(sessions{sesh}(end-1:end));
                if session == 1 ||  session == 10
                    for d = 1:numel(dom)
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).trained = analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).abstract;
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).untrained = analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).words;
                        results.(subjects{sub}).(sessions{sesh}).(dom{d}).trained = results.(subjects{sub}).(sessions{sesh}).(dom{d}).abstract;
                        results.(subjects{sub}).(sessions{sesh}).(dom{d}).untrained = results.(subjects{sub}).(sessions{sesh}).(dom{d}).words;
                    end
                else % sessions 2-9
                    analysis.(subjects{sub}).(sessions{sesh}).perception.trained = analysis.(subjects{sub}).(sessions{sesh}).perception.abstract;
                    results.(subjects{sub}).(sessions{sesh}).perception.trained = results.(subjects{sub}).(sessions{sesh}).perception.abstract;
                end
            end
        end
    elseif strcmp(char(results.(subjects{sub}).stimCond),'words')
        sessions = fieldnames(results.(subjects{sub}));
        for sesh = 1:numel(sessions)
            if strncmp(sessions{sesh},'session',7)
                session = str2double(sessions{sesh}(end-1:end));
                if session == 1 ||  session == 10
                    for d = 1:numel(dom)
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).trained = analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).words;
                        analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).untrained = analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).abstract;
                        results.(subjects{sub}).(sessions{sesh}).(dom{d}).trained = results.(subjects{sub}).(sessions{sesh}).(dom{d}).words;
                        results.(subjects{sub}).(sessions{sesh}).(dom{d}).untrained = results.(subjects{sub}).(sessions{sesh}).(dom{d}).abstract;
                    end
                else % sessions 2-9
                    analysis.(subjects{sub}).(sessions{sesh}).perception.trained = analysis.(subjects{sub}).(sessions{sesh}).perception.words;
                    results.(subjects{sub}).(sessions{sesh}).perception.trained = results.(subjects{sub}).(sessions{sesh}).perception.words;
                end
            end
        end
    else
        error('No stimCond defined for subject %s', subjects{sub});
    end
end


end

