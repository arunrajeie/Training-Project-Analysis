function [ pC ] = analyzePercentCorrect( results )
%GROUPDA Runs group d' analysis and optionally plots and exports the figures

dom = {'perception', 'memory'};
stim = {'abstract', 'words'};
for d = 1:numel(dom)
    for s = 1:numel(stim)
        pC.(dom{d}).(stim{s}).raw = [];
    end
end

subjects = fieldnames(results);
% Concatenate raw data
for sub = 1:numel(subjects)
    sessions = fieldnames(results.(subjects{sub}));
    for sesh = 1:numel(sessions)
        if strncmp(sessions{sesh},'session',7)
            session = str2double(sessions{sesh}(end-1:end));
            if session == 1 || session == 10
                for d = 1:numel(dom)
                    for s = 1:numel(stim)
                        pC.(dom{d}).(stim{s}).raw = vertcat(pC.(dom{d}).(stim{s}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).percentCorrect);
                    end
                end
            else % Sessions 2-9
                if strcmp(results.(subjects{sub}).stimCond, 'abstract')
                    pC.perception.abstract.raw = vertcat(pC.perception.abstract.raw, results.(subjects{sub}).(sessions{sesh}).perception.abstract.percentCorrect);
                elseif strcmp(results.(subjects{sub}).stimCond, 'words')
                    pC.perception.words.raw = vertcat(pC.perception.words.raw, results.(subjects{sub}).(sessions{sesh}).perception.words.percentCorrect);
                end
            end
        end
    end    
end

% Take mean and standard error
for d = 1:numel(dom)
    for s = 1:numel(stim)
        pC.(dom{d}).(stim{s}).mean = nanmean(pC.(dom{d}).(stim{s}).raw);
        pC.(dom{d}).(stim{s}).sem = nanstd(pC.(dom{d}).(stim{s}).raw)/sqrt(length(pC.(dom{d}).(stim{s}).raw));
    end
end

end