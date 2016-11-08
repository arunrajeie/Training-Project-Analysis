function [ results, supp ] = analyzeSupplementary( data, results )
%ANALYZESUPPLEMENTARY Performs analysis on supplementary data and saves
%individual subject values

ages = [53, 36, 29, 43, 36, 34, 29, 36, 29, 29, 23, 44, 25, 40, 32, 50, 47, 50, 29, 44, 53, 31, 34, 26, 32, 27, 55, 58, 34, 47,...
    27, 32, 46, 58, 53, 39, 31, 20, 32, 25, 26, 33, 26, 48, 33, 49, 41, 57, 56, 64, 38, 49, 62, 29, 34, 32, 30 23, 40, 26, 31, 27];
genders = {'m', 'f', 'f', 'f', 'm', 'f', 'm', 'f', 'm', 'f', 'f', 'f', 'f', 'f', 'm', 'f', 'm', 'm', 'f', 'f', 'f', 'f', 'm', 'm', 'm',...
    'm', 'f', 'f', 'm', 'm', 'm', 'f', 'm', 'f', 'f', 'f', 'f', 'f', 'm', 'm', 'f', 'm', 'm', 'f', 'm', 'm', 'm', 'f', 'm', 'f', 'f', 'f', 'f',...
    'm', 'm', 'f', 'm', 'f', 'f', 'f', 'f', 'm'};
bonuses = [14.6, 15.3, 20.9, 19.5, 12.8, 19.5, 18.1, 19.5, 18.8, 18.1, 20.9, 18.1, 16.7, 21.6, 18.1, 17.4, 20.9, 16.7, 18.8, 18.8, 14.6, 16,...
    16, 20.9, 14.6, 16.7, 17.4, 18.1, 20.9, 16, 20.9, 16.7, 18.1, 20.8, 20.9, 20.9, 16, 19.5, 16, 18.1, 16, 15.3, 20.9, 16.7, 19.5, 16.7, ...
    18.1, 14.6, 15.3, 17.4, 18.8, 18.8, 20.9, 18.8, 18.8, 20.2, 20.9, 18.8, 16.7, 18.8, 18.8, 17.4];
% stimConds = {'w', 'a', 'a', 'a', 'w', 'w', 'a', 'a', 'w', 'w', 'w', 'a', 'a', 'w', 'w', 'w', 'a', 'a', 'w', 'w', 'a', 'w', 'a', 'w', 'a', 'w',...
%     'w', 'w', 'a', 'w', 'a', 'w', 'w', 'a', 'w', 'w', 'w', 'w', 'a', 'a', 'a', 'a', 'w', 'w', 'w', 'a', 'a', 'a', 'w', 'w', 'a', 'a', 'a', 'w',...
%     'a', 'w', 'a', 'w', 'w', 'a', 'a', 'a'};
worktimes = [];
durations = [11, 10, 10, 10, 16, 13, 12, 19, 16, 14, 24, 22, 20, 18, 24, 10, 10, 17, 10, 14, 31, 9, 19, 10, 19, 31, 12, 15, 10, 25,...
    19, 16, 11, 10, 18, 14, 17, 10, 27, 12, 14, 11, 19, 35, 12, 15, 13, 13, 10, 16, 11, 16, 14, 12, 12, 23, 10, 13, 12, 17, 10, 12]; 

groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    age.(groups{g}).raw = [];
    gender.(groups{g}).raw = {};
    bonus.(groups{g}).raw = [];
    stimCond.(groups{g}).raw = {};
    worktime.(groups{g}).raw = [];
    duration.(groups{g}).raw = [];
end

subjects = fieldnames(results);
for sub = 1:numel(subjects)
    group = sprintf('group_%d', results.(subjects{sub}).group);
    results.(subjects{sub}).age = ages(sub);
    results.(subjects{sub}).gender = genders{sub};
    results.(subjects{sub}).bonus = bonuses(sub);
    results.(subjects{sub}).worktime = 0;
    for sesh = 1:10
        session = sprintf('session_%.2d',sesh);
        results.(subjects{sub}).worktime = results.(subjects{sub}).worktime + data.(subjects{sub}).(session).time_elapsed(end)/60000;
    end
    results.(subjects{sub}).duration = durations(sub);
    age.(group).raw = vertcat(age.(group).raw, results.(subjects{sub}).age);
    gender.(group).raw{end+1} = results.(subjects{sub}).gender;
    bonus.(group).raw = vertcat(bonus.(group).raw, results.(subjects{sub}).bonus);
    stimCond.(group).raw{end+1} = results.(subjects{sub}).stimCond;
    worktimes = vertcat(worktimes, results.(subjects{sub}).worktime);
    worktime.(group).raw = vertcat(worktime.(group).raw, results.(subjects{sub}).worktime);
    duration.(group).raw = vertcat(duration.(group).raw, results.(subjects{sub}).duration);
end

for g = 1:numel(groups)
    age.(groups{g}).mean = mean(age.(groups{g}).raw);
    age.(groups{g}).sem = std(age.(groups{g}).raw)/sqrt(numel(age.(groups{g}).raw));
    
    bonus.(groups{g}).mean = mean(bonus.(groups{g}).raw);
    bonus.(groups{g}).sem = std(bonus.(groups{g}).raw)/sqrt(numel(bonus.(groups{g}).raw));
    
    stimCond.(groups{g}).abstract = 0;
    stimCond.(groups{g}).words = 0;
    for sub = 1:numel(stimCond.(groups{g}).raw)
        if strcmp(stimCond.(groups{g}).raw{sub},'abstract')
            stimCond.(groups{g}).abstract = stimCond.(groups{g}).abstract + 1;
        elseif strcmp(stimCond.(groups{g}).raw{sub},'words')
            stimCond.(groups{g}).words = stimCond.(groups{g}).words + 1;
        end
    end
    
    worktime.(groups{g}).mean = mean(worktime.(groups{g}).raw);
    worktime.(groups{g}).sem = std(worktime.(groups{g}).raw)/sqrt(numel(worktime.(groups{g}).raw));
    
    duration.(groups{g}).mean = mean(duration.(groups{g}).raw);
    duration.(groups{g}).sem = std(duration.(groups{g}).raw)/sqrt(numel(duration.(groups{g}).raw));
end

% T-tests
[age.h, age.p, age.ci, age.stats] = ttest2(age.group_1.raw, age.group_2.raw);

[bonus.h, bonus.p, bonus.ci, bonus.stats] = ttest2(bonus.group_1.raw, bonus.group_2.raw);

[worktime.h, worktime.p, worktime.ci, worktime.stats] = ttest2(worktime.group_1.raw, worktime.group_2.raw);

[duration.h, duration.p, duration.ci, duration.stats] = ttest2(duration.group_1.raw, duration.group_2.raw);

supp = struct('age',age,'gender',gender,'bonus',bonus,'stimCond',stimCond,'worktime',worktime,'duration',duration);

end

