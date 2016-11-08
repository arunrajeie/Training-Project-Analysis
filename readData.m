function [ data ] = readData( subjects, sessions )
%READDATA Reads data and creates summary statistics for specified subject(s) and session(s)
%   Input 1: array of N subject's whose data you wish to read
%       e.g. subjects = 1:100;
%       e.g. subjects = [1, 3, 8, 23, 34];
%   Input 2: array of M sessions of data you wish to read
%       e.g. sessions = 1:10;
%       e.g. sessions = 2:9;
%       e.g. sessions = [1,10];
%
%   Output: data structure containing sub structures for each 
%   subject and session. Within each session, there are sub structures only 
%   for the columns of the CSV file. If only one subject/session is
%   entered, then the data structure directly references the columns of the
%   CSV file without having subject and session as a sub structure.

% set paths for code we need
rootdir = fileparts(which('readData')); % Directory in which we put readData.m
addpath(rootdir); % Add rootdir path
dataLocation = fullfile(rootdir,'data'); % where our data comes from
fs = filesep;

% Read data from CSV
numSubjects = length(subjects);
numSessions = length(sessions);
if numSubjects == 1 && numSessions == 1 % Don't create sub structures for 1 sub/session
    fName = sprintf('subject%.3d_session%.2d.csv',subjects,sessions);
    dataPath = [dataLocation fs fName];
    data = readSubjectSession(dataPath);
else
    for sub = 1:numSubjects
        subject = sprintf('subject_%.3d',subjects(sub));
        for sesh = 1:numSessions
            session = sprintf('session_%.2d',sessions(sesh));
            fName = sprintf('subject%.3d_session%.2d.csv',subjects(sub),sessions(sesh));
            dataPath = [dataLocation fs fName];
            data.(subject).(session) = readSubjectSession(dataPath);
        end
    end
end

end

