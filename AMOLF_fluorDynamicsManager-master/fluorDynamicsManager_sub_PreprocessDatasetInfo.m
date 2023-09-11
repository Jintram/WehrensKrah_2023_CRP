
%% This script determines dataset parameters
% Path parameters are extracted from the main XLS configuration file, 
% then some key parameters are determined and the config file (also xls) 
% is loaded, to establish more parameters that are necessary to run the
% analysis and/or load the data/plots.
%
% Input that is required:
% - dataIdx
% - ..
%
% Parameters that are determined by this script are:
% 
% rootDir
% movieDirPluspositionDir
% configFileName
% identifierInXLS
% fullexperimentalpath
% movieDate
% theMovieDateOnly
% movieName
% dataFileName
% theDirectoryWithThePlots
% dateDir
% ... (maybe some more that I forgot to mention)



%% else, process
rootDir=allXLSdata{dataIdx,1}
movieDirPluspositionDir=allXLSdata{dataIdx,2}
configFileName=allXLSdata{dataIdx,3}
identifierInXLS=allXLSdata{dataIdx,4}
% Test dataset
% ===
%{
% basedir
stringp1  = 'F:\EXPERIMENTAL_DATA_2014-2015_m1\';
% experimental data and analysis file path
stringp2  ='2015-06-12_CRP_asc852-853_plasmids\pos1crop\'
% configfile file name (should be located in date dir)
stringp3 = 'configFileMadeLater_2015-06-12_pos1.xlsx'
% ID to recognize dataset
stringp4 = 'CRP_plasmids_WT_rCRP'; 
%}

%% Now convert to more convenient names
% ===
fullexperimentalpath  =[rootDir movieDirPluspositionDir];

%if ~(fullpathstringposition(end)=='\') 
%    fullpathstringposition(end+1)='\';
%end

% Extract information from the path, note that this relies on some
% conventions.
% ===
% Examine the full path
slashLocations = strfind(fullexperimentalpath,'\');
% Get date movie + name of experiment (referred to as movieDate)
movieDate = fullexperimentalpath(slashLocations(end-2)+1:slashLocations(end-1)-1)
% Get the date by itself
theMovieDateOnly   = movieDate(1:10)
% Get the position name 
movieName = fullexperimentalpath(slashLocations(end-1)+1:slashLocations(end)-1)

% Get the full path to the experiment ('date dir')
dateDir = fullexperimentalpath(1:slashLocations(end-1))
% Get the path to the directory that holds the experimental dir ('root dir')
rootDir = rootDir % should equal fullexperimentalpath(1:slashLocations(end-2))
% Get the config file name
configFileName           = configFileName;

%% Set up ourSettings struct required to run analysis

% Clear ourSettings and other config variables
clear -global settings ourSettings p 

% Clear output parameters
clear -global schnitzcells s_rm output

% Set up minimal set of variables in ourSettings struct
ourSettings.mypathname=         dateDir;
ourSettings.myconfigfilename=   configFileName
ourSettings.configfilepath = [ourSettings.mypathname ourSettings.myconfigfilename];

% Load the configfile
runsections='reloadfile'; % 'reloadfile' is better than 'loadfile' since it doesn't ask user where file is located
Schnitzcells_masterscript

%%

% Update some more variables
if isempty(strfind(upper(ourSettings.ASCnumber),'ASC'))
        ourSettings.ASCnumber=['asc' num2str(ourSettings.ASCnumber)];
end
ourSettings.myID = [ourSettings.myGroupID '_' ourSettings.ASCnumber];

% determine location of datafile
dataFileName = [dateDir 'outputandsettings_v2_' theMovieDateOnly '_' movieName '.mat'];
    % note that Schnitzcells_masterscript also writes to this file
    
% Note that the directory might have been moved or renamed since it was
% analyzed, so we're forcing to use the information from the master excel
% file over the information from the config excel file
ourSettings.rootDir   = rootDir;
ourSettings.movieDate = movieDate;
    % note that p.movieDate = ourSettings.movieDate(1:10)

% Also determine the name of the directory where plots are stored
theDirectoryWithThePlots = [ourSettings.mypathname theMovieDateOnly  '_' movieName '_' ourSettings.myID  '\'];




    