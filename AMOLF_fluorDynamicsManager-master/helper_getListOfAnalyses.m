



% Dataminer script


%% CRP project data directories

%% 14/15, microscope 1 
baseDir = 'F:\EXPERIMENTAL_DATA_2014-2015_m1\';
experimentalDirList = {
'2015-06-12_CRP_asc852-853_plasmids',...
'2015-06-13_CRP_asc854-855_plasmids',...
'2015-07-27_CRP_asc854-855_plasmids',...
'2015-10-20_CRP_asc841-842_plasmids',...
'2015-12-02_CRP_asc893-941_plasmids',...
'2015-12-09_CRP_asc893-941_plasmids_baddata_canremove'};

%% 14/15, microscope 2
baseDir = 'F:\EXPERIMENTAL_DATA_2014-2015_m2\';
experimentalDirList = {
'2015-09-15_CRP_asc893-894_plasmids',...
'2015-09-16_CRP_asc893-894_plasmids',...
'2015-09-17_CRP_asc893-894_plasmids',...
'2015-12-09_CRP_asc893-941_plasmids'};

%% 16, microscope 1
baseDir = 'G:\EXPERIMENTAL_DATA_2016\';
experimentalDirList = {
'2016-06-21_CRP_asc893_941_cAMP60uM',...
'2016-06-28_CRP_asc941_cAMP300uM_baddata',...
'2016-06-28_CRP_asc941_cAMP300uM_crashed_baddata',...
'2016-06-30_CRP_asc893-941_cAMP300uM',...
'2016-09-02_CRP_asc990lac_badpad_baddata',...
'2016-09-02_CRP_asc990m9lac',...
'2016-09-06_CRP_asc990lac_extrapressure',...
'2016-09-06_CRP_asc990lac_lesspressure',...
'2016-09-28_CRP_M9xyl_asc990',...
'2016-10-19_CRP_asc1004_cAMP800',...
'2016-10-19_CRP_asc1004_cAMP800_m2',...
'2016-11-01_CRP_cAMPchromoLOW80_ascASC1004',...
...'2016-11-02_CRP_snaps_CRPstar',...
'2016-11-03_CRP_cAMPchromoLOW80_asc1004',...
'2016-11-08_CRP_cAMPchromoHIGH5000_asc1004',...
'2016-11-17_CRP_cAMPchromoXyl_asc990_TOANALYZE',...
'2016-12-08_CRP_asc990_lac',...
'2016-12-13_CRP_asc990_galactose',...
'2016-12-17_CRP_asc990_glucose_TOANALYZE',...
'2016-12-21_CRP_asc1004_glucose_225uM-cAMP_TOANALYZE'};

%% 17, microscope 1
baseDir = 'H:\EXPERIMENTAL_DATA_2017\';
experimentalDirList = {
'2017-01-06_CRP_asc990_aKG_lac',...
'2017-01-12_CRP_asc990_mothermachine_lacgluc',...
'2017-02-02_CRP_asc990_OAA_pulsing_3',...
'2017-02-05_CRP_aKG_pulsing_notsuredata',...
'2017-02-10_CRP_OAA_pulsing_asc1004_500uM-cAMP',...
'2017-03-22_CRP_asc1004_cAMP_pulsing'};


%% Identify directories

count=0;

% Set experimentalDirList and baseDir above
positionStruct=struct; 
for dirIdx = 1:numel(experimentalDirList)
    
    experimentalDir = [experimentalDirList{dirIdx} '\'];

    experimentalDirContents = dir([baseDir experimentalDir]);
    
    % Create a run to identify crop dirs
    for ii=1:numel(experimentalDirContents)

        % go over directories
        if experimentalDirContents(ii).isdir & experimentalDirContents(ii).name(1)~='.'
           % check whether there's directories that have a crop suffix (i.e.
           % longer name)
           if (experimentalDirContents(ii).name(1:3) == 'pos') & (numel(experimentalDirContents(ii).name) > 5)

               disp(['Identified ' experimentalDirContents(ii).name ' as possible data directory.']);
               count=count+1;

               posDirName = experimentalDirContents(ii).name;

               positionStruct(count).baseDir = baseDir;
               positionStruct(count).fullExperimentalDir =  [experimentalDir posDirName '\'];

               % It's hard to automatically identify whether this folder has a fully analyzed dataset, so this will have to be done manually..                 

           end
        end

    end
end

%% Create tab delimited list (can be pasted to Excel)

disp('*** OUTPUT ***');
for ii = 1:numel(positionStruct)
    disp([positionStruct(ii).baseDir 9  positionStruct(ii).fullExperimentalDir]);
end
disp('*** END OUTPUT ***');











