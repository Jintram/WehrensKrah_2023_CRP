
% Notes;
%
% Set masterExcelFileLocation first (see below)
%
% see also
% \\storage01\data\AMOLF\users\wehrens\ZZ_EXPERIMENTAL_DATA\MICROSCOPE_EXPERIMENTS_shortcuts\CRP_plasmid_data
% also for a .docx file with a list of plasmid data.
%
% Optional parameters for this script:
% - NOSAVEPLEASE

% SET THIS PARAMETER:
% 
%{
% CRP dataset
masterExcelFileLocation = '\\storage01\data\AMOLF\users\wehrens\ZZ_EXPERIMENTAL_DATA\MICROSCOPE_OVERVIEW_AND_FIGURES\_projectfile_CRP.xlsx';
OUTPUTFOLDER = 'U:\THESIS\Thesis\ChapterX_CRP\Figures\matlabExport\';
%}
%{
% RIBOSOMAL DATA
masterExcelFileLocation = '\\storage01\data\AMOLF\users\wehrens\ZZ_EXPERIMENTAL_DATA\MICROSCOPE_OVERVIEW_AND_FIGURES\_projectfile_Ribosomes.xlsx';
OUTPUTFOLDER = '\\storage01\data\AMOLF\users\wehrens\THESIS\Thesis\ChapterX_Ribosomes\Figures\MatlabExport\CCs\';
%}




if ~exist('FIGURESVISIBLE','var')
    FIGURESVISIBLE='on'; % choose 'off' to not have figures visible
end

%% A flashy new script to handle all my data

% Important parameter names in this script are
% - collectedOutput (stores output of multiple datasets into one struct)
% - processedOutput (holds mean values and single values per group)
% - fieldNameDictPerGroup (auxiliary output holds additional info per group per dataset)


% Documentation (outdated?)
% 
% DATA STORAGE ===
% The general idea is that there is an excel sheet that contains the names
% of the folders that hold the data.
% 
% Note sure necessary ---> 
% There will also be a column with identifiers, that describe 
% similar experiments such that data can be 
% grouped later.
% <-----
% 
% DATA&ANALYSIS STORAGE POLICY ===
% Generally, analyzed data will be stored in .mat files (or otherwise), 
% preferably in the directory with the data.
% Overview plots will be stored in an output directory.
% DATA STRUCTURE
% ===
% Experimental datasets are stored in folders that are labeled with the
% date the experiment has taken place, plus an additional description
% string (the "date" folder.
% Each experiment has different subfolders, where positions that were
% observed under the microscope are stored. Each analysis is performed on
% these "sub"experiments. Each of these sub experiments also has its own
% Excel configuration file, which allow Matlab to process the data. These
% configuration files are stored in the experimental date folder.
% The experimental folders are organized in "base" directories, that are
% organized per year.
%
% ANALYSIS === 
% Optimally, the script, when ran, would loop over all datasets and check
% whether analyses were already performed beforehand. Per dataset, if this
% is not the case, the analysis would be run again.
% This could be implemented by saving the analysis result per dataset to
% the folder itself (with a name specific to this script) and later
% checking whether that file already exists.
%
% 


%% Load Excel file to get information on directories

if ~exist('masterExcelFileLocation','var')
    error('masterExcelFileLocation was not set, set this parameter first');
end

[ndata, theTextData, allXLSdata] = xlsread(masterExcelFileLocation,'dataset_list','B16:E200')
nrDataLines=size(allXLSdata,1);

disp('Excel file with overview of data loaded..');

%% Define which datasets should be plotted together and how
% Sets of data that you want to plot are defined by a few key parameters, 

%% Execute any of these script before executing (the rest of) this script
%{
% FOR THE CRP datasets
z_plotssets_plasmids1
z_plotssets_plasmids2
z_plotssets_chromoCRPs70_1
z_plotssets_chromoCRPs70_prime
z_plotssets_chromo_misc

% FOR THE ribosome datasets
z_plotssets_riboM9steady
%}

%%

%{
z_plotssets_testplots
%}

    
%% Show some colors
%{
figure; hold on;
for i=1:size(somecolors,1)
    plot(i,1,'o','Color',somecolors(i,:),'MarkerSize',10,'MarkerFaceColor',somecolors(i,:))
end
%}

%% Running analyses if necessary

fluorDynamicsManager_sub_runAnalyses

%% Create parameter that has indices of datasets grouped

% First find out which lines of the xls file match the groups of plots that we want
applicableIndices = {};
% Go over the different plot groups
for plotGroupIdx = 1:numel(IDENTIFIERSTOPLOT)
    % Now find all lines in the xls file of which the identifier match any
    % of the IDENTIFIERSTOPLOT{plotGroupIdx}
    applicableIndices{plotGroupIdx} = ...
    find(...
        cellfun(@(currentIdentifierXls) any( ...
            cellfun(@(currentIdentifierToPlot) strcmp(currentIdentifierToPlot,currentIdentifierXls), IDENTIFIERSTOPLOT{plotGroupIdx})...
            ),  {allXLSdata{:,4}})...
    );
end

%% Now go and fetch all the corresponding filepaths, first for single params
figurePaths = struct;

for groupIdx = 1:numel(applicableIndices)
    for plotIdx = 1:numel(applicableIndices{groupIdx}) 
        
       
        %% 
        
        % Set appropriate index
        dataIdx =  applicableIndices{groupIdx}(plotIdx);
        
        % Load configuration file and pre-process dataset info
        fluorDynamicsManager_sub_PreprocessDatasetInfo
        
        % determine which fields are of interest
        fluorDynamicsManager_sub_GetInterestingFieldsForThisDataSet
            % i.e. get parameterOfInterestList
        
        for paramIdx = 1:numel(parameterOfInterestList)
            % Get the path to the PDF and put into figure list
            figFilename = ['FIG_PDF_' parameterOfInterestList{paramIdx} '_small.fig'];
            completeFigurePath = [theDirectoryWithThePlots figFilename];
            figurePaths.PDF{paramIdx}{groupIdx}{plotIdx}=completeFigurePath;
            % Get the path to the branches plot
            figFilename = ['FIG_branches_' parameterOfInterestList{paramIdx} '_small.fig'];
            completeFigurePath = [theDirectoryWithThePlots figFilename];
            figurePaths.branches{paramIdx}{groupIdx}{plotIdx}=completeFigurePath;
            % Get the path to the CV over time plot
            figFilename = ['FIG_CVovertime_' parameterOfInterestList{paramIdx} '_small.fig'];
            completeFigurePath = [theDirectoryWithThePlots figFilename];
            figurePaths.CV{paramIdx}{groupIdx}{plotIdx}=completeFigurePath;    
        end
        figurePaths.singleOrDual.CV=1;
        figurePaths.singleOrDual.PDF=1;
        figurePaths.singleOrDual.branches=1;
        
        
    end
end

disp('==');
allParamString=cell2mat(arrayfun(@(x) [10 '- ' parameterOfInterestList{x}],1:numel(parameterOfInterestList),'UniformOutput',0));
disp(['Parameters that have been defined in last run: ' allParamString]);

allDualParamString=cell2mat(arrayfun(@(x) [10 '- ' parameterOfInterestDoubleCombinatorialListString{x}],1:numel(parameterOfInterestDoubleCombinatorialListString),'UniformOutput',0));
disp(['Parameters pairs that have been defined in last run: ' allDualParamString]);
disp('==');

%% The same can be done for cross-correlations
%figurePaths = struct;

fieldNameDictPerGroup = {};

for groupIdx = 1:numel(applicableIndices)
    for plotIdx = 1:numel(applicableIndices{groupIdx}) 
        
       
        %% 
        
        % Set appropriate index
        dataIdx =  applicableIndices{groupIdx}(plotIdx);
        
        % Load configuration file and pre-process dataset info
        fluorDynamicsManager_sub_PreprocessDatasetInfo
        
        % determine which fields are of interest
        fluorDynamicsManager_sub_GetInterestingFieldsForThisDataSet
            % i.e. get parameterOfInterestList, parameterOfInterestDoubleCombinatorialList
        
        % Save parameternames in convenient reference ("dictionary like") struct
        fieldNameDictPerGroup{groupIdx}{plotIdx}.paramNames = paramNames;
            
        for paramIdx = 1:numel(parameterOfInterestDoubleCombinatorialListString)
            
            %%
            
            % Get the path                
            figFilename = ['FIG_crosscorrs_' parameterOfInterestDoubleCombinatorialListString{paramIdx} '_small.fig'];
            completeFigurePath = [theDirectoryWithThePlots figFilename];

            % Organize this into the figure list
            figurePaths.CCs{paramIdx}{groupIdx}{plotIdx}=completeFigurePath;
            figurePaths.singleOrDual.CCs=2;
            
        end
        
    end
end

%% Then place the figures neatly tiled into a figure

% Let's go over the figurePaths
namesOfDataGathered = fieldnames(figurePaths);

for plotTypeIdx = 1:numel(namesOfDataGathered)
    
    % Note there's one field that's just for administration, skip that one
    if strcmp(namesOfDataGathered{plotTypeIdx},'singleOrDual'), continue, end
    
    %% Make (multiple) combined plots
    
    plotType = namesOfDataGathered{plotTypeIdx};    
    % Get info on whether this is single or dual parameter (set)
    singleOrDual = figurePaths.singleOrDual.(plotType);
    
    % Note that the # of plots in figurePaths.(plotType) should match the
    % parameter lists in either parameterOfInterestList or 
    % parameterOfInterestDoubleCombinatorialList
    if singleOrDual == 1
        if numel(figurePaths.(plotType)) ~= numel(parameterOfInterestList)
            error('Consistency check failed. Check integrity of figurePaths.');
        end
    elseif singleOrDual == 2
        if numel(figurePaths.(plotType)) ~= numel(parameterOfInterestDoubleCombinatorialListString)
            error('Consistency check failed. Check integrity of figurePaths.');
        end
    end
        
    %
    for paramIdx=1:numel(figurePaths.(plotType))        
        % Note that paramIdx will map to field names in
        % in parameterOfInterestList or 
        % parameterOfInterestDoubleCombinatorialList
        
        %% Now create for each parameter (or combination of parameters) an overview plot        

        % Determine the size of our subplot figure first
        nrPanels      = cellfun(@(x) numel(x), applicableIndices);
        totalNrPanels = sum(nrPanels);
        xNrPanels = min(ceil(sqrt(totalNrPanels)), max(nrPanels));
        yNrPanels = sum(ceil(nrPanels./xNrPanels));

        % Determine (y) size of figure
        theHeight = min(yNrPanels*4.8, 19.2);

        % Make figure
        if exist('h1','var'), if ishandle(h1), close(h1), end, end % close if handle already existed
        h1=figure('Visible',FIGURESVISIBLE); clf; hold on; % make new
        MW_makeplotlookbetter(8,[],[12.8 theHeight]./2,1); 

        % Necessary later
        groupLabels = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

        %
        panelLine=0;
        panelNr=0;
        for groupIdx = 1:numel(applicableIndices)
            panelLine=panelLine+1;
            panelColumn=0;
            for plotIdx = 1:numel(applicableIndices{groupIdx}) 

                %%

                
                % Determine which panel to use
                panelColumn=panelColumn+1;
                if panelColumn>xNrPanels
                    panelColumn=1; panelLine=panelLine+1;
                end
                panelNr=(panelLine-1)*xNrPanels+panelColumn

                % Get the path for the current figure
                currentFigure = figurePaths.(plotType){paramIdx}{groupIdx}{plotIdx};

                % Skip if figure doesn't exist
                if ~exist(currentFigure,'file')
                    
                end
                
                % Load and get handle for current figure
                hCurrentFig = openfig(currentFigure);
                set(0,'CurrentFigure',hCurrentFig);
                set(hCurrentFig,'Visible','off');
                ax1 = gca;
                fig1 = get(ax1,'children');

                % Switch to figure with multiple panels
                set(0,'CurrentFigure',h1);
                s1=subplot(yNrPanels,xNrPanels,panelNr);        
                copyobj(fig1,s1);

                % XLim and YLim are not copied, so copy manually
                xlim(ax1.XLim);
                ylim(ax1.YLim);

                % Turn of axes ticks if desired
                %set(gca,'YTickLabel',[]); %set(gca,'XTickLabel',[]); 

                % give a title if it is the first of a group
                if plotIdx==1
                    t=title(groupLabels(groupIdx));
                    set(t, 'horizontalAlignment', 'left');
                    set(t, 'units', 'normalized');
                    t1 = get(t, 'position');
                    %set(t, 'position', [0 t1(2) t1(3)]);
                    set(t, 'position', [0-.1*xNrPanels t1(2) t1(3)]); % This was a bit trial and error
                end

                % Set limits to axes if desired        
                if singleOrDual == 1
                    if exist('CUSTOMLIMITSPERPARAMGROUP','var')
                        if ~any(isnan(CUSTOMLIMITSPERPARAMGROUP{paramIdx}(1:2)))
                            xlim(CUSTOMLIMITSPERPARAMGROUP{paramIdx}(1:2));
                        end
                        if ~any(isnan(CUSTOMLIMITSPERPARAMGROUP{paramIdx}(3:4)))
                            ylim(CUSTOMLIMITSPERPARAMGROUP{paramIdx}(3:4));
                        end
                    end
                elseif singleOrDual == 2
                    if exist('CUSTOMLIMITSPERPARAMGROUPCCS','var')            
                    if ~any(isnan(CUSTOMLIMITSPERPARAMGROUPCCS{paramIdx}(1:2)))
                        xlim(CUSTOMLIMITSPERPARAMGROUPCCS{paramIdx}(1:2));
                    end
                    if ~any(isnan(CUSTOMLIMITSPERPARAMGROUPCCS{paramIdx}(3:4)))
                        ylim(CUSTOMLIMITSPERPARAMGROUPCCS{paramIdx}(3:4));
                    end
                    end
                end

            end

        end

        set(0,'CurrentFigure',h1);

        % Also give overview of the different conditions:
        disp(['===' 10 'FIGURE LEGEND']);
        for idIdx = 1:size(IDENTIFIERSTOPLOT,2)
            if singleOrDual == 1
                disp([groupLabels(idIdx) ': ' IDENTIFIERSTOPLOT{idIdx}{:} ', ' parameterOfInterestList{paramIdx}]); 
            elseif singleOrDual == 2
                disp([groupLabels(idIdx) ': ' IDENTIFIERSTOPLOT{idIdx}{:} ', ' parameterOfInterestDoubleCombinatorialListString{paramIdx}]); 
            end
        end
        disp('===');

        % Now give the plots some axes labels etc
        if singleOrDual == 1
            parameterName = parameterOfInterestList{paramIdx};
        elseif singleOrDual == 2
            parameterName  = parameterOfInterestDoubleCombinatorialListString{paramIdx};
            parameterName1 = parameterOfInterestDoubleCombinatorialList{paramIdx}{1};
            parameterName2 = parameterOfInterestDoubleCombinatorialList{paramIdx}{2};
         end
        
        % Find out human readable equivalent of (first) parameter name
        if parameterName(1) == 'd'
            humanReadableName = 'Fluorophore production (a.u./min)';
            shortName = ['production ' upper(parameterName(2))];
        elseif parameterName(1) == 'm'
            humanReadableName = 'Growth rate (doublings/hr)';
            shortName = 'growth';
        elseif any(strcmp(upper(parameterName(1)),{'G','R','C','Y'}))        
            humanReadableName = 'Fluorophore concentration (a.u./pixel)';
            shortName = ['concentration ' upper(parameterName(1))];
        end
        
        % And for second if applicable
        if singleOrDual == 2
            % Note that #1 is not necessary because we use the first letter
            % of the string in which the two are combined
            if parameterName2(1) == 'd'
                humanReadableName2 = 'Fluorophore production (a.u./min)';
                shortName2 = ['production ' upper(parameterName2(2))];
            elseif parameterName2(1) == 'm'
                humanReadableName2 = 'Growth rate (doublings/hr)';
                shortName2 = 'growth';
            elseif any(strcmp(upper(parameterName2(1)),{'G','R','C','Y'}))        
                humanReadableName2 = 'Fluorophore concentration (a.u./pixel)';
                shortName2 = ['concentration ' upper(parameterName2(1))];
            end
        end
        
        % Now make corresponding plot names
        if strcmp(plotType,'PDF')            
            myTitle = '';
            myxlabel= humanReadableName;
            myylabel= 'Probability (normalized)';
        elseif strcmp(plotType,'branches')
            myTitle = '';
            myxlabel= 'Time (minutes)';
            myylabel= humanReadableName; 
        elseif strcmp(plotType,'CV')
            myTitle = '';
            myxlabel= 'Time (minutes)';
            myylabel= ['CV of ' shortName];
        elseif strcmp(plotType,'CCs')
            myTitle = '';
            myxlabel=   'Delay (hrs)';
            myylabel=   ['R(' shortName ',' shortName2 ')']; % XXX
        end
        subtitle_mw(myTitle,myxlabel,myylabel);
        
        %% Now save the plot               
        
        fileName = [GROUPNAME '_' plotType '_' parameterName];
        
        if ~exist('NOSAVEPLEASE')
            saveas(h1,[OUTPUTFOLDER 'tif_' fileName '.tif']);
            saveas(h1,[OUTPUTFOLDER 'fig_' fileName '.fig']);
            saveas(h1,[OUTPUTFOLDER 'svg_' fileName '.svg']);
            saveas(h1,[OUTPUTFOLDER 'pdf_' fileName '.pdf']);
        end
        
    end
    
end
    
disp('Section all done');

%%%%%%%%%%%%%%%
%SINGLEORDUALNAMES = {'single','CC'};
%SINGLEPLOTTYPES   = {'branches' 'PDF' 'CVovertime'};
%{
% dual or single, i.e. single parameter, or comparison of two parameters
% like for CC
for singleOrDualPlot = 1:2  

    if singleOrDualPlot==1
        paramIndices = 1:numel(parameterOfInterestList);        
    elseif singleOrDualPlot==2
        paramIndices = 1:numel(parameterOfInterestDoubleCombinatorialList);                
    else, error('Not recognized');
    end
    
    for paramIdx=paramIndices    

        if singleOrDualPlot == 1
            plotType = SINGLEPLOTTYPES{paramIdx};
        elseif singleOrDualPlot == 2
            plotType = 'CCs'; % CCs branches PDF CVovertime
        end
        
        %%%
        % CODE WAS CUT HERE
        %%%
        
    end
        
end
%}

%% Now collect things from "output" struct that is saved per dataset.
% Part 1 of processing
% Put everything in collectedOutput

PARAMETERNAMESTOCOLLECT = {'growthMean','rateMean','concentrationMean','CV','CC'};
NAMESTOSTORETHEMWITH    = {'Growth','Production','Concentration','CV','CC'};
    % Note 1: when changing the # of parameters that are collected here,
    % note that they are all collected into a subplot overview, so be
    % careful!
    % Note 2: parameters are renamed according to NAMESTOSTORETHEMWITH and
    % stored in collectedOutput. CV should not be renamed as that it is
    % assumed later it'll still have that name. Names are also used for y
    % plot labels later.

collectedOutput=struct;

for groupIdx = 1:numel(applicableIndices)
    for plotIdx = 1:numel(applicableIndices{groupIdx}) 


        %% 

        clear output
        
        %collectedOutput.(currentParameterNameToCollect) = {{}};

        % Set appropriate index
        dataIdx =  applicableIndices{groupIdx}(plotIdx);

        % Load configuration file and pre-process dataset info
        fluorDynamicsManager_sub_PreprocessDatasetInfo
            % also determines dataFileName

        load(dataFileName,'output');

        for paramIdx = 1:numel(PARAMETERNAMESTOCOLLECT)

            currentParameterNameToCollect=PARAMETERNAMESTOCOLLECT{paramIdx};

            if ~iscell(output.(currentParameterNameToCollect))
                collectedOutput.(NAMESTOSTORETHEMWITH{paramIdx}){groupIdx}{plotIdx} = ...
                    output.(currentParameterNameToCollect);
            else                
                    for cellIdx = 1:numel(output.(currentParameterNameToCollect))
                        collectedOutput.([NAMESTOSTORETHEMWITH{paramIdx} '_' upper(fluorValues{cellIdx})]){groupIdx}{plotIdx} = ...
                            output.(currentParameterNameToCollect){cellIdx};
                    end
            end
            
        end

    end
end
    
disp('Went over ouput files and collected values.');

%% Process data on simple fields and put in processed output
% Part 2 of processing
% Part 1/2 of putting things in processedOutput

processedOutput = struct;

parametersToCollect = fieldnames(collectedOutput);
%{'growthMean','rateMean','concentrationMean'};

for paramIdx = 1:numel(parametersToCollect)
    
    
    %% if multiple cells are available then just label them _1, _2 etcetera    
    
    currentParameterNameToCollect=parametersToCollect{paramIdx};
    
    % Add an exclusions for things we want to analyze later and not here
    % (Like CV, CC)
    if any(cellfun(@(fieldName) strcmp(currentParameterNameToCollect,fieldName) ,{'CV','CC'}))
        continue
    end
    
    %%
    processedOutput.(currentParameterNameToCollect).labelLocations = []; 
    processedOutput.(currentParameterNameToCollect).technicalReplicateNrs = [];
    processedOutput.(currentParameterNameToCollect).allMeans = []; 
    processedOutput.(currentParameterNameToCollect).allSEMs = [];
    processedOutput.(currentParameterNameToCollect).allValues = {};
    for groupIdx = 1:numel(applicableIndices)

        % Get values
        allCurrentValues = [collectedOutput.(currentParameterNameToCollect){groupIdx}{:}];
            % Note that the last {:} collects data from the multiple plots
            % in the group

        % Mean, SEM, etc
        currentMean = mean(allCurrentValues);
        currentSEM = std(allCurrentValues)/sqrt(numel(allCurrentValues));
        currentTechnicalReplicateNr = numel(allCurrentValues);

        % Remove SEMs based on 1 replicate
        if currentTechnicalReplicateNr==1
            currentSEM=NaN;
        end

        processedOutput.(currentParameterNameToCollect).labelLocations(end+1) = groupIdx;
        processedOutput.(currentParameterNameToCollect).technicalReplicateNrs(end+1) = currentTechnicalReplicateNr;
        processedOutput.(currentParameterNameToCollect).allMeans(end+1)= currentMean;
        processedOutput.(currentParameterNameToCollect).allSEMs(end+1) = currentSEM;
        processedOutput.(currentParameterNameToCollect).allValues{end+1} = allCurrentValues;
    end

end
    
%{
% print values on top
for groupIdx = 1:numel(applicableIndices)
    myValueString=sprintf('%0.2f',allMeans(groupIdx));
    text(groupIdx, myYlim(2),myValueString,'HorizontalAlignment','center');
end

labelNames = arrayfun(@(x) IDENTIFIERSTOPLOT{x}{:}, 1:numel(IDENTIFIERSTOPLOT),'UniformOutput', 0);

set(gca,'XTick',labelLocations,'XTickLabel',labelNames,'TickLabelInterpreter','None');
xtickangle(90);

ylabel('Growth rate [dbl/hr]');
MW_makeplotlookbetter(10);
%}

%% Now process data for the slightly more complicated CV
% Part 3 of processing
% Part 2/2 of putting things in processedOutput

% note that nrFluorColors is determined by
% fluorDynamicsManager_sub_GetInterestingFieldsForThisDataSet, and this
% section assumes (1) this has been run before, and (2) they're all the
% same for your datasets.

% set the search terms, e.g.: mu, dG, G5
searchTermsForWhatToPlot = {'mu'};
whatToPlotCommonName = {'growth'};
for fluorIdx = 1:nrFluorColors
    currentFluor = upper(fluorValues{fluorIdx});
    searchTermsForWhatToPlot{end+1} = ['^d' currentFluor];
    whatToPlotCommonName{end+1} = ['production_' currentFluor];
    searchTermsForWhatToPlot{end+1} = ['^' currentFluor '5']; 
    whatToPlotCommonName{end+1} = ['concentration_' currentFluor];
end

% Now gather all the CVs per parameter
for paramIdx=1:numel(searchTermsForWhatToPlot)
    
    currentSearchTermForWhatToPlot = searchTermsForWhatToPlot{paramIdx};
    currentwhatToPlotCommonName = whatToPlotCommonName{paramIdx};

    currentFieldToMake = ['CV_' currentwhatToPlotCommonName];
    
    processedOutput.(currentFieldToMake).labelLocations=[];
    processedOutput.(currentFieldToMake).technicalReplicateNrs = [];
    processedOutput.(currentFieldToMake).allMeans = []; 
    processedOutput.(currentFieldToMake).allSEMs = [];
    processedOutput.(currentFieldToMake).allValues = {};
    for groupIdx = 1:numel(applicableIndices)

        % Gather the values    
        allCurrentValues = [];
        for plotIdx=1:numel(collectedOutput.CV{groupIdx})
            % First determine how the value is called here..
            theFieldNames = fieldnames(collectedOutput.CV{groupIdx}{plotIdx});
            theFieldIdx   = find(arrayfun(@(x) any(regexp(theFieldNames{x},currentSearchTermForWhatToPlot)), 1:numel(theFieldNames)));
            whatToPlotName = theFieldNames{theFieldIdx};
            
            % Either take the mean of the last 10 values            
            % allCurrentValues(end+1) = collectedOutput.CV{groupIdx}{plotIdx}.(whatToPlotName).meanCoefficientOfVariationLast10;
            
            % Or calculate the weighed average
            %times  = collectedOutput.CV{groupIdx}{plotIdx}.(whatToPlotName).time; 
            values = collectedOutput.CV{groupIdx}{plotIdx}.(whatToPlotName).coefficientOfVariationOverTime;
            nrvalues = collectedOutput.CV{groupIdx}{plotIdx}.(whatToPlotName).numberOfValues;
            allCurrentValues(end+1) = sum(values.*nrvalues)/sum(nrvalues);
        end
        %{ 
        % Get values assuming all fields have the same name
        allCurrentValues = arrayfun(@(x) collectedOutput.CV{groupIdx}{x}.(whatToPlotName).meanCoefficientOfVariationLast10, 1:numel(collectedOutput.CV{groupIdx}) );
            % note that CV is determined in a complicated way
        %}

        % Mean, SEM, etc
        currentMean = mean(allCurrentValues);
        currentSEM = std(allCurrentValues)/sqrt(numel(allCurrentValues));
        currentTechnicalReplicateNr = numel(allCurrentValues);

        % Remove SEMs based on 1 replicate
        if currentTechnicalReplicateNr==1
            currentSEM=NaN;
        end

        % Store for later plotting    
        processedOutput.(currentFieldToMake).labelLocations(end+1) = groupIdx;
        processedOutput.(currentFieldToMake).technicalReplicateNrs(end+1) = currentTechnicalReplicateNr;    
        processedOutput.(currentFieldToMake).allMeans(end+1)= currentMean;
        processedOutput.(currentFieldToMake).allSEMs(end+1) = currentSEM;
        processedOutput.(currentFieldToMake).allValues{end+1} = allCurrentValues;
    end

end
    
plotType = 'CV'; % mu, CV


%{
% Find the field name based on simple search term
theFieldNames = fieldnames(collectedOutput.CV{1}{1});
%theFieldIdx   = find(arrayfun(@(x) any(strfind(theFieldNames{x},searchTermForWhatToPlot)), 1:numel(theFieldNames)));
theFieldIdx   = find(arrayfun(@(x) any(regexp(theFieldNames{x},currentSearchTermForWhatToPlot)), 1:numel(theFieldNames)));
whatToPlotName = theFieldNames{theFieldIdx};
%}

%% General way of plotting data per group..
% Now plot (part 4)

%plotType = 'CV'; % mu, CV

someColors = linspecer(numel(IDENTIFIERSTOPLOT));

myOverviewParametersToPlot = fieldnames(processedOutput);

h1=figure; clf; hold on;

numelMyOverviewParametersToPlot=numel(myOverviewParametersToPlot);
if nrFluorColors==1
    subplotXNr=3;
    sizeY=6.4;
    FontSize=10;
    myLineWidth=2;
else
    subplotXNr=3;
    sizeY=4.5;
    %sizeY=4.8;
    FontSize=8;
    myLineWidth=1;
end
subplotYNr = ceil(numelMyOverviewParametersToPlot/subplotXNr);
for paramIdx = 1:numelMyOverviewParametersToPlot

    %
            
    plotType = myOverviewParametersToPlot{paramIdx};
        
    % Now make an overview subplot
    %subtightplot(subplotYNr,3,paramIdx,0.01,0,0)
    if numel(IDENTIFIERSTOPLOT)>4  
        subtightplot(subplotYNr,subplotXNr,paramIdx,[0.2,0.05],[0.2 0.05],[0.1 0.1]); hold on;
    elseif nrFluorColors>1
        subtightplot(subplotYNr,subplotXNr,paramIdx,[0.12,0.12],[0.10 0.05],[0.1 0.1]); hold on;
    else
        subplot(subplotYNr,subplotXNr,paramIdx); hold on;
    end
    
    bar(1:numel(applicableIndices), processedOutput.(plotType).allMeans,'FaceColor',[.7 .7 .7]);
    %errorbar(1:numel(applicableIndices), allMeans, allSEMs,'k','LineStyle','none','LineWidth',2);

    myYlim=[0 max([processedOutput.(plotType).allMeans+processedOutput.(plotType).allSEMs, [processedOutput.(plotType).allValues{:}]])*2];

    % put printed values and single datapoints
    for groupIdx = 1:numel(applicableIndices)

        % single data points
        plot(ones(1,numel(processedOutput.(plotType).allValues{groupIdx})).*groupIdx,...
                processedOutput.(plotType).allValues{groupIdx},'o',...
                'LineWidth',myLineWidth,...
                'Color',someColors(groupIdx,:));%,...
                %'MarkerFaceColor',someColors(groupIdx,:)); 

        % write values on top
        %myValueString=sprintf('%0.2f',processedOutput.(plotType).allMeans(groupIdx));
        myValueString=sprintf('%0.2g',processedOutput.(plotType).allMeans(groupIdx));
        
        % Horizontal text at the top of the plot
        %text(groupIdx, myYlim(2),myValueString,'HorizontalAlignment','center');        
        
        % Vertical text at the bottom of each bar (and values)
        %yvalue = max([processedOutput.(plotType).allMeans(groupIdx), processedOutput.(plotType).allValues{groupIdx}]);        
        yvalue  = myYlim(2)*0.55;
        text(groupIdx, yvalue,myValueString,'HorizontalAlignment','left','Rotation',90);
    end

    % cosmetics
    ylim(myYlim);
    
    % get xtick labels
    if exist('HUMANREADABLENAMESFORGROUPS','var')
        % simply HUMANREADABLENAMESFORGROUPS, but replace comma for enter
        %labelNames = cellfun(@(value) strrep(value,', ',10),HUMANREADABLENAMESFORGROUPS,'UniformOutput',0);
        labelNames = HUMANREADABLENAMESFORGROUPS;
    else
        labelNames = arrayfun(@(x) IDENTIFIERSTOPLOT{x}{:}, 1:numel(IDENTIFIERSTOPLOT),'UniformOutput', 0);
    end
    % set xtick labels
    set(gca,'XTick',processedOutput.(plotType).labelLocations,'XTickLabel',labelNames,'TickLabelInterpreter','None');
    xtickangle(90);

    MW_makeplotlookbetter(10);

    if strcmp(plotType,'mu')
        ylabel('Growth rate [dbl/hr]');
    elseif strcmp(plotType,plotType)
        ylabel(strrep(plotType,'_',' '));%,'Interpreter','none');
    else
        ylabel('[unkown parameter]');
    end

    MW_makeplotlookbetter(FontSize,[],[12.8 subplotYNr*sizeY]/2,1)
    
end


fileName = [GROUPNAME '_overview_means'];

if ~exist('NOSAVEPLEASE')
    saveas(h1,[OUTPUTFOLDER 'tif_' fileName '.tif']);
    saveas(h1,[OUTPUTFOLDER 'fig_' fileName '.fig']);
    saveas(h1,[OUTPUTFOLDER 'svg_' fileName '.svg']);
    saveas(h1,[OUTPUTFOLDER 'pdf_' fileName '.pdf']);
end

%plotting CV
%output.CV.(theFieldNames).meanCoefficientOfVariationLast10

%% Make a modified plot, combining some of the parameters together

h2=figure(2); clf; hold on;

alldatax=[];alldatay=[];
for paramIdx = 1:numelMyOverviewParametersToPlot
    %myYlim=[0 max([processedOutput.(plotType).allMeans+processedOutput.(plotType).allSEMs, [processedOutput.(plotType).allValues{:}]])*2];

    % put printed values and single datapoints
    for groupIdx = 1:numel(applicableIndices)

        %thedts = myTimeBetweenShots(applicableIndices{groupIdx});
            % proved not necessary
                
        xdata = processedOutput.('Production_Y').allValues{groupIdx};
        ydata = processedOutput.('Concentration_Y').allValues{groupIdx} .* log(2)/60.*processedOutput.('Growth').allValues{groupIdx};
        
        % single data points
        plot(xdata,...
                ydata,...
                'o',...
                'LineWidth',myLineWidth,...
                'Color',someColors(groupIdx,:),...
                'MarkerFaceColor',someColors(groupIdx,:));
           
       alldatax = [alldatax xdata];
       alldatay = [alldatay ydata];
        
    end
    
    allvalues = [alldatax alldatay];
    plot([min(allvalues),max(allvalues)],[min(allvalues),max(allvalues)],'k-');
    
    % fit 1
    p=polyfit(alldatax,alldatay,1);
    xToFit = linspace(min(allvalues),max(allvalues),3);
    plot(xToFit,xToFit*p(1)+p(2),'k--')
    
    % fit 2
    offset = mean(alldatay-alldatax);
    plot([min(allvalues),max(allvalues)],[min(allvalues),max(allvalues)]+offset,'k:');
end

xlabel('Production Y [a.u./min]');
ylabel(['Concentration  times growth rate' 10 '[a.u./(px*min]']);

%% Make a text files as base for caption

fileID = fopen([OUTPUTFOLDER 'caption_' GROUPNAME '.txt'],'w');

fprintf(fileID,['===' 10 'FIGURE LEGEND' 10]);
for idIdx = 1:size(IDENTIFIERSTOPLOT,2)
    if singleOrDual == 1
        fprintf(fileID,[groupLabels(idIdx) ': ' IDENTIFIERSTOPLOT{idIdx}{:} ', ' parameterOfInterestList{paramIdx} 10]); 
    elseif singleOrDual == 2
        fprintf(fileID,[groupLabels(idIdx) ': ' IDENTIFIERSTOPLOT{idIdx}{:} ', ' parameterOfInterestDoubleCombinatorialListString{paramIdx} 10]); 
    end
end
fprintf(fileID,['===' 10]);

fclose(fileID);

%% Plot CV lines separately but color for group

%{
searchTermForWhatToPlot = 'mu'; % e.g.: mu, dG, G5

% Find the field name based on simple search term
theFieldNames = fieldnames(collectedOutput.CV{1}{1});
theFieldIdx   = find(arrayfun(@(x) any(strfind(theFieldNames{x},searchTermForWhatToPlot)), 1:numel(theFieldNames)));
whatToPlotName = theFieldNames{theFieldIdx};

% Now plot
h1=figure; clf; hold on;

% cosmetics
someColors = linspecer(numel(applicableIndices));

labelLocations = []; technicalReplicateNrs = [];
allMeans = []; allSEMs = [];
myLegendLines = []; myAxes = [];
allWeighedAveragesCV={};
for groupIdx = 1:numel(applicableIndices)

    allWeighedAveragesCV{groupIdx}=[];
    for lineIdx=1:numel(collectedOutput.CV{groupIdx})
    
        % Get values
        currentValuesToPlot = ...
            collectedOutput.CV{groupIdx}{lineIdx}.(whatToPlotName).coefficientOfVariationOverTime;
        currentTimePoints = ...
            collectedOutput.CV{groupIdx}{lineIdx}.(whatToPlotName).time;
        
            
        ax=subplot(numel(applicableIndices),1,groupIdx); hold on;
        l=plot(currentTimePoints,currentValuesToPlot,'Color',someColors(groupIdx,:),'LineWidth',2);
        
        % additionally, we could calculate a weighed average..
        currentNumberOfValues = ...
            collectedOutput.CV{groupIdx}{lineIdx}.(whatToPlotName).numberOfValues;
        currentWeighedAverage = sum(currentValuesToPlot.*currentNumberOfValues)/sum(currentNumberOfValues);
        
        allWeighedAveragesCV{groupIdx}(end+1) = currentWeighedAverage;
    end
    
    myLegendLines(end+1) = l;
    myAxes(end+1)=ax;       
    
end

ylabel(plotType);
xlabel('Time (hrs)');
%}

%% Now create legend from this figure

%{
labelNames = arrayfun(@(x) IDENTIFIERSTOPLOT{x}{:}, 1:numel(IDENTIFIERSTOPLOT),'UniformOutput', 0);
legendHandle=legend(myLegendLines,labelNames,'Interpreter','None','Location','NorthOutside');

saveLegendToImage(h1, legendHandle, myAxes);


%{
% Create separate legend plot
hLegend=figure(); clf; hold on;
myLegendLinesPrime=[];
for groupIdx = 1:numel(applicableIndices)
    lprime=plot([groupIdx,groupIdx],[1,2],'Color',someColors(groupIdx,:),'LineWidth',2);
    myLegendLinesPrime(end+1)=lprime;    
end
labelNames = arrayfun(@(x) IDENTIFIERSTOPLOT{x}{:}, 1:numel(IDENTIFIERSTOPLOT),'UniformOutput', 0);
legendHandle=legend(myLegendLinesPrime,labelNames,'Interpreter','None','Location','NorthOutside');

set(myLegendLinesPrime, 'visible', 'off');
set(gca, 'visible', 'off');

%hLegend=figure();
%copyobj(legendHandle,hLegend);
%}
%}

%% Now plot weighed averages

%{

someColors = linspecer(numel(applicableIndices));

h1=figure; clf; hold on;
for groupIdx = 1:numel(applicableIndices)
    
    bar(groupIdx,mean(allWeighedAveragesCV{groupIdx}),'FaceColor',[.7 .7 .7],'EdgeColor','none');
    
    plot(ones(1,numel(allWeighedAveragesCV{groupIdx})).*groupIdx,...
            allWeighedAveragesCV{groupIdx},'o',...
            'LineWidth',2,...
            'Color',someColors(groupIdx,:));            
end

labelNames = arrayfun(@(x) IDENTIFIERSTOPLOT{x}{:}, 1:numel(IDENTIFIERSTOPLOT),'UniformOutput', 0);
set(gca,'XTick',1:numel(applicableIndices),'XTickLabel',labelNames,'TickLabelInterpreter','None');
xtickangle(90);

ylabel(plotType);

%}

%% Plot CV against growth rate

try
    
h1=figure; clf; hold on;
for groupIdx = 1:numel(applicableIndices)
    
    plot([collectedOutput.(currentParameterNameToCollect){groupIdx}{:}],...        
         collectedOutput.CV{groupIdx}{:},...
            'o',...
            'LineWidth',2,...
            'Color',someColors(groupIdx,:)); 
         
end

xlabel('Growth rate [dbl/hr]');
ylabel('CV');

catch
    
    warning('Note: this does not work since there are multiple CV fields');
    % try processedOutput, but think about what to plot
    
end

%% Now plot combined CCs, both p-mu and C-mu, for both fluor colors

nrOfGroups = numel(IDENTIFIERSTOPLOT);
colorsPerGroup = linspecer(nrOfGroups);

groupLabels='ABCDEFGHIJKLMNOPQRSTUVWXYZ';

set_fieldNames = ...
    {{'concentration','muWithConcentration'},...
     {'rate','muWithRate'}};        
set_linecolors = ...
                {colorsPerGroup,
                zeros(nrOfGroups,3)};

for fluorIdx = 1:nrFluorColors
    %%
    h1=figure('Visible',FIGURESVISIBLE); clf; 
    
    set_fluorIndices = {[fluorIdx,fluorIdx],...
                        [fluorIdx,fluorIdx]};              

    totalNrPanels=nrOfGroups;
    yNrPanels = ceil(sqrt(totalNrPanels));
    xNrPanels = ceil(totalNrPanels./yNrPanels);
    if xNrPanels==1, xNrPanels=yNrPanels; yNrPanels=1; end % avoid one ugly outcome
        
    meanLineHandleCollection=cell(1,nrOfGroups);
    for setIdx = 1:numel(set_fieldNames)

        fieldNames = set_fieldNames{setIdx};
        fluorIndices = set_fluorIndices{setIdx};
        linecolors = set_linecolors{setIdx};

        for groupIdx = 1:nrOfGroups

            subplot(yNrPanels,xNrPanels,groupIdx); hold on;

            collectedCCs = ...
                    collectedOutput.CC{groupIdx};
            fieldNameDict = ...
                    fieldNameDictPerGroup{groupIdx};

            [meanLineHandle,individualLineHandles] = ...
                fluorDynamicsManager_sub_meanCCforGroup(h1,collectedCCs,fieldNameDict,fieldNames,fluorIndices,linecolors(groupIdx,:))

            meanLineHandleCollection{groupIdx}(end+1) = meanLineHandle;

        end    

    end

    % cosmetics
    sizeY = 19.2/3;
    for groupIdx = 1:nrOfGroups
        subplot(yNrPanels,xNrPanels,groupIdx);    

        uistack(meanLineHandleCollection{groupIdx},'top');        

        t=title([groupLabels(groupIdx) ': ' HUMANREADABLENAMESFORGROUPS{groupIdx}]);
        set(t, 'horizontalAlignment', 'left');
        set(t, 'units', 'normalized');
        t1 = get(t, 'position');
        %set(t, 'position', [0 t1(2) t1(3)]);
        set(t, 'position', [0-.1 t1(2) t1(3)]); % This was a bit trial and error

        subtitle_mw('','Delay (hrs)','Correlation')

        MW_makeplotlookbetter(10,[],[12.8 yNrPanels*sizeY]/2,1)
    end

    fileName = [GROUPNAME '_overview_correlations_CmuPmu_' upper(fluorValues{fluorIdx})];

    if ~exist('NOSAVEPLEASE')
        saveas(h1,[OUTPUTFOLDER 'tif_' fileName '.tif']);
        saveas(h1,[OUTPUTFOLDER 'svg_' fileName '.svg']);
        saveas(h1,[OUTPUTFOLDER 'fig_' fileName '.fig']);
        saveas(h1,[OUTPUTFOLDER 'pdf_' fileName '.pdf']);
    end
    
end

%% Now plot combined autocorrelations, combine for groups, split for params

nrOfGroups = numel(IDENTIFIERSTOPLOT);
colorsPerGroup = linspecer(nrOfGroups);

groupLabels='ABCDEFGHIJKLMNOPQRSTUVWXYZ';

% sets of input in cells
set_fieldNames = ...
    {{'muWithConcentration','muWithConcentration'}}
set_hrnames = {'Growth rate'};
set_fluorIndices = {[fluorIdx,fluorIdx]};                        
for fluorIdx = 1:nrFluorColors
     set_fieldNames{end+1} = {'concentration','concentration'};     
     set_hrnames{end+1} = ['Concentration ' upper(fluorValues{fluorIdx})];
     set_fluorIndices{end+1} = [fluorIdx,fluorIdx]; 
     
     set_fieldNames{end+1} = {'rate','rate'};               
     set_hrnames{end+1} = ['Production ' upper(fluorValues{fluorIdx})];
     set_fluorIndices{end+1} = [fluorIdx,fluorIdx]; 
end


totalNrPanels=numel(set_fieldNames);
xNrPanels = ceil(sqrt(totalNrPanels));
yNrPanels = sum(ceil(totalNrPanels./xNrPanels));         

% Plot the figure
h1=figure('Visible',FIGURESVISIBLE); clf; 

meanLineHandleCollection=cell(1,numel(set_fieldNames));
myAxes=[];
for setIdx = 1:numel(set_fieldNames)

    %%
    fieldNames = set_fieldNames{setIdx};
    fluorIndices = set_fluorIndices{setIdx};

    for groupIdx = 1:nrOfGroups

        myAxes(end+1) = subplot(xNrPanels,yNrPanels,setIdx); hold on;

        collectedCCs = ...
                collectedOutput.CC{groupIdx};
        fieldNameDict = ...
                fieldNameDictPerGroup{groupIdx};

        [meanLineHandle,individualLineHandles] = ...
            fluorDynamicsManager_sub_meanCCforGroup(h1,collectedCCs,fieldNameDict,fieldNames,fluorIndices,colorsPerGroup(groupIdx,:))

        meanLineHandleCollection{setIdx}(end+1) = meanLineHandle;

    end    

end

% cosmetics
subtitle_mw('','Delay (hrs)','Autocorrelation')

sizeY = 19.2/3;
for setIdx = 1:numel(set_fieldNames)
    
    subplot(xNrPanels,yNrPanels,setIdx);    

    uistack(meanLineHandleCollection{setIdx},'top');        

    t=title([groupLabels(setIdx) ': ' set_hrnames{setIdx}]);
    set(t, 'horizontalAlignment', 'left');
    set(t, 'units', 'normalized');
    t1 = get(t, 'position');
    %set(t, 'position', [0 t1(2) t1(3)]);
    set(t, 'position', [0-.1 t1(2) t1(3)]); % This was a bit trial and error    

    MW_makeplotlookbetter(10,[],[12.8 yNrPanels*sizeY]/2,1)
end

fileName = [GROUPNAME '_overview_autocorrelations'];%_' set_fieldNames{setIdx}{1} '_' upper(fluorValues{fluorIdx})];
if ~exist('NOSAVEPLEASE')
    saveas(h1,[OUTPUTFOLDER 'tif_' fileName '.tif']);
    saveas(h1,[OUTPUTFOLDER 'svg_' fileName '.svg']);
    saveas(h1,[OUTPUTFOLDER 'fig_' fileName '.fig']);
    saveas(h1,[OUTPUTFOLDER 'pdf_' fileName '.pdf']);
end

%% now make and save a separate legend
legendHandle = legend(meanLineHandleCollection{1}, HUMANREADABLENAMESFORGROUPS);
saveLegendToImage(h1, legendHandle, myAxes);
fileName = [GROUPNAME '_overview_autocorrelations_legend'];

if ~exist('NOSAVEPLEASE')
    saveas(h1,[OUTPUTFOLDER 'tif_' fileName '.tif']);
    saveas(h1,[OUTPUTFOLDER 'svg_' fileName '.svg']);
    saveas(h1,[OUTPUTFOLDER 'fig_' fileName '.fig']);
    saveas(h1,[OUTPUTFOLDER 'pdf_' fileName '.pdf']);
end

%% Make scatter plots

HUMANREADABLENAMESFORGROUPSpathString = ...
    strrep(strrep(HUMANREADABLENAMESFORGROUPS,' ','_'),',','');
    

% Note that some previous sections (like gathering the file paths and
% cross-corr paths) should be execute before

% Set up which cases we want to run this analysis for
casesFluorIdxs={};
casesFluorLetters ={};
casesFieldNames={};
for fluorIdx=1:nrFluorColors
    % concentration vs. growth
    casesFluorIdxs{end+1} = [fluorIdx fluorIdx];
    casesFluorLetters{end+1} = {upper(fluorValues{casesFluorIdxs{end}(1)}) upper(fluorValues{casesFluorIdxs{end}(2)})};
    casesFieldNames{end+1} = {'concentration', 'muWithConcentration'};
    
    % production (=rate) vs. growth
    casesFluorIdxs{end+1} = [fluorIdx fluorIdx];
    casesFluorLetters{end+1} = {upper(fluorValues{casesFluorIdxs{end}(1)}) upper(fluorValues{casesFluorIdxs{end}(2)})};
    casesFieldNames{end+1} = {'rate', 'muWithRate'};
    
    % concentration vs. production (fluor1->fluor1)
    casesFluorIdxs{end+1} = [fluorIdx fluorIdx];
    casesFluorLetters{end+1} = {upper(fluorValues{casesFluorIdxs{end}(1)}) upper(fluorValues{casesFluorIdxs{end}(2)})};
    casesFieldNames{end+1} = {'concentrationWithRate', 'rate'};
    
    for fluorIdx2 = (fluorIdx+1):nrFluorColors
        % production vs. production
        casesFluorIdxs{end+1} = [fluorIdx2 fluorIdx];
        casesFluorLetters{end+1} = {upper(fluorValues{casesFluorIdxs{end}(1)}) upper(fluorValues{casesFluorIdxs{end}(2)})};
        casesFieldNames{end+1} = {'rate', 'rate'};
    
        % concentration vs. concentration
        casesFluorIdxs{end+1} = [fluorIdx2 fluorIdx];
        casesFluorLetters{end+1} = {upper(fluorValues{casesFluorIdxs{end}(1)}) upper(fluorValues{casesFluorIdxs{end}(2)})};
        casesFieldNames{end+1} = {'concentration', 'concentration'};
        
        % concentration vs. production (fluor1->fluor2)
        casesFluorIdxs{end+1} = [fluorIdx fluorIdx2];
        casesFluorLetters{end+1} = {upper(fluorValues{casesFluorIdxs{end}(1)}) upper(fluorValues{casesFluorIdxs{end}(2)})};
        casesFieldNames{end+1} = {'concentrationWithRate', 'rate'};
        
        % concentration vs. production (fluor2->fluor1)
        casesFluorIdxs{end+1} = [fluorIdx2 fluorIdx];
        casesFluorLetters{end+1} = {upper(fluorValues{casesFluorIdxs{end}(1)}) upper(fluorValues{casesFluorIdxs{end}(2)})};
        casesFieldNames{end+1} = {'concentrationWithRate', 'rate'};
    end
    
end

% run over different cases to plot
% This is a bit redundant with the schnitzcells loading, but it was the
% most lazy way to do it...
for caseIdx = 1:numel(casesFieldNames)

    % get information on the case we're going to plot now
    fluorIndices=casesFluorIdxs{caseIdx};
    currentFluorLetters = casesFluorLetters{caseIdx};
    currentFieldNames = casesFieldNames{caseIdx};

    % close figures from previous round if applicable
    if exist('hScatter','var'), if ishandle(hScatter), close(hScatter); end, end
    
    % create figures
    hScatter=[];
    for groupIdx = 1:numel(applicableIndices)
        hScatter(groupIdx)=figure('Visible',FIGURESVISIBLE); clf; hold on; 
    end

    % get some colors
    somecolors = linspecer(numel(applicableIndices));

    % go over groups
    meanLineHandles={};contourHandles={};averagedYHandles={};
    suggxlims1=[]; suggxlims2=[]; suggylims1=[]; suggylims2=[];
    for groupIdx = 1:numel(applicableIndices)

        %% plot all plots within that group (for this case)
        hCurrent = hScatter(groupIdx);

        meanLineHandles{groupIdx}=[];
        contourHandles{groupIdx}=[];
        averagedYHandles{groupIdx}=[];    
        for plotIdx = 1:numel(applicableIndices{groupIdx}) 


            %% 

            clear output

            %collectedOutput.(currentParameterNameToCollect) = {{}};

            % Set appropriate index
            dataIdx =  applicableIndices{groupIdx}(plotIdx);

            % Load configuration file and pre-process dataset info
            fluorDynamicsManager_sub_PreprocessDatasetInfo
                % also determines dataFileName

            load(dataFileName,'s_rm');

            %%                
            lineColor = somecolors(groupIdx,:);

            theRawfieldNames = {fieldNameDictPerGroup{groupIdx}{plotIdx}.paramNames.(currentFieldNames{1}){fluorIndices(1)},...
                                fieldNameDictPerGroup{groupIdx}{plotIdx}.paramNames.(currentFieldNames{2}){fluorIndices(2)}};

            % actually make a plot
            [meanLineHandle,contourHandle,individualLineHandles,averagedYHandle,suggxlim,suggylim] = fluorDynamicsManager_sub_scatterCloudsForGroup(hCurrent,s_rm,theRawfieldNames,fluorIndices,lineColor)        
            meanLineHandles{groupIdx}(end+1)=meanLineHandle;
            contourHandles{groupIdx}(end+1)=contourHandle;
            if ~isempty(averagedYHandle)
                averagedYHandles{groupIdx}(end+1)=averagedYHandle;                
            end

            suggxlims1(end+1)=suggxlim(1);
            suggxlims2(end+1)=suggxlim(2);
            suggylims1(end+1)=suggylim(1);
            suggylims2(end+1)=suggylim(2);
            
        end

        xlim([min(suggxlims1) max(suggxlims2)]);
        ylim([min(suggylims1) max(suggylims2)]);        
        
        uistack(contourHandles{groupIdx},'top');        
        uistack(meanLineHandles{groupIdx},'top');
        uistack(averagedYHandles{groupIdx},'top');

        set(averagedYHandle,'LineWidth',1);
        set(contourHandles{groupIdx},'LineWidth',1);

        % Some cosmetics
        if strcmp(currentFieldNames{1}, 'concentration')
            xlabel(['Concentration ' currentFluorLetters{1} ' (a.u.)']);
        elseif strcmp(currentFieldNames{1}, 'rate')
            xlabel(['Production ' currentFluorLetters{1} ' (a.u.)']);
        end
        
        if strcmp(currentFieldNames{2}, 'concentration')
            ylabel(['Concentration ' currentFluorLetters{2} ' (a.u.)']);
        elseif strcmp(currentFieldNames{2}, 'rate')
            ylabel(['Production ' currentFluorLetters{2} ' (a.u.)']);
        else
            % assume it's growth if neither concentration or rate 
            ylabel('Growth rate (doublings/hr)');
        end
        

        MW_makeplotlookbetter(10,[],[12.8/2, 19.2/3]/2,1);

        % Save the plots
        if ~exist('NOSAVEPLEASE')
            fileName = [GROUPNAME '_scatters_' HUMANREADABLENAMESFORGROUPSpathString{groupIdx} '_' currentFieldNames{2} currentFluorLetters{2} '_' currentFieldNames{1} currentFluorLetters{1} ];%_' set_fieldNames{setIdx}{1} '_' upper(fluorValues{fluorIdx})];        
            saveas(hCurrent,[OUTPUTFOLDER 'tif_' fileName '.tif']);
            saveas(hCurrent,[OUTPUTFOLDER 'svg_' fileName '.svg']);
            saveas(hCurrent,[OUTPUTFOLDER 'fig_' fileName '.fig']);
            saveas(hCurrent,[OUTPUTFOLDER 'pdf_' fileName '.pdf']);
        end


    end
    
end

disp('Done making all scatters..');

%% Special case: production rate divided by mu, vs concentration
%
% Note that some previous sections (like gathering the file paths and
% cross-corr paths) should be execute before

HUMANREADABLENAMESFORGROUPSpathString = ...
    strrep(strrep(HUMANREADABLENAMESFORGROUPS,' ','_'),',','');
    
% run over different cases to plot
% This is a bit redundant with the schnitzcells loading, but it was the
% most lazy way to do it...
for fluorIdx = 1:nrFluorColors
        
    %% get information on the case we're going to plot now
    fluorIndices = [fluorIdx fluorIdx];
    currentFluorLetters = [fluorValues(fluorIdx) fluorValues(fluorIdx)];

    % close figures from previous round if applicable
    if exist('hScatter','var'), if ishandle(hScatter), close(hScatter); end, end
    
    %% create figures
    hScatter=[];
    for groupIdx = 1:numel(applicableIndices)
        hScatter(groupIdx)=figure('Visible',FIGURESVISIBLE); clf; hold on; 
    end

    % get some colors
    somecolors = linspecer(numel(applicableIndices));

    %% go over groups
    meanLineHandles={};contourHandles={};averagedYHandles={};
    suggxlims1=[]; suggxlims2=[]; suggylims1=[]; suggylims2=[];
    for groupIdx = 1:numel(applicableIndices)

        %% plot all plots within that group (for this case)
        hCurrent = hScatter(groupIdx);

        meanLineHandles{groupIdx}=[];
        contourHandles{groupIdx}=[];
        averagedYHandles{groupIdx}=[];    
        for plotIdx = 1:numel(applicableIndices{groupIdx}) 


            %% 

            clear output

            %collectedOutput.(currentParameterNameToCollect) = {{}};

            % Set appropriate index
            dataIdx =  applicableIndices{groupIdx}(plotIdx);

            % Load configuration file and pre-process dataset info
            fluorDynamicsManager_sub_PreprocessDatasetInfo
                % also determines dataFileName

            load(dataFileName,'s_rm');

            % Get the production rate field
            rateFieldName = fieldNameDictPerGroup{groupIdx}{plotIdx}.paramNames.rate{fluorIdx};
            % Get the mu field
            muFieldName = fieldNameDictPerGroup{groupIdx}{plotIdx}.paramNames.muWithRate{fluorIdx};
            % Get the concentration field
            concFieldName = fieldNameDictPerGroup{groupIdx}{plotIdx}.paramNames.concentrationWithRate{fluorIdx};
            
            %% Now create a temporary field where prod. rate is divided by mu
            
            for schnitzIdx=1:numel(s_rm)
                s_rm(schnitzIdx).tempFieldProdDivMu = ...
                    s_rm(schnitzIdx).(rateFieldName) ./ ...
                    (log(2)/60*s_rm(schnitzIdx).(muFieldName));
            end
            
            %%                
            lineColor = somecolors(groupIdx,:);

            theRawfieldNames = {'tempFieldProdDivMu',...
                                concFieldName};

            % actually make a plot
            [meanLineHandle,contourHandle,individualLineHandles,averagedYHandle,suggxlim,suggylim] = fluorDynamicsManager_sub_scatterCloudsForGroup(hCurrent,s_rm,theRawfieldNames,fluorIndices,lineColor);
            meanLineHandles{groupIdx}(end+1)=meanLineHandle;
            contourHandles{groupIdx}(end+1)=contourHandle;
            if ~isempty(averagedYHandle)
                averagedYHandles{groupIdx}(end+1)=averagedYHandle;                
            end

            suggxlims1(end+1)=suggxlim(1);
            suggxlims2(end+1)=suggxlim(2);
            suggylims1(end+1)=suggylim(1);
            suggylims2(end+1)=suggylim(2);
            
        end

        %% Some cosmetics
        
        xlim([min(suggxlims1) max(suggxlims2)]);
        ylim([min(suggylims1) max(suggylims2)]);        
        
        uistack(contourHandles{groupIdx},'top');        
        uistack(meanLineHandles{groupIdx},'top');
        uistack(averagedYHandles{groupIdx},'top');

        set(averagedYHandle,'LineWidth',1);
        set(contourHandles{groupIdx},'LineWidth',1);

        % labels
        xlabel(['Production ' currentFluorLetters{1} ' / growth (a.u/area)']);
        ylabel(['Concentration ' currentFluorLetters{2} ' (a.u. / area)']);
        
        % fonts, size
        MW_makeplotlookbetter(10,[],[12.8/2, 19.2/3]/2,1);

        % Save the plots
        if ~exist('NOSAVEPLEASE')
            fileName = [GROUPNAME '_scatters_' HUMANREADABLENAMESFORGROUPSpathString{groupIdx} '_prod' currentFluorLetters{1} 'DivGrowth_vs_Concentration' currentFluorLetters{2} ];%_' set_fieldNames{setIdx}{1} '_' upper(fluorValues{fluorIdx})];
            saveas(hCurrent,[OUTPUTFOLDER 'tif_' fileName '.tif']);
            saveas(hCurrent,[OUTPUTFOLDER 'svg_' fileName '.svg']);
            saveas(hCurrent,[OUTPUTFOLDER 'fig_' fileName '.fig']);
            saveas(hCurrent,[OUTPUTFOLDER 'pdf_' fileName '.pdf']);
        end


    end
    
end

disp('Done making all special scatters..');



%% 

%{
for i= 1:numel(identifiers)
disp(identifiers{i});
end
%}

%%


%% 

disp('All done, hurray!!');



































