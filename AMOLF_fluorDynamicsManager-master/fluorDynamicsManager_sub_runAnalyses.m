%% Running analyses if necessary
%
% Note: if the script crashes (e.g. a problem with setting axis limits)
% this might be due to NaN values in the data (check e.g. CorrData). This
% can be resolved by using MW_helper_schnitzcell_terror_counting to
% identify cells that have NaN values. After checking these errors are
% non-consequential for the data set (e.g. because it are the last
% schnitzes in the set that only have 1 frame and no offspring), they can
% simply be added to the "bad schnitzes" list in the excel config file of
% that data set. Open MW_GUI_schnitzcells to easily edit the config file.
%
% OPTIONAL PARAMETERS TO SET:
% - NOSELECTIONLOOKATALLANALYSES=1
% - RERUNANALYSIS=1

% Run over datafiles to run analysis if this is necessary
for dataIdx= 1:nrDataLines
    
    %% skip if empty    
    if isempty(allXLSdata{dataIdx,1}) | isnan(allXLSdata{dataIdx,1})
        continue
    end
    
    %% Load configuration file and pre-process dataset info
    fluorDynamicsManager_sub_PreprocessDatasetInfo        
        % some parameters that are determined are:
        % ourSettings, dateDir, theMovieDateOnly, movieName, dataFileName,
        % identifierInXLS
    
    %% Now the idea is to check for the existence of already analyzed data    
    
    % Now determine if we want to re-run the analysis
    % (I.e. if we want to go over all, or if dataset is mentioned in any of 
    % the to-plot identifiers.)
    if exist('NOSELECTIONLOOKATALLANALYSES','var') | any(cellfun(@(x) strcmp(identifierInXLS,x), [IDENTIFIERSTOPLOT{:}]))
    
        % And run the analysis if there is not already data
        if exist(dataFileName,'file') & ~exist('RERUNANALYSIS','var')

            disp(['Dataset [' theMovieDateOnly ', ' movieName '] was already analyzed before; I did not repeat analysis.']);

        else

            %% Otherwise run the analysis
            disp(['Dataset [' theMovieDateOnly ', ' movieName '] was not analyzed before (or user defined re-run); running analysis now.']);

            disp('Taking 5 seconds pause before starting..');
            pause(5);

            % ===
            % Update the production field names
            warning('Overwriting production field name');
            ourSettings.fluorDerivativeFieldName = 'dX5_divAreaPx_cycCor'; % this is the updated rate field            
            % overwrite the PLOTSCATTER setting
            ourSettings.PLOTSCATTER = 0;
            % Don't show figures
            FIGUREVISIBLE='off';      
            % Re-evaluate the bad schnitzes (in case the Excel file was
            % updated)
            ourSettings.alreadyRemovedInMatFile = 0;
            % tell master script schnitzcells data should be loaded from directory
            LOADDATASETATEND = 1;
            % Note:
            % also ourSettings.fitTimeCrosscorr is an important parameter
            % that decides which part of the branchData is taken.

            % run the masterscript analysis part
            runsections='makeoutputfull'; 
            Schnitzcells_masterscript

            % A number of (invisible) figures will have been
            % generated, so close all figures here to avoid accumulation..
            close all;
            
            % The masterscript will save the data
            disp('Analysis done and saved.');
        end
       
    end
        
end
disp('Section done, all datasets were inspected to see if analyis was necessary..');