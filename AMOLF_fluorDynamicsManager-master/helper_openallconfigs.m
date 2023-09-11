
% Helper that opens all of the config files
%
% Run loading excel section and load appropriate dataset parameters first.

for dataIdx= 1:nrDataLines
    
    %% skip if empty    
    if isempty(allXLSdata{dataIdx,1}) | isnan(allXLSdata{dataIdx,1})
        continue
    end
    
    %% Load configuration file and pre-process dataset info
    fluorDynamicsManager_sub_PreprocessDatasetInfo      
    
    winopen([dateDir configFileName]);
    
end

disp('All config files are now open');