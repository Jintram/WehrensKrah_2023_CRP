%% Open applicable excel files

%%
for groupIdx = 1:numel(applicableIndices)
       
    
    %%
    for plotIdx = 1:numel(applicableIndices{groupIdx}) 


        %% 

        %clear output
        
        %collectedOutput.(currentParameterNameToCollect) = {{}};

        % Set appropriate index
        dataIdx =  applicableIndices{groupIdx}(plotIdx);

        % Load configuration file and pre-process dataset info
        fluorDynamicsManager_sub_PreprocessDatasetInfo
            % also determines dataFileName
            
        winopen([dateDir configFileName]);
            
    end
    
end