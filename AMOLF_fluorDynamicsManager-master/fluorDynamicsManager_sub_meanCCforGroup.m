
function [meanLineHandle,individualLineHandles] = fluorDynamicsManager_sub_meanCCforGroup(h1,collectedCCs,fieldNameDict,fieldNames,fluorIndices,lineColor)
%%
% The idea is that this function receives a cell of cross-correlation
% functions that are determined, and then determines the average CC of that
% group and plots it, plus the different CCs that contribute to it in a
% lighter color.
% 
% Also required is a struct that allows this function to look up the
% parameternames per group.
%
% These functionalities are provided for by giving e.g. as input
%
% collectedOutput.CC{1}
% auxOutput.CC{1}
% fluorNrs=[1 1]
% fieldNames={'concentration', 'muWithConcentration'}
%
% Note that the cross-correlations that are asked for here should already
% be calculated beforehand (i.e. should exist).


%% Parameters

MARKER = '';
LINEWIDTH = 1;
FONTSIZE=15;
SHOWLEGEND=0;

% Focus on given figure
set(0,'CurrentFigure',h1);

        
%% We now collect the requested correlation functions into the
% parameters R and tau (we don't need the errorbars of the individual
% CCs)
collectedTau={}; collectedR={};
for ccIndex = 1:numel(collectedCCs)
    fieldName1 = fieldNameDict{ccIndex}.paramNames.(fieldNames{1}){fluorIndices(1)};
    fieldName2 = fieldNameDict{ccIndex}.paramNames.(fieldNames{2}){fluorIndices(2)};
    collectedTau{ccIndex} = collectedCCs{ccIndex}.([fieldName1 '_' fieldName2]).CorrData(:,1)';
    collectedR{ccIndex} = collectedCCs{ccIndex}.([fieldName1 '_' fieldName2]).CorrData(:,2)';
end
    

%% Plot separate ones
figure(h1); hold on;

mintaulist = []; maxtaulist = [];
individualLineHandles = [];
for i = 1:numel(collectedTau)
    % Plot
    %lighterColor = currentPlotColor+.5; lighterColor(lighterColor>1);
    lighterColor = (lineColor+[2,2,2])/3;
    l=plot(collectedTau{i},collectedR{i},['-' MARKER],'LineWidth',LINEWIDTH,'Color', lighterColor );% ,...
    if ~isempty(MARKER), set(l, 'MarkerSize', 10); end
                    %'MarkerEdgeColor', 'k');

    % Some stats for plot options
    mintaulist(end+1) = min(collectedTau{i});
    maxtaulist(end+1) = max(collectedTau{i});

    individualLineHandles(end+1) = l;
end


% Limits
if exist('USERXLIM','var')
    myxlimvalues = USERXLIM;
else
    myxlimvalues=[min(mintaulist), max(maxtaulist)];
end

xlim(myxlimvalues);    

if strcmp(fieldNames{1},fieldNames{2})
    % autocorrelation limits
    ylim([-.5 1]);
else
    % correlation limits
    ylim([-.5 .5]);
end

% Plot axes
plot(myxlimvalues,[0,0],'k-');
plot([0,0],[-1,1],'-k');

% Setting x-ticks
%{
if myxlimvalues(2)>6, tickDistance = 2; else, tickDistance=1; end
myxtickvalues=ceil(myxlimvalues(1)):tickDistance:floor(myxlimvalues(2));
set(gca,'XTick',myxtickvalues);       
%}

% Titles, fonts etc.
%xlabel('Delay (hrs)');
%ylabel(['Correlation'], 'Interpreter', 'none');
%Set all fontsizes

%% Determine binning
% ===   
% Previous determination of bins, was not optimal
% count taus for each dataset
tauCounts = cellfun(@(x) numel(x), collectedTau);
% determine which function has fewest datapoints
[nrDataPoints, indexLowestCount] = min(tauCounts);
% determine corresponding bin width using that one
binWidth = (max(collectedTau{indexLowestCount})-min(collectedTau{indexLowestCount}))/(nrDataPoints-1);
% determine window width from max and min values all taus
collectedTauPile = collectedTau{:};
windowBoundaries = [min(collectedTauPile)-binWidth*.5, max(collectedTauPile)+binWidth*.5];
windowWidth = windowBoundaries(2)-windowBoundaries(1);
% now create lin space
myBins = linspace(windowBoundaries(1),windowBoundaries(2),windowWidth/binWidth+1);
binCenters = myBins(2:end)-binWidth/2;    

% actual averaging
[meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins]=binnedaveraging(collectedTau,collectedR,myBins);

% plot meann line
meanLineHandle=plot(binCenters,meanValuesForBins,'k-','LineWidth',LINEWIDTH*2,'Color',lineColor);


%%
% Plot mean growth rate as bar
warning('TODO: also make growth rate bar');
%{
output.hrsPerDoublingMean = 1/mean([crossCorrData(idxsCrossCorrData).growthrateMean]);
plot([0,output.hrsPerDoublingMean],[.95,.95],'k-','LineWidth',10);
%}













