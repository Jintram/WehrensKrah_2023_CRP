function [meanLineHandle,contourHandle,individualLineHandle,averagedYHandle, suggxlim, suggylim] = fluorDynamicsManager_sub_scatterCloudsForGroup(h1,schnitzcells,theRawfieldNames,fluorIndices,lineColor)

%% Now do the same thing again but make scatter plots..
% -------------------------------------------------------------------------

    
NRCONTOURLINES=2;  
MARKERSIZE=2;
ALPHA=.25;

% Focus on given figure
set(0,'CurrentFigure',h1);
hold on;

%% Process & plot the data

% Get data
datax=[schnitzcells.(theRawfieldNames{1})];
datay=[schnitzcells.(theRawfieldNames{2})];

% Select data that is OK
if isfield(schnitzcells,'useForPlot')
    selectedIndices=find([schnitzcells.useForPlot]);
    disp([num2str(numel(selectedIndices)) ' will be selected out of ' num2str(numel(datax)) 'values.'])
    
    datax=datax(selectedIndices);
    datay=datay(selectedIndices);
    
end

% Plot individual points
%individualLineHandle = plot(datax,datay,'.','Color',lineColor);
individualLineHandle = scatter(datax,datay,MARKERSIZE^2,'MarkerEdgeColor',lineColor,'MarkerFaceColor',lineColor,'MarkerFaceAlpha',ALPHA,'MarkerEdgeAlpha',ALPHA);

% Calculate & plot contour (from kde)
[bandwidth,density,X,Y] = kde2d([datax;datay]');      
[C, contourHandle] = contour(X,Y,density,NRCONTOURLINES,'-k','LineWidth',2);

% Unfortunately contour data is hard to use to set limits because there can
% be some outliers
%{
Cx=C(1,:); Cy=C(2,:);
minCx=min(Cx); maxCx=max(Cx); minCy=min(Cy); maxCy=max(Cy);
dCx=maxCx-minCx; dCy=maxCy-minCy;
xlim([minCx-dCx/2 maxCx+dCx/2]);
ylim([minCy-dCy/2 maxCy+dCy/2]);
%}

% Get binned line
edges=linspace(min(datax),max(datax),100);
[meanValuesForBins, binCenters,stdValuesForBins,stdErrValuesForBins, counts,binnedValues,medianValuesForBins]=binnedaveraging({datax},{datay},edges)
selectionIdcs=counts>10;
averagedYHandle=plot(binCenters(selectionIdcs),meanValuesForBins(selectionIdcs),'-k');

% Calculate the mean
theScatterMean = [mean(datax(~isnan(datax))), mean(datay(~isnan(datay)))];  

% Plot mean
meanLineHandle = plot(theScatterMean(1),theScatterMean(2),'o');
set(meanLineHandle,'LineWidth',3,'MarkerFaceColor',lineColor,'MarkerEdgeColor','k','MarkerSize',5);

% Get the limits
sortedxdata=sort(datax);
sortedydata=sort(datay);
leftBoundx=sortedxdata(round(.1*numel(sortedxdata))); rightBoundx=sortedxdata(round(.9*numel(sortedxdata)));
leftBoundy=sortedydata(round(.1*numel(sortedydata))); rightBoundy=sortedydata(round(.9*numel(sortedydata)));
dx=rightBoundx-leftBoundx;
dy=rightBoundy-leftBoundy;
suggxlim=[leftBoundx-dx rightBoundx+dx];
suggylim=[leftBoundy-dy rightBoundy+dy];

disp('Plotting scatter done..');

    
end