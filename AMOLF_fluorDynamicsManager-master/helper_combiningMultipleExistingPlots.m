


%% Putting existing figures into subplot

hNew = figure();
SIZE=[12.8,4];
set(hNew,'Units','centimeters','Position',[1,1,SIZE*2]);
s1 = subplot(1,2,1); %create and get handle to the subplot axes
MW_makeplotlookbetter(8*2);
s2 = subplot(1,2,2);
MW_makeplotlookbetter(8*2);

h1=openfig('\\storage01\data\AMOLF\users\wehrens\Latex3\Misc\testFigure.fig');
%figure(h1);
ax1=gca;
fig1 = get(ax1,'children'); %get handle to all the children in the figure
fig2 = get(ax1,'children');

copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
copyobj(fig2,s2);

% export the new figure
set(hNew,'RendererMode','manual','Renderer','Painters');
saveas(hNew, '\\storage01\data\AMOLF\users\wehrens\Latex3\Misc\testAutoCombinedFigs.fig');
saveas(hNew, '\\storage01\data\AMOLF\users\wehrens\Latex3\Misc\testAutoCombinedFigs.svg');
