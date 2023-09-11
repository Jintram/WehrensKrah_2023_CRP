clear; clc;


%% This notebook assumes that only the PRIME data is used (comparable microscope settings)
MUT_LOW = dlmread("MUT_LOW.dat"); MUT_LOW(MUT_LOW == 0) = NaN;
MUT_HIGH = dlmread("MUT_HIGH.dat"); MUT_HIGH(MUT_HIGH == 0) = NaN;
MUT_OPT = dlmread("MUT_OPT.dat"); MUT_OPT(MUT_OPT == 0) = NaN; 

WT = dlmread("WT.dat"); WT(WT == 0) = NaN;

%% determine mean per experiment
meanLOW = mean(MUT_LOW,2,'omitnan');
meanHIGH = mean(MUT_HIGH,2,'omitnan');
meanOPT = mean(MUT_OPT,2, 'omitnan');

meanWT = mean(WT,2,'omitnan');

%% determine slopes for optimal mutant and WT
dWT = WT./meanWT;
dMUT = MUT_OPT./meanOPT;

WTslopesy_mu = []; WTslopesh_mu = []; % initialize WT and MUT slope lists
MUTslopesy_mu = []; MUTslopesh_mu = [];

for iexp = 0:length(dWT(:,1))/3-1
    lmWT = fitlm(dWT(iexp*3+1,~isnan(dWT(iexp*3+1,:))),dWT(iexp*3+3,~isnan(dWT(iexp*3+3,:))));   
    WTslopesy_mu = [WTslopesy_mu, lmWT.Coefficients{2,1}];
    
    lmWT = fitlm(dWT(iexp*3+2,~isnan(dWT(iexp*3+2,:))),dWT(iexp*3+3,~isnan(dWT(iexp*3+3,:))));
    WTslopesh_mu = [WTslopesh_mu, lmWT.Coefficients{2,1}];
end


for iexp = 0:length(dMUT(:,1))/3-1
    lmWT = fitlm(dMUT(iexp*3+1,~isnan(dMUT(iexp*3+1,:))),dMUT(iexp*3+3,~isnan(dMUT(iexp*3+3,:))));   
    MUTslopesy_mu = [MUTslopesy_mu, lmWT.Coefficients{2,1}];
    
    lmWT = fitlm(dMUT(iexp*3+2,~isnan(dMUT(iexp*3+2,:))),dMUT(iexp*3+3,~isnan(dMUT(iexp*3+3,:))));
    MUTslopesh_mu = [MUTslopesh_mu, lmWT.Coefficients{2,1}];
end

x = [WTslopesy_mu,  WTslopesh_mu,  MUTslopesy_mu, MUTslopesh_mu];
x2 = {WTslopesy_mu, WTslopesh_mu,  MUTslopesy_mu, MUTslopesh_mu}; xCenter = 1:numel(x2); 

gWTY = repmat({"WT-Y"},6,1); gWT = repmat({"WT"},6,1);
gMUTY=repmat({"MUT-Y"},4,1); gMUT=repmat({"MUT"},4,1);
gWTH = repmat({"WT-H"},6,1);
gMUTH = repmat({"MUT-H"},4,1);
g = [gWTY; gWTH; gMUTY; gMUTH]; 


MUT_OPT = MUT_OPT(1:9,:); % Store only PRIME data sets (first 3) and continue only with this set
WT_prime  = WT(1:3,:);  % Store only PRIME data sets (last one) and continue.

stdLOW = 1.5*nanstd(reshape(MUT_LOW, 3,2*length(MUT_LOW')),0,2);
stdHIGH = 1.5*nanstd(reshape(MUT_HIGH, 3,2*length(MUT_HIGH')),0,2);
stdOPT = 1.5*nanstd(reshape(MUT_OPT,3,3*length(MUT_OPT')),0,2);

fullmeanLOW = nanmean(reshape(MUT_LOW, 3,2*length(MUT_LOW'))');
fullmeanHIGH = nanmean(reshape(MUT_HIGH, 3,2*length(MUT_HIGH'))');
fullmeanOPT = nanmean(reshape(MUT_OPT,3,3*length(MUT_OPT'))');

figure;
subplot(3,1,1);
boxplot(x,g); hold on;
for i = 1:numel(x2)
    plot(xCenter(i), x2{i}, 'ko','linewidth', 1); 
end
ylabel("Slope \phi - \lambda")


subplot(3,1,2);
boxplot([ meanWT(3:3:end); meanOPT(3:3:end)],[gWT;gMUT]);  hold on;
plot(1, meanWT(3:3:end), 'ko','linewidth', 1); hold on;
plot(2, meanOPT(3:3:end), 'ko','linewidth', 1); hold off;
ylim([0.4,0.9]);
ylabel("Mean growth rate per experiment");

subplot(3,1,3);
boxplot([ meanWT(1:3:end); meanWT(2:3:end);  meanOPT(1:3:end);  meanOPT(2:3:end)],g);  hold on;
plot(1, meanWT(1:3:end),'ko','linewidth', 1); hold on;
    plot(1, mean(WT_prime(1,:)), 'k*','linewidth',3); hold on;
plot(2, meanWT(2:3:end),'ko','linewidth', 1); hold on;
   plot(2, mean(WT_prime(2,:)), 'k*','linewidth',3); hold on;
plot(3, meanOPT(1:3:end),'ko','linewidth', 1); hold on;
    plot(3, meanOPT(1:3:9),'k*','linewidth', 3); hold on;
plot(4, meanOPT(2:3:end),'ko','linewidth', 1); hold on;
    plot(4, meanOPT(2:3:10),'k*','linewidth', 3); hold off;
ylabel("Mean concentrations"); ylim([0, 500]);

%% Plotting data sets

% plot colors
RGBWTY = [30,144,255]/256;
RGBWTH = [135,206,250]/256;

RGBLOWY = [255,0,0]/256;
RGBLOWH = [255,99,71]/256;

RGBMUTY = [60,179,113]/256 ;
RGBMUTH = [144,238,144]/256;

RGBHIGHY = [255,140,0]/256;
RGBHIGHH = [255,165,0]/256;

% These values (a = -4.38*10^(-7); b = 2.55*10^(-4); c = -0.0276;) were raported in Martijns thesis, but here we fit new
% values:
f = figure;

scatter(MUT_LOW(1,:),MUT_LOW(3,:),20,RGBLOWY,'.','MarkerEdgeAlpha',.3); hold on;  ylim([0,1.5]); xlim([0, 500]);
scatter(MUT_OPT(1,:),MUT_OPT(3,:),20,RGBMUTY,'.','MarkerEdgeAlpha',.3); hold on; 
scatter(MUT_HIGH(1,:),MUT_HIGH(3,:),20,RGBHIGHY,'.','MarkerEdgeAlpha',.3); hold on; 
scatter(MUT_LOW(4,:),MUT_LOW(6,:),20,RGBLOWY,'.','MarkerEdgeAlpha',.3); hold on; 
scatter(MUT_OPT(4,:),MUT_OPT(6,:),20,RGBMUTY,'.','MarkerEdgeAlpha',.3); hold on; 
scatter(MUT_OPT(7,:),MUT_OPT(9,:),20,RGBMUTY,'.','MarkerEdgeAlpha',.3); hold on; 
scatter(MUT_HIGH(1,:),MUT_HIGH(3,:),20,RGBHIGHY,'.','MarkerEdgeAlpha',.3); hold on; 
scatter(MUT_HIGH(4,:),MUT_HIGH(6,:),20,RGBHIGHY,'.','MarkerEdgeAlpha',.3); hold on;
pbaspect([0.9 1 1]); 


LOWfit = reshape(MUT_LOW, 3, 2*length(MUT_LOW));
lowfit = fitlm(LOWfit(1,~isnan(LOWfit(1,:))),LOWfit(3,~isnan(LOWfit(3,:))));
fplot(@(x) lowfit.Coefficients{1,1} +lowfit.Coefficients{2,1}*x, [fullmeanLOW(1) - stdLOW(1),fullmeanLOW(1) + stdLOW(1)],':k', 'Linewidth',2); hold on;

OPTfit = reshape(MUT_OPT, 3, 3*length(MUT_OPT));
optfit = fitlm(OPTfit(1,~isnan(OPTfit(1,:))),OPTfit(3,~isnan(OPTfit(3,:))));
fplot(@(x) optfit.Coefficients{1,1} +optfit.Coefficients{2,1}*x, [fullmeanOPT(1) - stdOPT(1),fullmeanOPT(1) + stdOPT(1)],':k', 'Linewidth',2); hold on;

HIGHfit = reshape(MUT_HIGH, 3, 2*length(MUT_HIGH));
highfit = fitlm(HIGHfit(1,~isnan(HIGHfit(1,:))),HIGHfit(3,~isnan(HIGHfit(3,:))));
fplot(@(x) highfit.Coefficients{1,1} + highfit.Coefficients{2,1}*x, [fullmeanHIGH(1) - stdHIGH(1),fullmeanHIGH(1) + stdHIGH(1)],':k', 'Linewidth',2); hold on;


plot([fullmeanLOW(1), fullmeanOPT(1), fullmeanHIGH(1)], [fullmeanLOW(3),fullmeanOPT(3),fullmeanHIGH(3)],'.k','Markersize',15); hold on;
pY = polyfit([fullmeanLOW(1), fullmeanOPT(1), fullmeanHIGH(1)],[fullmeanLOW(3),fullmeanOPT(3),fullmeanHIGH(3)],2);
fplot(@(x) pY(1).*x.^2 + pY(2).*x + pY(3), [150,390],'-k');
%fplot(@(x) 60/log(2)*(a.*x.^2 + b.*x + c), [150 420],'-k'); hold off;  ylim([0,1.5]); xlim([0, 500]);
%legend("Low cAMP", "Optimal cAMP","High cAMP",'Location','northwest', 'FontSize',12); 
pbaspect([0.9 1 1]);
box on; xlabel("\phi_Y (a.u.)"); ylabel("\lambda (h^{-1})");

clear LOWFIT lowfit OPTfit optfit HIGHfit highfit
figure; ylim([0,1.5]); xlim([0, 500]);
scatter(MUT_LOW(2,:),MUT_LOW(3,:),20,RGBLOWH,'.','MarkerEdgeAlpha',.3); hold on; ylim([0,1.5]); xlim([0, 500]);
scatter(MUT_LOW(5,:),MUT_LOW(6,:),20,RGBLOWH,'.','MarkerEdgeAlpha',.3); hold on; 
scatter(MUT_OPT(2,:),MUT_OPT(3,:),20,RGBMUTH,'.','MarkerEdgeAlpha',.3); hold on; 
scatter(MUT_OPT(5,:),MUT_OPT(6,:),20,RGBMUTH,'.','MarkerEdgeAlpha',.3); hold on; 
scatter(MUT_OPT(8,:),MUT_OPT(9,:),20,RGBMUTH,'.','MarkerEdgeAlpha',.3); hold on; 
scatter(MUT_HIGH(2,:),MUT_HIGH(3,:),20,RGBHIGHH,'.','MarkerEdgeAlpha',.3); hold on; 
scatter(MUT_HIGH(5,:),MUT_HIGH(6,:),20,RGBHIGHH,'.','MarkerEdgeAlpha',.3); hold on;
pbaspect([0.9 1 1]); 

LOWfit = reshape(MUT_LOW, 3, 2*length(MUT_LOW));
lowfit = fitlm(LOWfit(2,~isnan(LOWfit(1,:))),LOWfit(3,~isnan(LOWfit(3,:))));
fplot(@(x) lowfit.Coefficients{1,1} +lowfit.Coefficients{2,1}*x, [fullmeanLOW(2) - stdLOW(2),fullmeanLOW(2)+ stdLOW(2)],':', 'color', [0.5,0.5,0.5], 'Linewidth',2); hold on;

OPTfit = reshape(MUT_OPT, 3, 3*length(MUT_OPT));
optfit = fitlm(OPTfit(2,~isnan(OPTfit(1,:))),OPTfit(3,~isnan(OPTfit(3,:))));
fplot(@(x) optfit.Coefficients{1,1} +optfit.Coefficients{2,1}*x, [fullmeanOPT(2) - stdOPT(2),fullmeanOPT(2) + stdOPT(2)],':', 'color', [0.5,0.5,0.5], 'Linewidth',2); hold on;

HIGHfit = reshape(MUT_HIGH, 3, 2*length(MUT_HIGH));
highfit = fitlm(HIGHfit(2,~isnan(HIGHfit(1,:))),HIGHfit(3,~isnan(HIGHfit(3,:))));
fplot(@(x) highfit.Coefficients{1,1} + highfit.Coefficients{2,1}*x, [fullmeanHIGH(2) - stdHIGH(2),fullmeanHIGH(2) + stdHIGH(2)],':', 'color', [0.5,0.5,0.5], 'Linewidth',2); hold on;

plot([fullmeanLOW(2), fullmeanOPT(2), fullmeanHIGH(2)], [fullmeanLOW(3),fullmeanOPT(3),fullmeanHIGH(3)],'.','Markersize',15,'color', [0.5,0.5,0.5]); hold on;
ylim([0,1.5]); xlim([0, 500]);
fplot(@(x) pY(1).*(445-x).^2 + pY(2).*(445-x) + pY(3), [50 300],'-', 'color', [0.5,0.5,0.5]); hold off;
pbaspect([0.9 1 1]);
box on; xlabel("\phi_H (a.u.)"); ylabel("\lambda (h^{-1})");


% Significance test (Welch's t-test):
% slopes H:
x = WTslopesh_mu; y = MUTslopesh_mu;
%t= (mean(x) - mean(y))/sqrt(var(x)/length(x)+ var(y)/length(y))
%df = (var(x)/length(x)+ var(y)/length(y))^2/(var(x)^2/(length(x)^2*(length(x)-1))+ var(y)^2/(length(y)^2*(length(y)-1)))
%p_slopeH = 2*(1-tcdf(abs(t), floor(df)));
[h,p_slopeH] = ttest2(x,y, 'Vartype','unequal');

% growth rates:
x = meanWT(3:3:end); y= meanOPT(3:3:end);
%t = (mean(x) - mean(y))/sqrt(var(x)/length(x)+ var(y)/length(y))
%df = (var(x)/length(x)+ var(y)/length(y))^2/(var(x)^2/(length(x)^2*(length(x)-1))+ var(y)^2/(length(y)^2*(length(y)-1)))
%p_mus = 2*(1-tcdf(abs(t), floor(df)));

[h,p_mus] = ttest2(x,y, 'Vartype','unequal');
% --> look up t statistics in t table.
