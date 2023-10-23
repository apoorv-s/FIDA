%%% Final Plots %%%

clear
clc
addpath('/Users/apoorv/Desktop/Git Repos/InviscidBurgersDA/MATLAB_ProfWeiKang/chebfun');

%% Plot 1 data

%%% Initial Condition

alpha = pi;
beta = 2;
gamma = pi;
icDomain = [-4, -1, 1, 4];
tmpVar = chebfun('tmpVar', icDomain);
g = chebfun({sin(alpha*tmpVar), beta, sin(gamma*tmpVar)}, ...
                    icDomain, 'splitting', 'on');
xStepsIC = linspace(-4, 4, 1000);
uVals = g(xStepsIC);
yLimitIC = [-1.2, 2.2];

%%% solution color plot
load("Run5/ex2_200_artDyn_jackKnife_aux.mat", 'observationMat', 'uAnalytic');
load("Run5/runData.mat");

xObs = [];
tObs = [];
for i = 1:length(observationMat)
    tempArr = cell2mat(observationMat{i});
    tempTArr = repmat(tSteps(i), 1, length(tempArr));
    xObs = [xObs, tempArr];
    tObs = [tObs, tempTArr];
end



%% Figure 1 tile

t = tiledlayout(1, 2);
% t.TileSpacing = 'compact';

nexttile
hold all
hFit = line(xStepsIC, uVals);
hData1 = line([-1, -1], [yLimitIC(1), 0]);
hData2 = line([1, 1], [yLimitIC(1), 0]);
hData3 = line([-4, -1], [beta, beta]);

set(hFit, 'Color', [0 0 .5], 'LineWidth', 2)
set(hData1, 'LineStyle', '--', 'Marker', '.', 'LineWidth', 2)
set(hData2, 'LineStyle', '--', 'Marker', '.', 'LineWidth', 2)
set(hData3, 'LineStyle', '--', 'Marker', '.', 'LineWidth', 2)

% hXLabel = xlabel('$x$', Interpreter='latex');
% hYLabel = ylabel('$u_{in}(x)$', Interpreter='latex');

% hText1 = text(-3.8, beta + 0.05, "$\beta = 2$", Interpreter="latex");


set(gca, 'FontName', 'Helvetica', 'ylim', yLimitIC, 'xlim', [-4, 4], 'FontSize', 16, 'LineWidth', 3)
% set([hXLabel, hYLabel, hText1], 'FontSize', 25)

nexttile
hold all
[tt, xx] = meshgrid(tSteps, xSteps);

uColorPlt = pcolor(xx, tt, uAnalytic);
uColorPlt.EdgeColor = 'none';
uColorPlt.FaceColor = 'interp';
uColorPlt.LineWidth = 0.00001; 

obsSkPosScatter = scatter(xObs, tObs, 50, 'filled', 'MarkerFaceColor','red');

% hXLabel = xlabel('$x$', Interpreter='latex');
% hYLabel = ylabel('$t$', Interpreter='latex');

legend(obsSkPosScatter, 'observed shock position')
set(legend, 'color', 'white', 'box', 'on', 'edgecolor', 'white')
cbar = colorbar;
% cbar.Label.String = '$u(x, t)$';
% cbar.Label.Interpreter = 'latex';

set(gca, 'xlim', [-4, 4], 'FontSize', 16, 'LineWidth', 6)
% set([cbar.Label, hXLabel, hYLabel], 'FontSize', 25)

%% Data for plot 2 and 3

load("Run5/runData.mat")
load("Run5/ex2_200_artDyn_jackKnife.mat")
% load("Run7/ex2_200_05_artDyn_jackKnife.mat")
data{1} = resampledParVecMat;
load("Run5/ex2_500_artDyn_jackKnife.mat")
% load("Run7/ex2_500_05_artDyn_jackKnife.mat")
data{2} = resampledParVecMat;
load("Run5/ex2_1000_artDyn_jackKnife.mat")
% load("Run7/ex2_1000_05_artDyn_jackKnife.mat")
data{3} = resampledParVecMat;

nData = 3;
DataVals = [200, 500, 1000]; % length(parIndices) = nData
nPars = 4;
parIndices = [1, 2, 3, 8]; % length(parIndices) = nPars
nTimeSteps = length(tSteps);
parOriginal = parametersOriginal(parIndices);
meanValues = zeros(nData, nPars, nTimeSteps);
stdValues = zeros(nData, nPars, nTimeSteps);
Quant95Values = zeros(nData, nPars, nTimeSteps);
Quant05Values = zeros(nData, nPars, nTimeSteps);
Quant95ValuesDiff = zeros(nData, nPars, nTimeSteps);
Quant05ValuesDiff = zeros(nData, nPars, nTimeSteps);
for i = 1:nData
    [meanValues(i, :, :), stdValues(i, :, :), ...
        Quant95Values(i, :, :), Quant05Values(i, :, :)] = getDataStats(data{i}, parIndices);
    Quant05ValuesDiff(i, :, :) = abs(Quant05Values(i, :, :) - meanValues(i, :, :));
    Quant95ValuesDiff(i, :, :) = abs(Quant95Values(i, :, :) - meanValues(i, :, :));
end
sumAllParsStd = zeros(nData, nTimeSteps);
for i = 1:nData
    sumAllParsStd(i, :) = sum(squeeze(stdValues(i, :, :)));
end
sumStdCutOff = 0.20000;
tCutOffId = find(sumAllParsStd(1, :) < sumStdCutOff, 1 );
% tCutOffId = 20;
tCutOffValue = tSteps(tCutOffId);

%% Plot 2

figure
hold all
yLimit = [0, 1];
hLine1 = line(tSteps, squeeze(stdValues(1, 1, :)), 'LineStyle', '-.', 'Color', 'r', 'LineWidth', 3);
hLine2 = line(tSteps, squeeze(stdValues(1, 2, :)), 'LineStyle', '--', 'Color', 'b', 'LineWidth', 3);
hLine3 = line(tSteps, squeeze(stdValues(1, 3, :)), 'LineStyle', '-.', 'Color', '#77AC30', 'LineWidth', 3);
hLine4 = line(tSteps, squeeze(stdValues(1, 4, :)), 'LineStyle', '--', 'Color', 'magenta', 'LineWidth', 3);
hLineSum = line(tSteps, sumAllParsStd(1, :), 'LineWidth', 3);

hLineCutOff = line([tSteps(1), tSteps(length(tSteps))], ...
    [sumStdCutOff, sumStdCutOff], "Color", 'black', 'LineStyle', ':');
hLineCutOff.LineWidth = 3;

hLineTimeCutOff = line([tCutOffValue, tCutOffValue], yLimit);
hLineTimeCutOff.LineStyle = ':';
hLineTimeCutOff.LineWidth = 3;
hLineTimeCutOff.Color = 'black';

% xLabel = xlabel("$t$", Interpreter="latex");
% yLabel = ylabel("$std(\sigma)$", Interpreter="latex");
ylim(yLimit)
lg = legend([hLine1, hLine2, hLine3, hLine4, hLineSum], "$\sigma_{\alpha}$", "$\sigma_{\beta}$", ...
        "$\sigma_{\gamma}$", "$\sigma_{\lambda}$", "$\sigma_{sum}$", 'Interpreter', 'latex');
set(lg, 'color', 'white', 'box', 'on', 'edgecolor', 'white', 'FontSize', 22)

% txt1 = text(0.02,0.98,strcat('$nParticles = $', string(DataVals(1))),'Units', 'Normalized', 'VerticalAlignment', 'Top',Interpreter='latex');
% txt2 = text(0.6,0.98,strcat('$t_{Off}$ = ', string(tCutOffValue)),'Units', 'Normalized', 'VerticalAlignment', 'Top', Interpreter='latex');

set(gca, 'FontSize', 16, 'LineWidth', 3)
% set([lg, txt1, txt2], 'FontSize', 25)
% set([xLabel, yLabel], 'FontSize', 25)

%% Plot 3
yLimits = {[3, 3.5], [1.8, 2.5], [3, 3.5], [1.6, 2.2]};
parName = {'$\alpha$', '$\beta$', '$\gamma$', '$\lambda$'};
tIndex = tCutOffId;
xLimits = [tSteps(1), tSteps(tIndex)];
parIdVarying = [1, 2, 3, 4];

figure
hold all
for id = 1:length(parIdVarying)
    parId = parIdVarying(id);
    yLimit = yLimits{parId};
    %
    dataId = 1;
    areaData = squeeze([Quant05Values(dataId, parId, 1:tIndex); ...
                   Quant95Values(dataId, parId, 1:tIndex) - Quant05Values(dataId, parId, 1:tIndex)])';
    
    splt1 = subplot(length(parIdVarying), 3, 3*(id - 1) + 1);
    hold on
    hLine1 = line(tSteps(1:tIndex), squeeze(meanValues(dataId, parId, 1:tIndex)));
    hLine1.LineWidth = 2;
    hLine2 = line([tSteps(1), tSteps(tIndex)], [parOriginal(parId), parOriginal(parId)], 'LineStyle', '-.', 'Color', 'r');
    hLine2.LineWidth = 1.5;
    hArea1 = area(tSteps(1:tIndex), areaData, BaseValue=yLimit(1), LineStyle=":");
    hArea1(1).FaceAlpha = 0;
    hArea1(2).FaceAlpha = 0.2;
    hArea1(2).FaceColor = [0.4660 0.6740 0.1880];
    
%     hXLabel = xlabel('$t$', Interpreter='latex');
%     hYLabel = ylabel(parName{parId}, Interpreter="latex");
    set(gca, 'FontName', 'Helvetica', 'FontSize', 16, 'ylim', yLimit, 'xlim', xLimits)
%     set([hXLabel, hYLabel], 'FontSize', 25)
    set(gca, 'Box', 'off', 'LineWidth', 2)
    legend([hLine1, hLine2, hArea1(2)], 'Mean value', 'True value', '5-95 percentile')
    set(legend, 'color', 'white', 'box', 'on', 'edgecolor', 'white', 'FontSize', 14)
%     text(0.02,0.98,'$nParticles = 200$','Units', 'Normalized', 'VerticalAlignment', 'Top', Interpreter='latex', FontSize=16)
    


    %
    dataId = 2;
    areaData = squeeze([Quant05Values(dataId, parId, 1:tIndex); ...
                   Quant95Values(dataId, parId, 1:tIndex) - Quant05Values(dataId, parId, 1:tIndex)])';
    
    splt2 = subplot(length(parIdVarying), 3, 3*(id - 1) + 2);
    hold on
    hLine1 = line(tSteps(1:tIndex), squeeze(meanValues(dataId, parId, 1:tIndex)));
    hLine1.LineWidth = 2;
    hLine2 = line([tSteps(1), tSteps(tIndex)], [parOriginal(parId), parOriginal(parId)], 'LineStyle', '-.', 'Color', 'r');
    hLine2.LineWidth = 1.5;
    hArea1 = area(tSteps(1:tIndex), areaData, BaseValue=yLimit(1), LineStyle=":");
    hArea1(1).FaceAlpha = 0;
    hArea1(2).FaceAlpha = 0.2;
    hArea1(2).FaceColor = [0.4660 0.6740 0.1880];
    
%     hXLabel = xlabel('$t$', Interpreter='latex');
%     hYLabel = ylabel(parName{parId}, Interpreter="latex");
    set(gca, 'FontName', 'Helvetica', 'FontSize', 16, 'ylim', yLimit, 'xlim', xLimits)
%     set([hXLabel, hYLabel], 'FontSize', 25)
    set(gca, 'Box', 'off', 'LineWidth', 2)
    legend([hLine1, hLine2, hArea1(2)], 'Mean value', 'True value', '5-95 percentile')
    set(legend, 'color', 'white', 'box', 'on', 'edgecolor', 'white', 'FontSize', 14)
%     text(0.02,0.98,'$nParticles = 200$','Units', 'Normalized', 'VerticalAlignment', 'Top', Interpreter='latex', FontSize=16)
     
    dataId = 3;
        areaData = squeeze([Quant05Values(dataId, parId, 1:tIndex); ...
                   Quant95Values(dataId, parId, 1:tIndex) - Quant05Values(dataId, parId, 1:tIndex)])';
    
    splt3 = subplot(length(parIdVarying), 3, 3*(id - 1) + 3);
    hold on
    hLine1 = line(tSteps(1:tIndex), squeeze(meanValues(dataId, parId, 1:tIndex)));
    hLine1.LineWidth = 2;
    hLine2 = line([tSteps(1), tSteps(tIndex)], [parOriginal(parId), parOriginal(parId)], 'LineStyle', '-.', 'Color', 'r');
    hLine2.LineWidth = 1.5;
    hArea1 = area(tSteps(1:tIndex), areaData, BaseValue=yLimit(1), LineStyle=":");
    hArea1(1).FaceAlpha = 0;
    hArea1(2).FaceAlpha = 0.2;
    hArea1(2).FaceColor = [0.4660 0.6740 0.1880];
    
%     hXLabel = xlabel('$t$', Interpreter='latex');
%     hYLabel = ylabel(parName{parId}, Interpreter="latex");
    set(gca, 'FontName', 'Helvetica', 'FontSize', 16, 'ylim', yLimit, 'xlim', xLimits)
%     set([hXLabel, hYLabel], 'FontSize', 25)
    set(gca, 'Box', 'off', 'LineWidth', 2)
    legend([hLine1, hLine2, hArea1(2)], 'Mean value', 'True value', '5-95 percentile')
    set(legend, 'color', 'white', 'box', 'on', 'edgecolor', 'white', 'FontSize', 14)
%     text(0.02,0.98,'$nParticles = 200$','Units', 'Normalized', 'VerticalAlignment', 'Top', Interpreter='latex', FontSize=16)

end


%% functions

function e = rmse(vec, mean)
    mean = double(mean);
    e = sqrt(sum(power(vec - mean, 2))/length(vec));
end

function [meanVal, stdVal, quant95ValDiff, quant05ValDiff] = getDataStats(data, parIndices)
    nPars = length(parIndices);
    nTimeSteps = size(data, 3);
    meanVal = zeros(nPars, nTimeSteps);
    stdVal = zeros(nPars, nTimeSteps);
    quant05ValDiff = zeros(nPars, nTimeSteps);
    quant95ValDiff = zeros(nPars, nTimeSteps);


    for i = 1:nPars
        meanVal(i, :) = mean(squeeze(data(parIndices(i), :, :)), 1);
        stdVal(i, :) = std(squeeze(data(parIndices(i), :, :)), 0, 1);
        quant05ValDiff(i, :) = quantile(squeeze(data(parIndices(i), :, :)), 0.05, 1);
        quant95ValDiff(i, :) = quantile(squeeze(data(parIndices(i), :, :)), 0.95, 1);
    end
end