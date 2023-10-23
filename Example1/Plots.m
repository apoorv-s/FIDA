%%% Final Plots %%%
clear
clc

%% Plot 1 Data

%%% Initial Condition
u_l = 2;
u_r = 1;
x_r = 1;
lambda = 2;
beta = 2/lambda;
alpha = (u_l - u_r)/x_r;
xStepsIC = linspace(-1, 5, 1000);
uValsIC = zeros(size(xStepsIC));
for n = 1:length(xStepsIC)
    if(xStepsIC(n) < 0)
        uValsIC(n) = u_l;
    elseif(0 <= xStepsIC(n) && xStepsIC(n) <= x_r)
        uValsIC(n) = (u_l - alpha*xStepsIC(n));
    else
        uValsIC(n) = u_r;
    end
end

%%% solution color plot
load('Run7/ex1_200_06-02-23_artificialDynamics_jackKnife.mat');
tStepsObs = tSteps;
load('analyticSolution.mat');

%% Figure 1 tile

t = tiledlayout(1, 2);
% t.TileSpacing = 'compact';

nexttile
hold all
hFit = line(xStepsIC, uValsIC);
hData1 = line([0, 0], [0, u_l]);
hData2 = line([x_r, x_r], [0, u_r]);

set(hFit, 'Color', [0 0 .5], 'LineWidth', 2)
set(hData1, 'LineStyle', '--', 'Marker', '.', 'LineWidth', 3)
set(hData2, 'LineStyle', '--', 'Marker', '.', 'LineWidth', 3)

% hXLabel = xlabel('$x$', Interpreter='latex');
% hYLabel = ylabel('$u_{in}(x)$', Interpreter='latex');
% 
% hText = text(0, u_l + 0.05, "$u_l = 2$", Interpreter="latex");
% hText2 = text(1, u_r + 0.05, "$u_r = 1$", Interpreter="latex");
% hText3 = text(1,  0.85, "$x_r = 1$", Interpreter="latex");

set(gca, 'FontName', 'Helvetica', 'ylim', [0.8, 2.2], 'xlim', [-1, 4], 'FontSize', 16, 'linewidth', 3);
% set([hXLabel, hYLabel, hText, hText2, hText3], 'FontSize', 20)



nexttile
hold all
[xx, tt] = meshgrid(xSteps, tSteps);

uColorPlt = pcolor(xx, tt, uAnalytic);
uColorPlt.EdgeColor = 'none';
uColorPlt.FaceColor = 'interp';

obsSkPosScatter = scatter(obsShockPos(1, 2:end), tStepsObs(1, 2:end), 50, 'filled', 'MarkerFaceColor','red');
% 
% hYLabelColorPlot = ylabel('$t$', Interpreter='latex');
% hXLabelColorPlot = xlabel('$x$', Interpreter='latex');

legend(obsSkPosScatter, 'observed shock position')
set(legend, 'color', 'white', 'box', 'on', 'edgecolor', 'white')
cbar = colorbar;
% cbar.Label.String = 'u(x, t)';
% cbar.Label.Interpreter = 'latex';

set(gca, 'xlim', [-1, 8], 'FontSize', 16, 'linewidth', 6)
% set([hXLabelColorPlot, hYLabelColorPlot], 'FontSize', 20)


%% Plot 2 ans 3 Data
clc

load("Run9/ex1_200_14-06-23_artificialDynamics_jackKnife.mat")
data{1} = permute(resampledParVecMat,[3 2 1]);
load("Run9/ex1_500_14-06-23_artificialDynamics_jackKnife.mat")
data{2} = permute(resampledParVecMat,[3 2 1]);
load("Run9/ex1_1000_14-06-23_artificialDynamics_jackKnife.mat")
data{3} = permute(resampledParVecMat,[3 2 1]);

nData = 3;
DataVals = [200, 500, 1000]; % length(parIndices) = nData
parIndices = [1, 2, 3, 4]; % length(parIndices) = nPars
nTimeSteps = size(data{1}, 3);
parOriginal = parametersOriginal(parIndices);
meanValues = zeros(nData, nParameters, nTimeSteps);
stdValues = zeros(nData, nParameters, nTimeSteps);
Quant95Values = zeros(nData, nParameters, nTimeSteps);
Quant05Values = zeros(nData, nParameters, nTimeSteps);
Quant95ValuesDiff = zeros(nData, nParameters, nTimeSteps);
Quant05ValuesDiff = zeros(nData, nParameters, nTimeSteps);

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
sumStdCutOff = 0.1;
tCutOffId = find(sumAllParsStd(1, :) < sumStdCutOff, 1 );
tCutOffValue = tSteps(tCutOffId);

%% Plot 2

figure
hold all
yLimit = [0, 0.3];
hLine1 = line(tSteps, squeeze(stdValues(1, 1, :)), "LineStyle", "-.", "Color", 'r', "LineWidth", 3);
hLine2 = line(tSteps, squeeze(stdValues(1, 4, :)), "LineStyle", "--", "Color", '#77AC30', "LineWidth", 3);
hLine3 = line([tCutOffValue, tCutOffValue], yLimit);
hLineSum = line(tSteps, sumAllParsStd(1, :), "LineWidth", 2);

hLineCutOff = line([tSteps(1), tSteps(length(tSteps))], ...
    [sumStdCutOff, sumStdCutOff], "Color", 'black', 'LineStyle', ':');
hLineCutOff.LineWidth = 3;

hLine3.LineStyle = ':';
hLine3.LineWidth = 3;


% xLabel = xlabel("$t$", Interpreter="latex");
% yLabel = ylabel("$std(\sigma)$", Interpreter="latex");
ylim(yLimit)

lg = legend([hLine1, hLine2, hLineSum], "$\sigma_{\lambda}$", "$\sigma_{x_r}$", "$\sigma_{sum} = \sigma_{\lambda} + \sigma_{x_r}$", 'Interpreter', 'latex');
set(lg, 'color', 'white', 'box', 'on', 'edgecolor', 'white', 'FontSize', 22)
% txt1 = text(0.02,0.98,strcat('nParticles = ', string(DataVals(1))),'Units', 'Normalized', 'VerticalAlignment', 'Top', Interpreter='latex');
% txt2 = text(0.8,0.4,strcat('$t_{Off}$ = ', string(tCutOffValue)),'Units', 'Normalized', 'VerticalAlignment', 'Top', Interpreter='latex');

set(gca, 'FontSize', 16, 'linewidth', 3)
% set([xLabel, yLabel], 'FontSize', 25)

%% UKF Data
% load("UKFResults2.mat")
load("ex1UKFFromPy.mat")
estimatedParsUKF = estimatedState';
estimatedParsStdUKF = estimatedStateStd';
%% Plot Pars

yLimits = {[1.9, 2.4], [1.5, 2.5], [0.5, 1.5], [0.6, 1.25]};
% yLimitsUKF = {[0.5, 3.8], [0.5, 1.8]};
yLimitsUKF = {[1.9, 2.4], [0.6, 1.25]};
parName = {'$\lambda$', '$u_l$', '$u_r$', '$x_r$'};
tIndex = tCutOffId;
xLimits = [tSteps(1), tSteps(tIndex)];
parIdVarying = [1, 4];

zScore = 1.644854;

%% Plot 3

ukfPlotCounter = 1;
figure
hold all
for id = 1:length(parIdVarying)
    parId = parIdVarying(id);
    yLimit = yLimits{parId};
    yLimitUKF = yLimitsUKF{ukfPlotCounter};
    %
    dataId = 1;
    areaData = squeeze([Quant05Values(dataId, parId, 1:tIndex); ...
                   Quant95Values(dataId, parId, 1:tIndex) - Quant05Values(dataId, parId, 1:tIndex)])';
    
    splt1 = subplot(length(parIdVarying), 4, 4*(id - 1) + 1);
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
    set(legend, 'color', 'white', 'box', 'on', 'edgecolor', 'white', 'FontSize', 16)
%     text(0.02,0.98,'$nParticles = 200$','Units', 'Normalized', 'VerticalAlignment', 'Top', Interpreter='latex', FontSize=16)
    
    
    %
    dataId = 2;
    areaData = squeeze([Quant05Values(dataId, parId, 1:tIndex); ...
                   Quant95Values(dataId, parId, 1:tIndex) - Quant05Values(dataId, parId, 1:tIndex)])';
    
    splt2 = subplot(length(parIdVarying), 4, 4*(id - 1) + 2);
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
    set(legend, 'color', 'white', 'box', 'on', 'edgecolor', 'white', 'FontSize', 16)
%     text(0.02,0.98,'$nParticles = 200$','Units', 'Normalized', 'VerticalAlignment', 'Top', Interpreter='latex', FontSize=16)
    
    %
    dataId = 3;
    areaData = squeeze([Quant05Values(dataId, parId, 1:tIndex); ...
                   Quant95Values(dataId, parId, 1:tIndex) - Quant05Values(dataId, parId, 1:tIndex)])';
    
    splt3 = subplot(length(parIdVarying), 4, 4*(id - 1) + 3);
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
    set(legend, 'color', 'white', 'box', 'on', 'edgecolor', 'white', 'FontSize', 16)
%     text(0.02,0.98,'$nParticles = 200$','Units', 'Normalized', 'VerticalAlignment', 'Top', Interpreter='latex', FontSize=16)

%UKF Plot
    areaData = squeeze([(estimatedParsUKF(1:tIndex, ukfPlotCounter) - zScore*estimatedParsStdUKF(1:tIndex, ukfPlotCounter))' ; ...
                   (2*zScore*estimatedParsStdUKF(1:tIndex, ukfPlotCounter))'])';
    
    splt4 = subplot(length(parIdVarying), 4, 4*(id - 1) + 4);
    hold on
    hLine1 = line(tSteps(1:tIndex), squeeze(estimatedParsUKF(1:tIndex, ukfPlotCounter)));
    hLine1.LineWidth = 2;
    hLine2 = line([tSteps(1), tSteps(tIndex)], [parOriginal(parId), parOriginal(parId)], 'LineStyle', '-.', 'Color', 'r');
    hLine2.LineWidth = 1.5;
    hArea1 = area(tSteps(1:tIndex), areaData, BaseValue=yLimitUKF(ukfPlotCounter), LineStyle=":");
    hArea1(1).FaceAlpha = 0;
    hArea1(2).FaceAlpha = 0.2;
    hArea1(2).FaceColor = [0.4660 0.6740 0.1880];
%     hXLabel = xlabel('$t$', Interpreter='latex');
%     hYLabel = ylabel(parName{parId}, Interpreter="latex");

    set(gca, 'FontName', 'Helvetica', 'FontSize', 16, 'ylim', yLimitUKF, 'xlim', xLimits)
%     set([hXLabel, hYLabel], 'FontSize', 25)
    set(gca, 'Box', 'off', 'LineWidth', 2)
    legend([hLine1, hLine2, hArea1(2)], 'Mean value', 'True value', '5-95 percentile')
    set(legend, 'color', 'white', 'box', 'on', 'edgecolor', 'white', 'FontSize', 16)
%     text(0.02,0.98,'$nParticles = 200$','Units', 'Normalized', 'VerticalAlignment', 'Top', Interpreter='latex', FontSize=16)
    ukfPlotCounter = ukfPlotCounter + 1;

end

%%

% yLimitsUKF = {[0.5, 3.8], [0.5, 1.8]};
yLimitsUKF = {[1.9, 2.4], [0.6, 1.25]};

ukfPlotCounter = 1;
figure
hold all
yLimitUKF = yLimitsUKF{ukfPlotCounter};
areaData = squeeze([(estimatedParsUKF(1:tIndex, ukfPlotCounter) - estimatedParsStdUKF(1:tIndex, ukfPlotCounter))' ; ...
               (2*estimatedParsStdUKF(1:tIndex, ukfPlotCounter))'])';

splt1 = subplot(1, 2, 1);
hold on
hLine1 = line(tSteps(1:tIndex), squeeze(estimatedParsUKF(1:tIndex, ukfPlotCounter)));
hLine1.LineWidth = 2;
hLine2 = line([tSteps(1), tSteps(tIndex)], [parOriginal(1), parOriginal(1)], 'LineStyle', '-.', 'Color', 'r');
hLine2.LineWidth = 1.5;
hArea1 = area(tSteps(1:tIndex), areaData, BaseValue=yLimitUKF(ukfPlotCounter), LineStyle=":");
hArea1(1).FaceAlpha = 0;
hArea1(2).FaceAlpha = 0.2;
hArea1(2).FaceColor = [0.4660 0.6740 0.1880];
%     hXLabel = xlabel('$t$', Interpreter='latex');
%     hYLabel = ylabel(parName{parId}, Interpreter="latex");

set(gca, 'FontName', 'Helvetica', 'FontSize', 16, 'ylim', yLimitUKF, 'xlim', xLimits)
%     set([hXLabel, hYLabel], 'FontSize', 25)
set(gca, 'Box', 'off', 'LineWidth', 2)
legend([hLine1, hLine2, hArea1(2)], 'Mean value', 'True value', '2 std')
set(legend, 'color', 'white', 'box', 'on', 'edgecolor', 'white', 'FontSize', 16)
%     text(0.02,0.98,'$nParticles = 200$','Units', 'Normalized', 'VerticalAlignment', 'Top', Interpreter='latex', FontSize=16)
ukfPlotCounter = ukfPlotCounter + 1;

yLimitUKF = yLimitsUKF{ukfPlotCounter};
areaData = squeeze([(estimatedParsUKF(1:tIndex, ukfPlotCounter) - estimatedParsStdUKF(1:tIndex, ukfPlotCounter))' ; ...
               (2*estimatedParsStdUKF(1:tIndex, ukfPlotCounter))'])';

splt2 = subplot(1, 2, 2);
hold on
hLine1 = line(tSteps(1:tIndex), squeeze(estimatedParsUKF(1:tIndex, ukfPlotCounter)));
hLine1.LineWidth = 2;
hLine2 = line([tSteps(1), tSteps(tIndex)], [parOriginal(4), parOriginal(4)], 'LineStyle', '-.', 'Color', 'r');
hLine2.LineWidth = 1.5;
hArea1 = area(tSteps(1:tIndex), areaData, BaseValue=yLimitUKF(ukfPlotCounter), LineStyle=":");
hArea1(1).FaceAlpha = 0;
hArea1(2).FaceAlpha = 0.2;
hArea1(2).FaceColor = [0.4660 0.6740 0.1880];
%     hXLabel = xlabel('$t$', Interpreter='latex');
%     hYLabel = ylabel(parName{parId}, Interpreter="latex");

set(gca, 'FontName', 'Helvetica', 'FontSize', 16, 'ylim', yLimitUKF, 'xlim', xLimits)
%     set([hXLabel, hYLabel], 'FontSize', 25)
set(gca, 'Box', 'off', 'LineWidth', 2)
legend([hLine1, hLine2, hArea1(2)], 'Mean value', 'True value', 'mean \pm \sigma')
set(legend, 'color', 'white', 'box', 'on', 'edgecolor', 'white', 'FontSize', 16)
%     text(0.02,0.98,'$nParticles = 200$','Units', 'Normalized', 'VerticalAlignment', 'Top', Interpreter='latex', FontSize=16)

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