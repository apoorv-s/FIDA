clear
clc

%%  Plot 1
load("Run2/ex3_200_01-02-23_artificialDynamics_jackKnife.mat");
load("Run2/ex3_200_01-02-23_auxilliary_artificialDynamics_jackKnife.mat", "xx", "yy");
load("Run2/analyticParticles.mat");

xSteps = xx(:, 1);
ySteps = yy(1, :);

figure
hold all

for i = 1:6
    subplot(2, 3, i)
    hold on
%     tString = strcat('$h(\mathbf{x}, t = ', string(tSteps(i)), ')$');

    uMat = reshape(analyticParticles(i, 1, :, :), size(xx));
    uColorPlt = imagesc(xx(:, 1), yy(1, :), uMat);

    if(i ~= 1)
        origin = scatter([0, 0], [0, 0], 'filled');
        rLine = line([0, obsShockPos(i)/sqrt(2)], [0, obsShockPos(i)/sqrt(2)]);
        rLine.LineWidth = 2.5;
        rLine.Color = 'red';
        rLine.LineStyle = '-';
        legend(rLine, 'Observed radius', Location='southeast');
        
    end

    cbar = colorbar;
%     cbar.Label.String = tString;
%     cbar.Label.Interpreter = 'latex';


%     hXLabel = xlabel('$x$', Interpreter='latex');
%     hYLabel = ylabel('$y$', Interpreter='latex');
    
%     set(gca, 'FontName', 'Helvetica')
%     set([hXLabel, hYLabel, cbar.Label], 'FontSize', 25)
    set(gca, 'xlim', [-2.5, 2.5], 'ylim', [-2.5, 2.5], 'FontSize', 16, 'LineWidth', 3)

end

%% Plot 2 Data

load("Run2/ex3_200_01-02-23_artificialDynamics_jackKnife.mat")
data{1} = permute(resampledParVecMat, [3, 2, 1]);
load("Run2/ex3_500_01-02-23_artificialDynamics_jackKnife.mat")
data{2} = permute(resampledParVecMat, [3, 2, 1]);
load("Run2/ex3_1000_02-02-23_artificialDynamics_jackKnife.mat")
data{3} = permute(resampledParVecMat, [3, 2, 1]);

nData = 3;
DataVals = [200, 500, 1000]; % length(parIndices) = nData
parIndices = [1, 2, 3, 4, 5, 6]; % length(parIndices) = nPars
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

%% plot 2

yLimits = {[18, 29], [0.5, 1.5], [0.2, 0.7], [-0.2, 0.2], [-0.2, 0.2], [6, 12]};
parName = {'$h_{in}$', '$h_{out}$', '$r$', '$x0$', '$y0$', '$g$'};
tIndex = nTimeSteps;
xLimits = [tSteps(1), tSteps(tIndex)];
parIdVarying = [1, 6];

figure
hold all
for id = 1:length(parIdVarying)
    parId = parIdVarying(id);
    yLimit = yLimits{parId};

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
    set(legend, 'color', 'white', 'box', 'on', 'edgecolor', 'white', 'FontSize', 16)
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
    set(legend, 'color', 'white', 'box', 'on', 'edgecolor', 'white', 'FontSize', 16)
%     text(0.02,0.98,'$nParticles = 200$','Units', 'Normalized', 'VerticalAlignment', 'Top', Interpreter='latex', FontSize=16)
    
    %
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
    set(legend, 'color', 'white', 'box', 'on', 'edgecolor', 'white', 'FontSize', 16)
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
