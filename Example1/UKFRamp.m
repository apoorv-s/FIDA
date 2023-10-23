clear
clc

global time
time = 0;

addpath("MatlabUKF/")

load("Run9/ex1_200_14-06-23_artificialDynamics_jackKnife.mat", "tSteps", "obsShockPos", ...
    "artificialDynamicsCov", "parameterOrder", "parametersOriginal", "parametersRange", "obsNoiseStd", "")

uncertainParIndices = [1, 4];

trueState = parametersOriginal(uncertainParIndices);
parameterRangesMeanValue = mean(parametersRange, 2);
%%

% Your initial state guess at time k, utilizing measurements up to time k-1: xhat[k|k-1]
initialStateGuess = [parameterRangesMeanValue(uncertainParIndices(1));
                       parameterRangesMeanValue(uncertainParIndices(2))]; % xhat[k|k-1]
% Construct the filter
ukf = unscentedKalmanFilter(...
    @ForwardModel,... % State transition function
    @ObsOperator,... % Measurement function
    initialStateGuess,...
    'HasAdditiveMeasurementNoise',true);

%%
R = obsNoiseStd^2; % Variance of the measurement noise v[k]
ukf.MeasurementNoise = R;


% T = 0.05; % [s] Filter sample time
% tSteps = 0:T:5;

% rng(1); % Fix the random number generator for reproducible results


Nsteps = numel(tSteps); % Number of time steps

% yTrue = zeros(Nsteps,1);
% yMeas = zeros(Nsteps,1); % sqrt(R): Standard deviation of noise

yMeas = obsShockPos;

xCorrectedUKF = zeros(Nsteps,2); % Corrected state estimates
PCorrected = zeros(Nsteps,2,2); % Corrected state estimation error covariances
e = zeros(Nsteps,1); % Residuals (or innovations)

xCorrectedUKF(1, :) = initialStateGuess;
PCorrected(1, :, :) = diag([0.25, 0.25]);

%%
% yTrue = zeros(length(tSteps),1);
% for k=1:length(tSteps)
%     time = tSteps(k);
%     yTrue(k) = ObsOperator(trueState);
% end
% 
% figure
% hold all
% 
% line(tSteps, obsShockPos)
% line(tSteps, yTrue)

%%
for k=2:Nsteps
    % Let k denote the current time.
    
    

    time = tSteps(k);
%     yTrue(k) = ObsOperator(trueState);
%     yMeas(k) = yTrue(k) + normrnd(0, sqrt(R));

    ukf.ProcessNoise = diag([0.025, 0.025])/(20*time);
    ukf.MeasurementNoise = R;

    % Residuals (or innovations): Measured output - Predicted output
    e(k) = yMeas(k) - ObsOperator(ukf.State); % ukf.State is x[k|k-1] at this point
    % Incorporate the measurements at time k into the state estimates by
    % using the "correct" command. This updates the State and StateCovariance
    % properties of the filter to contain x[k|k] and P[k|k]. These values
    % are also produced as the output of the "correct" command.
    [xCorrectedUKF(k,:), PCorrected(k,:,:)] = correct(ukf,yMeas(k));
    % Predict the states at next time step, k+1. This updates the State and
    % StateCovariance properties of the filter to contain x[k+1|k] and
    % P[k+1|k]. These will be utilized by the filter at the next time step.
    predict(ukf);
end

%%

estimatedParsUKF = xCorrectedUKF;
estimatedParsStdUKF = zeros(size(estimatedParsUKF));
estimatedParsStdUKF(:, 1) = sqrt(PCorrected(:,1,1));
estimatedParsStdUKF(:, 2) = sqrt(PCorrected(:,2,2));


save("UKFResults2.mat", "estimatedParsStdUKF", "estimatedParsUKF")


%%
figure();
subplot(2,1,1);
plot(tSteps,repelem(trueState(1), Nsteps),tSteps,xCorrectedUKF(:,1));
legend('True','UKF estimate')
% ylim([-2.6 2.6]);
ylabel('x_1');
subplot(2,1,2);
plot(tSteps,repelem(trueState(2), Nsteps),tSteps,xCorrectedUKF(:,2));
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('x_2');

%%
figure();
plot(tSteps, e);
xlabel('Time [s]');
ylabel('Residual (or innovation)');

%%
mean(e)

%%
[xe,xeLags] = xcorr(e,'coeff'); % 'coeff': normalize by the value at zero lag
% Only plot non-negative lags
idx = xeLags>=0;
figure();
plot(xeLags(idx),xe(idx));
xlabel('Lags');
ylabel('Normalized correlation');
title('Autocorrelation of residuals (innovation)');

%%
eStates = double(trueState)-xCorrectedUKF;
figure();
subplot(2,1,1);
plot(tSteps,eStates(:,1),...               % Error for the first state
    tSteps, sqrt(PCorrected(:,1,1)),'r', ... % 1-sigma upper-bound
    tSteps, -sqrt(PCorrected(:,1,1)),'r');   % 1-sigma lower-bound
xlabel('Time [s]');
ylabel('Error for state 1');
title('State estimation errors');
subplot(2,1,2);
plot(tSteps,eStates(:,2),...               % Error for the second state
    tSteps,sqrt(PCorrected(:,2,2)),'r', ...  % 1-sigma upper-bound
    tSteps,-sqrt(PCorrected(:,2,2)),'r');    % 1-sigma lower-bound
xlabel('Time [s]');
ylabel('Error for state 2');
legend('State estimate','1-sigma uncertainty bound',...
    'Location','Best');

%%
[xeStates1,xeStatesLags1] = xcorr(eStates(:,1),'coeff'); % 'coeff': normalize by the value at zero lag
[xeStates2,xeStatesLags2] = xcorr(eStates(:,2),'coeff'); % 'coeff'
% Only plot non-negative lags
idx = xeStatesLags1>=0;
figure();
subplot(2,1,1);
plot(xeStatesLags1(idx),xeStates1(idx));
xlabel('Lags');
ylabel('For state 1');
title('Normalized autocorrelation of state estimation error');
subplot(2,1,2);
plot(xeStatesLags2(idx),xeStates2(idx));
xlabel('Lags');
ylabel('For state 2');