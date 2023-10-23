function SbatchFunction()

%%% Particle Filter: MATLAB Implementation
% Case with multiple shocks: IC1 :
%
% chebfun({sin(tempAlpha*tmpVar), tempBeta, sin(tempGamma*tmpVar)}, ...
%                     icDomain, 'splitting', 'on');

clear
clc

nParticlesVec = [200, 500, 1000];

addpath('/home/darve/apoorv1/InviscidBurgersDA/matlabLab/chebfun');
addpath('/home/darve/apoorv1/InviscidBurgersDA/ArtificialDynamics/Example2/supportFunctions');

nParameters = 8;

tDomain = [0, 2];
delta_t = 0.1;
tSteps = tDomain(1):delta_t:tDomain(2);

xDomain = [-4, 4];
delta_x = 0.01;
xSteps = xDomain(1):delta_x:xDomain(2);
xSteps = xSteps';

alpha_original = pi; % par 1
beta_original = 2; % par 2
gamma_original = pi; % par 3

a_original = -4; % par 4
b_original = -1; % par 5
c_original = 1; % par 6
d_original = 4; % par 7

lambda_original = 2; % par 8

alpha_range = [3, 3.5];
beta_range = [beta_original - 0.1, beta_original + 0.4];
gamma_range = [3, 3.5];

a_range = [a_original - 0, a_original + 0];
b_range = [b_original - 0, b_original + 0];
c_range = [c_original - 0, c_original + 0];
d_range = [d_original - 0, d_original + 0];

lambda_range = [lambda_original - 0.4, lambda_original + 0.1];

parametersOriginal = [alpha_original; beta_original; gamma_original; a_original; b_original; c_original; d_original; lambda_original];
paramatersRange = [alpha_range; beta_range; gamma_range; a_range; b_range; c_range; d_range; lambda_range];

initialParameterVec = initialParameterVecGenerator(max(nParticlesVec), nParameters, paramatersRange);

save("runData.mat", 'nParticlesVec', "initialParameterVec", 'nParameters', 'xSteps', ...
        'tSteps', 'parametersOriginal', 'paramatersRange');

for nParticleIndex = 1:length(nParticlesVec)

    nParticles = nParticlesVec(nParticleIndex);
    disp("nParticles")
    disp(nParticles)

    propagatedParVecMat = zeros(nParameters, nParticles, length(tSteps));
    resampledParVecMat = zeros(nParameters, nParticles, length(tSteps));
    particlesFreqMat = zeros(nParticles, length(tSteps));
    likelihoodMat = zeros(nParticles, length(tSteps));
    posteriorMat = zeros(nParticles, length(tSteps));
    estimatedParMatWeighed = zeros(nParameters, length(tSteps));
    estimatedParMatResampled = zeros(nParameters, length(tSteps));
    observationMat = cell(1, length(tSteps));
    propagatedStateMat = cell(nParticles, length(tSteps));
    resampledStateMat = cell(nParticles, length(tSteps));
    trueParticles = cell(1, length(tSteps));
    uWeightedVecMat = zeros(length(xSteps), length(tSteps));
    uResampledVecMat = zeros(length(xSteps), length(tSteps));
    uWeightedMat = zeros(length(xSteps), nParticles, length(tSteps));
    uResampledMat = zeros(length(xSteps), nParticles, length(tSteps));

    uAnalytic = zeros(length(xSteps), length(tSteps));

    
    propagatedParVecMat(:, :, 1) = initialParameterVec(:, 1:nParticles);
    resampledParVecMat(:, :, 1) = propagatedParVecMat(:, :, 1);
    particlesFreqMat(:, 1) = ones(nParticles, 1);
    posteriorMat(:, 1) = ones(nParticles, 1)./nParticles;
    estimatedParMatWeighed(:, 1) = propagatedParVecMat(:, :, 1)*posteriorMat(:, 1);
    estimatedParMatResampled(:, 1) = estimatedParMatWeighed(:, 1);
    
    [propagatedStateMat(:, 1), uWeightedVecMat(:, 1)] = initialParticlesGeneratorIC1(nParticles, ...
        propagatedParVecMat(:, :, 1), xSteps);
    resampledStateMat(:, 1) = propagatedStateMat(:, 1);
    uResampledVecMat(:, 1) = uWeightedVecMat(:, 1);
    trueParticles{1} = IBSolverIC1(parametersOriginal, 0);
    uAnalytic(:, 1) = trueParticles{1}(xSteps);
    
    obsNoiseMean = 0;
    obsNoiseStd = 0.1;
    likelihoodStd = obsNoiseStd;
    artificialDynamicsMean = zeros(nParameters, 1);
    artificialDynamicsCov = diag(power([alpha_range(2) - alpha_range(1), ...
                                        beta_range(2) - beta_range(1), ...
                                        gamma_range(2) - gamma_range(1), ...
                                        a_range(2) - a_range(1), ...
                                        b_range(2) - b_range(2), ...
                                        c_range(2) - c_range(1), ...
                                        d_range(2) - d_range(1), ...
                                        lambda_range(2) - lambda_range(1)], 2))/100;

    disp("artificialDynamicsCov")
    disp(artificialDynamicsCov)

    forwardModelParameters = {artificialDynamicsMean, artificialDynamicsCov};
    observationParameters = [obsNoiseMean, obsNoiseStd];
    likelihoodParameters = likelihoodStd;
    
    
    for i = 2:length(tSteps)
        time = tSteps(i);
        disp("time")
        disp(time)
        trueParticles{i} = IBSolverIC1(parametersOriginal, tSteps(i));
        uAnalytic(:, i) = trueParticles{i}(xSteps);
        observationMat{i} = observationGeneratorIC1(trueParticles{i}, ...
                                                    observationParameters);

        [propagatedStateMat(:, i), propagatedParVecMat(:, :, i),...
            likelihoodMat(:, i), uWeightedVecMat(:, i), uWeightedMat(:, :, i)] = ParticleFilter(nParticles, ...
            resampledParVecMat(:, :, i - 1), time, observationMat{i}, forwardModelParameters, likelihoodParameters, xSteps);

        estimatedParMatWeighed(:, i) = propagatedParVecMat(:, :, i)*likelihoodMat(:, i);

        [resampledStateMat(:, i), resampledParVecMat(:, :, i),...
            posteriorMat(:, i), particlesFreqMat(:, i), uResampledVecMat(:, i),...
            uResampledMat(:, :, i)] = Resampling(nParticles, propagatedParVecMat(:, :, i), ...
            propagatedStateMat(:, i), posteriorMat(:, i - 1), likelihoodMat(:, i), uWeightedMat(:, :, i));

        estimatedParMatResampled(:, i) = resampledParVecMat(:, :, i)*posteriorMat(:, i);
    end
    

    
    filename1 = strcat('ex2_', int2str(nParticles), '_artDyn_jackKnife.mat');
    filename2 = strcat('ex2_', int2str(nParticles), '_artDyn_jackKnife_aux.mat');
    
    save(filename1, 'nParticles', 'resampledParVecMat', 'posteriorMat', "estimatedParMatResampled", ...
        "uResampledVecMat", 'tSteps');
    save(filename2, 'propagatedParVecMat', 'particlesFreqMat', 'likelihoodMat', 'estimatedParMatWeighed', ...
        'observationMat', 'propagatedStateMat', 'resampledStateMat', 'trueParticles', 'uWeightedVecMat', ...
        'uWeightedMat', 'uResampledMat', 'uAnalytic');

end

end