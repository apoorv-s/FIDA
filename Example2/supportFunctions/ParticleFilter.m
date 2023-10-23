function [currentStateVec, currentParameterVec, likelihoodVec,...
    uEstimateWeighed, uPropagatedMat] = ParticleFilter(nParticles, prevParametersVec, time, observations, ...
    forwardModelParameters, likelihoodParameters, xSteps)
%  input: nParticles (number of particles), timeStep, prevParametersVec (matrix comprising of nParticles 
%         previous parameter vectors as rows), observation (observations  
%         for the current step), forwardModelParameters, (parameters for the artificial dynamics, covMat),
%         likelihoodParameters (parameters (standard deviation) for the likelihood function)
%  Output: currentStateVec (chebfun objects for the current state), currentParameterVec (perturbed parameters for the 
%          current state), likelihoodVec (likelihood for particles, normalized)

likelihoodVec = zeros(nParticles, 1);
currentStateVec = cell(nParticles, 1);
uPropagatedMat = zeros(length(xSteps), nParticles);

currentParameterVec = prevParametersVec + mvnrnd(forwardModelParameters{1}, forwardModelParameters{2}/(time*10), nParticles)';

parfor i = 1:nParticles
% for i = 1:nParticles
    disp(i)
    [tempParticle, ] = IBSolverIC1(currentParameterVec(:, i), time);
    currentStateVec{i} = tempParticle;
    likelihoodVec(i) = particleLikelihoodIC1(tempParticle, observations, likelihoodParameters);
    uPropagatedMat(:, i) = tempParticle(xSteps);
end

likelihoodVec = likelihoodVec/sum(likelihoodVec);
uEstimateWeighed = uPropagatedMat*likelihoodVec;

end

