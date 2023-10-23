function [resampledStateVec, resampledParameterMat, newPosteriorVec,...
    frequencyVec, uEstimateResampled, uResampledMat] = Resampling(nParticles, currentParameterMat, ...
    currentStateVec, priorVec, likelihoodVec, uWeighedMat)
% input: nParticles (number of particles), currentStateVectors (matrix comprising of nParticles 
%        current state vectors as rows, parameters in present case), priorVec (describing   
%        for the current step), likelihoodParameters (parameters for the likelihood function)
% Output: likelihoodVec (likelihood for particles, not normalized)

posteriorVec = priorVec.*likelihoodVec;
posteriorVec = posteriorVec/sum(posteriorVec);

if(round(sum(posteriorVec), 2) ~= round(1, 2))
        disp("error 1 in resampling");
end

cumProb = cumsum(posteriorVec);

resampledStateVec = cell(size(currentStateVec));
resampledParameterMat = zeros(size(currentParameterMat));
uResampledMat = zeros(size(uWeighedMat));
frequencyVec = zeros(nParticles, 1);

uSeed = rand/nParticles;
k = 1;
for j = 1:nParticles
    tempRandomU = uSeed + ((j - 1)/nParticles);
    while(tempRandomU > cumProb(k))
        k = k + 1;
        if(k > nParticles)
            disp("resampling error")
        end
    end
    
    resampledParameterMat(:, j) = currentParameterMat(:, k);
    resampledStateVec{j} = currentStateVec{k};
    uResampledMat(:, j) = uWeighedMat(:, k);
    frequencyVec(k) = frequencyVec(k) + 1;
end

newPosteriorVec = ones(size(posteriorVec))/nParticles;
uEstimateResampled = uResampledMat*newPosteriorVec;

end