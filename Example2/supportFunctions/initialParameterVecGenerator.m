function parameterVec = initialParameterVecGenerator(nParticles, nParameters, parametersRange)
% Input: nParticles (number of particles), nParameters (number of
%        parameters), parametersRange (range of parameters)
% Output: parameterVec (nParticles realization of parameter vectors from uniform distribution)


parameterVec = zeros(nParameters, nParticles);

for i = 1:nParameters
    parameterVec(i, :) = ((parametersRange(i, 2) - parametersRange(i, 1))*rand(nParticles, 1) ...
                           + parametersRange(i, 1));
end

end