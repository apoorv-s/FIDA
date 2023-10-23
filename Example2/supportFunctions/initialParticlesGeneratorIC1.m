function [particles, uVec] = initialParticlesGeneratorIC1(nParticles, parametersVec, xSteps)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

uVec = zeros(length(xSteps), 1);

particles = cell(nParticles, 1);
for i = 1:nParticles
    tempAlpha = parametersVec(1, i);
    tempBeta = parametersVec(2, i);
    tempGamma = parametersVec(3, i);
    tempIcDomain = [parametersVec(4, i), parametersVec(5, i), parametersVec(6, i), parametersVec(7, i)];
    tmpVar = chebfun('tmpVar', tempIcDomain);
    tmpChebfun = chebfun({sin(tempAlpha*tmpVar), tempBeta, sin(tempGamma*tmpVar)}, ...
                    tempIcDomain, 'splitting', 'on');
    particles{i} = tmpChebfun;
    uVec = uVec + tmpChebfun(xSteps);
end

uVec = uVec/nParticles;

end