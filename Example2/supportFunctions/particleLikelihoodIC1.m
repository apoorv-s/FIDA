function likelihood = particleLikelihoodIC1(particle, obsShockPos, observationStd)
% Gives likelihood of a given particle given the observed data and the
% particle

obsShockPos = obsShockPos{1};
particleShockPos = shockPosFromParticleIC1(particle);
likelihood = 1;

for i = 1:length(obsShockPos)
    likelihood = likelihood*max(normpdf(particleShockPos, obsShockPos(i), observationStd));
end

if(likelihood == 1)
    disp("error in likelihood")
    disp(obsShockPos)
    disp(particleShockPos)
end

end