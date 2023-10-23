function obsShockPos = observationGeneratorIC1(particle, observationParameters)
% Input: particle (true particle), observationParameters (noise mean and standard deviation)
% Output: obsShockPos (Noisy shock position in a cell)


shockPos = shockPosFromParticleIC1(particle);

shockPos = shockPos + normrnd(observationParameters(1), observationParameters(2), size(shockPos));

obsShockPos = {shockPos};

end