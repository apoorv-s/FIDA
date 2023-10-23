function shockPos = shockPosFromParticleIC1(particle)
% Finds the position of the shock given the particle (chebfun)

jumpPoints = particle.ends;
jumpPoints = jumpPoints(2:length(jumpPoints) - 1);

currState_Prime = diff(particle);
currState_dPrime = diff(currState_Prime);

prime_jumpPoints = currState_Prime(jumpPoints);
dPrime_jumpPoints = currState_dPrime(jumpPoints);

shockPos = [];

for i = 1:length(jumpPoints)
    if isinf(dPrime_jumpPoints(i)) || isinf(prime_jumpPoints(i))
        shockPos = [shockPos, jumpPoints(i)];
    end
end


end