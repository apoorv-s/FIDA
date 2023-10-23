import numpy as np
from clawpack.riemann.shallow_roe_with_efix_2D_constants import depth, x_momentum, y_momentum, num_eqn
import clawpack.petclaw as pyclaw
from clawpack import riemann
import scipy.stats as stats
import scipy.io
from datetime import datetime


# ## Model functions


def shallowWaterEqnPropagation_2d(parametersVec, timeStep, forwardModelParameters):

    from clawpack.riemann.shallow_roe_with_efix_2D_constants import depth, x_momentum, y_momentum, num_eqn

    delta_t = forwardModelParameters[0];
    nInterStep = forwardModelParameters[1];
    num_cells = forwardModelParameters[2];
    xDomain = forwardModelParameters[3];
    yDomain = forwardModelParameters[4];

    x = pyclaw.Dimension(xDomain[0], xDomain[1], num_cells[0], name='x');
    y = pyclaw.Dimension(yDomain[0], yDomain[1], num_cells[1], name='y');

    domain = pyclaw.Domain([x, y]);

    state = pyclaw.State(domain, num_eqn);
    
    H_in = parametersVec[0];
    H_out = parametersVec[1];
    radius = parametersVec[2];
    x0 = parametersVec[3];
    y0 = parametersVec[4];
    grav = parametersVec[5];
    xx = forwardModelParameters[5];
    yy = forwardModelParameters[6];
    
    state.problem_data['grav'] = grav;

    r = np.sqrt((xx-x0)**2 + (yy-y0)**2);
    state.q[depth, :, :] = H_in*(r<=radius) + H_out*(r>radius);
    state.q[x_momentum, :, :] = 0.;
    state.q[y_momentum, :, :] = 0.;

    solver=pyclaw.ClawSolver2D(riemann.shallow_roe_with_efix_2D);
    solver.limiters = pyclaw.limiters.tvd.MC;
    solver.dimensional_split=1

    solver.all_bcs = pyclaw.BC.extrap;

    claw = pyclaw.Controller();
    claw.keep_copy = True;
    claw.solver = solver;
    claw.solution = pyclaw.Solution(state, domain);
    claw.tfinal = timeStep;
    claw.num_output_times = nInterStep*int(timeStep/delta_t);
    claw.output_format = None;

    data = claw.run();
    print("******************")
    return claw.frames[-1].q;



def shockRadiusExact(parameterVector, timeStep, forwardModelParameters):
    
    frame = shallowWaterEqnPropagation_2d(parameterVector, timeStep, forwardModelParameters);
    hFrame = frame[0,...]; #extract frame corresponding to the depth
    
    
    cutOffParameter = 0.01;
    (vx,vy) = np.gradient(hFrame)
    vs = np.sqrt(vx**2 + vy**2)
    shockFrame = vs > cutOffParameter;
    shockFrame = 1*shockFrame; # shock Frame
    
    xx = forwardModelParameters[5];
    yy = forwardModelParameters[6];
    
    nonZeroIndices = np.nonzero(shockFrame);
    xNonZero = xx[nonZeroIndices];
    yNonZero = yy[nonZeroIndices];
    radiusValues = np.sqrt(xNonZero**2 + yNonZero**2);
    
    return max(radiusValues)



def generate_XY_grid(xDomain, yDomain, numCells):
    x = pyclaw.Dimension(xDomain[0], xDomain[1], num_cells[0], name='x');
    y = pyclaw.Dimension(yDomain[0], yDomain[1], num_cells[1], name='y');

    domain = pyclaw.Domain([x, y]);

    state = pyclaw.State(domain, 1);
    xx, yy = state.grid.p_centers;
    
    return [xx, yy];


# ## PF helper functions




def observedShockRadius(parameterVector, timeStep, forwardModelParameters, observationParameters):
    frame = shallowWaterEqnPropagation_2d(parameterVector, timeStep, forwardModelParameters);
    hFrame = frame[0,...]; #extract frame corresponding to the depth
    cutOffParameter = 0.01;
    (vx,vy) = np.gradient(hFrame)
    vs = np.sqrt(vx**2 + vy**2)
    shockFrame = vs > cutOffParameter;
    shockFrame = 1*shockFrame; # shock Frame
    
    xx = forwardModelParameters[5];
    yy = forwardModelParameters[6];
    
    obsNoiseMean = observationParameters[0];
    obsNoiseStd = observationParameters[1];
    nonZeroIndices = np.nonzero(shockFrame);
    
    xNonZero = xx[nonZeroIndices];
    yNonZero = yy[nonZeroIndices];
    
    xNonZero = xNonZero + np.random.normal(obsNoiseMean, obsNoiseStd, np.shape(xNonZero));
    yNonZero = yNonZero + np.random.normal(obsNoiseMean, obsNoiseStd, np.shape(yNonZero));
    
    radiusValues = np.sqrt(xNonZero**2 + yNonZero**2);
    
    return max(radiusValues)





def initialParametersVecGenerator(nParticles, nParameters, parametersRange):
    parametersVec = np.zeros((nParticles, nParameters));
    for i in range(0, nParticles):
        for j in range(0, nParameters):
            parametersVec[i, j] = np.random.uniform(parametersRange[j][0], parametersRange[j][1]);

    return parametersVec


# ## Particle Filter and Resampling



def ParticleFilter(nParticles, timeStep, prevParameterVectors, observation, forwardModelParameters,
                   forwardParModelParameters, likelihoodParameters):
#     input: nParticles (number of particles), timeStep, prevStateVectors (matrix comprising of nParticles 
#            previous state vectors as rows, parameters in present case), observation (observations  
#            for the current step), forwardModelParameters, (parameters for the artificial dynamics, covMat),
#            likelihoodParameters (parameters (standard deviation) for the likelihood function)
#     Output: likelihoodVec (likelihood for particles, normalized)
    
    # Propagating the states(artificial dynamics added to each parameter)
    currentParameterVectors = prevParameterVectors + np.random.multivariate_normal(forwardParModelParameters[0], 
                                                                forwardParModelParameters[1]/(timeStep*20),
                                                                nParticles);
    # Estimating likelihood 
    likelihoodVec = np.zeros((nParticles, ));
    for i in range(0, nParticles):
        tempParticleObs = shockRadiusExact(currentParameterVectors[i, :], timeStep, 
                                           forwardModelParameters);
        likelihoodVec[i] = stats.norm(observation, likelihoodParameters[0]).pdf(tempParticleObs);
        
    likelihoodVec = likelihoodVec/np.sum(likelihoodVec);
        
    
    return [currentParameterVectors, likelihoodVec]



def ResamplingWithArtificialDynamics(nParticles, currentParameterVectors, priorVec, likelihoodVec):
#     input: nParticles (number of particles), currentStateVectors (matrix comprising of nParticles 
#            current state vectors as rows, parameters in present case), priorVec (describing   
#            for the current step), likelihoodVec (Likelihood vector from PF)
#     Output: resampledStates (states after resampling), newPosteriorVec (uniform pdf after resampling),
#            frequencyVec (vector indicating how many time a particle is resampled)
    
    posteriorVec = priorVec*likelihoodVec;
    posteriorVec = posteriorVec/np.sum(posteriorVec); # Normalizing the pdf
    
    if(np.round(sum(posteriorVec), 2) != np.round(1, 2)):
        print("error 1 in resampling");

    cumProb = np.cumsum(posteriorVec);
    resampledParameterVectors = np.zeros(np.shape(currentParameterVectors));
    frequencyVec = np.zeros((nParticles, ));

    uSeed = np.random.uniform()/nParticles;
    k = 0;
    for j in range(0, nParticles):
        tempRandomU = uSeed + (j/nParticles);
        while(tempRandomU > cumProb[k]):
            k = k + 1;
            if(k >= nParticles):
                print("resampling error")
        
        resampledParameterVectors[j, :] = currentParameterVectors[k, :];
        frequencyVec[k] = frequencyVec[k] + 1;

    newPosteriorVec = np.ones(posteriorVec.size)/nParticles;
    
    return [resampledParameterVectors, newPosteriorVec, frequencyVec]



nParticlesVec = [200, 500, 1000];
nParameters = 6;

tDomain = [0, 0.3];
delta_t = 0.05;
tSteps = np.arange(tDomain[0], tDomain[1], delta_t);
tSteps = np.round(tSteps, 4);

xDomain = [-2.5, 2.5];
yDomain = [-2.5, 2.5];

num_cells = (250, 250);
particleSize = (num_eqn, num_cells[0], num_cells[1]);

hInOriginal = 25;
hOutOriginal = 1;
radiusOriginal = 0.5;
xPosOriginal = 0;
yPosOriginal = 0;
gravOriginal = 10;

hIn_range = [hInOriginal - 7, hInOriginal + 2];
hOut_range = [hOutOriginal - 0, hOutOriginal + 0];
radius_range = [radiusOriginal - 0, radiusOriginal + 0];
xPos_range = [xPosOriginal - 0, xPosOriginal + 0];
yPos_range = [yPosOriginal - 0, yPosOriginal + 0];
grav_range = [gravOriginal - 4, gravOriginal + 1];

parameterOrder = ["hIn", "hOut", "radius", "xPos", "yPos", "grav"];
parametersOriginal = [hInOriginal, hOutOriginal, radiusOriginal, xPosOriginal, yPosOriginal, gravOriginal];
parametersRange = [hIn_range, hOut_range, radius_range, xPos_range, yPos_range, grav_range];

initialParameterVec = initialParametersVecGenerator(max(nParticlesVec), nParameters, parametersRange);

nInterSteps = 10;
obsNoiseMean = 0;
obsNoiseStd = 0.1;
likelihoodStd = obsNoiseStd;

artificialDynamicsMean = np.zeros((nParameters, ));
artificialDynamicsCov = np.diag([np.square((hIn_range[1] - hIn_range[0])/10),
                                 np.square((hOut_range[1] - hOut_range[0])/10),
                                 np.square((radius_range[1] - radius_range[0])/10),
                                 np.square((xPos_range[1] - xPos_range[0])/10),
                                 np.square((yPos_range[1] - yPos_range[0])/10),
                                 np.square((grav_range[1] - grav_range[0])/10)]);

forwardParModelParameters = [artificialDynamicsMean, artificialDynamicsCov];
observationParameters = [obsNoiseMean, obsNoiseStd];
likelihoodParameters = [likelihoodStd];

[xx, yy] = generate_XY_grid(xDomain, yDomain, num_cells);

forwardModelParameters = [delta_t, nInterSteps, num_cells, xDomain, yDomain, xx, yy];



for i in range(0, len(nParticlesVec)):
    nParticles = nParticlesVec[i];
    print(nParticles)

    propagatedParVecMat = np.zeros((len(tSteps), nParticles, nParameters));
    resampledParVecMat = np.zeros((len(tSteps), nParticles, nParameters));
    particlesFreqMat = np.zeros((len(tSteps), nParticles));
    estimatedParMatWeighed = np.zeros((len(tSteps), nParameters));
    estimatedParMatResampled = np.zeros((len(tSteps), nParameters));
    likelihoodMat = np.zeros((len(tSteps), nParticles));
    posteriorMat = np.zeros((len(tSteps), nParticles));
    obsShockPos = np.zeros((len(tSteps), ));

    propagatedParVecMat[0,...] = initialParameterVec[0:nParticles,...];
    resampledParVecMat[0,...] = propagatedParVecMat[0,...];
    particlesFreqMat[0,...] = np.ones((1, nParticles));
    posteriorMat[0,...] = np.ones((1, nParticles))/nParticles;

    estimatedParMatWeighed[0, :] = np.matmul(np.transpose(propagatedParVecMat[0,...]), posteriorMat[0,...]);
    estimatedParMatResampled[0, :] = estimatedParMatWeighed[0, :];
    
    for i in range(1, len(tSteps)):
        timeStep = tSteps[i];
        print(timeStep);

        obsShockPos[i] = observedShockRadius(parametersOriginal, timeStep, forwardModelParameters,
                                             observationParameters);

        [propagatedParVecMat[i,...],
         likelihoodMat[i,...]] = ParticleFilter(nParticles, timeStep,
                                                resampledParVecMat[i - 1,...], obsShockPos[i],
                                                forwardModelParameters, forwardParModelParameters,
                                                likelihoodParameters);

        estimatedParMatWeighed[i, :] = np.matmul(np.transpose(propagatedParVecMat[i,...]), likelihoodMat[i,...]);

        [resampledParVecMat[i,...], posteriorMat[i,...],
         particlesFreqMat[i,...]] = ResamplingWithArtificialDynamics(nParticles, 
                                                                     propagatedParVecMat[i,...],
                                                                     posteriorMat[i - 1,...], likelihoodMat[i,...]);
        estimatedParMatResampled[i, :] = np.matmul(np.transpose(resampledParVecMat[i,...]), posteriorMat[i,...]);
        
    filename1 = 'ex3_' + str(nParticles) + '_' + datetime.today().strftime('%d-%m-%y') + '_artificialDynamics_jackKnife.mat';
    filename2 = 'ex3_' + str(nParticles) + '_' + datetime.today().strftime('%d-%m-%y') + '_auxilliary_artificialDynamics_jackKnife.mat';

    scipy.io.savemat(filename1, dict(nParameters = nParameters, nParticles = nParticles, tSteps = tSteps, 
                                    parametersRange = parametersRange, parametersOriginal = parametersOriginal,
                                    resampledParVecMat = resampledParVecMat, posteriorMat = posteriorMat, 
                                    estimatedParMatResampled = estimatedParMatResampled, 
                                    observationParameters = observationParameters, parameterOrder = parameterOrder, 
                                    obsShockPos = obsShockPos, likelihoodStd = likelihoodStd, 
                                    artificialDynamicsCov = artificialDynamicsCov));


    scipy.io.savemat(filename2, dict(propagatedParVecMat = propagatedParVecMat, particlesFreqMat = particlesFreqMat,
                                    estimatedParMatWeighed = estimatedParMatWeighed, likelihoodMat = likelihoodMat,
                                    xDomain = xDomain, yDomain = yDomain, num_cells = num_cells,
                                    nInterSteps = nInterSteps, forwardModelParameters = forwardModelParameters, xx = xx, yy = yy)); 
    






