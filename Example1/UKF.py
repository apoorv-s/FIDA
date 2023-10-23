import numpy as np
import numpy.linalg

class Ramp:
    def __init__(self, lmd:float, xR:float) -> None:
        self.lmd = lmd
        self.uL = 2
        self.uR = 1
        self.xR = xR
        
        self.beta = 2.0/self.lmd
        self.alpha = (self.uL-self.uR)/(self.xR)
        self.tBreak = 1/(self.beta*self.alpha)
        
    def AnalyticSolution(self, x:float, t:float) -> float:
        u = np.inf
        if t < self.tBreak:
            if x < self.beta*self.uL*t:
                u = self.uL
            elif self.beta*self.uL*t <= x < self.xR + self.beta*self.uR*t:
                u = (self.uL - self.alpha*x)/(1 - self.alpha*self.beta*t)
            else:
                u = self.uR
        else:
            skPos = self.GetShockPosition(t)
            if x < skPos:
                u = self.uL
            elif x == skPos:
                u = (self.uL + self.uR)/2
            else:
                u = self.uR
        return u
            
    def AnalyticalSolutionVectorized(self, xSteps:np.ndarray, t:float) -> np.ndarray:
        uVec = np.zeros(len(xSteps))
        for i in range(len(xSteps)):
            uVec[i] = self.AnalyticSolution(xSteps[i], t)
        return uVec
            
    def GetShockPositionForDA(self, t:float) -> float:
        if t < self.tBreak:
            skPos = np.round((self.beta*self.uL*t + self.xR + self.beta*self.uR*t)/2, 4)
        if t >= self.tBreak:
            sigma = (self.uL + self.uR)/self.lmd
            skPos = sigma*t + (self.xR/2)
        return skPos
    
def GenerateNoisyShockObservation(rampObj:Ramp, t:float, noiseMean:float, noiseStd:float) -> float:
    return rampObj.GetShockPositionForDA(t) + np.random.normal(noiseMean, noiseStd)

def GetSigmaPoints(xEstimate, stateCov, w0):
    lenState = len(xEstimate)
    nSigmaPoints = 2*lenState + 1
    
    # Numpy cholesky of the for A = L*L.T, so columns of L to be used, (Unscented Filtering and Nonlinear Estimation, Julier and Uhlmann, 2004)
    sqrtCovX = np.linalg.cholesky(stateCov*(lenState/(1 - w0)))
    
    sigmaPoints = np.zeros((nSigmaPoints, lenState))
    weights = np.zeros(nSigmaPoints)
    
    sigmaPoints[0] = xEstimate
    weights[0] = w0
    
    for i in range(0, lenState):
        sigmaPoints[i + 1] = xEstimate + sqrtCovX[:, i]
        weights[i + 1] = (1 - w0)/(2*lenState)
        sigmaPoints[i + lenState + 1] = xEstimate - sqrtCovX[:, i]
        weights[i + lenState + 1] = (1 - w0)/(2*lenState)
    
    return sigmaPoints, weights


# def GetSigmaPoints2(xEstimate, stateCov):
#     lenState = len(xEstimate)
#     nSigmaPoints = 2*lenState
#     sqrtCovX = np.linalg.cholesky(stateCov*(lenState))
    
#     sigmaPoints = np.zeros((nSigmaPoints, lenState))
#     sigmaWeights = np.ones(nSigmaPoints)*(1/(2*lenState))

#     for i in range(0, lenState):
#         sigmaPoints[i] = xEstimate + sqrtCovX[:, i]
#         sigmaPoints[i + lenState] = xEstimate - sqrtCovX[:, i]
    
#     return sigmaPoints, sigmaWeights

def GetSigmaPoints(stateEstimate, stateCov, alpha, beta, kappa):
    nStates = len(stateEstimate)
    nSigmaPoints = 2*nStates + 1
    sigmaPointsMat = np.zeros((nStates, nSigmaPoints))
    sigmaWeightsMean = np.zeros(nSigmaPoints)
    sigmaWeightsCov = np.zeros(nSigmaPoints)
    
    lambdaForSigmaPoints = (alpha**2)*(nStates + kappa) - nStates
    
    # Numpy cholesky of the for A = L*L.T, so columns of L to be used, (Unscented Filtering and Nonlinear Estimation, Julier and Uhlmann, 2004)
    sqrtStateCov = np.linalg.cholesky(stateCov*(nStates + lambdaForSigmaPoints))
    
    sigmaPointsMat[:, 0] = stateEstimate
    sigmaWeightsMean[0] = lambdaForSigmaPoints/(nStates + lambdaForSigmaPoints)
    sigmaWeightsCov[0] = (lambdaForSigmaPoints/(nStates + lambdaForSigmaPoints)) + (1 - (alpha**2) + beta)
    
    for i in range(nStates):
        sigmaPointsMat[:, i + 1] = stateEstimate + sqrtStateCov[:, i]
        sigmaWeightsMean[i + 1] = 1/(2*(nStates + lambdaForSigmaPoints))
        sigmaWeightsCov[i + 1] = 1/(2*(nStates + lambdaForSigmaPoints))
        
        sigmaPointsMat[:, i + nStates + 1] = stateEstimate - sqrtStateCov[:, i]
        sigmaWeightsMean[i + nStates + 1] = 1/(2*(nStates + lambdaForSigmaPoints))
        sigmaWeightsCov[i + nStates + 1] = 1/(2*(nStates + lambdaForSigmaPoints))
        
    return sigmaPointsMat, sigmaWeightsMean, sigmaWeightsCov
    

def GetWeightedMean(pointsMat, weights):
    # pointsMat : len(weight) x nStates
    return np.array(np.average(pointsMat, axis = 0, weights = weights))

def GetWeighteCov(pointsMat, weights):
    # pointsMat : len(weight) x nStates
    mean = GetWeightedMean(pointsMat, weights)
    cov = np.zeros((mean.size, mean.size))
    for i in range(len(pointsMat)):
        cov = cov + weights[i]*np.outer(pointsMat[i] - mean, pointsMat[i] - mean)
    return cov
    
def GetCrossCov(pointsMat1, pointsMat2, weights):
    # pointsMat1,2 : len(weight) x nStates
    mean1 = GetWeightedMean(pointsMat1, weights)
    mean2 = GetWeightedMean(pointsMat2, weights)
    cov = np.zeros((mean1.size, mean2.size))
    for i in range(len(weights)):
        cov = cov + weights[i]*np.outer(pointsMat1[i] - mean1, pointsMat2[i] - mean2)
    return cov

def ForwardModel(x):
    return x + np.random.normal(0, 0.1, x.shape)

def ObservationModel(x, time):
    Ramp(*x).GetShockPositionForDA(time)
    
def UnscentedKF(prevStateEstimate, prevStateCov, forwardModel, observationModel, observation, time, w0 = 0.5):
    
    ## 1. Generate Sigma Points 
    sigmaPoints, sigmaWeights = GetSigmaPoints(prevStateEstimate, prevStateCov, w0)
    
    ## 2. Instantiating each point through process model
    transformedSigmaPoints = forwardModel(sigmaPoints)
    
    ## 3. Predicted Mean
    predictedMean = GetWeightedMean(transformedSigmaPoints, sigmaWeights)

    ## 4. Predicted Covariance
    predictedCov = GetWeighteCov(transformedSigmaPoints, sigmaWeights)

    ## 5. Instatiating through the observation Model
    obsFromSigmaPoints = np.zeros(len(sigmaWeights))
    for iSigmaWeight in range(len(sigmaWeights)):
        obsFromSigmaPoints[iSigmaWeight] = observationModel(sigmaPoints[iSigmaWeight], time)

    ## 6. Predicted Observation
    predictedObs = GetWeightedMean(obsFromSigmaPoints, sigmaWeights)

    ## 7. Innovation Covariance
    innovationCov = GetWeighteCov(obsFromSigmaPoints, sigmaWeights)

    ## 8. Cross Covariance Matrix
    crossCovMat = GetCrossCov(transformedSigmaPoints, obsFromSigmaPoints, sigmaWeights)

    ## 9. Update
    W = np.dot(crossCovMat, np.linalg.inv(innovationCov))
    updatedEstimate =  predictedMean + np.dot(W, (observation - predictedObs)).reshape(predictedMean.shape)
    updatedCov = predictedCov - np.dot(np.dot(W, innovationCov), W.T)
    
    return updatedEstimate, predictedMean, updatedCov, predictedCov






