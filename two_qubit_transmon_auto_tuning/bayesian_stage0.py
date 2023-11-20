import numpy as np
from matplotlib import pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
import scipy

## Bayesian Optimisation for a random Objective function-----example

#Define objective function - this gives the real underlying evaluations
def objective(x, noise=0.1):
 n = np.random.normal(loc=0, scale=noise)
 out:float = x**2 * np.sin(5 * np.pi * x)**6.0 + n
 return out

# surrogate or approximation for the objective function
def surrogate(model, X):
 return model.predict(X, return_std=True)

#plot real observations vs surrogate function
def plot(X, y, model):
 plt.scatter(X, y)
 Xsamples = np.asarray(np.arange(0, 1, 0.001))
 Xsamples = Xsamples.reshape(len(Xsamples), 1)
 ysamples, _ = surrogate(model, Xsamples)
 plt.plot(Xsamples, ysamples)
 plt.show()

# probability of improvement acquisition function
def acquisition(X, Xsamples, model):
 # calculate the best surrogate score found so far
 yhat, _ = surrogate(model, X)
 best = max(yhat)
 # calculate mean and stdev via surrogate function
 mu, std = surrogate(model, Xsamples)
#  mu = mu[:, 0]
 # calculate the probability of improvement
 probs = scipy.stats.norm.cdf((mu - best) / (std+1E-9))
 return probs

# optimize the acquisition function
def opt_acquisition(X, y, model):
 # random search, generate random samples
 Xsamples = np.random.random(100)
 Xsamples = Xsamples.reshape(len(Xsamples), 1)
 # calculate the acquisition function for each sample
 scores = acquisition(X, Xsamples, model)
 # locate the index of the largest scores
 ix = np.argmax(scores)
 return Xsamples[ix, 0]

# %% 

#Define sample points 
X = np.random.random(100)
y = np.asarray([objective(x) for x in X])

X = X.reshape([len(X), 1])
y = y.reshape([len(y), 1])
# f = [objective(x,noise=0) for x in X]

#define model
model = GaussianProcessRegressor()
model.fit(X, y)
# plot(X,y,model)

#define optimization process

for i in range(100):
 x = opt_acquisition(X,y,model)
 actual = objective(x)
 est, _ = surrogate(model, [[x]])
#  print('>x=%.3f, f()=%3f, actual=%.3f' % (x, est, actual))
 X = np.vstack((X, [[x]]))
 y = np.vstack((y, [[actual]]))
 # update the model
 model.fit(X, y)

# plot all samples and the final surrogate function
plot(X, y, model)
plot(X, [objective(x,0) for x in X], model)
# best result
ix = np.argmax(y)
print('Best Result: x=%.3f, y=%.3f' % (X[ix], y[ix]))





