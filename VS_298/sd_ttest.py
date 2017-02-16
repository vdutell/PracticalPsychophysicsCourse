import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import stats

rseed = 298
np.random.seed(rseed) # fixes the random numbers generated for the trial data across python runs

mean_pop = .5 # mean of population
std_pop = 1 # standard deviation of population
ntrials = 16 # number of trials in an experimental sample
nruns = 100000 # number of simulations

se_pop = std_pop/np.sqrt(ntrials) # standard error of the population  

#trial_data = np.random.normal(mean_pop, std_pop, ntrials) # generates normal data with pop stats
trial_data = np.array([-0.3968, 1.1858, 0.1122, 1.1531, 1.6512, 0.5867, 0.5016, -0.6992, 0.3109,
                     0.4996, 1.3534, 1.2568, 0.1775, 1.1165, -0.2890,  2.0588])
null_data = np.zeros(ntrials) # null data of all zeros

mean_samp = np.mean(trial_data) # mean of the trial samples taken from the population
std_samp = np.std(trial_data) # standard deviation of the trial samples taken from the population 

se_samp = std_samp/np.sqrt(ntrials) # standard error of the trial samples taken from the population

print 'number of trials per experimental sample: ' + str(ntrials)
print 'mean of the population: ' + str(mean_pop)
print 'standard deviation of the population: ' + str(mean_pop)
print 'standard error of the population: ' + str(se_pop) + '\n'

print 'mean of the trial sample: ' + str(mean_samp)
print 'standard deviation of the trial sample: ' + str(std_samp)
print 'standard error of the trial sample: ' + str(se_samp) + '\n'

# run t-test on the trial data with the null hypothesis being the null data
[t_score, p_val] = scipy.stats.ttest_ind(trial_data, null_data)
print 't score of the trial sample: ' + str(t_score) + ' with corresponding probability value of: ' + str(p_val) + '\n'

sim_data = np.random.normal(mean_pop, std_pop, (nruns,ntrials))
sim_null_data = np.zeros(sim_data.shape)

mean_sim = np.mean(sim_data)
std_sim = np.std(sim_data)

print 'number of runs in the simulation: ' + str(nruns)
print 'mean of all simulation trials: ' + str(mean_sim)
print 'standard deviation of all simulation trials: ' + str(std_sim) + '\n' 

# run t-test on the trial data with the null hypothesis being the null data, nruns number of times
[t_score_sim, p_val_sim] = scipy.stats.ttest_ind(sim_data, sim_null_data, axis=1)
print 't scores of the simulation: ' + str(t_score_sim) + ' with corresponding probability values of: ' + str(p_val_sim) + '\n'

mean_ts = np.mean(t_score_sim) # mean of all the t-scores of the simulation
mean_ps = np.mean(p_val_sim) # mean of all the prob vals of the simulation

std_ts = np.std(t_score_sim) # standard deviation of all the t-scores of the simulation
std_ps = np.std(p_val_sim) # standard deviation of all the prob vals of the simulation 

se_ts = np.std(t_score_sim)/np.sqrt(nruns) # mean of all the t-scores of the simulation
se_ps = np.std(p_val_sim)/np.sqrt(nruns) # mean of all the prob vals of the simulation 

print 'mean of all the t-scores of the simulation: ' + str(mean_ts)
print 'mean of all the two-tailed probability values of the simulation: ' + str(mean_ps) + '\n'
print 'standard deviation of all the t-scores of the simulation: ' + str(std_ts)
print 'standard deviation of all the two-tailed probability values of the simulation: ' + str(std_ps) + '\n'
print 'standard error of all the t-scores of the simulation: ' + str(se_ts)
print 'standard error of all the two-tailed probability values of the simulation: ' + str(se_ps)

f = plt.figure()
ax1 = f.add_subplot(121)
ax1.hist(t_score_sim, bins=1000)
ax1.set_title('Histogram of\n the t-values from\n all simulation runs')
ax2 = f.add_subplot(122)
ax2.hist(p_val_sim, bins=1000)
ax2.set_title('Histogram of\n the probability values\n from all simulation runs')
plt.show()
