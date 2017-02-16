from pylab import *

n_subj = 10
xs = array([1, 5, 10, 20, 50]).astype(float)
ys = array([1, 5, 10, 20, 50]).astype(float)
sd = array([.1, .5, 1, 2, 5]).astype(float)

xs_inv = 1.0/xs
sd_inv = 1.0/sd

se_m = sd/sqrt(n_subj)
se_sd = sd/sqrt(2*n_subj)

min_sd = sd-se_sd
max_sd = sd+se_sd

f = figure()
ax1 = f.add_subplot(111)
ax1.set_xlim([0,60]); ax1.set_ylim([0,60])
ax1.errorbar(xs,ys,yerr=min_sd, color='b')
ax1.errorbar(xs,ys,yerr=sd, color='g')
ax1.errorbar(xs,ys,yerr=max_sd, color='r')
show()


### how to plot where the linearly increasing expected SD lies,
### and where real SDs could be but how we can interpret above w/ sampled SD
## still TODO above, w/e

mean = 0
std_dev = 10.0
subj = 10

vals = normal(mean, std_dev, (subj,1000))

a = rand(10000,1000)*100
