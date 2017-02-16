from pylab import *
from scipy.optimize import curve_fit

def zahn_func1(cond, p1):
    exp = zeros(len(cond))
    for i in range(len(cond)):
        exp[i] = p1
    return exp

def zahn_func2(cond, p1, p2):
    exp = zeros(len(cond))
    for i in range(len(cond)):
        exp[i] = p1
        if cond[i] == 1 or cond[i] == 2:
            exp[i] = exp[i] + p2
        elif cond[i] == 3 or cond[i] == 4:
            exp[i] = exp[i] - p2
    return exp

def zahn_func3(cond, p1, p2, p3):
    exp = zeros(len(cond))
    for i in range(len(cond)):
        exp[i] = p1
        if cond[i] == 1 or cond[i] == 2:
            exp[i] = exp[i] + p2
        elif cond[i] == 3 or cond[i] == 4:
            exp[i] = exp[i] - p2
        if cond[i] == 1 or cond[i] == 3:
            exp[i] = exp[i] + p3
        elif cond[i] == 2 or cond[i] == 4:
            exp[i] = exp[i] - p3
    return exp

def zahn_func4(cond, p1, p2, p3, p4):
    exp = zeros(len(cond))
    if iter == 4:
        exp[i] = p1
        if cond[i] == 1 or cond[i] == 2:
            exp[i] = exp[i] + p2
        elif cond[i] == 3 or cond[i] == 4:
            exp[i] = exp[i] - p2
        if cond[i] == 1 or cond[i] == 3:
            exp[i] = exp[i] + p3
        elif cond[i] == 2 or cond[i] == 4:
            exp[i] = exp[i] - p3
        if cond[i] == 1 or cond[i] == 4:
            exp[i] = exp[i] + p4
        elif cond[i] == 2 or cond[i] == 3:
            exp[i] = exp[i] - p4
    return exp

print 'Unbalanced subjects create hidden interactions in the data\n'
print 'Data are from Zahn\'s article with male/female salaries and edu data\n'

data_2d = array([[24, 26, 25, 24, 27, 24, 27, 23, 15, 17, 20, 16, 25, 29, 27, 19, 18, 21, 20, 21, 22, 19],
                 [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4]])

# in the second row, 1 = F with deg, 2 = F w/o deg, 3 = M w/ deg, 4 = M w/o deg
data = data_2d[0,:]
cond = data_2d[1,:]
n_subj = len(data)

# start by calculating the standard deviations of the data, assuming each person comes from the same distribution
f_deg_inds = cond == 1
f_nodeg_inds = cond == 2
m_deg_inds = cond == 3
m_nodeg_inds = cond == 4

means = [average(data[f_deg_inds]), average(data[f_nodeg_inds]),
         average(data[m_deg_inds]), average(data[m_nodeg_inds])]
stddevs = [std(data[f_deg_inds]), std(data[f_nodeg_inds]),
           std(data[m_deg_inds]), std(data[m_nodeg_inds])]

expected = zeros(n_subj)
expected[f_deg_inds] = means[0]
expected[f_nodeg_inds] = means[1]
expected[m_deg_inds] = means[2]
expected[m_nodeg_inds] = means[3]

df = n_subj - 4 # fit 4 parameters for calculating SD with degrees of freedom
SD = sqrt(sum(((data - expected)**2)/df))

print 'SD = ' + str(SD) + ' with 4 means and ' + str(n_subj) + '-4 = ' + str(df) + ' degrees of freedom'

'''
f_inds = ((f_deg_inds)|(f_nodeg_inds))
m_inds = ((m_deg_inds)|(m_nodeg_inds))
nodeg_inds = ((f_nodeg_inds)|(m_nodeg_inds))
deg_inds = ((f_deg_inds)|(m_deg_inds))

design_matrix = ones((4, n_subj))
design_matrix[1, m_inds] = -1 # males
design_matrix[2, nodeg_inds] = -1 # low edu
design_matrix[3, :] = design_matrix[1,:]*design_matrix[2,:] # interactions
'''


popt, pcov = curve_fit(lambda cond, p1, p2: zahn_func4(cond, p1, p2, 0, 0), # function to fit to the data
                       cond, # x data
                       data, # y data
                       p0 = [1, 1], # number of parameters change relative to the zahn function number (ie, zahn2 = [1 1])
                       sigma = SD) # standard deviation to scale the Least Squares weights

#exp = zahn_func(cond,popt)



print popt
print pcov
