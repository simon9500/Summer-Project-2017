import sys
import numpy as np
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# This program takes in input files of diffrate/time and the integrated rate/energy threshold
# and samples multiple neutrino start times, which are saved in a histogram.

# Read names of output file and input file from command line
if len(sys.argv) != 3:
    print "Wrong number of arguments."
    print "Usage: " + sys.argv[0] + " <(diffrate) input file> <(energy threshold) input file>"
    quit()
else:
    infileName1 = sys.argv[1]
    infileName2 = sys.argv[2]

# Open input file for reading
infile1 = open(infileName1, "r")
infile2 = open(infileName2, "r")

# Initialise arrays
time_val_list = []
diffrate_list = []
intrate_list = []
enthreshold_list = []

# Scaling factors for total events
energy_threshold = 0. # Energy threshold of the detector in keV
detector_mass = 7. # Detector mass in tonnes (7 for LZ)
SN_dist_unit = 3.0857e+20 # 10kpc in m - all readings are relative to this distance and scale as 1/r^2
SN_dist_actual = 3.0857e+20 # 1kpc in m is 3.0857e+19
SN_dist_event_scaling = 1./(SN_dist_actual/SN_dist_unit)**2 # Scaling factor for events by distance SN is from earth

# Fill arrays
for line in infile1:
    tokens = line.split()
    time_val_list.append(float(tokens[0]))
    diffrate_list.append(float(tokens[1]))

for line in infile2:
    tokens = line.split(",")
    enthreshold_list.append(float(tokens[0]))
    intrate_list.append(float(tokens[1]))

# Close files
infile1.close()
infile2.close()

# Convert lists to numpy arrays
time_val_list = np.array(time_val_list)
diffrate_list = np.array(diffrate_list)
intrate_list = np.array(intrate_list)
enthreshold_list = np.array(enthreshold_list)

# Fit threshold data with an exponential
def exp_func(x, a, b, c):
    return a * np.exp(-b * x) + c

popt, pcov = curve_fit(exp_func, enthreshold_list, intrate_list)

total_events = np.sum(diffrate_list)/1000. # Total integral of the file
# Scale by SN distance, detector mass and the threshold energy of the detector
# Scaling by mass is done by multiplying by the mass in tonnes of the detector
# The total events scale off as 1/r^2, where r is the SN distance
# The scaling by energy threshold is done by finding the value of the fitted function at the energy
# threshold of the detector and scaling by this value over the value of the function at 0 energy threshold
total_events_scaled = SN_dist_event_scaling*detector_mass*exp_func(energy_threshold, *popt)
print "Total integral of file = %f" %total_events
print "Total number of events = %d" %int(total_events_scaled)
start_time_samples = [] # Array to store all start times
num_start_time = 1000 # Number of start times to add to histogram
num_events = 1000 #int(total_events_scaled)
outputfile_name = 'histstarttime_' + str(num_events) + 'events_' + str(num_start_time) + 'ent.png'

# Running integral (over full 10s range) array - to save computation in the loop
running_integral_10s = np.empty(len(diffrate_list))
integral_sofar = 0.
for m in range(len(diffrate_list)):
    integral_sofar += diffrate_list[m]/1000.
    running_integral_10s[m] = integral_sofar

for k in range(num_start_time):
    start_time = 20. # Set start time to large number for the purpose of calculating a minimum
    for i in range(num_events):
        random_sample = random.random()*total_events
        for j in range(len(time_val_list)):
            if running_integral_10s[j] >= random_sample:
                # Interpolate between the interval endpoints
                # Current integral is the integral up to the start of the interval
                # tvals and dr_interp are the interpolated data points over the 1ms interval
                if j == 0:
                    tvals = np.linspace(0., time_val_list[j], 101)
                    dr_interp = np.interp(tvals, [0., time_val_list[j]],
                                         [diffrate_list[j], diffrate_list[j]])
                    current_integral = 0.
                else:
                    tvals = np.linspace(time_val_list[j - 1], time_val_list[j], 101)
                    dr_interp = np.interp(tvals, [time_val_list[j - 1], time_val_list[j]],
                                         [diffrate_list[j - 1], diffrate_list[j]])
                    current_integral = running_integral_10s[j - 1]
                # Run through integral and find the time value of the random sample to a precision of 0.01ms
                for n in range(1,len(tvals)):
                    if current_integral >= random_sample:
                        first_encounter = tvals[n]
                        break
                    current_integral += dr_interp[n] / 100000.
                break
        if first_encounter < start_time:
            start_time = first_encounter
    print "%d: Start time = %fms" % (k+1, start_time * 1000)
    start_time_samples.append(start_time)

start_time_samples = np.array(start_time_samples)

"""
plt.plot(enthreshold_list, intrate_list, 'b-', label='data')
plt.plot(enthreshold_list, exp_func(enthreshold_list, *popt), 'r-', label='fit')
plt.xlabel('Energy threshold')
plt.ylabel('Integrated rate')
plt.axis([0, 10, 0, 30])
plt.legend()
plt.savefig('thresholdfit.png', bbox_inches='tight')
plt.show()

n1, bins1, patches1 = plt.hist(time_samples, bins = 100, range = (0,10), normed = True)

plt.xlabel('Time')
plt.ylabel('Differential rate')
plt.title('Histogram of Differential rate')
plt.axis([0, 10, 0, 1])
plt.savefig('histdiffrate.png', bbox_inches='tight')
plt.show()
"""

n2, bins2, patches2 = plt.hist(start_time_samples, bins = 1000, range = (0,0.1))

plt.xlabel('Start time')
plt.ylabel('Num of Occurences')
plt.title('Histogram of Start times')
plt.axis([0, 0.1, 0, num_start_time])
plt.savefig(outputfile_name, bbox_inches='tight')


