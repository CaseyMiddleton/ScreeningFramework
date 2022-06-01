import numpy as np
from functions import *
from parameters import *
import matplotlib as plt
from prettyplotlib import *
import pandas as pd

''' Global Parameters '''
# File paths
fig_path = r'/Users/caseymiddleton/env/ActiveProjects/TestingFramework/figs/'
data_path = r'/Users/caseymiddleton/env/ActiveProjects/TestingFramework/TestingFrameworkRepo/data/'
# Time series data parameters
start = 0
stop = 10
step = 0.0001
numels = int((stop - start)/step + 1)
t_vals = np.linspace(start,stop,numels)
dt = t_vals[1]-t_vals[0]
# Sensitivity threshold
sensitivity_threshold = 0

def create_plot(freq):
    ''' Create a plot that compares testing every five days using a scheduled vs random probability '''
    # Read in data
    random_file = data_path + "simulation_reduction_random_testing.csv"
    random_pop_reduction = pd.read_csv(random_file)
    scheduled_file = data_path + "simulation_reduction_scheduled_testing.csv"
    scheduled_pop_reduction = pd.read_csv(scheduled_file)

    avg_random = np.mean(random_pop_reduction)
    avg_scheduled = np.mean(scheduled_pop_reduction)
    x_random = ["random" for i in range(0,len(random_pop_reduction))]
    x_scheduled = ["scheduled" for i in range(0,len(scheduled_pop_reduction))]

    # include reduction predicted by model
    mod_prediction = population_reduction_model(freq)

    fig, ax1 = plt.subplots()

    ax1.scatter(x_random, random_pop_reduction, color=DARK_COLOR, alpha=0.2, label="Random Testing")
    ax1.scatter("random", avg_random, color=ACCENT_COLOR_1)
    ax1.scatter(x_scheduled, scheduled_pop_reduction, color=DARK_COLOR, alpha=0.2, label="Scheduled Testing")
    ax1.scatter("scheduled", avg_scheduled, color=ACCENT_COLOR_1)
    ax1.plot(["random","scheduled"],[mod_prediction,mod_prediction], color=ACCENT_COLOR_1, label="Model Prediction", linestyle="dashed")
    ax1.set_ylim(0,1.1)

    ax1.set_xlabel("Testing Scenario")
    ax1.set_ylabel("Reduction in Transmission Potential")
    hide_right_top_axis(ax1)
    finalize(ax1)

    plt.tight_layout()
    plt.show()

def generate_data(freq):
    ''' Create a plot that compares testing every five days using a schedule vs random probability '''
    n = 100

    random_pop_reduction = population_reduction_random(n,freq)
    scheduled_pop_reduction = population_reduction_scheduled(n,freq)

    fname = data_path + "simulation_reduction_random_testing" + ".csv"
    np.savetxt(fname, random_pop_reduction, delimiter=',')
    fname = data_path + "simulation_reduction_scheduled_testing" + ".csv"
    np.savetxt(fname, scheduled_pop_reduction, delimiter=',')

def population_reduction_model(freq):
    VL = generate_example_trajectory()
    p = example_params_det()
    p['sensitivity_threshold'] = 0
    p['Q'] = freq
    det = np.array([has_results(t,p['Tlatent'],p['Tpeak'],p['Ypeak'],p['Tclear'],p['Ydetect'],p['sensitivity_threshold'],p['F'],p['delay']) for t in t_vals])
    tst = np.array([prob_test(t,p['Q']) for t in t_vals])
    return R_reduction(VL,det,tst,t_vals)

def population_reduction_random(n,freq):
    reduction = [ individual_reduction_random(freq) for i in range(0,n) ]
    return reduction

def population_reduction_scheduled(n,freq):
    reduction = [ individual_reduction_scheduled(freq) for i in range(0,n) ]
    return reduction

def individual_reduction_random(freq):
    ''' Returns the reduction in outward infectiousness for a random testing schedule '''
    VL = generate_example_trajectory()
    tested = generate_testing_random(freq)
    return individual_reduction(VL,tested)

def individual_reduction_scheduled(freq):
    ''' Returns the reduction in outward infectiousness for a scheduled testing frequency, freq '''
    VL = generate_example_trajectory()
    tested = generate_testing_schedule(freq)
    return individual_reduction(VL,tested)

def individual_reduction(VL,schedule):
    ''' Returns the reduction in outward infectiousness for a any testing schedule and VL '''
    time_detected = calculate_time_detected(VL,schedule)
    inf_without_testing = get_total_inf(VL,t_vals)[-1]
    inf_with_testing = get_total_inf(VL[0:time_detected],t_vals[0:time_detected])[-1]
    return ((inf_without_testing - inf_with_testing)/inf_without_testing)

def calculate_time_detected(VL,test_schedule):
    ''' Returns the day of the first positive test result '''
    for i in range(0,len(VL)):
        if test_schedule[i] == 1:
            if VL[i] > sensitivity_threshold:
                return i
    # if never detected, return end of trajectory
    return len(VL)

def generate_example_trajectory():
    ''' Generates a det trajectory of example viral load '''
    p = example_params_det()
    inf = np.array([get_inf(t,**p) for t in t_vals])
    return inf

def generate_testing_schedule(freq):
    ''' Generates a random testing time for an individual according to testing frequency freq '''
    first_test = np.random.randint(low=0,high=(freq/dt+dt)) # index of first testing day in t_vals
    test_schedule = []
    test_day = first_test
    for i in range (0,len(t_vals)):
        if i == test_day:
            test_schedule.append(1)
            test_day += freq/dt
        else:
            test_schedule.append(0)
    return test_schedule

def generate_testing_random(freq):
    ''' Generates random testing for an individual according to a daily probability of being tested  '''
    prob_test = dt/freq
    test_schedule = [ (np.random.random() <= prob_test) for t in t_vals ]
    return test_schedule

def main():
    freq = 5
    generate_data(freq)
    create_plot(freq)

if __name__ == '__main__':
    main()
