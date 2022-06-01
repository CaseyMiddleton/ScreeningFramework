from scipy.stats import gamma
import numpy as np
from functions import *
import os

def get_VL(t,Tlatent,Tpeak,Ypeak,Tclear,Ydetect):
    '''
    Computes the scalar value of a hinge function with control points:
        (Tlatent, Ydetect) : time spent between exposure and log Ydetect VL
        (Tpeak, Ypeak) : time between log 3 VL and peak infection, VL at peak
        (Tclear, Ydetect) : time between peak and crossing log Ydetect VL
    Returns zero whenever (1) t<T3 or (2) hinge(t) is negative.
    '''
    if t < Tlatent:
        return 0
    if t < Tpeak:
        VL = (t-Tlatent) * (Ypeak-Ydetect)/(Tpeak-Tlatent) + Ydetect
    else:
        VL = np.max([(t-Tpeak) * (Ydetect-Ypeak)/(Tclear-Tpeak) + Ypeak,0])
    return VL

def get_inf(t,Tlatent,Tpeak,Ypeak,Tclear,Ydetect):
    '''
    Computes the scalar value of infectiousness with control points:
        (Tlatent, Ydetect) : time spent between exposure and log Ydetect VL
        (Tpeak, Ypeak) : time between log 3 VL and peak infection, VL at peak
        (Tclear, Ydetect) : time between peak and crossing log Ydetect VL
    Returns infectious viral load at time t.
    '''
    inf = get_VL(t,Tlatent,Tpeak,Ypeak,Tclear,Ydetect) - Ydetect
    return np.max([inf,0])

def has_results(t,Tlatent,Tpeak,Ypeak,Tclear,Ydetect,sensitivity_threshold,failure_rate=0,delay=0,participation=1,compliance=1,frequency=1):
    '''
    Computes the probability of recieving results (including test + turnaround time) at time t given the following parameters:
        (Tlatent, Ydetect), (Tpeak, Ypeak), (Tclear, Ydetect) : VL params (see above)
        sensitivity_threshold : log10 VL threshold at which the virus can be detected by given test
        failure_rate : proportion of tests with VL >= sensitivity_threshold that return false negatives
        delay : time between testing and return of results (in days)
        participation : proportion of individuals who participate in testing, ever
        compliance : proportion of scheduled testa actually administered to individuals that participate in testing
        frequency : testing frequency (days). This parameter is unused in this function, but passed so as to keep all testing parameters in one dictionary.
    Returns scalar probability of recieving a positive test result at time t.
    '''
    VL = get_VL(t-delay,Tlatent,Tpeak,Ypeak,Tclear,Ydetect)
    if (t - delay >= 0) & (VL >= sensitivity_threshold):
        return (1 - failure_rate)*compliance
    return 0

def prob_test(t,frequency,**kwargs):
    '''
    Computes the probability of being tested at time t, assuming a uniform distribution of testing every Q days.
    Returns probability of being tested at time t.
    '''
    return 1/frequency

def create_testing_program(pathogen,**kwargs):
    param_func = pathogen + "_test_params"
    params = globals()[param_func]()
    for param_name in kwargs.keys():
        params[str(param_name)] = kwargs[param_name]
    return params

# Virus Parameter values ----------------------------------------------------------------

def COVID_wt_params_stoch():
    '''
    Defines a dictionary of stochastically drawn parameter values for COVID wt infectiousness.
    '''
    params = {}
    params['Tlatent'] = np.random.random()+2.5;
    params['Tpeak'] = gamma.rvs(1.5,loc=0.5);
    while params['Tpeak'] > 3:
        params['Tpeak'] = gamma.rvs(1.5,loc=0.5);
    params['Tpeak'] = params['Tpeak'] + params['Tlatent']
    params['Ypeak'] = np.random.uniform(7,11);
    params['Tclear'] = params['Tpeak'] + np.random.uniform(6.1,14.4);
    params['Ydetect'] = 3;
    #params['name'] = 'COVID_wt';
    return params

def COVID_wt_params_det():
    '''
    Defines a dictionary of deterministic parameter values for COVID wt infectiousness.
    '''
    params = {}
    params['Tlatent'] = 3;
    params['Tpeak'] = params['Tlatent'] + 2;
    params['Ypeak'] = 8;
    params['Tclear'] = params['Tpeak'] + 9;
    params['Ydetect'] = 3;
    #params['name'] = 'COVID_wt';
    return params

def COVID_delta_params_det():
    '''
    Defines a dictionary of deterministic parameter values for COVID delta infectiousness.
    '''
    params = {}
    params['Tlatent'] = 3
    params['Tpeak'] = params['Tlatent'] + 3
    params['Ypeak'] = 9.3
    params['Tclear'] = params['Tpeak'] + 7
    params['Ydetect'] = 3
    params['sensitivity_threshold'] = 3
    #params['name'] = 'COVID_delta';
    return params

def RSV_params_stoch():
    '''
    Defines a dictionary of stochastically drawn parameter values for RSV infectiousness.
    '''
    params = {}
    params['Tlatent'] = np.random.uniform(2,5);
    params['Tpeak'] = params['Tlatent'] + np.random.uniform(1,2) #4-7 days after challenge
    params['Ypeak'] = np.random.uniform(0.4,5.0);
    params['Tclear'] = params['Tpeak'] + np.random.uniform(1,3); # 8-10 days after challenge
    params['Ydetect'] = 0.9;
    #params['name'] = 'RSV';
    return params

def RSV_params_det():
    '''
    Defines a dictionary of deterministic parameter values for RSV infectiousness.
    '''
    params = {}
    params['Tlatent'] = 3;
    params['Tpeak'] = params['Tlatent'] + 2;
    params['Ypeak'] = 2.2;
    params['Tclear'] = params['Tpeak'] + 3;
    params['Ydetect'] = 0.9;
    #params['name'] = 'RSV';
    return params

def FLU_params_stoch():
    '''
    Defines a dictionary of stochastic parameter values for FLU infectiousness.
    '''
    params = {}
    params['Tlatent'] = np.random.uniform(1,4);
    params['Tpeak'] = params['Tlatent'] + np.random.uniform(1,3)
    params['Ypeak'] = np.random.uniform(4,9);
    params['Tclear'] = params['Tpeak'] + np.random.uniform(3,5);
    params['Ydetect'] = 0.7;
    #params['name'] = 'FLU';
    return params

def FLU_params_det():
    '''
    Defines a dictionary of deterministic parameter values for FLU infectiousness.
    '''
    params = {}
    params['Tlatent'] = 2;
    params['Tpeak'] = params['Tlatent'] + 2;
    params['Ypeak'] = 8;
    params['Tclear'] = params['Tpeak'] + 4;
    params['Ydetect'] = 0.7;
    #params['name'] = 'FLU';
    return params


# Testing parameters ---------------------------------------------------------
'''
Defines a dictionary of testing parameters
'''
def COVID_test_params():
    params = {}
    params['sensitivity_threshold'] = 3;
    params['frequency'] = 1;
    params['delay'] = 0;
    params['failure_rate'] = 0;
    return params

def RSV_test_params():
    params = {}
    params['sensitivity_threshold'] = 2; #0.5-3
    params['frequency'] = 1;
    params['delay'] = 0;
    params['failure_rate'] = 0;
    return params

def FLU_test_params():
    params = {}
    params['sensitivity_threshold'] = 3;
    params['frequency'] = 1;
    params['delay'] = 0;
    params['failure_rate'] = 0;
    return params

# Symptom parameters ---------------------------------------------------------
def COVID_is_feverish(t,dt,Tpeak):
	'''
	Computes the probability of detection of fever for COVID given the following parameters:
		t, dt : current time step and step size
		Tpeak : time of peak
	Returns scalar probability of detection at time t.
	'''
	fever_onset = Tpeak
	fever_offset = Tpeak + 2
	if t >= fever_onset and t < fever_offset:
		return 1
	else:
		return 0

def FLU_is_feverish(t,dt,Tpeak):
	'''
	Computes the probability of detection of fever for FLU given the following parameters:
		t, dt : current time step and step size
		Tpeak : time of peak
	Returns scalar probability of detection at time t.
	'''
	fever_onset = 2 # days after infection
	fever_offset = 5 # days after infection
	if t >= fever_onset and t < fever_offset:
		return 1
	else:
		return 0

def RSV_is_feverish(t,dt,Tpeak):
	'''
	Computes the probability of detection of fever for RSV given the following parameters:
		t, dt : current time step and step size
		Tpeak : time of peak
	Returns scalar probability of detection at time t.
	'''
	fever_onset = 6 # days after infection
	fever_offset = 8 # days after infection
	if t >= fever_onset and t < fever_offset:
		return 1
	else:
		return 0

def COVID_is_symptomatic(t,dt,Tpeak):
	'''
	Computes the probability of symptoms given the following parameters:
		t, dt : current time step and step size
		Tpeak : time of peak
	Returns scalar probability of symptoms at time t.
	'''
	symp_onset = Tpeak # onset at peak
	symp_offset = Tpeak + 1.5 # last 1.5 days
	if t >= symp_onset and t < symp_offset:
		return 1
	else:
		return 0

# Test Parameter values ----------------------------------------------------------------

def example_params_det():
    '''
    Defines a dictionary of deterministic parameter values for examplevirus infectiousness.
    '''
    params = {}
    params['Tlatent'] = 0;
    params['Tpeak'] = params['Tlatent'] + 5;
    params['Ypeak'] = 10;
    params['Tclear'] = params['Tpeak'] + 5;
    params['Ydetect'] = 0;
    #params['name'] = 'Examplevirus';
    return params

def example_test_params():
    params = {}
    params['sensitivity_threshold'] = 0; # perfect testing scenario
    params['frequency'] = 1;
    params['delay'] = 0;
    params['failure_rate'] = 0;
    return params

def example_is_symptomatic(t,Tpeak,**kwargs):
    '''
    Computes the probability of symptoms given the following parameters:
    t : current time step
    Tpeak : time of peak
    Returns scalar probability of symptoms at time t.
    '''
    symp_onset = 3 # onset at day 3
    symp_offset = 6 # last 3 days
    if t >= symp_onset and t < symp_offset:
        return 1
    else:
        return 0
