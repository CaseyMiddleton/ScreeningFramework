import numpy as np
import importlib
import pandas as pd

# change directory
import sys
sys.path.insert(1,r'/Users/caseymiddleton/Documents/ActiveProjects/TestingTheory/TestingFrameworkRepo/py')

from functions import *
from parameters import *

data_path = r'/Users/caseymiddleton/Documents/ActiveProjects/TestingTheory/TestingFrameworkRepo/data'

''' Time series data parameters '''
start = 0
stop = 15
step = 0.0001
numels = int((stop - start)/step + 1)
t_vals = np.linspace(start,stop,numels)
dt = t_vals[1]-t_vals[0]

def get_reduction(L,Q,pathogen):
    ''' Calculates the predicted reduction in infectiousness
            for a given set of parameters:
        L : Limit of detection
        Q : Days between tests
        pathogen : COVID_wt, COVID_delta, RSV, or FLU
    '''
    inf_params,test_params = get_params(L,Q,pathogen)
    inf = np.array([get_inf(t,**inf_params) for t in t_vals])
    is_detectable = np.array([has_results(t,**inf_params,**test_params) for t in t_vals])
    pr_test = np.array([prob_test(t,**test_params) for t in t_vals])

    reduction = R_reduction(inf,is_detectable,pr_test,t_vals)
    return reduction

def get_params(L,Q,pathogen):
    '''
    get_params() defines the needed parameter values for each infectiousness and testing function.
        L : Limit of detection
        Q : Days between tests
        pathogen : COVID_wt, COVID_delta, RSV, FLU, or example
    '''
    inf_params = globals()[pathogen+"_params_det"]()
    test_params = create_testing_program(pathogen,sensitivity_threshold=L,frequency=D)
    return inf_params,test_params


# ------------------------------------------------------------------------------------------------------------
# ---------------------------------- Examplevirus figures ----------------------------------------------------
# ------------------------------------------------------------------------------------------------------------








# ------------------------------------------------------------------------------------------------------------
# ----------------------------------- Manuscript figures -----------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

def line_Q_LOD(pathogen):
    '''
    Line plot comparing sensitivity (L) to testing frequency (Q) for any pathogen
         pathogen : COVID_wt, COVID_delta, RSV, or FLU
    '''
    Qvals = np.arange(1,31,1);
    Lvals = [3,5,6,7]
    lims = {}
    for L in Lvals:
        reds = []
        for Q in Qvals:
            # calculate reduction values and store in list
            reduction = get_reduction(L,Q,pathogen)
            reds.append(reduction)
        # store list of reduction values in dictionary for L
        key_val = "reduction_LOD"+str(L)
        lims[key_val] = reds
    #save data
    df = pd.DataFrame.from_dict(lims)
    print(df)
    filename = data_path + "line_Q_LOD_" + pathogen + ".csv"
    df.to_csv(filename, index = False, header=True)

def generate_stochastic_trajectories(pathogen,n):
    '''
    Generates n stochastic trajectories for the pathogen specified
    '''
    inf_traj_draws = [];
    for i in range(1,n+1):
        p = globals()[pathogen+"_params_stoch"]()
        inf_traj_draws.append(np.array([get_inf(t,**p) for t in t_vals]))
        #is_detectable = np.array([results(t,p['Tlatent'],p['Tpeak'],p['Ypeak'],p['Tclear'],p['Ydetect'],p['sensitivity_threshold']) for t in t_vals])
    # write data to file
    inf_traj_draws = np.array(inf_traj_draws)
    fname = data_path + "inf_trajectories_" + pathogen + ".csv"
    np.savetxt(fname, inf_traj_draws, delimiter=',')
