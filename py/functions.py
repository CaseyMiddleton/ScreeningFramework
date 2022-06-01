import numpy as np

def rect_rule(y,t):
    '''
    This function adds up the area of rectangles
    defined by the points in vectors (t,y)

    Returns a scalar
    '''
    dt = t[1]-t[0]
    area = np.sum(y)*dt
    return area

def cumulative_rect_rule(y,t):
    '''
    This function adds up the area of rectangles
    defined by the points in (y,t)

    Returns a vector of the total area so far as we
    integrate from left to right.
    '''
    dt = t[1]-t[0]
    areas = np.cumsum(y)*dt
    return areas


def R_reduction(infectiousness,detectability,testing,t,**kwargs):
    '''
    This function calculates the reduction in transmission as
    a result of:
        infectiousness: vector of infectious viral load
        detectability: pr(detection | tested and true positive) at time t
        testing: pr(tested) at time t

    Returns a scalar in [0,1] of reduction in R
    '''
    total_inf = get_total_inf(infectiousness,t)[-1]
    resid_inf = get_resid_inf(infectiousness,detectability,testing,t,**kwargs)[-1]
    reduction = resid_inf/total_inf
    return 1-reduction

def R_reduction_time_series(infectiousness,detectability,testing,t,**kwargs):
    '''
    This function calculates the reduction in transmission as
    a result of:
        infectiousness: vector of infectious viral load
        detectability: pr(detection | tested and true positive) at time t
        testing: pr(tested) at time t

    Returns a scalar in [0,1] of reduction in R
    '''
    total_inf = get_total_inf(infectiousness,t)
    resid_inf = get_resid_inf(infectiousness,detectability,testing,t,**kwargs)
    reduction = resid_inf/total_inf
    return 1-reduction

def get_total_inf(my_inf,t):
    '''
    This function calculates infectiousness using the rectangle rule as a function of:
        infectiousness: vector of infectious viral load

    Returns a vector of the cumulative area under the infectious curve at each time point in t
    '''
    area = cumulative_rect_rule(my_inf,t)
    return area

def get_resid_inf(my_inf,my_det,my_tes,t,participation=1):
    '''
    This function calculates the infectiousness using the rectangle rule as
    a function of:
        infectiousness: vector of infectious viral load
        detectability: pr(detection | tested and true positive) at time t
        testing: pr(tested) at time t
        t: vector of time values
        participation : proportion of individuals who participate in testing
    Returns a vector for the likely cumulative infectiousness at each time point in t
    '''
    # begin by calculating cumulative probability of detection at time t (inner intergral)
    prob_det_at_t = np.multiply(my_det,my_tes)
    cumulative_prob_det = np.minimum(cumulative_rect_rule(prob_det_at_t,t), np.ones(len(my_det)))

    # now calculate the probability that the individual is *not* detected at time t (1 - inner integral)
    cumulative_prob_not_det = [1 - prob for prob in cumulative_prob_det]

    # now calculate residual infectiousness (outer integral)
    integrand = np.multiply(my_inf,cumulative_prob_not_det)
    tested_inf = cumulative_rect_rule(integrand,t)

    # incorporate participation
    untested_inf = get_total_inf(my_inf,t)
    resid_inf = (1-participation)*untested_inf + participation*tested_inf

    return resid_inf


## Self isolation due to symptoms

def R_reduction_with_symp(infectiousness,detectability,testing,symptoms,t,prob_symp,prob_isolate,**kwargs):
    '''
    This function calculates the reduction in transmission as
    a result of:
        infectiousness: vector of infectious viral load
        detectability: pr(detection | tested and true positive) at time t
        testing: pr(tested) at time t
        symptoms: pr(symptomatic) at time t for individuals showing symptoms during the course of infection
        prob_symp : p(symptomatic | infected)
        prob_isolate : p(isolate | symptomatic)
    Returns a scalar in [0,1] of reduction in R
    '''
    total_inf = get_total_inf_with_symp(infectiousness,symptoms,t,prob_symp,prob_isolate)[-1]
    resid_inf = get_resid_inf_with_symp(infectiousness,detectability,testing,symptoms,t,prob_symp,prob_isolate,**kwargs)[-1]
    reduction = resid_inf/total_inf
    return 1-reduction

def R_reduction_with_symp_time_series(infectiousness,detectability,testing,symptoms,t,prob_symp,prob_isolate,**kwargs):
    '''
    This function calculates the reduction in transmission as
    a result of:
        infectiousness: vector of infectious viral load
        detectability: pr(detection | tested and true positive) at time t
        testing: pr(tested) at time t
        symptoms: pr(symptomatic) at time t for individuals showing symptoms during the course of infection
        prob_symp : p(symptomatic | infected)
        prob_isolate : p(isolate | symptomatic)
    Returns a scalar in [0,1] of reduction in R
    '''
    total_inf = get_total_inf_with_symp(infectiousness,symptoms,t,prob_symp,prob_isolate,**kwargs)
    resid_inf = get_resid_inf_with_symp(infectiousness,detectability,testing,symptoms,t,prob_symp,prob_isolate,**kwargs)
    reduction = resid_inf/total_inf
    return 1-reduction

def get_total_inf_with_symp(my_inf,my_symp,t,prob_symp,prob_isolate):
    '''
    This function calculates the total infectiousness assuming some individuals self isolate with symptoms
        my_symp : returns 0 when asymptomatic, 1 when symptomatic throughout course of infection
        prob_symp : p(symptomatic | infected)
        prob_isolate : p(isolate | symptomatic)
    '''
    asymp_inf = my_inf * (1-prob_symp)
    symp_inf = my_inf * prob_symp * (1 - prob_isolate*my_symp)
    area = cumulative_rect_rule(asymp_inf + symp_inf,t)
    return area

def get_resid_inf_with_symp(my_inf,my_det,my_tes,my_symp,t,prob_symp,prob_isolate,participation=1):
    '''
    This function calculates the residual infectiousness assuming some individuals self isolate with symptoms
    '''
    my_max_det = participation*np.ones(len(my_det))
    # begin by calculating cumulative probability of detection at time t (inner intergral)
    prob_det_at_t = np.multiply(my_det,my_tes)
    cumulative_prob_det = np.minimum(cumulative_rect_rule(prob_det_at_t,t), my_max_det)

    # now calculate the probability that the individual is *not* detected at time t (1 - inner integral)
    cumulative_prob_not_det = [1 - (1-prob_symp*prob_isolate)*prob for prob in cumulative_prob_det]

    # now calculate residual infectiousness (outer integral)
    integrand = np.multiply(my_inf,cumulative_prob_not_det)
    resid_inf = cumulative_rect_rule(integrand,t)
    return resid_inf

## Multiple screening strategies

def R_reduction_multiscreen(infectiousness,detectability,testing,symptoms,t,prob_symp,prob_isolate,**kwargs):
    '''
    This function calculates the reduction in transmission as
    a result of:
        infectiousness: vector of infectious viral load
        detectability: list of [pr(detection | tested and true positive) at time t] for each test
        testing: list of [pr(tested) at time t] for each test
        symptoms: pr(symptomatic) at time t for individuals showing symptoms during the course of infection
        prob_symp : p(symptomatic | infected)
        prob_isolate : p(isolate | symptomatic)
    Returns a scalar in [0,1] of reduction in R
    '''
    total_inf = get_total_inf_with_symp(infectiousness,symptoms,t,prob_symp,prob_isolate)[-1]
    resid_inf = get_resid_inf_multiscreen(infectiousness,detectability,testing,symptoms,t,prob_symp,prob_isolate,**kwargs)[-1]
    reduction = resid_inf/total_inf
    return 1-reduction

def R_reduction_multiscreen_time_series(infectiousness,detectability,testing,symptoms,t,prob_symp,prob_isolate,**kwargs):
    '''
    This function calculates the reduction in transmission as
    a result of:
        infectiousness: vector of infectious viral load
        detectability: list of [pr(detection | tested and true positive) at time t] for each test
        testing: list of [pr(tested) at time t] for each test
        symptoms: pr(symptomatic) at time t for individuals showing symptoms during the course of infection
        prob_symp : p(symptomatic | infected)
        prob_isolate : p(isolate | symptomatic)
    Returns a scalar in [0,1] of reduction in R
    '''
    total_inf = get_total_inf_with_symp(infectiousness,symptoms,t,prob_symp,prob_isolate)
    resid_inf = get_resid_inf_multiscreen(infectiousness,detectability,testing,symptoms,t,prob_symp,prob_isolate,**kwargs)
    reduction = resid_inf/total_inf
    return 1-reduction

def get_resid_inf_multiscreen(my_inf,my_det,my_tes,my_symp,t,prob_symp=0,prob_isolate=0,participation=1):
    '''
    This function calculates the residual infectiousness assuming some individuals self isolate with symptoms
    '''
    # ensure proper number of arguments are passed
    assert(len(my_det) == len(my_tes))

    my_max_det = participation*np.ones(len(my_det[0]))
    # First, calculate probability of detection (inner integral)
    summed_prob_detection = np.zeros(len(t))
    summed_prob_double_overlap = np.zeros(len(t))
    prob_total_overlap = np.ones(len(t))

    # Probability of detection by any individual test (additive)
    for ii in range(0,len(my_det)):
        prob_det_at_t = np.multiply(my_det[ii],my_tes[ii])
        cumulative_prob_det = np.minimum(cumulative_rect_rule(prob_det_at_t,t), my_max_det)
        summed_prob_detection = np.minimum(np.add(summed_prob_detection,cumulative_prob_det), np.ones(len(my_det[0])))
        # Probability of being detected by two tests (subtractive)
        if ii < len(my_det) - 1:
            prob_det_at_t_test2 = np.multiply(my_det[ii+1],my_tes[ii+1])
            cumulative_prob_det_test2 = np.minimum(cumulative_rect_rule(prob_det_at_t_test2,t), my_max_det)
            prob_double_overlap = np.multiply(cumulative_prob_det,cumulative_prob_det_test2)
            summed_prob_double_overlap = np.minimum(np.add(summed_prob_double_overlap,prob_double_overlap), np.ones(len(my_det[0])))
        # Probability of all tests overlapping (additive)
        prob_total_overlap = np.multiply(prob_total_overlap,cumulative_prob_det)

    # now calculate the probability that the individual is *not* detected at time t (1 - inner integral)
    probs = np.subtract(summed_prob_detection,summed_prob_double_overlap)
    probs = np.add(probs,np.array(prob_total_overlap)*(len(my_det)-1))

    # print(probs)
    cumulative_prob_not_det = [1 - (1-prob_symp*prob_isolate)*prob for prob in probs]

    # now calculate residual infectiousness (outer integral)
    integrand = np.multiply(my_inf,cumulative_prob_not_det)
    resid_inf = cumulative_rect_rule(integrand,t)
    return resid_inf


## Serial Interval stuff

def calc_serial_int(exp_inf_notesting, exp_inf,t):
    '''
    This function calculates the expected probability of transmission (serial interval) at time t given
        exp_inf_notesting : the expected infectiousness at time t without testing control measures
        exp_inf_notesting : the expected infectiousness at time t with or without testing control measures
    Returns a vector of the likelihood of transmission at time t
    '''
    # convert expected outward infectiousness to probability function
    prob_transmits = exp_inf / max(exp_inf_notesting)
    exp_value = []
    for i in range(0,len(t)-1):
        exp_value.append(rect_rule(t[i:i+2] * prob_transmits[i:i+2],t[i:i+2]))
    return exp_value
