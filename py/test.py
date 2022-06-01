import numpy as np
from functions import *
from parameters import *
from sample_simulations import *

def test1():
    '''
    rect_rule calculates the area under a line
    Test that area under line y = x+1 on interval [0,1] = 1.5
    '''
    # define a line
    t = np.linspace(0,1,1001)
    m = 1
    b = 1
    y = m*t+b
    true_area = 1.5
    area = rect_rule(y,t)
    residual = np.abs(true_area - area)
    assert(residual < 0.01)

test1()

def test2():
    '''
    cumulative_rect_rule returns vector of the total area so far as we integrate from left to right.
    Test that area under a horizontal line equals the sum of previous rectangles.
    '''
    t = np.linspace(0,1,101)
    y = np.ones(101)
    true_areas = np.arange(1,102)*0.01
    areas = cumulative_rect_rule(y,t)
    assert((areas==true_areas).all())

test2()

def test3():
    '''
    R_reduction returns the a scalar value of the reduction in R from testing.
    Test that approximately no testing gives approximately no reduction in R.
    '''
    t_vals = np.linspace(0,10,1001)
    predicted_reduction = 0

    inf_params = example_params_det()
    test_params = create_testing_program("example",frequency=1000000)

    inf = np.array([get_inf(t,**inf_params) for t in t_vals])
    is_detectable = np.array([has_results(t,**inf_params,**test_params) for t in t_vals])
    pr_test = np.array([prob_test(t,**test_params) for t in t_vals])

    actual_reduction = R_reduction(inf,is_detectable,pr_test,t_vals)
    residual = np.abs(actual_reduction - predicted_reduction)
    assert(residual < 0.001)

test3()

def test4():
    '''
    Test that approximately constant testing gives approximately 100% reduction in R.
    '''
    t_vals = np.linspace(0,10,1001)
    predicted_reduction = 1

    inf_params = example_params_det()
    test_params = create_testing_program("example",frequency=0.0001)

    inf = np.array([get_inf(t,**inf_params) for t in t_vals])
    is_detectable = np.array([has_results(t,**inf_params,**test_params) for t in t_vals])
    pr_test = np.array([prob_test(t,**test_params) for t in t_vals])

    actual_reduction = R_reduction(inf,is_detectable,pr_test,t_vals)
    residual = np.abs(actual_reduction - predicted_reduction)
    assert(residual < 0.001)

test4()

def test5():
    '''
    Test that approximately constant testing with a 5 day turnaround time gives approximately 50% reduction in R.
    This is due to the symmetry of the viral growth and decay around day 5 for examplevirus.
    '''
    t_vals = np.linspace(0,10,1001)
    predicted_reduction = 0.5

    inf_params = example_params_det()
    test_params = create_testing_program("example",frequency=0.0001,delay=5)

    inf = np.array([get_inf(t,**inf_params) for t in t_vals])
    is_detectable = np.array([has_results(t,**inf_params,**test_params) for t in t_vals])
    pr_test = np.array([prob_test(t,**test_params) for t in t_vals])

    actual_reduction = R_reduction(inf,is_detectable,pr_test,t_vals)
    residual = np.abs(actual_reduction - predicted_reduction)
    assert(residual < 0.001)

test5()

def test6():
    '''
    Test that near constant testing with the limit of detection ~= peak viral load gives 50% reduction in R
    '''
    t_vals = np.linspace(0,10,1001)
    predicted_reduction = 0.5

    inf_params = example_params_det()
    test_params = create_testing_program("example",frequency=.001,sensitivity_threshold=inf_params['Ypeak'] - 0.01)

    inf = np.array([get_inf(t,**inf_params) for t in t_vals])
    is_detectable = np.array([has_results(t,**inf_params,**test_params) for t in t_vals])
    pr_test = np.array([prob_test(t,**test_params) for t in t_vals])

    actual_reduction = R_reduction(inf,is_detectable,pr_test,t_vals)
    residual = np.abs(actual_reduction - predicted_reduction)
    assert(residual < 0.001)

test6()

def test7():
    '''
    Test that near constant testing with a perfect test at 50% participation leads to ~50% reduction in R
    '''
    t_vals = np.linspace(0,10,1001)
    predicted_reduction = 0.5

    inf_params = example_params_det()
    test_params = create_testing_program("example",frequency=.0001)

    reductions = []
    for i in range(0,100):
        # 50% of individuals have a 100% failure rate, representing non-participation
        if i%2 ==0:
            test_params['failure_rate'] = 1
        else:
            test_params['failure_rate'] = 0
        inf = np.array([get_inf(t,**inf_params) for t in t_vals])
        is_detectable = np.array([has_results(t,**inf_params,**test_params) for t in t_vals])
        pr_test = np.array([prob_test(t,**test_params) for t in t_vals])

        reductions.append(R_reduction(inf,is_detectable,pr_test,t_vals))

    avg_reduction = np.mean(reductions)
    residual = np.abs(avg_reduction - predicted_reduction)
    assert(residual < 0.001)

test7()

def test8():
    '''
    Test that simulations for example virus with detection at day 5 produces 50% reduction in infectiousness
    '''
    # Time series data parameters
    start = 0
    stop = 10
    step = 0.0001
    numels = int((stop - start)/step + 1)
    t_vals = np.linspace(start,stop,numels)
    dt = t_vals[1]-t_vals[0]

    VL = generate_example_trajectory()
    tested = np.zeros(len(VL))
    tested[int(5/dt)] = 1 # tested at day 5
    reduction = individual_reduction(VL,tested)
    residual = np.abs(reduction - 0.5)
    assert(residual < 0.001)

test8()

def test9():
    '''
    Test that under a 5-day testing schedule, example virus has at least a
        50% reduction in infectiousness for scheduled testing
    '''
    n = 10
    freq = 5
    reductions = population_reduction_scheduled(n,freq)
    for r in reductions:
        assert(r >= 0.5)

test9()

def test10():
    ''' test multiple screening strategies '''
    # if one strategy is fully effective and the other is never effective,
        # predict 100% reduction
    # get trajectory
    inf_params = example_params_det()
    vl = np.array([get_VL(t,**inf_params) for t in t_vals])
    inf = np.array([get_inf(t,**inf_params) for t in t_vals])
    symp = np.zeros(len(vl)) # no symptoms
    prop_symp = 0; prop_isolate = 0;

    test_freq = 0.001

    # test 1 parameters : viral detection test
    test1_LOD = 0
    test1_params = create_testing_program("example",sensitivity_threshold=test1_LOD,frequency=test_freq)
    test1_det = np.array([has_results(t,**inf_params,**test1_params) for t in t_vals])
    tst1 = np.array([prob_test(t,**test1_params) for t in t_vals])
    # test 2 parameters : symptom screening
    test2_LOD = np.max(vl) + 1
    test2_params = create_testing_program("example",sensitivity_threshold=test2_LOD,frequency=test_freq)
    test2_det = np.array([has_results(t,**inf_params,**test2_params) for t in t_vals])
    tst2 = np.array([prob_test(t,**test2_params) for t in t_vals])

    # get R_reduction
    det = [test1_det,test2_det]
    tst = [tst1,tst2]
    reduction = R_reduction_multiscreen(inf,det,tst,symp,t_vals,prop_symp,prop_isolate)
    expected_reduction = 1
    error = np.abs(reduction - expected_reduction)
    assert(error <= 0.01)

    # Reduction should never be greater than 1
    # test 2 parameters : symptom screening
    test2_LOD = 0
    test2_params = create_testing_program("example",sensitivity_threshold=test2_LOD,frequency=test_freq)
    test2_det = np.array([has_results(t,**inf_params,**test2_params) for t in t_vals])
    tst2 = np.array([prob_test(t,**test2_params) for t in t_vals])

    # get R_reduction
    det = [test1_det,test2_det]
    tst = [tst1,tst2]
    reduction = R_reduction_multiscreen(inf,det,tst,symp,t_vals,prop_symp,prop_isolate)
    assert(reduction <= 1)

test10()

def test11():
    ''' test participation '''
    # if participation = 0, expect no reduction in infectiousness
    # get trajectory
    inf_params = example_params_det()
    vl = np.array([get_VL(t,**inf_params) for t in t_vals])
    inf = np.array([get_inf(t,**inf_params) for t in t_vals])

    test_freq = 0.001

    # test parameters : perfect viral detection test
    test_LOD = 0
    test_params = create_testing_program("example",sensitivity_threshold=test_LOD,frequency=test_freq)
    det = np.array([has_results(t,**inf_params,**test_params) for t in t_vals])
    tst = np.array([prob_test(t,**test_params) for t in t_vals])

    reduction = R_reduction(inf,det,tst,t_vals,participation=0)
    expected = 0
    error = np.abs(reduction - expected)
    assert(error <= 0.01)

    # if using a perfect test with 100% participation, expect 100% reduction
    reduction = R_reduction(inf,det,tst,t_vals,participation=1)
    expected = 1
    error = np.abs(reduction - expected)
    assert(error <= 0.01)

    # if using a perfect test wth participatin p, expect p reduction (ex p = 50%)
    reduction = R_reduction(inf,det,tst,t_vals,participation=0.5)
    expected = 0.5
    error = np.abs(reduction - expected)
    assert(error <= 0.01)

test11()
