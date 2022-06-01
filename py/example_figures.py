import importlib
import pandas as pd

# colors
from palettable.cartocolors.sequential import Mint_5 as COVID_colors
from palettable.colorbrewer.sequential import Reds_5 as delta_colors
from palettable.cartocolors.sequential import Burg_5 as Flu_colors
from palettable.colorbrewer.sequential import YlOrBr_5 as RSV_colors #YeOrRd_5
from palettable.cartocolors.qualitative import Safe_10 as diag_colors
my_purp = "#0A014F"
light_pink = "#F29CA3"
med_pink = "#D4215D"
dark_pink = "#710627"
black = "#042A2B"

# change directory
import sys
sys.path.insert(1,r'/Users/caseymiddleton/Documents/ActiveProjects/TestingTheory/TestingFrameworkRepo/py')

from functions import *
from parameters import *
from prettyplotlib import *

fig_path = r'/Users/caseymiddleton/Documents/ActiveProjects/TestingTheory/figs/'
data_path = r'/Users/caseymiddleton/Documents/ActiveProjects/TestingTheory/TestingFrameworkRepo/data/'

# ------------------------------------------------------------------------------------------------------------
# ---------------------------------- Figure Drivers ----------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
def main():
    # examplevirus_trajectory()
    # Fig1_generate_data()
    Draw_Fig1()
    # Fig3_generate_data()
    Draw_Fig3()
    # Fig5_generate_data()
    Draw_Fig5()
    # Fig6_generate_data()
    Draw_Fig6()
    # Fig7_generate_data()
    Draw_Fig7()

# ------------------------------------------------------------------------------------------------------------
# ---------------------------------- Figure Functions ----------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

''' Create time series data '''
start = 0
stop = 10
step = 0.01
numels = int((stop - start)/step + 1)
t_vals = np.linspace(start,stop,numels)
dt = t_vals[1]-t_vals[0]

def examplevirus_trajectory():
    '''
    Plot trajectory of examplevirus with and without testing
    '''
    inf_params = example_params_det()
    test_params = create_testing_program("example",sensitivity_threshold=2,frequency=7,delay=0)
    LOD = np.ones(len(t_vals))*test_params['sensitivity_threshold'] # Limit of Detection vector
    inf_threshold = np.ones(len(t_vals))*inf_params['Ydetect'] # Threshold for where virus becomes infectious

    # get trajectory
    vl = np.array([get_VL(t,**inf_params) for t in t_vals])
    inf = np.array([get_inf(t,**inf_params) for t in t_vals])
    det = np.array([has_results(t,**inf_params,**test_params) for t in t_vals])
    tst = np.array([prob_test(t,**test_params) for t in t_vals])

    # find infectiousness with testing
    cumulative_pr_detected = np.minimum(cumulative_rect_rule(det*tst,t_vals),1)
    removed = inf * cumulative_pr_detected
    testing = vl - removed

    # get pr daily testing
    p_test = np.array([prob_test(t,**test_params) * has_results(t,**inf_params,**test_params) for t in t_vals])
    p_tested = np.minimum(cumulative_rect_rule(p_test,t_vals),1)

    # Find first point of detectability
    first_detectable_time = 0
    for i in range(0,len(t_vals)):
        if first_detectable_time == 0:
            if inf[i] >= test_params['sensitivity_threshold']:
                first_detectable_time = t_vals[i]
        else:
            break

    #Plot results
    plt.rcParams['hatch.linewidth'] = 2  # hatch linewidth
    fig, (ax1,ax2) = plt.subplots(nrows=1,ncols=2,figsize=(12,4))

    # top left panel ------------------------------------------------------------------------------------------------
    ax1.fill_between(t_vals, vl, color="lightgray", alpha = 0.5)
    ax1.fill_between(t_vals, vl, inf_threshold, where = (vl > inf_threshold), color=LIGHT_COLOR, alpha = 1)
    ax1.plot(t_vals, vl, color=DARK_COLOR, label="Infectiousness", linewidth=3)
    # limit of detection
    #ax1.plot(t_vals, LOD, color = "white", linewidth=2)
    ax1.plot(t_vals, LOD, color = ALMOST_BLACK, linewidth=2, linestyle = "dashed")
    ax1.plot(first_detectable_time*np.ones(test_params['sensitivity_threshold']+1), range(0,test_params['sensitivity_threshold']+1),\
        color = ALMOST_BLACK, linewidth=1, linestyle = "dashed")
    ax1.set_ylim(0,10.2)
    ax1.set_xlim(0,10.5)

    ax1.set_xlabel("Time (days)")
    ax1.set_ylabel("Viral Concentration")
    hide_right_top_axis(ax1)
    no_ticks(ax1)
    finalize(ax1)

    # probability of daily testing
    ax1_right = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    #ax2.set_ylabel('Probability of Test', color=diag_colors.hex_colors[8])  # we already handled the x-label with ax1
    ax1_right.plot(t_vals, p_tested, color=dark_pink, linewidth=2)
    ax1_right.tick_params(axis='y', labelcolor=dark_pink)
    ax1_right.set_ylim(0,1.05)
    hide_right_top_axis(ax1_right)
    no_ticks(ax1_right)

    ax1.text(6.5, test_params['sensitivity_threshold'], r'Limit of Detection', ha='center', va='center', fontsize=LABEL_SIZE-1, color=ALMOST_BLACK, \
        alpha=1, backgroundcolor=LIGHT_COLOR)
    ax1.text(8.8, 8.7, 'Cumulative\n Probability\n of Detection', ha='center', va='center', fontsize=LABEL_SIZE-1, color=dark_pink, alpha=1)
    ax1.text(first_detectable_time, -0.5, 'First Detectable Instance', ha='center', va='center', fontsize=LABEL_SIZE-3, color=ALMOST_BLACK, alpha=1)

    # top right panel ------------------------------------------------------------------------------------------------
    ax2.fill_between(t_vals, vl, color="lightgray", alpha = 0.5)
    ax2.fill_between(t_vals, vl, inf_threshold, where = (vl > inf_threshold), color=LIGHT_COLOR, alpha = 1)
    ax2.fill_between(t_vals, testing, inf_threshold, where = (testing > inf_threshold), hatch='\\', \
        zorder=2, edgecolor=med_pink, facecolor=LIGHT_COLOR, alpha = 0.6)
    ax2.plot(t_vals, vl, color=DARK_COLOR, label="Infectiousness", linewidth=3)
    ax2.plot(t_vals, testing, color=med_pink, label="With testing", linewidth=3)

    ax2.text(3.5, 1.5, r'E(Infectiousness | testing)', ha='center', va='center', fontsize=LABEL_SIZE-1, color=med_pink, \
        alpha=1, backgroundcolor=LIGHT_COLOR)
    ax2.text(5, 6, r'E(Infectiousness)', ha='center', va='center', fontsize=LABEL_SIZE-1, color=DARK_COLOR, \
        alpha=1, backgroundcolor=LIGHT_COLOR)

    ax2.set_xlabel("Time (days)")
    hide_right_top_axis(ax2)
    no_ticks(ax2)
    finalize(ax2)

    plt.tight_layout()
    save_string = fig_path + "example_model_diagram_2"
    plt.savefig(save_string)
    plt.show()


def Fig1_generate_data():
    ''' Compare frequent testing with a low sensitivity test to infrequent testing with a highly sensitive test '''
    # High sensitivity test
    # get trajectory
    inf_params = example_params_det()
    vl = np.array([get_VL(t,**inf_params) for t in t_vals])
    inf = np.array([get_inf(t,**inf_params) for t in t_vals])

    header = ["low","high","1 day TAT","3 day TAT"]
    results = []
    test_frequencies = [1,3.5,7,14]
    test_sensitivities = [4,0,0,0]
    delays = [0,0,1,3]
    for freq in test_frequencies:
        scenario_results = []
        for ii in range(0,len(test_sensitivities)):
            test_params = create_testing_program("example",sensitivity_threshold=test_sensitivities[ii],frequency=freq,delay=delays[ii])
            # get R_reduction
            det = np.array([has_results(t,**inf_params,**test_params) for t in t_vals])
            tst = np.array([prob_test(t,**test_params) for t in t_vals])
            reduction = R_reduction(inf,det,tst,t_vals)
            scenario_results.append(reduction)
        results.append(scenario_results)

    # Save data
    fname = data_path + "example_fig1_data.csv"
    data = pd.DataFrame(results, columns=header)
    data.to_csv(fname, index=False)


def Draw_Fig1():
    data_file = data_path + "example_fig1_data.csv"
    data= pd.read_csv(data_file)

    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,4))

    bar_labs = ["Low Sensitivity","High Sens. - No TAT", "High Sens. - 1 Day TAT", "High Sens. - 3 Day TAT"]
    xlabels = ["Daily","Twice Weekly","Weekly","Biweekly"]
    x = np.arange(len(xlabels))  # label locations
    width = 0.20                 # width of the bars: can also be len(x) sequence

    ax.bar(x-3*width/2, data["low"], width, label=bar_labs[0], color = ALMOST_BLACK, edgecolor = ALMOST_BLACK)
    ax.bar(x-width/2, data["high"], width, label=bar_labs[1], color = dark_pink, edgecolor = ALMOST_BLACK)
    ax.bar(x+width/2, data["1 day TAT"], width, label=bar_labs[2], color = med_pink, edgecolor = ALMOST_BLACK)
    ax.bar(x+3*width/2, data["3 day TAT"], width, label=bar_labs[3], color = light_pink, edgecolor = ALMOST_BLACK)

    ax.set_xlabel('Testing Frequency')
    ax.set_ylabel('Transmission Reduction ($R/R_0$)')
    ax.set_xticks(x)
    ax.set_xticklabels(xlabels)
    ax.legend(frameon=False, fontsize=LABEL_SIZE-2, loc=(.58, .67))
    finalize(ax)
    plt.tight_layout()

    fname = fig_path + "example_Fig1.png"
    plt.savefig(fname)
    plt.show()

def Fig3_generate_data():
    ''' Expected serial interaval figure '''
    inf_params = example_params_det()
    test_params = create_testing_program("example",sensitivity_threshold=2,frequency=3.5,delay=0)

    # get trajectory
    vl = np.array([get_VL(t,**inf_params) for t in t_vals])
    inf = np.array([get_inf(t,**inf_params) for t in t_vals])
    det = np.array([has_results(t,**inf_params,**test_params) for t in t_vals])
    tst = np.array([prob_test(t,**test_params) for t in t_vals])

    # infectiousness with testing
    cumulative_pr_detected = np.minimum(cumulative_rect_rule(det*tst,t_vals),1)
    removed = inf * cumulative_pr_detected
    resid = inf - removed

    SI_notesting = calc_serial_int(inf,inf,t_vals)
    SI_withtesting = calc_serial_int(inf,resid,t_vals)

    d = {"No Testing": SI_notesting, "Testing": SI_withtesting}
    data = pd.DataFrame(data = d)
    fname = data_path + "Example_fig3_data.csv"
    data.to_csv(fname, index=False)


def Draw_Fig3():
    data_file = data_path + "example_fig3_data.csv"
    data= pd.read_csv(data_file)

    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,4))
    x_vals = t_vals[0:len(t_vals)-1] # remove last t value

    ax.plot(x_vals, data["No Testing"], label="No Testing", color = ALMOST_BLACK)
    ax.plot(x_vals, data["Testing"], label="Testing", color = dark_pink)

    ax.text(1, 0.04, r'No Testing', ha='left', va='center', fontsize=LABEL_SIZE-1, color=ALMOST_BLACK, alpha=0.8)
    ax.text(2.1, 0.01, r'Testing', ha='left', va='center', fontsize=LABEL_SIZE-1, color=dark_pink, alpha=0.8)

    ax.set_xlabel('Serial Interval (days)')
    ax.set_ylabel('Likelihood')
    #ax.legend(frameon=False, fontsize=LABEL_SIZE-2, loc=(0.2, .6))
    finalize(ax)
    plt.tight_layout()

    fname = fig_path + "example_Fig3.png"
    plt.savefig(fname)
    plt.show()

def Fig5_generate_data():
    ''' Compare R/R0 with and without symptomatic self-isolation '''
    # Symptomatic self isolation assumptions
    prop_symp = 0.5
    prop_isolate = 0.5

    # get trajectory
    inf_params = example_params_det()
    vl = np.array([get_VL(t,**inf_params) for t in t_vals])
    inf = np.array([get_inf(t,**inf_params) for t in t_vals])
    symp = np.array([example_is_symptomatic(t,**inf_params) for t in t_vals])

    header = ["low","high"]
    results = []
    test_frequencies = [1,3.5,7,14]
    test_sensitivities = [4,0]
    for freq in test_frequencies:
        scenario_results = []
        for ii in range(0,len(test_sensitivities)):
            test_params = create_testing_program("example",sensitivity_threshold=test_sensitivities[ii],frequency=freq)
            # get R_reduction
            det = np.array([has_results(t,**inf_params,**test_params) for t in t_vals])
            tst = np.array([prob_test(t,**test_params) for t in t_vals])
            reduction = R_reduction_with_symp(inf,det,tst,symp,t_vals,prop_symp,prop_isolate)
            scenario_results.append(reduction)
        results.append(scenario_results)

    # Save data
    fname = data_path + "example_fig5_data.csv"
    data = pd.DataFrame(results, columns=header)
    data.to_csv(fname, index=False)

def Draw_Fig5():
    asymp_data_file = data_path + "example_fig1_data.csv"
    symp_data_file = data_path + "example_fig5_data.csv"
    asymp_data = pd.read_csv(asymp_data_file)
    symp_data = pd.read_csv(symp_data_file)

    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,4))

    #bar_labs = ["Low Sensitivity","High Sens."]
    xlabels = ["Daily","Twice Weekly","Weekly","Biweekly"]
    symp_screening = ["With Self Isolation","W/O Self Isolation"]
    x = np.arange(len(xlabels))  # label locations
    width = 0.20                 # width of the bars: can also be len(x) sequence

    plt.rcParams['hatch.color'] = "white"  # hatch color
    ax.bar(x-3*width/2, asymp_data["low"], width, label="Low Sens. w/o Self Isolation", color = LIGHT_COLOR, edgecolor = ALMOST_BLACK)
    ax.bar(x-width/2, symp_data["low"], width, label="Low Sens. with Self Isolation", color = LIGHT_COLOR, edgecolor = ALMOST_BLACK, hatch='\\')
    ax.bar(x+width/2, asymp_data["high"], width, label="High Sens. w/o Self Isolation", color = dark_pink, edgecolor = ALMOST_BLACK)
    ax.bar(x+3*width/2, symp_data["high"], width, label="High Sens. with Self Isolation", color = dark_pink, edgecolor = ALMOST_BLACK, hatch='\\')

    ax.set_xlabel('Testing Frequency')
    ax.set_ylabel('Transmission Reduction ($R/R_0$)')
    ax.set_xticks(x)
    ax.set_xticklabels(xlabels)
    ax.legend(frameon=False, fontsize=LABEL_SIZE-2, loc=(.48, .67))
    finalize(ax)
    plt.tight_layout()

    fname = fig_path + "example_Fig5.png"
    plt.savefig(fname)
    plt.show()


def Fig6_generate_data():
    ''' Multiple screening strategies '''
    # Symptomatic self isolation assumptions
    prop_symp = 0.5
    prop_isolate = 0.5

    # get trajectory
    inf_params = example_params_det()
    vl = np.array([get_VL(t,**inf_params) for t in t_vals])
    inf = np.array([get_inf(t,**inf_params) for t in t_vals])

    # test 1 parameters : viral detection test
    test1_LOD = 4
    test1_params = create_testing_program("example",sensitivity_threshold=test1_LOD)
    test1_det = np.array([has_results(t,**inf_params,**test1_params) for t in t_vals])

    # test 2 parameters : symptom screening
    test2_det = np.array([example_is_symptomatic(t,**inf_params) for t in t_vals])
    symp = test2_det

    results = []
    test_frequencies = np.arange(start=1,stop=15,step=1)
    header = test_frequencies
    for freq_test1 in test_frequencies:
        test1_params = create_testing_program("example",frequency=freq_test1)
        tst1 = np.array([prob_test(t,**test1_params) for t in t_vals])
        scenario_results = []
        for freq_test2 in test_frequencies:
            test2_params = create_testing_program("example",frequency=freq_test2)
            tst2 = np.array([prob_test(t,**test2_params) for t in t_vals])
            # get R_reduction
            det = [test1_det,test2_det]
            tst = [tst1,tst2]
            reduction = R_reduction_multiscreen(inf,det,tst,symp,t_vals,prop_symp,prop_isolate)
            scenario_results.append(reduction)
        results.append(scenario_results)

    # Save data
    fname = data_path + "example_fig6_data.csv"
    # rows = temperature screening, columns = viral screening
    data = pd.DataFrame(results, columns=header)
    data.to_csv(fname, index=False)

def Draw_Fig6():
    data_file = data_path + "example_fig6_data.csv"
    data = pd.read_csv(data_file)

    test_frequencies = np.arange(start=1,stop=17,step=2) - 2

    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,4))

    plt.imshow(data)

    ax.set_xlabel('Temperature Screening Frequency (days)')
    ax.set_ylabel('Viral Screening Frequency (days)')
    ax.set_xticklabels(test_frequencies)
    ax.set_yticklabels(test_frequencies)

    cbar = plt.colorbar()
    cbar.ax.get_yaxis().labelpad = 10
    cbar.ax.set_ylabel('Testing Effectiveness')
    plt.tight_layout()

    fname = fig_path + "example_Fig6.png"
    plt.savefig(fname)
    plt.show()

def Fig7_generate_data():
    ''' TE for various freq, sens, and TAT params '''
    # get trajectory
    inf_params = example_params_det()
    vl = np.array([get_VL(t,**inf_params) for t in t_vals])
    inf = np.array([get_inf(t,**inf_params) for t in t_vals])

    ''' Testing frequency '''
    test_LOD = 2

    df = pd.DataFrame()
    test_frequencies = np.arange(start=1,stop=31,step=1)
    participations = [1, 1, 0.5, 0.5]
    compliances = [1, 0.5, 1, 0.5]
    for ii in range(0,len(participations)):
        scenario_results = []
        for freq_test in test_frequencies:
            test_params = create_testing_program("example",frequency=freq_test,compliance=compliances[ii],sensitivity_threshold=test_LOD)
            tst = np.array([prob_test(t,**test_params) for t in t_vals])
            det = np.array([has_results(t,**inf_params,**test_params) for t in t_vals])
            reduction = R_reduction(inf,det,tst,t_vals,participation=participations[ii])
            scenario_results.append(reduction)
        df_string = "part" + str(participations[ii]) + "_comp" + str(compliances[ii])
        df[df_string] = scenario_results

    # Save data
    fname = data_path + "example_fig7_freq.csv"
    df.to_csv(fname, index=False)

    ''' Test sensitivity '''
    freq = 1 # daily testing with no turnaround time

    df = pd.DataFrame()
    test_sensitivities = np.arange(start=0,stop=np.max(vl),step=1)
    participations = [1, 1, 0.5, 0.5]
    compliances = [1, 0.5, 1, 0.5]
    for ii in range(0,len(participations)):
        scenario_results = []
        for sens in test_sensitivities:
            test_params = create_testing_program("example",frequency=freq,compliance=compliances[ii],sensitivity_threshold=sens)
            tst = np.array([prob_test(t,**test_params) for t in t_vals])
            det = np.array([has_results(t,**inf_params,**test_params) for t in t_vals])
            reduction = R_reduction(inf,det,tst,t_vals,participation=participations[ii])
            scenario_results.append(reduction)
        df_string = "part" + str(participations[ii]) + "_comp" + str(compliances[ii])
        df[df_string] = scenario_results
    # Save data
    fname = data_path + "example_fig7_sens.csv"
    df.to_csv(fname, index=False)

    ''' Test turnaround time '''
    freq = 1
    sens = 2

    df = pd.DataFrame()
    test_TATs = np.arange(start=0,stop=7.5,step=0.5)
    participations = [1, 1, 0.5, 0.5]
    compliances = [1, 0.5, 1, 0.5]
    for ii in range(0,len(participations)):
        scenario_results = []
        for TAT in test_TATs:
            test_params = create_testing_program("example",frequency=freq,compliance=compliances[ii],sensitivity_threshold=sens,delay=TAT)
            tst = np.array([prob_test(t,**test_params) for t in t_vals])
            det = np.array([has_results(t,**inf_params,**test_params) for t in t_vals])
            reduction = R_reduction(inf,det,tst,t_vals,participation=participations[ii])
            scenario_results.append(reduction)
        df_string = "part" + str(participations[ii]) + "_comp" + str(compliances[ii])
        df[df_string] = scenario_results
    # Save data
    fname = data_path + "example_fig7_TAT.csv"
    df.to_csv(fname, index=False)

def Draw_Fig7():
    # read in data
    data_file = data_path + "example_fig7_freq.csv"
    freq_data = pd.read_csv(data_file)
    data_file = data_path + "example_fig7_sens.csv"
    sens_data = pd.read_csv(data_file)
    data_file = data_path + "example_fig7_TAT.csv"
    TAT_data = pd.read_csv(data_file)

    participations = [1, 1, 0.5, 0.5]
    compliances = [1, 0.5, 1, 0.5]
    frequencies = np.arange(start=1,stop=31,step=1)
    test_sensitivities = np.arange(start=0,stop=10,step=1)
    test_TATs = np.arange(start=0,stop=7.5,step=0.5)

    colors = [ALMOST_BLACK, dark_pink, med_pink, light_pink]

    fig,(ax1,ax2,ax3) = plt.subplots(nrows=1,ncols=3,figsize=(12,4))

    for ii in range(0,len(participations)):
        # get data
        df_string = "part" + str(participations[ii]) + "_comp" + str(compliances[ii])
        freq_plot = freq_data[df_string]
        sens_plot = sens_data[df_string]
        TAT_plot = TAT_data[df_string]
        # specify legend
        if participations[ii] == 1:
            if compliances[ii] == 1:
                lab = "High participation, high compliance"
            else:
                lab = "High participation, low compliance"
        else:
            if compliances[ii] == 1:
                lab = "Low participation, high compliance"
            else:
                lab = "Low participation, low compliance"
        # plot
        ax1.plot(frequencies,freq_plot,color=colors[ii],label=lab)
        ax2.plot(test_sensitivities,sens_plot,color=colors[ii])
        ax3.plot(test_TATs,TAT_plot,color=colors[ii])

    ax1.set_xlabel('Test Frequency (days)')
    ax1.set_ylabel('Testing Effectiveness (TE)')
    ax1.legend()
    finalize(ax1)

    ax2.get_yaxis().set_ticks([])
    ax2.set_xlabel('Test Sensitivity $(\log_{10}$ cp RNA/mL)')
    finalize(ax2)
    ax3.get_yaxis().set_ticks([])
    ax3.set_xlabel('Test Turnaround Time (days)')
    finalize(ax3)

    # ax.text(2, 1, 'High participation, high compliance', ha='left', va='center', fontsize=LABEL_SIZE-1, color=ALMOST_BLACK)
    # ax.plot([3,10],[freq_data['part1_comp0.5'][2], freq_data['part1_comp0.5'][2]],color=dark_pink)
    # ax.text(10, freq_data['part1_comp0.5'][2], 'High part., low comp.', ha='left', va='center', fontsize=LABEL_SIZE-1, color=dark_pink)
    # ax.text(0.1, 0.6, 'Low part.,\n high comp.', ha='left', va='center', fontsize=LABEL_SIZE-1, color=med_pink)
    # ax.text(5, 0.05, 'Low part., low comp.', ha='left', va='center', fontsize=LABEL_SIZE-1, color=light_pink)

    plt.tight_layout()

    fname = fig_path + "example_Fig7.png"
    plt.savefig(fname)
    plt.show()



if __name__ == "__main__":
    main()
