import data_generators as dg
from parameters import *

#pathogen = "COVID_wt"
#dg.line_Q_LOD(pathogen)

pathogen = "COVID_wt"
dg.generate_stochastic_trajectories(pathogen,n=1)
