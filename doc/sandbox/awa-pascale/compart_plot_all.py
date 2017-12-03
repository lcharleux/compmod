# Plots all the simulations in the workdir
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from compart_settings import settings

exp = pd.read_csv(settings["expdir"] + settings["experiment"], 
                  delim_whitespace = True)
exp.stress *= settings['exp_stress_factor']
exp.strain *= settings['exp_strain_factor']
simpath = [p for p in os.listdir(settings["workdir"]) if p.endswith(".xls") and p.startswith(settings["label"])]


fig = plt.figure()
for sp in simpath:
  sim = pd.read_excel(settings["workdir"] + sp)
  plt.plot(sim.strain, sim.stress, 'r-', label = sp, linewidth = 2.)
plt.plot(exp.strain, exp.stress, 'b-', label = 'Experiment', linewidth = 2.)
plt.xlabel('Tensile Strain, $\epsilon$')
plt.ylabel(' Tensile Stress $\sigma$')
plt.grid()
plt.legend(loc = "best")
plt.savefig(settings['workdir'] + settings['label'] + '_all_stress_strain.pdf')

