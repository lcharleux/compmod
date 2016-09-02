import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from compart_settings import settings

df = pd.read_csv(settings["workdir"] + settings["label"] + ".csv")
exp = pd.read_csv(settings["expdir"] + settings["experiment"], 
                  delim_whitespace = True,
                  comment = "#")
exp.stress *= settings['exp_stress_factor']
exp.strain *= settings['exp_strain_factor']

fig = plt.figure()

plt.plot(df.strain, df.stress, 'r-', label = 'Simulation', linewidth = 2.)
plt.plot(exp.strain, exp.stress, 'b-', label = 'Experiment', linewidth = 2.)
plt.xlabel('Tensile Strain, $\epsilon$')
plt.ylabel(' Tensile Stress $\sigma$')
plt.grid()
plt.legend(loc = "best")
plt.savefig(settings['workdir'] + settings['label'] + '_exp.pdf')
