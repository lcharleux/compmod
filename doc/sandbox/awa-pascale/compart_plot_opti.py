import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from compart_settings import settings

exp = pd.read_csv(settings["expdir"] + settings["experiment"], 
                  delim_whitespace = True)
exp.stress *= settings['exp_stress_factor']
exp.strain *= settings['exp_strain_factor']
opti = pd.read_csv(settings["workdir"] + settings["label"] + "_opti_results.csv", index_col = "sim_id")
sim = pd.read_csv(settings["workdir"] + settings["label"] + "_{0}.csv".format(int(opti.index.max())))


fig = plt.figure()
plt.plot(sim.strain, sim.stress, 'r-', label = 'Simulation', linewidth = 2.)
plt.plot(exp.strain, exp.stress, 'b-', label = 'Experiment', linewidth = 2.)
plt.xlabel('Tensile Strain, $\epsilon$')
plt.ylabel(' Tensile Stress $\sigma$')
plt.grid()
plt.legend(loc = "best")
plt.savefig(settings['workdir'] + settings['label'] + '_opti_stress_strain.pdf')

fig = plt.figure()
ax = fig.add_subplot(3,1,1)
plt.plot(opti.index, opti.sy_mean, 'or-',linewidth = 2.)
plt.xlabel("Simulation id")
plt.ylabel('Mean Yield Stress, $\sigma_{ym}$')
plt.grid()
ax = fig.add_subplot(3,1,2)
plt.plot(opti.index, opti.n_bil, 'or-',linewidth = 2.)
plt.xlabel("Simulation id")
plt.ylabel('Hardening Slope, $n$')
plt.grid()
ax = fig.add_subplot(3,1,3)
plt.plot(opti.index, opti.sigma_sat, 'or-',linewidth = 2.)
plt.xlabel("Simulation id")
plt.ylabel('Staturation Stress, $\sigma_{sat}$')
plt.grid()
plt.tight_layout()
plt.savefig(settings['workdir'] + settings['label'] + '_opti_convergence.pdf')
