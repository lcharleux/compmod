import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from compart_settings import settings

df = pd.read_csv(settings["workdir"] + settings["label"] + ".csv")

fig = plt.figure()

sp1 = fig.add_subplot(2, 1, 1)
plt.plot(df.disp, df.force, 'ro-')
plt.xlabel('Displacement, $U$')
plt.ylabel('Force, $F$')
plt.grid()
sp1 = fig.add_subplot(2, 1, 2)
plt.plot(df.strain, df.stress, 'ro-', label = 'simulation curve', linewidth = 2.)
plt.xlabel('Tensile Strain, $\epsilon$')
plt.ylabel(' Tensile Stress $\sigma$')
plt.grid()
plt.savefig(settings['workdir'] + settings['label'] + '.pdf')
