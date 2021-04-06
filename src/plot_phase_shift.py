import numpy as np
import matplotlib.pyplot as plt

e = np.linspace(-10,10,50)
es = 0.1
gamma = 0.01
phase = np.pi/2.0 + np.atan((e-es)/gamma)
plt.plot(e,phase)
plt.show()