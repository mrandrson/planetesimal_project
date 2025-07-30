import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

dat = pd.read_csv('positionscentralmass.csv')

x = dat['x']
y = dat['y']
r = np.sqrt(x**2+y**2)
plt.plot(x, y)
plt.show()

