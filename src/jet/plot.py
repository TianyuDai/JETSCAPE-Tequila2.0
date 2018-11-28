import matplotlib.pyplot as plt
import numpy as np
data = np.loadtxt("LLBT/rate1d_table_muperp1.393422gg.dat")
x = data.T[0]
y = data.T[1]
yp = [np.exp(y[i])/(x[i]*x[i]+1) for i in range(len(x))]
print(x[yp.index(max(yp))])
plt.plot(x, yp)
plt.axvline(x=0)
plt.show()
