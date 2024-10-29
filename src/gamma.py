
import numpy as np
from scipy.stats import gamma
import scipy.special as sc
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 1)
md"""
# a=39
"""
A=np.array(np.linspace(19.5,39,39))
chisq=[0, 1, 10, 100,1000]
for a in A:
	# P=sc.gammainc(a,chisq)
	# np.sum(P)
	# ax.plot(chisq,a,label=f'P(z,a)')
	mean, var, skew, kurt = gamma.stats(a, moments='mvsk')
	x = np.linspace(gamma.ppf(0.01, a),
	                gamma.ppf(0.99, a), 100)
	ax.plot(x, gamma.pdf(x, a),
	       'r-', lw=5, alpha=0.6, label=f'gamma pdf {np.round(a)}')
# 	rv = gamma(a)
	# ax.plot(x, rv.pdf(x), 'k-', lw=2)
# 	vals = gamma.ppf([0.001, 0.5, 0.999], a)

plt.legend()
plt.show()
# print(A