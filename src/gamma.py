
import numpy as np
from scipy.stats import gamma
import scipy.special as sc
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2)

A=np.array(np.linspace(19.5,39,10))
chisq=[1,10,100,1000]
for a in range(1,len(A)):
	P=sc.gammainc(A[a],chisq);rv=gamma(A[a])
	ax[0].plot(chisq,P,label=f'N={np.round(A[a])}')
	#mean, var, skew, kurt = gamma.stats(a, moments='mvsk')
	#x = np.linspace(gamma.ppf(0.01, A[a]),gamma.ppf(0.99,A[a]), 100)
	#ax[1].plot(x, gamma.pdf(x, A[a]),lw=5, alpha=0.6, label=f'N={np.round(A[a])}')
	if min(P)>0:
        	ax[1].plot(chisq,np.log(P))
                
	#ax[1].plot(x, rv.pdf(x), lw=2)
# 	vals = gamma.ppf([0.001, 0.5, 0.999], a)
ax[0].set_xlabel("chisq")#;ax[1].set_xlabel('gamma')
ax[0].set_ylabel(f"P(z,a)");ax[1].set_ylabel('gamma pdf')
plt.legend()
plt.show()
#ax[1].title('Gamma ')
plt.tight_layout()
# print(A
