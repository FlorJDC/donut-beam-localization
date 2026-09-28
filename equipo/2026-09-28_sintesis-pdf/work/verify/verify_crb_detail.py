# -*- coding: utf-8 -*-
import numpy as np
exec(open('verify_pminflux_B.py').read().split('# --- MLE propio')[0].split('print("== 2. CRB ratios")')[0])
print("ideal disk mean, excl centre:", c0[(xd**2+yd**2)>0].mean())
r = np.hypot(xd, yd)
for tau in (1,2,3,4,5):
    M = Mmat(tau); cm = crb(xd, yd, M)
    c0L = c0.copy(); c0L[r==0] = 33.63411105379955
    print(tau, "%.5f  (limit-centre ideal) %.5f" % (cm.mean()/c0.mean(), cm.mean()/c0L.mean()))
