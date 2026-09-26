# -*- coding: utf-8 -*-
import numpy as np, sys
sys.path.insert(0, '..')
from vfisher import crb
from scipy.optimize import minimize_scalar
for ped in ('gauss', 'const'):
    print('pedestal', ped)
    for eps in (0.002, 0.01, 0.05, 0.15):
        f = lambda L: crb([0., 0.], L, 300., 100, eps=eps, ped=ped)
        LL = np.linspace(1, 400, 800); v = [f(L) for L in LL]; i = int(np.argmin(v))
        r = minimize_scalar(f, bounds=(LL[max(i-2,0)], LL[i+2]), method='bounded', options={'xatol': 1e-8})
        print('  eps=%g L_opt=%.3f CRB=%.4f  L/(fwhm sqrt eps)=%.3f' % (eps, r.x, r.fun, r.x/(300*np.sqrt(eps))))
# eps=0 monotone?
LL = np.linspace(0.5, 600, 2000)
v0 = np.array([crb([0., 0.], L, 300., 100) for L in LL])
vl = np.array([crb([1e-6*L, 0.], L, 300., 100) for L in LL])
print('eps=0 point: monotone increasing?', np.all(np.diff(v0) > 0), v0[:3], ' limit:', np.all(np.diff(vl) > 0))
# other fwhm scaling check
for fw in (200., 400.):
    for eps in (0.002, 0.01, 0.05):
        f = lambda L: crb([0., 0.], L, fw, 100, eps=eps, ped='gauss')
        r = minimize_scalar(f, bounds=(1, 2*fw), method='bounded', options={'xatol': 1e-8})
        print('fwhm', fw, 'eps', eps, 'Lopt/(fwhm sqrt eps)=%.3f' % (r.x/(fw*np.sqrt(eps))))
