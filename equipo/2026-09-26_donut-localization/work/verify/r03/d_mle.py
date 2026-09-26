# -*- coding: utf-8 -*-
import json, os, sys, time
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(H, '..'))
from mymle import tcp, pvec, mle, sig_c, boot_se
from vfisher import crb
ROOT = os.path.abspath(os.path.join(H, *['..']*5))
PN = json.load(open(os.path.join(ROOT, 'data', 'paper_numbers.json')))
L = 50.; C = tcp(L); SEED = 777031
res = {}
def cmp(k, v, se):
    t = PN[k]; tse = t.get('se', 0) or 0; z = (v-t['value'])/np.hypot(se, tse)
    res[k] = dict(mine=v, se=se, theirs=t['value'], their_se=tse, z=z)
    print('%-42s mine=%.4f±%.4f theirs=%.4f±%.4f z=%+.2f' % (k, v, se, t['value'], tse, z), flush=True)
def run(r0, N, sbr, R, reps, seed, ncoarse=60):
    rng = np.random.default_rng(seed)
    n = rng.multinomial(N, pvec(np.array(r0, float), C, sbr), size=reps)
    return mle(n, C, sbr, R, ncoarse=ncoarse)-np.array(r0, float), n
t = time.time()
e, _ = run([0, 0], 100, 10., 50., 10000, SEED)
s = sig_c(e); se = boot_se(e, 1000, 1); c10 = crb([0, 0], L, 300., 100, sbr=10)
cmp('mle_sigma_center_sbr10_nm', s, se); cmp('mle_efficiency_center', s/c10, se/c10)
print('  bias', e.mean(0), ' gaussian SE', s/np.sqrt(2*len(e))/c10, ' max|e|', np.abs(e).max(), time.time()-t)
e, n = run([0, 0], 100, np.inf, 50., 10000, SEED+1)
s = sig_c(e); se = boot_se(e, 1000, 2); lim = crb([1e-7, 0], L, 300., 100); s27 = crb([0, 0], L, 300., 100)
cmp('mle_nobg_sigma_center_nm', s, se); cmp('mle_nobg_sigma_over_crb_limit_center', s/lim, se/lim); cmp('mle_nobg_sigma_over_s27_center', s/s27, se/s27)
e, _ = run([2, 0], 100, np.inf, 50., 40000, SEED+2)
cmp('mle_nobg_bias_x_r2_nm', e[:, 0].mean(), e[:, 0].std()/np.sqrt(len(e)))
for N in (100, 1000):
    t = time.time()
    e, _ = run([15, 0], N, 10., 100., 10000, SEED+3+N, ncoarse=160)
    s = sig_c(e); se = boot_se(e, 1000, 3); cr = crb([15, 0], L, 300., N, sbr=10)
    cmp('mle_sbr10_sigma_over_crb_x15_N%d' % N, s/cr, se/cr)
    d = np.sqrt((e**2).sum(1))
    print('  N=%d crb=%.4f  frac |err|>5crb: %.4f  robust sigma (1.4826 MAD per axis)/crb: %.3f  max|err| %.1f  bias %s  [%.0fs]' %
          (N, cr, np.mean(d > 5*cr), np.mean(1.4826*np.median(np.abs(e-np.median(e, 0)), 0))/cr, d.max(), np.round(e.mean(0), 3), time.time()-t))
    if N == 100:
        # sensitivity to disk radius (outlier tail is protocol dependent?)
        for Rd in (50., 75.):
            e2, _ = run([15, 0], N, 10., Rd, 10000, SEED+3+N, ncoarse=int(1.6*Rd))
            print('   disk R=%g: sigma/crb=%.3f±%.3f' % (Rd, sig_c(e2)/cr, boot_se(e2, 300, 4)/cr))
# LMS no background, centre (Eq. S49-S50 linear estimator)
rng = np.random.default_rng(SEED+9); N = 100
n = rng.multinomial(N, pvec(np.zeros(2), C, np.inf), size=10000)/N
g = 1-L**2*np.log(2)/300.**2
e = -(n[:, :3]@C[:3])/g
s = sig_c(e); se = boot_se(e, 1000, 5); s27 = crb([0, 0], L, 300., N); lim = crb([1e-7, 0], L, 300., N)
cmp('lms_nobg_sigma_over_s27_center', s/s27, se/s27); cmp('lms_nobg_sigma_over_crb_limit_center', s/lim, se/lim)
json.dump(res, open(os.path.join(H, 'd_mle.json'), 'w'), indent=1)
