# -*- coding: utf-8 -*-
"""Verifier r03: own iterative MINFLUX (protocol from the docstring only): emitters uniform in disk
L0/4, TCP_0 at origin, fixed geometric L 150->25 (4 it.), equal photon split (remainder to last),
MLE in disk 0.75 L_k around c_k, re-centring, final = last estimate. Own seeds."""
import json, os, sys, time
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(H, '..'))
from mymle import tcp, pvec, mle, sig_c
ROOT = os.path.abspath(os.path.join(H, *['..']*5))
PN = json.load(open(os.path.join(ROOT, 'data', 'paper_numbers.json')))

def split(N, K=4):
    Nk = np.full(K, N//K); Nk[-1] += N-Nk.sum(); return Nk

def run(Ntot, reps, seed, Ls=None, sbr=np.inf, recenter=True, rdisk=37.5):
    rng = np.random.default_rng(seed)
    if Ls is None: Ls = 150.*(25./150.)**(np.arange(4)/3)
    rr = rdisk*np.sqrt(rng.random(reps)); th = 2*np.pi*rng.random(reps)
    x = np.stack([rr*np.cos(th), rr*np.sin(th)], 1); c = np.zeros_like(x)
    Nk = split(Ntot, len(Ls)); errs = []
    for L, nk in zip(Ls, Nk):
        C = tcp(L)
        n = rng.multinomial(nk, pvec(x-c, C, sbr))
        est = c + mle(n, C, sbr, 0.75*L)
        errs.append(est-x)
        if recenter: c = est
    return errs

def bse(e, nb, rng):
    R = len(e); return np.array([sig_c(e[rng.integers(0, R, R)]) for _ in range(nb)])

if __name__ == '__main__':
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
    NL = [250, 500, 1000, 2000, 4000, 8000]
    brng = np.random.default_rng(99)
    sig, se, B, out = [], [], [], {}
    for i, N in enumerate(NL):
        t = time.time()
        e = run(N, reps, 5150+i)[-1]
        s = sig_c(e); b = bse(e, 1000, brng)
        sig.append(s); se.append(b.std(ddof=1)); B.append(b)
        k = 'iterative_sigma_N%d_nm' % N; tv, ts = PN[k]['value'], PN[k]['se']
        print('N=%5d sigma=%.4f±%.4f (gauss %.4f)  theirs %.4f±%.4f  z=%+.2f  ratio_cam=%.4f  kurt=%.1f [%.0fs]' %
              (N, s, se[-1], s/(2*np.sqrt(reps)), tv, ts, (s-tv)/np.hypot(se[-1], ts), s/(100/np.sqrt(N)),
               np.mean(((e-e.mean(0))/e.std(0))**4), time.time()-t), flush=True)
        out[N] = dict(sigma=s, se=se[-1])
    NL = np.array(NL, float); sig = np.array(sig); B = np.stack(B, 1)
    for name, m in (('all', NL > 0), ('N_ge_500', NL >= 500)):
        sl = np.polyfit(np.log(NL[m]), np.log(sig[m]), 1)[0]
        bs = np.array([np.polyfit(np.log(NL[m]), np.log(B[j, m]), 1)[0] for j in range(B.shape[0])])
        lo, hi = np.percentile(bs, [2.5, 97.5]); tv = PN['iterative_slope_%s' % name]
        print('slope %-9s = %.4f ± %.4f  CI95 [%.4f, %.4f]   theirs %.4f ± %.4f [%.4f, %.4f]  z=%+.2f' %
              (name, sl, bs.std(ddof=1), lo, hi, tv['value'], tv['se'], PN['iterative_slope_%s_ci95_lo' % name]['value'],
               PN['iterative_slope_%s_ci95_hi' % name]['value'], (sl-tv['value'])/np.hypot(bs.std(ddof=1), tv['se'])))
        out['slope_'+name] = dict(v=sl, se=bs.std(ddof=1), ci=[lo, hi])
    r = sig/(100/np.sqrt(NL)); print('ratio range %.4f..%.4f (theirs %.4f..%.4f)' % (r.min(), r.max(), PN['iterative_ratio_min']['value'], PN['iterative_ratio_max']['value']))
    # extras at N=1000
    for tag, kw, key in (('SBR=10', dict(sbr=10.), 'iterative_sbr10_sigma_nm'),
                         ('adaptive [150,25,25,25]', dict(Ls=np.array([150., 25, 25, 25])), 'iterative_adaptive_sigma_nm'),
                         ('recenter=False', dict(recenter=False), 'iterative_no_recenter_sigma_nm_ARTEFACT')):
        e = run(1000, reps, 6100+len(tag), **kw)
        s = sig_c(e[-1]); b = bse(e[-1], 500, brng).std(ddof=1); tv = PN[key]
        print('%-24s sigma=%.4f±%.4f theirs %.4f±%.4f z=%+.2f   per-iter %s' % (tag, s, b, tv['value'], tv['se'], (s-tv['value'])/np.hypot(b, tv['se']),
              np.round([sig_c(q) for q in e], 3)))
    # adaptive iteration-0 stats (L=150, 250 photons)
    e0 = run(250, reps, 6200, Ls=np.array([150.]))[0]
    d = np.sqrt((e0**2).sum(1)); f = np.mean(d > 12.5)
    print('iter0: sigma=%.4f (theirs %.4f±%.4f)  frac |err|>12.5 = %.4f±%.4f (theirs %.4f±%.4f)  frac>3*3.47=%.4f' %
          (sig_c(e0), PN['adaptive_iter0_sigma_nm']['value'], PN['adaptive_iter0_sigma_nm']['se'], f, np.sqrt(f*(1-f)/reps),
           PN['adaptive_frac_outside_next_radius']['value'], PN['adaptive_frac_outside_next_radius']['se'], np.mean(d > 3*3.46995)))
    json.dump({str(k): v for k, v in out.items()}, open(os.path.join(H, 'e_iter.json'), 'w'), indent=1, default=float)
