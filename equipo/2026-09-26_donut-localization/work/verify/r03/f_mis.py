# -*- coding: utf-8 -*-
"""Verifier r03: own misalignment MC. L=100, N=500, SBR=10, every zero (4) displaced by delta in an
independent uniform direction; counts from the true model; honest MLE (true centres) vs naive (ideal
TCP), both in disk 0.75 L around the origin. P patterns x R reps. Own seed."""
import json, os, sys, time
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(H, '..'))
from mymle import tcp, pvec, mle
ROOT = os.path.abspath(os.path.join(H, *['..']*5))
PN = json.load(open(os.path.join(ROOT, 'data', 'paper_numbers.json')))
L, N, SBR = 100., 500, 10.
P = int(sys.argv[1]) if len(sys.argv) > 1 else 400; R = int(sys.argv[2]) if len(sys.argv) > 2 else 200
C0 = tcp(L); D = (0., 2., 5., 10.); POS = {'center': np.array([0., 0.]), 'Lq': np.array([L/4, 0.])}
dirs = 2*np.pi*np.random.default_rng(31337).random((P, 4))
rng = np.random.default_rng(4242)
res = {}
t0 = time.time()
for pn, r0 in POS.items():
    for d in D:
        st = {m: dict(b=[], s=[], ms=[]) for m in ('naive', 'honest')}
        same = True
        for j in range(P):
            C = C0 + d*np.stack([np.cos(dirs[j]), np.sin(dirs[j])], 1)
            n = rng.multinomial(N, pvec(r0, C, SBR), size=R)
            for m, Cm in (('naive', C0), ('honest', C)):
                e = mle(n, Cm, SBR, 0.75*L, ncoarse=60) - r0
                st[m]['b'].append(np.linalg.norm(e.mean(0))); st[m]['s'].append(np.sqrt(e.var(0, ddof=1).sum()/2))
                st[m]['ms'].append((e**2).sum(1).mean()/2)
                if m == 'naive': en = e
                else: same &= bool(np.array_equal(en, e))
        for m in st:
            b, s = np.array(st[m]['b']), np.array(st[m]['s'])
            res[(m, pn, d)] = dict(bias=b.mean(), bias_se=b.std(ddof=1)/np.sqrt(P), sigma=s.mean(), sigma_se=s.std(ddof=1)/np.sqrt(P),
                                   rmse=np.sqrt(np.mean(st[m]['ms'])), floor=s.mean()*np.sqrt(np.pi/(2*R)))
            kb = 'misalignment_%s_bias_abs_%s_d%d_nm' % (m, pn, d); ks = 'misalignment_%s_sigma_%s_d%d_nm' % (m, pn, d)
            kr = 'misalignment_%s_rmse_%s_d%d_nm' % (m, pn, d)
            q = res[(m, pn, d)]
            tb, ts = PN[kb], PN[ks]
            print('%-6s %-6s d=%4.1f |bias|=%.3f±%.3f (theirs %.3f±%.3f, z=%+.2f) sigma=%.3f±%.4f (theirs %.3f±%.4f, z=%+.2f) rmse=%.3f (theirs %.3f) floor=%.3f' %
                  (m, pn, d, q['bias'], q['bias_se'], tb['value'], tb['se'], (q['bias']-tb['value'])/np.hypot(q['bias_se'], tb['se']),
                   q['sigma'], q['sigma_se'], ts['value'], ts['se'], (q['sigma']-ts['value'])/np.hypot(q['sigma_se'], ts['se']),
                   q['rmse'], PN[kr]['value'], q['floor']), flush=True)
        print('   delta=%g identical naive/honest: %s   [%.0fs]' % (d, same, time.time()-t0), flush=True)
for pn in POS:
    dd = np.array(D[1:]); b = np.array([res[('naive', pn, d)]['bias'] for d in dd]); s = np.array([res[('naive', pn, d)]['bias_se'] for d in dd])
    w = 1/s**2; sl = (w*dd*b).sum()/(w*dd*dd).sum(); se = 1/np.sqrt((w*dd*dd).sum())
    k = 'misalignment_naive_bias_over_delta_%s' % pn; t = PN[k]
    print('slope naive |bias|/delta %-6s = %.4f ± %.4f  theirs %.4f ± %.4f  z=%+.2f   (unweighted per-delta ratios %s)' %
          (pn, sl, se, t['value'], t['se'], (sl-t['value'])/np.hypot(se, t['se']), np.round(b/dd, 3)))
json.dump({'|'.join(map(str, k)): v for k, v in res.items()}, open(os.path.join(H, 'f_mis_%dx%d.json' % (P, R)), 'w'), indent=1, default=float)
