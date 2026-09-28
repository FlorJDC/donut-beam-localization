# -*- coding: utf-8 -*-
import numpy as np
_src = open('verify_pminflux_B.py').read(); exec(_src.split('print("== 3. sesgo')[0])
exec(_src.split('# --- 4. flickering')[1].split('NR = 2000')[0])
# rejilla desplazada 1.25 nm (no cae en ceros de las donas)
g = np.arange(-50+1.25, 50, 2.5); X, Y = np.meshgrid(g, g); m = X**2+Y**2 <= 2500
xs_, ys_ = X[m], Y[m]; cs = crb(xs_, ys_)
print("shifted grid n", xs_.size, [round(crb(xs_, ys_, Mmat(t)).mean()/cs.mean(), 5) for t in (1,2,3,4,5)])
# fina 0.5 nm desplazada
g = np.arange(-50+0.25, 50, 0.5); X, Y = np.meshgrid(g, g); m = X**2+Y**2 <= 2500
xs_, ys_ = X[m], Y[m]; cs = crb(xs_, ys_)
print("fine grid n", xs_.size, [round(crb(xs_, ys_, Mmat(t)).mean()/cs.mean(), 5) for t in (1,2,3,4,5)])
# sesgo en malla de 11 puntos
xs = np.linspace(0, 50, 11)
for tau, order in [(3,(0,1,2,3)),(3,(0,2,1,3)),(5,(0,1,2,3)),(5,(0,2,1,3))]:
    M = Mmat(tau, order); e = mle(p_of(xs, 0*xs)@M.T, p_of) - np.c_[xs, 0*xs]
    print(tau, order, np.round(np.hypot(e[:,0], e[:,1]), 3))
# signo de la fuga: con orden A el foton del slot s cae en s+1 (no s-1)
M = Mmat(5.0); print("M[1,0] (0->1) %.4f  M[3,0] (0->3) %.4f" % (M[1,0], M[3,0]))

ivs = telegraph(6000); p0 = p_of(np.array(0.), np.array(0.))
W = np.array([w_int(iv) for iv in ivs]); ok = (W*p0).sum(1) > 0
f = (W[ok]*p0)/(W[ok]*p0).sum(1, keepdims=True); e = mle(f, p_of)
print('interleaved n=6000 pop sigma_fl %.4f nm; rel weight std %.5f' % (np.sqrt(0.5*np.mean((e**2).sum(1))), (W[ok]/W[ok].mean(1,keepdims=True)).std()))
# N-independencia (SimuFLUX): seq5 con N=1000 multinomial
ivs = telegraph(1500); W = np.array([w_seq(iv, 5) for iv in ivs]); ok = (W*p0).sum(1) > 0
qq = (W[ok]*p0)/(W[ok]*p0).sum(1, keepdims=True)
ep = mle(qq, p_of); sfl = np.sqrt(0.5*np.mean((ep**2).sum(1)))
for N in (100, 1000):
    c = np.array([rng.multinomial(N, q) for q in qq]).astype(float); c0_ = rng.multinomial(N, p0, size=len(qq)).astype(float)
    e = mle(c, p_of); e0 = mle(c0_, p_of); v = 0.5*e.var(0).sum(); v0 = 0.5*e0.var(0).sum()
    print('seq5 N=%d STD %.3f ctrl %.3f sfl_from_N %.3f  pop %.3f' % (N, np.sqrt(v), np.sqrt(v0), np.sqrt(max(v-v0,0)), sfl))
