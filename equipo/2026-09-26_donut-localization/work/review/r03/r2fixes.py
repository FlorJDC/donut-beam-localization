import sys, time, warnings, numpy as np
sys.path.insert(0, "src")
from donutloc import fisher, photons, patterns, beams, experiments as ex, vectorial as v
beam = beams.make_beam("donut", fwhm=300.0)
p = photons.make_model(patterns.tcp_centers(50.0), beam)
# 2: crb limit with N array
print("crb limit N arr", fisher.crb(p, [0, 0], np.array([100, 400]), zero_policy="limit"),
      [fisher.crb_limit(p, n) for n in (100, 400)], fisher.crb(p, [0, 0], np.array([100, 400])))
# extra: scalar r grid + scalar N
print("grid", fisher.crb(p, np.array([[0., 0.], [1., 0.]]), 100, zero_policy="limit"))
# N (2,) vs r (2,2) -> broadcast (2,)
print("paired", fisher.crb(p, np.array([[0., 0.], [3., 1.]]), np.array([100, 400]), zero_policy="limit"))
# 1: misalignment
t = time.time()
m = ex.misalignment_study([10.0], L=100.0, N=500, sbr=10, n_patterns=20, n_rep=200, positions=[(0.0, 0.0)])
for k in ("honest", "naive"):
    print(k, "sigma", m[k]["sigma_mean"][0, 0], "se", m[k]["sigma_se"][0, 0], "within", m[k]["sigma_se_within"][0, 0])
print("t", time.time() - t)
# 4: linear vectorial default
o = dict(polarization="linear")
bh = v.make_vectorial_beam(rho_max=300.0, **o)
be = v.make_vectorial_beam(mode="exact", **o)
c = patterns.tcp_centers(50.0)
ph, pe = photons.make_model(c, bh), photons.make_model(c, be)
print("lin crb(7,3)", float(fisher.crb(ph, np.array([7., 3.]), 100)), float(fisher.crb(pe, np.array([7., 3.]), 100)))
print("lin limit", fisher.crb_limit(ph, 100), fisher.crb_limit(pe, 100))
# pol_angle rotated
bh2 = v.make_vectorial_beam(rho_max=300.0, polarization="linear", pol_angle=0.7)
be2 = v.make_vectorial_beam(mode="exact", polarization="linear", pol_angle=0.7)
pts = np.random.default_rng(0).uniform(-200, 200, (40, 2))
print("pol_angle maxrel", np.max(np.abs(bh2(pts[:,0],pts[:,1]) - be2(pts[:,0],pts[:,1]))/np.max(be2(pts[:,0],pts[:,1]))))
# charge 2 linear
bh3 = v.make_vectorial_beam(rho_max=300.0, polarization="linear", charge=2)
be3 = v.make_vectorial_beam(mode="exact", polarization="linear", charge=2)
print("charge2 maxrel", np.max(np.abs(bh3(pts[:,0],pts[:,1]) - be3(pts[:,0],pts[:,1]))/np.max(be3(pts[:,0],pts[:,1]))))
# eps with linear
bh4 = v.make_vectorial_beam(rho_max=300.0, polarization="linear", eps=0.01)
be4 = v.make_vectorial_beam(mode="exact", polarization="linear", eps=0.01)
print("eps maxrel", np.max(np.abs(bh4(pts[:,0],pts[:,1]) - be4(pts[:,0],pts[:,1]))))
print("near 0", bh(np.array([0., 0.3, 1.0]), np.array([0., 0.2, 0.0])), be(np.array([0., 0.3, 1.0]), np.array([0., 0.2, 0.0])))
