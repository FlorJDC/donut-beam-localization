# -*- coding: utf-8 -*-
"""Newton polish (complex-step gradient, FD Hessian) of the noise-free naive MLE, eps=0.002."""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcore import tcp, probs
cen = tcp(50.0); FW = 300.0
for x in (10.0, 20.0):
    r = np.array([x, 0.0]); n = 100 * probs(r, cen, FW, eps=0.002, sbr=10.0)
    def grad(v):
        g = []
        for j in range(2):
            vv = np.asarray(v, complex).copy(); vv[j] += 1e-20j
            g.append(np.imag(np.sum(n * np.log(probs(vv, cen, FW, sbr=10.0)))) / 1e-20)
        return np.array(g)
    v = r.copy()
    for it in range(30):
        g = grad(v); h = 1e-4
        H = np.array([(grad(v + h * np.eye(2)[j]) - grad(v - h * np.eye(2)[j])) / (2 * h) for j in range(2)])
        v = v - np.linalg.solve(H, g)
    print(x, v - r, np.hypot(*(v - r)), "grad", grad(v), "N_eq", None)
