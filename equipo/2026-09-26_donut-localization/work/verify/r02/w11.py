# -*- coding: utf-8 -*-
import numpy as np
from w10 import run, sig
Ns = [250, 500, 1000, 2000, 4000, 8000]; S = []
for i, N in enumerate(Ns):
    Ls, errs = run(800, Ntot=N, seed=100+i)
    s, se = sig(errs[-1]); S.append(s)
    print(N, '%.4f ± %.4f' % (s, se), 'ratio cam %.4f' % (s/(100/np.sqrt(N))))
print('slope', np.polyfit(np.log(Ns), np.log(S), 1)[0])
