import numpy as np
from scipy.optimize import minimize
exec(open('/home/user/donut-beam-localization/equipo/2026-09-28_final-review/work/verify/indep.py').read().split('def crb')[0])
L=100.;N=500.;K=4
b=lam(0,0,L).sum()/(K*5)
def mu(x,y,bb): s=lam(x,y,L)+bb; return s/s.sum()
for x in [10.,50.]:
    pt=mu(x,0,b)
    f=lambda v: -np.sum(pt*np.log(mu(v[0],v[1],0)+1e-300))
    r=minimize(f,[x,0.5],method='Nelder-Mead',options={'xatol':1e-8,'fatol':1e-14,'maxiter':5000})
    # CRB b known / free (params x,y,b) multinomial
    def J(fun,v,h=1e-5):
        return np.array([(fun(*(v+h*e))-fun(*(v-h*e)))/(2*h) for e in np.eye(len(v))])
    F3=lambda v: N*(J(lambda *w: mu(w[0],w[1],w[2]),v)/pt)@J(lambda *w: mu(w[0],w[1],w[2]),v).T
    Fm=F3(np.array([x,0,b]))
    C2=np.linalg.inv(Fm[:2,:2]); C3=np.linalg.inv(Fm)
    print(x,"bias",r.x[0]-x,r.x[1],"CRBknown",np.sqrt((C2[0,0]+C2[1,1])/2),"CRBfree",np.sqrt((C3[0,0]+C3[1,1])/2))
