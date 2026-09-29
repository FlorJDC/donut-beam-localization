import numpy as np
from scipy.optimize import minimize_scalar, least_squares
fw=300.; a=4*np.log(2)/fw**2
def I(r2,eps=0.0): return (1-eps)*np.e*a*r2*np.exp(-a*r2)+eps  # eps: pedestal fraction of max (max of e*a r2 exp(-a r2) is 1)
def beams(L): 
    ang=np.deg2rad([90,210,330]); return np.vstack([np.c_[L/2*np.cos(ang),L/2*np.sin(ang)],[0,0]])
def lam(x,y,L,eps=0.0,b=0.0):
    B=beams(L); return I((x-B[:,0])**2+(y-B[:,1])**2,eps)
def crb(x,y,L,N,eps=0.0,sbr=0.0,exclude_zero=False,h=1e-4):
    # multinomial on p with background: counts ~ signal + bkg; SBR at centre defined via signal sum at center
    def p(x,y):
        s=lam(x,y,L,eps); 
        if sbr>0:
            s0=lam(0,0,L,eps).sum(); bk=s0/sbr/4; s=s+bk
        return s/s.sum()
    p0=p(x,y); dx=(p(x+h,y)-p(x-h,y))/2/h; dy=(p(x,y+h)-p(x,y-h))/2/h
    m=p0>1e-14 if exclude_zero else p0>0
    F=N*np.array([[np.sum(dx[m]**2/p0[m]),np.sum(dx[m]*dy[m]/p0[m])],[0,np.sum(dy[m]**2/p0[m])]]); F[1,0]=F[0,1]
    C=np.linalg.inv(F); return np.sqrt((C[0,0]+C[1,1])/2)
print("CRB centre limit (r=1e-3):",crb(1e-3,0,50,100,h=1e-5))
print("CRB centre pointwise S27:",crb(0,0,50,100,exclude_zero=True,h=1e-3))
print("ratio",crb(0,0,50,100,exclude_zero=True,h=1e-3)/crb(1e-3,0,50,100,h=1e-5),"sqrt5/2",np.sqrt(5)/2)
print("SBR10 centre:",crb(0,0,50,100,sbr=10,h=1e-3),crb(1e-3,0,50,100,sbr=10,h=1e-5))
for eps in [0.002,0.01,0.05]:
    r=minimize_scalar(lambda L: crb(0,0,L,100,eps=eps,h=1e-3),bounds=(2,150),method='bounded',options={'xatol':1e-3})
    print("eps",eps,"Lopt",r.x,"CRB",r.fun)
for tau in [4.0]:
    q=np.exp(-12.5/tau); print("M col:",[(1-q)*q**m/(1-q**4) for m in range(4)])
# L=100 mode1 p's
p=lam(25,0,100); print("p at x=L/4 (x-axis) ",p/p.sum()); print("ring radius",np.sqrt(1/a),np.sqrt(1/a)/fw)
