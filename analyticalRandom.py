import numpy as np

###analytical RR
def RR(N,wp2):
    boxsize=2000.0
    rpbins=np.logspace(wp2['logrpmin'][0],wp2['logrpmax'][0],wp2['nrpbins'][0]+1)
    pibins=np.linspace(0,wp2['pimax'][0],wp2['npibins'][0]+1)
    boxvec=np.array([boxsize,boxsize,boxsize])
    V = np.pi * np.outer(rpbins**2, (2.*pibins))
    dV = np.diff(np.diff(V, axis=0), axis=1)
    dVd=dV/boxvec.prod()
    r1r2=N**2*dVd
    return r1r2

###analytical RRR
def sphereOverlapVolume(d,R,r):
    V=0
    if(r>R):
        r,R=R,r
    if(d<(R+r)):
        if(d>(R-r)):
            V=(np.pi*(R+r-d)*(R+r-d)*(d*d+2.0*d*r-3.0*r*r+2.0*d*R+6.0*r*R-3.0*R*R))/(12.0*d)
        else:
            V=(4.0*np.pi/3.0)*r*r*r
    return V

def crossSectionVolume(r1,r2,r3):
    V_oo = sphereOverlapVolume(r1,r3+0.5*delta_r,r2+0.5*delta_r)
    V_oi = sphereOverlapVolume(r1,r3+0.5*delta_r,r2-0.5*delta_r)
    V_io = sphereOverlapVolume(r1,r3-0.5*delta_r,r2+0.5*delta_r)
    V_ii = sphereOverlapVolume(r1,r3-0.5*delta_r,r2-0.5*delta_r)
    return V_oo-V_oi-V_io+V_ii

def getPermutations(r1,r2,r3):
    perm=1
    if(r1!=r2 and r1!=r3 and r2!=r3):
        perm=6
    elif((r1==r2 and r1!=r3) or (r1==r3 and r1!=r2) or (r2==r3 and r2!=r1)):
        perm=3
    return perm

def sphericalShellVolume(r):
    r_o=r+0.5*delta_r
    r_i=r-0.5*delta_r
    return 4*np.pi*(r_o*r_o*r_o-r_i*r_i*r_i)/3.0

def gaussQuadCrossSection(r1,r2,r3):
    result=0.0
    for i in range(len(w)):
        r_1=r1+0.5*delta_r*x[i]
        result+=0.5*delta_r*w[i]*crossSectionVolume(r_1,r2,r3)*r_1*r_1
    return result

def RRR(Nrand,Nshell,rbin,V_box):
    RRR_analy=[]
    nbar=Nrand/V_box
    for i in range(Nshell):
        for j in range(i,Nshell):
            for k in range(j,Nshell):
                if(rbin[k]<=(rbin[i]+rbin[j])):
                    ind=k+Nshell*(j+Nshell*i)
                    V=gaussQuadCrossSection(rbin[i],rbin[j],rbin[k])
                    n_perm=getPermutations(rbin[i],rbin[j],rbin[k])
                    RRR_analy.append(4*np.pi*n_perm*nbar*nbar*V*Nrand)
    return np.array(RRR_analy)
