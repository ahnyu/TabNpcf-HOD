from scipy import special
import numpy as np
import paircounts
import hodmodel as hm



###read HH HP PP
def read_XX(path,root,npoints,nbins):
    XX=np.empty((npoints,nbins,nbins))
    for i in range(nbins):
        for j in range(i+1):
            XX[:,i,j]=np.loadtxt(path+root+'_%02d' %i + '_%02d.dat' %j, usecols=(4))
    for i in range(nbins):
        for j in range(i+1,nbins):
            XX[:,i,j]=XX[:,j,i]
    return XX

def read_XY(path,root,npoints,nbins):
    XY=np.empty((npoints,nbins,nbins))
    for i in range(nbins):
        for j in range(nbins):
            XY[:,i,j]=np.loadtxt(path+root+'_%02d' %i + '_%02d.dat' %j, usecols=(4))
    return XY

###read HHH HPP PHH PPP
def read_XXX(path,root,npoints,nbins):
    XXX=np.empty((npoints,nbins,nbins,nbins))
    for i in range(nbins):
        for j in range(i+1):
            for k in range(j+1):
                XXX[:,i,j,k]=np.loadtxt(path+root+'_%02d' %(i) +'_%02d' %(j) + '_%02d.dat' %(k),usecols=(3))
    for i in range(nbins):
        for j in range(i+1,nbins):
            for k in range(i+1):
                XXX[:,i,j,k]=XXX[:,j,i,k]
    for i in range(nbins):
        for j in range(i+1,nbins):
            for k in range(i+1,nbins):
                XXX[:,i,j,k]=XXX[:,j,k,i]
    for i in range(nbins):
        for j in range(i+1):
            for k in range(j+1,nbins):
                XXX[:,i,j,k]=XXX[:,i,k,j]
    return XXX

def read_XYY(path,root,npoints,nbins):
    XYY=np.zeros((npoints,nbins,nbins,nbins))
    for i in range(nbins):
        for j in range(nbins):
            for k in range(j+1):
                XYY[:,i,j,k]=np.loadtxt(path+root+'_%02d' %(i) +'_%02d' %(j) + '_%02d.dat' %(k),usecols=(3))
    for i in range(nbins):
        for j in range(nbins):
            for k in range(j+1,nbins):
                XYY[:,i,j,k]=XYY[:,i,k,j]
    return XYY


#DD_weight
def DD_elg(hod,HH,HP,PP,bins,tab,model_e):
    DD_weight_HH=np.zeros((bins[0],bins[1]))
    DD_weight_HP=np.zeros((bins[0],bins[1]))    
    DD_weight_PP=np.zeros((bins[0],bins[1]))    
    DD_weight=np.zeros((bins[0],bins[1]))
    
    bin_mid=tab[:,0]
    num_halo=tab[:,1]
    num_part=tab[:,2]
    num_tab=len(bin_mid)

    w=hm.weight_elg(bin_mid,hod,model_e)
         
    for i in range(num_tab):
        wi_halo=w[0][i]
        wi_part=w[1][i]*num_halo[i]/num_part[i]
        for j in range(num_tab):
            wj_halo=w[0][j]
            wj_part=w[1][j]*num_halo[j]/num_part[j]
            
            tmp_HH=HH[:,i,j]
            tmp_HP=HP[:,i,j]            
            tmp_PP=PP[:,i,j]            
            
            DD_weight_HH+=(tmp_HH*wi_halo*wj_halo).reshape(bins[0],bins[1])
            DD_weight_HP+=(tmp_HP*wi_halo*wj_part).reshape(bins[0],bins[1]) 
            DD_weight_PP+=(tmp_PP*wi_part*wj_part).reshape(bins[0],bins[1])
    DD_weight=DD_weight_HH+2*DD_weight_HP+DD_weight_PP
    return DD_weight

def DD_qso(hod,HH,HP,PP,bins,tab,model_q):
    DD_weight_HH=np.zeros((bins[0],bins[1]))
    DD_weight_HP=np.zeros((bins[0],bins[1]))    
    DD_weight_PP=np.zeros((bins[0],bins[1]))    
    DD_weight=np.zeros((bins[0],bins[1]))
    
    bin_mid=tab[:,0]
    num_halo=tab[:,1]
    num_part=tab[:,2]
    num_tab=len(bin_mid)

    w=hm.weight_qso(bin_mid,hod,model_q)
    
    for i in range(num_tab):
        wi_halo=w[0][i]
        wi_part=w[1][i]*num_halo[i]/num_part[i]
        for j in range(num_tab):
            wj_halo=w[0][j]
            wj_part=w[1][j]*num_halo[j]/num_part[j]
            
            tmp_HH=HH[:,i,j]
            tmp_HP=HP[:,i,j]            
            tmp_PP=PP[:,i,j]            
            
            DD_weight_HH+=(tmp_HH*wi_halo*wj_halo).reshape(bins[0],bins[1])
            DD_weight_HP+=(tmp_HP*wi_halo*wj_part).reshape(bins[0],bins[1]) 
            DD_weight_PP+=(tmp_PP*wi_part*wj_part).reshape(bins[0],bins[1])
    DD_weight=DD_weight_HH+2*DD_weight_HP+DD_weight_PP
    return DD_weight

def DD_lrg(hod,HH,HP,PP,bins,tab,model_l):
    DD_weight_HH=np.zeros((bins[0],bins[1]))
    DD_weight_HP=np.zeros((bins[0],bins[1]))    
    DD_weight_PP=np.zeros((bins[0],bins[1]))    
    DD_weight=np.zeros((bins[0],bins[1]))
    
    bin_mid=tab[:,0]
    num_halo=tab[:,1]
    num_part=tab[:,2]
    num_tab=len(bin_mid)

    w=hm.weight_lrg(bin_mid,hod,model_l)
    
    for i in range(num_tab):
        wi_halo=w[0][i]
        wi_part=w[1][i]*num_halo[i]/num_part[i]
        for j in range(num_tab):
            wj_halo=w[0][j]
            wj_part=w[1][j]*num_halo[j]/num_part[j]
            
            tmp_HH=HH[:,i,j]
            tmp_HP=HP[:,i,j]            
            tmp_PP=PP[:,i,j]
            
            DD_weight_HH+=(tmp_HH*wi_halo*wj_halo).reshape(bins[0],bins[1])
            DD_weight_HP+=(tmp_HP*wi_halo*wj_halo*wj_part).reshape(bins[0],bins[1])
            DD_weight_PP+=(tmp_PP*wi_halo*wj_halo*wi_part*wj_part).reshape(bins[0],bins[1])

            
    DD_weight=DD_weight_HH+2*DD_weight_HP+DD_weight_PP
    return DD_weight

def DD_lxe(hod_lrg,hod_elg,HH,HP,PP,bins,tab,model_l,model_e):
    DD_weight_HH=np.zeros((bins[0],bins[1]))
    DD_weight_HP=np.zeros((bins[0],bins[1]))    
    DD_weight_PP=np.zeros((bins[0],bins[1]))    
    DD_weight=np.zeros((bins[0],bins[1]))
    
    bin_mid=tab[:,0]
    num_halo=tab[:,1]
    num_part=tab[:,2]
    num_tab=len(bin_mid)

    wl=hm.weight_lrg(bin_mid,hod_lrg,model_l)
    we=hm.weight_elg(bin_mid,hod_elg,model_e)
    
    for i in range(num_tab):
        wi_halo_l=wl[0][i]
        wi_part_l=wl[1][i]*num_halo[i]/num_part[i]
        wi_halo_e=we[0][i]
        wi_part_e=we[1][i]*num_halo[i]/num_part[i]
        for j in range(num_tab):
            wj_halo_l=wl[0][j]
            wj_part_l=wl[1][j]*num_halo[j]/num_part[j]
            wj_halo_e=we[0][j]
            wj_part_e=we[1][j]*num_halo[j]/num_part[j]
            
            tmp_HH=HH[:,i,j]
            tmp_HP=HP[:,i,j]            
            tmp_PP=PP[:,i,j]
            
            ##Halo-Halo
            DD_weight_HH+=0.5*(tmp_HH*(wi_halo_l*wj_halo_e+wi_halo_e*wj_halo_l)).reshape(bins[0],bins[1])
            
            ##Halo-Part
            DD_weight_HP+=(tmp_HP*(wi_halo_e*wj_halo_l*wj_part_l+wi_halo_l*wj_part_e)).reshape(bins[0],bins[1])
            
            ##Part-Part
            DD_weight_PP+=0.5*(tmp_PP*(wi_part_e*wj_halo_l*wj_part_l+wj_part_e*wi_halo_l*wi_part_l)).reshape(bins[0],bins[1])
            
    DD_weight=DD_weight_HH+DD_weight_HP+DD_weight_PP
    return DD_weight

def DD_lxq(hod_lrg,hod_qso,HH,HP,PP,bins,tab,model_l,model_q):

    DD_weight_HH=np.zeros((bins[0],bins[1]))
    DD_weight_HP=np.zeros((bins[0],bins[1]))    
    DD_weight_PP=np.zeros((bins[0],bins[1]))    
    DD_weight=np.zeros((bins[0],bins[1]))
    
    bin_mid=tab[:,0]
    num_halo=tab[:,1]
    num_part=tab[:,2]
    num_tab=len(bin_mid)

    wl=hm.weight_lrg(bin_mid,hod_lrg,model_l)
    wq=hm.weight_qso(bin_mid,hod_qso,model_q)
    
    for i in range(num_tab):
        wi_halo_l=wl[0][i]
        wi_part_l=wl[1][i]*num_halo[i]/num_part[i]
        wi_halo_q=wq[0][i]
        wi_part_q=wq[1][i]*num_halo[i]/num_part[i]
        for j in range(num_tab):
            wj_halo_l=wl[0][j]
            wj_part_l=wl[1][j]*num_halo[j]/num_part[j]
            wj_halo_q=wq[0][j]
            wj_part_q=wq[1][j]*num_halo[j]/num_part[j]
            
            tmp_HH=HH[:,i,j]
            tmp_HP=HP[:,i,j]            
            tmp_PP=PP[:,i,j]
            
            ##Halo-Halo
            DD_weight_HH+=0.5*(tmp_HH*(wi_halo_l*wj_halo_q+wi_halo_q*wj_halo_l)).reshape(bins[0],bins[1])
            
            ##Halo-Part
            DD_weight_HP+=(tmp_HP*(wi_halo_q*wj_halo_l*wj_part_l+wi_halo_l*wj_part_q)).reshape(bins[0],bins[1])
            
            ##Part-Part
            DD_weight_PP+=0.5*(tmp_PP*(wi_part_q*wj_halo_l*wj_part_l+wj_part_q*wi_halo_l*wi_part_l)).reshape(bins[0],bins[1])
            
    DD_weight=DD_weight_HH+DD_weight_HP+DD_weight_PP
    return DD_weight

def DD_exq(hod_elg,hod_qso,HH,HP,PP,bins,tab,model_e,model_q):
    
    DD_weight_HH=np.zeros((bins[0],bins[1]))
    DD_weight_HP=np.zeros((bins[0],bins[1]))    
    DD_weight_PP=np.zeros((bins[0],bins[1]))    
    DD_weight=np.zeros((bins[0],bins[1]))
    
    bin_mid=tab[:,0]
    num_halo=tab[:,1]
    num_part=tab[:,2]
    num_tab=len(bin_mid)

    we=hm.weight_elg(bin_mid,hod_elg,model_e)
    wq=hm.weight_qso(bin_mid,hod_qso,model_q)
    
    for i in range(num_tab):
        wi_halo_e=we[0][i]
        wi_part_e=we[1][i]*num_halo[i]/num_part[i]
        wi_halo_q=wq[0][i]
        wi_part_q=wq[1][i]*num_halo[i]/num_part[i]
        for j in range(num_tab):
            wj_halo_e=we[0][j]
            wj_part_e=we[1][j]*num_halo[j]/num_part[j]
            wj_halo_q=wq[0][j]
            wj_part_q=wq[1][j]*num_halo[j]/num_part[j]
            
            tmp_HH=HH[:,i,j]
            tmp_HP=HP[:,i,j]            
            tmp_PP=PP[:,i,j]
            
            ##Halo-Halo
            DD_weight_HH+=0.5*(tmp_HH*(wi_halo_e*wj_halo_q+wi_halo_q*wj_halo_e)).reshape(bins[0],bins[1])
            
            ##Halo-Part
            DD_weight_HP+=(tmp_HP*(wi_halo_q*wj_part_e+wi_halo_e*wj_part_q)).reshape(bins[0],bins[1])
            
            ##Part-Part
            DD_weight_PP+=0.5*(tmp_PP*(wi_part_q*wj_part_e+wj_part_q*wi_part_e)).reshape(bins[0],bins[1])
            
    DD_weight=DD_weight_HH+DD_weight_HP+DD_weight_PP
    return DD_weight


def wp2(DDw,RRw,bins):
    nrpbins=bins[0]
    npibins=bins[1]
    xirppi_weight=DDw/RRw - 1
    wpweight=np.zeros(nrpbins)
    dpi=1.0
    for rpbin in range(nrpbins):
        for piind in range(npibins):
            wpweight[rpbin]=wpweight[rpbin]+xirppi_weight[rpbin][piind]*dpi*2
    return wpweight

###DDD_weight
def DDD(hod,HHH,HPP,PHH,PPP,xi3,fit_common,massbin):

    nsbins=xi3['nsbins'][0]
    num_massbin_halo=massbin[0]
    num_massbin_part=massbin[1]
    massbin_mid=massbin[2]
    nbins=fit_common['nmassbins']


    DDD_weight_HHH=np.zeros(nsbins)
    DDD_weight_PHH=np.zeros(nsbins)
    DDD_weight_HPP=np.zeros(nsbins)
    DDD_weight_PPP=np.zeros(nsbins)
    DDD_weight=np.zeros(nsbins)

    w=hodmodel.weight(massbin_mid,hod)
    for i in range(nbins):
        wi_halo=w[0][i]
        wi_part=w[1][i]*num_massbin_halo[i]/num_massbin_part[i]

        for j in range(nbins):
            wj_halo=w[0][j]
            wj_part=w[1][j]*num_massbin_halo[j]/num_massbin_part[j]

            for k in range(nbins):
                wk_halo=w[0][k]
                wk_part=w[1][k]*num_massbin_halo[k]/num_massbin_part[k]

                tmp_HHH=HHH[:,i,j,k]
                tmp_PHH=PHH[:,i,j,k]
                tmp_HPP=HPP[:,i,j,k]
                tmp_PPP=PPP[:,i,j,k]

                DDD_weight_HHH+=tmp_HHH*wi_halo*wj_halo*wk_halo
                DDD_weight_PHH+=tmp_PHH*wi_part*wj_halo*wk_halo
                DDD_weight_HPP+=tmp_HPP*wi_halo*wj_part*wk_part
                DDD_weight_PPP+=tmp_PPP*wi_part*wj_part*wk_part

    DDD_weight=DDD_weight_HHH+3*DDD_weight_PHH+3*DDD_weight_HPP+DDD_weight_PPP

    return DDD_weight

def xi3(DDDw,RRRw):
    return DDDw/RRRw-1
