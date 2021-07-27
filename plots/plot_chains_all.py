import emcee
import numpy as np
import corner
filename='../test.h5'
reader=emcee.backends.HDFBackend(filename)
#tar=reader.get_autocorr_time()
burn=3000
samples_all=reader.get_chain(flat=True)
samples_burn=reader.get_chain(flat=True,discard=burn)
log_pro_burn=reader.get_log_prob(discard=burn,flat=True)
log_pro_all=reader.get_log_prob(flat=True)
#labels=[r'$log_{M_{cut}}^{lrg}$',r'$log_{M_1}^{lrg}$',r'$\sigma^{lrg}$',r'$\alpha^{lrg}$',r'$\kappa^{lrg}$',r'$p_{max}^{lrg}$',r'$log_{M_{cut}}^{elg}$',r'$log_{M_1}^{elg}$',r'$\sigma^{elg}$',r'$\alpha^{elg}$',r'$\kappa^{elg}$',r'$p_{max}^{elg}$',r'$\gamma^{elg}$',r'$Q^{elg}$',r'$log_{M_{cut}}^{qso}$',r'$log_{M_1}^{qso}$',r'$\sigma^{qso}$',r'$\alpha^{qso}$',r'$\kappa^{qso}$',r'$p_{max}^{qso}$',r'$loglike$']
labels=[r'$log_{M_{cut}}^{lrg}$',r'$log_{M_1}^{lrg}$',r'$\sigma^{lrg}$',r'$\alpha^{lrg}$',r'$\kappa^{lrg}$',r'$p_{max}^{lrg}$',r'$log_{M_{cut}}^{elg}$',r'$log_{M_1}^{elg}$',r'$\sigma^{elg}$',r'$\alpha^{elg}$',r'$p_{max}^{elg}$',r'$\gamma^{elg}$',r'$Q^{elg}$',r'$log_{M_{cut}}^{qso}$',r'$log_{M_1}^{qso}$',r'$\sigma^{qso}$',r'$\alpha^{qso}$',r'$\kappa^{qso}$',r'$p_{max}^{qso}$',r'$loglike$']
#labels=[r'$log_{M_{cut}}^{lrg}$',r'$log_{M_1}^{lrg}$',r'$\sigma^{lrg}$',r'$\alpha^{lrg}$',r'$\kappa^{lrg}$',r'$p_{max}^{lrg}$',r'$log_{M_{cut}}^{elg}$',r'$log_{M_1}^{elg}$',r'$\sigma^{elg}$',r'$\alpha^{elg}$',r'$p_{max}^{elg}$',r'$\gamma^{elg}$',r'$Q^{elg}$',r'$log_{M_{cut}}^{qso}$',r'$log_{M_1}^{qso}$',r'$\sigma^{qso}$',r'$\alpha^{qso}$',r'$\kappa^{qso}$',r'$p_{max}^{qso}$']
all_all=np.concatenate((samples_all,log_pro_all[:,None]),axis=1)
all_burn=np.concatenate((samples_burn,log_pro_burn[:,None]),axis=1)
theta_max=samples_all[np.argmax(reader.get_log_prob(flat=True))]
theta_max_all=np.append(theta_max,[np.amax(reader.get_log_prob(flat=True))])
print(theta_max_all)
fig=corner.corner(all_burn,show_titles=True,labels=labels,plot_datapoints=True,quantiles=[0.025,0.16,0.5,0.84,0.975],truths=theta_max_all)
fig.savefig('test.png')
