import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
#print(np.__version__)


obs_spec = 'observed_spectrum4.txt'
syn_spec = 'synthetic_spectrum_5000-5500.spec'

# read data
w_obs, f_obs = np.loadtxt(obs_spec, unpack=True, comments='#')
w_syn, f_syn, e_syn = np.loadtxt(syn_spec, unpack=True, comments='#')

print('Wavelength range of your observed spectrum = [%f, %f]'%(w_obs.min(), w_obs.max()))
print('Wavelength range of your synthetic spectrum = [%f, %f]'%(w_syn.min(), w_syn.max()))
print('')
print('Enter observed wavelength range [wave_min, wave_max] for cross-correlation.')
print('It should be well within given synthetic wavelength range (+/- 15 Angstrom)')
w_min, w_max = eval(input('wave_min, wave_max = '))

# generate velocity range
v = np.linspace(-500.0, +500.0, 1001)	# km/s
c = 2.997e5	# km/s
n = len(v)
r_coef = np.empty(shape=n)		# empty array for cross-correlation coefficients

# mask the spectrum for given range
mask = np.where((w_obs >= w_min) & (w_obs <= w_max))

w_obs1 = w_obs[mask]
f_obs1 = f_obs[mask]

# cross-correlation
for i in range(n):
	w_syn1 = v[i]/c * w_syn + w_syn		#(w - w0) / w0 = v/c
	interp = interp1d(w_syn1, f_syn)
	w_syn2 = w_obs1
	f_syn2 = interp(w_syn2)
	crscor = pearsonr(f_obs1, f_syn2)
	r_coef[i] = crscor[0]
	#print(crscor)
	"""
	plt.plot(w_obs1, f_obs1, label='observed')
	plt.plot(w_syn2, f_syn2, label='synthetic')
	plt.legend()
	plt.show()
	plt.close()
	"""

# plotting ===== before RV correction =====
fig, ax = plt.subplots(3,1, figsize=(8,12))
ax[0].plot(w_obs1, f_obs1, color='r', label='observed (before RV correction)')
ax[0].plot(w_syn, f_syn, color='g', label='synthetic')
ax[0].legend()
ax[0].set_xlim(w_obs1.min(), w_obs1.max())
ax[0].set_xlabel('Wavelength')
ax[0].set_ylabel('flux')

# plotting ===== RV calculation =====
mask_rv = np.where(r_coef == r_coef.max())
radvel = v[mask_rv]
if (len(radvel)>1): radvel = radvel[0]
print('Radial Velocity = %6.2f km/s'%(radvel))
w_obs2 = -radvel/c * w_obs1 + w_obs1
f_obs2 = f_obs1

ax[1].plot(v, r_coef, color='k', label=str('RV = %6.2f'%(radvel)))
ax[1].legend()
ax[1].set_xlabel('velocity [km/s]')
ax[1].set_ylabel('cross-correlation coefficient')

# plotting ===== after RV correction =====
ax[2].plot(w_obs2, f_obs2, color='r', label='observed (after RV correction)')
ax[2].plot(w_syn, f_syn, color='g', label='synthetic')
ax[2].legend()
ax[2].set_xlim(w_obs1.min(), w_obs1.max())
ax[2].set_xlabel('Wavelength')
ax[2].set_ylabel('flux')


plt.show()
plt.close()
# END
