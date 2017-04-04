# -*- coding: utf-8 -*-
import numpy as np

########################## Constants ############################
Rsun = 6.957e10         # Radius of the Sun (cm)
Rj = 7149200000.0       # Radius of Jupiter (cm)
kB = 1.38064852e-16     # Boltzman's constant (in erg/K)
gamma = 0.57721         # Euler-Mascheroni constant (dimensionless)
mamu = 1.660539040e-24  # Atomic mass unit (grams)
mH = 1.00794            # Molecular weight of Hydrogen (amus)
mHe = 4.002602          # Same for Helium (amus)
mH2O = 18.01528         # Same for water (amus)
mCH4 = 16.0425          # Same for methane (amus)
mCO2 = 44.0095          # Same for carbon dioxyde (amus)
mCO = 28.0101           # Same for carbon monoxyde (amus)
mTiO = 63.8664          # Same for titaniun oxyde (amus)


#################################################################

####################### Planet inputs ###########################
R0 = 1.20*Rj            # Reference radius, in cm
P0 = 10.                # Reference pressure, in bars
Rs = 1.038*Rsun          # Stellar radii, in cm
g = 713.               # Planetary gravity, in cm/s^2
T = 1399.               # Temperature in the atmosphere, in K
sigma_cloud = 0.0       # If you have an oppacity for the clouds, add it here.
#################################################################

# Mixing ratios for T = 1250 K, solar composition HJ:
"""
H2 = 0.83613
He = 0.16228
H2O = 0.00096782
CH4 = 0.00041993
CO2 = 3.0772e-08
CO = 4.1632e-05
TiO = 3.9358e-08
"""
# Mixing ratios for T = 2087 K, solar composition HJ atmosphere:
H2 = 0.83593
He = 0.16193
H2O = 5.4885e-4
CH4 = 1.439e-6
CO2 = 6.1649e-8
CO = 4.5939e-4
TiO = 1.4848e-7
# Mixing raitos for T = 2087 K, C/O = 1 solar composition HJ atmosphere:
"""
H2 = 0.83596
He = 0.162
H2O = 3.9812e-05
CH4 = 4.1815e-05
CO2 = 9.4243e-09
CO = 0.00096819
TiO = 9.8804e-08
"""
# Mixing raitos for T = 2087 K, C/O = 0.1 solar composition HJ atmosphere:
"""
H2 = 0.83587
He = 0.16213
H2O = 0.00090784
CH4 = 1.9079e-07
CO2 = 2.237e-08
CO = 0.00010077
TiO = 1.4808e-07
"""
# Conver bars to cgs:
P0 = P0*1e6          # 1 bar = 1e6 g/cm/s^2

# Import opacities, index them:
Top = '2000'
min_nu = 5000.
max_nu = 30000.
nu = np.arange(min_nu,max_nu,0.1)
wavelength = (1./nu)
for molecule in ['H2O','CH4','CO2','CO','TiO']:
    # Import and put everything on a common wavenumber scale:
    exec "sigma_"+molecule+" = np.zeros(len(nu))"
    exec "nu_"+molecule+",csigma_"+molecule+" = np.loadtxt('cross_sections/"+molecule+"_"+Top+".sigma',unpack=True)"
    exec "c_min_wav = np.min(nu_"+molecule+")"
    exec "c_max_wav = np.max(nu_"+molecule+")"
    exec "idx = np.where((nu>=c_min_wav)&(nu<=c_max_wav))[0]"
    exec "sigma_"+molecule+"[idx] = csigma_"+molecule
    # Invert axis to follow wavelength indexing defined after this for loop:
    exec "sigma_"+molecule+" = sigma_"+molecule+"[::-1]"
    print "Loaded "+molecule+" cross-section."

# Invert axis in order to have increasing wavelengths:
nu = nu[::-1]
wavelength = wavelength[::-1]

# Compute total cross-section:
sigma = 0.0
for molecule in ['H2O','CH4','CO2','CO','TiO']:
    exec "sigma = "+molecule+"*sigma_"+molecule+" + sigma"

# Add cloud cross-section:
sigma = sigma + sigma_cloud

# Keep computations only where sigma != 0 (i.e., were there is 
# at least one opacity source):
idx = np.where(sigma != 0)[0]
wavelength = wavelength[idx]
sigma = sigma[idx]
print wavelength

# Calculate mean molecular mass (note we assume a molecular 
# hydrogen-dominated, H2, atmopshere):
m = (H2*(2.*mH) + He*mHe + H2O*mH2O + CH4*mCH4 + CO2*mCO2 + CO*mCO + TiO*mTiO)/(H2+He+H2O+CH4+CO2+CO+TiO)

print 'Mean molecular mass: ',m
# Compute scale-height:
H = (kB*T)/(m*mamu*g)

# Get reference optical depth:
tau0 = (P0*sigma/(kB*T)) * np.sqrt(2.*np.pi*H*R0)

# Get radius as a function of wavelength:
R = R0 + H*(gamma + np.log(tau0))

# Bin model:
idx = range(0,len(wavelength),100)
w_bin = np.zeros(len(idx))
depth_bin = np.zeros(len(idx))
model_depth = (R/Rs)**2
for k in range(len(idx)):
    if k != len(idx)-1:
        w_bin[k] = np.mean(wavelength[idx[k]:idx[k+1]])   
        depth_bin[k] = np.mean(model_depth[idx[k]:idx[k+1]])
    else:
        w_bin[k] = np.mean(wavelength[idx[k]:])
        depth_bin[k] = np.mean(model_depth[idx[k]:])

# Plot it:
from scipy.ndimage.filters import gaussian_filter1d
import matplotlib.pyplot as plt
plt.style.use('ggplot')
plt.plot(wavelength*1e4,R/Rj)
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('Radius ($R_J$)')
plt.show()

plt.plot(wavelength*1e4,model_depth*1e6,label='Full model')
plt.plot(wavelength*1e4,gaussian_filter1d(model_depth*1e6,50),label = 'Gaussian filter')
plt.plot(w_bin*1e4,depth_bin*1e6,label='Mean binning')
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('$(R_p/R_s)^2$ (ppm)')
plt.xlim([1.15,1.65])
plt.legend()
#plt.ylim([1.32,1.43])
plt.show()

# Save model:
f = open('transpec_cross_section_out.dat','w')
for i in range(len(w_bin)):
    f.write(str((w_bin[i]*1e4))+' '+str(depth_bin[i])+'\n')
f.close()
