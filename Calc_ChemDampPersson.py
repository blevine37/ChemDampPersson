import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

hbar    = 6.582119569e-16       # Plancks constant (eV.s) 
Ef      = 5.53e0                # Fermi energy in (eV). (CRC Handbook of Chemistry and Physics, 84th Ed, 2003-2004, Page 12-233)
delEaEf = 1.65e0                # Au-TiO2 Schottky barrier + (1/2)*Γ (eV). Paper (0.9-1.2 eV), old paper (1−1.5 eV). I use 1.25+0.5*0.8
Ea      = Ef + delEaEf          # TiO2 induced virtual state centered at Ea (eV)
spr     = 1.7712e0              # Observed position of Surface Plasmon Resonance of AuNR@TiO2 (eV). 700nm = 1.7712 eV 
Gama    = 0.8e0                 # Γ (eV) of TiO2. 1/2(FWHM of the TiO2 band). b/c 1/2 width is in conduction and 1/2 in valence band. 1/2(1.6 ev).
e       = 3.794733e0            # Constant (eV^0.5 Å^0.5), e^2 = (e'^2/4πε) = (8.98755x10^9Nm^2C^-2)*(1.6022177x10^-19C)^2 = 14.4 eV Å 
d       = 0.962e0               # Distance b/w dynamic image plane and COM of Orbital (Å) *(Read details below in Appendix A)  
neuF    = 1.40e16               # Fermi velocity of gold (Å/s). (CRC Handbook of Chemistry and Physics, 84th Ed, 2003-2004, Page 12-233)
Q       = 0.33e0                # Q for px or py orbital of Oxygen (See Persson Page 156)
n       = 0.059e0               # Carrier concentration 5.90×10^28 m^-3 = 0.059 Å^-3. (CRC Handbook of Chemistry and Physics, 84th Ed, 2003-2004, Page 12-233) 
na      = 0.1e0                 # No. of Adsorbate per unit surface area (Å^-2) [using Ag@CO value from Persson 1993]
epi0    = 60                    # Bulk dielectric constant of TiO2 matrix from 10.1016/S0040-6090(00)01027-0
r       = 152.78                # Radius (Å) of sphere having same mean free path as AuNR@TiO2 **(Read details below in Appendix B)
omgF    = Ef/hbar               # Fermi frequency i.e. Ef/hbar (s^-1)
hw      = np.linspace(0,8,1000) # Energy axis h.omega (eV)
 
#### Part-1: Calculate the tangential contribution to the surface plasmon width
# γ∥  = ΔA∥/R
# ΔA∥ = γ∥*R                      (Persson Page 158)
# γ∥  = (3/8)*(ν_F/R)*n_a*σ_diff  (Persson Eq. 15)
# ΔA∥ = (3/8)*ν_F*n_a*σ_diff
# ΔA∥ = (3/8)*ν_F*n_a*σ_0*J       (∵ σ_diff = σ_0*J Persson Eq. 18) 

# calculate σ_0 (Å^2) using Persson Eq. 19 (For factor 2, see Persson Page 158) 
sig0 = (64/(3*np.pi))*(omgF*Q/(n*neuF))*2

# Get J integral for plot using simplified Persson Eq. 17 and Eq. 16       
def integrandJ(E,Eaj,hwi):
    rhoe  = (1/np.pi)*((Gama/2)/((E-Eaj)**2 + (Gama/2)**2))
    rhoew = (1/np.pi)*((Gama/2)/((E+hwi-Eaj)**2 + (Gama/2)**2))
    return E*Gama*rhoew + (E+hwi)*Gama*rhoe 

Jint = []
for i in range(1,len(hw)): # Skip hw[0] to avoide divided by zero warning
    curJ = (np.pi/(4*Ef))*(1/hw[i])*quad(integrandJ, Ef-hw[i], Ef, args=(Ea,hw[i]))[0]    
    Jint.append(curJ)

# Get vlaue of J integral at SPR to use in Persson Eq. 18
JatSPR = (np.pi/(4*Ef))*(1/spr)*quad(integrandJ, Ef-spr, Ef, args=(Ea,spr))[0]

# Put values in Persson Eq. 18 to get ΔA∥ (Å/s)
delAParll = (3/8)*neuF*na*sig0*JatSPR

#Convert unit of ΔA∥ to eV.Å 
delAParll = delAParll*hbar 

# Divide by R to get Eq. 15 which gives γ∥ (eV)
gamalowerParll = delAParll/r

print("")
print(" Tangential Contribution to the Surface Plasmon Width:")
print(" +=============================+================+===========+")
print(" |         Parameter           |      Value     |   Units   |") 
print(" |=============================|================|===========|")
print(" {}{:12.5f}{}".format("|     σ0 (Persson Eq. 19)     | ",sig0,          "   |    Å^2    |"))
print(" {}{:12.5f}{}".format("|  J at SPR (Persson Eq. 17)  | ",JatSPR,        "   |    ---    |"))
print(" {}{:12.5f}{}".format("|            ΔA               | ",delAParll*1000,     "   |    meV.Å  |"))
print(" {}{:12.5f}{}".format("|     Δγ (Persson Eq. 15)     | ",gamalowerParll*1000,"   |    meV    |"))
print(" +=============================+================+===========+")

#### Part-2: Calculate the normal contribution to the surface plasmon width
# γ⊥  = ΔA⊥/R
# ΔA⊥ = γ⊥*R                            (Persson Page 158)
# γ⊥  = (16Pi/(1+2ε0))*(na/R)*εa*Im α⊥  (Persson Eq. 24)
# ΔA∥ = (16Pi/(1+2ε0))*na*εa*Im α⊥

# Get Im α⊥ integral for plot using Persson Eq. 25 and Eq. 16
def integrandalp(E,Eaj,hwi):
    rhoe  = (1/np.pi)*((Gama/2)/((E-Eaj)**2 + (Gama/2)**2))
    rhoew = (1/np.pi)*((Gama/2)/((E+hwi-Eaj)**2 + (Gama/2)**2))
    return rhoew*rhoe 

alpha = []
for i in range(len(hw)): # Here factor 2 account for both 2pi* orbitals (Persson Page 158) 
    ai = 2*(2*np.pi*(e*d)**2)*quad(integrandalp, Ef-hw[i], Ef, args=(Ea,hw[i]))[0]
    alpha.append(ai)

# Get vlaue of Im α⊥ integral at SPR to use in Persson Eq. 24
alpatSPR = 2*(2*np.pi*(e*d)**2)*quad(integrandalp, Ef-spr, Ef, args=(Ea,spr))[0]

# Put values in Persson Eq. 24 to get ΔA⊥ (Å/s)
delAnorml = (16*np.pi/(1+2*epi0))*na*spr*alpatSPR

# Divide by R to get Eq. 24 which gives γ∥ (eV)
gamalowerNorml = delAnorml/r

print("")
print(" Normal Contribution to the Surface Plasmon Width:")
print(" +=============================+================+===========+")
print(" |         Parameter           |      Value     |   Units   |") 
print(" |=============================|================|===========|")
print(" {}{:12.5f}{}".format("| Imα at SPR (Persson Eq. 25) | ",alpatSPR,        "   |    Å^3    |"))
print(" {}{:12.5f}{}".format("|            ΔA               | ",delAnorml*1000,     "   |   meV.Å   |"))
print(" {}{:12.5f}{}".format("|     Δγ (Persson Eq. 24)     | ",gamalowerNorml*1000,"   |    meV    |"))
print(" +=============================+================+===========+")
print("")
# print("{}{:12.8f}{}\n".format(" Total adsorbate induced contribution ΔA = ΔA∥+ΔA⊥ = ",delAParll+delAnorml," eV.Å"))
print(" Saving plots of J and Imα integrals in the current folder...")

plt.rc('font', family='serif', size=12)
fig1, ax1 = plt.subplots()
ax1.plot(hw[1:],Jint,color='r',label=str("{}{:4.3f}".format(r'$\mathrm{\epsilon_a} - \mathrm{\epsilon_F} = $',delEaEf)))
ax1.axvline(x = spr, color = 'black', linestyle="--",lw=1.0, label = 'Position of SPR')
ax1.set_xlim([0, 8])
ax1.set_xlabel("$\mathrm{\hbar\omega}$ (eV)",fontsize=12)
ax1.set_ylabel((str("J")),fontsize=12)
ax1.legend(frameon=False)
ax1.annotate('{} {:3.2f} {}'.format(r'$\mathrm{\epsilon_F}$ =',Ef,'eV'),xy=(4.83, 0.155))
ax1.annotate('{} {:3.2f} {}'.format(r'$\mathrm{\Gamma}$ =',Gama,'eV'),xy=(4.83, 0.135))
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())
plt.subplots_adjust(left=0.15)
plt.savefig("J_TiO2AuNR.png", format='png', dpi=300)
# plt.show()

fig2, ax2 = plt.subplots()
ax2.plot(hw,alpha,color='b',label=str("{}{:4.3f}".format(r'$\mathrm{\epsilon_a} - \mathrm{\epsilon_F} = $',delEaEf)))
ax2.axvline(x = spr, color = 'black', linestyle="--",lw=1.0, label = 'Position of SPR')
ax2.set_xlim([0, 8])
ax2.set_xlabel("$\mathrm{\hbar\omega}$ (eV)",fontsize=12)
ax2.set_ylabel((str("Im α ($\mathrm{\AA}^3$)")),fontsize=12)
ax2.legend(frameon=False)
ax2.annotate('{} {:3.2f} {}'.format(r'$\mathrm{\epsilon_F}$ =',Ef,'eV'),xy=(4.83, 2.4))
ax2.annotate('{} {:3.2f} {}'.format(r'$\mathrm{\Gamma}$ =',Gama,'eV'),xy=(4.83, 2.15))
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())
plt.subplots_adjust(left=0.15)
plt.savefig("Im_alpha_TiO2AuNR.png", format='png', dpi=300)
# plt.show()
print(" Done.\n")

a = input('Press a key to exit')
if a:
    exit(0)

### APPENDIX A ###
## * How I calculated the distance: ## 
# D1 = 2.040 Å (Distance between mean planes of two Au layers) 
# D2 = 1.982 Å (Distance between mean plane of outer Au layers and TiO2 layer adsorbed on Au)
# d = D2 - (D1/2) = 0.962 Å 

### APPENDIX B ###
## ** How I calculated the radius of the particle ##
# Reason: Our particle is rod, but Persson 1993 assumes that particle is spherical 
# SOlution: We calculated the mean free path of AuNR@TiO2 by using George C. Schatz 
# J. Chem. Phys. 119, 3926 (2003) method and then using that mean free path we shall 
# calculate rhe radius of a sphere which has same mean free path as AuNR@TiO2  
#
# According to Eq. 4.6 of Schatz 2003, effective mean free path (L_eff) of type II 
# (prolate) cylinders is: L_eff = 2d/(r+2) where d = Diameter, D = Height, r = d/D (aspect ratio)    
# AuNR@TiO2 has size 25 × 55 nm. So d = 25 nm, D = 55 nm, r = 0.4545 and L_eff = 20.3707 nm.  
# Now we shall calculate the radius of a sphere which has usng this L_eff we shall calculate 
# radius of sphere which has L_eff = 20.3707 nm. According to Schatz 2003 (Page 3929, Section III)
# for a sphere with radius R, L_eff = (4/3)R = 2d/3 => d = (3/2)L_eff = (3/2)* 20.3707 = 30.5561 nm
# So r = (30.5561/2) nm = 15.278 nm = 152.78 Å    
