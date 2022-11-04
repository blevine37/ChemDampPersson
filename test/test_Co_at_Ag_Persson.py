import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

hbar    = 6.582119569e-16       # Plancks constant (eV.s) 
Ef      = 5.49e0                # Fermi energy in (eV). (CRC Handbook of Chemistry and Physics, 84th Ed, 2003-2004, Page 12-233)
delEaEf = 3.000e0               # ΔE = Ea - Ef (eV): CO induced resonance or virtual state 
Ea      = Ef + delEaEf          # Ea = CO induced resonance or virtual state (eV) 
spr     = 3.41                  # Observed position of Surface Plasmon Resonance of Ag@CO (eV)
Gama    = 1.0e0                 # Gamma Line width (eV) of CO induced virtual state  
e       = 3.794733e0            # Constant (eV^0.5 Å^0.5), e^2 = (e'^2/4πε) = (8.98755x10^9Nm^2C^-2)*(1.6022177x10^-19C)^2 = 14.4 eV Å
d       = 0.5e0                 # Distance b/w dynamic image plane and center of mass of Orbital (Å)  
neuF    = 1.39e16               # Fermi velocity of Ag (Å/s). (CRC Handbook of Chemistry and Physics, 84th Ed, 2003-2004, Page 12-233)
Q       = 0.33e0                # Q for px or py orbital of Oxygen (See Persson Page 156)
n       = 0.0586e0              # Carrier concentration 5.86×10^22 cm^-3 = 0.0586 Å^-3. (CRC Handbook of Chemistry and Physics, 84th Ed, 2003-2004, Page 12-233) 
na      = 0.1                   # No. of Adsorbate per unit surface area (Å^-2) (See Persson Page 158)
epi0    = 2.0                   # Bulk dielectric constant of CO matrix (See Persson Page 158)
r       = 3                     # rs = 3 jellium 
omgF    = Ef/hbar               # Fermi frequency i.e. Ef/hbar (s^-1)
hw      = np.linspace(0,8,1000) # Energy axis h.omega (eV))
 
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

outfile = open("Output.txt",'w',encoding='utf-8')

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

outfile.write("\n")
outfile.write(" Tangential Contribution to the Surface Plasmon Width:\n")
outfile.write(" +=============================+================+===========+\n")
outfile.write(" |         Parameter           |      Value     |   Units   |\n") 
outfile.write(" |=============================|================|===========|\n")
outfile.write(" {}{:12.5f}{}\n".format("|     σ0 (Persson Eq. 19)     | ",sig0,          "   |    Å^2    |"))
outfile.write(" {}{:12.5f}{}\n".format("|  J at SPR (Persson Eq. 17)  | ",JatSPR,        "   |    ---    |"))
outfile.write(" {}{:12.5f}{}\n".format("|            ΔA               | ",delAParll*1000,     "   |    meV.Å  |"))
outfile.write(" {}{:12.5f}{}\n".format("|     Δγ (Persson Eq. 15)     | ",gamalowerParll*1000,"   |    meV    |"))
outfile.write(" +=============================+================+===========+\n")

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

outfile.write("\n")
outfile.write(" Normal Contribution to the Surface Plasmon Width:\n")
outfile.write(" +=============================+================+===========+\n")
outfile.write(" |         Parameter           |      Value     |   Units   |\n") 
outfile.write(" |=============================|================|===========|\n")
outfile.write(" {}{:12.5f}{}\n".format("| Imα at SPR (Persson Eq. 25) | ",alpatSPR,        "   |    Å^3    |"))
outfile.write(" {}{:12.5f}{}\n".format("|            ΔA               | ",delAnorml*1000,     "   |   meV.Å   |"))
outfile.write(" {}{:12.5f}{}\n".format("|     Δγ (Persson Eq. 24)     | ",gamalowerNorml*1000,"   |    meV    |"))
outfile.write(" +=============================+================+===========+\n")

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
ax1.annotate('{} {:3.2f} {}'.format(r'$\mathrm{\epsilon_F}$ =',Ef,'eV'),xy=(0.4, 0.138))
ax1.annotate('{} {:3.2f} {}'.format(r'$\mathrm{\Gamma}$ =',Gama,'eV'),xy=(0.4, 0.126))
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())
plt.subplots_adjust(left=0.15)
plt.savefig("J.png", format='png', dpi=300)
# plt.show()

fig2, ax2 = plt.subplots()
ax2.plot(hw,alpha,color='b',label=str("{}{:4.3f}".format(r'$\mathrm{\epsilon_a} - \mathrm{\epsilon_F} = $',delEaEf)))
ax2.axvline(x = spr, color = 'black', linestyle="--",lw=1.0, label = 'Position of SPR')
ax2.set_xlim([0, 8])
ax2.set_xlabel("$\mathrm{\hbar\omega}$ (eV)",fontsize=12)
ax2.set_ylabel((str("Im α ($\mathrm{\AA}^3$)")),fontsize=12)
ax2.legend(frameon=False)
ax2.annotate('{} {:3.2f} {}'.format(r'$\mathrm{\epsilon_F}$ =',Ef,'eV'),xy=(0.4, 0.35))
ax2.annotate('{} {:3.2f} {}'.format(r'$\mathrm{\Gamma}$ =',Gama,'eV'),xy=(0.4, 0.32))
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())
plt.subplots_adjust(left=0.15)
plt.savefig("Im_alpha.png", format='png', dpi=300)
# plt.show()
print(" Done.\n")

a = input('Press a key to exit')
if a:
    exit(0)