# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 15:38:16 2022

@author: ML
"""

import numpy as np
import matplotlib.pyplot as plt

#User defined
Eps0=1.8   #Dielectric constant of medium
EpsAu=6.9   #Dielectric constant of nanoparticle material (gold)
aNP=40     #Radius of nanoparticle (nm)
r12=1       #Distance between electric charges in dipole (dye/protein)
R=40     #Centre-to-centre distance between nanoparticle and biomolecule (nm)
Gamma=40    #Full-width at medium height of dipol absorption peak (nm)
WavelengthB=300 #wavelength of the biomolecule absorption band taking part in the interaction
WavelengthLSPR=655 #wavelength of the nanoparticle LSPR band taking part in the interaction
RelaxAu=0.07 #Dielectric function relaxation factor


#Initializing constants - to SI units
Eps0=Eps0*8.854*10**-12
EpsAu=EpsAu*8.854*10**-12
Planck=6.58*10**-16 #in eV
aNP=aNP*10**-9
R=R*10**-9
r12=r12*10**-10
OmegaB=(2*np.pi*3*10**8)/(WavelengthB*10**-9)   #frequency of protein ABS band
OmegaLSPR=(2*np.pi*3*10**8)/(WavelengthLSPR*10**-9)     #frequency of gold LSPR band
Gamma=Gamma*OmegaB*Planck/(2*WavelengthB)
Beta=(aNP**3)*((EpsAu-Eps0)/(EpsAu+2*Eps0)) #nanoparticle (spherical!) coefficient
e=1.6*10**-19
me=9.109*10**-31
c=3*10**8
r=np.array([r12,r12,r12])   #initialize vector
p=np.array([r12,r12,r12], dtype=complex)    #electric-dipole transition 
p[0]=complex(r12,me*OmegaB*r[0])
p[1]=complex(0,me*OmegaB*r[1])
p[2]=complex(0,me*OmegaB*r[2])
u12=e*r #electric dipole
# m12=np.real(-e/(2*me*c)*np.cross(r,p))
m12=(0.05*10**-10)*OmegaB*r12*(-e/(2*c)) #magnetic dipole
m12=np.array([m12,0,0])
P=np.array([[1-Beta/R**3,0,0],[0,1-Beta/R**3,0],[0,0,1+2*Beta/R**3]]) #Field enhancement factor
Gomega=-((e**2)/((OmegaB**2) * (me**2) * (R**6) * Eps0))*Beta*(p[0]**2 + p[1]**2 + 4*p[2]**2) #dipole limit, bands shifts and broadening

#Calculations of CD molecule
wavelength=np.arange(200.0,800.0)
omega=(2*np.pi*3*10**8)/(wavelength*10**-9)
A1=(8/3)*np.sqrt(Eps0) #const
A2=Gamma/(np.abs(Planck*omega-Planck*OmegaB+np.complex(0,Gamma))**2) #const
A3=np.real(np.dot(np.dot(P,u12),m12))
CDmolecule=A1*A2*A3*omega


#Calculations of CD NP
# r=np.array([R,R,R])
# p=np.array([R,R,R], dtype=complex)
# p[0]=complex(1,me*OmegaLSPR*r[0])
# p[1]=complex(0,me*OmegaLSPR*r[1])
# p[2]=complex(0,me*OmegaLSPR*r[2])
# u12=e*r
# m12=-e/(2*me*c)*np.cross(r,p)
# P=np.array([[1-Beta/R**3,0,0],[0,1-Beta/R**3,0],[0,0,1+2*Beta/R**3]])
B1=(8/9)*(np.abs(3*Eps0/(EpsAu+2*Eps0))**2)*((aNP**3)/(np.sqrt(Eps0)*R**3))
# B2=(OmegaLSPR**2)/(omega*(omega**2 + RelaxAu**2))
B2=EpsAu
B3=((np.dot(u12[0],m12[0])+np.dot(u12[1],m12[1])-2*np.dot(u12[2],m12[2]))/(np.abs(Planck*omega-Planck*OmegaLSPR+np.complex(0,Gamma))))
CDnp=B1*B2*B3*omega

CD=CDmolecule+CDnp
CD=(CDmolecule+CDnp)*(2*np.pi*(6*10**23)*(10**-4))/(0.23*(3*10**8)*np.sqrt(Eps0))

#plot
plt.figure()
# plt.plot(wavelength,CDmolecule)
# plt.plot(wavelength,CDnp)
plt.plot(wavelength,CD)
# plt.title('Protein band width= '+str(int(GammaW))+'nm')

#1 Brakuje odpowiednich stałych aby policzyć rzeczywiste wartosci sprzezenia  dipolowego
#2 Brakuje w modelu opisu anizotropowych nanoczastek ktore posiadaja tez plazmon poprzeczny

