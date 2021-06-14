import matplotlib.pyplot as plt
from scipy import special
import numpy
import cmath
from mpmath import coulombf


Bess = []
BessAbs = []
BessRe = []
BessIm = []

Coul = []
CoulAbs = []
CoulRe = []
CoulIm = []

Momentum = 100
RMAX = 10
NR = 10
FmToNu = 5.067731237e-3
AlphaFS = 0.0072973525664

for RAD in numpy.arange(0,RMAX,RMAX/NR):
    RHO = Momentum*RAD*FmToNu;
    l = 1
    z = complex(RHO,0.5);
    bz = special.jv(l,z)
    Bess.append(special.jv(1,RHO))
    BessAbs.append(abs(bz))
    BessRe.append(bz.real)
    BessIm.append(bz.imag)


    RedMass = 500;
    q1q2 = 1;
    eta = RedMass*q1q2*AlphaFS/Momentum;
    cz = coulombf(l,eta,z);
    Coul.append(coulombf(l,eta,RHO))
    CoulAbs.append(abs(cz))
    CoulRe.append(cz.real)
    CoulIm.append(cz.imag)
    print(cz)

plt.figure(1)
plt.plot(numpy.arange(0,RMAX,RMAX/NR),Bess,numpy.arange(0,RMAX,RMAX/NR),BessAbs,numpy.arange(0,RMAX,RMAX/NR),BessRe,numpy.arange(0,RMAX,RMAX/NR),BessIm)
plt.ylabel('some numbers')

plt.figure(2)
plt.plot(numpy.arange(0,RMAX,RMAX/NR),Coul,numpy.arange(0,RMAX,RMAX/NR),CoulAbs,numpy.arange(0,RMAX,RMAX/NR),CoulRe,numpy.arange(0,RMAX,RMAX/NR),CoulIm)
plt.ylabel('some numbers')


plt.show()
