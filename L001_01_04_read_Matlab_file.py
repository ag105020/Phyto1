'''
Created on Apr 26, 2019
Light limited growth rate of Croco
@author: keiin
'''

from pylab import *
from FigSetting2 import *
from Savefig3 import *
from Savetxt2 import *
from Genfromtxt import *
from scipy.io import loadmat

Mat = loadmat('..\\Data\\HOT_data.mat')

DepthMat = Mat['depth']
NO3Mat = Mat['no3']
PO4Mat = Mat['po4']

rcParams.update({'ytick.right': 'True'})
rcParams.update({'ytick.direction': 'in'})
rcParams.update({'xtick.top': 'True'})
rcParams.update({'xtick.direction': 'in'})

#==================================
# Main Funcation
#==================================
z = arange(0,150+1,1) #(m) depth
def MuZ(a0,b0,c0):
    z0 = 30         #(m) Depth of 1/e light of surface
    I0 = 1000       #(umol m-2 s-1) Surface light intensity
    I = I0*exp(-z/z0) #(umol m-2 s-1)
    Mu = zeros(size(z)) 
    Qc = 2.4*10**(-13)*10**15 #(fmol C cell-1) from Masuda-san's info + Azoto stoichiometry
    
    a = a0 * 24/Qc #a*Mu + b (Kei 215-15~16) 24 is to convert from h-1 to d-1
    b = b0 * 24/Qc
    c = c0 * 24/Qc #c*MU (Kei 215-15~17)
    
    Pmax = 7. #(fmol C cell-1 h-1)
    #I0 = 100 #(umol m-2 s-1)
    I0P = 100
    
    P = Pmax*(1-exp(-I/I0P))
    
    Mu1 = (P-b)/(a+c)
    Mu2 = P/(2*c)
    
    Mu[Mu1<Mu2] = Mu1[Mu1<Mu2]
    Mu[Mu1>=Mu2] = Mu2[Mu1>=Mu2]
    
    Mu[Mu<0] = 0

    if a0 == 35.373:
        figure(1)
        plot(Mu,z,color='b')

        figure(2)
        plot(I,z,color='orange')
    
    elif a0 == 17.708:
        figure(1)
        plot(Mu,z,color='#00DCDC',dashes=(5,2))
    
    zNan = copy(z).astype(float)
    zNan[Mu==0]=nan
    zMax = max(zNan)
    
    return zMax,Mu
        

a00 = 35.373  #From C:\Users\keiin\OneDrive\Desktop\figures\01\HeteroC00\HighRes100\04 BoxSum0.xlsx run:617 04 23 00
b00 = 11.184
c00 = 78.847

a01 = 31.623  #From C:\Users\keiin\OneDrive\Desktop\figures\01\HeteroC00\HighRes101\BoxSum.xlsx run:617 04 23 01
b01 = 11.184
c01 = 62.192

a02 = 17.708  #From C:\Users\keiin\OneDrive\Desktop\figures\01\HeteroC00\HighRes103\BoxSum.xlsx run:617 04 23 03
b02 = 22.369
c02 = 78.847

zMax0,Mu0 = MuZ(a00,b00,c00)
#MuZ(a01,b01,c01)
zMax2,Mu2 = MuZ(a02,b02,c02)
#===================================
# Supporting funcation
#===================================
def ud(FigNumber):
    figure(FigNumber)
    gca().invert_yaxis()

#========================
# Plot preparation
#========================

Ylabel = 'Depth (m)'

Savefolder = '02\\05 LightLimitationCroco'
Ylim = (0,130)
Yticks = arange(0,120+20,20)

DPI = 100
#=======================
# Plotting
#=======================

figure(1)
#---For legend---
rcParams.update({'legend.fancybox': False})
plot([],[],color='#00DCDC',dashes=(5,2),label='Homo.')
plot([],[],color='b',label='Hetero.')
legend(loc=4)
#----------------
ylabel(Ylabel)
xlabel('$\mathit{\mu}$ (d$^{-1}$)')
ylim(Ylim)
yticks(Yticks)
xlim(left=0,right=0.5)
ud(1)
Savefig3(Savefolder,'Growth rate',DPI)
savetxt('..\\Output\\Model.csv',vstack((z,Mu0,Mu2)),delimiter=",",fmt='%s')

figure(2)
ylabel(Ylabel)
xlabel('Light ($\mu$mol m$^{-2}$ s$^{-1}$)')
ylim(Ylim)
yticks(Yticks)
xlim(left=0)
ud(2)
Savefig3(Savefolder,'Light',DPI)

HOT = genfromtxt('..\\Data\\HOT.csv',delimiter=',')[:17].T #[:17] is for up to 150 (m)

figure(3)
errorbar(HOT[3],HOT[0],xerr=(HOT[4],HOT[4]),fmt='o-',color='r',markeredgecolor='k',elinewidth=1,capthick=2,capsize=7)

ylabel(Ylabel)
ylim(Ylim)
xlabel('PO$_4^{3-}$ ($\mu$mol kg$^{-1}$)')


PO4range = arange(-0.002,0.302+0.01,0.01)
y1 = zMax0*ones(size(PO4range))
y2 = zMax2*ones(size(PO4range))
fill_between(PO4range, y1, y2, where=y1 >= y2, facecolor='#FDD2D0',edgecolor = "none")

xlim(0.025,0.17)
ud(3)
Savefig3(Savefolder,'PO4hot',DPI)

figure(33)
plot(PO4Mat,DepthMat,'o',color='k',markersize=1.5)
plot(HOT[3],HOT[0],'--',color='r',markeredgecolor='k',markersize=5)
xlim(-0.002,0.302)
ylabel(Ylabel)
ylim(Ylim)
xlabel('NO$_3^{-}$ ($\mu$mol kg$^{-1}$)')
fill_between(PO4range, y1, y2, where=y1 >= y2, facecolor='#FDD2D0',edgecolor = "none")
gca().invert_yaxis()
Savefig3(Savefolder,'PO4hotDot',DPI)

figure(4)
errorbar(HOT[1],HOT[0],xerr=(HOT[2],HOT[2]),fmt='o-',color='r',markeredgecolor='k',elinewidth=1,capthick=2,capsize=7)

ylabel(Ylabel)
ylim(Ylim)
xlabel('NO$_3^{-}$ ($\mu$mol kg$^{-1}$)')

NO3range = arange(0,3.0+0.01,0.01)
y1 = zMax0*ones(size(NO3range))
y2 = zMax2*ones(size(NO3range))
fill_between(NO3range, y1, y2, where=y1 >= y2, facecolor='#FDD2D0',edgecolor = "none")

xlim(0.0,1.1)
ud(4)
Savefig3(Savefolder,'NO3hot',DPI)

figure(34)
plot(NO3Mat,DepthMat,'o',color='k',markersize=1.5)
plot(HOT[1],HOT[0],'--',color='r',markeredgecolor='k',markersize=5)
xlim(0,3)
ylabel(Ylabel)
ylim(Ylim)
xlabel('NO$_3^{-}$ ($\mu$mol kg$^{-1}$)')
fill_between(NO3range, y1, y2, where=y1 >= y2, facecolor='#FDD2D0',edgecolor = "none")
gca().invert_yaxis()
Savefig3(Savefolder,'NO3hotDot',DPI)

show()

