'''
Created on May 10, 2019

@author: keiin
'''

from pylab import *
from Savefig2 import *
from FigSetting2 import *
import matplotlib.ticker as ptick

#OOOOOOOOOOOOOOOOOOOOOOOOOOOO
# Data preparation
#OOOOOOOOOOOOOOOOOOOOOOOOOOOO
a = genfromtxt('C:\\Users\\keiin\\Google Drive\\0 Paper\\3 phytoplanktonmodel\\Data\\Jahn 2018\\Jahn2018Altogether.csv',delimiter=',')

MuLHC = a[:,1]
LHC = a[:,2]

MuRIB = a[:,4]
RIB = a[:,5]

MuCBM = a[:,7]
CBM = a[:,8]


#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# Line preparation (Linear interpolation)
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

MuFit = arange(0,0.12+0.01,0.01)
LHCfit = -1.7718*MuFit + 0.3152  #Jahn2018Altogether.xlsx
RIBfit = 0.9507*MuFit + 0.1012
CBMfit = 0.7206*MuFit + 0.0921

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# Preparing functrions
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO



def Dpi():
    return 600

def SaveFolder():
    return 'Phytoplankton model\\Jahn2018'

#OOOOOOOOOOOOOOO
# Figures
#OOOOOOOOOOOOOOO

ColorLHC = 'orange'
ColorRIB = '#00ff00'
ColorCBM = 'g'

# figure(1)
# plot(MuFit*24,LHCfit,color = ColorLHC)
# plot(MuLHC*24,LHC,'o',color = ColorLHC)
# ylim(ymin=0)
# 
# figure(2)
# plot(MuFit*24,RIBfit,color = ColorRIB)
# plot(MuRIB*24,RIB,'o',color = ColorRIB)
# ylim(ymin=0)
# 
# figure(3)
# plot(MuFit*24,CBMfit,color = ColorCBM)
# plot(MuCBM*24,CBM,'o',color = ColorCBM)
# ylim(ymin=0)

Dashes=(5,3)

figure(4)
plot(MuFit*24,LHCfit,dashes=Dashes,color = ColorLHC,label='LHC')
plot(MuLHC*24,LHC,'o',color = ColorLHC)

plot(MuFit*24,RIBfit,dashes=Dashes,color = ColorRIB,label='RIB')
plot(MuRIB*24,RIB,'o',color = ColorRIB)

plot(MuFit*24,CBMfit,dashes=Dashes,color = ColorCBM,label='CBM')
plot(MuCBM*24,CBM,'o',color = ColorCBM)

xlabel('$\mu_{max}$ (d$^{-1}$)')
ylabel('Protein mass fraction')
ylim(ymin=0)
title('Data: Synechocystis',y=1.02)
legend(loc=4,ncol=2,fontsize=20)
Savefig2(SaveFolder(),'DataPlot',Dpi())



# figure(1)
# xlabel('$\mu$ (d$^{-1}$)')
# ax1=gca(); ax2=ax1.twinx()
# ax1.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2 
# ax1.patch.set_visible(False) # hide the 'canvas' 
# 
# Ebar2(ax1,Carb,CarbMean,CarbSTDup,CarbSTDdown,'g','Carbohydrate')
# Ebar2(ax1,Protein,ProteinMean,ProteinSTDup,ProteinSTDdown,'r','Protein')
# Ebar2(ax2,D*nan,Ddata*nan,Ddata*nan,Ddata*nan,'r','Protein')
# Ebar2(ax2,D*nan,Ddata*nan,Ddata*nan,Ddata*nan,'g','Carb.')
# Ebar2(ax2,Lipid,LipidMean,LipidSTDup,LipidSTDdown,'b','Lipid')
# 
# title('Prochlorococcus',y=1.02)
# ax1.set_ylim(18,80)  #or (18,70)
# ax2.set_ylim(4,17.23) #or (4,15.1)
# ax1.set_ylabel('C allocation ($\%$)')
# ax2.set_ylabel('C allocation ($\%$): Lipid')
# 
# ax2.legend(loc=2,fontsize=23)
# 
# Savefig2(SaveFolder(),'DataPlot',Dpi())

show()