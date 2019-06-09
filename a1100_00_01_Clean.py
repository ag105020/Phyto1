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
a = genfromtxt('C:\\Users\\keiin\\Google Drive\\0 Paper\\3 phytoplanktonmodel\\Data\\Felcmanova 2017\\Felcmanova2017altogether.csv',delimiter=',')

Ddata = a[:,1]

LipidMean = a[:,2]
LipidSTDup = a[:,3]
LipidSTDdown = a[:,4]

ProteinMean = a[:,6]
ProteinSTDup = a[:,7]
ProteinSTDdown = a[:,8] 

CarbMean = a[:,10]
CarbSTDup = a[:,11]
CarbSTDdown = a[:,12]

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# Line preparation (Linear interpolation)
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

D = arange(0,0.5+0.1,0.1)
Lipid = -8.285*D + 9.8333  #Felcmanova2017altogether.xlsx
Protein = 74.264*D + 38.331
Carb = -58.401*D + 49.989

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# Preparing functrions
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

def Ebar2(ax,y,Mean,STDup,STDdown,Color,Label):
    ax.plot(D,y,dashes=(5,3),color=Color)
    ax.errorbar(Ddata,Mean,(STDdown,STDup),fmt='o'\
         ,color=Color,elinewidth=1,capthick=2,capsize=7,label=Label)
    xlim(0,0.4)
    xticks([0,0.1,0.2,0.3,0.4])

def Dpi():
    return 600

def SaveFolder():
    return 'Phytoplankton model\\Felcmanova2017'

#OOOOOOOOOOOOOOO
# Figures
#OOOOOOOOOOOOOOO


figure(1)
xlabel('$\mu$ (d$^{-1}$)')
ax1=gca(); ax2=ax1.twinx()
ax1.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2 
ax1.patch.set_visible(False) # hide the 'canvas' 

Ebar2(ax1,Carb,CarbMean,CarbSTDup,CarbSTDdown,'b','Carbohydrate')
Ebar2(ax1,Protein,ProteinMean,ProteinSTDup,ProteinSTDdown,'r','Protein')
Ebar2(ax2,D*nan,Ddata*nan,Ddata*nan,Ddata*nan,'r','Protein')
Ebar2(ax2,D*nan,Ddata*nan,Ddata*nan,Ddata*nan,'b','Carb.')
Ebar2(ax2,Lipid,LipidMean,LipidSTDup,LipidSTDdown,'g','Lipid')

title('Data: Prochlorococcus',y=1.02)
ax1.set_ylim(18,80)  #or (18,70)
ax2.set_ylim(4,17.23) #or (4,15.1)
ax1.set_ylabel('C allocation ($\%$)')
ax2.set_ylabel('C allocation ($\%$): Lipid')

ax2.legend(loc=2,fontsize=20)

Savefig2(SaveFolder(),'DataPlot',Dpi())

show()