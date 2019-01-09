'''
Created on Aug 15, 2018
Here I read Mu max files from three different runs (Healey85, Challup90, Sakshaug89)
and plot them with data.
I add Pchl model

Originally from a1000 00 02
Double axes from a822 04 02

@author: Keisuke
'''
from pylab import *
from Savefig2 import *
from FigSetting import *
import matplotlib.ticker as ptick


#OOOOOOOOOOOOOOOOOOOOOOOOO
# Define functions
#OOOOOOOOOOOOOOOOOOOOOOOOO
def Dpi():
    return 600

def SaveFolder():
    return 'Phytoplankton model\\N-limiting\\MuMax'

def ReadData(folder,name):
    return genfromtxt(folder+name,delimiter=',')

def CreateI(MuMax,step):  #create I array from Mumax with input of step
    return arange(step,size(MuMax)*step+step,step)
    
def Plot(FigNumber,x1,y1,x2,y2,x3,y3,Title):
    figure(FigNumber)
    plot(x1,y1);plot(x2,y2);plot(x3,y3)
    title(Title)

def Plot2(FigNumber,X,Y,C,Title,Label):
    figure(FigNumber)
    for i in arange(shape(X)[0]):
        plot(X[i],Y[i],color=C[i],label=Label[i])
        title(Title)

def Plot1(FigNumber,Title,xData,yData,cData,xModel,yModel,cModel,xPchl,yPchl,cPchl,Loc):
    figure(str(FigNumber),figsize=(8.5,6.5))
    ax1=gca()
    ax2=ax1.twinx()
    ln1 = ax1.plot(xData,yData,'o',color=cData,label='$\mu_{max}$ data',zorder=10)
    ln2 = ax1.plot(xModel,yModel,color=cModel,label='$\mu_{max}$ model',zorder=5)
    ln3 = ax2.plot(xPchl,yPchl,'--',color=cPchl,label='$P_{Chl}$ $\ $model',zorder=1)
    ax1.set_ylim(0,1.8)
    ax2.set_ylim(0,0.003)
    title(Title)
    ax1.set_xlabel('Irradiance \\textrm{\\greektext m}mol m$^{-2}$ s$^{-1}$')
    ax1.set_ylabel('$\mu_{max}$ (d$^{-1}$)')
    ax2.set_ylabel('$P_{Chl}$ (mol C h$^{-1}$ Chl mol C$^{-1}$)')
    ax2.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
    lns = ln1 + ln2 + ln3
    labs = [l.get_label() for l in lns]
    legend(lns, labs, loc=Loc,fontsize=27)
    xlim(0,200)
   
def DataPlot(X,Y,C):
    for i in arange(shape(X)[0]):
        plot(X[i],Y[i],'o',color=C[i])

#OOOOOOOOOOOOOOOOOOOOOOOOOOO
# Define Loading locations
#OOOOOOOOOOOOOOOOOOOOOOOOOOO
LoadFolderHealey85 = 'C:\\Users\\Keisuke\\Desktop\\figures\\Phytoplankton model\\N-limiting\\MuMax\\'
LoadFolderChalup90 = 'C:\\Users\\Keisuke\\Desktop\\figures\\Phytoplankton model\\N-limiting\\MuMax\\'
LoadFolderSakshaug89 = 'C:\\Users\\Keisuke\\Desktop\\figures\\Phytoplankton model\\N-limiting\\MuMax\\'

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# Load files (model results)
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOO
Step = 0.1

MuMaxModelHealey85 = ReadData(LoadFolderHealey85,'Healey85.csv')
IModelHealey85 = CreateI(MuMaxModelHealey85,Step) 

MuMaxModelChalup90 = ReadData(LoadFolderChalup90,'Chalup90.csv')
IModelChalup90 = CreateI(MuMaxModelChalup90,Step) 

MuMaxModelSakshaug89 = ReadData(LoadFolderSakshaug89,'Sakshaug89.csv')
IModelSakshaug89 = CreateI(MuMaxModelSakshaug89,Step)

MuMaxModelTuple = (MuMaxModelHealey85,MuMaxModelChalup90,MuMaxModelSakshaug89)
IModelTuple = (IModelHealey85,IModelChalup90,IModelSakshaug89)


#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# Preparing Data for plotting
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

IDataHealey85 = array([12,22,38,62,144])
MuMaxDataHealey85 = genfromtxt('Healey_data_Dmax_for_AW_method_N_limited.csv',delimiter=',')

IDataChalup90 = (62.96,188.66)
MuMaxDataChalup90 = (0.57,1.14)

IDataSakshaug89 = (12.03,70.56,99.44)
MuMaxDataSakshaug89 = (0.53,1.2,1.4)

IDataTuple = (IDataHealey85,IDataChalup90,IDataSakshaug89)
MuMaxDataTuple = (MuMaxDataHealey85,MuMaxDataChalup90,MuMaxDataSakshaug89)

#OOOOOOOOOOOOOOOOOOOOOOOOOOO
#Preparing Pchl
#OOOOOOOOOOOOOOOOOOOOOOOOOOO
I = arange(0.1,200+0.1,0.1)

#Healey 85
PmaxH85 = 0.0029241122800652097
OTH85 = 0.010874458727716233
PchlH85 = PmaxH85*(1-exp(-OTH85*I))

#Chalup90
PmaxC90 = 0.0025 #from a822 20 00
OTC90 = 0.010874458727716233/2 #from a822 20 00
PchlC90 = PmaxC90*(1-exp(-OTC90*I))

#Sakshaug89
PmaxS89 = 0.005 #from a822 30 00
OTS89 = 0.010874458727716233/1.8 #from a822 30 00
PchlS89 = PmaxS89*(1-exp(-OTS89*I))

#OOOOOOOOOOOOOOOOOOOOOOOOOOO
# Plot
#OOOOOOOOOOOOOOOOOOOOOOOOOOO

Plot1('Healey85','$S.$ $linearis$',IDataHealey85,MuMaxDataHealey85,'b',IModelHealey85,MuMaxModelHealey85,'b',I,PchlH85,'#00FFFF','lower right')
Savefig2(SaveFolder(),'MuMaxH85',Dpi())

Plot1('Chalup90','$P.$ $lutheri$',IDataChalup90,MuMaxDataChalup90,'g',IModelChalup90,MuMaxModelChalup90,'g',I,PchlC90,'#00FF00','upper left')
Savefig2(SaveFolder(),'MuMaxC90',Dpi())

Plot1('Sakshaug89','$S.$ $costatum$',IDataSakshaug89,MuMaxDataSakshaug89,'#E50209',IModelSakshaug89,MuMaxModelSakshaug89,'#E50209',I,PchlS89,'pink','lower right')
Savefig2(SaveFolder(),'MuMaxS89',Dpi())

show()
