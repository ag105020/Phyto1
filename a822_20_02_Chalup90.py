'''
Created on Sep 29, 2016
This one matches the values of a800_04_28_00
@author: Keisuke
'''

from pylab import *
from Dmax_computation15 import *
from Dmax_approximation01 import *
from PlotSetting01 import *
from Savetxt import *

rc('text', usetex=True)
rc('font', family='serif')
rc('font', serif='Times New Roman')

Lightintensity=arange(0.1,200+0.1,0.1)
What_is_limiting=1

m=4.78224745056655E-19         #(mol C s-1 cell-1) maintenance carbonhydrate consumption (idea from 172-7)
Pmax0=7                     #(g C /(g Chl h) Maximum production rate per chlorophyll (around 6 by Cullen 1990)
Mchl=893.49             #(g / mol chlorophyll) mollar mass of chlorophyll
Pmax=0.00209835928190267      #(mol C s-1 mol chl-1) carbon fixing rate (156-10) (156-15) for unit conversion)

OT=0.0065661503085636
Ynphoto_chl=3.0056463754677          #((molN cell-1)/(molC chl cell-1)) the stoichiometric ratio for cell photosynthetic enzyme (Rubisco etc.) nitrogen to chlorophyll (193-25)
Cnbiosynth=2.53367759328165E-10         #(molN cell-1 s) Constant for varible part of biosynthesis protein nitrogen (193-37)
Nconst_protein=1.18274979779916E-15    #(molN cell-1) Constant protein pool in nitrogen (193-25)
Nstore_max=2.91679384515998E-15        #(molN cell-1) Constant protein pool in nitrogen (193-25)
Cnrna_variable=6212.59249917364       #(s) Constant for Variable part of RNA (193-26)
Ypthylakoid_chl=0.0281633095303638        #((molP cell-1)/(molC chl cell-1)) the shoichiometric ratio for cell phosphorus in thylakoid membrane to chlorophyll (193-26)
Pconst_other=5.44534485638617E-17            #(molP cell-1) Constant part of phosphorus (193-26) * This includes ATP ADP, Phospholipid and DNA RNA
Qp_max=25.26/(3.097e16)                                                                              #(molP cell-1) total phosphorus content in the cell (193-26)
Cessential=3.57361118517991E-14         #(molC cell-1) essential carbon (lipid membrane, etc.) *8.33e-14/10 is 10%

Dmax0=dmax_comp(Lightintensity,What_is_limiting,m,Pmax,OT,Ynphoto_chl,Cnbiosynth,Nconst_protein,Cnrna_variable,\
        Ypthylakoid_chl,Pconst_other,Qp_max,Cessential,Nstore_max)

Dmax_appro0=dmax_appro(Lightintensity,What_is_limiting,m,Pmax,OT,Ynphoto_chl,Cnbiosynth,Nconst_protein,Cnrna_variable,\
        Ypthylakoid_chl,Pconst_other,Qp_max,Cessential,Nstore_max)

Dmax=Dmax0.Dmax
ChltoCplot=Dmax0.ChltoCplot
NtoCplot=Dmax0.NtoCplot
PtoCplot=Dmax0.PtoCplot
BiomassC=Dmax0.BiomassC

Dmax_appro=Dmax_appro0.Dmax

Dmax[Dmax<0]=0
Dmax_appro[Dmax_appro<0]=0

I=Lightintensity
Pchl=Pmax*(1-exp(-OT*I))

A=800
B=0.1
Dmax_PI=Pchl*A-B
Dmax_PI[Dmax_PI<0]=0

fig=figure(1,figsize=(8,6.5))
ax1=fig.add_subplot(111)
#ax2=ax1.twinx()

Data_light=(62.96,188.66)
Dmaxdata=(0.57,1.14)

ax1.plot(Data_light,Dmaxdata,'o',color='cyan',label='Data',zorder=4)
ax1.plot(Lightintensity,Dmax*86400,label='Model',zorder=3)
#ax2.plot(Lightintensity,Pchl,'--',color='#00FF00',label='$P_{Chl}$',zorder=1)
ax1.set_xlabel('Irradiance \\textrm{\\greektext m}mol m$^{-2}$ s$^{-1}$')
ax1.set_ylabel('$\mu_{max}$ (d$^{-1}$)')
#ax2.set_ylabel('$P_{Chl}$ (mol C h$^{-1}$ Chl mol C$^{-1}$)')
#ax2.set_ylim(0,0.003)

ax1.set_ylim(0,1.8)
#ax1.plot(Lightintensity,Pchl/0.003*1.8,'--',color='#00FF00',label='$P_{Chl}$')

ax1.legend(loc=4,fontsize=25,borderaxespad=0.5,labelspacing=0.3)
Whatislimiting='N-limiting'
savefig("C:\\Users\\Keisuke\\desktop\\figures\\Phytoplankton model\\"+Whatislimiting+"\\MuMax\\Chalup90.png")
Savetxt(Dmax*86400,"Phytoplankton model\\"+Whatislimiting+"\\MuMax","Chalup90")
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#Plot free parameters
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
from decimal import Decimal
def sci(Numb):
    Numb='%.2E' % Decimal(Numb)
    return Numb

print("m",sci(m))
print("Pmax",sci(Pmax))
print("Apho",sci(OT))
print("Ynphoto_chl",sci(Ynphoto_chl))
print("Abio",sci(Cnbiosynth))
print("Cother_protein",sci(Nconst_protein))
print("Arna",sci(Cnrna_variable))
print("Ythylakoid_chl_P",sci(Ypthylakoid_chl))
print("Pconst_other",sci(Pconst_other))
print("Nstore_max",sci(Nstore_max))
print("Cessential",sci(Cessential))
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# 
# 
# figure(2)
# plot(Dmax,BiomassC)

show()