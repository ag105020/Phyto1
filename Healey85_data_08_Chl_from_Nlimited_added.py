'''
Created on May 18, 2014

@author: Keisuke
'''

from pylab import *


class Healey85C:
    def __init__(self,Lightintensity,plotcolor,DL,L,Dchl,chl,D_PtoC,PtoC,D_NtoC,NtoC, D_NtoP, NtoP,DChl_Nlimited,Chl_Nlimited,i):
        self.Lightintensity=Lightintensity
        self.plotcolor=plotcolor
        self.DL=DL                  #dilution rate
        self.L=L                    #carbon concentration
        self.Dchl=Dchl              #(-d) dilution rate to chlorophyll to Carbon
        self.chl=chl                #(ug chlorophyll / mg C) chlorophyll to Carbon
        self.D_PtoC=D_PtoC          #(-d) dilution rate for the data of organic P to C ratio
        self.PtoC=PtoC              #(ug P / mg C) organic P to C ratio
        self.D_NtoC=D_NtoC          #(-d) dilution rate for the data of orgenic N to C ratio
        self.NtoC=NtoC              #(ug N / mg C) organic N to C ratio
        self.D_NtoP=D_NtoP          #(-d) dilution rate for the data of organic N to P ratio
        self.NtoP=NtoP              #(ug N / ug P) organic N to P ratio
        self.DChl_Nlimited=DChl_Nlimited  #(-d) dilution rate for chlorophyll to Carbon
        self.Chl_Nlimited=Chl_Nlimited  #(-d) dilution rate for chlorophyll to Carbon
        self.i=i
        
def healey85C():
    data0=genfromtxt('Healey_1985_H_Nlimited_Chl_added.csv', delimiter=',')
    nodata=zeros(data0.shape[0])+nan
    Lightintensity=[12,22,38,62,144]
    plotcolor=['red','green','blue', 'cyan','purple' ]
    I=arange(0,5,1)
    
    healey_12=Healey85C(Lightintensity[0],plotcolor[0], data0[:,1], data0[:,2], data0[:,18], data0[:,19], data0[:,35], data0[:,36], data0[:,52], data0[:,53], data0[:,69], data0[:,70], data0[:,86], data0[:,87],I[0])
    healey_22=Healey85C(Lightintensity[1],plotcolor[1], data0[:,4], data0[:,5], data0[:,21], data0[:,22], data0[:,38], data0[:,39], data0[:,55], data0[:,56], data0[:,72], data0[:,73], data0[:,89], data0[:,90],I[1])
    healey_38=Healey85C(Lightintensity[2],plotcolor[2], data0[:,7], data0[:,8], data0[:,24], data0[:,25], data0[:,41], data0[:,42], data0[:,58], data0[:,59], data0[:,75], data0[:,76], data0[:,92], data0[:,93],I[2])
    healey_62=Healey85C(Lightintensity[3],plotcolor[3], data0[:,10], data0[:,11], data0[:,27], data0[:,28], data0[:,44], data0[:,45], data0[:,61], data0[:,62], data0[:,78], data0[:,79], data0[:,95], data0[:,96],I[3])
    healey_144=Healey85C(Lightintensity[4],plotcolor[4], data0[:,13], data0[:,14], data0[:,30], data0[:,31], data0[:,47], data0[:,48], data0[:,64], data0[:,65], data0[:,81], data0[:,82], data0[:,98], data0[:,99],I[4])   
#     
#     healey_chl_12=Healey85C(Lightintensity[0],plotcolor[0], data0[:,18], data0[:,19])
#     healey_chl_22=Healey85C(Lightintensity[1],plotcolor[1], data0[:,21], data0[:,22])
#     healey_chl_38=Healey85C(Lightintensity[2],plotcolor[2], data0[:,24], data0[:,25])
#     healey_chl_62=Healey85C(Lightintensity[3],plotcolor[3], data0[:,27], data0[:,28])
#     healey_chl_144=Healey85C(Lightintensity[4],plotcolor[4], data0[:,30], data0[:,31])   
#     
    
    return (healey_12,healey_22,healey_38,healey_62,healey_144)#,healey_chl_12,healey_chl_22,healey_chl_38,healey_chl_62,healey_chl_144)#DL12, L12, DL22, L22, DL38, L38, DL62, L62, DL144, L144

# a=healey85C()
# print(a)
# import csv
# with open ('TEST.csv') as f:
#     reader=csv.reader(f)
#     for column in reader:
#         print(column)
# def Healey85C():    
# 
#     a=genfromtxt('Healey_1985_C.csv', delimiter=',')
# #    print(a)        
# #savetxt("RR.csv", a, delimiter=",",fmt='%2.8f')
# 
# 
#     DL12=a[:,1]
#     L12=a[:,2]
#     DL22=a[:,4]
#     L22=a[:,5]
#     DL38=a[:,7]
#     L38=a[:,8]
#     DL62=a[:,10]
#     L62=a[:,11]
#     DL144=a[:,13]
#     L144=a[:,14]
#     
    
    




# 
# b=a[:,1]
# c=a[:,2]
# d=a[:,4]
# e=a[:,5]
# 
# figure(1)
# plot(b,c,'ro')
# plot(d,e,'bo')

# show()