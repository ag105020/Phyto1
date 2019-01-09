'''
Created on Mar 6, 2018

@author: Keisuke
'''

from pylab import *

class Chalup90:
    def __init__(self,Lightintensity,plotcolor,Dchl,Chl,Dnc,NC,i):
        self.Lightintensity=Lightintensity
        self.plotcolor=plotcolor
        self.Dchl=Dchl
        self.Chl=Chl
        self.Dnc=Dnc
        self.NC=NC
        self.i=i
    
def chalup90():
    data0=genfromtxt('Chalup90.csv', delimiter=',')
    Lightintensity=[62.96,188.66]
    plotcolor=['red','blue']
    I=arange(0,2,1)
    chalup1=Chalup90(Lightintensity[0], plotcolor[0], data0[:,1], data0[:,2], data0[:,8], data0[:,9],I[0])
    chalup2=Chalup90(Lightintensity[1], plotcolor[1], data0[:,4], data0[:,5], data0[:,11], data0[:,12],I[1])
    
    return(chalup1,chalup2)