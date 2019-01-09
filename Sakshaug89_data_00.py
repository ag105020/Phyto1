'''
Created on Mar 6, 2018

@author: Keisuke
'''

from pylab import *

class Data:
    def __init__(self,Lightintensity,plotcolor,Dchl,Chl,Dnc,NC,i):
        self.Lightintensity=Lightintensity
        self.plotcolor=plotcolor
        self.Dchl=Dchl
        self.Chl=Chl
        self.Dnc=Dnc
        self.NC=NC
        self.i=i
    
def data():
    data0=genfromtxt('Sakshaug89.csv', delimiter=',')
    Lightintensity=[12.03,70.56,99.44,1202.78]
    plotcolor=['red','green','blue','cyan']
    I=arange(0,4,1)
    d0=Data(Lightintensity[0], plotcolor[0], data0[:,1], data0[:,2], data0[:,14], data0[:,15],I[0])
    d1=Data(Lightintensity[1], plotcolor[1], data0[:,4], data0[:,5], data0[:,17], data0[:,18],I[1])
    d2=Data(Lightintensity[2], plotcolor[2], data0[:,7], data0[:,8], data0[:,20], data0[:,21],I[2])
    d3=Data(Lightintensity[3], plotcolor[3], data0[:,10], data0[:,11], data0[:,23], data0[:,24],I[3])
    return(d0,d1,d2,d3)