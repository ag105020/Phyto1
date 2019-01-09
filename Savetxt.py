'''
Created on May 18, 2014

@author: Keisuke
'''
from pylab import * 

def Savetxt(parameter,savefolder,name):
    First_part="C:\\Users\\Keisuke\\Desktop\\figures\\"
    Second_part=savefolder+"\\"+name
    Last_part=".csv"
    savetxt(First_part+Second_part+Last_part, parameter, delimiter=",",fmt='%.8e')