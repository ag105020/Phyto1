#import pylab as mpl

from pylab import *

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# Plot
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

rcParams.update({'font.size': 30,
                 'lines.markersize':10,
                 'lines.markeredgewidth':0.5})
rcParams.update({'xtick.major.pad': 15})
rcParams.update({'xtick.major.pad': 15})
rcParams.update({'font.serif': 'Times New Roman'})
rcParams.update({'figure.autolayout': True})
rcParams['figure.figsize']=8,6.5
rcParams.update({'figure.facecolor':'W'})
rcParams.update({'lines.linewidth':2.5})   

#for typing mu (micro) 
rcParams.update({'text.usetex':True})   #to call real latex
rcParams.update({'text.latex.preamble':['\\usepackage[greek,english]{babel}']})

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO