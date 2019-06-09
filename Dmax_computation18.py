'''
Created on Sep 29, 2016

@author: Keisuke
'''

from pylab import *
from af001_energy_calculation import *
from Solver_2D import *

class Dmax_comp:
    def __init__(self,Dmax,ChltoCplot,NtoCplot,NtoPplot,PtoCplot,BiomassC,Prna,Pthylakoid,Pdna,Pstore,Limitation,Cphoto,Cbiosynth,Cconst_protein):
        self.Dmax=Dmax
        self.ChltoCplot=ChltoCplot
        self.NtoCplot=NtoCplot
        self.NtoPplot=NtoPplot
        self.PtoCplot=PtoCplot
        self.BiomassC=BiomassC
        self.Prna=Prna
        self.Pthylakoid=Pthylakoid
        self.Pdna=Pdna
        self.Pstore=Pstore
        self.Limitation=Limitation
        self.Cphoto=Cphoto
        self.Cbiosynth=Cbiosynth
        self.Cconst_protein=Cconst_protein

def dmax_comp(Lightintensity,What_is_limiting,m,Pmax,OT,Ynphoto_chl,Cnbiosynth,Nconst_protein,Cnrna_variable,\
        Ypthylakoid_chl,Pconst_other,Qp_max,Cessential,Nstore_max):
    
    I=Lightintensity
    
    if What_is_limiting==0:
    #For P-limiting case
        Pin=0.002  #(mol/m3) Phosphorus concentration in the incoming medium (Healey 1985)
        Nin=0.2    #(mol/m3) Nitrate concentration in the incoming medium (Healey 1985)
        Qc=1.00*10**(-12)/12      #(molC/cell) biomass C per cell (196-18)(from Healey 1985)
    
    elif What_is_limiting==1:
    #For N-limiting case
        Pin=0.02  #(mol/m3) Phosphorus concentration in the incoming medium (Healey 1985)
        Nin=0.05    #(mol/m3) Nitrate concentration in the incoming medium (Healey 1985)
        Qc=10**(-12)/12      #(molC/cell) biomass C per cell (196-18)(from Healey 1985)
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Model preparation
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Other phytoplankton parameters
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #ooooooooooooooooooooooooooooooooooooooooooooooooooo
    #From the previous work (mainly from a800_04_27)
    #ooooooooooooooooooooooooooooooooooooooooooooooooooo
    
    E3=evalue()
    E=E3.E
    
    Molar_mass_DNA_AT_average=307.47        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    Molar_mass_DNA_CG_average=307.97        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    
    Molar_mass_RNA_AT_average=316.47        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    Molar_mass_RNA_CG_average=323.97        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    
    #================================
    #E coli
    #================================
    CG_Ecoli=0.506          #(dimensionless) from [http://www.ncbi.nlm.nih.gov/genome/167 (accessed 06/18/2016)]
    AT_Ecoli=1-CG_Ecoli     #(dimensionless) 
    
    Molar_mass_DNA_Ecoli=Molar_mass_DNA_AT_average*CG_Ecoli+Molar_mass_DNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of DNA unit
    Molar_mass_RNA_Ecoli=Molar_mass_RNA_AT_average*CG_Ecoli+Molar_mass_RNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of RNA unit
    
    RNA_DNA_mass_ratio=20/7.6   #(ug/ug) Bremer and Dennis 1996
    RNA_DNA_mass_ratio=17.844/6.5239  #(ug/ug) from values ad D=0 "07 Bremer and Dennis 1996 data plot.xlsx"
    
    RNA_DNA_molar_ratio=RNA_DNA_mass_ratio/Molar_mass_RNA_Ecoli*Molar_mass_DNA_Ecoli    #(mol mol-1)
    #================================
    #Stoichiometric parameters
    #================================
    YcyanoC_N=2                             #(molC molN) C/N molar ratio of cyanophycin
    YpgC_P=40                           #(molC molP) C/P molar ratio of PG: Phosphatidyl glycerol (assuming C16 fatty acids (since actually mostly C16 (Huflejt et al., 1990)
    
    CG=0.563                  #GC%    [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 12/14/2018)]
    YnucacidP_N=1/(3.5*(1-CG)+4*CG)               #(molP molN-1) P/N molar ratio of RNA (193-26) values (193-28) excel file "08 N to P ratio in DNA and RNA.xlsx"
    
    YdnaC_N=3.5*(1-CG)+2.5*CG       #(molC molN-1) C/N molar ratio of dna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)
    YrnaC_N=3.25*(1-CG)+2.5*CG      #(molC molN-1) C/N molar ratio of rna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)
    
    DNAmb=2.1269                   #(Mb) Megabase pair of synechococcus DNA in mega (million) base pairs [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 12/14/2018)]
    Avogadro=6.022*10**23           #(molecules mol-1) Avogadro constant
    Pdna_const=DNAmb*2*10**6/Avogadro                #(molP cell-1) Constant part of DNA in phosphorus 
    Prna_const=Pdna_const*RNA_DNA_molar_ratio       #(molP cell-1) Constant part of RNA in phosphorus
    #* Make sure to multiply by 2 as they are base PAIRs"
    Ndna_const=Pdna_const/YnucacidP_N      #(molN cell-1) Constant part of DNA in nitrogen
    Nrna_const=Ndna_const*RNA_DNA_molar_ratio   #(molN cell-1) Constatn part of RNA in nitrogen
    Ndna=Ndna_const    #(molN cell-1) DNA in nitrogen (here assuming constant)
    Pdna=Ndna*YnucacidP_N   #(mol P cell-1) DNA in phosphorus
    
    YphotoFe_N=0.001636364  #(molFe molN-1) Fe/N ratio in photosystem iron
    
    #=======================================
    #Calculation of carbon usage (195-16)
    #=======================================
    CNprotein=4.49   #(molC molN) the ratio of C to N in protein (derived from Brown 1991) calculation in "13 Amino acid composition of different phytoplankton.xlsx"
    YdnaC_N=3.5*(1-CG)+2.5*CG       #(molC molN-1) C/N molar ratio of dna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)
    YrnaC_N=3.25*(1-CG)+2.5*CG      #(molC molN-1) C/N molar ratio of rna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)
    YpgC_P=40                           #(molC molP) C/P molar ratio of PG: Phosphatidyl glycerol (assuming C16 fatty acids (since actually mostly C16 (Huflejt et al., 1990)
    
    #ooooooooooooooooooooooooooooooooooo
    #Photosynthetic parameters
    #ooooooooooooooooooooooooooooooooooo
    Pmax0=7                     #(g C /(g Chl h) Maximum production rate per chlorophyll (around 6 by Cullen 1990)
    Mchl=893.49             #(g / mol chlorophyll) mollar mass of chlorophyll
    Pchl=Pmax*(1-exp(-OT*I))  #(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)(193-25)
    
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    # Preparation Kei 206-4
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    A=((1+E)*Qc*Ynphoto_chl)/Pchl+Cnbiosynth
    B=Nconst_protein+(m*Ynphoto_chl)/Pchl
    L=((1+E)*Qc*Ypthylakoid_chl)/Pchl
    M=(m*Ypthylakoid_chl)/Pchl
    N=YrnaC_N*A*Cnrna_variable
    O=(1+E)*Qc/Pchl+CNprotein*A+YrnaC_N*B*Cnrna_variable+YpgC_P*L
    P=m/Pchl+CNprotein*B+YrnaC_N*Nrna_const+YpgC_P*M+Cessential+Ndna*YdnaC_N-Qc
    
    if size(I)==1:
        Dmax0=solver_2D(float(N),float(O),float(P))
    else:
        Dmax0=solver_2D(N,O,P)
    
    Dmax=Dmax0.X
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Computation of elemental ratios
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ls=Dmax*Qc
    D=Dmax
    Chl=((1+E)*ls+m)/Pchl       #(molC chl cell-1) Chlrophyll concentration (193-25) 
    Nphoto=Chl*Ynphoto_chl  #(molN cell-1) Photosynthesis related protein nitrogen (193-25)
    Nbiosynth=D*Cnbiosynth             #(molN cell-1) various part of biosynthesis related protein in N (193-37)
    Nprotein=Nphoto+Nconst_protein+Nbiosynth    #(molN cell-1) All the proteins in N (193-26)
    Nrna_variable=Nprotein*D*Cnrna_variable        #(molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
    Ndna_variable=0*D                     #(molN cell-1) While Bremer and Dennis shows increasing trend, Parrott 1980 shows decreasing trend.
    Nchl=Chl*4/55                           #(molN cell-1) Chlorophyll nitrogen (actually almost negligiable)
    Pthylakoid=Chl*Ypthylakoid_chl          #(molP cel-1) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)
    Prna_variable=Nrna_variable*YnucacidP_N     #(molP cell-1) variable part of phosphorus in RNA (193-26)
    Pdna_variable=Ndna_variable*YnucacidP_N     #(molP cell-1) variable part of phosphorus in DNA (193-26)

    #=================================
    #Total calculation
    #=================================
    Qn_max=Nprotein+Nrna_variable+Nrna_const+Ndna_variable+Ndna_const+Nchl+Nstore_max      #(molN cell-1)  nitrogen content in the cell (193-26)                                 #(molN cell-1) total phosphorus content in the cell (193-26)
    Qp_max=25.26/(3.097e16)                                                                              #(molP cell-1) total phosphorus content in the cell (193-26)
    Qn_min=Nprotein+Nrna_variable+Nrna_const+Ndna_variable+Ndna_const+Nchl            #(molN cell-1) total nitrogen in the cell without storage
    Qp_min=Pconst_other+Pthylakoid+Prna_variable+Prna_const+Pdna_variable+Pdna_const      #(molP cell-1) total phosphorus in the cell without storage

    #=================================
    #Vector preparation
    #=================================
    Nstore=zeros(size(I))
    X=zeros(size(I))
    Qn_test=zeros(size(I))
    Qp_test=copy(X)
    Qp=copy(X)
    Qn=copy(X)
    Pstore=copy(X)
    Limitation=copy(X)
    #=================================
    #Population calculation
    #=================================
    Xn_max=Nin/Qn_min
    Xp_max=Pin/Qp_min
    U=arange(0,size(I),1)
    for i in U:
        if Xn_max[i]>Xp_max[i]:
            X[i]=Xp_max[i]
            Qp[i]=Qp_min[i]
            Qn_test[i]=Nin/X[i]
            Limitation[i]=0
            if Qn_test[i]<Qn_max[i]:
                Qn[i]=Qn_test[i]
                Nstore[i]=Qn_test[i]-Nprotein[i]-Nrna_variable[i]-Nrna_const-Ndna_variable[i]-Ndna_const-Nchl[i]  #(molN cell-1) Nitrogen storage in the cell
            else:
                Qn[i]=Qn_max[i]
                Nstore[i]=Nstore_max
        else:
            X[i]=Xn_max[i]
            Qn[i]=Qn_min[i]
            Qp_test[i]=Pin/X[i]
            Limitation[i]=1
            if Qp_test[i]<Qp_max:
                Qp[i]=Qp_test[i]
            else:
                Qp[i]=Qp_max
            Pstore[i]=Qp[i]-Pconst_other-Pthylakoid[i]-Prna_variable[i]-Prna_const-Pdna_variable[i]-Pdna_const   #(molP cell-1) Stored phosphorus in the cell
    

    #-------------------
    #C protein
    #-------------------
    Cphoto=Nphoto*CNprotein     #(molC cell-1) carbon in photosystem protein (195-16)
    Cbiosynth=Nbiosynth*CNprotein   #(molC cell-1) carbon in biosynthesis protein (195-16)
    Cconst_protein=Nconst_protein*CNprotein*ones(size(I))  #(molC cell-1) carbon in other protein assumed constant (195-16)

    #----------------------
    #C DNA RNA
    #----------------------
    Crna_const=Nrna_const*YrnaC_N       #(molC cell-1) carbon in variable part of RNA (195-16)
    Crna_variable=Nrna_variable*YrnaC_N     #(molC cell-1) carbon in variable part of RNA (195-16)
    
    Cdna_const=Ndna_const*YdnaC_N       #(molC cell-1) carbon in constant part of DNA (195-16)
    Cdna_variable=Ndna_variable*YdnaC_N     #(molC cell-1) carbon in variable part of DNA (195-16)
    
    Cnstore=Nstore*YcyanoC_N        #(molC cell-1) carbon in nitrogen storage (cyanophycin)
    CthylakoidPG=Pthylakoid*YpgC_P           #(molC cell-1) carbon in PG (phosphatidyl glycerol) in thylakoid membranes
    
    
    #---------------------------------------------------
    #C other: Here revised to include Nstore reduction
    #---------------------------------------------------
    Cchl=Chl
    #print(Qc)
    Cother_without_Nstore=Qc-Cphoto-Cbiosynth-Cconst_protein-Cchl\
            -Crna_const-Crna_variable-Cdna_const-Cdna_variable\
            -Cessential-CthylakoidPG
    
    Cother_with_full_Nstore=Qc-Cphoto-Cbiosynth-Cconst_protein-Cchl\
            -Crna_const-Crna_variable-Cdna_const-Cdna_variable\
            -Cessential-Cnstore-CthylakoidPG
            
#     Cother=Cother_with_full_Nstore            
#     
#     Nstore_reduce=logical_and(Cother_without_Nstore>0, Cother_with_full_Nstore<0)
#     
#     Cnstore[Nstore_reduce]=Cother_without_Nstore[Nstore_reduce]
#     Nstore0=copy(Nstore)
#     Nstore[Nstore_reduce]=Cother_without_Nstore[Nstore_reduce]/YcyanoC_N
#     Cother[Nstore_reduce]=0
    Qn=Qn-Nstore #+Nstore[Nstore_reduce]
    
    
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For plotting 1
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    BiomassC=12*X*Qc   #(mg C L-1) Biomass concentration
    NtoCplot=Qn/Qc*14*10**6/(12*10**3)    #(ug N / mg C) biomass N to C ratio (164-20)
    PtoCplot=Qp/Qc*30.97*10**6/(12*10**3)    #(ug P / mg C) biomass P to C ratio (164-20)
    NtoPplot=Qn/Qp*14*10**6/(30.97*10**6)        #(ug N /ug P) biomass N to P ratio (164-20)
    ChltoC0=Chl/Qc         #(mol C chl mol C -1) Chlorophyll to carbon ratio
    Mchl=893.49             #(g / mol chlorophyll) mollar mass of chlorophyll
    ChltoCplot=ChltoC0/12/1000*Mchl/55*10**6     #(ug chlorophyll a mg C-1) (see 157-36 for conversion)
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For stack plot
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Nrna=Nrna_variable+Nrna_const
    Prna=Nrna*YnucacidP_N
    Pthylakoid=Pthylakoid
    
    
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@





    X1=Dmax_comp(Dmax,ChltoCplot,NtoCplot,NtoPplot,PtoCplot,BiomassC,Prna,Pthylakoid,Pdna,Pstore,Limitation,Cphoto,Cbiosynth,Cconst_protein)


 
    
    return(X1)