'''
Created on Jan 23, 2016

Note:

About a803s
Simplify the model 
*putting DNA and RNA together
*putting togehter the constant pool
*removing chl N and so on.

Revision note: 193-37

note:
This is a good one. This can be possible for thesis.
This is very similar to a800 05 12 06 except we used same Qc for both N and P limited
@author: Keisuke
'''

from pylab import *
from af001_energy_calculation import *
from Savetxt import *



rc('text', usetex=True)
rc('font', family='serif')
rc('font', serif='Times New Roman')

What_is_limiting=0    #0: P-limiting  1:N-limiting

DPI = 600

#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#Function beging here
#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
def kkI(healey85data,What_is_limiting,m,Pmax,OT,Ynphoto_chl,Cnbiosynth,Nconst_protein,Nstore_max,Cnrna_variable,Ypthylakoid_chl,Pconst_other,Qp_max,Cessential):   #this function calculate for the same irradiance

    I=healey85data.Lightintensity
    plotcolor=healey85data.plotcolor

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Parameters
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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

    
    E3=evalue()
    E=E3.E
    
    R=1.63*10**(-6)
    h=0.0007174       #used in the next equation
    
    ep3=h      #diffusivity of glycolipid layer compared to alginate layer (132-23)
    Ra=0         #(dimensionless) the ratio of alginate layer to the cell radius
    La=R*Ra         #(m) thickness of alginate layer (132-23)
    
    Rg=0.1           #ratio of glycolipid layer to the cell radius
    Lg=R*Rg             #(m) thickness of glycolipid layer
    x0=1/h*(1/R-1/(R+Lg))+ep3/h*1/(R+Lg)+(1/(R+Lg+La))*(1-ep3/h)  #effect of
    r5=1/(x0*R)         #diffusivity efficiency of cell membrane
    V=4/3*pi*R**3   #Volume of the cell (m3/cell)
    C=6
    
    Ddmax=1.6
    #Dmax=4
    Dstep=0.001
    Dd=arange(Dstep,Ddmax+Dstep,Dstep)       #(h-1) growth rate
    U=arange(0,Ddmax/Dstep,1)
    D=Dd/(3600*24)
    #Qc=10**(-12)/12      #(molC/cell) biomass C per cell (196-18)(from Healey 1985)
    rho=Qc/V             #(molC m-3) biomass C per volume in the cell
    Do20=2.12*10**(-9)  #Diffusion coefficient of O2 in the water (m2/s)
    Do2=Do20*r5         #Diffusion coefficient of O2 in the water (m2/s)
    Dch0=6.728*10**(-10)  #Diffusion coefficient of glucose in the water (m2/s)
    Dch=Dch0*r5          #Diffusion coefficient of glucose in the water (m2/s)
    Dnh40=1.98*10**(-9)  #Diffusion coefficient of ammonium in the water (m2/s)
    Dnh4=Dnh40*r5        #Diffusion coefficient of ammonium in the water (m2/s)
    Dno30=1.9*10**(-9)  #Diffusion coefficient of nitrate in the water (m2/s) (Li 1974 at 25C)
    Dno3=Dnh40*r5        #Diffusion coefficient of nitrate in the water (m2/s)
    Dh2po40=8.46*10**(-10) #Diffusion coefficient of H2PO4- in water at 25C (m2 s-1) (Li 1974)
    Dhpo40=7.34*10**(-10) #Diffusion coefficient of HPO42- in water at 25C (m2 s-1) (Li 1974)
    Dp0=(Dh2po40+Dhpo40)/2   #Take the average as at around Ph 7.5 (healey 1985 paper' initial condition) They exist about half and half
    Dp=Dp0*r5           #diffusion coefficient of P through membrane
    Mchl=893.49             #(g / mol chlorophyll) mollar mass of chlorophyll
    
    #==============================
    #New parameters
    #==============================
    
    #------------------------------
    #Photosynthesis
    #------------------------------

    Pchl=Pmax*(1-exp(-OT*I)) #(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)(193-25)
    
    ls=D*Qc                    #(molC s-1) Biomass synthesis rate (193-25)
    #------------------------------
    
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Key parameters: parameterization ideas -> Kei 193-28
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    #================================
    #Constant parameters
    #================================

    #================================
    
    #================================
    #Y intercept parameters
    #================================

    #================================
    
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
    
    CG=0.563                   #GC% not CG but I started with CG so I stick with it; it does not matter as "AT GC".   [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]
    YnucacidP_N=1/(3.5*(1-CG)+4*CG)               #(molP molN-1) P/N molar ratio of RNA (193-26) values (193-28) excel file "08 N to P ratio in DNA and RNA.xlsx"
    
    YdnaC_N=3.5*(1-CG)+2.5*CG       #(molC molN-1) C/N molar ratio of dna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)
    YrnaC_N=3.25*(1-CG)+2.5*CG      #(molC molN-1) C/N molar ratio of rna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)
    
    print('N:P ',1/YnucacidP_N)
    print('C:P (DNA)',YdnaC_N/YnucacidP_N)
    print('C:P (RNA)',YrnaC_N/YnucacidP_N)
    
    DNAmb=2.1269                   #(Mb) Megabase pair of synechococcus DNA in mega (million) base pairs [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]
    Avogadro=6.022*10**23           #(molecules mol-1) Avogadro constant
    Pdna_const=DNAmb*2*10**6/Avogadro                #(molP cell-1) Constant part of DNA in phosphorus 
    Prna_const=Pdna_const*RNA_DNA_molar_ratio       #(molP cell-1) Constant part of RNA in phosphorus
    #* Make sure to multiply by 2 as they are base PAIRs"
    Ndna_const=Pdna_const/YnucacidP_N      #(molN cell-1) Constant part of DNA in nitrogen
    Nrna_const=Ndna_const*RNA_DNA_molar_ratio   #(molN cell-1) Constatn part of RNA in phosphorus
    
    
    #print(Ndna_const*Avogadro,Nconst_protein*Avogadro)
    #YdnaN_P=1/3.75               #(molP molN-1) P/N molar ratio of DNA (193-26) values (193-28) excel file "08 N to P ratio in DNA and RNA.xlsx"
    
    
    Ynrnaconst_dnarnaconst=1/2            #(dimensionless) N molar ratio of RNA (constant part) to DNA + RNA (constant part) (193-33) reffering to around p.112 of Biology of Prokyariotes
    Yndnaconst_dnarnaconst=1-Ynrnaconst_dnarnaconst     #(dimensionless) N molar ratio of DNA (constant part) to DNA + RNA (constant part)  (193-33) refering to around p.112 of Biology of Prokyariotes
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Calculation
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    Chl=((1+E)*ls+m)/Pchl       #(molC chl cell-1) Chlrophyll concentration (193-25) 
    Nphoto=Chl*Ynphoto_chl  #(molN cell-1) Photosynthesis related protein nitrogen (193-25)
    Nbiosynth=D*Cnbiosynth             #(molN cell-1) various part of biosynthesis related protein in N (193-37)
    Nprotein=Nphoto+Nconst_protein+Nbiosynth    #(molN cell-1) All the proteins in N (193-26)
    Nrna_variable=Nprotein*D*Cnrna_variable        #(molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
    Ndna_variable=Ndna_const*Dd/1.2*(18.3-7.6)/7.6        #(molN cell-1) variable part of nitrogen in DNA (193-26) Increasing ratio based on Bremmer 1996
    Ndna_variable=0*Dd                     #(molN cell-1) While Bremer and Dennis shows increasing trend, Parrott 1980 shows decreasing trend. 
   
    #print(Dd/1.2*(18.3-7.6)/7.6)
    #Ndnarna=Nconst_dnarna+Nnucacid_variable
    Nchl=Chl*4/55                           #(molN cell-1) Chlorophyll nitrogen (actually almost negligiable)
    Pthylakoid=Chl*Ypthylakoid_chl          #(molP cel-1) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)
    Prna_variable=Nrna_variable*YnucacidP_N     #(molP cell-1) variable part of phosphorus in RNA (193-26)
    Pdna_variable=Ndna_variable*YnucacidP_N     #(molP cell-1) variable part of phosphorus in DNA (193-26)
    
  
    #=================================
    #Constant part
    #=================================
 #   Nconst_dna=Nconst_dnarna*Yndnaconst_dnarnaconst     #(molN cell-1) Constant part of DNA in nitorgen (193-33)
 #   Nconst_rna=Nconst_dnarna*Ynrnaconst_dnarnaconst     #(molN cell-1) Constant part of RNA in nitrogen (193-33)
 #   Pconst_dna=Nconst_dna*YdnaN_P   #(molP cell-1) phosphorus in constant part of DNA (193-33)
 #   Pconst_rna=Nconst_rna*YrnaN_P   #(molP cell-1) phosphorus in constant part of RNA (193-33)
    #=================================
    

    
    
    #=================================
    #Total calculation
    #=================================
    Qn_max=Nprotein+Nrna_variable+Nrna_const+Ndna_variable+Ndna_const+Nchl+Nstore_max      #(molN cell-1)  nitrogen content in the cell (193-26)                                 #(molN cell-1) total phosphorus content in the cell (193-26)                                                                             #(molP cell-1) total phosphorus content in the cell (193-26)
    Qn_min=Nprotein+Nrna_variable+Nrna_const+Ndna_variable+Ndna_const+Nchl            #(molN cell-1) total nitrogen in the cell without storage
    Qp_min=Pconst_other+Pthylakoid+Prna_variable+Prna_const+Pdna_variable+Pdna_const      #(molP cell-1) total phosphorus in the cell without storage
    
    
    #=================================
    #Vector preparation
    #=================================
    Nstore=zeros(size(Dd))
    X=zeros(size(Dd))
    Qn_test=zeros(size(Dd))
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
    
    #X=Pin/Qp        #(cells m-3) Number density of cells (193-27)
   # x=Nin/Qn        #(cells m-3) number density of cells (193-27)
    
    #
    
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
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For plotting 2  (calculation of dna, rna, etc)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

   # Ndna=Nconst_dna+Ndna_variable       #(molN cell-1) nitorgen in DNA (193-33)
   # Nrna=Nconst_rna+Nrna_variable       #(molN cell-1) nitrogen in RNA (193-33)
    
   # Pdna=Ndna*YdnaN_P           #(molP cell-1) phosphorus in DNA (193-33)
   # Prna=Nrna*YrnaN_P           #(molP cell-1) phosphorus in RNA (193-33)
    
    #Pconst_other=Pconst-Pconst_dna-Pconst_rna      #(molP cell-1) phosphorus in other parts (ex. phospholipid in outer membrane, ATP, ADP, etc. (assuming constant) (193-33)

    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For plotting 3 (unit adjustment)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Nunit=1/Qc#*14*10**6/(12*10**3)          #((ug N / mgC)/(molN cell-1) unit conversion term (164-20)
    Punit=1/Qc#*30.97*10**6/(12*10**3)       #((ug P / mgC)/(molP cell-1) unit conversion term (164-20)
    Numbertoarray=ones(size(Dd))            #(dimensionless) Number to array converter
    
    NunitData = 14*10**6/(12*10**3)
    PunitData = 30.97*10**6/(12*10**3)
    
    
    #=======================================
    #Calculation of carbon usage (195-16)
    #=======================================
    CNprotein=4.49   #(molC molN) the ratio of C to N in protein (derived from Brown 1991) calculation in "13 Amino acid composition of different phytoplankton.xlsx"
    
    #-------------------
    #C protein
    #-------------------
    Cphoto=Nphoto*CNprotein     #(molC cell-1) carbon in photosystem protein (195-16)
    Cbiosynth=Nbiosynth*CNprotein   #(molC cell-1) carbon in biosynthesis protein (195-16)
    Cconst_protein=Nconst_protein*CNprotein  #(molC cell-1) carbon in other protein assumed constant (195-16)
 
    #----------------------
    #C chlorophyll
    #----------------------
    Cchl=Chl                    #(molC cell-1) carbon in chlorophyll (195-16)
    
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

    #print(Qc)
    Cother_without_Nstore=Qc-Cphoto-Cbiosynth-Cconst_protein-Cchl\
            -Crna_const-Crna_variable-Cdna_const-Cdna_variable\
            -Cessential-CthylakoidPG
    
    Cother_with_full_Nstore=Qc-Cphoto-Cbiosynth-Cconst_protein-Cchl\
            -Crna_const-Crna_variable-Cdna_const-Cdna_variable\
            -Cessential-Cnstore-CthylakoidPG
            
    Cother=Cother_with_full_Nstore            
    
    Nstore_reduce=logical_and(Cother_without_Nstore>0, Cother_with_full_Nstore<0)
    
    Cnstore[Nstore_reduce]=Cother_without_Nstore[Nstore_reduce]
    Nstore0=copy(Nstore)
    Nstore[Nstore_reduce]=Cother_without_Nstore[Nstore_reduce]/YcyanoC_N
    Cother[Nstore_reduce]=0
    Qn[Nstore_reduce]=Qn[Nstore_reduce]+Nstore[Nstore_reduce]-Nstore0[Nstore_reduce]
    
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For plotting 1
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    BiomassC=12*X*Qc   #(mg C L-1) Biomass concentration
    NtoCplot=Qn/Qc#*14*10**6/(12*10**3)    #(mol N / mol C) biomass N to C ratio (164-20)
    
    PtoCplot=Qp/Qc#*30.97*10**6/(12*10**3)    #(mol P / mol C) biomass P to C ratio (164-20)
    NtoPplot=Qn/Qp#*14*10**6/(30.97*10**6)        #(mol N /mol P) biomass N to P ratio (164-20)
    ChltoC0=Chl/Qc         #(mol C chl mol C -1) Chlorophyll to carbon ratio
    Mchl=893.49             #(g / mol chlorophyll) mollar mass of chlorophyll
    ChltoCplot=ChltoC0/12/1000*Mchl/55*10**6     #(ug chlorophyll a mg C-1) (see 157-36 for conversion)
    

    Dd[Cother<0]=nan
    
    #OOOOOOOOOOOOOOOOOOOOOOOOO
    #Print parameters
    #OOOOOOOOOOOOOOOOOOOOOOOOO
    #Updated values
    #RNA_DNA_mass_ratio
    #Ndna_variable=0       
    #CG=0.5755
    #OOOOOOOOOOOOOOOOOOOOOOOOO
#     print("E",E)
#     print("Qc",Qc) 
#     print("YnucacidP_N",YnucacidP_N)
#     print("YpgC_P",YpgC_P)
#     print("YcyanoC_N",YcyanoC_N)
#     print("YdnaC_N",YdnaC_N)
#     print("YrnaC_N",YrnaC_N)
#    print("Qpmax",25/(3.097e16))
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#Plot free parameters #Here units adjusted for paper
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    from decimal import Decimal
    def sci(Numb):
        Numb='%.2E' % Decimal(Numb)
        return Numb
    
    print("m",sci(m/Qc*86400))
    print("Pmax",sci(Pmax*86400))
    print("Apho",sci(OT))
    print("Ynphoto_chl",sci(Ynphoto_chl*CNprotein))
    print("Abio",sci(Cnbiosynth/Qc*CNprotein/86400))
    print("Cother_protein",sci(Nconst_protein/Qc*CNprotein))
    print("Arna",sci(Cnrna_variable/CNprotein/86400*YnucacidP_N))
    print("Ythylakoid_chl_P",sci(Ypthylakoid_chl))
    print("Pconst_other",sci(Pconst_other/Qc))
    print("Nstore_max",sci(Nstore_max/Qc))
    print("Cessential",sci(Cessential/Qc))
    
    #Other parameters
    print("E",sci(E))
    print("Cdna",sci(Cdna_const/Qc))
    print("Prnamin",sci(Prna_const/Qc))
    print("QpMax",sci(Qp_max/Qc))
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    
    #==================================
    #For thesis revision
    #==================================
#     print("Nrna_const",Nrna_const)
#     print("Ndna_const",Ndna_const)
#     print("YnucacidP_N",YnucacidP_N)
#     print("YnucacidP_N",YnucacidP_N)
#     print("YrnaC_N",YrnaC_N)
#     print("YdnaC_N",YdnaC_N)
    
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #5.Plot
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
    
    lowlim=0
    highlim=Ddmax+1e-5
    step=0.2
     
    #==================================
    #Plot control * 1=on other=off
    #==================================
    Plot_C=1
    Plot_NC=1
    Plot_PC=1
    Plot_NP=1
    Plot_Chl=1
    Plot_NC_stack=1
    Plot_PC_stack=1
    Plot_C_stack=1
    Plot_C_stack_small=1

    xmax0=1.60000001
    if healey85data.i==0:
        plotcolor='#902682'
    if healey85data.i==1:
        plotcolor='#3953A9'
    if healey85data.i==2:
        plotcolor='#05A457'
    if healey85data.i==3:
        plotcolor='#FDC014'
    if healey85data.i==4:
        plotcolor='#EE2C3B'
    #print(BiomassC)
    if Plot_C==1:
        figure(3)
        plot(Dd,BiomassC, color=plotcolor,label=str(I)+'\\textrm{\\greektext m}E m$^{-2}$s$^{-1}$')
        if healey85data.i == 3 and What_is_limiting == 0:
            plot(healey85data.DL, healey85data.L,'o',color=plotcolor,zorder=0)
        else:
            plot(healey85data.DL, healey85data.L,'o',color=plotcolor)
          
        xlabel('$\mu$ (d$^{-1}$)')
        ylabel('mg C L$^{-1}$')
        #title('Carbon concentration')
        ylim(ymin=0)
        if What_is_limiting==0:
            ylim(ymax=25)
        elif What_is_limiting==1:
            ylim(ymax=10)
        xlim(xmax=xmax0)
    #    legend(loc='upper right', borderaxespad=0.6, fontsize=20)#, fancybox=True)
     
    if Plot_PC==1:    
        figure(5)
        ord=100+healey85data.i
        plot(healey85data.D_PtoC, healey85data.PtoC/PunitData,'o' ,color=plotcolor,zorder=10000)

        line=plot(Dd,PtoCplot, color=plotcolor,label=str(I)+'\\textrm{\\greektext m}E m$^{-2}$s$^{-1}$',zorder=ord)
        
        if What_is_limiting==1:#ax.set_linewidth(10)
            if healey85data.i==0:
                line[0].set_linewidth(20)
            if healey85data.i==1:
                line[0].set_linewidth(16)
            if healey85data.i==2:
                line[0].set_linewidth(12)
            if healey85data.i==3:
                line[0].set_linewidth(8)
            if healey85data.i==4:
                line[0].set_linewidth(4)
            
        
        xlabel('$\mu$ (d$^{-1}$)')
        ylabel('mol P mol C$^{-1}$')
        #title('P-C ratio')
    #     ylim(ymin=0)
    #    legend(loc='upper left', borderaxespad=0.6, fontsize=18, ncol=2)
        ylim(ymin=0)
        if What_is_limiting==0:
            ylim(ymax=0.010)  
        elif What_is_limiting==1:
            ylim(ymax=0.014)
            
        xlim(xmax=xmax0)
 
    if Plot_NC==1:
        figure(6)
        plot(Dd, NtoCplot, color=plotcolor,label=str(I)+'\\textrm{\\greektext m}E m$^{-2}$s$^{-1}$')
        if healey85data.i == 3 and What_is_limiting == 0:
            plot(healey85data.D_NtoC, healey85data.NtoC/NunitData,'o' ,color=plotcolor,zorder = 0)
        else:
            plot(healey85data.D_NtoC, healey85data.NtoC/NunitData,'o' ,color=plotcolor)
        
        xlabel('$\mu$ (d$^{-1}$)')
        ylabel('mol N mol C$^{-1}$')
        ylim(ymin=0)
        xlim(xmax=xmax0)
        if What_is_limiting==0:
            ylim(ymax=0.3000000001)  
        elif What_is_limiting==1:
            ylim(ymax=0.3000000001)

#########################
# Here New NP plot
#########################
    if Plot_NP==1:
        figure(10)
        #------Determining Zorder factor depending on what is limiting-------
        if What_is_limiting==0:
            ZorderFactor = 1
        elif What_is_limiting==1:
            ZorderFactor = -1
        #--------------------------------------------------------------------
        plot(Dd, NtoPplot, color=plotcolor,label=str(I)+'\\textrm{\\greektext m}E m$^{-2}$s$^{-1}$',zorder=1000-I*ZorderFactor)
        plot(healey85data.D_NtoP, healey85data.NtoP/NunitData*PunitData,'o' ,color=plotcolor,zorder=2000-I*ZorderFactor)
        xlabel('$\mu$ (d$^{-1}$)')
        ylabel('mol N mol P$^{-1}$')
        ylim(ymin=0)
        xlim(xmax=xmax0)
        if What_is_limiting==0:
            ylim(ymax=140)  
        elif What_is_limiting==1:
            ylim(ymax=140)
    #    
###########################

    if Plot_Chl==1:
        figure(7)
        plot(Dd,ChltoCplot,color=plotcolor,label=str(I)+'\\textrm{\\greektext m}E m$^{-2}$s$^{-1}$')
        #  plot(Dh,Dh,color=plotcolor)
        plot(healey85data.Dchl, healey85data.chl,'o' ,color=plotcolor)
        plot(healey85data.DChl_Nlimited, healey85data.Chl_Nlimited,'D' ,color=plotcolor)
           
        xlabel('$\mu$ (d$^{-1}$)')
        ylabel('\\textrm{\\greektext m}g chlorophyll a mg C$^{-1}$')
         # title('Chlorophyll a concentration')
        xlim(xmax=xmax0)
        ylim([0,55])
        yticks(arange(0,56,10))
         
    return

#AAAAAAAAAAAAAAAAAAAAAAAAAAAAA
# Main part
#AAAAAAAAAAAAAAAAAAAAAAAAAAAAA


if What_is_limiting==0:
    #for P limiting data
    from Healey85_data_08_Chl_from_Nlimited_added import healey85C
elif What_is_limiting==1:
    #for N limiting data
    from Healey85_data_09_Nlimited_case import healey85C


#=============================
# Parameter sets
#=============================

m=3.79146798299876E-19         #(mol C s-1 cell-1) maintenance carbonhydrate consumption (idea from 172-7)
Pmax=0.00320513285659728

OT=0.00863364097132997
Ynphoto_chl=3.56099164557551          #((molN cell-1)/(molC chl cell-1)) the stoichiometric ratio for cell photosynthetic enzyme (Rubisco etc.) nitrogen to chlorophyll (193-25)
Cnbiosynth=4.34728279914354E-10        #(molN cell-1 s) Constant for varible part of biosynthesis protein nitrogen (193-37)
Nconst_protein=4.45336898828389E-15     #(molN cell-1) Constant protein pool in nitrogen (193-25)
Nstore_max=2.91679384515998E-15         #(molN cell-1) Constant protein pool in nitrogen (193-25)
Cnrna_variable=6212.59249917364        #(s) Constant for Variable part of RNA (193-26)
Ypthylakoid_chl=0.0281633095303638        #((molP cell-1)/(molC chl cell-1)) the shoichiometric ratio for cell phosphorus in thylakoid membrane to chlorophyll (193-26)
Pconst_other=5.44534485638617E-17               #(molP cell-1) Constant part of phosphorus (193-26) * This includes ATP ADP, Phospholipid and DNA RNA
Qp_max=25.26/(3.097e16)                                                                              #(molP cell-1) total phosphorus content in the cell (193-26)
Cessential=1.51786753491048E-15          #(molC cell-1) essential carbon (lipid membrane, etc.) *8.33e-14/10 is 10%

#==============================

healey85data=healey85C()
for a in healey85data:
    kkI(a,What_is_limiting,m,Pmax,OT,Ynphoto_chl,Cnbiosynth,Nconst_protein,Nstore_max,Cnrna_variable,Ypthylakoid_chl,Pconst_other,Qp_max,Cessential)


#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

from Dmax_computation17 import *
from PlotSetting01 import *

rc('text', usetex=True)
rc('font', family='serif')
rc('font', serif='Times New Roman')

Lightintensity=arange(0.1,500+0.1,0.1)

Dmax0=dmax_comp(Lightintensity,What_is_limiting,m,Pmax,OT,Ynphoto_chl,Cnbiosynth,Nconst_protein,Cnrna_variable,\
        Ypthylakoid_chl,Pconst_other,Qp_max,Cessential,Nstore_max)

Dmax=Dmax0.Dmax
ChltoCplot=Dmax0.ChltoCplot
NtoCplot=Dmax0.NtoCplot
NtoPplot=Dmax0.NtoPplot
PtoCplot=Dmax0.PtoCplot
BiomassC=Dmax0.BiomassC
Prna=Dmax0.Prna
Pthylakoid=Dmax0.Pthylakoid
Pdna=Dmax0.Pdna
Pstore=Dmax0.Pstore
Limitation=Dmax0.Limitation

ChltoCplot[Dmax<0]=nan
NtoCplot[Dmax<0]=nan
NtoPplot[Dmax<0]=nan
PtoCplot[Dmax<0]=nan
BiomassC[Dmax<0]=nan
Prna[Dmax<0]=nan
Pthylakoid[Dmax<0]=nan
Pstore[Dmax<0]=nan
Limitation[Dmax<0]=nan

Dmax[Dmax<0]=0

plottxt=vstack((Dmax,PtoCplot))
savetxt("CPplot.csv", transpose(plottxt), delimiter=",",fmt='%8e')

figure(0)
plot(Lightintensity,Dmax*86400,label='Model')
xlim(0,200.000001)
if What_is_limiting==1:
    Savetxt(Dmax*86400,'Phytoplankton model\\N-limiting\\MuMax','Healey85')


if What_is_limiting==0:
    Dmaxdata=genfromtxt('Healey_data_Dmax_for_AW_method_P_limited.csv',delimiter=',')
if What_is_limiting==1:
    Dmaxdata=genfromtxt('Healey_data_Dmax_for_AW_method_N_limited.csv',delimiter=',')
    
#============================================
#What is limiting: output folder controll
#============================================
if What_is_limiting==0:
    Whatislimiting="P-limiting"
elif What_is_limiting==1:
    Whatislimiting="N-limiting"
#============================================    

NunitData = 14*10**6/(12*10**3)
PunitData = 30.97*10**6/(12*10**3)

Data_light=array([12,22,38,62,144])
plot(Data_light,Dmaxdata,'o',color='cyan',label='Data')
xlabel('Irradiance \\textrm{\\greektext m}mol m$^{-2}$ s$^{-1}$')
ylabel('$\mu_{max}$ (d$^{-1}$)')
legend(loc=4,fontsize=30,borderaxespad=1)
savefig("C:\\Users\\Keisuke\\desktop\\figures\\Phytoplankton model\\"+Whatislimiting+"\\Healey85\\Mumax.png",dpi=DPI)

figure(3)
plot(Dmax*86400,BiomassC,':',color='black')
savefig("C:\\Users\\Keisuke\\desktop\\figures\\Phytoplankton model\\"+Whatislimiting+"\\Healey85\\C concentration.png",dpi=DPI)

figure(5)
plot(Dmax*86400,PtoCplot/PunitData,':',color='black',zorder=1000)
savefig("C:\\Users\\Keisuke\\desktop\\figures\\Phytoplankton model\\"+Whatislimiting+"\\Healey85\\P to C.png",dpi=DPI)

figure(6)
plot(Dmax*86400,NtoCplot/NunitData,':',color='black',zorder=-1)
savefig("C:\\Users\\Keisuke\\desktop\\figures\\Phytoplankton model\\"+Whatislimiting+"\\Healey85\\N to C.png",dpi=DPI)

figure(7)
plot(Dmax*86400,ChltoCplot,':',color='black',zorder=-1)
savefig("C:\\Users\\Keisuke\\desktop\\figures\\Phytoplankton model\\"+Whatislimiting+"\\Healey85\\Chl to C.png",dpi=DPI)

figure(8)
stackplot(Dmax*86400,Pconst_other*ones(size(Dmax)),Pdna*ones(size(Dmax)),Pthylakoid,Prna,Pstore)

figure(9)
plot(Dmax*86400,Limitation)

figure(10)
plot(Dmax*86400,NtoPplot/NunitData*PunitData,':',color='black')
savefig("C:\\Users\\Keisuke\\desktop\\figures\\Phytoplankton model\\"+Whatislimiting+"\\Healey85\\N to P.png",dpi=DPI)

show()
    
    