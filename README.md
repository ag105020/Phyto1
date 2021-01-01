# Code Policy
Please state “Cell Flux Model” and “Keisuke Inomura” in the acknowledgement when your
publication includes the results based on the original/revised code. Or you may consider
including Keisuke Inomura as a co-author depending on the contribution. In either case, the
publication must cite the following paper:
* Inomura K, Omta AW, Talmy D, Bragg J, Deutsch C, Follows MJ. 2020. A Mechanistic model of macromolecular allocation, elemental stoichiometry, and growth rate in phytoplankton. Frontiers in Microbiology 11:1–22.

(Paper downloaded from https://www.inomura.com/papers)

Keisuke Inomura (University of Rhode Island)
inomura@uri.edu

# Phyto1
Phytoplankton model 1

The code is in python 3. Basic python modules such as NumPy and Matplotlib are necessary.

Preparation;
1. Download all the files and put it in one folder.
2. Arrange path for each file

Ploting figures;
Note that each file has a flag for N limiting and P limiting;
"What_is_limiting=1  #0: P-limiting  1:N-limiting".
Change values (0 or 1) according to the purpose.

Fig. 4; Run a800_05_12_07

Fig. 5, 6; Run a800_05_12_06 (for subplot A),
        Run a800_04_43_05 (for subplot B)

Fig. 7; Run a800_04_33_90

Fig. 8; First, run a800_05_12_06, a822_20_02, and a822_30_02_10 to create output files.
        Then, run a1000_01_05 for plotting them.

Fig. 9; Run a800_05_12_07

Fig. 10; Run a800_04_33_92 (for subplot A),
         Run 1100_00_01 (for subplot B),
         Run a800_05_12_07_01 (for subplot C),
         Run a1101_00_00 (for subplot D).

Fig. S2; Run a800_04_60_11_01 and a800_04_70_11_10_01 for the figures in the top and bottom rows, respectively.

Fig. S3; Run a800_04_43_57

Fig. S4; Run a800_05_12_07 (for subplot i),
         Run a800_04_43_06 (for other subplots)
         
Fig. S5; Run a800_04_33_89 (For subplot A),
         Run a800_05_12_07 (For subplot B)
         
Fig. S6; Run a800_04_60_11_01 (For subplot A),
         Run a800_07_11_10_01 (For subplot B)
         
Fig. S7; Run a800_05_12_07

