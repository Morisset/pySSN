This directory is where the fortran program is. It is easily compiled with:

gfortran XSSN_Phyat.f -o XSSN_Phyat.exe

For the code to run, you will need here another directory "data" containing additional files in separate directories:

* in d1: the 162 Hummer and Storey files, e.g. r2b0005.d. 
* in d2: the 23 Porter files e.g. 23000_1_0_50_Results.txt
* in d3: the 144 Storey file e.g. lines_3.0_3.0_dr 
* in d4: the file HeI_Porter13.dat

