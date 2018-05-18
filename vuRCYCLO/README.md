#### CODE NAME:<br>
   vuRCYCLO<br><br>


####  AUTHOR:<br>
  ZHAO Deng<br>
  email: zhao.deng@foxmail.com<br><br>


 #### ADVISER: <br>
   R. E. Waltz<br><br>
   
#### PURPOSE: <br>
  To plot the data from rCYCLO code, and analyze the physics.<br><br>
   
   
####  HOW TO START:<br>
 1. Make sure rCYCLO output files needed for plotting and all the vuRCYCLO .pro files are in the same directory.<br>
 2. Compile cgplot first.<br>
 3. Compile vuRCYCLO<br>
 4. Run vuRCYCLO<br><br>
 
####  CODE STRUCTURE:<br>
 * Main plots:<br>       To plot the basic information of a case, including: total Chi vs time, total_D vs time,
                   and growth rate (and frequency) vs ky of gyrokinetics(GK), cyclokinetics in Fourier harmonic
                   representation (CKinFH), and cyclokinetics in cyclotron harmonic representation (CKinCH) 
                   respectively.<br>
 
 * More details:<br>     To plot more details e.g. the 3D plot of Chi, the distribution function in velocity space.<br>
 * Colorful plot:<br>    To plot the colorful contour figures of growth rate and frequency, and the eddy in (x,y) space
                   evolution in time. Pay attention to that the 'Main plot' and 'More details' applies the default.  
                   system color, however, the color table used for colorful plot is given by the file 'color.dat'.
                   So there may be a color display problem when operation switch between them. To solve this problem,
                   you need to restart idl.<br>
                   
                    
 * Movie:<br>            Disable, the eddy movie is plotted in 'Colorful plot'<br>
 * ColorMap: <br>        To plot the color table for the Colorful plot<br><br>
 

####  INPUT FILES: <br>
*  inputRCYCLO.txt   -- in order to read the grid parameters and physics variables.<br>
*  inputvuRCYCLO.txt -- in order to read the control parameters for IDL plotting.<br><br>
 
 
#### CONTROL VARIABLES INTRODUCE:<br>
 The variables in inputvuRCYCLO.txt is used for control the idl plot.<br><br>

*****************************************************************************************************************
####  OTHER SPECIFIED FIGURES:<br>
 
 **How to plot the frequency spectrum of Chi (|Phi|^2 is the same) new method (on 2014.7.11)** <br>

 1. Apply the after '2014.7.10.f90' code, and set 'verify_NL=2'<br>

 (a). You can get the 'Phi_kOMG.txt' and 'GOMG.txt' from nt=0 to nt=ntmax, by set Stopnt=0 in 'inputrCYCLO.txt'.<br>

 (b). Or you can restart from any old running restart point. Then 'Phi_kOMG.txt' and 'GOMG.txt' will be from nt=Stopnt <br>
 to nt=ntmax. Remember to set Stopnt in 'inputrCYCLO.txt' be the restart point.<br>
 
 PS: If you use a different output_step value as the original one, then remeber to set fortran code text some 'output_step' <br>
 equal to the smaller one to prevent output/input error. <br>
 
 2. Use the 'Chi_vs_nnew.pro' code control=2, Begintime='the begin ploting point, 'jump' is jumping some time points in <br>
 order to reduce the calculation. (better be 1).Then IDL will output the freqency result in file 'y.txt' (output 'y.txt'<br>
 will save time for repeating operations).  <br>
 3. Set control=2.1 or 2.2 to continue to plot Chi vs omega or Integral_Chi vs omega.<br>

 4. As an code verification, remember to recover Parseval's Theorem.<br><br>


 **A fast way to get 'y.txt' file on cluster:** <br>
 Since the file for freq spectrum plot is very large, the transfor of the file to a PC is very terrible.<br>
 Thus, we offer a way to calculate the short 'y.txt' file on cluster. Following by:<br>
 1.Put a suit of idl code in to the cluster, most importantly, including 'idl2014.7.12_rCYCLO.pro'<br>
 
 2.When Stopnt=0, set the 'Begintime' in the 'inputvuRCYCLO.txt' file to be the  begin plotting point, and set ntmax be the end<br>
 plotting point. When Stopnt>0, ploting will start from Stopnt to ntmax.<br>
 
 3.run 'idl2014.7.12_rCYCLO.pro', and click the 'Energy(DW, ZF, GAM)', then will generate 'y.txt' file.<br>


*****************************************************************************************************************
