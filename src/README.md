 #### CODE NAME:<br>
   rCYCLO<br><br>


 #### AUTHOR:<br>
  ZHAO Deng<br>
  email: zhao.deng@foxmail.com<br><br>


  #### ADVISER:<br>
   R. E. Waltz<br><br>

  #### SCRIPT:<br>
   * rCYCLO.f90 ---- the running script.<br>
   * inputRCYCLO.txt ---- the input file.<br><br>

  #### PURPOSE:<br>
rCYCLO code is developed in order to explore the explore the missing transport near L-mode edge problems[1].<br><br>


  #### HOW TO RUN THE CODE:<br>
 Before started, you should make sure there is an available lapack library.<br>
 1. Compile rCYCLO code, remember to connected the lapack lib. e.g.<br>
    For NERSC: ftn rCYCLO.f90<br>
    If install the lapack by yourself: mpif90 2014.5.14_rCYCLO.o -L../lib -llapack -ltmglib -lblas<br>
    lib is the directory where you put the '*.a' compiled lapack library source files.<br>
 2. Submit the job, the number of processors should follow the rules in the annotation of GK_FK_CK variable below.<br><br>


  #### INTRODUCE:<br>
 	Gyrokinetic simulations of L-mode near edge tokamak plasmas with the GYRO code underpredict both
 the transport and the turbulence levels by 5 to 10 fold[1], which suggest either some important
 mechanism is missing from current gyrokinetic codes like GYRO or the gyrokinetic approximation itself
 is breaking down. It is known that GYRO drift-kinetic simulations with gyro-averaging suppressed
 recover most of the missing transport[2]. With these motivations, we developed a flux tube nonlinear
 cyclokinetic[3] code rCYCLO with the parallel motion and variation suppressed. rCYCLO dynamically follows
 the high frequency ion gyro-phase motion (with no averaging) which is nonlinearly coupled into the low frequency
 drift-waves thereby interrupting and possibly suppressing the gyro-averaging. By comparison with the corresponding
 gyrokinetic simulations, we can test the conditions for the breakdown of gyrokinetics. rCYCLO nonlinearly
 couples grad-B driven ion temperature gradient (ITG) modes and collisional fluid electron drift modes to
 ion cyclotron (IC) modes.<br><br>

   rCYCLO code includes four independent parts, controlled by the parameter GK_FK_CK=0, 1, 2, and -2.<br>
 GK_FK_CK=0  means to choose gyrokinetics.<br>
 GK_FK_CK=1  means to choose cyclokinetics in Fourier harmonic representation (CKinFH).<br>
 GK_FK_CK=2  means to choose cyclokinetics in cyclotron harmonic representation (CKinCH).<br>
 GK_FK_CK=-2 is designed for deep parallelized CKinCH, since CKinCH is extremely expensive.<br><br>


*************************************************************************************************************;


 The input parameters should be set in rCYCLO input file: inputRCYCLO.txt .<br>
 The control variables are introduced below, other parameters are explained in input file or Ref [Z-W]<br>

 restart<br>
              * rCYCLO code has the ability to restart. The code will back up the basic information at each<br>
              * backup point. After the case running finished or stopped, we can restart rCYCLO either from<br>
              * the last backup point or rollback to the point before last backup point.<br>
              * restart=0 for run a case from time=0;<br>
              * restart=1 for restart from the last backup point;<br>
              * restart=2 for restart and rollback to the point before last backup point.<br><br>
 GK_FK_CK<br>
              * 0 for GK; 1 for CKinFH; 2 for CKinCH ; -2 for deep parallelized CKinCH<br>
              * Since the parallelization of rCYCLO is entirely relied on the dimensions, the number of the<br>
              * job processors is already determined when the grid variables are fixed:<br>
              * If GK_FK_CK=0. The number of processors should be: N_mu+1<br>
              * If GK_FK_CK=1. The number of processors should be: (N_mu+1)*(2*N_FT-1)<br>
              * If GK_FK_CK=2. The number of processors should be: (N_mu+1)*(2*N_CY-1)<br>
              * If GK_FK_CK=3. The number of processors should be: (2*pmax-1)*(N_mu+1)*(2*N_CY-1)<br><br>
 CDW <br>
              * rCYCLO offers two electron descriptions which is controlled by CDW parameter:<br>
              * CDW=0 for i*delta electron;<br>
              * CDW=1 for collisional drift wave (CDW) electron<br><br>
 muDtype<br>
              * muDtype is used for choosing mu-derivative operator.<br>
              * muDtype=0 is recommended since it is the one been proved to conserve the incremental entropy [Z-W].<br><br>
 mugridtype<br>
              * mugridtype is used for choosing the mu grid and the corresponding weight. Mugridtype=1 (equal<br>
              * mod weight grid) is recommended since it is most efficient. 0 for Gauss-Legendre grid. 2 for<br>
              * equal chi grid...<br><br>
 output_step<br>
              * output interval of diagnosing<br><br>
 backup_num<br>
              * The number of backups<br><br>
 Const_NL<br>
              * The coefficient in front of NL term, normally  Const_NL=1<br><br>
 Const_L<br>
              * The coefficient in front of Linear term, normally  Const_L=1<br><br>
 verify_L<br>
              * verify_L=0 is defult value, for a high performence running<br>
              * verify_L=1 for linear verification and time diagnose of each subrutine;<br>
              * verify_L=2 for time initial growth rate convergence test, after each time step all the variables will <br>
              *            be divided by a certain number: gammax(k)*Cnorm.<br><br>
 verify_NL<br>
              * verify_NL=0 is defult value, for a high performence running and only output the basic informations<br>
              * for plotting total Chi vs time, total_D vs time, and growth rate (and frequency) vs ky<br>
              * verify_NL=1 nonlinear verification and outputing more variables;<br>
              * verify_NL=2 for give the output files Phi_kOMG.txt and GOMG.txt used for plotting freq spectrum.<br><br>


*************************************************************************************************************;


 ** The rCYCLO offers a visualized program in order to plot and analyze the data, named vuRCYCLO.<br>

 ** The standard example cases of GK, CKinFH and CKinCH are given in the folder of vuRCYCLO.<br>


*************************************************************************************************************;

  #### REFERENCES:<br>
 [1] C. Holland, A.E. White, et al., Phys. Plasmas 16, 052301 (2009)<br>
 [2] R.E. Waltz, BAPS Series II, Vol. 57, No. 12, (2012) p. 105, DI3-2<br>
 [3] R. E. Waltz and Zhao Deng, Phys. Plasmas 20, 012507 (2013)<br>
