!
! CODE NAME:
!   rCYCLO
!
!
! AUTHOR:
!  ZHAO Deng
!  email: zhao.deng@foxmail.com
!
!
! ADVISER:
!   R. E. Waltz
!
!
! PURPOSE:
!	We develop rCYCLO in order to explore the explore the missing transport near L-mode edge problems[1].
!
!
! HOW TO START:
! Before started, you should make sure there is an available lapack library.
! 1. Compile rCYCLO code, remember to connected the lapack lib. e.g.
!    For NERSC: ftn rCYCLO.f90
!    If install the lapack by yourself: mpif90 2014.5.14_rCYCLO.o -L../lib -llapack -ltmglib -lblas
!    lib is the directory where you put the '*.a' compiled lapack library source files.
! 2. Submit the job, the number of processors should follow the rules in the annotation of GK_FK_CK variable below.
!
!
! INTRODUCE:
! 	Gyrokinetic simulations of L-mode near edge tokamak plasmas with the GYRO code underpredict both
! the transport and the turbulence levels by 5 to 10 fold[1], which suggest either some important
! mechanism is missing from current gyrokinetic codes like GYRO or the gyrokinetic approximation itself
! is breaking down. It is known that GYRO drift-kinetic simulations with gyro-averaging suppressed
! recover most of the missing transport[2]. With these motivations, we developed a flux tube nonlinear
! cyclokinetic[3] code rCYCLO with the parallel motion and variation suppressed. rCYCLO dynamically follows
! the high frequency ion gyro-phase motion (with no averaging) which is nonlinearly coupled into the low frequency
! drift-waves thereby interrupting and possibly suppressing the gyro-averaging. By comparison with the corresponding
! gyrokinetic simulations, we can test the conditions for the breakdown of gyrokinetics. rCYCLO nonlinearly
! couples grad-B driven ion temperature gradient (ITG) modes and collisional fluid electron drift modes to
! ion cyclotron (IC) modes.
!
!   rCYCLO code includes four independent parts, controlled by the parameter GK_FK_CK=0, 1, 2, and -2.
! GK_FK_CK=0  means to choose gyrokinetics.
! GK_FK_CK=1  means to choose cyclokinetics in Fourier harmonic representation (CKinFH).
! GK_FK_CK=2  means to choose cyclokinetics in cyclotron harmonic representation (CKinCH).
! GK_FK_CK=-2 is designed for deep parallelized CKinCH, since CKinCH is extremely expensive.
!
!
!*************************************************************************************************************;
!
!
! The input parameters should be set in rCYCLO input file: inputRCYCLO.txt .
! The control variables are introduced below, other parameters are explained in input file or Ref [Z-W]
!
! restart      * rCYCLO code has the ability to restart. The code will back up the basic information at each
!              * backup point. After the case running finished or stopped, we can restart rCYCLO either from
!              * the last backup point or rollback to the point before last backup point.
!              * restart=0 for run a case from time=0;
!              * restart=1 for restart from the last backup point;
!              * restart=2 for restart and rollback to the point before last backup point.
! GK_FK_CK     * 0 for GK; 1 for CKinFH; 2 for CKinCH ; -2 for deep parallelized CKinCH
!              * Since the parallelization of rCYCLO is entirely relied on the dimensions, the number of the
!              * job processors is already determined when the grid variables are fixed:
!              * If GK_FK_CK=0. The number of processors should be: N_mu+1
!              * If GK_FK_CK=1. The number of processors should be: (N_mu+1)*(2*N_FT-1)
!              * If GK_FK_CK=2. The number of processors should be: (N_mu+1)*(2*N_CY-1)
!              * If GK_FK_CK=3. The number of processors should be: (2*pmax-1)*(N_mu+1)*(2*N_CY-1)
! CDW          * rCYCLO offers two electron descriptions which is controlled by CDW parameter:
!              * CDW=0 for i*delta electron;
!              * CDW=1 for collisional drift wave (CDW) electron
! muDtype      * muDtype is used for choosing mu-derivative operator.
!              * muDtype=0 is recommended since it is the one been proved to conserve the incremental entropy [Z-W].
! mugridtype   * mugridtype is used for choosing the mu grid and the corresponding weight. Mugridtype=1 (equal
!              * mod weight grid) is recommended since it is most efficient. 0 for Gauss-Legendre grid. 2 for
!              * equal chi grid...
! output_step  * output interval of diagnosing
! backup_num   * The number of backups
! Const_NL     * The coefficient in front of NL term, normally  Const_NL=1
! Const_L      * The coefficient in front of Linear term, normally  Const_L=1
! verify_L     * verify_L=0 is defult value, for a high performence running
!              * verify_L=1 for linear verification and time diagnose of each subrutine;
!              * verify_L=2 for time initial growth rate convergence test, after each time step all the variables will 
!              *            be divided by a certain number: gammax(k)*Cnorm.
! verify_NL    * verify_NL=0 is defult value, for a high performence running and only output the basic informations
!              * for plotting total Chi vs time, total_D vs time, and growth rate (and frequency) vs ky
!              * verify_NL=1 nonlinear verification and outputing more variables;
!              * verify_NL=2 for give the output files Phi_kOMG.txt and GOMG.txt used for plotting freq spectrum.
!
!
!*************************************************************************************************************;
!
!
! ** The rCYCLO offers a visualized program in order to plot and analyze the data, named vuRCYCLO.
!
! ** The standard example cases of GK, CKinFH and CKinCH are given in the folder of vuRCYCLO.
!
!
!*************************************************************************************************************;
!
! REFERENCES:
! [1] C. Holland, A.E. White, et al., Phys. Plasmas 16, 052301 (2009)
! [2] R.E. Waltz, BAPS Series II, Vol. 57, No. 12, (2012) p. 105, DI3-2
! [3] R. E. Waltz and Zhao Deng, Phys. Plasmas 20, 012507 (2013)
!
!*************************************************************************************************************;
! REVISE TIME & CONTENT:
! 2014-8-6: add annotation
!
!*************************************************************************************************************;
!
!
Module parameter

implicit none

! numerics
  integer,parameter :: singlep=selected_real_kind(6), doublep=selected_real_kind(12)
  real,parameter :: PI=3.1415926536,SMALL=1.0e-10

!!##Control_variables
  integer :: restart         !# 0 new start; 1 restart without rollback; 2 restart rollback
  integer :: GK_FK_CK        !# 0 GK ion; 1 CKinFH ion; 2 CKinCH ion
  integer :: CDW             !# 0 i*delta electron; 1 CDW electron
  integer :: muDtype         !# 0 Waltz's; 1 ZD's; 3 M12 control; 4 equal space
  integer :: mugridtype      !# 0 Guass-Legendre; 1 equal space; 2 equal Chi; 3 M12; 4 equal Mod upd
  integer :: output_step     !# Diagnose output interval
  integer :: backup_num      !# The number of backups
  real :: Const_NL           !# The coefficient in front of Nonlinear term
  real :: Const_L            !# The coefficient in front of Linear term
  integer :: verify_L        !# 1 linear verification; others do not
  integer :: verify_NL       !# 1 nonlinear verification; others do not
  integer :: stdout          !# output to display;to file "2D_tubulence_stdout.out"

!!##Grid_variables
  integer :: ntmax           !# Number of total time steps
  real :: tstep              !# time step
  real ::  kxmax, kymax      !# Max Kx, Ky
  integer :: pmax, nmax      !# Max num of Kx, Ky grid
  real :: mumax              !# Max mu for Guass-Legendre grid type
  integer :: N_mu            !# Number of mu grid
  integer :: N_CY            !# Max harmonic number of CKinCH
  integer :: N_FT            !# Max harmonic number of CKinFH

!!##Physics_variables
  real :: Omega_star         !# Ion relative cyclotron frequency
  real :: a_Ln               !# Density gradient
  real :: a_LTi              !# Ion temperature gradient
  real :: lambda_k
  real :: lambda_n           !# Drift wave inertia
  real :: lambda_0           !# Zonal flow inertia
  real :: lambda_D           !# Debye length
  real :: AlphaA             !# CDW adiabaticity
  real :: delta_1            !# Electron non-adiabatic response
  real :: mu_HK              !# Damping of high-k modes
  real :: mu_LK              !# Damping of low-k modes
  real :: F_k_int            !# Initial value
  real :: Epsilon            !# aspect ratio
  real :: ratio_TiTe         !# Ratio of Ti/Te
  real :: G                  !# Drift-kinetic coefficient, added in Bessel function

!!##Trial_variables
  integer :: isotropic       !# 0 real Physics; 1 delta_k and Omiga* isotropic
  integer :: Add_GAM         !# 1  add GAM; 0 doesn't
  real :: F_GAM_int          !# GAM initial value
  real,dimension(:,:),allocatable :: nu_DWcode
  real,dimension(:),allocatable :: nu_ZFcode
  real :: Omega_GAM          !# GAM frequency
  real :: eta
  real :: nu_DW              !# Const damping to all modes
  real :: nu_ZF              !# Damping of Zonal modes
  real :: CDWid              !# CDW i*delta coefficient
  real :: q                  !#
  real :: Beta               !#
  real :: NLi                !#
  real :: NLe                !#
  real :: kpar               !#
  real :: uD                 !#
  real :: gamE               !#
  real :: gamIC              !# The ion cyclotron (IC) growth rate.
  real :: D_IC
  real,dimension(:,:),allocatable :: gamIC_k
  integer :: CHpola
  integer :: Stopnt

!!Other parameters;
  integer :: p,n,a,i,j,j_p,j_m,l,l_p,l_m,jL,lL
  real FM_integ
  integer :: ntstart
  integer :: ntmaxstep
  integer :: ntstarttemp
  integer :: ntmaxtemp
  integer :: ntmaxcount
  integer :: nt_Odd
  integer :: nt_Even
  integer :: ierr
  real :: MDmu12

!cyclo-kinetic
  integer :: a_p, a_m
  integer :: jp,lp,lpp
  integer :: p1, n1, p2, n2, ap,ip
!4D Fourier harmonic
  integer :: kk
  integer :: kmax
  integer,dimension(:),allocatable :: kk_ngrid
  real,dimension(:),allocatable :: k
  real,dimension(:),allocatable :: energy_k
  real,dimension(:),allocatable ::  energy_kk
  complex(singlep) :: NLkk1
! Trial_variables:
  integer :: Ccos
  real :: Alpha
  real :: k0
  real :: A0
  !varialble for diagnose
!  real :: total_energy_zhu
  integer ::dir_dia_count
  integer ::indir_dia_count
  !control parameter for output
  integer :: nt
  integer :: ki,kimax
  DOUBLE PRECISION :: tempOMG
  integer :: NLcount
  real :: NLt1,NLtb1,NLte1
  real :: Loopt,Looptb,Loopte
  real :: Mt,Mtb,Mte, Mtm1,Mtm2,Mt1,Mt2,Mt3
  real :: DBt,DBte,DBtb
  real :: CtNL,CtM,CtNLb,CtNLe,CtMb,CtMe
  real :: GKMt,GKMtb,GKMte,GKNLt,GKNLtb,GKNLte
  real :: loade, loadb
  real :: n_n0
  complex(singlep) :: n_n0com
  real :: Cnorm

!=================================================================================================================================
!Public
 ! varialble array
  real,dimension(:),allocatable :: kx, k1x, k2x
  real,dimension(:),allocatable :: ky, k1y, k2y
  real,dimension(:,:),allocatable :: delta_k
  ! varialble for integratioin of mu
  real,dimension(:,:),allocatable:: Integ_den
  real,dimension(:,:),allocatable:: signk

!math function
  complex(singlep),dimension(:,:),allocatable:: Integ_num,Integ_numcy  !MPI
  ! arrays mu_point and w_point of length N_mu, containing the abscissas and weights of
  !  the Gauss-Legendre n-point quadrature formula
  real,dimension(:),allocatable :: mu_Point,w_point
  real,dimension(:,:,:),allocatable:: krho
  real,dimension(:),allocatable:: rho  ! krho=k*rho
  real,dimension(:),allocatable ::FM   !Maxwell distribution -- FM(mu)

!=================================================================================================================================
!GK
  real :: total_D
  real :: total_D_3G
  real :: total_Chi
  complex(singlep),dimension(:,:),allocatable::phi_k
  complex(singlep),dimension(:,:),allocatable:: F_k,nonlinear_term !MPI
  complex(singlep),dimension(:,:,:),allocatable:: F_kmu
  complex(singlep),dimension(:,:),allocatable::H_m_delta
  complex(singlep),dimension(:,:,:),allocatable::Omega_matr
  complex(singlep),dimension(:,:,:),allocatable :: F_kGKrec
  complex(singlep),dimension(:,:),allocatable :: F_GAMrec
  complex(singlep),dimension(:),allocatable::F_GAM,F_GAMa  !MPI
  real(singlep),dimension(:,:,:),allocatable ::Fk_re_rt,Fk_im_rt !MPI
  real(singlep),dimension(:,:),allocatable ::F_GAM_re_rt,F_GAM_im_rt !MPI
  real,dimension(:,:),allocatable :: D_kx_ky,D_3G_kx_ky,Chi_kx_ky
  complex(doublep),dimension(:,:,:,:),allocatable:: R_Matrix,M_Matrix
  complex(singlep),dimension(:,:,:),allocatable:: NLSa
  complex(singlep),dimension(:),allocatable:: Source
  real(singlep),dimension(:,:),allocatable :: energy_kx_ky,entropy_kx_ky
  real(singlep),dimension(:),allocatable :: ene_GAM,Chi_mu,Fk_mu
  complex(singlep),dimension(:,:),allocatable:: F_GAMdia
  complex(singlep),dimension(:),allocatable:: Fk_mucom

!=================================================================================================================================
!CKinFH
  !intermediate result
  complex(singlep),dimension(:,:),allocatable::phi_km,phi_kp  !applied in nonlinear terms to store some intermediate result
  complex(singlep),dimension(:,:),allocatable::Cphi_m,Cphi_p
  real,dimension(:,:),allocatable:: MDmu
!4D Fourier harmonic
  complex(singlep),dimension(:,:),allocatable::phi_kft
 !complex(singlep),dimension(:,:),allocatable:: F_kft,G_kft
  complex(singlep),dimension(:,:,:),allocatable:: F_kfta, NLSfta
  complex(singlep),dimension(:),allocatable:: Sourceft
  complex(singlep),dimension(:,:),allocatable:: NLSft,NLSmuft,NLSalft
  complex(singlep),dimension(:,:),allocatable:: N_NLSmuft,N_NLSalft
  complex(singlep),dimension(:,:),allocatable:: lamb_m_delta
  real(singlep),dimension(:,:,:),allocatable:: F_kft_re_rt,F_kft_im_rt
  complex(doublep),dimension(:,:,:),allocatable:: R_Matrixft,M_Matrixft
  complex(singlep),dimension(:,:),allocatable:: muFftp, muFftm, alFftp, alFftm
  complex(singlep),dimension(:,:),allocatable:: N_muFftp, N_muFftm, N_alFftp, N_alFftm
  complex(singlep),dimension(:,:,:),allocatable:: Omega_matrft
  complex(singlep),dimension(:,:),allocatable:: Omegaft
  complex(singlep),dimension(:,:,:),allocatable:: Omegakmft
  complex(doublep),dimension(:,:),allocatable:: Matrtemp
  complex(doublep),dimension(:,:),allocatable:: MatrtempREV
  complex(doublep),dimension(:,:),allocatable:: MatrtempOUT
  complex(singlep),dimension(:,:,:),allocatable:: Fkkmft
  complex(singlep),dimension(:),allocatable::G_m_deltaft
  complex(singlep),dimension(:,:),allocatable:: Fkmyid
  real,dimension(:,:),allocatable :: GamFunft
  real(singlep),dimension(:,:),allocatable :: D_kft,Chi_kft
  real(singlep),dimension(:),allocatable :: Chi_muft,Fk_ftmu,Fk_ftmuT
  complex(singlep),dimension(:),allocatable:: Fk_ftmucom, Fk_ftcomT
  complex(singlep),dimension(:,:),allocatable:: F_kftbef
  real,dimension(:,:),allocatable:: Rate
  real,dimension(:),allocatable:: Ratave,delRat
  real,dimension(:,:),allocatable:: gammax

!=================================================================================================================================
!CKinCH
!cyclo-kinetic
  real,dimension(:,:,:,:),allocatable ::Jn
  real,dimension(:,:,:),allocatable :: GamFun
  real,dimension(:,:),allocatable :: Omega_d,Omega_nT
  complex(singlep),dimension(:,:),allocatable::G_m_delta
  complex(singlep),dimension(:,:,:,:,:,:),allocatable::Delta_nnp
  complex(singlep),dimension(:,:),allocatable:: expibp,expibm
  complex(singlep),dimension(:,:,:),allocatable:: Omega_matrcy
  complex(singlep),dimension(:),allocatable:: Sourcecy
  complex(singlep),dimension(:,:),allocatable:: Omegacy
!functions
  complex(singlep),dimension(:,:),allocatable::phi_kcy
  complex(singlep),dimension(:,:,:),allocatable:: F_kcy,G_kcy,F_kcya
  real(singlep),dimension(:,:,:),allocatable:: F_kcy_re_rt,F_kcy_im_rt
  complex(doublep),dimension(:,:,:),allocatable:: R_Matrixcy,M_Matrixcy
  complex(singlep),dimension(:,:,:),allocatable:: Omegakmcy
  complex(singlep),dimension(:,:,:),allocatable:: Fkkmcy
  complex(singlep),dimension(:,:),allocatable:: NLScy
  complex(singlep),dimension(:,:,:),allocatable:: NLScya
  complex(singlep),dimension(:,:,:),allocatable:: NLScyatemp
  complex(singlep),dimension(:,:),allocatable:: Fkmyidcy
  complex(singlep),dimension(:,:,:,:,:),allocatable:: NL1
  complex(singlep),dimension(:,:,:,:,:,:),allocatable:: NL2
  real(singlep),dimension(:,:),allocatable :: D_kcy,Chi_kcy
  complex(singlep),dimension(:,:),allocatable:: GOMG
  real(singlep),dimension(:),allocatable :: Chi_mucy,Fk_cymuT
  complex(singlep),dimension(:),allocatable:: Fk_cycomT
!=================================================================================================================================
!CKinCH deeply paralized
  complex(singlep),dimension(:,:,:,:),allocatable::Delta_nnpdp
  complex(singlep),dimension(:,:,:,:),allocatable:: NL1dp
  complex(singlep),dimension(:,:,:,:,:),allocatable:: NL2dp
  complex(singlep),dimension(:),allocatable:: NLScydp
  complex(singlep),dimension(:,:),allocatable:: NLScyadp
  complex(singlep),dimension(:,:),allocatable:: NLScyadptemp

!=================================================================================================================================
end Module parameter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main   !MPI

  use mpi

  use parameter
  implicit none
  integer myid, numprocs,rc !MPI
  real ctime0,ctime1,ctime
  character(len=79):: read_filename,plot_name_temp
  integer ::error
!  logical alive
  integer,parameter:: fileid=1001
  character(8):: date

 !MPI initial
 call MPI_INIT( ierr )
 call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
 call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
 print *, 'Process ', myid, ' of ', numprocs, ' is alive'

   stdout=1
  !read input data from input files
   read_filename='inputRCYCLO.txt'
  open(fileid,file=read_filename)
    plot_name_temp='aaa'
    Do While (.Not.(plot_name_temp == '##Control_variables'))
    read(fileid,"(A19)" ,iostat=error)plot_name_temp
    plot_name_temp=Trim(plot_name_temp)
    End Do
    read(fileid,"(8x,i2)" ,iostat=error)restart
    read(fileid,"(9x,i2)" ,iostat=error)GK_FK_CK
    read(fileid,"(4x,i2)" ,iostat=error)CDW
    read(fileid,"(8x,i2)" ,iostat=error)muDtype
    read(fileid,"(11x,i2)" ,iostat=error)mugridtype
    read(fileid,"(12x,i6)" ,iostat=error)output_step
    read(fileid,"(11x,i4)" ,iostat=error)backup_num
    read(fileid,"(9x,f6.3)" ,iostat=error)Const_NL
    read(fileid,"(8x,f6.3)" ,iostat=error)Const_L
    read(fileid,"(9x,i2)" ,iostat=error) verify_L
    read(fileid,"(10x,i2)" ,iostat=error) verify_NL
    read(fileid,"(7x,i10)" ,iostat=error)Stopnt


    plot_name_temp='aaa'
    Do While (.Not.(plot_name_temp == '##Grid_variables'))
    read(fileid,"(A16)" ,iostat=error)plot_name_temp
    plot_name_temp=Trim(plot_name_temp)
    End Do
    read(fileid,"(6x,i10)" ,iostat=error)ntmax
    read(fileid,"(6x,f12.9)" ,iostat=error)tstep
    read(fileid,"(6x,f10.5,/,6x,f10.5)" ,iostat=error) kxmax, kymax
    read(fileid,"(5x,i5,/,5x,i5)" ,iostat=error)pmax,nmax
    read(fileid,"(6x,f8.3,/,5x,i5)" ,iostat=error)mumax,N_mu
    read(fileid,"(5x,i6)" ,iostat=error)N_FT
    read(fileid,"(5x,i6)" ,iostat=error)N_CY


    plot_name_temp='aaa'
    Do While (.Not.(plot_name_temp == '##Physics_variables'))
    read(fileid,"(A19)" ,iostat=error)plot_name_temp
    plot_name_temp=Trim(plot_name_temp)
    End Do
    read(fileid,"(11x,f7.4)" ,iostat=error)Omega_star
    read(fileid,"(5x,f11.6,/,6x,f11.6)" ,iostat=error)a_Ln,a_LTi
    read(fileid,"(9x,f6.3,/,9x,f6.3,/,9x,f8.5)" ,iostat=error)lambda_n,lambda_0,lambda_D
    read(fileid,"(7x,f9.3)" ,iostat=error)AlphaA
    read(fileid,"(8x,f6.3)" ,iostat=error)delta_1
    read(fileid,"(6x,f6.3)" ,iostat=error)mu_HK
    read(fileid,"(6x,f7.4)" ,iostat=error)mu_LK
    read(fileid,"(8x,f12.7)" ,iostat=error)F_k_int
    read(fileid,"(8x,f10.7)" ,iostat=error)Epsilon
    read(fileid,"(11x,f7.3)" ,iostat=error)ratio_TiTe
    read(fileid,"(2x,f7.4)" ,iostat=error)G
    read(fileid,"(6x,f9.4)" ,iostat=error)gamIC
    read(fileid,"(5x,f9.4)" ,iostat=error)D_IC    


    plot_name_temp='aaa'
    Do While (.Not.(plot_name_temp == '##Trial_variables'))
    read(fileid,"(A17)" ,iostat=error)plot_name_temp
    plot_name_temp=Trim(plot_name_temp)
    End Do
    read(fileid,"(10x,i2)" ,iostat=error)isotropic
    read(fileid,"(8x,i2)" ,iostat=error)Add_GAM
    read(fileid,"(10x,f12.7)" ,iostat=error)F_GAM_int
    read(fileid,"(10x,f10.6)" ,iostat=error)Omega_GAM
    read(fileid,"(4x,f6.3)" ,iostat=error)eta
    read(fileid,"(6x,f7.4)" ,iostat=error)nu_DW
    read(fileid,"(6x,f7.4)" ,iostat=error)nu_ZF
    read(fileid,"(6x,f9.6)" ,iostat=error)CDWid
    read(fileid,"(5x,i6)" ,iostat=error)Ccos
    read(fileid,"(3x,f10.5)" ,iostat=error)A0
    read(fileid,"(3x,f9.4)" ,iostat=error)k0
    read(fileid,"(6x,f10.5)" ,iostat=error)Alpha
    read(fileid,"(5x,f9.4)" ,iostat=error)Beta
    read(fileid,"(4x,f9.4)" ,iostat=error)NLi
    read(fileid,"(4x,f9.4)" ,iostat=error)NLe
    read(fileid,"(7x,f9.4)" ,iostat=error)MDmu12
    read(fileid,"(5x,f7.4)" ,iostat=error)kpar
    read(fileid,"(3x,f9.4)" ,iostat=error)uD
    read(fileid,"(5x,f9.4)" ,iostat=error)gamE
    read(fileid,"(7x,i2)" ,iostat=error)CHpola
    read(fileid,"(6x,f12.4)" ,iostat=error)Cnorm
  close(fileid)


  If(myid==0)then
  !print output data on screen
    write(*,*)" "
    write(*,*)"##Control_variables:"
    write(*,*)"restart=",restart
    write(*,*)"GK_FK_CK=",GK_FK_CK
    write(*,*)"CDW=",CDW
    write(*,*)"muDtype=",muDtype
    write(*,*)"mugridtype=",mugridtype
    write(*,*)"output_step=",output_step
    write(*,*)"backup_num=",backup_num
    write(*,*)"Const_NL=",Const_NL
    write(*,*)"Const_L=",Const_L
    write(*,*)"verify_L=",verify_L
    write(*,*)"verify_NL=",verify_NL
    write(*,*)"Stopnt=",Stopnt
    write(*,*)" "
    write(*,*)" "


    write(*,*)"##Grid_variables:"
    write(*,*)"ntmax=",ntmax
    write(*,*)"tstep=",tstep
    write(*,*)"kxmax=",kxmax
    write(*,*)"kymax=",kymax
    write(*,*)"pmax=",pmax
    write(*,*)"nmax=",nmax
    write(*,*)"mumax=",mumax
    write(*,*)"N_mu=",N_mu
    write(*,*)"N_FT=",N_FT
    write(*,*)"N_CY=",N_CY
    write(*,*)" "
    write(*,*)" "


    write(*,*)"##Physics_variables:"
    write(*,*)"Omega_star=",Omega_star
    write(*,*)"a_Ln=",a_Ln
    write(*,*)"a_LTi=",a_LTi
    write(*,*)"lambda_n=",lambda_n
    write(*,*)"lambda_0=",lambda_0
    write(*,*)"lambda_D=",lambda_D
    write(*,*)"AlphaA=",AlphaA
    write(*,*)"delta_1=",delta_1
    write(*,*)"mu_HK=",mu_HK
    write(*,*)"mu_HK=",mu_LK
    write(*,*)"F_k_int=",F_k_int
    write(*,*)"Epsilon=",Epsilon
    write(*,*)"ratio_TiTe=",ratio_TiTe
    write(*,*)"G=",G
    write(*,*)"gamIC=",gamIC
    write(*,*)"D_IC=",D_IC    
    write(*,*)" "
    write(*,*)" "


    write(*,*)"##Trial_variables:"
    write(*,*)"isotropic=",isotropic
    write(*,*)"Add_GAM=",Add_GAM
    write(*,*)"F_GAM_int=",F_GAM_int
    write(*,*)'Omega_GAM=',Omega_GAM
    write(*,*)"eta=",eta
    write(*,*)"nu_DW=",nu_DW
    write(*,*)"nu_ZF=",nu_ZF
    write(*,*) 'CDWid=',CDWid
    write(*,*)"Ccos=",Ccos
    write(*,*)"A0=",A0
    write(*,*)"k0=",k0
    write(*,*)"Alpha=",Alpha
    write(*,*)"Beta=",Beta
    write(*,*)"NLi=",NLi
    write(*,*)"NLe=",NLe
    write(*,*)"MDmu12=",MDmu12
    write(*,*)"kpar=",kpar
    write(*,*)"uD=",uD
    write(*,*)"gamE=",gamE
    write(*,*)"CHpola=",CHpola
    write(*,*)"CHpola=",CHpola
    write(*,*)"Cnorm=",Cnorm
    write(*,*)" "
  Endif

  IF(myid==0)then
    If(GK_FK_CK==0)then
      if(numprocs==N_mu+1)then
      else
      write(*,*)"!!!"
      write(*,*)"In Gyro-Kinetic case: The number of processors should be equal to N_mu+1"
      write(*,*)" "
      stop
      endif
     Else if(GK_FK_CK==1)then
      if(numprocs==(N_mu+1)*(2*N_FT-1))then
      else
      write(*,*)"!!!"
      write(*,*)"In CKinFH case: The number of processors should be equal to (N_mu+1)*(2*N_FT-1)"
      write(*,*)" "
      stop
      endif
     Else if(GK_FK_CK==2)then
      if(numprocs==(N_mu+1)*(2*N_CY-1))then
      else
      write(*,*)"!!!"
      write(*,*)"In CKinCH case: The number of processors should be equal to (N_mu+1)*(2*N_CY-1)"
      write(*,*)" "
      stop
      endif
     Else if(GK_FK_CK== -2)then
      if(numprocs==(2*pmax-1)*(N_mu+1)*(2*N_CY-1))then
      else
      write(*,*)"!!!"
      write(*,*)"In this case: The number of processors should be equal to (2*pmax-1)*(N_mu+1)*(2*N_CY-1)"
      write(*,*)" "
      stop
      endif
     Else
      write(*,*)"!!!"
      write(*,*)"parameter 'GK_FK_CK' should only be 0, 1, 2 or -2."
      write(*,*)" "
      stop
     Endif
  Endif

   !We set the backup_num below:(backup_num is the pieces of the devivied ntmax,
   ! we could avoid the super large variable which may exceed the limit)
  if(backup_num.EQ.0)then
    backup_num=1 !defalt value of backup_num
    if(ntmax .GT. 20*1000)then
    backup_num=ntmax/1000
    endif
    do j=1,20
      if((pmax .GE. j*10).and.(pmax.LE.(j+1)*10))then
      backup_num=backup_num*j*10
      endif
      if((nmax .GE. j*10).and.(nmax.LE.(j+1)*10))then
      backup_num=backup_num*j*10
      endif
    enddo
    if((ntmax.LE.100).and.((pmax.GT.10).or.(nmax.GT.10)))then
    backup_num=ntmax
    endif
  endif


  IF(myid==0)then
  call CPU_TIME(ctime0)
  IF(verify_L==1)then
   call CPU_TIME(loadb)
  Endif
  ENdif

  call load(myid)

  IF(myid==0)then
  IF(verify_L==1)then
   call CPU_TIME(loade)
   write(*,*)"Load time= ",loade-loadb
  Endif
  Endif

!------------------------------------------------------------------------------------------------------------
  if(restart==0)then
  call diagnose(myid)  !output the initial value of energy, D, power at nt=0
  endif


!the loops of each ntmax pieces
DO  ntmaxcount=1,backup_num
    if(restart==0)then
    ntstarttemp=(ntmaxcount-1)*ntmaxstep     !ntstarttemp is the starting point of each loop
    ntmaxtemp=ntmaxcount*ntmaxstep           !ntmaxtemp is the ending point of each loop
    else if(restart==1.or.restart==2)then
    ntstarttemp=(ntmaxcount-1)*ntmaxstep+ntstart
    ntmaxtemp=ntmaxcount*ntmaxstep+ntstart
    endif

    IF(ntmaxcount==backup_num)then
    ntmaxtemp=ntmax
    Endif
    !the first bin have an extra zero point.
    !take backup_num=10 for example: the first bin from 0 to 10, the second bin from 11 to 20.
    !We put this plus 1 after load because of a more conviniet to write the judgement:ntstarttemp==ntstart

    ntstarttemp=ntstarttemp+1

    do nt=ntstarttemp,ntmaxtemp
         call nonlinear(myid)    ! calculate non-linear term
           if(mod(nt,output_step)==0)then  !use F_k, phi_k, nonlinear_term to caculate the growthrate and frequency of
             call indir_diagnose(myid)                        !time when mod(nt,output_step)==0.
           endif
         call motion(myid)  ! calculate the final resault---the value of N_k
           if(mod(nt,output_step)==0)then  !every output_step make an output
             call diagnose(myid)
           endif
!           if(nt==ntmax)then !caculate the indir_diagnose at the last time point.
!             call nonlinear(myid)
!             call indir_diagnose(myid)
!           endif
    enddo

    nt=nt-1
   call backup(myid)
    nt=nt+1

enddo

!------------------------------------------------------------------------------------------------------------

  if(myid==0 .and. verify_L==1)then
  call CPU_TIME(Loopte)
  Loopt=Loopte-Looptb
   write(*,*)"Loopt= ",Loopt
   if(GK_FK_CK==0)then
   open(stdout+103,file="cpu_timeGK.txt",status='old',POSITION='APPEND')
   else IF(GK_FK_CK==1)then
   open(stdout+103,file="cpu_timeFK.txt",status='old',POSITION='APPEND')
   else IF(GK_FK_CK==2 .or. GK_FK_CK== -2)then
   open(stdout+103,file="cpu_timeCK.txt",status='old',POSITION='APPEND')
   endif
   write(stdout+103,*)"Loopt= ",Loopt
   close(stdout+103)
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  IF(myid==0)then
! obtain ending CPU time
  call CPU_TIME(ctime1)
  call DATE_AND_TIME(DATE=date)
  ctime=ctime1-ctime0
  IF(GK_FK_CK==0)then
  open(stdout+103,file="cpu_timeGK.txt",status='old',POSITION='APPEND')
  Else IF(GK_FK_CK==1)then
  open(stdout+103,file="cpu_timeFK.txt",status='old',POSITION='APPEND')
  Else IF(GK_FK_CK==2 .or. GK_FK_CK== -2)then
  open(stdout+103,file="cpu_timeCK.txt",status='old',POSITION='APPEND')
  Endif
  write(stdout+103,*)'Program CPU TIME=', ctime
  write(stdout+103,*)"Date=",date
  close(stdout+103)
  Endif



  write(*,*)"End Process ",myid," ...."

  call MPI_FINALIZE(rc)


end   !MPI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




Subroutine load(myid)

  use mpi
  use parameter
  implicit none
  integer,intent(in) :: myid
  integer info
  real BESSJ,BESSI
  real b
  integer imax
  integer nc
  real,ALLOCATABLE,DIMENSION(:) :: mu_pointtemp,w_pointtemp,mu_bar
  real,ALLOCATABLE,DIMENSION(:,:) :: MDmum
  real,ALLOCATABLE,DIMENSION(:,:) :: MDmup
  real(doublep) veritemp,veritempd
  real alp,alpp, Area
  real,ALLOCATABLE,DIMENSION(:,:) :: random
  real,ALLOCATABLE,DIMENSION(:) :: randomGAM

  complex*16, ALLOCATABLE, DIMENSION(:,:) :: Matrix
  complex*16, ALLOCATABLE, DIMENSION(:) :: eigenO
  complex*16, ALLOCATABLE, DIMENSION(:) :: WORK
  INTEGER :: nmatr
  complex*16  DUMMY(1,1)

 !cyclo-kinetic
  integer :: nmatrcy
  complex*16, ALLOCATABLE, DIMENSION(:,:) :: Matrixcy
  complex*16, ALLOCATABLE, DIMENSION(:) :: eigenOcy
  complex*16, ALLOCATABLE, DIMENSION(:) :: WORKcy
  character(len=30):: Omegacyname

  !4D Fourier harmonic
  integer :: jj
  integer :: nmatrft
  complex*16, ALLOCATABLE, DIMENSION(:,:) :: Matrixft
  complex*16, ALLOCATABLE, DIMENSION(:) :: eigenOft
  complex*16, ALLOCATABLE, DIMENSION(:) :: WORKft
  character(len=30):: Omegaftname
  integer, ALLOCATABLE, DIMENSION(:) :: IPIV

!=================================================================================================================================
!Public
  allocate(MDmum(1:N_mu,1:N_mu),MDmup(1:N_mu,1:N_mu))
  allocate(random(-pmax+1:pmax-1,-nmax+1:nmax-1),randomGAM(-pmax+1:pmax-1))
  allocate(Jn(-pmax+1:pmax-1,-nmax+1:nmax-1,1:N_mu,min(-2*N_CY,-N_FT):max(2*N_CY,N_FT)))
  allocate(kx(-pmax+1:pmax-1),k1x(-pmax+1:pmax-1),k2x(-pmax+1:pmax-1))
  allocate(ky(-nmax+1:nmax-1),k1y(-nmax+1:nmax-1),k2y(-nmax+1:nmax-1))
  allocate(delta_k(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(nu_DWcode(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(nu_ZFcode(-nmax+1:nmax-1))
  allocate(krho(-pmax+1:pmax-1,-nmax+1:nmax-1,1:N_mu),rho(1:N_mu))
  allocate(FM(1:N_mu))
  allocate(signk(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(mu_point(1:N_mu),w_point(1:N_mu))
  allocate(mu_bar(1:N_mu-1))
  allocate(mu_pointtemp(1:N_mu-1),w_pointtemp(1:N_mu-1))
  allocate(Integ_den(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(Integ_num(-pmax+1:pmax-1,0:nmax-1))  !MPI
  allocate(Omega_nT(-nmax+1:nmax-1,1:N_mu))

IF(GK_FK_CK==1 .or. GK_FK_CK==2 .or. GK_FK_CK==-2)then
  allocate(Omega_d(-nmax+1:nmax-1,1:N_mu))
  allocate(expibp(-pmax+1:pmax-1,-nmax+1:nmax-1),&
           expibm(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(MDmu(1:N_mu,1:N_mu))
  allocate(gamIC_k(-pmax+1:pmax-1,0:nmax-1))
  !intermediate varables
EndIf
!=================================================================================================================================
!GK
IF(GK_FK_CK==0)then
  !distribute memory to each allocated variables
  allocate(phi_k(-pmax+1:pmax-1,-nmax+1:nmax-1),F_k(-pmax+1:pmax-1,-nmax+1:nmax-1)) !MPI
 ! allocate(F_ka(-pmax+1:pmax-1,-nmax+1:nmax-1),phi_ka(-pmax+1:pmax-1,-nmax+1:nmax-1)) !MPI
  allocate(F_kmu(-pmax+1:pmax-1,-nmax+1:nmax-1,0:N_mu))
  allocate(nonlinear_term(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(D_kx_ky(-pmax+1:pmax-1,-nmax+1:nmax-1)&
           ,D_3G_kx_ky(-pmax+1:pmax-1,-nmax+1:nmax-1)&
           ,Chi_kx_ky(-pmax+1:pmax-1,-nmax+1:nmax-1))
 ! allocate(F_kGKrec(-pmax+1:pmax-1,-nmax+1:nmax-1,0:N_mu))
  allocate(F_GAMrec(-pmax+1:pmax-1,1:N_mu))

  allocate(energy_kk(1:(pmax-1)**2+(nmax-1)**2))
  allocate(kk_ngrid(1:(pmax-1)**2+(nmax-1)**2))
  IF(Add_GAM==1)then
  allocate(F_GAM(-pmax+1:pmax-1),F_GAMa(-pmax+1:pmax-1))
  allocate(F_GAM_re_rt(-pmax+1:pmax-1,1:N_mu),&
            F_GAM_im_rt(-pmax+1:pmax-1,1:N_mu))
  F_GAM=(0.0,0.0)
  Endif

  if(myid==0)then
  allocate(H_m_delta(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(NLSa(-pmax+1:pmax-1,-nmax+1:nmax-1,0:N_mu))
  allocate(Source(0:N_mu))
  allocate(Fk_re_rt(-pmax+1:pmax-1,-nmax+1:nmax-1,0:N_mu),&
            Fk_im_rt(-pmax+1:pmax-1,-nmax+1:nmax-1,0:N_mu))
  allocate(energy_kx_ky(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(entropy_kx_ky(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(ene_GAM(-pmax+1:pmax-1))
  allocate(Chi_mu(1:N_mu))
  allocate(Fk_mu(1:N_mu))
  allocate(F_GAMdia(-pmax+1:pmax-1,1:N_mu))
  allocate(Fk_mucom(1:N_mu))

  open(stdout+103,file="cpu_timeGK.txt",status='replace')
  write(stdout+103,*)"cpu_timeGK.txt"
  close(stdout+103)
  endif

Endif
!=================================================================================================================================
!CKinFH
IF(GK_FK_CK==1)then
!4D Fourier harmonic
  allocate(phi_kft(-pmax+1:pmax-1,-nmax+1:nmax-1))
  !allocate(F_kft(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(F_kfta(-pmax+1:pmax-1,-nmax+1:nmax-1,0:(N_mu+1)*(2*N_FT-1)-1))
  allocate(NLSfta(-pmax+1:pmax-1,0:nmax-1,0:(N_mu+1)*(2*N_FT-1)-1))
 ! allocate(G_kft(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(muFftp(-pmax+1:pmax-1,-nmax+1:nmax-1),&
           muFftm(-pmax+1:pmax-1,-nmax+1:nmax-1),&
  		   alFftp(-pmax+1:pmax-1,-nmax+1:nmax-1),&
  		   alFftm(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(N_muFftp(-pmax+1:pmax-1,-nmax+1:nmax-1),&
           N_muFftm(-pmax+1:pmax-1,-nmax+1:nmax-1),&
  		   N_alFftp(-pmax+1:pmax-1,-nmax+1:nmax-1),&
  		   N_alFftm(-pmax+1:pmax-1,-nmax+1:nmax-1))
  !Math
  if( mod((2*pmax-1)*nmax,(N_mu+1)*(2*N_FT-1))==0 )then
  kimax=((2*pmax-1)*nmax) / ((N_mu+1)*(2*N_FT-1))
  else
  kimax=((2*pmax-1)*nmax) / ((N_mu+1)*(2*N_FT-1)) +1
  endif
  allocate(NLSft(-pmax+1:pmax-1,0:nmax-1))
  allocate(NLSmuft(-pmax+1:pmax-1,0:nmax-1))
  allocate(NLSalft(-pmax+1:pmax-1,0:nmax-1))
  allocate(N_NLSmuft(-pmax+1:pmax-1,0:nmax-1))
  allocate(N_NLSalft(-pmax+1:pmax-1,0:nmax-1))
  allocate(Sourceft(0:(N_mu+1)*(2*N_FT-1)-1))
  allocate(GamFunft(0:kimax-1,-N_CY+1:N_CY-1))
  allocate(G_m_deltaft(0:kimax-1))
  allocate(Fkmyid(0:kimax-1,0:(N_mu+1)*(2*N_FT-1)-1))
  allocate(phi_km(-pmax+1:pmax-1,-nmax+1:nmax-1),phi_kp(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(Cphi_m(-pmax+1:pmax-1,-nmax+1:nmax-1),Cphi_p(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(lamb_m_delta(-pmax+1:pmax-1,-nmax+1:nmax-1))

 if(myid==0)then
  allocate(F_kft_re_rt(-pmax+1:pmax-1,-nmax+1:nmax-1,0:(N_mu+1)*(2*N_FT-1)-1),&
            F_kft_im_rt(-pmax+1:pmax-1,-nmax+1:nmax-1,0:(N_mu+1)*(2*N_FT-1)-1))
  allocate(Fkkmft(0:kimax-1,0:(N_mu+1)*(2*N_FT-1)-1,0:(N_mu+1)*(2*N_FT-1)-1))
  allocate(D_kft(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(Chi_kft(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(Chi_muft(1:N_mu))
  allocate(Fk_ftmu(1:N_mu))
  allocate(Fk_ftmuT(1:N_mu))
  allocate(Fk_ftmucom(1:N_mu))
  allocate(Fk_ftcomT(1:N_mu))
  if(verify_L==2)then
  allocate(F_kftbef(-nmax+1:nmax-1,0:(N_mu+1)*(2*N_FT-1)-1))
  allocate(Rate(1:nmax-1,0:N_mu*(2*N_FT-1)-1))
  allocate(Ratave(1:nmax-1))
  allocate(delRat(1:nmax-1))
  allocate(gammax(-pmax+1:pmax-1,-nmax+1:nmax-1))
  endif

  open(stdout+103,file="cpu_timeFK.txt",status='replace')
  write(stdout+103,*)"cpu_timeFK.txt"
  close(stdout+103)
 endif

  !4D Fourier harmonic
  if(N_FT < 10)then
  write(Omegaftname,"(i1)")N_FT
  else if(N_FT <100)then
  write(Omegaftname,"(i2)")N_FT
  endif
  Omegaftname=Trim('Omega_matrftN_FT' // Trim(Omegaftname) // '.txt')

Endif

!=================================================================================================================================
!CKinCH
IF(GK_FK_CK==2)then
!cyclo-kinetic
  allocate(Delta_nnp(-N_CY+1:N_CY-1,-pmax+1:pmax-1,-nmax+1:nmax-1,-pmax+1:pmax-1,0:nmax-1,-N_CY+1:N_CY-1))
!F_k and Phi_k
  allocate(phi_kcy(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(F_kcya(0:(N_mu+1)*(2*N_CY-1)-1,-pmax+1:pmax-1,-nmax+1:nmax-1))
!Math
  allocate(NLScy(-pmax+1:pmax-1,0:nmax-1))
  allocate(G_m_delta(-pmax+1:pmax-1,0:nmax-1))
  allocate(Sourcecy(0:(N_mu+1)*(2*N_CY-1)-1))
  allocate(NLScya(0:(N_mu+1)*(2*N_CY-1)-1,-pmax+1:pmax-1,0:nmax-1))

  i=mod(myid,N_mu+1)
  a=myid/(N_mu+1)-(N_CY-1)
  IF(i .gt. 0)then
    allocate(NL1(-N_CY+1:N_CY-1,-pmax+1:pmax-1,-nmax+1:nmax-1,-pmax+1:pmax-1,0:nmax-1))
    allocate(NL2(1:N_mu,-N_CY+1:N_CY-1,-pmax+1:pmax-1,-nmax+1:nmax-1,-pmax+1:pmax-1,0:nmax-1))
  Endif

  if( mod((2*pmax-1)*nmax,(N_mu+1)*(2*N_CY-1)-1)==0 )then
  kimax=((2*pmax-1)*nmax) / ((N_mu+1)*(2*N_CY-1)-1)
  else
  kimax=((2*pmax-1)*nmax) / ((N_mu+1)*(2*N_CY-1)-1) +1
  endif
  allocate(Fkmyidcy(0:(N_mu+1)*(2*N_CY-1)-1,0:kimax-1))

 if(myid==0)then
  allocate(G_kcy(0:(N_mu+1)*(2*N_CY-1)-1,-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(F_kcy_re_rt(0:(N_mu+1)*(2*N_CY-1)-1,-pmax+1:pmax-1,-nmax+1:nmax-1),&
            F_kcy_im_rt(0:(N_mu+1)*(2*N_CY-1)-1,-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(GamFun(-N_CY+1:N_CY-1,-pmax+1:pmax-1,0:nmax-1))
  allocate(Fkkmcy(0:(N_mu+1)*(2*N_CY-1)-1,0:kimax-1,0:(N_mu+1)*(2*N_CY-1)-1))
  allocate(Integ_numcy(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(NLScyatemp(-pmax+1:pmax-1,0:nmax-1,0:(N_mu+1)*(2*N_CY-1)-1))
  allocate(D_kcy(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(Chi_kcy(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(Chi_mucy(1:N_mu))
  allocate(Fk_cymuT(1:N_mu))
  allocate(Fk_cycomT(1:N_mu))

  open(stdout+103,file="cpu_timeCK.txt",status='replace')
  write(stdout+103,*)"cpu_timeCK.txt"
  close(stdout+103)

 endif

  !cyclo-kinetic
  if(N_CY < 10)then
  write(Omegacyname,"(i1)")N_CY
  else if(N_CY <100)then
  write(Omegacyname,"(i2)")N_CY
  endif
  Omegacyname=Trim('Omega_matrcyN_CY' // Trim(Omegacyname) // '.txt')

EndIf

call MPI_Barrier(MPI_COMM_WORLD, ierr)
!=================================================================================================================================
!CKinCH deeply paralized
IF(GK_FK_CK== -2)then
!cyclo-kinetic
  allocate(Delta_nnpdp(-N_CY+1:N_CY-1,-pmax+1:pmax-1,-nmax+1:nmax-1,0:nmax-1))
!F_k and Phi_k
  allocate(phi_kcy(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(F_kcya(0:(N_mu+1)*(2*N_CY-1)-1,-pmax+1:pmax-1,-nmax+1:nmax-1))
!Math
  allocate(NLScydp(0:nmax-1))
  allocate(G_m_delta(-pmax+1:pmax-1,0:nmax-1))
  allocate(Sourcecy(0:(N_mu+1)*(2*N_CY-1)-1))
  allocate(NLScyadp(0:(2*pmax-1)*(N_mu+1)*(2*N_CY-1)-1,0:nmax-1))

  i=mod(myid,N_mu+1)
  a=mod((myid-i)/(N_mu+1),(2*N_CY-1))-(N_CY-1)
  IF(i .gt. 0)then
    allocate(NL1dp(-N_CY+1:N_CY-1,-pmax+1:pmax-1,-nmax+1:nmax-1,0:nmax-1))
    allocate(NL2dp(1:N_mu,-N_CY+1:N_CY-1,-pmax+1:pmax-1,-nmax+1:nmax-1,0:nmax-1))
  Endif

  if( mod((2*pmax-1)*nmax,(2*pmax-1)*(N_mu+1)*(2*N_CY-1)-1)==0 )then
  kimax=((2*pmax-1)*nmax) / ((2*pmax-1)*(N_mu+1)*(2*N_CY-1)-1)
  else
  kimax=((2*pmax-1)*nmax) / ((2*pmax-1)*(N_mu+1)*(2*N_CY-1)-1) +1
  endif
  allocate(Fkmyidcy(0:(N_mu+1)*(2*N_CY-1)-1,0:kimax-1))

  if(myid==0)then
  allocate(G_kcy(0:(N_mu+1)*(2*N_CY-1)-1,-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(F_kcy_re_rt(0:(N_mu+1)*(2*N_CY-1)-1,-pmax+1:pmax-1,-nmax+1:nmax-1),&
            F_kcy_im_rt(0:(N_mu+1)*(2*N_CY-1)-1,-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(GamFun(-N_CY+1:N_CY-1,-pmax+1:pmax-1,0:nmax-1))
  allocate(Fkkmcy(0:(N_mu+1)*(2*N_CY-1)-1,0:kimax-1,0:(2*pmax-1)*(N_mu+1)*(2*N_CY-1)-1))
  allocate(Integ_numcy(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(NLScyadptemp(0:nmax-1,0:(2*pmax-1)*(N_mu+1)*(2*N_CY-1)-1))
  allocate(D_kcy(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(Chi_kcy(-pmax+1:pmax-1,-nmax+1:nmax-1))
  allocate(Chi_mucy(1:N_mu))
  allocate(Fk_cymuT(1:N_mu))
  allocate(Fk_cycomT(1:N_mu))

  open(stdout+103,file="cpu_timeCK.txt",status='replace')
  write(stdout+103,*)"cpu_timeCK.txt"
  close(stdout+103)

  endif

  !cyclo-kinetic
  if(N_CY < 10)then
  write(Omegacyname,"(i1)")N_CY
  else if(N_CY <100)then
  write(Omegacyname,"(i2)")N_CY
  endif
  Omegacyname=Trim('Omega_matrcyN_CY' // Trim(Omegacyname) // '.txt')

EndIf
!=================================================================================================================================
!Public
 !set the mu grids with subrutine gauss_legendre()
  mu_point(:)=0
  w_point(:)=0
  mu_pointtemp(:)=0
  w_pointtemp(:)=0

If(mugridtype==0)then
  call gauss_legendre(0.0,mumax,mu_pointtemp,w_pointtemp,N_mu-1)
  mu_point(1:N_mu-1)=mu_pointtemp(:)
  w_point(1:N_mu-1)=w_pointtemp(:)
  do i=1,N_mu-1
    w_point(i)=w_point(i)*2*PI
  enddo
  mu_point(N_mu)=mumax
  w_point(N_mu)=2*PI
Else if(mugridtype==1)then
  do i=1,N_mu-1
  mu_bar(i)=-ALOG(1.0-1.0*i/N_mu)
  enddo
  mu_point(1)=(mu_bar(1)+0)/2.0
  do i=2,N_mu-1
  mu_point(i)=(mu_bar(i)+mu_bar(i-1))/2.0
  enddo
  mu_point(N_mu)=-ALOG(1.0/N_mu)
  w_point(1)=(mu_bar(1)-0)*2*PI
  do i=2,N_mu-1
  w_point(i)=(mu_bar(i)-mu_bar(i-1))*2*PI
  enddo
  w_point(N_mu)=2*PI

Else if(mugridtype==2)then
  alp=1.5
  Do i=1,N_mu-1
    Do
    alpp=-Alog((1.0-1.0*i/N_mu)/(alp+1.0))
      if(abs(alp-alpp) .le. 0.0000001)exit
    alp=alpp
    Enddo
  mu_bar(i)=alpp
  Enddo

  mu_point(N_mu)= mu_bar(N_mu-1)
  mu_point(1)=(mu_bar(1)+0)/2.0
  do i=2,N_mu-1
  mu_point(i)=(mu_bar(i)+mu_bar(i-1))/2.0
  enddo

  w_point(1)=mu_bar(1)*2*PI
  do i=2,N_mu-1
  w_point(i)=(mu_bar(i)-mu_bar(i-1))*2*PI
  enddo
  w_point(N_mu)=2*PI

Else if(mugridtype==3)then  !construct grid. weight comes from Area/height
  do i=1,N_mu-1
  mu_bar(i)=-ALOG(1.0-1.0*i/N_mu)
  enddo
  mu_point(1)=(mu_bar(1)+0)/2.0
  mu_point(2)=(mu_bar(2)+mu_bar(1))/2.0

  w_point(1)=2.0*PI/(N_mu*exp(-mu_point(1)))
  w_point(2)=2.0*PI/(N_mu*exp(-mu_point(2)))

  mu_point(3)=mu_point(1)+w_point(2)/(w_point(1)*MDmu12)
  w_point(3)=2.0*PI/(N_mu*exp(-mu_point(3)))

  do j=3, N_mu-1
  mu_point(j+1)=(mu_point(j)-mu_point(j-2))*w_point(j)/w_point(j-1)+mu_point(j-1)
  w_point(j+1)=2.0*PI/(N_mu*exp(-mu_point(j+1)))
  enddo

Else if(mugridtype==4)then
  Area=(1-exp(-mumax))/(N_mu-1)
  do i=1,N_mu-1
  mu_bar(i)=-ALOG(1.0-i*Area)
  enddo
  mu_point(1)=(mu_bar(1)+0)/2.0
  do i=2,N_mu-1
  mu_point(i)=(mu_bar(i)+mu_bar(i-1))/2.0
  enddo
  mu_point(N_mu)=mumax

  do i=1,N_mu-1
  w_point(i)=Area/exp(-mu_point(i))*2.0*PI
  enddo
  w_point(N_mu)=2*PI


Else if(mugridtype==5)then   !construct grid. weight comes from space.

  do i=1,N_mu-1
  mu_bar(i)=-ALOG(1.0-1.0*i/N_mu)
  enddo
  mu_point(1)=(mu_bar(1)+0)/2.0
  mu_point(2)=(mu_bar(2)+mu_bar(2-1))/2.0

  w_point(1)=mu_bar(1)*2*PI
  w_point(2)=(mu_bar(2)-mu_bar(2-1))*2*PI

  mu_point(3)=mu_point(1)+w_point(2)/(w_point(1)*MDmu12)
  w_point(3)=(2*mu_point(3)*2*PI-w_point(1)-w_point(2))

  do j=3, N_mu-1
  mu_point(j+1)=(mu_point(j)-mu_point(j-2))*w_point(j)/w_point(j-1)+mu_point(j-1)
  w_point(j+1)=2*mu_point(j+1)*2*PI
   do i=1,j
   w_point(j+1)=w_point(j+1)-w_point(i)
   enddo
  enddo

  w_point(N_mu)=0
  do i=1,N_mu-1
  w_point(N_mu)=w_point(N_mu)+w_point(i)
  enddo
  w_point(N_mu)=exp(-w_point(N_mu)/(2*PI))/exp(-mu_point(N_mu))*2*PI  !Area/height


Else if(mugridtype==6)then
  do i=1,N_mu
  mu_point(i)=(i-0.5)*Mumax/N_mu
  w_point(i)=2.0*PI*Mumax/N_mu
  enddo

EndIf

  !gauss_legendre(x1,x2,x,w,n)
  If(myid==1)then
    write(*,*)"mu_point="
    do i=1,N_mu
    write(*,*)mu_point(i)
    enddo
    write(*,*)" "
    write(*,*)"w_point="
    do i=1,N_mu
    write(*,*)w_point(i)
    enddo
    write(*,*)" "
    write(*,*)"Sum_WGuass-Leg(1:N-1)=",sum(w_pointtemp(1:N_mu-1))
    write(*,*)" "
    write(*,*)"Sum_W(1:N-1)=",sum(w_point(1:N_mu-1))
    write(*,*)"Sum_W(1:N)=",sum(w_point(1:N_mu))
    write(*,*)"----------------------- "
   Endif

!Set the initial value of k, kx, ky
  do p=-pmax+1,pmax-1     !  kxmax-kxmax/pmax =< kx <=kxmax+kxmax/pmax ; kymax-kymax/pmax =< ky <= kymax-kymax/pmax
    kx(p)=kxmax*real(p)/real(pmax-1)
    k1x(p)=kxmax*real(p)/real(pmax-1)
    k2x(p)=kxmax*real(p)/real(pmax-1)
  enddo
  do n=-nmax+1,nmax-1
    ky(n)=kymax*real(n)/real(nmax-1)
    k1y(n)=kymax*real(n)/real(nmax-1)
    k2y(n)=kymax*real(n)/real(nmax-1)
  enddo

!calculate Maxwell distribution FM(mu)
  do i=1,N_mu
  FM(i)=(0.5/PI)*exp((-1.0)*mu_point(i))
  enddo

  if(myid==1)then
  FM_integ=0
  do i=1,N_mu
     FM_integ = FM_integ+w_point(i)*FM(i)
  enddo
  write(*,*)' '
  write(*,*)'Maxwell integration is',FM_integ
  write(*,*)' '
  FM_integ=0
  do i=1,N_mu
     FM_integ = FM_integ+w_point(i)*mu_point(i)*exp(-mu_point(i)**2)
  enddo
  write(*,*)'Test1=',FM_integ
  write(*,*)' '
  FM_integ=0
  do i=1,N_mu
     FM_integ = FM_integ+w_point(i)*mu_point(i)*exp(-mu_point(i))
  enddo
  write(*,*)'Test2=',FM_integ
  write(*,*)' '
  Endif

  !calculate k*rho

  do n=-nmax+1,nmax-1
  do p=-pmax+1,pmax-1
    do i=1,N_mu
       rho(i)=sqrt(2*ratio_TiTe*mu_point(i))
       krho(p,n,i)=sqrt(kx(p)**2+ky(n)**2)*rho(i)*G
    enddo
  enddo
  enddo

  do n=-nmax+1,nmax-1
  do p=-pmax+1,pmax-1
    if(n/=0) then
    signk(p,n)=sign(1.0,ky(n))
    else if (n==0)then
    signk(p,n)=sign(1.0,kx(p))
    endif
  enddo
  enddo
  !get the value of delta_k, and we should first defined signk
If(isotropic==0)then
  do n=-nmax+1,nmax-1
  do p=-pmax+1,pmax-1
    delta_k(p,n)=delta_1*10/Omega_star*ky(n)/(1.0+eta*kx(p)**2)
    if(ky(n).EQ.0.0)then
    nu_DWcode(p,n)=0
    else
    nu_DWcode(p,n)=nu_DW
    endif
    if(ky(n).EQ.0.0 )then
    nu_ZFcode(n)=nu_ZF
    else
    nu_ZFcode(n)=0
    endif
  enddo
  enddo

Else if(isotropic==1)then
 ! isotropic delta_k scheme1
 signk(:,:)=0
  do n=0,nmax-1
  do p=-pmax+1,pmax-1
    if(n/=0) then
    signk(p,n)=sign(1.0,ky(n))
    else if (n==0)then
    signk(p,n)=sign(1.0,kx(p))
    end if
    delta_k(p,n)=delta_1*sqrt(kx(p)**2 + ky(n)**2)*signk(p,n)/(1.0+eta*(kx(p)**2 + ky(n)**2))
  enddo
  enddo
  delta_k(0,0)=1
Else if(isotropic==2)then
 ! isotropic delta_k scheme2
  signk(:,:)=0
  do n=0,nmax-1
  do p=-pmax+1,pmax-1
   if (kx(p)==0.0 .and. ky(n)==0.0)then
   signk(p,n)=1.0
   elseif(ky(n) .GE. 0)then
   signk(p,n)=sign(1.0,cos(Ccos*acos(kx(p)/sqrt(kx(p)**2+ky(n)**2))))
   elseif (ky(n) .LT. 0)then
   signk(p,n)=sign(1.0,cos(Ccos*(PI+acos(kx(p)/sqrt(kx(p)**2+ky(n)**2)))))
   endif
   delta_k(p,n)=delta_1*sqrt(kx(p)**2 + ky(n)**2)*signk(p,n)/(1.0+eta*(kx(p)**2 + ky(n)**2))
  write(*,*)'(',p,n,')',signk(p,n)
  enddo
  enddo
  delta_k(0,0)=1
Endif

     do p=-pmax+1,pmax-1
       do n=1,nmax-1
       delta_k(-p,-n)= -1.0*delta_k(p,n)
       enddo
     enddo

  do n=-nmax+1,nmax-1
  do p=-pmax+1,pmax-1
    do a=0,max(2*N_CY,N_FT)
    do i=1,N_mu
    Jn(p,n,i,a)=BESSJ(a,krho(p,n,i))
    enddo
    enddo
  enddo
  enddo

  do n=-nmax+1,nmax-1
  do p=-pmax+1,pmax-1
    do a=1,max(2*N_CY,N_FT)! must to calculate the a<0 value by this way, since BESSJ() get wrong value with a<0
    do i=1,N_mu
      Jn(p,n,i,-a)=Jn(p,n,i,a)*((-1)**a)
    enddo
    enddo
  enddo
  enddo


 If(myid==0 .and. verify_L==1)then
 open(300,file="Jn.txt",status='replace')
 do a=min(-2*N_CY,-N_FT),max(2*N_CY,N_FT)
  do i=1,N_mu
    do p=-pmax+1,pmax-1
      do n=-nmax+1,nmax-1
      write(300,*)"Jn(",krho(p,n,i),",",a,")=",Jn(p,n,i,a)
      enddo
     enddo
   enddo
 enddo
 close(300)
 Endif
!load the initial value of F_k
  F_k_int=F_k_int/((2*pmax-1)*(2*nmax-1))
  F_GAM_int=F_GAM_int/(2*pmax-1)

  dir_dia_count=1  !dir_dia_count is the times of call diagnose
  indir_dia_count=1  !indir_dia_count is the times of call indir_diagnose

   do i=1,N_mu
    do n=-nmax+1,nmax-1
    Omega_nT(n,i)=ky(n)*(a_Ln+a_LTi*(mu_point(i)-1))
    enddo
   enddo

IF(GK_FK_CK==1 .or. GK_FK_CK==2 .or. GK_FK_CK== -2)then
  do n=-nmax+1,nmax-1
  do p=-pmax+1,pmax-1
    IF(n==0 .and. p==0)then
    expibp(p,n)=(1,0)
    expibm(p,n)=(1,0)
    else
    expibp(p,n)=(kx(p)/sqrt(kx(p)**2+ky(n)**2) + (0,1.0)*ky(n)/sqrt(kx(p)**2+ky(n)**2))
    expibm(p,n)=(kx(p)/sqrt(kx(p)**2+ky(n)**2) - (0,1.0)*ky(n)/sqrt(kx(p)**2+ky(n)**2))
    endif
  enddo
  enddo

   do i=1,N_mu
    do n=-nmax+1,nmax-1
    Omega_d(n,i)=ky(n)*ratio_TiTe*mu_point(i)*Epsilon
    enddo
   enddo

  do n=0,nmax-1
  do p=-pmax+1,pmax-1
    gamIC_k(p,n)=-1.778*gamIC*(ky(n)**2)+2.667*gamIC*sqrt(ky(n)**2)
  enddo
  enddo

!-mu_LK/((abs(ky(n))+0.03)*(kx(p)**2+ky(n)**2)) &


  MDmu(:,:)=0
  MDmup(:,:)=0
  MDmum(:,:)=0
  if(muDtype==0 .or. muDtype==2 .or. muDtype==-2)then
  do j=1,N_mu
    do jp=1,N_mu
    if(jp == j+1 .and. j>1)then
    MDmup(j,jp)=1/(mu_point(j+1)-mu_point(j-1))
    elseif(jp == j+1 .and. j==1 ) then
    MDmup(1,2)=1.0/(mu_point(2)-mu_point(1))

    !elseif(jp == j+1 .and. j==1 ) then
    !MDmup(1,2)=2.0/(mu_point(2)-mu_point(1))
    !elseif(jp == 1 .and. j==1 ) then
    !MDmup(1,1)= -2.0/(mu_point(2)-mu_point(1))
    endif

    if(jp == j-1 .and. j<N_mu) then
    MDmum(j,jp)= -1/(mu_point(j+1)-mu_point(j-1))
    elseif(jp == j-1 .and. j==N_mu ) then
    MDmum(N_mu,N_mu-1)=-1/(mu_point(N_mu)-mu_point(N_mu-1))

    !elseif(jp == j-1 .and. j==N_mu ) then
    !MDmum(N_mu,N_mu-1)= -2.0/(mu_point(N_mu)-mu_point(N_mu-1))
    !elseif(jp == N_mu .and. j==N_mu ) then
    !MDmu(N_mu,N_mu)= 2.0/(mu_point(N_mu)-mu_point(N_mu-1))
    endif
    enddo
  enddo

  do j=1,N_mu
    do jp=1,N_mu
    if(jp == j-1)then
    MDmup(j,jp)= -(w_point(jp))/(w_point(j))*MDmup(jp,j)
    endif
    if(jp == j+1)then
    MDmum(j,jp)= -(w_point(jp))/(w_point(j))*MDmum(jp,j)
    endif
    enddo
  enddo
  MDmu(:,:)=MDmup(:,:)/2+MDmum(:,:)/2

  else if(muDtype==1)then
  do j=1,N_mu
    do jp=1,N_mu
    if(jp == j+1 .and. j>1 .and. j<N_mu)then
    MDmu(j,jp)=1/(mu_point(j+1)-mu_point(j-1))
    elseif(jp == j+1 .and. j==1 ) then
    MDmup(1,2)=1/(mu_point(j+1)-mu_point(j))
    elseif(jp == 1 .and. j==1 ) then
    MDmu(1,1)=-1/(mu_point(j+1)-mu_point(j))
    endif

    if(jp == j-1 .and. j>1 .and. j<N_mu)then
    MDmu(j,jp)=-1/(mu_point(j+1)-mu_point(j-1))
    elseif(jp ==N_mu .and. j==N_mu ) then
    MDmu(N_mu,N_mu)=1/(mu_point(j)-mu_point(j-1))
    elseif(jp == N_mu-1 .and. j==N_mu ) then
    MDmu(N_mu,N_mu-1)=-1/(mu_point(j)-mu_point(j-1))
    endif
    enddo
  enddo

  else if(muDtype==3)then
  MDmu(:,:)=0
  MDmu(1,2)=MDmu12
 ! MDmu(1,1)=(-FM(1)-MDmu(1,2)*FM(2))/FM(1)
  do j=2,N_mu-1
  MDmu(j,j-1)=-w_point(j-1)/w_point(j)*MDmu(j-1,j)
  MDmu(j,j+1)=-MDmu(j,j-1)
  enddo
  j=N_mu
  MDmu(j,j-1)=-w_point(j-1)/w_point(j)*MDmu(j-1,j)

  else if(muDtype==4)then
  MDmu(:,:)=0
  do j=2,N_mu-1
  MDmu(j,j+1)=1.0/(2.0*Mumax/N_mu)
  MDmu(j,j-1)=-1.0/(2.0*Mumax/N_mu)
  enddo
  MDmu(1,2)=1.0/(Mumax/N_mu)
  MDmu(N_mu,N_mu-1)=-1.0/(Mumax/N_mu)

  else if(muDtype==5)then
  MDmu(:,:)=0
  MDmu(1,2)=MDmu12
  do j=2,N_mu-1
  MDmu(j,j-1)=-w_point(j-1)/w_point(j)*MDmu(j-1,j)
  MDmu(j,j+1)=-MDmu(j,j-1)
  enddo
  j=N_mu
  MDmu(j,j-1)=-w_point(j-1)/w_point(j)*MDmu(j-1,j)
  endif

  if(myid==0 .and. verify_NL==1 )then

  veritemp=0
  veritempd=0
  do j=1,N_mu
    do jp=1,N_mu
    veritemp=veritemp+abs((w_point(j))*MDmu(j,jp)+(w_point(jp))*MDmu(jp,j))
    veritempd=veritempd+abs((w_point(j))*MDmu(j,jp)-(w_point(jp))*MDmu(jp,j))
    enddo
  enddo
  write(*,*)'Sum_j,jp[abs(Cj*Mjjp+Cjp*Mjpj)]/Sum_j,jp [abs(Cj*Mjjp-Cjp*Mjpj)]=',veritemp/veritempd
  endif

EndIF

!=================================================================================================================================
!GK
IF(GK_FK_CK==0)then
 !Calculate the grids for energy spectrum Energy(k)

  kmax=0
  kk_ngrid(:)=0
  do p=-pmax+1,pmax-1
    do n=-nmax+1,nmax-1
    kk=p**2+n**2
    if(kk/=0)then
    kk_ngrid(kk)=kk_ngrid(kk)+1
    endif
    enddo
  enddo

  do kk=1,(pmax-1)**2+(nmax-1)**2
    if(kk_ngrid(kk)/=0)then
     kmax=kmax+1
    endif
  enddo

  allocate(k(1:kmax))
  allocate(energy_k(1:kmax))

  i=0
  do kk=1,(pmax-1)**2+(nmax-1)**2
    if(kk_ngrid(kk)/=0)then
      i=i+1
      k(i)=sqrt(real(kk))*kxmax/pmax
    endif
  enddo

  open(stdout+1,file="k.txt",status='replace')
  write(stdout+1,*)kmax
  write(stdout+1,99)k
  close(stdout+1)

  F_k=(0.0,0.0)
  phi_k=(0.0,0.0)
 !set the value of the start point of nt in one run of this 3D code
  if(restart==0)then    !restart = 0 no restart; =1 restart without rollback; =2 restart and roll back
      nt_Odd=0
      nt_Even=0
      ntstart=0
      ntmaxstep=ntmax/backup_num    !ntmaxstep is the lenth of the ntmax pieces
      nt=0  !initial value of nt when run the code at the first time
    IF(myid==0)then
      open(101,file="nt_Odd.txt",status='replace')
      write(101,*)0
      close(101)
      open(101,file="nt_Even.txt",status='replace')
      write(101,*)0
      close(101)

     !We should set the initial value of F_k=F_k_int, when we run this 2Dcode at the first time
     Do i=0,N_mu
       if(i==0)then
       !F_kmu(:,:,i)=cmplx(F_k_int,F_k_int)
       call random_number(random)
       F_kmu(:,:,i)=F_k_int*exp((0,1.0)*2*PI*random(:,:))
       else
       F_kmu(:,:,i)=FM(i)*F_k_int*exp((0,1.0)*2*PI*random(:,:))
       !F_kmu(:,:,i)=FM(i)*cmplx(F_k_int,F_k_int)
       endif
     Enddo

       do p=-pmax+1,pmax-1
         do n=0,nmax-1
      F_kmu(-p,-n,:)=conjg(F_kmu(p,n,:))
         enddo
       enddo
       F_kmu(0,0,:)=(0,0)

    Endif


   IF(Add_GAM==1)then
    if(myid.ne.0)then
     call random_number(randomGAM)
      F_GAM(:)=FM(myid)*F_GAM_int*exp((0,1.0)*2*PI*randomGAM(:))
     ! F_GAM(:)=FM(myid)*cmplx(F_GAM_int,F_GAM_int)
       do p=0,pmax-1
      F_GAM(-p)=conjg(F_GAM(p))
       enddo
       F_GAM(0)=(0,0)
     endif
   Endif
     !!Add a initial pulse to the modes near k0
     !!the form of the pulse is: A0*exp(-(k-k0)^2/Alpha), where k=(kxmax/pmax)*sqrt(real(p**2+n**2))
       do p=-pmax+1,pmax-1
         do n=-nmax+1,nmax-1
         F_k(p,n)=F_k(p,n)+F_k(p,n)*A0*exp((-1.0*((kxmax/pmax)*sqrt(real(p**2+n**2))-k0)**2)/Alpha)
         enddo
       enddo

  else if(restart==1)then
    IF(myid==0)then
      open(101,file="nt_Odd.txt")
      read(101,*)nt_Odd
      close(101)
      open(101,file="nt_Even.txt")
      read(101,*)nt_Even
      close(101)
      if(ntmax .LT. max(nt_Odd,nt_Even))then
      write(*,*)'You should set restart=0 or set the ntmax much bigger for restart!!!'
      write(*,*)'New ntmax must larger than old ntmax of last run when you try to restart from backup point'
      stop
      endif
      nt=max(nt_Odd,nt_Even)
      if(nt==0)then
      write(*,*)"myid=",myid,"has exited the program!   Since the code has never made backup at last run."
      stop
      elseif(nt .GT. 0)then
      write(*,*)'continue the code from',(nt)*tstep,'s',' nt=',nt
      ntmaxstep=(ntmax-nt)/backup_num
      ntstart=nt
      endif
    Endif
  else if(restart==2)then
    IF(myid==0)then
      open(101,file="nt_Odd.txt")
      read(101,*)nt_Odd
      close(101)
      open(101,file="nt_Even.txt")
      read(101,*)nt_Even
      close(101)
      if(ntmax .LT. min(nt_Odd,nt_Even))then
      write(*,*)'You should set restart=0 or set the ntmax much bigger for restart!!!'
      write(*,*)'!!!new ntmax must larger than backup point of last run when you try to restart from backup point'
      stop
      endif
      nt=min(nt_Odd,nt_Even)
      if(nt==0)then
        write(*,*)"Exit the program!!!   Since the code has only made backup once, you should not roll back when you restart!"
      stop
      elseif(nt .GT. 0)then
      write(*,*)'Roll back!  We continue the code from',(nt)*tstep,'s' ,' nt=',nt
      ntmaxstep=(ntmax-nt)/backup_num
      ntstart=nt
      endif
    Endif
  else
      write(*,*)"restartGK=",restart
      write(*,*)"myid=",myid,'Has exited the program!    The value of restart must be 0 ,1 OR 2'
      stop
  endif

   call MPI_BCAST(ntmaxstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process
   call MPI_BCAST(ntstart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process

    If(restart==1 .or.restart==2 )then
     IF(myid==0)then
      !We should set the initial value of F_k from the restart point
      if(nt_Odd==ntstart)then
      open(1000,file='Fk_re_Odd')
      open(1001,file='Fk_im_Odd')
      elseif(nt_Even==ntstart)then
      open(1000,file='Fk_re_Even')
      open(1001,file='Fk_im_Even')
      endif
      read(1000,99)Fk_re_rt
      read(1001,99)Fk_im_rt
      close(1000)
      close(1001)

      F_kmu=cmplx(Fk_re_rt(:,:,:),Fk_im_rt(:,:,:))

      IF(Add_GAM==1)then
      if(nt_Odd==ntstart)then
      open(1004,file='FGAM_re_Odd')
      open(1005,file='FGAM_im_Odd')
      elseif(nt_Even==ntstart)then
      open(1004,file='FGAM_re_Even')
      open(1005,file='FGAM_im_Even')
      endif
      read(1004,99)F_GAM_re_rt
      read(1005,99)F_GAM_im_rt
      close(1004)
      close(1005)
      Endif
     Endif


     IF(Add_GAM==1)then
     call MPI_BCAST(F_GAM_re_rt,(2*pmax-1)*(N_mu),MPI_REAL,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process
     call MPI_BCAST(F_GAM_im_rt,(2*pmax-1)*(N_mu),MPI_REAL,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process
     if(myid.ne.0)then
     F_GAM=cmplx(F_GAM_re_rt(:,myid),F_GAM_im_rt(:,myid))
     endif
     Endif
  Endif

    IF(Add_GAM==1)then
    if(myid.ne.0)then
    do p=-pmax+1,pmax-1
      F_k(p,0)=F_k(p,0)+F_GAM(p)
    enddo
    endif
    Endif

    IF(myid==0)then

    do p=-pmax+1,pmax-1
      do n=0,nmax-1
        if(n/=0)then
        lambda_k=lambda_n
        else
        lambda_k=lambda_0
        end if
        if(n==0.and.p==0)then
        Integ_den(0,0)=0
        H_m_delta(p,n)=(1,0)  !used for calculate the denometor of the possion equation
        else
        Integ_den(p,n)=0
        do i=1,N_mu
          Integ_den(p,n) = Integ_den(p,n)+w_point(i)*FM(i)*(1-Jn(p,n,i,0)**2)/(G*G)
        enddo
        !calculate the denometor of the possion equation
        H_m_delta(p,n)=Integ_den(p,n)/ratio_TiTe+(1-CDW)*(lambda_k-(0,1.0)*delta_k(p,n))+&
                       CDW*(lambda_D**2)*(kx(p)**2+ky(n)**2)
        endif
       enddo
     enddo

     do p=-pmax+1,pmax-1
       do n=1,nmax-1
       H_m_delta(-p,-n)=H_m_delta(p,n)
       enddo
     enddo


    do p=-pmax+1,pmax-1
      do n=0,nmax-1
       if(n==0.and.p==0)then
       phi_k(0,0)=(0,0)!the real part of (0,0)mode do not change, the image part go to 0 because of conjugation!!
       else
        Integ_num(p,n)=(0,0)
        do i=1,N_mu
        Integ_num(p,n)=Integ_num(p,n)+w_point(i)*Jn(p,n,i,0)*F_kmu(p,n,i)
        enddo
        i=0
        Integ_num(p,n)=Integ_num(p,n)-CDW*F_kmu(p,n,i)

        phi_k(p,n)=Integ_num(p,n)/H_m_delta(p,n)    !!!***===
       endif
      enddo
    enddo

     do p=-pmax+1,pmax-1
       do n=1,nmax-1
       phi_k(-p,-n)=conjg(phi_k(p,n))
       enddo
     enddo

    Endif

  call MPI_BCAST(phi_k,(2*pmax-1)*(2*nmax-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process
  call MPI_BCAST(F_kmu(:,:,:),(2*pmax-1)*(2*nmax-1)*(N_mu+1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process


!----------------------------------------------------------------------------------------------
  IF(myid==0)then
    nmatr=N_mu  !the dimentioin of matrix
    ALLOCATE(Matrix(0:nmatr,0:nmatr))
    ALLOCATE(eigenO(0:nmatr))
    ALLOCATE(WORK(1:4*(nmatr+1)))
    Allocate(Omega_matr(-pmax+1:pmax-1,0:nmax-1,0:N_mu))
    Allocate(R_Matrix(-pmax+1:pmax-1,0:nmax-1,0:N_mu,0:N_mu),&
           M_Matrix(-pmax+1:pmax-1,0:nmax-1,0:N_mu,0:N_mu))
  ALLOCATE(IPIV(0:nmatr))
  ALLOCATE(Matrtemp(0:N_mu,0:N_mu))
  ALLOCATE(MatrtempREV(0:N_mu,0:N_mu))
  ALLOCATE(MatrtempOUT(0:N_mu,0:N_mu))

  !Calculate the linear Omega theoriticaly
  M_Matrix=(0,0)
  do p=-pmax+1,pmax-1
    do n=0,nmax-1
    If(p==0 .and. n==0)then
    Omega_matr(p,n,:)=(0,0)-(0,1.0)*nu_ZFcode(0)
    Else
      Matrix(:,:)=(0,0)
      do  jp=1,N_mu
        do  j=1,N_mu
        If(jp==j)then
        Matrix(j,jp)=Matrix(j,jp)+(0,1.0)*ky(n)*ratio_TiTe*Epsilon*mu_point(j)   !Only the variables comes from Phi_k must use the column subscript i, others use the row  subscript j

        Matrix(j,jp)=Matrix(j,jp)- mu_HK*(kx(p)**2+ky(n)**2)**2&
                      -nu_DWcode(p,n)&
                      -nu_ZFcode(n)&
                      !-mu_LK/((abs(ky(n))+0.03)*(kx(p)**2+ky(n)**2))
                      -mu_LK/(kx(p)**2+ky(n)**2)
        Endif
        Matrix(j,jp)=Matrix(j,jp)+ (0,1.0)*ky(n)*FM(j)*Jn(p,n,j,0)*(Epsilon*mu_point(j)-a_Ln-a_LTi*(mu_point(j)-1))&
                    *w_point(jp)*Jn(p,n,jp,0)/H_m_delta(p,n)
        enddo
      enddo
      jp=0
      do j=1,N_mu
      Matrix(j,jp)=Matrix(j,jp)+(0,1.0)*ky(n)*FM(j)*Jn(p,n,j,0)*(Epsilon*mu_point(j)-a_Ln-a_LTi*(mu_point(j)-1))&
                  *(-CDW)/H_m_delta(p,n)
      enddo

      !CDW electron dynamics in the following
      j=0
      jp=0
      Matrix(j,jp)=Matrix(j,jp)- ((0,1.0)*ky(n)*Epsilon + AlphaA)&
                            + (AlphaA*(1.0-(0,1.0)*delta_k(p,n)*CDWid) + (0,1.0)*ky(n)*(Epsilon-a_Ln))*(-CDW)/H_m_delta(p,n)

      Matrix(j,jp)=Matrix(j,jp)- mu_HK*(kx(p)**2+ky(n)**2)**2&
                      -nu_DWcode(p,n)&
                      -nu_ZFcode(n)&
                      !-mu_LK/((abs(ky(n))+0.03)*(kx(p)**2+ky(n)**2)) &
                      -mu_LK/(kx(p)**2+ky(n)**2)&
                      -(0,1.0)*kpar*uD


      do jp=1,N_mu
      Matrix(j,jp)=Matrix(j,jp)+ (AlphaA + (0,1.0)*ky(n)*(Epsilon-a_Ln))*w_point(jp)*Jn(p,n,jp,0)/H_m_delta(p,n)
      enddo

      M_Matrix(p,n,:,:)=Matrix(:,:)

 call ZGEEV('N', 'N', nmatr+1, Matrix, nmatr+1, eigenO, DUMMY, nmatr+1, DUMMY, nmatr+1, WORK, nmatr*2+2, WORK, info)
!    eigenO=1
    Omega_matr(p,n,:)=eigenO(:)/(0,-1)
    Endif
    enddo
  enddo

!Implicit time forward method for 4D Fourier kinetics
   IF(muDtype==2 .or. muDtype==-2)then
    if(Const_L==0.0)then
    M_Matrix(:,:,:,:)=(0,0)
    endif
   Else

    if(Const_L==0.0)then
    M_Matrix(:,:,:,:)=(0,0)
    endif

    do p=-pmax+1,pmax-1
      do n=0,nmax-1
         R_Matrix(p,n,:,:)=M_Matrix(p,n,:,:)*(-tstep/2)
         do l=0,N_mu
         R_Matrix(p,n,l,l)=R_Matrix(p,n,l,l)+1
         enddo
       enddo
     enddo

    do p=-pmax+1,pmax-1
      do n=0,nmax-1
         M_Matrix(p,n,:,:)=M_Matrix(p,n,:,:)*(tstep/2)
         do l=0,N_mu
         M_Matrix(p,n,l,l)=M_Matrix(p,n,l,l)+1
         enddo
       enddo
     enddo

    do p=-pmax+1,pmax-1
      do n=0,nmax-1
      Matrtemp(:,:)=R_Matrix(p,n,:,:)
      call zgetrf( nmatr+1, nmatr+1,Matrtemp , nmatr+1, IPIV, INFO)
      if(info/=0) write(0,*) 'Error occured in zgetrf!'
      call zgetri(nmatr+1,Matrtemp,nmatr+1, IPIV, WORK, nmatr+1, info ) !Be careful of the sixth parameter, should be set to n but not -1.
      if(info/=0) write(0,*) 'Error occured in zgetri!'
      R_Matrix(p,n,:,:)=Matrtemp(:,:)
      enddo
     enddo

   Endif


  If(restart==0)then
      open(stdout+29,file="Omega_matr.txt",status='replace')
   elseif(restart==1)then
      open(stdout+29,file="Omega_matr.txt",status='old',POSITION='APPEND')
   elseif(restart==2)then
      open(stdout+29,file="Omega_matr.txt",status='old')
   endif

    If(verify_L==1)then
      IF(nmax .GT. 3)then
      open(stdout+35,file="Omega_matr_3.txt",status='replace')
      write(stdout+35,*)'Omega_matr'
      write(stdout+35,*)"kx=","0    ; ","ky=",3*kymax/(nmax-1)
      do i=1,N_mu
        write(stdout+35,"('(',e14.7,',',e14.7,')')")&
        real(Omega_matr(0,3,i)),aimag(Omega_matr(0,3,i))
      ! write(stdout+35,"('(',i4,',',i4,')','=(',e14.7,',',e14.7,')')")&
      ! i,a,real(Omega_matrcy(pmax-1,nmax,l)),aimag(Omega_matrcy(pmax-1,nmax,l))
      enddo
      close(stdout+35)
      Endif
    Endif

     do  i=0,N_mu
      do n=0,nmax-1
       do p=-pmax+1,pmax-1
      write(stdout+29,"('(',e14.7,',',e14.7,')')")real(Omega_matr(p,n,i)),aimag(Omega_matr(p,n,i))
         enddo
       enddo
     enddo
     close(stdout+29)

  deallocate(Matrix)
  deallocate(eigenO)
  deallocate(WORK)
  deallocate(Omega_matr)
  deallocate(IPIV)
  deallocate(Matrtemp)
  deallocate(MatrtempREV)
  deallocate(MatrtempOUT)

  Endif


    IF(Add_GAM==1)then
    F_k(:,0)=F_k(:,0)-F_GAM(:)
    Endif

!----------------------------------------------------------------------------------------------
  IF(myid==0)then
  deallocate(Fk_re_rt)
  deallocate(Fk_im_rt)
  Endif
  IF(Add_GAM==1)then
  deallocate(F_GAM_re_rt)
  deallocate(F_GAM_im_rt)
  endif

  if(verify_L==1)then
  if(myid==0)then
  call CPU_TIME(Looptb)
  GKMt=0
  endif
  GKNLt=0
  endif

Endif


!=================================================================================================================================
!CKinFH
IF(GK_FK_CK==1)then

  do p=-pmax+1,pmax-1
    do n=-nmax+1,nmax-1
    Cphi_p(p,n)=(0,1.0)/ratio_TiTe*sqrt(kx(p)**2+ky(n)**2)/2*expibp(p,n)
    Cphi_m(p,n)=(0,1.0)/ratio_TiTe*sqrt(kx(p)**2+ky(n)**2)/2*expibm(p,n)
    enddo
  enddo

    do p=-pmax+1,pmax-1
      do n=-nmax+1,nmax-1
        if(n/=0)then
        lambda_k=lambda_n
        else
        lambda_k=lambda_0
        end if
        if(n==0.and.p==0)then
        lamb_m_delta(p,n)=(1,0)
        else
        lamb_m_delta(p,n)=(1-CDW)*(lambda_k-(0,1.0)*delta_k(p,n))+CDW*(lambda_D**2)*(kx(p)**2+ky(n)**2)
        endif
       enddo
     enddo

!---------------------------------------------------------------------------------
 !set the value of the start point of nt in one run of this 3D code
  if(restart==0)then    !restart = 0 no restart; =1 restart without rollback; =2 restart and roll back
      nt_Odd=0
      nt_Even=0
      ntstart=0
      ntmaxstep=ntmax/backup_num    !ntmaxstep is the lenth of the ntmax pieces
      nt=0  !initial value of nt when run the code at the first time
    IF(myid==0)then
      open(101,file="nt_Odd.txt",status='replace')
      write(101,*)0
      close(101)
      open(101,file="nt_Even.txt",status='replace')
      write(101,*)0
      close(101)
    !4D Fourier harmonic
    do l=0,(N_mu+1)*(2*N_FT-1)-1
      i=l/(2*N_FT-1)
      if(i==0)then
      call random_number(random)
      F_kfta(:,:,l)=F_k_int*exp((0,1.0)*2*PI*random(:,:))
      !F_kfta(:,:,l)=cmplx(F_k_int,F_k_int)
      else
      !F_kfta(:,:,l)=FM(i)*cmplx(F_k_int,F_k_int)
      F_kfta(:,:,l)=FM(i)*F_k_int*exp((0,1.0)*2*PI*random(:,:))
      endif
    enddo

    do p=-pmax+1,pmax-1
     do n=0,nmax-1       !Be careful you must start with n=0 (not n=1) to make sure the conjagation rule
      do l=0,(N_mu+1)*(2*N_FT-1)-1
     i=l/(2*N_FT-1)
     a=mod(l,2*N_FT-1)-(N_FT-1)
     j=i*(2*N_FT-1)-a+(N_FT-1)
     F_kfta(-p,-n,j)=conjg(F_kfta(p,n,l))
      enddo
     enddo
    enddo
    F_kfta(0,0,:)=(0,0)
   Endif

  else if(restart==1)then
    IF(myid==0)then
      open(101,file="nt_Odd.txt")
      read(101,*)nt_Odd
      close(101)
      open(101,file="nt_Even.txt")
      read(101,*)nt_Even
      close(101)
      if(ntmax .LT. max(nt_Odd,nt_Even))then
      write(*,*)'You should set restart=0 or set the ntmax much bigger for restart!!!'
      write(*,*)'New ntmax must larger than old ntmax of last run when you try to restart from backup point'
      stop
      endif
      nt=max(nt_Odd,nt_Even)
      if(nt==0)then
      write(*,*)"!! Exit the program!   Since the code has never made backup at last run."
      stop
      elseif(nt .GT. 0)then
      write(*,*)'continue the code from',(nt)*tstep,'s' ,' nt=',nt
      ntmaxstep=(ntmax-nt)/backup_num
      ntstart=nt
      endif
    Endif
  else if(restart==2)then
    IF(myid==0)then
      open(101,file="nt_Odd.txt")
      read(101,*)nt_Odd
      close(101)
      open(101,file="nt_Even.txt")
      read(101,*)nt_Even
      close(101)
      if(ntmax .LT. min(nt_Odd,nt_Even))then
      write(*,*)'You should set restart=0 or set the ntmax much bigger for restart!!!'
      write(*,*)'!!!new ntmax must larger than backup point of last run when you try to restart from backup point'
      stop
      endif
      nt=min(nt_Odd,nt_Even)
      if(nt==0)then
        write(*,*)"Exit the program!!!   Since the code has only made backup once, you should not roll back when you restart!"
      stop
      elseif(nt .GT. 0)then
      write(*,*)'Roll back!  We continue the code from',(nt)*tstep,'s' ,' nt=',nt
      ntmaxstep=(ntmax-nt)/backup_num
      ntstart=nt
      endif
    Endif
  else
    IF(myid==0)then
      write(*,*)"restartFK=",restart
      write(*,*)'!! Exit the program!    The value of restart must be 0 ,1 OR 2'
      stop
    Endif
  endif

   call MPI_BCAST(ntmaxstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process
   call MPI_BCAST(ntstart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process

  If(restart==1 .or.restart==2 )then
    IF(myid==0)then
  !We should set the initial value of F_kft from the restart point
  !4D Fourier harmonic
      if(nt_Odd==ntstart)then
      open(1008,file='Fkft_re_Odd')
      open(1009,file='Fkft_im_Odd')
      elseif(nt_Even==ntstart)then
      open(1008,file='Fkft_re_Even')
      open(1009,file='Fkft_im_Even')
      endif
      read(1008,*)F_kft_re_rt
      read(1009,*)F_kft_im_rt
      close(1008)
      close(1009)
    F_kfta(:,:,:)=cmplx(F_kft_re_rt(:,:,:),F_kft_im_rt(:,:,:))
    Endif
  Endif

  If(myid==0)then
     do n=0,nmax-1
        if(n/=0)then
        lambda_k=lambda_n
        else
        lambda_k=lambda_0
        endif
     do p=-pmax+1,pmax-1
        if(n==0.and.p==0)then
        phi_kft(0,0)=(0,0)  !4D Fourier harmonic
        else
        !4D Fourier harmonic
        phi_kft(p,n)=(0,0)
        do i=1,N_mu
        a=0
        l=i*(2*N_FT-1)+a+(N_FT-1)
        phi_kft(p,n)=phi_kft(p,n)+w_point(i)/Beta*F_kfta(p,n,l)
        enddo
        a=0
        i=0
        l=i*(2*N_FT-1)+a+(N_FT-1)
        phi_kft(p,n)=phi_kft(p,n)-CDW*F_kfta(p,n,l)

        phi_kft(p,n)=phi_kft(p,n)/lamb_m_delta(p,n)
        endif
     enddo
     enddo

     do n=1,nmax-1
     do p=-pmax+1,pmax-1
       phi_kft(-p,-n)=conjg(phi_kft(p,n))  !4D Fourier harmonic
     enddo
     enddo

  Endif

  call MPI_BCAST(phi_kft,(2*pmax-1)*(2*nmax-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process
  call MPI_BCAST(F_kfta,(2*pmax-1)*(2*nmax-1)*(N_mu+1)*(2*N_FT-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process

!-------------------------------------------------------------------------------------

 !IF(myid==0)then
  !4D Fourier harmonic
  nmatrft=(N_mu+1)*(2*N_FT-1)
  IF(myid==0)then
  ALLOCATE(Omega_matrft(-pmax+1:pmax-1,0:nmax-1,0:(N_mu+1)*(2*N_FT-1)-1))
  ALLOCATE(Omegakmft(0:kimax-1,0:(N_mu+1)*(2*N_FT-1)-1,0:(N_mu+1)*(2*N_FT-1)-1))
  Endif
  ALLOCATE(Matrixft(0:nmatrft-1,0:nmatrft-1))
  ALLOCATE(eigenOft(0:nmatrft-1))
  ALLOCATE(WORKft(0:2*nmatrft-1))
  ALLOCATE(IPIV(0:nmatrft-1))
  ALLOCATE(Omegaft(0:kimax-1,0:(N_mu+1)*(2*N_FT-1)-1))
  ALLOCATE(Matrtemp(0:(N_mu+1)*(2*N_FT-1)-1,0:(N_mu+1)*(2*N_FT-1)-1))
  ALLOCATE(MatrtempREV(0:(N_mu+1)*(2*N_FT-1)-1,0:(N_mu+1)*(2*N_FT-1)-1))
  ALLOCATE(MatrtempOUT(0:(N_mu+1)*(2*N_FT-1)-1,0:(N_mu+1)*(2*N_FT-1)-1))
  allocate(R_Matrixft(0:kimax-1,0:(N_mu+1)*(2*N_FT-1)-1,0:(N_mu+1)*(2*N_FT-1)-1),&
           M_Matrixft(0:kimax-1,0:(N_mu+1)*(2*N_FT-1)-1,0:(N_mu+1)*(2*N_FT-1)-1))

    M_Matrixft=(0,0)
    R_Matrixft=(0,0)
  !4D Fourier harmonic with the exp(-i*Alpha) Fourier transform
   do ki=0,kimax-1
   j=myid*kimax+ki
   p=mod(j,(2*pmax-1))-(pmax-1)
   n=j/(2*pmax-1)

    Matrixft(:,:)=(0,0)
    If(p==0 .and. n==0)then
    Omegaft(ki,:)=(0,0)-(0,1.0)*nu_ZFcode(0)
    Else if (n .LT. nmax)then
      do a=-N_FT+1,N_FT-1
      do i=1,N_mu
        j=i*(2*N_FT-1)+a+(N_FT-1)
        jp=j
        a_p=a+1
        a_m=a-1
        if(a==-N_FT+1)then      !cyclic condition of Fourier harmonic
        a_m=N_FT-1
        else if(a==N_FT-1)then  !cyclic condition of Fourier harmonic
        a_p=-N_FT+1
        endif
        Matrixft(j,jp)=Matrixft(j,jp)+(0,1.0)*(Omega_d(n,i)+a*Omega_star)
        j_p=i*(2*N_FT-1)+a_p+(N_FT-1)
        Matrixft(j,j_p)=Matrixft(j,j_p)-(0,1.0)*krho(p,n,i)*Omega_star/2.0&
                                                    *expibp(p,n)  !jp=(a+N_FT-1+1+(i-1)*(2*N_FT-1))
        j_m=i*(2*N_FT-1)+a_m+(N_FT-1)
        Matrixft(j,j_m)=Matrixft(j,j_m)-(0,1.0)*krho(p,n,i)*Omega_star/2.0&
                                                    *expibm(p,n)
        Matrixft(jp,jp)=Matrixft(jp,jp)-mu_HK*(kx(p)**2+ky(n)**2)**2&
                      -nu_DWcode(p,n)&
                      -nu_ZFcode(n)&
                      !-mu_LK/((abs(ky(n))+0.03)*(kx(p)**2+ky(n)**2))
                      -mu_LK/(kx(p)**2+ky(n)**2)

      enddo
      enddo

      do a=-N_FT+1,N_FT-1
      do i=1,N_mu                   !for row, which gives a and i
      j=i*(2*N_FT-1)+a+(N_FT-1)
        do ip=1,N_mu
        ap=0
        jp=ip*(2*N_FT-1)+ap+(N_FT-1)                         !for column, while ap=0
        if(a==0)then                !as delta_n=0 limit the row to a=0 ,row is decide by both j and a
        Matrixft(j,jp)=Matrixft(j,jp) +(0,1.0)*(Omega_d(n,i)/ratio_TiTe+a*Omega_star/(ratio_TiTe)&  !go through the column whose a=0 from i=1 to N_mu
                          -Omega_nT(n,i))*FM(i)*w_point(ip)/lamb_m_delta(p,n)
        else if(a==-1)then          !as delta_n=-1 limit the row to a=-1
        Matrixft(j,jp)=Matrixft(j,jp)-(0,1.0)*krho(p,n,i)*Omega_star/(ratio_TiTe*2)*FM(i)*w_point(ip)/lamb_m_delta(p,n)&
                         *expibp(p,n)
        else if(a==1)then           !as delta_n=+1 limit the row to a=+1
        Matrixft(j,jp)=Matrixft(j,jp)-(0,1.0)*krho(p,n,i)*Omega_star/(ratio_TiTe*2)*FM(i)*w_point(ip)/lamb_m_delta(p,n)&
                         *expibm(p,n)
        endif
        enddo

        ip=0                !for column, while ap=0
        ap=0
        jp=ip*(2*N_FT-1)+ap+(N_FT-1)
        if(a==0)then                !as delta_n=0 limit the row to a=0 ,row is decide by both j and a
        Matrixft(j,jp)=Matrixft(j,jp)+(0,1.0)*(Omega_d(n,i)/ratio_TiTe+a*Omega_star/(ratio_TiTe)&
                          -Omega_nT(n,i))*FM(i)*(-CDW)*Beta/lamb_m_delta(p,n)
        else if(a==-1)then          !as delta_n=-1 limit the row to a=-1
        Matrixft(j,jp)=Matrixft(j,jp)-(0,1.0)*krho(p,n,i)*Omega_star/(ratio_TiTe*2)*FM(i)&
                         *expibp(p,n)*(-CDW)*Beta/lamb_m_delta(p,n)
        else if(a==1)then           !as delta_n=+1 limit the row to a=+1
        Matrixft(j,jp)=Matrixft(j,jp)-(0,1.0)*krho(p,n,i)*Omega_star/(ratio_TiTe*2)*FM(i)&
                         *expibm(p,n)*(-CDW)*Beta/lamb_m_delta(p,n)
        endif

        do nc=-1,1
        !do nc=-N_FT+1,N_FT-1
        do ap=-N_FT+1,N_FT-1
        ip=i
        jp=ip*(2*N_FT-1)+ap+(N_FT-1)
        Matrixft(j,jp)=Matrixft(j,jp)+(expibm(p,n)**a)*gamIC_k(p,n)*(expibp(p,n)**ap)*(nc*Jn(p,n,i,nc-a)*Jn(p,n,i,nc-ap))

        Enddo
        Enddo

      enddo
      enddo

      !CDW electron dynamics in the following on dynamics in the following
      a=0
      i=0
      j=i*(2*N_FT-1)+a+(N_FT-1)
      ap=0
      ip=0
      jp=ip*(2*N_FT-1)+ap+(N_FT-1)
      Matrixft(j,jp)=Matrixft(j,jp)- ((0,1.0)*ky(n)*Epsilon + AlphaA)&
                            + (AlphaA*(1.0-(0,1.0)*delta_k(p,n)*CDWid) + (0,1.0)*ky(n)*(Epsilon-a_Ln))*(-CDW)/lamb_m_delta(p,n)

      Matrixft(j,jp)=Matrixft(j,jp)-mu_HK*(kx(p)**2+ky(n)**2)**2&
                      -nu_DWcode(p,n)&
                      -nu_ZFcode(n)&
                      !-mu_LK/((abs(ky(n))+0.03)*(kx(p)**2+ky(n)**2)) &
                      -mu_LK/(kx(p)**2+ky(n)**2)&
                      -(0,1.0)*kpar*uD

      Matrixft(j,j)=Matrixft(j,j)+signk(p,n)*gamE*(expibp(p,n)**a)*(Jn(p,n,i,1-a)-Jn(p,n,i,-1-a))

      ap=0
      do ip=1,N_mu
      jp=ip*(2*N_FT-1)+ap+(N_FT-1)
      Matrixft(j,jp)=Matrixft(j,jp)+ (AlphaA + (0,1.0)*ky(n)*(Epsilon-a_Ln))*w_point(ip)/Beta/lamb_m_delta(p,n)
      enddo


    M_Matrixft(ki,:,:)=Matrixft(:,:)

 call ZGEEV('N', 'N', nmatrft, Matrixft, nmatrft, eigenOft, DUMMY, nmatrft, DUMMY, nmatrft, WORKft, nmatrft*2, WORKft, info)
    Omegaft(ki,:)=eigenOft/(0,-1)
    Endif

   !Implicit time forward method for 4D Fourier kinetics
     IF(muDtype==2 .or. muDtype==-2)then
      if(Const_L==0.0)then
      M_Matrixft=(0,0)
      endif
     Else

      R_Matrixft(ki,:,:)=M_Matrixft(ki,:,:)*(-tstep/2)
      do l=0,(N_mu+1)*(2*N_FT-1)-1
      R_Matrixft(ki,l,l)=R_Matrixft(ki,l,l)+1
      enddo

      M_Matrixft(ki,:,:)=M_Matrixft(ki,:,:)*(tstep/2)
      do l=0,(N_mu+1)*(2*N_FT-1)-1
      i=l/(2*N_FT-1)
      a=mod(l,2*N_FT-1)-(N_FT-1)
      M_Matrixft(ki,l,l)=M_Matrixft(ki,l,l)+1
      if((i .eq. 0) .and. (a .ne. 0))then
      M_Matrixft(ki,l,l)=(0,0)
      endif
      enddo

      Matrtemp(:,:)=R_Matrixft(ki,:,:)
      call zgetrf( nmatrft, nmatrft,Matrtemp , nmatrft, IPIV, INFO)
      if(info/=0) write(0,*) 'Error occured in zgetrf!'
      call zgetri(nmatrft,Matrtemp,nmatrft, IPIV, WORKft, nmatrft, info ) !Be careful of the sixth parameter, should be set to n but not -1.
      if(info/=0) write(0,*) 'Error occured in zgetri!'
      R_Matrixft(ki,:,:)=Matrtemp(:,:)

      if(Const_L==0.0)then
      M_Matrixft(ki,:,:)=(0,0)
      R_Matrixft(ki,:,:)=(0,0)
      do l=0,(N_mu+1)*(2*N_FT-1)-1
      M_Matrixft(ki,l,l)=(1,0)
      R_Matrixft(ki,l,l)=(1,0)
      enddo
      do l=0,(N_mu+1)*(2*N_FT-1)-1
      i=l/(2*N_FT-1)
      a=mod(l,2*N_FT-1)-(N_FT-1)
      if((i .eq. 0) .and. (a .ne. 0))then
      M_Matrixft(ki,l,l)=(0,0)
      R_Matrixft(ki,l,l)=(0,0)
      endif
      enddo
      endif

   Endif
  enddo

!  Endif

  call MPI_GATHER(Omegaft(:,:),kimax*(N_mu+1)*(2*N_FT-1),MPI_COMPLEX,&
                  Omegakmft(:,:,:),kimax*(N_mu+1)*(2*N_FT-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)

    If(myid==0)then
      do ki=0,kimax-1
       do i=0,(N_mu+1)*(2*N_FT-1)-1
        j=i*kimax+ki
        p=mod(j,(2*pmax-1))-(pmax-1)
        n=j/(2*pmax-1)
        if (n .LT. nmax)then
        Omega_matrft(p,n,:) =Omegakmft(ki,:,i)
        endif
       enddo
      enddo

      open(stdout+40,file="Omega_matrft.txt",status='replace')
      do  i=0,(N_mu+1)*(2*N_FT-1)-1
       do n=0,nmax-1
        do p=-pmax+1,pmax-1
      write(stdout+40,"('(',e14.7,',',e14.7,')')")real(Omega_matrft(p,n,i)),aimag(Omega_matrft(p,n,i))
        enddo
       enddo
      enddo
      close(stdout+40)


      do  i=0,(N_mu+1)*(2*N_FT-1)-1
       do n=0,nmax-1
        do p=-pmax+1,pmax-1
        tempOMG=tempOMG+aimag(Omega_matrft(p,n,i)) /1000000.0
        enddo
       enddo
      enddo
      write(*,*)''
      write(*,*)"tempOMG=",tempOMG
      write(*,*)"sumOMG=",sum(aimag(Omega_matrft(:,:,:)))/1000000.0
      write(*,*)''
    Endif

  if(verify_L==2)then
  If(myid==0)then
    do n=-nmax+1,nmax-1
    do p=-pmax+1,pmax-1
    tempOMG=maxval(aimag(Omega_matrft(p,n,:)))
    gammax(p,n)=(1.0+tstep*tempOMG/2.0)/(1.0-tstep*tempOMG/2.0)
    enddo
    enddo
  Endif
  endif

  deallocate(Matrixft)
  deallocate(eigenOft)
  deallocate(WORKft)
  deallocate(IPIV)
  deallocate(Matrtemp)
  deallocate(Omegaft)

  IF(myid==0)then
  !4D Fourier harmonic
  deallocate(MatrtempREV)
  deallocate(MatrtempOUT)
  Endif

  IF(myid==0)then
  deallocate(Omega_matrft)
  deallocate(Omegakmft)
  deallocate(F_kft_re_rt)
  deallocate(F_kft_im_rt)
  Endif


  IF(myid==0)then
  IF(verify_NL==1)then
  !open(stdout+58,file="NLSft.txt",status='replace')
  Endif
  Endif

!  call MPI_Barrier(MPI_COMM_WORLD, ierr)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Calculate Omegacy in FK subrutine since idl plot Omegaft need Omegacy
 IF(restart==0)then
  do ki=0,kimax-1
   j=myid*kimax+ki
   p=mod(j,(2*pmax-1))-(pmax-1)
   n=j/(2*pmax-1)
   if (n .LT. nmax)then
   do a=-N_CY+1,N_CY-1
     !b=ratio_TiTe*(kx(p)**2+ky(n)**2)
     GamFunft(ki,a)=0
     do i=1,N_mu
     GamFunft(ki,a)=GamFunft(ki,a)+(Jn(p,n,i,a)**2)*FM(i)*w_point(i)
     enddo
   enddo

        if(n/=0)then
        lambda_k=lambda_n
        else
        lambda_k=lambda_0
        end if
        if(n==0.and.p==0)then
        G_m_deltaft(ki)=(1,0)
        else
        G_m_deltaft(ki)=(0,0)   !set the integral variables to zero before integral
        do a=-N_CY+1,N_CY-1
          G_m_deltaft(ki)=G_m_deltaft(ki)+GamFunft(ki,a)
        enddo
        G_m_deltaft(ki)=CHpola*(1-G_m_deltaft(ki))/ratio_TiTe+ &
                      (1-CDW)*(lambda_k-(0,1.0)*delta_k(p,n))+CDW*(lambda_D**2)*(kx(p)**2+ky(n)**2)
        endif

  endif
  enddo

 !cyclo-kinetic
  IF(myid==0)then
  ALLOCATE(Omega_matrcy(-pmax+1:pmax-1,0:nmax-1,0:(N_mu+1)*(2*N_CY-1)-1))
  ALLOCATE(Omegakmcy(0:kimax-1,0:(N_mu+1)*(2*N_CY-1)-1,0:(N_mu+1)*(2*N_FT-1)-1))
  Endif
  nmatrcy=(N_mu+1)*(2*N_CY-1)
 ! open(stdout+301,file="nmatrcy.txt",status='replace')
  ALLOCATE(Matrixcy(0:nmatrcy-1,0:nmatrcy-1))
  ALLOCATE(eigenOcy(0:nmatrcy-1))
  ALLOCATE(WORKcy(0:2*nmatrcy-1))
  ALLOCATE(IPIV(1:nmatrcy-1))
  ALLOCATE(Omegacy(0:kimax-1,0:(N_mu+1)*(2*N_CY-1)-1))

!cyclo-kinetic
   do ki=0,kimax-1
   j=myid*kimax+ki
   p=mod(j,(2*pmax-1))-(pmax-1)
   n=j/(2*pmax-1)
   If(p==0 .and. n==0)then
    Omegacy(ki,:)=(0,0)-(0,1.0)*nu_ZFcode(0)
  !  Omega_matrcy(0,0,:)=(0,0)-(0,1.0)*nu_ZFcode(0)
   Else if (n .LT. nmax)then
      Matrixcy=(0,0)
      do a=-N_CY+1,N_CY-1
      do i=1,N_mu
      j=i*(2*N_CY-1)+a+(N_CY-1)
        do ap=-N_CY+1,N_CY-1
        do ip=1,N_mu
        jp=ip*(2*N_CY-1)+ap+(N_CY-1)
        If(j==jp)then
        Matrixcy(j,jp)=Matrixcy(j,jp)&
                       +(0,1.0)*(Omega_d(n,i)+a*Omega_star)&
                       +(0,1.0)*(Omega_d(n,i)/ratio_TiTe+a*Omega_star/(ratio_TiTe)-Omega_nT(n,i))&
                       *FM(i)*Jn(p,n,i,a)*w_point(ip)*Jn(p,n,ip,ap)/G_m_deltaft(ki)

        Matrixcy(j,jp)=Matrixcy(j,jp)-mu_HK*(kx(p)**2+ky(n)**2)**2&
                          -nu_DWcode(p,n)&
                          -nu_ZFcode(n)&
                          -mu_LK/(kx(p)**2+ky(n)**2)

          If(a==1 .or.a==-1)then
          Matrixcy(j,jp)=Matrixcy(j,jp)+a*gamIC_k(p,n)
          Endif

          If(a .ne. 0)then
          Matrixcy(j,jp)=Matrixcy(j,jp)+  D_IC*(&
                           mu_HK*(kx(p)**2+ky(n)**2)**2&
                          +nu_DWcode(p,n)&
                          +nu_ZFcode(n)&
                          +mu_LK/(kx(p)**2+ky(n)**2) )          
          Endif
          
        Else
        Matrixcy(j,jp)=Matrixcy(j,jp)&
                       +(0,1.0)*(Omega_d(n,i)/ratio_TiTe+a*Omega_star/(ratio_TiTe)-Omega_nT(n,i))&
                       *FM(i)*Jn(p,n,i,a)*w_point(ip)*Jn(p,n,ip,ap)/G_m_deltaft(ki)
        Endif
        enddo
        enddo

        ip=0
        ap=0
        jp=ip*(2*N_CY-1)+ap+(N_CY-1)
        Matrixcy(j,jp)=Matrixcy(j,jp)&
                       +(0,1.0)*(Omega_d(n,i)/ratio_TiTe+a*Omega_star/(ratio_TiTe)-Omega_nT(n,i))&
                       *FM(i)*Jn(p,n,i,a)*(-CDW)/G_m_deltaft(ki)
      enddo
      enddo

      !CDW electron dynamics in the following
      a=0
      i=0
      j=i*(2*N_CY-1)+a+(N_CY-1)
      ap=0
      ip=0
      jp=ip*(2*N_CY-1)+ap+(N_CY-1)
      Matrixcy(j,jp)=Matrixcy(j,jp)- ((0,1.0)*ky(n)*Epsilon + AlphaA)&
                            + (AlphaA*(1.0-(0,1.0)*delta_k(p,n)*CDWid) + (0,1.0)*ky(n)*(Epsilon-a_Ln))*(-CDW)/G_m_deltaft(ki)

      Matrixcy(j,jp)=Matrixcy(j,jp)-mu_HK*(kx(p)**2+ky(n)**2)**2&
                          -nu_DWcode(p,n)&
                          -nu_ZFcode(n)&
                          !-mu_LK/((abs(ky(n))+0.03)*(kx(p)**2+ky(n)**2)) &
                          -mu_LK/(kx(p)**2+ky(n)**2)&
                          -(0,1.0)*kpar*uD


      do ap=-N_CY+1,N_CY-1
      do ip=1,N_mu
      jp=ip*(2*N_CY-1)+ap+(N_CY-1)
      Matrixcy(j,jp)=Matrixcy(j,jp)+ (AlphaA + (0,1.0)*ky(n)*(Epsilon-a_Ln))*w_point(ip)*Jn(p,n,ip,ap)/G_m_deltaft(ki)
      enddo
      enddo

 call ZGEEV('N', 'N', nmatrcy, Matrixcy, nmatrcy, eigenOcy, DUMMY, nmatrcy, DUMMY, nmatrcy, WORKcy, nmatrcy*2, WORKcy, info)
!    eigenOcy=1
    Omegacy(ki,:)=eigenOcy/(0,-1)
   ! write(stdout+35,*)"eigenOcy",eigenOcy
  Endif
 enddo
!------------------------------------------------------------------------
  call MPI_GATHER(Omegacy(:,:),kimax*(N_mu+1)*(2*N_CY-1),MPI_COMPLEX,&
                  Omegakmcy(:,:,:),kimax*(N_mu+1)*(2*N_CY-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)

    If(myid==0)then
      do ki=0,kimax-1
       do i=0,(N_mu+1)*(2*N_FT-1)-1
        j=i*kimax+ki
        p=mod(j,(2*pmax-1))-(pmax-1)
        n=j/(2*pmax-1)
        if (n .LT. nmax)then
        Omega_matrcy(p,n,:) =Omegakmcy(ki,:,i)
        endif
       enddo
      enddo

      open(stdout+30,file="Omega_matrcy.txt",status='replace')
      do  i=0,(N_mu+1)*(2*N_CY-1)-1
       do n=0,nmax-1
        do p=-pmax+1,pmax-1
      write(stdout+30,"('(',e14.7,',',e14.7,')')")real(Omega_matrcy(p,n,i)),aimag(Omega_matrcy(p,n,i))
        enddo
       enddo
      enddo
      close(stdout+30)
    endif
 !cyclo-kinetic
  deallocate(Matrixcy)
  deallocate(eigenOcy)
  deallocate(WORKcy)
  deallocate(IPIV)
  deallocate(Omegacy)
  IF(myid==0)then
  deallocate(Omega_matrcy)
  deallocate(Omegakmcy)
  endif

 Endif    !End the "IF(restart==0)then", also means end CK Matrix solver in FK
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  deallocate(GamFunft)
  deallocate(G_m_deltaft)

  if(verify_L==1)then
  NLt1=0
  Mt=0
  Mt1=0
  Mt2=0
  Mt3=0
  CtNL=0
  CtM=0
  NLcount=0
  if(myid==0)then
  call CPU_TIME(Looptb)
  DBt=0
  endif
  endif

Endif
!=================================================================================================================================
!CKinCH
IF(GK_FK_CK==2)then

 if(myid==0)then
 Delta_nnp=(0,0)
 do a=-N_CY+1,N_CY-1
  do n=0,nmax-1
  do p=-pmax+1,pmax-1
    do n1=-nmax+1,nmax-1
    do p1=-pmax+1,pmax-1     !sum_k1
      do ap=-N_CY+1,N_CY-1
      p2=p-p1
      n2=n-n1
      if((p2.GE.-pmax+1).AND.(p2.LE.pmax-1).AND.(n2.GE.-nmax+1).AND.(n2.LE.nmax-1))then
      Delta_nnp(ap,p1,n1,p,n,a)=(expibm(p1,n1)**a)*&
                                  (expibp(p,n)**a)*&
                                  (expibp(p1,n1)**ap)*&
                                  (expibm(p2,n2)**ap)
      endif
      enddo
      enddo
    enddo
    enddo
  enddo
  enddo

   GamFun=0
   do n=0,nmax-1
   do p=-pmax+1,pmax-1
     do a=-N_CY+1,N_CY-1
     do i=1,N_mu
       GamFun(a,p,n)=GamFun(a,p,n)+(Jn(p,n,i,a)**2)*FM(i)*w_point(i)
     enddo
     enddo
   enddo
   enddo

    do n=0,nmax-1
    do p=-pmax+1,pmax-1
        if(n/=0)then
        lambda_k=lambda_n
        else
        lambda_k=lambda_0
        end if
        if(n==0.and.p==0)then
        G_m_delta(p,n)=(1,0)
        else
        G_m_delta(p,n)=(0,0)   !set the integral variables to zero before integral
        do a=-N_CY+1,N_CY-1
          G_m_delta(p,n)=G_m_delta(p,n)+GamFun(a,p,n)
        enddo
        G_m_delta(p,n)=CHpola*(1-G_m_delta(p,n))/ratio_TiTe+ &
                      (1-CDW)*(lambda_k-(0,1.0)*delta_k(p,n))+CDW*(lambda_D**2)*(kx(p)**2+ky(n)**2)
        endif
     enddo
     enddo
 endif

 call MPI_BCAST(Delta_nnp,nmax*(2*nmax-1)*((2*pmax-1)*(2*N_CY-1))**2,MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(G_m_delta,nmax*(2*pmax-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)


i=mod(myid,N_mu+1)
a=myid/(N_mu+1)-(N_CY-1)
IF(i .gt. 0)then
 NL1=(0,0)
 NL2=(0,0)
 If(Const_NL .ne. 0) then
  do n=0,nmax-1
  do p=-pmax+1,pmax-1
    do n1=-nmax+1,nmax-1
    do p1=-pmax+1,pmax-1     !sum_k1
      do ap=-N_CY+1,N_CY-1
      p2=p-p1
      n2=n-n1
      !-----------------------------------------------------------------------
        NL1(ap,p1,n1,p,n)=&
            +2.0*(k1x(p1)*k2y(n2)-k2x(p2)*k1y(n1))*Jn(p1,n1,i,a-ap)*Delta_nnp(ap,p1,n1,p,n,a)&

            +sqrt(kx(p1)**2+ky(n1)**2)/rho(i)*(Jn(p1,n1,i,a-ap+1)-Jn(p1,n1,i,a-ap-1))/(2.0*(0,1.0))&
                      *ap*Delta_nnp(ap,p1,n1,p,n,a)&

            +sqrt(kx(p1)**2+ky(n1)**2)/rho(i)*(Jn(p1,n1,i,a-ap+1)-Jn(p1,n1,i,a-ap-1))/(2.0*(0,1.0))&
                      *a*Delta_nnp(ap,p1,n1,p,n,a)

!            +sqrt(kx(p1)**2+ky(n1)**2)/rho(i)*((ap-1)*Jn(p1,n1,i,a-ap+1)-(ap+1)*Jn(p1,n1,i,a-ap-1))/(2.0*(0,1.0))&
!                        *Delta_nnp(ap,p1,n1,p,n,a)&
!
!            +(kx(p1)**2+ky(n1)**2)*(Jn(p1,n1,i,a-ap+2)-Jn(p1,n1,i,a-ap-2))/(4.0*(0,1.0))&
!                        *Delta_nnp(ap,p1,n1,p,n,a)

          do ip=1,N_mu
          NL2(ip,ap,p1,n1,p,n)=&
              +sqrt(kx(p1)**2+ky(n1)**2)/rho(i)*(Jn(p1,n1,i,a-ap+1)+Jn(p1,n1,i,a-ap-1))/(2.0*(0,1.0))&
                 *2.0*mu_point(i)*MDmu(i,ip)*Delta_nnp(ap,p1,n1,p,n,a)&

              +sqrt(kx(p1)**2+ky(n1)**2)*rho(ip)/ratio_TiTe*Delta_nnp(ap,p1,n1,p,n,a)/(2.0*(0,1.0))&
                 *MDmu(i,ip)*(Jn(p1,n1,ip,a-ap+1)+Jn(p1,n1,ip,a-ap-1))
          enddo
      !-----------------------------------------------------------------------
      enddo
    enddo
    enddo
  enddo
  enddo
 Endif
Endif

Endif
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!CKinCH deeply paralized
IF(GK_FK_CK== -2)then

 if(myid==0)then
   GamFun=0
   do n=0,nmax-1
   do p=-pmax+1,pmax-1
     do a=-N_CY+1,N_CY-1
     do i=1,N_mu
       GamFun(a,p,n)=GamFun(a,p,n)+(Jn(p,n,i,a)**2)*FM(i)*w_point(i)
     enddo
     enddo
   enddo
   enddo

    do n=0,nmax-1
    do p=-pmax+1,pmax-1
        if(n/=0)then
        lambda_k=lambda_n
        else
        lambda_k=lambda_0
        end if
        if(n==0.and.p==0)then
        G_m_delta(p,n)=(1,0)
        else
        G_m_delta(p,n)=(0,0)   !set the integral variables to zero before integral
        do a=-N_CY+1,N_CY-1
          G_m_delta(p,n)=G_m_delta(p,n)+GamFun(a,p,n)
        enddo
        G_m_delta(p,n)=CHpola*(1-G_m_delta(p,n))/ratio_TiTe+ &
                      (1-CDW)*(lambda_k-(0,1.0)*delta_k(p,n))+CDW*(lambda_D**2)*(kx(p)**2+ky(n)**2)
        endif
     enddo
     enddo
 endif

 call MPI_BCAST(G_m_delta,nmax*(2*pmax-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)


  i=mod(myid,N_mu+1)
  a=mod((myid-i)/(N_mu+1),(2*N_CY-1))-(N_CY-1)
  p=myid/((N_mu+1)*(2*N_CY-1))-(pmax-1)
  Delta_nnpdp=(0,0)
  do n=0,nmax-1
    do n1=-nmax+1,nmax-1
    do p1=-pmax+1,pmax-1     !sum_k1
      do ap=-N_CY+1,N_CY-1
      p2=p-p1
      n2=n-n1
      if((p2.GE.-pmax+1).AND.(p2.LE.pmax-1).AND.(n2.GE.-nmax+1).AND.(n2.LE.nmax-1))then
      Delta_nnpdp(ap,p1,n1,n)=(expibm(p1,n1)**a)*&
                                  (expibp(p,n)**a)*&
                                  (expibp(p1,n1)**ap)*&
                                  (expibm(p2,n2)**ap)
      endif
      enddo
    enddo
    enddo
  enddo

 IF(i .gt. 0)then
 NL1dp=(0,0)
 NL2dp=(0,0)
 If(Const_NL .ne. 0) then
  do n=0,nmax-1
    do n1=-nmax+1,nmax-1
    do p1=-pmax+1,pmax-1     !sum_k1
      do ap=-N_CY+1,N_CY-1
      p2=p-p1
      n2=n-n1
      !-----------------------------------------------------------------------
        NL1dp(ap,p1,n1,n)=&
            +2.0*(k1x(p1)*k2y(n2)-k2x(p2)*k1y(n1))*Jn(p1,n1,i,a-ap)*Delta_nnpdp(ap,p1,n1,n)&

            +sqrt(kx(p1)**2+ky(n1)**2)/rho(i)*(Jn(p1,n1,i,a-ap+1)-Jn(p1,n1,i,a-ap-1))/(2.0*(0,1.0))&
                      *ap*Delta_nnpdp(ap,p1,n1,n)&

            +sqrt(kx(p1)**2+ky(n1)**2)/rho(i)*(Jn(p1,n1,i,a-ap+1)-Jn(p1,n1,i,a-ap-1))/(2.0*(0,1.0))&
                      *a*Delta_nnpdp(ap,p1,n1,n)

!            +sqrt(kx(p1)**2+ky(n1)**2)/rho(i)*((ap-1)*Jn(p1,n1,i,a-ap+1)-(ap+1)*Jn(p1,n1,i,a-ap-1))/(2.0*(0,1.0))&
!                        *Delta_nnpdp(ap,p1,n1,n)&
!
!            +(kx(p1)**2+ky(n1)**2)*(Jn(p1,n1,i,a-ap+2)-Jn(p1,n1,i,a-ap-2))/(4.0*(0,1.0))&
!                        *Delta_nnpdp(ap,p1,n1,n)

          do ip=1,N_mu
          NL2dp(ip,ap,p1,n1,n)=&
              +sqrt(kx(p1)**2+ky(n1)**2)/rho(i)*(Jn(p1,n1,i,a-ap+1)+Jn(p1,n1,i,a-ap-1))/(2.0*(0,1.0))&
                 *2.0*mu_point(i)*MDmu(i,ip)*Delta_nnpdp(ap,p1,n1,n)&

              +sqrt(kx(p1)**2+ky(n1)**2)*rho(ip)/ratio_TiTe*Delta_nnpdp(ap,p1,n1,n)/(2.0*(0,1.0))&
                 *MDmu(i,ip)*(Jn(p1,n1,ip,a-ap+1)+Jn(p1,n1,ip,a-ap-1))
          enddo
      !------------------------------------------------------------------------
      enddo
    enddo
    enddo
  enddo
 Endif
Endif

Endif
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!CKinCH
IF(GK_FK_CK==2 .or. GK_FK_CK== -2)then

!---------------------------------------------------------------------------------
 !set the value of the start point of nt in one run of this 3D code
  if(restart==0)then    !restart = 0 no restart; =1 restart without rollback; =2 restart and roll back
      nt_Odd=0
      nt_Even=0
      ntstart=0
      ntmaxstep=ntmax/backup_num    !ntmaxstep is the lenth of the ntmax pieces
      nt=0  !initial value of nt when run the code at the first time
    IF(myid==0)then
      open(101,file="nt_Odd.txt",status='replace')
      write(101,*)0
      close(101)
      open(101,file="nt_Even.txt",status='replace')
      write(101,*)0
      close(101)
!--------------------------------
    !call random_number(random)
    !random=0.6
    !do l=0,(N_mu+1)*(2*N_CY-1)-1
    !  i=mod(l,N_mu+1)
    !  if(i==0)then
    !  F_kcya(l,:,:)=F_k_int*exp((0,1.0)*2*PI*random(:,:))
    !  else
    !  F_kcya(l,:,:)=FM(i)*F_k_int*exp((0,1.0)*2*PI*random(:,:))
    !  endif
    !enddo
!--------------------------------
    do l=0,(N_mu+1)*(2*N_CY-1)-1
      i=mod(l,N_mu+1)
      if(i==0)then
      call random_number(random)
      F_kcya(l,:,:)=F_k_int*exp((0,1.0)*2*PI*random(:,:))
     ! F_kcya(l,:,:)=cmplx(F_k_int,F_k_int)
      else
      F_kcya(l,:,:)=FM(i)*F_k_int*exp((0,1.0)*2*PI*random(:,:))
      !F_kcya(l,:,:)=FM(i)*cmplx(F_k_int,F_k_int)
      endif
    enddo


    do n=0,nmax-1       !Be careful you must start with n=0 (not n=1) to make sure the conjagation rule
    do p=-pmax+1,pmax-1
      do l=0,(N_mu+1)*(2*N_CY-1)-1
      i=mod(l,N_mu+1)
      a=l/(N_mu+1)-(N_CY-1)
      j=(-a+N_CY-1)*(N_mu+1)+i
      F_kcya(j,-p,-n)=conjg(F_kcya(l,p,n))*((-1)**a)
      enddo
    enddo
    enddo
    F_kcya(:,0,0)=(0,0)

    do a=-N_CY+1,N_CY-1
    i=0
    j=(a+N_CY-1)*(N_mu+1)+i
    if(a .ne. 0)then
    F_kcya(j,:,:)=(0,0)
    endif
    enddo

   Endif

  else if(restart==1)then
    IF(myid==0)then
      open(101,file="nt_Odd.txt")
      read(101,*)nt_Odd
      close(101)
      open(101,file="nt_Even.txt")
      read(101,*)nt_Even
      close(101)
      if(ntmax .LT. max(nt_Odd,nt_Even))then
      write(*,*)'You should set restart=0 or set the ntmax much bigger for restart!!!'
      write(*,*)'New ntmax must larger than old ntmax of last run when you try to restart from backup point'
      stop
      endif
      nt=max(nt_Odd,nt_Even)
      if(nt==0)then
      write(*,*)"!! Exit the program!   Since the code has never made backup at last run."
      stop
      elseif(nt .GT. 0)then
      write(*,*)'continue the code from',(nt)*tstep,'s' ,' nt=',nt
      ntmaxstep=(ntmax-nt)/backup_num
      ntstart=nt
      endif
    Endif
  else if(restart==2)then
    IF(myid==0)then
      open(101,file="nt_Odd.txt")
      read(101,*)nt_Odd
      close(101)
      open(101,file="nt_Even.txt")
      read(101,*)nt_Even
      close(101)
      if(ntmax .LT. min(nt_Odd,nt_Even))then
      write(*,*)'You should set restart=0 or set the ntmax much bigger for restart!!!'
      write(*,*)'!!!new ntmax must larger than backup point of last run when you try to restart from backup point'
      stop
      endif
      nt=min(nt_Odd,nt_Even)
      if(nt==0)then
        write(*,*)"Exit the program!!!   Since the code has only made backup once, you should not roll back when you restart!"
      stop
      elseif(nt .GT. 0)then
      write(*,*)'Roll back!  We continue the code from',(nt)*tstep,'s' ,' nt=',nt
      ntmaxstep=(ntmax-nt)/backup_num
      ntstart=nt
      endif
    Endif
  else
    IF(myid==0)then
      write(*,*)"restartFK=",restart
      write(*,*)'!! Exit the program!    The value of restart must be 0 ,1 OR 2'
      stop
    Endif
  endif

   call MPI_BCAST(ntmaxstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process
   call MPI_BCAST(ntstart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process

  If(restart==1 .or.restart==2 )then
    IF(myid==0)then
  !We should set the initial value of F_kft from the restart point
  !4D Fourier harmonic
      if(nt_Odd==ntstart)then
      open(1006,file='Fkcy_re_Odd')
      open(1007,file='Fkcy_im_Odd')
      elseif(nt_Even==ntstart)then
      open(1006,file='Fkcy_re_Even')
      open(1007,file='Fkcy_im_Even')
      endif
      read(1006,99)F_kcy_re_rt!F_kcy_real_restart
      read(1007,99)F_kcy_im_rt!F_kcy_imag_restart
      close(1006)
      close(1007)
      F_kcya(:,:,:)=cmplx(F_kcy_re_rt(:,:,:),F_kcy_im_rt(:,:,:))
    Endif

  Endif

  If(myid==0)then
    do n=0,nmax-1
        if(n/=0)then
        lambda_k=lambda_n
        else
        lambda_k=lambda_0
        endif
    do p=-pmax+1,pmax-1
        if(n==0.and.p==0)then
        phi_kcy(0,0)=(0,0)  !cyclo-kinetic
        else
        !cyclo-kinetic
        phi_kcy(p,n)=(0,0)
        do a=-N_CY+1,N_CY-1
          do i=1,N_mu
          l=(a+N_CY-1)*(N_mu+1)+i
          phi_kcy(p,n)=phi_kcy(p,n)+w_point(i)*Jn(p,n,i,a)*F_kcya(l,p,n)
          enddo
        enddo
        i=0
        a=0
        l=(a+N_CY-1)*(N_mu+1)+i
        phi_kcy(p,n)=phi_kcy(p,n)-CDW*F_kcya(l,p,n)

        phi_kcy(p,n)=phi_kcy(p,n)/G_m_delta(p,n)
        endif
     enddo
     enddo

     do n=1,nmax-1
     do p=-pmax+1,pmax-1
       phi_kcy(-p,-n)=conjg(phi_kcy(p,n))
     enddo
     enddo

  Endif

  call MPI_BCAST(phi_kcy,(2*pmax-1)*(2*nmax-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process
  call MPI_BCAST(F_kcya,(2*pmax-1)*(2*nmax-1)*(N_mu+1)*(2*N_CY-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process

  nmatrcy=(N_mu+1)*(2*N_CY-1)

  IF(GK_FK_CK==2 .and. myid==0)then
  ALLOCATE(Omegakmcy(0:(N_mu+1)*(2*N_CY-1)-1,0:kimax-1,0:(N_mu+1)*(2*N_CY-1)-1))
  Endif
  IF(GK_FK_CK== -2 .and. myid==0)then
  ALLOCATE(Omegakmcy(0:(N_mu+1)*(2*N_CY-1)-1,0:kimax-1,0:(2*pmax-1)*(N_mu+1)*(2*N_CY-1)-1))
  Endif

  IF(myid==0)then
  ALLOCATE(Omega_matrcy(0:(N_mu+1)*(2*N_CY-1)-1,-pmax+1:pmax-1,0:nmax-1))
  Endif
  ALLOCATE(Matrixcy(0:nmatrcy-1,0:nmatrcy-1))
  ALLOCATE(eigenOcy(0:nmatrcy-1))
  ALLOCATE(WORKcy(0:2*nmatrcy-1))
  ALLOCATE(IPIV(1:nmatrcy-1))
  ALLOCATE(Omegacy(0:(N_mu+1)*(2*N_CY-1)-1,0:kimax-1))
  ALLOCATE(Matrtemp(0:(N_mu+1)*(2*N_CY-1)-1,0:(N_mu+1)*(2*N_CY-1)-1))
  allocate(R_Matrixcy(0:(N_mu+1)*(2*N_CY-1)-1,0:(N_mu+1)*(2*N_CY-1)-1,0:kimax-1),&
           M_Matrixcy(0:(N_mu+1)*(2*N_CY-1)-1,0:(N_mu+1)*(2*N_CY-1)-1,0:kimax-1))

    M_Matrixcy=(0,0)
    R_Matrixcy=(0,0)
!cyclo-kinetic
  do ki=0,kimax-1
   j=(myid-1)*kimax+ki
    If(j .GE. 0)then
    p=mod(j,(2*pmax-1))-(pmax-1)
    n=j/(2*pmax-1)
    Endif

    Matrixcy=(0,0)
    If(p==0 .and. n==0)then
    Omegacy(:,ki)=(0,0)-(0,1.0)*nu_ZFcode(0)
    Else if (n .LT. nmax)then
      do a=-N_CY+1,N_CY-1
      do i=1,N_mu
        j=(a+N_CY-1)*(N_mu+1)+i
        do ap=-N_CY+1,N_CY-1
        do ip=1,N_mu
        jp=(ap+N_CY-1)*(N_mu+1)+ip
        If(j==jp)then
        Matrixcy(jp,j)=Matrixcy(jp,j)&
                       +(0,1.0)*(Omega_d(n,i)+a*Omega_star)&
                       +(0,1.0)*(Omega_d(n,i)/ratio_TiTe+a*Omega_star/(ratio_TiTe)-Omega_nT(n,i))&
                       *FM(i)*Jn(p,n,i,a)*w_point(ip)*Jn(p,n,ip,ap)/G_m_delta(p,n)

        Matrixcy(jp,j)=Matrixcy(jp,j)-mu_HK*(kx(p)**2+ky(n)**2)**2&
                          -nu_DWcode(p,n)&
                          -nu_ZFcode(n)&
                          -mu_LK/(kx(p)**2+ky(n)**2)

          If(a==1 .or. a==-1)then
          Matrixcy(j,jp)=Matrixcy(j,jp)+a*gamIC_k(p,n)
          Endif
          
          If(a .ne. 0)then
          Matrixcy(j,jp)=Matrixcy(j,jp)+  D_IC*(&
                           mu_HK*(kx(p)**2+ky(n)**2)**2&
                          +nu_DWcode(p,n)&
                          +nu_ZFcode(n)&
                          +mu_LK/(kx(p)**2+ky(n)**2) )          
          Endif
          
        Else
        Matrixcy(jp,j)=Matrixcy(jp,j)&
                       +(0,1.0)*(Omega_d(n,i)/ratio_TiTe+a*Omega_star/(ratio_TiTe)-Omega_nT(n,i))&
                       *FM(i)*Jn(p,n,i,a)*w_point(ip)*Jn(p,n,ip,ap)/G_m_delta(p,n)
        Endif
        enddo
        enddo

        ip=0
        ap=0
        jp=(ap+N_CY-1)*(N_mu+1)+ip
        Matrixcy(jp,j)=Matrixcy(jp,j)&
                       +(0,1.0)*(Omega_d(n,i)/ratio_TiTe+a*Omega_star/(ratio_TiTe)-Omega_nT(n,i))&
                       *FM(i)*Jn(p,n,i,a)*(-CDW)/G_m_delta(p,n)
      enddo
      enddo

      !CDW electron dynamics in the following
      a=0
      i=0
      j=(a+N_CY-1)*(N_mu+1)+i
      ap=0
      ip=0
      jp=(ap+N_CY-1)*(N_mu+1)+ip
      Matrixcy(jp,j)=Matrixcy(jp,j)- ((0,1.0)*ky(n)*Epsilon + AlphaA)&
                    + (AlphaA*(1.0-(0,1.0)*delta_k(p,n)*CDWid) + (0,1.0)*ky(n)*(Epsilon-a_Ln))*(-CDW)/G_m_delta(p,n)

      Matrixcy(jp,j)=Matrixcy(jp,j)-mu_HK*(kx(p)**2+ky(n)**2)**2&
                          -nu_DWcode(p,n)&
                          -nu_ZFcode(n)&
                          !-mu_LK/((abs(ky(n))+0.03)*(kx(p)**2+ky(n)**2)) &
                          -mu_LK/(kx(p)**2+ky(n)**2)&
                          -(0,1.0)*kpar*uD

      a=0
      i=0
      j=(a+N_CY-1)*(N_mu+1)+i
      do ap=-N_CY+1,N_CY-1
      do ip=1,N_mu
      jp=(ap+N_CY-1)*(N_mu+1)+ip
      Matrixcy(jp,j)=Matrixcy(jp,j)+ (AlphaA + (0,1.0)*ky(n)*(Epsilon-a_Ln))*w_point(ip)*Jn(p,n,ip,ap)/G_m_delta(p,n)
      enddo
      enddo

      M_Matrixcy(:,:,ki)=Matrixcy(:,:)

 call ZGEEV('N', 'N', nmatrcy, Matrixcy, nmatrcy, eigenOcy, DUMMY, nmatrcy, DUMMY, nmatrcy, WORKcy, nmatrcy*2, WORKcy, info)
  Omegacy(:,ki)=eigenOcy/(0,-1)
  Endif

!------------------------------------------------------------------------
!Implicit time forward method for Cyclo-kinetic
   IF(muDtype==2 .or. muDtype==-2)then
      if(Const_L==0.0)then
      M_Matrixcy=(0,0)
      endif
   Else

      R_Matrixcy(:,:,ki)=M_Matrixcy(:,:,ki)*(-tstep/2)
      do l=0,(N_mu+1)*(2*N_CY-1)-1
      R_Matrixcy(l,l,ki)=R_Matrixcy(l,l,ki)+1
      enddo

      M_Matrixcy(:,:,ki)=M_Matrixcy(:,:,ki)*(tstep/2)
      do l=0,(N_mu+1)*(2*N_CY-1)-1
      a=l/(N_mu+1)-(N_CY-1)
      i=mod(l,N_mu+1)
      M_Matrixcy(l,l,ki)=M_Matrixcy(l,l,ki)+1
      if((i .eq. 0) .and. (a .ne. 0))then
      M_Matrixcy(l,l,ki)=(0,0)
      endif
      enddo

      Matrtemp(:,:)=R_Matrixcy(:,:,ki)
      call zgetrf( nmatrcy, nmatrcy,Matrtemp , nmatrcy, IPIV, INFO)
      if(info/=0) write(0,*) 'Error occured in zgetrf!'
      call zgetri(nmatrcy,Matrtemp,nmatrcy, IPIV, WORKcy, nmatrcy, info ) !Be careful of the sixth parameter, should be set to n but not -1.
      if(info/=0) write(0,*) 'Error occured in zgetri!'
      R_Matrixcy(:,:,ki)=Matrtemp(:,:)

      if(Const_L==0.0)then
      M_Matrixcy(:,:,ki)=(0,0)
      R_Matrixcy(:,:,ki)=(0,0)
      do l=0,(N_mu+1)*(2*N_CY-1)-1
      M_Matrixcy(l,l,ki)=(1,0)
      R_Matrixcy(l,l,ki)=(1,0)
      enddo
      do l=0,(N_mu+1)*(2*N_CY-1)-1
      a=l/(N_mu+1)-(N_CY-1)
      i=mod(l,N_mu+1)
      if((i .eq. 0) .and. (a .ne. 0))then
      M_Matrixcy(l,l,ki)=(0,0)
      R_Matrixcy(l,l,ki)=(0,0)
      endif
      enddo
      endif


   Endif
 enddo

  call MPI_GATHER(Omegacy(:,:),kimax*(N_mu+1)*(2*N_CY-1),MPI_COMPLEX,&
                  Omegakmcy(:,:,:),kimax*(N_mu+1)*(2*N_CY-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)

   If(myid==0)then

     IF(GK_FK_CK==2)then
     imax=(N_mu+1)*(2*N_CY-1)-1
     ENDIF
     IF(GK_FK_CK== -2)then
     imax=(2*pmax-1)*(N_mu+1)*(2*N_CY-1)-1
     ENDIF
      do i=0,imax
      do ki=0,kimax-1
      j=(i-1)*kimax+ki
      If(j .GE. 0)then
      p=mod(j,(2*pmax-1))-(pmax-1)
      n=j/(2*pmax-1)
        if (n .LT. nmax)then
        Omega_matrcy(:,p,n) =Omegakmcy(:,ki,i)
        endif
      Endif
      enddo
      enddo


      open(stdout+40,file="Omega_matrcy.txt",status='replace')

      do  i=0,(N_mu+1)*(2*N_CY-1)-1
      do n=0,nmax-1
      do p=-pmax+1,pmax-1
      write(stdout+40,"('(',e14.7,',',e14.7,')')")real(Omega_matrcy(i,p,n)),aimag(Omega_matrcy(i,p,n))
      enddo
      enddo
      enddo
      close(stdout+40)
   Endif

!------------------------------------------------------------------------
 IF(GK_FK_CK==2)then
  deallocate(Matrixcy)
  deallocate(eigenOcy)
  deallocate(WORKcy)
  deallocate(IPIV)
  deallocate(Matrtemp)
  deallocate(Omegacy)
  deallocate(Delta_nnp)
  IF(myid==0)then
  deallocate(Omega_matrcy)
  deallocate(Omegakmcy)
  deallocate(F_kcy_re_rt)
  deallocate(F_kcy_im_rt)
  deallocate(GamFun)
  Endif
 ENDIF


 IF(GK_FK_CK== -2)then
  deallocate(Matrixcy)
  deallocate(eigenOcy)
  deallocate(WORKcy)
  deallocate(IPIV)
  deallocate(Matrtemp)
  deallocate(Omegacy)
  deallocate(Delta_nnpdp)
  IF(myid==0)then
  deallocate(Omega_matrcy)
  deallocate(Omegakmcy)
  deallocate(F_kcy_re_rt)
  deallocate(F_kcy_im_rt)
  deallocate(GamFun)
  Endif
 ENDIF

  IF(myid==0)then
  IF(verify_L==1)then
  open(stdout+57,file="NLScy.txt",status='replace')
  Endif
  Endif

 if(verify_L==1)then
  NLt1=0
  Mt=0
  Mt1=0
  Mt2=0
  Mt3=0
  CtNL=0
  CtM=0
  NLcount=0
  if(myid==0)then
  call CPU_TIME(Looptb)
  endif
  IF(myid==0)then
  DBt=0
  Endif
 endif

Endif
!=================================================================================================================================
  !Deallocate Public variables

  deallocate(MDmum)
  deallocate(MDmup)
  deallocate(random)
  deallocate(randomGAM)
  deallocate(krho)
  deallocate(nu_DWcode)
  deallocate(nu_ZFcode)
  deallocate(signk)
  deallocate(mu_bar)
  deallocate(mu_pointtemp)
  deallocate(w_pointtemp)
  deallocate(Integ_den)
  deallocate(Omega_nT)

  IF(GK_FK_CK==1 .or. GK_FK_CK==2 .or. GK_FK_CK== -2)then
  deallocate(Omega_d)
  deallocate(expibp)
  deallocate(expibm)
  deallocate(gamIC_k)
  deallocate(delta_k)
  Endif

  IF(GK_FK_CK==1)then
  deallocate(Jn)
  ENDIF
  IF(GK_FK_CK==2 .or. GK_FK_CK== -2)then
  if(myid .ne. 0)then
  deallocate(Jn)
  endif
  ENDIF

  99 format(e14.7)
end Subroutine load
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine nonlinear(myid)

  use mpi
  use parameter
  implicit none
  integer,intent(in) :: myid

!=================================================================================================================================
!GK
IF(GK_FK_CK==0)then
  nonlinear_term(:,:)=(0,0)

  if(myid==1 .and. verify_L==1)then
  call cpu_time(GKNLtb)
  endif

   IF(Add_GAM==1)then
  ! F_k(:,0)=F_k(:,0)+F_GAM(:)
   endif

  if(Const_NL==0) then
  else
  do p=-pmax+1,pmax-1    !the loop of nonlinear term,kx all and 0=<ky<=kymax
    do n=0,nmax-1   !2011.6.15 do the whole loop of the nonliear_term, the image part of verify1 get to 0
        do p1=-pmax+1,pmax-1     !sum_k1
          do n1=-nmax+1,nmax-1
          p2=p-p1
          n2=n-n1
          if((p2.GE.-pmax+1).AND.(p2.LE.pmax-1).AND.(n2.GE.-nmax+1).AND.(n2.LE.nmax-1))then
          IF(myid==0)then
          nonlinear_term(p,n)=nonlinear_term(p,n)+&
                 NLe*(k1x(p1)*k2y(n2)-k2x(p2)*k1y(n1))*phi_k(p1,n1)*F_kmu(p2,n2,myid)
          Else
          nonlinear_term(p,n)=nonlinear_term(p,n)+&
                 NLi*(k1x(p1)*k2y(n2)-k2x(p2)*k1y(n1))*phi_k(p1,n1)*F_kmu(p2,n2,myid)*Jn(p1,n1,myid,0)
          Endif
          endif
          enddo
        enddo
    enddo
  enddo

  do p=-pmax+1,pmax-1
     do n=1,nmax-1
     nonlinear_term(-p,-n)=conjg(nonlinear_term(p,n))
     enddo
  enddo
  endif

  if(myid==1 .and. verify_L==1)then
  call cpu_time(GKNLte)
  GKNLt=GKNLt+ (GKNLte-GKNLtb)
   if(nt==ntmax)then
   write(*,*)"GKNLt= ",GKNLt
   open(stdout+103,file="cpu_timeGK.txt",status='old',POSITION='APPEND')
   write(stdout+103,*)"GKNLt= ",GKNLt
   close(stdout+103)
   endif
  endif

  call MPI_GATHER(nonlinear_term(:,:),(2*pmax-1)*(2*nmax-1),MPI_COMPLEX,&
                  NLSa(:,:,:),(2*pmax-1)*(2*nmax-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)

  IF(Add_GAM==1)then
 ! F_k(:,0)=F_k(:,0)-F_GAM(:)
  Endif

Endif

!=================================================================================================================================
!CKinFH
IF(GK_FK_CK==1)then

  if(verify_L==1)then
  call cpu_time(NLtb1)
  endif

  NLSft(:,:)=(0,0)
  If(Const_NL==0) then
  Else
!Fourier kinetics

 i=myid/(2*N_FT-1)
 a=mod(myid,2*N_FT-1)-(N_FT-1)
 IF(i==0)then
   if(a==0)then
   do n=0,nmax-1   !2011.6.15 do the whole loop of the nonliear_term, the image part of verify1 get to 0
   do p=-pmax+1,pmax-1    !the loop of nonlinear term,kx all and 0=<ky<=kymax
     do n1=-nmax+1,nmax-1
     do p1=-pmax+1,pmax-1     !sum_k1
       p2=p-p1
       n2=n-n1
       if((p2.GE.-pmax+1).AND.(p2.LE.pmax-1).AND.(n2.GE.-nmax+1).AND.(n2.LE.nmax-1))then
      ! write(*,*)"p=",p,"n=",n,"p1=",p1,"n1=",n1,"p2=",p2,"n2=",n2
       NLSft(p,n)=NLSft(p,n)+&
              NLe*(k1x(p1)*k2y(n2)-k2x(p2)*k1y(n1))*phi_kft(p1,n1)*F_kfta(p2,n2,myid)
       endif
       enddo
       enddo
   enddo
   enddo
   endif

 Else If(i.GE.1)then
  phi_kp(:,:)=(0,0)
  phi_km(:,:)=(0,0)
  do n=-nmax+1,nmax-1
  do p=-pmax+1,pmax-1
     phi_kp(p,n)=Cphi_p(p,n)*phi_kft(p,n)
     phi_km(p,n)=Cphi_m(p,n)*phi_kft(p,n)
  enddo
  enddo

   a_p=a+1
   a_m=a-1
   if(a==-N_FT+1)then      !cyclic condition of Fourier harmonic
   a_m=N_FT-1
   else if(a==N_FT-1)then  !cyclic condition of Fourier harmonic
   a_p=-N_FT+1
   endif
   l_p=i*(2*N_FT-1)+a_p+(N_FT-1)
   l_m=i*(2*N_FT-1)+a_m+(N_FT-1)
 !--------------------------------------------------------------------------------
 !Preparing for Conservative form of nonlinear term
   muFftp=(0,0)
   muFftm=(0,0)
   alFftp=(0,0)
   alFftm=(0,0)
  do n=-nmax+1,nmax-1
  do p=-pmax+1,pmax-1
     alFftp(p,n)=a/rho(i)*F_kfta(p,n,l_p)*ratio_TiTe
     alFftm(p,n)=a/rho(i)*F_kfta(p,n,l_m)*ratio_TiTe
  enddo
  enddo

  do ip=1,N_mu
    do n=-nmax+1,nmax-1
    do p=-pmax+1,pmax-1
    muFftp(p,n)=muFftp(p,n) + MDmu(i,ip)*rho(ip)*F_kfta(p,n,ip*(2*N_FT-1)+a_p+(N_FT-1))
    muFftm(p,n)=muFftm(p,n) + MDmu(i,ip)*rho(ip)*F_kfta(p,n,ip*(2*N_FT-1)+a_m+(N_FT-1))
    enddo
    enddo
  enddo
 !--------------------------------------------------------------------------------
 !Preparing for Non-Conservative form of nonlinear term
  N_muFftp=(0,0)
  N_muFftm=(0,0)
  N_alFftp=(0,0)
  N_alFftm=(0,0)
  do n=-nmax+1,nmax-1
  do p=-pmax+1,pmax-1
    N_alFftp(p,n)=a_p/rho(i)*F_kfta(p,n,l_p)*ratio_TiTe
    N_alFftm(p,n)=a_m/rho(i)*F_kfta(p,n,l_m)*ratio_TiTe
  enddo
  enddo

!--------------------------------------------------
  if(myid==10 .and. verify_L==1)then
  open(stdout+57,file="MDmu.txt",status='replace')
  write(stdout+57,*)MDmu
  close(stdout+57)
  endif
!--------------------------------------------------

  do n=-nmax+1,nmax-1
  do p=-pmax+1,pmax-1
    do ip=1,N_mu
    N_muFftp(p,n)=N_muFftp(p,n) + MDmu(i,ip)*F_kfta(p,n,ip*(2*N_FT-1)+a_p+(N_FT-1))
    N_muFftm(p,n)=N_muFftm(p,n) + MDmu(i,ip)*F_kfta(p,n,ip*(2*N_FT-1)+a_m+(N_FT-1))
    enddo
    N_muFftp(p,n)=N_muFftp(p,n)*rho(i)
    N_muFftm(p,n)=N_muFftm(p,n)*rho(i)
  enddo
  enddo
 !--------------------------------------------------------------------------------
 !Calculating the nonlinear terms
 ! NLSft(:,:)=(0,0)
  NLSmuft(:,:)=(0,0)
  NLSalft(:,:)=(0,0)
  N_NLSmuft(:,:)=(0,0)
  N_NLSalft(:,:)=(0,0)

  do n=0,nmax-1   !2011.6.15 do the whole loop of the nonliear_term, the image part of verify1 get to 0
  do p=-pmax+1,pmax-1    !the loop of nonlinear term,kx all and 0=<ky<=kymax
    do n1=-nmax+1,nmax-1
    do p1=-pmax+1,pmax-1     !sum_k1
      p2=p-p1
      n2=n-n1
      if((p2.GE.-pmax+1).AND.(p2.LE.pmax-1).AND.(n2.GE.-nmax+1).AND.(n2.LE.nmax-1))then

     ! NLSft(p,n)=NLSft(p,n)&
     !               +phi_km(p1,n1)*muFftm(p2,n2) + phi_kp(p1,n1)*muFftp(p2,n2)&
     !               -phi_km(p1,n1)*alFftm(p2,n2) + phi_kp(p1,n1)*alFftp(p2,n2)

      NLSmuft(p,n)=NLSmuft(p,n)+&
                      phi_km(p1,n1)*muFftm(p2,n2) + phi_kp(p1,n1)*muFftp(p2,n2)

      NLSalft(p,n)=NLSalft(p,n)&
                      -phi_km(p1,n1)*alFftm(p2,n2)+ phi_kp(p1,n1)*alFftp(p2,n2)


      N_NLSmuft(p,n)=N_NLSmuft(p,n)+&
                      phi_km(p1,n1)*N_muFftm(p2,n2)+ phi_kp(p1,n1)*N_muFftp(p2,n2)

      N_NLSalft(p,n)=N_NLSalft(p,n)&
                      -phi_km(p1,n1)*N_alFftm(p2,n2) + phi_kp(p1,n1)*N_alFftp(p2,n2)
      endif
    enddo
    enddo
  enddo
  enddo

  NLSft(:,:)= NLi*0.5*( NLSmuft(:,:)+NLSalft(:,:) + N_NLSmuft(:,:)+ N_NLSalft(:,:) )

  Endif

  Endif


  if(verify_L==1)then
  call cpu_time(NLte1)
  NLt1=NLt1+ (NLte1-NLtb1)
   if(nt==ntmax)then
   NLcount=NLcount+1
     if(NLcount==2)then
     i=myid/(2*N_FT-1)
     a=mod(myid,2*N_FT-1)-(N_FT-1)
     write(*,*)"i,a,myid=",i,a,myid,"  NLt=",NLt1
     endif
   endif
  endif

  if(verify_L==1)then
  call cpu_time(CtNLb)
  endif

  call MPI_GATHER(NLSft(:,:),(2*pmax-1)*nmax,MPI_COMPLEX,&
                  NLSfta(:,:,:),(2*pmax-1)*nmax,MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(NLSfta,(2*pmax-1)*nmax*(N_mu+1)*(2*N_FT-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)

!--------------------------------------------------
  if(myid==0 .and. verify_L==1)then
  open(stdout+57,file="NLSfta.txt",status='replace')
  write(stdout+57,*)NLSfta
  close(stdout+57)
  endif
!--------------------------------------------------
  if(verify_L==1)then
  call cpu_time(CtNLe)
  CtNL=CtNL+ (CtNLe-CtNLb)
   if(nt==ntmax)then
     if(NLcount==2)then
     i=myid/(2*N_FT-1)
     a=mod(myid,2*N_FT-1)-(N_FT-1)
     write(*,*)"i,a,myid=",i,a,myid,"  CtNL=",CtNL
     endif
   endif
  endif

Endif

!=================================================================================================================================
!CKinCH
IF(GK_FK_CK==2)then

  if(verify_L==1)then
  call cpu_time(NLtb1)
  endif

  NLScy=(0,0)
  if(myid==0)then
  NLScya(:,:,:)=(0,0)
  endif

 If(Const_NL==0) then
 Else
 i=mod(myid,N_mu+1)
 a=myid/(N_mu+1)-(N_CY-1)
 IF(i==0)then
   if(a==0)then
   do n=0,nmax-1   !2011.6.15 do the whole loop of the nonliear_term, the image part of verify1 get to 0
   do p=-pmax+1,pmax-1    !the loop of nonlinear term,kx all and 0=<ky<=kymax
     do n1=-nmax+1,nmax-1
     do p1=-pmax+1,pmax-1
       p2=p-p1
       n2=n-n1
       if((p2.GE.-pmax+1).AND.(p2.LE.pmax-1).AND.(n2.GE.-nmax+1).AND.(n2.LE.nmax-1))then
       NLScy(p,n)=NLScy(p,n)+&
              (k1x(p1)*k2y(n2)-k2x(p2)*k1y(n1))*Phi_kcy(p1,n1)*F_kcya(myid,p2,n2)
       endif
     enddo
     enddo
   enddo
   enddo
   NLScy(:,:)=NLScy(:,:)*NLe
   endif
 Else

  do n=0,nmax-1
  do p=-pmax+1,pmax-1
    do n1=-nmax+1,nmax-1
    do p1=-pmax+1,pmax-1     !sum_k1
    p2=p-p1
    n2=n-n1
      if((p2.GE.-pmax+1).AND.(p2.LE.pmax-1).AND.(n2.GE.-nmax+1).AND.(n2.LE.nmax-1))then
      NLkk1=(0,0)
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        do ap=-N_CY+1,N_CY-1
        lp=(ap+N_CY-1)*(N_mu+1)+i
          NLkk1=NLkk1+ F_kcya(lp,p2,n2)*NL1(ap,p1,n1,p,n)
        enddo
      !-----------------------------------------------------------------------
        do ap=-N_CY+1,N_CY-1
        do ip=1,N_mu
          lpp=(ap+N_CY-1)*(N_mu+1)+ip
          NLkk1=NLkk1+ F_kcya(lpp,p2,n2)*NL2(ip,ap,p1,n1,p,n)
        enddo
        enddo
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      NLScy(p,n)=NLScy(p,n)+Phi_kcy(p1,n1)*NLkk1
      endif
    enddo
    enddo
  enddo
  enddo
  NLScy=NLScy/2.0
 Endif
 Endif


  if(verify_L==1)then
  call cpu_time(NLte1)
  NLt1=NLt1+ (NLte1-NLtb1)
   if(nt==ntmax)then
   NLcount=NLcount+1
     if(NLcount==2)then
     i=mod(myid,N_mu+1)
     a=myid/(N_mu+1)-(N_CY-1)
     write(*,*)"i,a,myid=",i,a,myid,"  NLt=",NLt1
     endif
   endif
  endif

  if(verify_L==1)then
  call cpu_time(CtNLb)
  endif

  call MPI_GATHER(NLScy(:,:),(2*pmax-1)*nmax,MPI_COMPLEX,&
                 NLScyatemp(:,:,:),(2*pmax-1)*nmax,MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)
  If(myid==0)then
    do n=0,nmax-1
    do p=-pmax+1,pmax-1
    NLScya(:,p,n)=NLScyatemp(p,n,:)
    enddo
    enddo
  Endif
  call MPI_BCAST(NLScya,(2*pmax-1)*nmax*(N_mu+1)*(2*N_CY-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)

  if(verify_L==1)then
  call cpu_time(CtNLe)
  CtNL=CtNL+ (CtNLe-CtNLb)
   if(nt==ntmax)then
     if(NLcount==2)then
     i=mod(myid,N_mu+1)
     a=myid/(N_mu+1)-(N_CY-1)
     write(*,*)"i,a,myid=",i,a,myid,"  CtNL=",CtNL
     endif
   endif
  endif

  IF(verify_L==1)then
  IF(myid==0)then
  do n=0,nmax-1
  do p=-pmax+1,pmax-1
    do a=-N_CY+1,N_CY-1
    do i=0,N_mu
     l=(a+N_CY-1)*(N_mu+1)+i
     write(stdout+57,*)'i,a,p,n=',i,a,p,n
     write(stdout+57,*)&
     real(NLScya(l,p,n)),aimag(NLScya(l,p,n))
    enddo
    enddo
  enddo
  enddo
  EndIF
  Endif

Endif
!=================================================================================================================================
!CKinCH deeply paralized
IF(GK_FK_CK== -2)then

  if(verify_L==1)then
  call cpu_time(NLtb1)
  endif

  NLScydp=(0,0)
  if(myid==0)then
  NLScyadp=(0,0)
  endif

 If(Const_NL==0) then
 Else
  i=mod(myid,N_mu+1)
  a=mod((myid-i)/(N_mu+1),(2*N_CY-1))-(N_CY-1)
  p=myid/((N_mu+1)*(2*N_CY-1))-(pmax-1)
 IF(i==0)then
   if(a==0)then
   l=(a+N_CY-1)*(N_mu+1)+i
   do n=0,nmax-1   !2011.6.15 do the whole loop of the nonliear_term, the image part of verify1 get to 0
     do n1=-nmax+1,nmax-1
     do p1=-pmax+1,pmax-1
       p2=p-p1
       n2=n-n1
       if((p2.GE.-pmax+1).AND.(p2.LE.pmax-1).AND.(n2.GE.-nmax+1).AND.(n2.LE.nmax-1))then
       NLScydp(n)=NLScydp(n)+&
              (k1x(p1)*k2y(n2)-k2x(p2)*k1y(n1))*Phi_kcy(p1,n1)*F_kcya(l,p2,n2)
       endif
     enddo
     enddo
   enddo
   NLScydp(:)=NLScydp(:)*NLe
   endif
 Else

  do n=0,nmax-1
    do n1=-nmax+1,nmax-1
    do p1=-pmax+1,pmax-1     !sum_k1
    p2=p-p1
    n2=n-n1
      if((p2.GE.-pmax+1).AND.(p2.LE.pmax-1).AND.(n2.GE.-nmax+1).AND.(n2.LE.nmax-1))then
      NLkk1=(0,0)
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        do ap=-N_CY+1,N_CY-1
        lp=(ap+N_CY-1)*(N_mu+1)+i
          NLkk1=NLkk1+ F_kcya(lp,p2,n2)*NL1dp(ap,p1,n1,n)
        enddo
      !-----------------------------------------------------------------------
        do ap=-N_CY+1,N_CY-1
        do ip=1,N_mu
          lpp=(ap+N_CY-1)*(N_mu+1)+ip
          NLkk1=NLkk1+ F_kcya(lpp,p2,n2)*NL2dp(ip,ap,p1,n1,n)
        enddo
        enddo
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      NLScydp(n)=NLScydp(n)+Phi_kcy(p1,n1)*NLkk1
      endif
    enddo
    enddo
  enddo
  NLScydp=NLScydp/2.0
 Endif
 Endif

  if(verify_L==1)then
  call cpu_time(NLte1)
  NLt1=NLt1+ (NLte1-NLtb1)
   if(nt==ntmax)then
   NLcount=NLcount+1
     if(NLcount==2)then
     write(*,*)"myid=",myid,"  NLt=",NLt1
     endif
   endif
  endif

  if(verify_L==1)then
  call cpu_time(CtNLb)
  endif

  call MPI_GATHER(NLScydp(:),nmax,MPI_COMPLEX,&
                 NLScyadptemp(:,:),nmax,MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)
  If(myid==0)then
    do n=0,nmax-1
    NLScyadp(:,n)=NLScyadptemp(n,:)
    enddo
  Endif
  call MPI_BCAST(NLScyadp,nmax*(2*pmax-1)*(N_mu+1)*(2*N_CY-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)

  if(verify_L==1)then
  call cpu_time(CtNLe)
  CtNL=CtNL+ (CtNLe-CtNLb)
   if(nt==ntmax)then
     if(NLcount==2)then
     write(*,*)"myid=",myid,"  CtNL=",CtNL
     endif
   endif
  endif

  IF(verify_L==1)then
  IF(myid==0)then
  do n=0,nmax-1
  do p=-pmax+1,pmax-1
    do a=-N_CY+1,N_CY-1
    do i=0,N_mu
    jL=(p+pmax-1)*(N_mu+1)*(2*N_CY-1)+(a+N_CY-1)*(N_mu+1)+i
    write(stdout+57,*)'i,a,p,n=',i,a,p,n
     write(stdout+57,*)&
     real(NLScyadp(jL,n)),aimag(NLScyadp(jL,n))
     enddo
  enddo
  enddo
  enddo
  EndIF
  Endif


Endif
!=================================================================================================================================
end Subroutine nonlinear
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




Subroutine motion(myid)

  use mpi
  use parameter
  implicit none
  integer,intent(in) :: myid
  complex(singlep) :: F_ktemp

!=================================================================================================================================
!GK
IF(GK_FK_CK==0)then
   !store the value before step, as for diagnose


    F_k=(0,0)
    phi_k=(0,0)

    IF(Add_GAM==1)then
    F_GAMa(:) = F_GAM(:)
    F_GAM=(0,0)
    Endif

 If(myid==0)then

  if(verify_L==1)then
  call cpu_time(GKMtb)
  endif


 IF(abs(muDtype) .ne. 2)then
  do p=-pmax+1,pmax-1
    do n=0,nmax-1
      Source(:)=(0,0)
      do  l=0,N_mu
        do  j=0,N_mu
        Source(l)=Source(l)+M_Matrix(p,n,l,j)*F_kmu(p,n,j)
        enddo
      enddo
      F_kmu(p,n,:)=(0,0)
      do  l=0,N_mu
        do  j=0,N_mu
        F_kmu(p,n,l)=F_kmu(p,n,l)+R_Matrix(p,n,l,j)*(Source(j)+tstep*NLSa(p,n,j))
        enddo
      enddo
    enddo
  enddo
 Else if (muDtype==2)then
  do p=-pmax+1,pmax-1
    do n=0,nmax-1
      do  l=0,N_mu
      F_kmu(p,n,l)=F_kmu(p,n,l)+ tstep*NLSa(p,n,l)
      enddo
    enddo
  enddo
 Else if (muDtype==-2)then
  do p=-pmax+1,pmax-1
    do n=0,nmax-1
      do  l=0,N_mu
      F_ktemp=(0,0)
        do  j=0,N_mu
        F_ktemp=F_ktemp+M_Matrix(p,n,l,j)*F_kmu(p,n,j)
        enddo
      F_kmu(p,n,l)=F_kmu(p,n,l)+ tstep*(F_ktemp+NLSa(p,n,l))
      enddo
    enddo
  enddo
 Endif


       do p=-pmax+1,pmax-1
         do n=0,nmax-1
      F_kmu(-p,-n,:)=conjg(F_kmu(p,n,:))
         enddo
       enddo
       F_kmu(0,0,:)=(0,0)


    do p=-pmax+1,pmax-1
      do n=0,nmax-1
       if(n==0.and.p==0)then
       phi_k(0,0)=(0,0)!the real part of (0,0)mode do not change, the image part go to 0 because of conjugation!!
       else
        Integ_num(p,n)=(0,0)
        do i=1,N_mu
        Integ_num(p,n)=Integ_num(p,n)+w_point(i)*Jn(p,n,i,0)*F_kmu(p,n,i)
        enddo
        i=0
        Integ_num(p,n)=Integ_num(p,n)-CDW*F_kmu(p,n,i)

        phi_k(p,n)=Integ_num(p,n)/H_m_delta(p,n)    !!!***===
       endif
      enddo
    enddo

     do p=-pmax+1,pmax-1
       do n=0,nmax-1
       phi_k(-p,-n)=conjg(phi_k(p,n))
       enddo
     enddo

  if(verify_L==1)then
  call cpu_time(GKMte)
  GKMt=GKMt+ (GKMte-GKMtb)
   if(nt==ntmax)then
   write(*,*)"GK Motion time= ",GKMt
   open(stdout+103,file="cpu_timeGK.txt",status='old',POSITION='APPEND')
   write(stdout+103,*)"GK Motion time= ",GKMt
   close(stdout+103)
   endif
  endif

  Endif


  call MPI_BCAST(phi_k,(2*pmax-1)*(2*nmax-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process
  call MPI_BCAST(F_kmu(:,:,:),(2*pmax-1)*(2*nmax-1)*(N_mu+1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process

    IF(Add_GAM==1)then
    F_k(:,0)=F_k(:,0)-F_GAM(:)
    Endif


Endif

!=================================================================================================================================
!CKinFH
IF(GK_FK_CK==1)then
    !4D Fourier harmonic

  if(verify_L==1)then
  call cpu_time(Mtb)
  endif

   do ki=0,kimax-1
   j=myid*kimax+ki
   p=mod(j,(2*pmax-1))-(pmax-1)
   n=j/(2*pmax-1)
   if (n .LT. nmax)then
    IF(abs(muDtype) .ne. 2)then
     Sourceft(:)=(0,0)
      do  j=0,(N_mu+1)*(2*N_FT-1)-1
      do  l=0,(N_mu+1)*(2*N_FT-1)-1
        Sourceft(l)=Sourceft(l)+M_Matrixft(ki,l,j)*F_kfta(p,n,j)
      enddo
      enddo

   !   Sourceft(:)=MATMUL(R_Matrixft(ki,:,:),Sourceft(:)+tstep*NLSfta(p,n,:) )

      Fkmyid(ki,:)=(0,0)
      do  j=0,(N_mu+1)*(2*N_FT-1)-1
      do  l=0,(N_mu+1)*(2*N_FT-1)-1
        Fkmyid(ki,l)=Fkmyid(ki,l)+R_Matrixft(ki,l,j)*(Sourceft(j)+tstep*NLSfta(p,n,j))
      enddo
      enddo

      do a=-N_FT+1,N_FT-1
      i=0
      j=i*(2*N_FT-1)+a+(N_FT-1)
      if(a .ne. 0)then
      Fkmyid(ki,j)=(0,0)
      endif
      enddo

     Else if (muDtype==2)then
      Fkmyid(ki,:)=(0,0)
      Fkmyid(ki,:)=F_kfta(p,n,:)+ tstep*NLSfta(p,n,:)

      do a=-N_FT+1,N_FT-1
      i=0
      j=i*(2*N_FT-1)+a+(N_FT-1)
      if(a .ne. 0)then
      Fkmyid(ki,j)=(0,0)
      endif
      enddo

     Else if (muDtype==-2)then
      Fkmyid(ki,:)=(0,0)
      do  l=0,(N_mu+1)*(2*N_FT-1)-1
      F_ktemp=(0,0)
        do  j=0,(N_mu+1)*(2*N_FT-1)-1
        F_ktemp=F_ktemp+M_Matrixft(ki,l,j)*F_kfta(p,n,j)
        enddo
      Fkmyid(ki,l)=F_kfta(p,n,l)+ tstep*(F_ktemp+NLSfta(p,n,l))
      enddo

      do a=-N_FT+1,N_FT-1
      i=0
      j=i*(2*N_FT-1)+a+(N_FT-1)
      if(a .ne. 0)then
      Fkmyid(ki,j)=(0,0)
      endif
      enddo
     Endif

  endif
  enddo

  if(verify_L==1)then
  call cpu_time(Mtm1)
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  call MPI_GATHER(Fkmyid(:,:),kimax*(N_mu+1)*(2*N_FT-1),MPI_COMPLEX,&
                  Fkkmft(:,:,:),kimax*(N_mu+1)*(2*N_FT-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)

  if(verify_L==1)then
  call cpu_time(Mtm2)
  endif
  if(verify_L==2)then
    If(myid==0)then
    F_kftbef(:,:)=F_kfta(0,:,:)
    Endif
  endif


  F_kfta(:,:,:)=(0,0)
  If(myid==0)then
      do ki=0,kimax-1
      do i=0,(N_mu+1)*(2*N_FT-1)-1
        j=i*kimax+ki
        p=mod(j,(2*pmax-1))-(pmax-1)
        n=j/(2*pmax-1)
        if (n .LT. nmax)then
        F_kfta(p,n,:) =Fkkmft(ki,:,i)
        endif
      enddo
      enddo

     do l=0,(N_mu+1)*(2*N_FT-1)-1
     i=l/(2*N_FT-1)
     a=mod(l,2*N_FT-1)-(N_FT-1)
     j=i*(2*N_FT-1)-a+(N_FT-1)
       do n=0,nmax-1
       do p=-pmax+1,pmax-1
       F_kfta(-p,-n,j)=conjg(F_kfta(p,n,l))
       enddo
     enddo
     enddo
     F_kfta(0,0,:)=(0,0)

     do n=0,nmax-1
        if(n/=0)then
        lambda_k=lambda_n
        else
        lambda_k=lambda_0
        endif
     do p=-pmax+1,pmax-1
        if(n==0.and.p==0)then
        phi_kft(0,0)=(0,0)  !4D Fourier harmonic
        else
        !4D Fourier harmonic
        phi_kft(p,n)=(0,0)
        do i=1,N_mu
        a=0
        l=i*(2*N_FT-1)+a+(N_FT-1)
        phi_kft(p,n)=phi_kft(p,n)+w_point(i)/Beta*F_kfta(p,n,l)
        enddo
        a=0
        i=0
        l=i*(2*N_FT-1)+a+(N_FT-1)
        phi_kft(p,n)=phi_kft(p,n)-CDW*F_kfta(p,n,l)

        phi_kft(p,n)=phi_kft(p,n)/lamb_m_delta(p,n)
        endif
     enddo
     enddo

     do n=0,nmax-1
     do p=-pmax+1,pmax-1
       phi_kft(-p,-n)=conjg(phi_kft(p,n))  !4D Fourier harmonic
     enddo
     enddo
  Endif



  if(verify_L==1)then
  call cpu_time(Mte)
  Mt=Mt+ (Mte-Mtb)
  Mt1=Mt1+ (Mtm1-Mtb)
  Mt2=Mt2+ (Mtm2-Mtm1)
  Mt3=Mt3+ (Mte-Mtm2)
   if(nt==ntmax)then
   write(*,*)"myid= ",myid,"Mt=",Mt,":  Mt1=",Mt1," Mt2=",Mt2," Mt3=",Mt3
   endif
  call cpu_time(CtMb)
  endif

  if(verify_L==2)then
  If(myid==0)then
    Rate=(0,0)
    Ratave=(0,0)
    delRat=(0,0)
    do n=1,nmax-1
      do i=1,N_mu
      do a=-N_FT+1,N_FT-1
      j=i*(2*N_FT-1)+a+(N_FT-1)
      l=(i-1)*(2*N_FT-1)+a+(N_FT-1)
      Rate(n,l)=2.0*real((F_kfta(0,n,j)-F_kftbef(n,j))/(F_kfta(0,n,j)+F_kftbef(n,j))/tstep)
      enddo
      enddo
    enddo

    do n=1,nmax-1
    Ratave(n)=sum(Rate(n,:))/(N_mu*(2*N_FT-1))
    delRat(n)=sqrt(abs(Sum((Rate(n,:)-Ratave(n)))**2)/(N_mu*(2*N_FT-1)))
    enddo

    if(nt==1.and.restart==0)then
    open(stdout+60,file="Ratave.txt",status='replace')
    open(stdout+61,file="delRat.txt",status='replace')
    endif
    write(stdout+60,*) Ratave
    write(stdout+61,*) delRat

    do n=-nmax+1,nmax-1
    do p=-pmax+1,pmax-1
    phi_kft(p,n)=phi_kft(p,n)/gammax(p,n)/Cnorm
    F_kfta(p,n,:)=F_kfta(p,n,:)/gammax(p,n)/Cnorm
    enddo
    enddo

  Endif
  endif

  call MPI_BCAST(phi_kft,(2*pmax-1)*(2*nmax-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process call MPI_BCAST(F_kfta(:,:,:),(2*pmax-1)*(2*nmax-1)*(N_mu*(2*N_FT-1)),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process
  call MPI_BCAST(F_kfta(:,:,:),(2*pmax-1)*(2*nmax-1)*(N_mu+1)*(2*N_FT-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process


  if(verify_L==1)then
  call cpu_time(CtMe)
  CtM=CtM+ (CtMe-CtMb)
   if(nt==ntmax)then
   i=myid/(2*N_FT-1)
   a=mod(myid,2*N_FT-1)-(N_FT-1)
   write(*,*)"i,a,myid=",i,a,myid,"  CtM=",CtM
   endif
  endif

Endif

!=================================================================================================================================
!CKinCH
IF(GK_FK_CK==2)then

  if(verify_L==1)then
  call cpu_time(Mtb)
  endif

 do ki=0,kimax-1
 j=(myid-1)*kimax+ki
 IF(j .GE. 0)then
  p=mod(j,(2*pmax-1))-(pmax-1)
  n=j/(2*pmax-1)

  if (n .LT. nmax)then
  IF(abs(muDtype) .ne. 2)then
    Sourcecy(:)=(0,0)
      do  l=0,(N_mu+1)*(2*N_CY-1)-1
      do  j=0,(N_mu+1)*(2*N_CY-1)-1
        Sourcecy(l)=Sourcecy(l)+M_Matrixcy(j,l,ki)*F_kcya(j,p,n)
      enddo
      enddo

      Fkmyidcy(:,ki)=(0,0)
      do  l=0,(N_mu+1)*(2*N_CY-1)-1
      do  j=0,(N_mu+1)*(2*N_CY-1)-1
        Fkmyidcy(l,ki)=Fkmyidcy(l,ki)+R_Matrixcy(j,l,ki)*(Sourcecy(j)+tstep*NLScya(j,p,n))
      enddo
      enddo

      do a=-N_CY+1,N_CY-1
      i=0
      j=(a+N_CY-1)*(N_mu+1)+i
      if(a .ne. 0)then
      Fkmyidcy(j,ki)=(0,0)
      endif
      enddo

  Else if(muDtype==2)then
      Fkmyidcy(:,ki)=(0,0)
      Fkmyidcy(:,ki)=F_kcya(:,p,n)+ tstep*NLScya(:,p,n)

      do a=-N_CY+1,N_CY-1
      i=0
      j=(a+N_CY-1)*(N_mu+1)+i
      if(a .ne. 0)then
      Fkmyidcy(j,ki)=(0,0)
      endif
      enddo

  Else if(muDtype==-2)then
      Fkmyidcy(:,ki)=(0,0)
      do  l=0,(N_mu+1)*(2*N_CY-1)-1
      F_ktemp=(0,0)
        do  j=0,(N_mu+1)*(2*N_CY-1)-1
        F_ktemp=F_ktemp+M_Matrixcy(j,l,ki)*F_kcya(j,p,n)
        enddo
      Fkmyidcy(l,ki)=F_kcya(l,p,n)+ tstep*(F_ktemp+NLScya(l,p,n))
      enddo

      do a=-N_CY+1,N_CY-1
      i=0
      j=(a+N_CY-1)*(N_mu+1)+i
      if(a .ne. 0)then
      Fkmyidcy(j,ki)=(0,0)
      endif
      enddo
  Endif
  endif
 ENDIF
 enddo

  if(verify_L==1)then
  call cpu_time(Mtm1)
  endif

  call MPI_GATHER(Fkmyidcy(:,:),kimax*(N_mu+1)*(2*N_CY-1),MPI_COMPLEX,&
                  Fkkmcy(:,:,:),kimax*(N_mu+1)*(2*N_CY-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)

  if(verify_L==1)then
  call cpu_time(Mtm2)
  endif

  F_kcya(:,:,:)=(0,0)
  If(myid==0)then
    do i=0,(N_mu+1)*(2*N_CY-1)-1
    do ki=0,kimax-1
    j=(i-1)*kimax+ki
    If(j .GE. 0)then
      p=mod(j,(2*pmax-1))-(pmax-1)
      n=j/(2*pmax-1)
      if (n .LT. nmax)then
      F_kcya(:,p,n) =Fkkmcy(:,ki,i)
      endif
    Endif
    enddo
    enddo

    do n=0,nmax-1
    do p=-pmax+1,pmax-1
      do l=0,(N_mu+1)*(2*N_CY-1)-1
      i=mod(l,N_mu+1)
      a=l/(N_mu+1)-(N_CY-1)
      j=(-a+N_CY-1)*(N_mu+1)+i
      F_kcya(j,-p,-n)=conjg(F_kcya(l,p,n))*((-1)**a)
      enddo
    enddo
    enddo
    F_kcya(:,0,0)=(0,0)

    do n=0,nmax-1
        if(n/=0)then
        lambda_k=lambda_n
        else
        lambda_k=lambda_0
        endif
    do p=-pmax+1,pmax-1
        if(n==0.and.p==0)then
        phi_kcy(0,0)=(0,0)  !cyclo-kinetic
        else
        !cyclo-kinetic
        phi_kcy(p,n)=(0,0)
        do a=-N_CY+1,N_CY-1
          do i=1,N_mu
          l=(a+N_CY-1)*(N_mu+1)+i
          phi_kcy(p,n)=phi_kcy(p,n)+w_point(i)*Jn(p,n,i,a)*F_kcya(l,p,n)
          enddo
        enddo
        i=0
        a=0
        l=(a+N_CY-1)*(N_mu+1)+i
        phi_kcy(p,n)=phi_kcy(p,n)-CDW*F_kcya(l,p,n)

        phi_kcy(p,n)=phi_kcy(p,n)/G_m_delta(p,n)
        endif
     enddo
     enddo

     do n=0,nmax-1
     do p=-pmax+1,pmax-1
       phi_kcy(-p,-n)=conjg(phi_kcy(p,n))
     enddo
     enddo

 Endif

  if(verify_L==1)then
  call cpu_time(Mte)
  Mt=Mt+ (Mte-Mtb)
  Mt1=Mt1+ (Mtm1-Mtb)
  Mt2=Mt2+ (Mtm2-Mtm1)
  Mt3=Mt3+ (Mte-Mtm2)
   if(nt==ntmax)then
   write(*,*)"myid= ",myid,"Mt=",Mt
   write(*,*)"Mt1=",Mt1," Mt2=",Mt2," Mt3=",Mt3
   endif
  call cpu_time(CtMb)
  endif

  call MPI_BCAST(phi_kcy,(2*pmax-1)*(2*nmax-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process call MPI_BCAST(F_kfta(:,:,:),(2*pmax-1)*(2*nmax-1)*(N_mu*(2*N_FT-1)),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process
  call MPI_BCAST(F_kcya(:,:,:),(2*pmax-1)*(2*nmax-1)*(N_mu+1)*(2*N_CY-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process

  if(verify_L==1)then
  call cpu_time(CtMe)
  CtM=CtM+ (CtMe-CtMb)
   if(nt==ntmax)then
   write(*,*)"myid=",myid,"  CtM=",CtM
   endif
  endif

Endif
!=================================================================================================================================
!CKinCH deeply paralized
IF(GK_FK_CK== -2)then

  if(verify_L==1)then
  call cpu_time(Mtb)
  endif

 do ki=0,kimax-1
 j=(myid-1)*kimax+ki
 If(j .GE. 0)then
  p=mod(j,(2*pmax-1))-(pmax-1)
  n=j/(2*pmax-1)

  if (n .LT. nmax)then
  IF(muDtype .ne. 2)then
    Sourcecy(:)=(0,0)
      do  l=0,(N_mu+1)*(2*N_CY-1)-1
      do  j=0,(N_mu+1)*(2*N_CY-1)-1
        Sourcecy(l)=Sourcecy(l)+M_Matrixcy(j,l,ki)*F_kcya(j,p,n)
      enddo
      enddo

      Fkmyidcy(:,ki)=(0,0)
      do  l=0,(N_mu+1)*(2*N_CY-1)-1
        do  a=-N_CY+1,N_CY-1
        do  i=0,N_mu
       j=(a+N_CY-1)*(N_mu+1)+i
       jL=(p+pmax-1)*(N_mu+1)*(2*N_CY-1)+(a+N_CY-1)*(N_mu+1)+i
       Fkmyidcy(l,ki)=Fkmyidcy(l,ki)+R_Matrixcy(j,l,ki)*(Sourcecy(j)+tstep*NLScyadp(jL,n))
        enddo
        enddo
      enddo

      do a=-N_CY+1,N_CY-1
      i=0
      j=(a+N_CY-1)*(N_mu+1)+i
      if(a .ne. 0)then
      Fkmyidcy(j,ki)=(0,0)
      endif
      enddo

  Else if(muDtype==2)then
      Fkmyidcy(:,ki)=(0,0)
      do  a=-N_CY+1,N_CY-1
      do  i=0,N_mu
      j=(a+N_CY-1)*(N_mu+1)+i
      jL=(p+pmax-1)*(N_mu+1)*(2*N_CY-1)+(a+N_CY-1)*(N_mu+1)+i
      Fkmyidcy(j,ki)=F_kcya(j,p,n)+ tstep*NLScyadp(jL,n)
      enddo
      enddo

  Else if(muDtype==-2)then
      Fkmyidcy(:,ki)=(0,0)
      do  l=0,(N_mu+1)*(2*N_CY-1)-1
      F_ktemp=(0,0)
        do  j=0,(N_mu+1)*(2*N_CY-1)-1
        F_ktemp=F_ktemp+M_Matrixcy(j,l,ki)*F_kcya(j,p,n)
        enddo
        lL=(p+pmax-1)*(N_mu+1)*(2*N_CY-1)+l
        Fkmyidcy(l,ki)=F_kcya(l,p,n)+ tstep*(F_ktemp+NLScyadp(lL,n))
      enddo
  Endif
  endif
 Endif
 enddo

  if(verify_L==1)then
  call cpu_time(Mtm1)
  endif

  call MPI_GATHER(Fkmyidcy(:,:),kimax*(N_mu+1)*(2*N_CY-1),MPI_COMPLEX,&
                  Fkkmcy(:,:,:),kimax*(N_mu+1)*(2*N_CY-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)

  if(verify_L==1)then
  call cpu_time(Mtm2)
  endif

  F_kcya(:,:,:)=(0,0)
  If(myid==0)then
   do i=0,(2*pmax-1)*(N_mu+1)*(2*N_CY-1)-1
   do ki=0,kimax-1
   j=(i-1)*kimax+ki
   If(j .GE. 0)then
   p=mod(j,(2*pmax-1))-(pmax-1)
   n=j/(2*pmax-1)
     if (n .LT. nmax)then
     F_kcya(:,p,n) =Fkkmcy(:,ki,i)
     endif
   Endif
   enddo
   enddo

   do n=0,nmax-1
   do p=-pmax+1,pmax-1
     do l=0,(N_mu+1)*(2*N_CY-1)-1
     i=mod(l,N_mu+1)
     a=l/(N_mu+1)-(N_CY-1)
     j=(-a+N_CY-1)*(N_mu+1)+i
     F_kcya(j,-p,-n)=conjg(F_kcya(l,p,n))*((-1)**a)
     enddo
   enddo
   enddo
    F_kcya(:,0,0)=(0,0)

    do n=0,nmax-1
        if(n/=0)then
        lambda_k=lambda_n
        else
        lambda_k=lambda_0
        endif
    do p=-pmax+1,pmax-1
        if(n==0.and.p==0)then
        phi_kcy(0,0)=(0,0)  !cyclo-kinetic
        else
        !cyclo-kinetic
        phi_kcy(p,n)=(0,0)
        do a=-N_CY+1,N_CY-1
          do i=1,N_mu
          l=(a+N_CY-1)*(N_mu+1)+i
          phi_kcy(p,n)=phi_kcy(p,n)+w_point(i)*Jn(p,n,i,a)*F_kcya(l,p,n)
          enddo
        enddo
        i=0
        a=0
        l=(a+N_CY-1)*(N_mu+1)+i
        phi_kcy(p,n)=phi_kcy(p,n)-CDW*F_kcya(l,p,n)

        phi_kcy(p,n)=phi_kcy(p,n)/G_m_delta(p,n)
        endif
     enddo
     enddo

     do n=0,nmax-1
     do p=-pmax+1,pmax-1
       phi_kcy(-p,-n)=conjg(phi_kcy(p,n))
     enddo
     enddo

 Endif

  if(verify_L==1)then
  call cpu_time(Mte)
  Mt=Mt+ (Mte-Mtb)
  Mt1=Mt1+ (Mtm1-Mtb)
  Mt2=Mt2+ (Mtm2-Mtm1)
  Mt3=Mt3+ (Mte-Mtm2)
   if(nt==ntmax)then
   write(*,*)"myid= ",myid,"Mt=",Mt
   write(*,*)"Mt1=",Mt1," Mt2=",Mt2," Mt3=",Mt3
   endif
  call cpu_time(CtMb)
  endif

  call MPI_BCAST(phi_kcy,(2*pmax-1)*(2*nmax-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process call MPI_BCAST(F_kcya(:,:,:),(2*pmax-1)*(2*nmax-1)*(N_mu*(2*N_FT-1)),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process
  call MPI_BCAST(F_kcya(:,:,:),(2*pmax-1)*(2*nmax-1)*(N_mu+1)*(2*N_CY-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)  !0 mark of the root process

  if(verify_L==1)then
  call cpu_time(CtMe)
  CtM=CtM+ (CtMe-CtMb)
   if(nt==ntmax)then
   write(*,*)"myid=",myid,"  CtM=",CtM
   endif
  endif

Endif
!=================================================================================================================================

end Subroutine motion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Subroutine backup(myid)

  use mpi
  use parameter
  implicit none
  integer,intent(in) :: myid

!=================================================================================================================================
!GK
IF(GK_FK_CK==0)then

  IF(Add_GAM==1)then
  F_GAMrec(:,:)=(0,0)
  call MPI_GATHER(F_GAM(:),(2*pmax-1),MPI_COMPLEX,&
                  F_GAMrec(:,:),(2*pmax-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)
  Endif

 If(myid==0)then
  if(nt_Odd==min(nt_Odd,nt_Even))then
    open(1000,file='Fk_re_Odd',status='replace')
    open(1001,file='Fk_im_Odd',status='replace')
    IF(Add_GAM==1)then
    open(1004,file='FGAM_re_Odd',status='replace')
    open(1005,file='FGAM_im_Odd',status='replace')
    Endif
    open(1003,file='nt_Odd.txt',status='replace')
    write(1003,*)nt
    close(1003)
    nt_Odd=nt
  elseif(nt_Even==min(nt_Odd,nt_Even))then
    open(1000,file='Fk_re_Even',status='replace')
    open(1001,file='Fk_im_Even',status='replace')
    IF(Add_GAM==1)then
    open(1004,file='FGAM_re_Odd',status='replace')
    open(1005,file='FGAM_im_Odd',status='replace')
    Endif
    open(1003,file='nt_Even.txt',status='replace')
    write(1003,*)nt
    close(1003)
    nt_Even=nt
  endif

  write(1000,99) real(F_kmu)
  write(1001,99) aimag(F_kmu)

  close(1000)
  close(1001)

  IF(Add_GAM==1)then
  write(1004,99) real(F_GAMrec)
  write(1005,99) aimag(F_GAMrec)
  close(1004)
  close(1005)
  Endif

  !---------------------------------------------------------
  close(stdout+20)
  close(stdout+21)
  close(stdout+22)
  close(stdout+36)
  close(stdout+37)
  close(stdout+55)
  close(stdout+63) 
  open(stdout+20,file="total_D_3G.txt",status='old',POSITION='APPEND')
  open(stdout+21,file="Phi_k.txt",status='old',POSITION='APPEND')
  open(stdout+22,file="Chi_mu.txt",status='old',POSITION='APPEND')
  open(stdout+36,file="Chi.txt",status='old',POSITION='APPEND')
  open(stdout+37,file="total_Chi.txt",status='old',POSITION='APPEND')
  open(stdout+55,file="sqrn_n0.txt",status='old',POSITION='APPEND')
  open(stdout+63,file="E.txt",status='old',POSITION='APPEND') 
  !+++++++++++++++++++++++++++++++++++++++++++++++
  IF(verify_NL==1)then
  close(stdout+15)
  close(stdout+13)
  close(stdout+16)
  close(stdout+17)
  close(stdout+18)
  close(stdout+19)
  close(stdout+23)
  close(stdout+24)
  close(stdout+25)
  close(stdout+31)
  close(stdout+38)  
  close(stdout+32)  
  IF(Add_GAM==1)then
  close(stdout+26)
  close(stdout+27)
  close(stdout+28)
  Endif
  open(stdout+15,file="energy_kx_ky.txt",status='old',POSITION='APPEND')
  open(stdout+13,file="entropy_kx_ky.txt",status='old',POSITION='APPEND')
  open(stdout+16,file="D.txt",status='old',POSITION='APPEND')
  open(stdout+17,file="total_D.txt",status='old',POSITION='APPEND')
  open(stdout+18,file="Fk_mu.txt",status='old',POSITION='APPEND')
  open(stdout+19,file="D_3G.txt",status='old',POSITION='APPEND')
  open(stdout+23,file="tot_Entropy.txt",status='old',POSITION='APPEND')
  open(stdout+24,file="total_energy.txt",status='old',POSITION='APPEND')
  open(stdout+25,file="energy_k.txt",status='old',POSITION='APPEND')
  open(stdout+31,file="NLconser.txt",status='old',POSITION='APPEND')
  open(stdout+38,file="Phisqr.txt",status='old',POSITION='APPEND')
  open(stdout+32,file="ne_k.txt",status='old',POSITION='APPEND')
  IF(Add_GAM==1)then
  open(stdout+26,file="tot_ene_GAM.txt",status='old',POSITION='APPEND')
  open(stdout+27,file="ene_GAM.txt",status='old',POSITION='APPEND')
  open(stdout+28,file="F_GAM_1.txt",status='old',POSITION='APPEND')
  Endif
  Endif

  IF(verify_NL==2)then
  close(stdout+58)
  close(stdout+59)
  open(stdout+58,file="GOMG.txt",status='old',POSITION='APPEND')
  open(stdout+59,file="Phi_kOMG.txt",status='old',POSITION='APPEND')
  Endif
  !+++++++++++++++++++++++++++++++++++++++++++++++

 Endif

Endif

!=================================================================================================================================
!CKinFH
!4D Fourier harmonic
IF(GK_FK_CK==1)then

 If(myid==0)then
  if(verify_L==1)then
  call cpu_time(DBtb)
  endif

  if(nt_Odd==min(nt_Odd,nt_Even))then
    open(1008,file='Fkft_re_Odd',status='replace')
    open(1009,file='Fkft_im_Odd',status='replace')
    open(1003,file='nt_Odd.txt',status='replace')
    write(1003,*)nt
    close(1003)
    nt_Odd=nt
  elseif(nt_Even==min(nt_Odd,nt_Even))then
    open(1008,file='Fkft_re_Even',status='replace')
    open(1009,file='Fkft_im_Even',status='replace')
    open(1003,file='nt_Even.txt',status='replace')
    write(1003,*)nt
    close(1003)
    nt_Even=nt
  endif

  write(1008,*) real(F_kfta)
  write(1009,*) aimag(F_kfta)
  close(1008)
  close(1009)

  !---------------------------------------------------------
  close(stdout+43)
  close(stdout+45)
  close(stdout+46)
  close(stdout+47)
  close(stdout+55)
  close(stdout+63) 
  open(stdout+43,file="Phi_kft.txt",status='old',POSITION='APPEND')
  open(stdout+45,file="total_Dft.txt",status='old',POSITION='APPEND')
  open(stdout+46,file="Chi_kft.txt",status='old',POSITION='APPEND')
  open(stdout+47,file="total_Chift.txt",status='old',POSITION='APPEND')
  open(stdout+55,file="sqrn_n0ft.txt",status='old',POSITION='APPEND')
  open(stdout+63,file="Eft.txt",status='old',POSITION='APPEND') 
  !+++++++++++++++++++++++++++++++++++++++++++++++
  IF(verify_NL==1)then
  close(stdout+38)  
  close(stdout+44)
  close(stdout+48)
  close(stdout+49)
  close(stdout+51)
  close(stdout+52)
  close(stdout+74)
  open(stdout+38,file="Phisqrft.txt",status='old',POSITION='APPEND')
  open(stdout+44,file="D_kft.txt",status='old',POSITION='APPEND')
  open(stdout+48,file="tot_Entropyft.txt",status='old',POSITION='APPEND')
  open(stdout+49,file="Chi_muft.txt",status='old',POSITION='APPEND')
  open(stdout+51,file="Fk_ftmu.txt",status='old',POSITION='APPEND')
  open(stdout+52,file="Fk_ftmuT.txt",status='old',POSITION='APPEND')
  open(stdout+74,file="NLconserft.txt",status='old',POSITION='APPEND')
  Endif

  IF(verify_NL==2)then
  close(stdout+58)
  close(stdout+59)
  open(stdout+58,file="GOMG.txt",status='old',POSITION='APPEND')
  open(stdout+59,file="Phi_kOMG.txt",status='old',POSITION='APPEND')
  Endif

  if(verify_L==1)then
  call cpu_time(DBte)
  DBt=DBt+ (DBte-DBtb)
  endif

 Endif
  !---------------------------------------------------------

Endif
!=================================================================================================================================
!CKinCH
!cyclo-kinetic
IF(GK_FK_CK==2 .or. GK_FK_CK==- 2)then

 If(myid==0)then
  if(verify_L==1)then
  call cpu_time(DBtb)
  endif

  if(nt_Odd==min(nt_Odd,nt_Even))then
    open(1006,file='Fkcy_re_Odd',status='replace')
    open(1007,file='Fkcy_im_Odd',status='replace')
    open(1003,file='nt_Odd.txt',status='replace')
    write(1003,*)nt
    close(1003)
    nt_Odd=nt
  elseif(nt_Even==min(nt_Odd,nt_Even))then
    open(1006,file='Fkcy_re_Even',status='replace')
    open(1007,file='Fkcy_im_Even',status='replace')
    open(1003,file='nt_Even.txt',status='replace')
    write(1003,*)nt
    close(1003)
    nt_Even=nt
  endif

  write(1006,99) real(F_kcya)
  write(1007,99) aimag(F_kcya)
  close(1006)
  close(1007)

  !---------------------------------------------------------
  close(stdout+34)
  close(stdout+51)
  close(stdout+46)
  close(stdout+53)
  close(stdout+55)
  close(stdout+62)
  close(stdout+63)    
  open(stdout+34,file="Phi_kcy.txt",status='old',POSITION='APPEND')
  open(stdout+51,file="total_Dcy.txt",status='old',POSITION='APPEND')
  open(stdout+46,file="Chi_kcy.txt",status='old',POSITION='APPEND')
  open(stdout+53,file="total_Chicy.txt",status='old',POSITION='APPEND')
  open(stdout+55,file="sqrn_n0cy.txt",status='old',POSITION='APPEND') 
  open(stdout+62,file="ni_n0.txt",status='old',POSITION='APPEND')
  open(stdout+63,file="Ecy.txt",status='old',POSITION='APPEND')  
  !+++++++++++++++++++++++++++++++++++++++++++++++
  IF(verify_NL==1)then
  close(stdout+38)  
  close(stdout+44)
  close(stdout+49)
  close(stdout+52)
  close(stdout+54)
  close(stdout+56)
  close(stdout+74)
  open(stdout+38,file="Phisqrcy.txt",status='old',POSITION='APPEND')
  open(stdout+44,file="D_kcy.txt",status='old',POSITION='APPEND')
  open(stdout+49,file="Chi_mucy.txt",status='old',POSITION='APPEND')
  open(stdout+52,file="Fk_cymuT.txt",status='old',POSITION='APPEND')
  open(stdout+54,file="tot_Entropycy.txt",status='old',POSITION='APPEND')
  open(stdout+56,file="ne.txt",status='old',POSITION='APPEND') 
  open(stdout+74,file="NLconsercy.txt",status='old',POSITION='APPEND')
  Endif

  IF(verify_NL==2)then
  close(stdout+58)
  close(stdout+59)
  open(stdout+58,file="GOMG.txt",status='old',POSITION='APPEND')
  open(stdout+59,file="Phi_kOMG.txt",status='old',POSITION='APPEND')
  Endif

  if(verify_L==1)then
  call cpu_time(DBte)
  DBt=DBt+ (DBte-DBtb)
  endif

 Endif

Endif
!=================================================================================================================================

  If(myid==0)then
  write(*,*)'make a restart point at',' nt=',nt,'   time=',nt*tstep
  Endif

  99 format(e14.7)

end Subroutine backup


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine diagnose(myid)

  use mpi
  use parameter
  implicit none
  integer,intent(in) :: myid
  real readtemp
  real CNaN
  character(20)::str
  real a_Lntmp, a_LTitmp
  !GK
  real tot_Entropy
  real total_energy
  real Phisqr
  real tot_ene_GAM
  complex F_GAM_1,readtempcom
  !CKinFH
  real tot_Entropyft
  real total_Dft
  real total_Chift
  !CKinCH
  real tot_Entropycy
  real total_Dcy
  real total_Chicy
  real ni_n0,Phisks

    IF(a_Ln==0)then
    a_Lntmp=1
    else
    a_Lntmp=a_Ln
    endif

    IF(a_LTi==0)then
    a_LTitmp=1
    else
    a_LTitmp=a_LTi
    endif

    n_n0com=(0,0)
    n_n0=0
    Phisqr=0

!=================================================================================================================================
!GK
IF(GK_FK_CK==0)then
    total_D=0
    total_D_3G=0
    total_Chi=0
    D_kx_ky(:,:)=0
    D_3G_kx_ky(:,:)=0
    Chi_kx_ky(:,:)=0
    tot_Entropy=0
    total_energy=0
    energy_kx_ky(:,:)=0
    entropy_kx_ky(:,:)=0
    Chi_mu(:)=0
    Fk_mucom(:)=(0,0)
    Fk_mu(:)=0
    kk=0
    energy_kk(:)=0
    IF(add_GAM==1)then
    tot_ene_GAM=0
    ene_GAM(:)=0
    F_GAM_1=(0,0)
    Endif
    IF(verify_NL == 2)then
    allocate(GOMG(-pmax+1:pmax-1,-nmax+1:nmax-1))
    GOMG(:,:)=(0,0)
    Endif

  If(myid==0)then
     do n=0,nmax-1
     do p=-pmax+1,pmax-1
        do  i=1,N_mu    !gauss legendre integration
        Integ_num(p,n)=Integ_num(p,n)+w_point(i)*Jn(p,n,i,0)*F_kmu(p,n,i)
        enddo
        D_3G_kx_ky(p,n)=real((0,1.0)*ky(n)*conjg(phi_k(p,n))*Integ_num(p,n)/a_Lntmp)
        total_D_3G=total_D_3G + D_3G_kx_ky(p,n)

        Integ_num(p,n)=(0,0)
        do  i=1,N_mu    !gauss legendre integration
        Integ_num(p,n)=Integ_num(p,n)+w_point(i)*Jn(p,n,i,0)*F_kmu(p,n,i)*mu_point(i)
        enddo
        Chi_kx_ky(p,n)=real((0,1.0)*ky(n)*conjg(phi_k(p,n))*Integ_num(p,n)/a_LTitmp)

        if(n==0)then
        do  i=1,N_mu
        Chi_mu(i)=Chi_mu(i)+real(w_point(i)*Jn(p,n,i,0)*F_kmu(p,n,i)*mu_point(i)*(0,1.0)*ky(n)*conjg(phi_k(p,n))/a_LTitmp)
        enddo
        elseif (n>0)then
        do  i=1,N_mu
        Chi_mu(i)=Chi_mu(i)+2*real(w_point(i)*Jn(p,n,i,0)*F_kmu(p,n,i)*mu_point(i)*(0,1.0)*ky(n)*conjg(phi_k(p,n))/a_LTitmp)
        enddo
        endif

       IF(verify_NL == 2)then
        Integ_num(p,n)=(0,0)
        do  i=1,N_mu    !gauss legendre integration
        Integ_num(p,n)=Integ_num(p,n)+w_point(i)*Jn(p,n,i,0)*F_kmu(p,n,i)*mu_point(i)
        enddo
        GOMG(p,n)=Integ_num(p,n)
       Endif
     enddo
     enddo
    total_Chi=total_Chi + sum(Chi_mu(:))


     do n=-nmax+1,nmax-1
     do p=-pmax+1,pmax-1
        n_n0com=n_n0com + F_kmu(p,n,0)*conjg(F_kmu(p,n,0))/(Omega_star**2)
     enddo
     enddo
     n_n0=abs(n_n0com)

     Phisks=0
     do p=-pmax+1,pmax-1
      do n=-nmax+1,nmax-1
      Phisks=Phisks+(kx(p)**2+ky(n)**2)*conjg(phi_k(p,n))*phi_k(p,n)/(Omega_star**2)/2.0
      enddo
    enddo

    write(str,*)total_Chi
    if(Trim(ADJUSTL(str)) .eq. 'NaN')then
    write(*,*)"The total_Chi goes to NaN !!!"
    stop
    endif

     do n=1,nmax-1
     do p=-pmax+1,pmax-1
       D_3G_kx_ky(-p,-n)=D_3G_kx_ky(p,n)
       total_D_3G=total_D_3G+D_3G_kx_ky(-p,-n)
       Chi_kx_ky(-p,-n)=Chi_kx_ky(p,n)

       IF(verify_NL == 2)then
       GOMG(-p,-n)=GOMG(p,n)
       ENDIF
     enddo
     enddo


! output file of energy, diffusive coefficient etc.
    if(dir_dia_count==1.and.restart==0)then  !open new output files at the first time in a run
      open(stdout+20,file="total_D_3G.txt",status='replace')
      open(stdout+21,file="Phi_k.txt",status='replace')
      open(stdout+22,file="Chi_mu.txt",status='replace')
      open(stdout+36,file="Chi.txt",status='replace')
      open(stdout+37,file="total_Chi.txt",status='replace')
      open(stdout+55,file="sqrn_n0.txt",status='replace')
      open(stdout+63,file="E.txt",status='replace')  
    elseif(dir_dia_count==1.and.(restart==1.or.restart==2))then !open old output files at the beginning of the restart
      open(stdout+20,file="total_D_3G.txt",status='old')
      open(stdout+21,file="Phi_k.txt",status='old')
      open(stdout+22,file="Chi_mu.txt",status='old')
      open(stdout+36,file="Chi.txt",status='old')
      open(stdout+37,file="total_Chi.txt",status='old')
      open(stdout+55,file="sqrn_n0.txt",status='old')
      open(stdout+63,file="E.txt",status='old')  
    endif

    If(dir_dia_count==1.and.(restart==1.or.restart==2))then
      do i=1,nt/output_step
        read(stdout+20,*)readtemp
      enddo
     do i=1,(2*nmax-1)*(2*pmax-1)*nt/output_step
        read(stdout+21,*)readtempcom
     enddo
     do i=1,N_mu*nt/output_step
        read(stdout+22,*)readtemp
     enddo
      do i=1,(2*nmax-1)*(2*pmax-1)*nt/output_step
        read(stdout+36,*)readtemp
      enddo
      do i=1,nt/output_step
        read(stdout+37,*)readtemp
      enddo
      do i=1,nt/output_step
        read(stdout+55,*)readtemp
      enddo
      do i=1,nt/output_step
        read(stdout+63,*)readtemp
      enddo          
    Endif

      write(stdout+20,102) total_D_3G
      write(stdout+36,102) Chi_kx_ky
      write(stdout+37,102) total_Chi
      write(stdout+22,102) Chi_mu
      write(stdout+55,102) n_n0
      write(stdout+63,102) Phisks

      do n=-nmax+1,nmax-1
       do p=-pmax+1,pmax-1
       write(stdout+21,"('(',e14.7,',',e14.7,')')")real(Phi_k(p,n)),aimag(Phi_k(p,n))
       enddo
      enddo

  !+++++++++++++++++++++++++++++++++++++++++++++++
  IF(verify_NL==1)then
    IF(Add_GAM==1)then
    call MPI_GATHER(F_GAM(:),(2*pmax-1),MPI_COMPLEX,&
                  F_GAMdia(:,:),(2*pmax-1),MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)
    Endif

     do n=0,nmax-1
     do p=-pmax+1,pmax-1
        do  i=1,N_mu
        energy_kx_ky(p,n)=energy_kx_ky(p,n)+ w_point(i)*(F_kmu(p,n,i)*conjg(phi_k(p,n))+&
                              conjg(F_kmu(p,n,i))*phi_k(p,n))*Jn(p,n,i,0)/2.0
        entropy_kx_ky(p,n)=entropy_kx_ky(p,n)+ w_point(i)*conjg(F_kmu(p,n,i))*F_kmu(p,n,i)
        enddo
!        i=0
!        entropy_kx_ky(p,n)=entropy_kx_ky(p,n)+conjg(F_kmu(p,n,i))*F_kmu(p,n,i)*CDW

        total_energy=total_energy+energy_kx_ky(p,n)
        tot_Entropy=tot_Entropy+entropy_kx_ky(p,n)
        D_kx_ky(p,n)=ky(n)*delta_k(p,n)*conjg(phi_k(p,n))*phi_k(p,n)/a_Lntmp
        total_D=total_D + D_kx_ky(p,n)
     enddo
     enddo

     do n=1,nmax-1
     do p=-pmax+1,pmax-1
       energy_kx_ky(-p,-n)=energy_kx_ky(p,n)
       total_energy=total_energy+energy_kx_ky(-p,-n)

       D_kx_ky(-p,-n)=D_kx_ky(p,n)
       total_D=total_D+D_kx_ky(-p,-n)

       entropy_kx_ky(-p,-n)=entropy_kx_ky(p,n)
       tot_Entropy=tot_Entropy+entropy_kx_ky(-p,-n)
     enddo
     enddo

     do n=0,nmax-1
     do p=-pmax+1,pmax-1
        do  i=1,N_mu
        Fk_mucom(i)=Fk_mucom(i)+F_kmu(p,n,i)
        enddo
     enddo
     enddo

     do n=1,nmax-1
     do p=-pmax+1,pmax-1
       do  i=1,N_mu
       Fk_mucom(i)=Fk_mucom(i)+conjg(F_kmu(p,n,i))
       enddo
     enddo
     enddo
     Fk_mu(:)=2*PI*abs(Fk_mucom(:))

     do p=-pmax+1,pmax-1
      do n=-nmax+1,nmax-1
      Phisqr=Phisqr+conjg(phi_k(p,n))*phi_k(p,n)/(Omega_star**2)
      enddo
    enddo
    
  IF(Add_GAM==1)then
  do p=-pmax+1,pmax-1
    do  i=1,N_mu
    ene_GAM(p)=ene_GAM(p)+w_point(i)*(F_GAMdia(p,i)*conjg(phi_k(p,0))+&
                             conjg(F_GAMdia(p,i))*phi_k(p,0))*Jn(p,0,i,0)/2.0
    enddo
   ! ene_GAM(p)=ene_GAM(p)-conjg(phi_k(p,0))*phi_k(p,0)*Integ_den(p,0)/ratio_TiTe
    tot_ene_GAM=tot_ene_GAM+ ene_GAM(p)
  enddo

  total_energy=total_energy+tot_ene_GAM
  Endif

!others
    do p=-pmax+1,pmax-1
      do n=-nmax+1,nmax-1
      kk=p**2+n**2
      if(kk/=0)then
      energy_kk(kk)=energy_kk(kk)+energy_kx_ky(p,n)
      endif
      enddo
    enddo

    i=0
    do kk=1,(pmax-1)**2+(nmax-1)**2
    if(kk_ngrid(kk)/=0)then
    i=i+1
    energy_kk(kk)=(energy_kk(kk)/kk_ngrid(kk))*2*Pi*k(i)
    endif
    enddo

    i=0
    do kk=1,(pmax-1)**2+(nmax-1)**2
      if(kk_ngrid(kk)/=0)then
      i=i+1
      energy_k(i)=energy_kk(kk)
      endif
    enddo

    if(dir_dia_count==1.and.restart==0)then
      open(stdout+15,file="energy_kx_ky.txt",status='replace')
      open(stdout+13,file="entropy_kx_ky.txt",status='replace')
      open(stdout+16,file="D.txt",status='replace')
      open(stdout+18,file="Fk_mu.txt",status='replace')
      open(stdout+19,file="D_3G.txt",status='replace')
      open(stdout+17,file="total_D.txt",status='replace')
      open(stdout+23,file="tot_Entropy.txt",status='replace')
      open(stdout+24,file="total_energy.txt",status='replace')
      open(stdout+25,file="energy_k.txt",status='replace')
      open(stdout+38,file="Phisqr.txt",status='replace')
      open(stdout+32,file="ne_k.txt",status='replace')
      IF(Add_GAM==1)then
      open(stdout+26,file="tot_ene_GAM.txt",status='replace')
      open(stdout+27,file="ene_GAM.txt",status='replace')
      open(stdout+28,file="F_GAM_1.txt",status='replace')
      Endif
    elseif(dir_dia_count==1.and.(restart==1.or.restart==2))then
      open(stdout+15,file="energy_kx_ky.txt",status='old')
      open(stdout+13,file="entropy_kx_ky.txt",status='old')
      open(stdout+16,file="D.txt",status='old')
      open(stdout+17,file="total_D.txt",status='old')
      open(stdout+18,file="Fk_mu.txt",status='old')
      open(stdout+19,file="D_3G.txt",status='old')
      open(stdout+23,file="tot_Entropy.txt",status='old')
      open(stdout+24,file="total_energy.txt",status='old')
      open(stdout+25,file="energy_k.txt",status='old')
      open(stdout+38,file="Phisqr.txt",status='old')
      open(stdout+32,file="ne_k.txt",status='old')
      IF(Add_GAM==1)then
      open(stdout+26,file="tot_ene_GAM.txt",status='old')
      open(stdout+27,file="ene_GAM.txt",status='old')
      open(stdout+28,file="F_GAM_1.txt",status='old')
      !open(stdout+29,file="F_GAM_1_Amp.txt",status='old')
      Endif
    endif

    If(dir_dia_count==1.and.(restart==1.or.restart==2))then
      do i=1,(2*nmax-1)*(2*pmax-1)*nt/output_step
        read(stdout+15,*)readtemp
      enddo
      do i=1,(2*nmax-1)*(2*pmax-1)*nt/output_step
        read(stdout+13,*)readtemp
      enddo
      do i=1,(2*nmax-1)*(2*pmax-1)*nt/output_step
        read(stdout+16,*)readtemp
      enddo
      do i=1,nt/output_step
        read(stdout+17,*)readtemp
      enddo
     do i=1,N_mu*nt/output_step
        read(stdout+18,*)readtemp
     enddo
      do i=1,(2*nmax-1)*(2*pmax-1)*nt/output_step
        read(stdout+19,*)readtemp
      enddo
      do i=1,nt/output_step
        read(stdout+23,*)readtemp
      enddo

      do i=1,nt/output_step
        read(stdout+24,*)readtemp
      enddo
      do i=1,kmax*nt/output_step
        read(stdout+25,*)readtemp
      enddo

      do i=1,nt/output_step
        read(stdout+38,*)readtemp
      enddo
      
     do i=1,(2*nmax-1)*(2*pmax-1)*nt/output_step
        read(stdout+32,*)readtempcom
     enddo
      
      IF(Add_GAM==1)then
      do i=1,nt/output_step
        read(stdout+26,*)readtemp
      enddo
      do i=1,(2*pmax-1)*nt/output_step
        read(stdout+27,*)readtemp
      enddo
      do i=1,3*nt/output_step
        read(stdout+28,*)readtemp
      enddo
      Endif

    Endif

      write(stdout+15,102) energy_kx_ky
      write(stdout+13,102) entropy_kx_ky
      write(stdout+16,102) D_kx_ky
      write(stdout+17,102) total_D
      write(stdout+18,102) Fk_mu
      write(stdout+19,102) D_3G_kx_ky
      write(stdout+23,102) tot_Entropy
      write(stdout+24,102) total_energy
      write(stdout+25,102) energy_k
      write(stdout+38,102) Phisqr
      IF(Add_GAM==1)then
      write(stdout+26,102) tot_ene_GAM
      write(stdout+27,102) ene_GAM
      write(stdout+28,102) real(F_GAM_1)
      write(stdout+28,102) aimag(F_GAM_1)
      write(stdout+28,102) CABS(F_GAM_1)
      Endif

      do n=-nmax+1,nmax-1
       do p=-pmax+1,pmax-1
       write(stdout+32,"('(',e14.7,',',e14.7,')')")real(F_kmu(p,n,0)),aimag(F_kmu(p,n,0))
       enddo
      enddo
      
  Endif

  !-----------------------------------------------
  IF(verify_NL==2)then
    if(dir_dia_count==1.and.restart==0)then
    open(stdout+58,file="GOMG.txt",status='replace')
    open(stdout+59,file="Phi_kOMG.txt",status='replace')
    endif

    If(dir_dia_count==1.and.(restart==1.or.restart==2))then
     If(Stopnt==0)then
      open(stdout+58,file="GOMG.txt",status='old')
      open(stdout+59,file="Phi_kOMG.txt",status='old')
      do i=1,(2*nmax-1)*(2*pmax-1)*(nt)/output_step
        read(stdout+58,*)readtempcom
      enddo
      do i=1,(2*nmax-1)*(2*pmax-1)*(nt)/output_step
        read(stdout+59,*)readtempcom
      enddo
     ELSE
       if(nt==Stopnt .or. nt==Stopnt+1)then
        open(stdout+58,file="GOMG.txt",status='replace')
        open(stdout+59,file="Phi_kOMG.txt",status='replace')
       else
        open(stdout+58,file="GOMG.txt",status='old')
        open(stdout+59,file="Phi_kOMG.txt",status='old')

        do i=1,(2*nmax-1)*(2*pmax-1)*(nt-Stopnt)/output_step
        read(stdout+58,*)readtempcom
        enddo
        do i=1,(2*nmax-1)*(2*pmax-1)*(nt-Stopnt)/output_step
        read(stdout+59,*)readtempcom
        enddo
       endif
     ENDIF
    Endif


      do n=-nmax+1,nmax-1
       do p=-pmax+1,pmax-1
       write(stdout+58,"('(',e14.7,',',e14.7,')')")real(GOMG(p,n)),aimag(GOMG(p,n))
       enddo
      enddo


      do n=-nmax+1,nmax-1
       do p=-pmax+1,pmax-1
      write(stdout+59,"('(',e14.7,',',e14.7,')')")real(Phi_k(p,n)),aimag(Phi_k(p,n))
       enddo
      enddo

      deallocate(GOMG)
  Endif

  !+++++++++++++++++++++++++++++++++++++++++++++++

  if(dir_dia_count==1)then   !write the directory name of the output data at the first time in a run of this 3D code, and calculate linear disperson relation Omega

  open(stdout+990,file="data_file_name.txt",status='replace')
  open(stdout+991,file="cpu_timeGK.txt",status='old',POSITION='APPEND')

  If(Add_GAM==0)then
  write(stdout+990,&
"('iso',i2.2,'NL',f4.1,'nt',i6.6,'tsp',f9.6,'pmax',i2.2,'delta',f6.3,'mu',f6.3,'lambda',f5.2,'Omega_star',f7.2,'N_CY',i1)")&
 isotropic,Const_NL,ntmax,tstep,pmax,delta_1,mu_HK,lambda_0,Omega_star,N_CY!,F_k_int*((2*pmax-1)*(2*nmax-1)*(N_mu))
    write(stdout+991,*)"mu_point="
    do i=1,N_mu
    write(stdout+991,*)mu_point(i)
    enddo
    write(stdout+991,*)" "
    write(stdout+991,*)"w_point="
    do i=1,N_mu
    write(stdout+991,*)w_point(i)
    enddo
    write(stdout+991,*)" "
    write(stdout+991,*)"w_pointMod="
    do i=1,N_mu
    write(stdout+991,*)w_point(i)*FM(i)
    enddo
    write(stdout+991,*)" "

  FM_integ=0
  do i=1,N_mu
     FM_integ = FM_integ+w_point(i)*FM(i)
  enddo
  write(stdout+991,*)'Sum_Wmod=',FM_integ
  write(stdout+991,*)" "
  FM_integ=0
  do i=1,N_mu
     FM_integ = FM_integ+w_point(i)*mu_point(i)*exp(-mu_point(i)**2)
  enddo
  write(stdout+991,*)'Test1=',FM_integ
  write(stdout+991,*)' '
  FM_integ=0
  do i=1,N_mu
     FM_integ = FM_integ+w_point(i)*mu_point(i)*exp(-mu_point(i))
  enddo
  write(stdout+991,*)'Test2=',FM_integ
  write(stdout+991,*)' '
  write(stdout+991,*)"GK"
  write(stdout+991,"('CDW=',I2)")CDW
  write(stdout+991,"('mugridtype=',I2)")mugridtype
  if(mugridtype==3)then
  write(stdout+991,"('MDmu12=',f5.2)")MDmu12
  endif
  write(stdout+991,*)"Add_GAM=",Add_GAM
  write(stdout+991,*)"ntmax=",ntmax
  write(stdout+991,"('tstep=',f8.5)")tstep
  write(stdout+991,"('output_step=',I5)")output_step
  write(stdout+991,"('pmax=',I2)")pmax
  write(stdout+991,"('kxmax=',f6.3)")kxmax
  write(stdout+991,"('mumax=',f4.1)")mumax
  write(stdout+991,"('N_mu=',I2)")N_mu
  write(stdout+991,"('a_Ln=',f6.3)")a_Ln
  write(stdout+991,"('a_LTi=',f6.3)")a_LTi
  write(stdout+991,"('lambda_0=',f6.3)")lambda_0
  write(stdout+991,"('OMG=',f5.1)")Omega_star
  If(CDW==1)then
  write(stdout+991,"('lambda_D=',f6.3)")lambda_D
  write(stdout+991,"('AlphaA=',f8.3)")AlphaA
  Endif
  write(stdout+991,"('delta_1=',f6.3)")delta_1
  write(stdout+991,"('mu_HK=',f6.3)")mu_HK
  write(stdout+991,"('mu_LK=',f7.4)")mu_LK
  write(stdout+991,"('nu_DW=',f6.3)")nu_DW
  write(stdout+991,"('nu_ZF=',f6.3)")nu_ZF
  write(stdout+991,"('G=',f6.3)")G
  write(stdout+991,"('Epsilon=',f6.3)")Epsilon
  write(stdout+991,"('kpar=',f7.3)")kpar
  write(stdout+991,"('uD=',f8.4)")uD
  write(stdout+991,"('CDWid=',f5.2)")CDWid
  write(stdout+991,"('gamE=',f8.4)")gamE
  write(stdout+991,"('Beta=',f6.2)")Beta
  Else if (Add_GAM==1)then
  write(stdout+990,&
"('nt',i6.6,'tsp',f8.5,'pmax',i2.2,'delta',f6.3,'muCD',f6.3,'nuDW',f7.4,'nuZF',f7.4,'lambda',f5.2,'Omega_st',f7.4,'N_CY',i1)")&!,'Int',f8.4)")&
 ntmax,tstep,pmax,delta_1,mu_HK,nu_DW,nu_ZF,lambda_0,Omega_star,N_CY
    write(stdout+991,*)"mu_point="
    do i=1,N_mu
    write(stdout+991,*)mu_point(i)
    enddo
    write(stdout+991,*)" "
    write(stdout+991,*)"w_point="
    do i=1,N_mu
    write(stdout+991,*)w_point(i)
    enddo
    write(stdout+991,*)" "
    write(stdout+991,*)"w_pointMod="
    do i=1,N_mu
    write(stdout+991,*)w_point(i)*FM(i)
    enddo
    write(stdout+991,*)" "

  FM_integ=0
  do i=1,N_mu
     FM_integ = FM_integ+w_point(i)*FM(i)
  enddo
  write(stdout+991,*)'Sum_Wmod=',FM_integ
  write(stdout+991,*)" "
  FM_integ=0
  do i=1,N_mu
     FM_integ = FM_integ+w_point(i)*mu_point(i)*exp(-mu_point(i)**2)
  enddo
  write(stdout+991,*)'Test1=',FM_integ
  write(stdout+991,*)' '
  FM_integ=0
  do i=1,N_mu
     FM_integ = FM_integ+w_point(i)*mu_point(i)*exp(-mu_point(i))
  enddo
  write(stdout+991,*)'Test2=',FM_integ
  write(stdout+991,*)' '
  write(stdout+991,*)"GK"
  write(stdout+991,"('CDW=',I2)")CDW
  write(stdout+991,"('mugridtype=',I2)")mugridtype
  if(mugridtype==3)then
  write(stdout+991,"('MDmu12=',f5.2)")MDmu12
  endif
  write(stdout+991,*)"Add_GAM=",Add_GAM
  write(stdout+991,*)"ntmax=",ntmax
  write(stdout+991,"('tstep=',f8.5)")tstep
  write(stdout+991,"('output_step=',I5)")output_step
  write(stdout+991,"('pmax=',I2)")pmax
  write(stdout+991,"('kxmax=',f6.3)")kxmax
  write(stdout+991,"('mumax=',f4.1)")mumax
  write(stdout+991,"('N_mu=',I2)")N_mu
  write(stdout+991,"('a_Ln=',f6.3)")a_Ln
  write(stdout+991,"('a_LTi=',f6.3)")a_LTi
  write(stdout+991,"('lambda_0=',f6.3)")lambda_0
  write(stdout+991,"('OMG=',f5.1)")Omega_star
  If(CDW==1)then
  write(stdout+991,"('lambda_D=',f6.3)")lambda_D
  write(stdout+991,"('AlphaA=',f8.3)")AlphaA
  Endif
  write(stdout+991,"('delta_1=',f6.3)")delta_1
  write(stdout+991,"('mu_HK=',f6.3)")mu_HK
  write(stdout+991,"('mu_LK=',f7.4)")mu_LK
  write(stdout+991,"('nu_DW=',f6.3)")nu_DW
  write(stdout+991,"('nu_ZF=',f6.3)")nu_ZF
  write(stdout+991,"('G=',f6.3)")G
  write(stdout+991,"('Epsilon=',f6.3)")Epsilon
  write(stdout+991,"('uD=',f8.4)")uD
  write(stdout+991,"('CDWid=',f5.2)")CDWid
  write(stdout+991,"('gamE=',f8.4)")gamE
  write(stdout+991,"('Beta=',f6.2)")Beta
  write(stdout+991,"('CDWid=',f6.3)")CDWid
  Endif
  close(stdout+990)
  close(stdout+991)
  endif

 Endif

Endif

!=================================================================================================================================
!CKinFH
IF(GK_FK_CK==1)then
!Fourier kinetics
 If(myid==0)then

  if(verify_L==1)then
  call cpu_time(DBtb)
  endif

    tot_Entropyft=0
    total_Dft=0
    total_Chift=0
    D_kft=0
    Chi_kft=0
    Chi_muft=0
    Fk_ftmucom(:)=(0,0)
    Fk_ftmu(:)=0
    Fk_ftcomT(:)=(0,0)
    Fk_ftmuT(:)=0
    IF(verify_NL == 2)then
    allocate(GOMG(-pmax+1:pmax-1,-nmax+1:nmax-1))
    GOMG(:,:)=(0,0)
    Endif

   do n=0,nmax-1
   do p=-pmax+1,pmax-1
        Integ_num(p,n)=(0,0)  !set the integral variables to zero before integral
        do  i=1,N_mu    !gauss legendre integration
        a=0
        l=i*(2*N_FT-1)+a+(N_FT-1)
        Integ_num(p,n)=Integ_num(p,n)+w_point(i)*(F_kfta(p,n,l)+Phi_kft(p,n)*FM(i)/ratio_TiTe)
        enddo
        D_kft(p,n)=real((0,1.0)*ky(n)*conjg(Phi_kft(p,n))*Integ_num(p,n)/a_Lntmp)
        total_Dft=total_Dft+D_kft(p,n)

        Integ_num(p,n)=(0,0)  !set the integral variables to zero before integral
        do  i=1,N_mu    !gauss legendre integration
        a=0
        l=i*(2*N_FT-1)+a+(N_FT-1)
        Integ_num(p,n)=Integ_num(p,n)+w_point(i)*ratio_TiTe*mu_point(i)*(F_kfta(p,n,l)+Phi_kft(p,n)*FM(i)/ratio_TiTe)
        enddo
        Chi_kft(p,n)=real((0,1.0)*ky(n)*conjg(Phi_kft(p,n))*Integ_num(p,n)/a_LTitmp)

        IF(n==0)then
        do  i=1,N_mu
        a=0
        l=i*(2*N_FT-1)+a+(N_FT-1)
        Chi_muft(i)=Chi_muft(i)+real((0,1.0)*ky(n)*conjg(Phi_kft(p,n))*w_point(i)*ratio_TiTe*mu_point(i)*(F_kfta(p,n,l)&
        +Phi_kft(p,n)*FM(i)/ratio_TiTe))/a_LTitmp
        enddo
        elseif (n>0)then
        do  i=1,N_mu
        a=0
        l=i*(2*N_FT-1)+a+(N_FT-1)
        Chi_muft(i)=Chi_muft(i)+2*real((0,1.0)*ky(n)*conjg(Phi_kft(p,n))*w_point(i)*ratio_TiTe*mu_point(i)*(F_kfta(p,n,l)&
                                     +Phi_kft(p,n)*FM(i)/ratio_TiTe))/a_LTitmp
        enddo
        endif

       IF(verify_NL == 2)then
        Integ_num(p,n)=(0,0)  !set the integral variables to zero before integral
        do  i=1,N_mu    !gauss legendre integration
        a=0
        l=i*(2*N_FT-1)+a+(N_FT-1)
        Integ_num(p,n)=Integ_num(p,n)+w_point(i)*mu_point(i)*(F_kfta(p,n,l)+Phi_kft(p,n)*FM(i)/ratio_TiTe)
        enddo
        GOMG(p,n)=Integ_num(p,n)
       Endif
     enddo
     enddo
     total_Chift=sum(Chi_muft(:))


     do n=-nmax+1,nmax-1
     do p=-pmax+1,pmax-1
       i=0
       a=0
       l=i*(2*N_FT-1)+a+(N_FT-1)
       n_n0com=n_n0com + F_kfta(p,n,l)*conjg(F_kfta(p,n,l))/(Omega_star**2)
     enddo
     enddo
     n_n0=abs(n_n0com)

     Phisks=0
     do p=-pmax+1,pmax-1
      do n=-nmax+1,nmax-1
      Phisks=Phisks+(kx(p)**2+ky(n)**2)*conjg(Phi_kft(p,n))*Phi_kft(p,n)/(Omega_star**2)/2.0
      enddo
    enddo    

    IF(verify_NL == 2)then
    do p=-pmax+1,pmax-1
      do n=1,nmax-1
       GOMG(-p,-n)=conjg(GOMG(p,n))
      enddo
    enddo
    Endif

    write(str,*)total_Chift
    if(Trim(ADJUSTL(str)) .eq. 'NaN')then
    write(*,*)"The total_Chift goes to NaN !!!"
    stop
    endif


     do n=1,nmax-1
     do p=-pmax+1,pmax-1
       D_kft(-p,-n)=(D_kft(p,n))
       total_Dft=total_Dft+D_kft(-p,-n)

       Chi_kft(-p,-n)=(Chi_kft(p,n))
     enddo
     enddo

! output file of energy, diffusive coefficient etc.
  if(dir_dia_count==1.and.restart==0)then   !open new output files at the first time in a run
      open(stdout+43,file="Phi_kft.txt",status='replace')
      open(stdout+45,file="total_Dft.txt",status='replace')
      open(stdout+46,file="Chi_kft.txt",status='replace')
      open(stdout+47,file="total_Chift.txt",status='replace')
      open(stdout+55,file="sqrn_n0ft.txt",status='replace')
      open(stdout+63,file="Eft.txt",status='replace')   
    elseif(dir_dia_count==1.and.(restart==1.or.restart==2))then   !open old output files at the beginning of the restart
      open(stdout+43,file="Phi_kft.txt",status='old')
      open(stdout+45,file="total_Dft.txt",status='old')
      open(stdout+46,file="Chi_kft.txt",status='old')
      open(stdout+47,file="total_Chift.txt",status='old')
      open(stdout+55,file="sqrn_n0ft.txt",status='old')
      open(stdout+63,file="Eft.txt",status='old')
    endif

    If(dir_dia_count==1.and.(restart==1.or.restart==2))then
     do i=1,(2*nmax-1)*(2*pmax-1)*nt/output_step
        read(stdout+43,*)readtempcom
     enddo
      do i=1,nt/output_step
        read(stdout+45,*)readtemp
      enddo
      do i=1,(2*nmax-1)*(2*pmax-1)*nt/output_step
        read(stdout+46,*)readtemp
      enddo
      do i=1,nt/output_step
        read(stdout+47,*)readtemp
      enddo
      do i=1,nt/output_step
        read(stdout+55,*)readtemp
      enddo
      do i=1,nt/output_step
        read(stdout+63,*)readtemp
      enddo          
   Endif

      write(stdout+45,102) total_Dft
      write(stdout+46,102) Chi_kft
      write(stdout+47,102) total_Chift
      write(stdout+55,102) n_n0
      write(stdout+63,102) Phisks

      do n=-nmax+1,nmax-1
       do p=-pmax+1,pmax-1
       write(stdout+43,"('(',e14.7,',',e14.7,')')")real(Phi_kft(p,n)),aimag(Phi_kft(p,n))
       enddo
      enddo

  !+++++++++++++++++++++++++++++++++++++++++++++++
  IF(verify_NL==1)then
   do n=0,nmax-1
   do p=-pmax+1,pmax-1
    if(n==0)then
      do i=1,N_mu
      do a=-N_FT+1,N_FT-1
      l=i*(2*N_FT-1)+a+(N_FT-1)
      tot_Entropyft=tot_Entropyft+ w_point(i)*conjg(F_kfta(p,n,l))*F_kfta(p,n,l)
      enddo
      enddo
!      i=0
!      a=0
!      l=i*(2*N_FT-1)+a+(N_FT-1)
!      tot_Entropyft=tot_Entropyft+ conjg(F_kfta(p,n,l))*F_kfta(p,n,l)
    elseif (n>0)then
      do i=1,N_mu
      do a=-N_FT+1,N_FT-1
      l=i*(2*N_FT-1)+a+(N_FT-1)
      tot_Entropyft=tot_Entropyft+ 2*(w_point(i)*conjg(F_kfta(p,n,l))*F_kfta(p,n,l))
      enddo
      enddo
!      i=0
!      a=0
!      l=i*(2*N_FT-1)+a+(N_FT-1)
!      tot_Entropyft=tot_Entropyft+ 2*(conjg(F_kfta(p,n,l))*F_kfta(p,n,l))
    endif
   enddo
   enddo

     do n=0,nmax-1
     do p=-pmax+1,pmax-1
        do  i=1,N_mu
        a=0
        l=i*(2*N_FT-1)+a+(N_FT-1)
        Fk_ftmucom(i)=Fk_ftmucom(i)+F_kfta(p,n,l)
        enddo
     enddo
     enddo

     do n=1,nmax-1
     do p=-pmax+1,pmax-1
       do  i=1,N_mu
        a=0
        l=i*(2*N_FT-1)+a+(N_FT-1)
        Fk_ftmucom(i)=Fk_ftmucom(i)+conjg(F_kfta(p,n,l))
       enddo
     enddo
     enddo
    Fk_ftmu(:)=2*PI*abs(Fk_ftmucom(:))

   do n=0,nmax-1
   do p=-pmax+1,pmax-1
      do i=1,N_mu
      do a=-N_FT+1,N_FT-1
      l=i*(2*N_FT-1)+a+(N_FT-1)
      Fk_ftcomT(i)=Fk_ftcomT(i)+F_kfta(p,n,l)
      enddo
      enddo
   enddo
   enddo

     do n=1,nmax-1
     do p=-pmax+1,pmax-1
      do i=1,N_mu
      do a=-N_FT+1,N_FT-1
      j=i*(2*N_FT-1)-a+(N_FT-1)
      Fk_ftcomT(i)=Fk_ftcomT(i)+conjg(F_kfta(p,n,j))
      enddo
      enddo
     enddo
     enddo
    Fk_ftmuT(:)=2*PI*abs(Fk_ftcomT(:))

     do p=-pmax+1,pmax-1
      do n=-nmax+1,nmax-1
      Phisqr=Phisqr+conjg(Phi_kft(p,n))*Phi_kft(p,n)/(Omega_star**2)
      enddo
    enddo    

    if(dir_dia_count==1.and.restart==0)then
      open(stdout+38,file="Phisqrft.txt",status='replace')
      open(stdout+44,file="D_kft.txt",status='replace')
      open(stdout+48,file="tot_Entropyft.txt",status='replace')
      open(stdout+49,file="Chi_muft.txt",status='replace')
      open(stdout+51,file="Fk_ftmu.txt",status='replace')
      open(stdout+52,file="Fk_ftmuT.txt",status='replace')
    elseif(dir_dia_count==1.and.(restart==1.or.restart==2))then
      open(stdout+38,file="Phisqrft.txt",status='old')
      open(stdout+44,file="D_kft.txt",status='old')
      open(stdout+48,file="tot_Entropyft.txt",status='old')
      open(stdout+49,file="Chi_muft.txt",status='old')
      open(stdout+51,file="Fk_ftmu.txt",status='old')
      open(stdout+52,file="Fk_ftmuT.txt",status='old')
    endif

    If(dir_dia_count==1.and.(restart==1.or.restart==2))then
      do i=1,nt/output_step
        read(stdout+38,*)readtemp
      enddo      
      do i=1,(2*nmax-1)*(2*pmax-1)*nt/output_step
        read(stdout+44,*)readtemp
      enddo
      do i=1,nt/output_step
        read(stdout+48,*)readtemp
      enddo
      do i=1,N_mu*nt/output_step
        read(stdout+49,*)readtemp
      enddo
      do i=1,N_mu*nt/output_step
        read(stdout+51,*)readtemp
      enddo
      do i=1,N_mu*nt/output_step
        read(stdout+52,*)readtemp
      enddo
    Endif

      write(stdout+38,102) Phisqr
      write(stdout+44,102) D_kft
      write(stdout+48,102) tot_Entropyft
      write(stdout+49,102) Chi_muft
      write(stdout+51,102) Fk_ftmu
      write(stdout+52,102) Fk_ftmuT
  Endif

  !-----------------------------------------------
  IF(verify_NL==2)then
    if(dir_dia_count==1.and.restart==0)then
    open(stdout+58,file="GOMG.txt",status='replace')
    open(stdout+59,file="Phi_kOMG.txt",status='replace')
    endif

    If(dir_dia_count==1.and.(restart==1.or.restart==2))then
     If(Stopnt==0)then
      open(stdout+58,file="GOMG.txt",status='old')
      open(stdout+59,file="Phi_kOMG.txt",status='old')
      do i=1,(2*nmax-1)*(2*pmax-1)*(nt)/output_step
        read(stdout+58,*)readtempcom
      enddo
      do i=1,(2*nmax-1)*(2*pmax-1)*(nt)/output_step
        read(stdout+59,*)readtempcom
      enddo
     ELSE
       if(nt==Stopnt .or. nt==Stopnt+1)then
        open(stdout+58,file="GOMG.txt",status='replace')
        open(stdout+59,file="Phi_kOMG.txt",status='replace')
       else
        open(stdout+58,file="GOMG.txt",status='old')
        open(stdout+59,file="Phi_kOMG.txt",status='old')

        do i=1,(2*nmax-1)*(2*pmax-1)*(nt-Stopnt)/output_step
        read(stdout+58,*)readtempcom
        enddo
        do i=1,(2*nmax-1)*(2*pmax-1)*(nt-Stopnt)/output_step
        read(stdout+59,*)readtempcom
        enddo
       endif
     ENDIF
    Endif

      do n=-nmax+1,nmax-1
       do p=-pmax+1,pmax-1
       write(stdout+58,"('(',e14.7,',',e14.7,')')")real(GOMG(p,n)),aimag(GOMG(p,n))
       enddo
      enddo


      do n=-nmax+1,nmax-1
       do p=-pmax+1,pmax-1
      write(stdout+59,"('(',e14.7,',',e14.7,')')")real(Phi_kft(p,n)),aimag(Phi_kft(p,n))
       enddo
      enddo

      deallocate(GOMG)
  Endif
  !+++++++++++++++++++++++++++++++++++++++++++++++

  if(dir_dia_count==1)then   !write the directory name of the output data at the first time in a run of this 3D code, and calculate linear disperson relation Omega
  open(stdout+990,file="data_file_name.txt",status='replace')
  open(stdout+991,file="cpu_timeFK.txt",status='old',POSITION='APPEND')
  write(stdout+990,&
"('iso',i2.2,'NL',f4.1,'nt',i6.6,'tsp',f9.6,'pmax',i2.2,'delta',f6.3,'mu',f6.3,'lambda',f5.2,'Omega_star',f6.2,'N_CY',i1)")&
 isotropic,Const_NL,ntmax,tstep,pmax,delta_1,mu_HK,lambda_0,Omega_star,N_CY!,F_k_int*((2*pmax-1)*(2*nmax-1)*(N_mu))
    write(stdout+991,*)"mu_point="
    do i=1,N_mu
    write(stdout+991,*)mu_point(i)
    enddo
    write(stdout+991,*)" "
    write(stdout+991,*)"w_point="
    do i=1,N_mu
    write(stdout+991,*)w_point(i)
    enddo
    write(stdout+991,*)" "
    write(stdout+991,*)"w_pointMod="
    do i=1,N_mu
    write(stdout+991,*)w_point(i)*FM(i)
    enddo
    write(stdout+991,*)" "
  FM_integ=0
  do i=1,N_mu
     FM_integ = FM_integ+w_point(i)*mu_point(i)*exp(-mu_point(i)**2)
  enddo
  write(stdout+991,*)'Test1=',FM_integ
  write(stdout+991,*)' '
  FM_integ=0
  do i=1,N_mu
     FM_integ = FM_integ+w_point(i)*mu_point(i)*exp(-mu_point(i))
  enddo
  write(stdout+991,*)'Test2=',FM_integ
  write(stdout+991,*)' '
  write(stdout+991,*)"FK"
  write(stdout+991,"('CDW=',I2)")CDW
  write(stdout+991,"('mugridtype=',I2)")mugridtype
  if(mugridtype==3)then
  write(stdout+991,"('MDmu12=',f5.2)")MDmu12
  endif
  write(stdout+991,*)"ntmax=",ntmax
  write(stdout+991,"('tstep=',f10.7)")tstep
  write(stdout+991,"('output_step=',I5)")output_step
  write(stdout+991,"('pmax=',I2)")pmax
  write(stdout+991,"('kxmax=',f6.3)")kxmax
  write(stdout+991,"('mumax=',f4.1)")mumax
  write(stdout+991,"('N_mu=',I2)")N_mu
  write(stdout+991,"('a_Ln=',f6.3)")a_Ln
  write(stdout+991,"('a_LTi=',f6.3)")a_LTi
  write(stdout+991,"('lambda_0=',f6.3)")lambda_0
  If(CDW==1)then
  write(stdout+991,"('lambda_D=',f6.3)")lambda_D
  write(stdout+991,"('AlphaA=',f8.3)")AlphaA
  Endif
  write(stdout+991,"('delta_1=',f6.3)")delta_1
  write(stdout+991,"('mu_HK=',f6.3)")mu_HK
  write(stdout+991,"('mu_LK=',f7.4)")mu_LK
  write(stdout+991,"('nu_DW=',f6.3)")nu_DW
  write(stdout+991,"('nu_ZF=',f6.3)")nu_ZF
  write(stdout+991,"('G=',f6.3)")G
  write(stdout+991,"('Epsilon=',f6.3)")Epsilon
  write(stdout+991,"('OMG=',f5.1)")Omega_star
  write(stdout+991,"('kpar=',f7.3)")kpar
  write(stdout+991,"('uD=',f8.4)")uD
  write(stdout+991,"('CDWid=',f5.2)")CDWid
  write(stdout+991,"('gamE=',f8.4)")gamE
  write(stdout+991,"('gamIC=',f8.4)")gamIC
  write(stdout+991,"('Beta=',f6.2)")Beta
  write(stdout+991,"('N_CY=',I2)")N_CY
  write(stdout+991,"('N_FT=',I2)")N_FT

  close(stdout+990)
  close(stdout+991)
  endif
 Endif

 if(myid==0 .and. verify_L==1)then
  call cpu_time(DBte)
  DBt=DBt+ (DBte-DBtb)
   if(nt==ntmax)then
   write(*,*)"Dia&Back time= ",DBt
   open(stdout+103,file="cpu_timeFK.txt",status='old',POSITION='APPEND')
   write(stdout+103,*)"Dia&Back time= ",DBt
   close(stdout+103)
   endif
 endif

Endif

!=================================================================================================================================
!CKinCH
IF(GK_FK_CK==2 .or. GK_FK_CK== -2)then
! Cyclotron kinetics
 If(myid==0)then
    tot_Entropycy=0
    total_Dcy=0
    total_Chicy=0
    D_kcy=0
    Chi_kcy=0
    G_kcy=(0,0)
    Chi_muft=0

    if(verify_L==1)then
    call cpu_time(DBtb)
    endif

    IF(verify_NL == 2)then
    allocate(GOMG(-pmax+1:pmax-1,-nmax+1:nmax-1))
    GOMG=(0,0)
    Endif


    do n=-nmax+1,nmax-1
    do p=-pmax+1,pmax-1
        Integ_numcy(p,n)=(0,0)  !set the integral variables to zero before integral
        do a=-N_CY+1,N_CY-1
        do i=1,N_mu
        l=(a+N_CY-1)*(N_mu+1)+i
        G_kcy(l,p,n)=F_kcya(l,p,n)+Jn(p,n,i,a)*Phi_kcy(p,n)*FM(i)/ratio_TiTe
        Integ_numcy(p,n)=Integ_numcy(p,n)+w_point(i)*Jn(p,n,i,a)*G_kcy(l,p,n)
        enddo
        enddo
        D_kcy(p,n)=real((0,1.0)*ky(n)*conjg(Phi_kcy(p,n))*Integ_numcy(p,n)/a_Lntmp)
        total_Dcy=total_Dcy+D_kcy(p,n)

        Integ_numcy(p,n)=(0,0)  !set the integral variables to zero before integral
        do a=-N_CY+1,N_CY-1
        do i=1,N_mu
        l=(a+N_CY-1)*(N_mu+1)+i
        Integ_numcy(p,n)=Integ_numcy(p,n)+w_point(i)*Jn(p,n,i,a)*G_kcy(l,p,n)*ratio_TiTe*mu_point(i)
        enddo
        enddo
        Chi_kcy(p,n)=real((0,1.0)*ky(n)*conjg(Phi_kcy(p,n))*Integ_numcy(p,n)/a_LTitmp)
        total_Chicy=total_Chicy+Chi_kcy(p,n)

      IF(verify_NL == 2)then
        Integ_numcy(p,n)=(0,0)  !set the integral variables to zero before integral
        do a=-N_CY+1,N_CY-1
        do i=1,N_mu    !gauss legendre integration
        l=(a+N_CY-1)*(N_mu+1)+i
        Integ_numcy(p,n)=Integ_numcy(p,n)+w_point(i)*mu_point(i)*Jn(p,n,i,a)*&
                         (F_kcya(l,p,n)+Jn(p,n,i,a)*Phi_kcy(p,n)*FM(i)/ratio_TiTe)
        enddo
        enddo
        GOMG(p,n)=Integ_numcy(p,n)
      Endif
    enddo
    enddo

     do n=-nmax+1,nmax-1
     do p=-pmax+1,pmax-1
       i=0
       a=0
       l=(a+N_CY-1)*(N_mu+1)+i
       n_n0com=n_n0com + F_kcya(l,p,n)*conjg(F_kcya(l,p,n))/(Omega_star**2)
     enddo
     enddo
     n_n0=abs(n_n0com)
     

     ni_n0=0
     do n=-nmax+1,nmax-1
     do p=-pmax+1,pmax-1
        n_n0com=(0,0)
        do a=-N_CY+1,N_CY-1
        do i=1,N_mu
        l=(a+N_CY-1)*(N_mu+1)+i
        n_n0com=n_n0com + Jn(p,n,i,a)*F_kcya(l,p,n)/Omega_star
        enddo
        enddo
        ni_n0=ni_n0+n_n0com*conjg(n_n0com)
     enddo
     enddo
     
     Phisks=0
     do p=-pmax+1,pmax-1
      do n=-nmax+1,nmax-1
      Phisks=Phisks+(kx(p)**2+ky(n)**2)*conjg(Phi_kcy(p,n))*Phi_kcy(p,n)/(Omega_star**2)/2.0
      enddo
    enddo
  
    write(str,*)total_Chicy
    if(Trim(ADJUSTL(str)) .eq. 'NaN')then
    write(*,*)"The total_Chicy goes to NaN !!!"
    stop
    endif

! output file of energy, diffusive coefficient etc.
    if(dir_dia_count==1.and.restart==0)then   !open new output files at the first time in a run
      open(stdout+34,file="Phi_kcy.txt",status='replace')
      open(stdout+51,file="total_Dcy.txt",status='replace')
      open(stdout+46,file="Chi_kcy.txt",status='replace')
      open(stdout+53,file="total_Chicy.txt",status='replace')
      open(stdout+55,file="sqrn_n0cy.txt",status='replace')
      open(stdout+62,file="ni_n0.txt",status='replace')
      open(stdout+63,file="Ecy.txt",status='replace')           
    elseif(dir_dia_count==1.and.(restart==1.or.restart==2))then  !open old output files at the beginning of the restart
      open(stdout+34,file="Phi_kcy.txt",status='old')
      open(stdout+51,file="total_Dcy.txt",status='old')
      open(stdout+46,file="Chi_kcy.txt",status='old')
      open(stdout+53,file="total_Chicy.txt",status='old')
      open(stdout+55,file="sqrn_n0cy.txt",status='old')
      open(stdout+62,file="ni_n0.txt",status='old')
      open(stdout+63,file="Ecy.txt",status='old')         
    endif

    If(dir_dia_count==1.and.(restart==1.or.restart==2))then
     do i=1,(2*nmax-1)*(2*pmax-1)*nt/output_step
        read(stdout+34,*)readtempcom
     enddo
      do i=1,nt/output_step
        read(stdout+51,*)readtemp
      enddo
      do i=1,(2*nmax-1)*(2*pmax-1)*nt/output_step
        read(stdout+46,*)readtemp
      enddo
      do i=1,nt/output_step
        read(stdout+53,*)readtemp
      enddo
      do i=1,nt/output_step
        read(stdout+55,*)readtemp
      enddo
      do i=1,nt/output_step
        read(stdout+62,*)readtemp
      enddo
      do i=1,nt/output_step
        read(stdout+63,*)readtemp
      enddo            
    Endif

      write(stdout+51,102) total_Dcy
      write(stdout+46,102) Chi_kcy
      write(stdout+53,102) total_Chicy
      write(stdout+55,102) n_n0
      write(stdout+62,102) ni_n0
      write(stdout+63,102) Phisks
      
      do n=-nmax+1,nmax-1
      do p=-pmax+1,pmax-1
       write(stdout+34,"('(',e14.7,',',e14.7,')')")real(Phi_kcy(p,n)),aimag(Phi_kcy(p,n))
      enddo
      enddo

  !+++++++++++++++++++++++++++++++++++++++++++++++
  IF(verify_NL==1)then

    Fk_cycomT(:)=(0,0)
    Fk_cymuT(:)=0
    Chi_mucy(:)=0

   do n=0,nmax-1
   do p=-pmax+1,pmax-1
      if(n==0)then
      do a=-N_CY+1,N_CY-1
      do i=1,N_mu
      l=(a+N_CY-1)*(N_mu+1)+i
      tot_Entropycy=tot_Entropycy+ w_point(i)*conjg(F_kcya(l,p,n))*F_kcya(l,p,n)
      enddo
      enddo
!      i=0
!      a=0
!      l=(a+N_CY-1)*(N_mu+1)+i
!      tot_Entropycy=tot_Entropycy+ conjg(F_kcya(l,p,n))*F_kcya(l,p,n)
      elseif (n>0)then
      do a=-N_CY+1,N_CY-1
      do i=1,N_mu
      l=(a+N_CY-1)*(N_mu+1)+i
      tot_Entropycy=tot_Entropycy+ 2*w_point(i)*conjg(F_kcya(l,p,n))*F_kcya(l,p,n)
      enddo
      enddo
!      i=0
!      a=0
!      l=(a+N_CY-1)*(N_mu+1)+i
!      tot_Entropycy=tot_Entropycy+ 2*conjg(F_kcya(l,p,n))*F_kcya(l,p,n)
      endif
   enddo
   enddo

    do n=-nmax+1,nmax-1
    do p=-pmax+1,pmax-1
      do a=-N_CY+1,N_CY-1
      do i=1,N_mu
      l=(a+N_CY-1)*(N_mu+1)+i
      Chi_mucy(i)=Chi_mucy(i)+real((0,1.0)*ky(n)*conjg(Phi_kcy(p,n))*w_point(i)*Jn(p,n,i,a)&
                             *G_kcy(l,p,n)*ratio_TiTe*mu_point(i)/a_LTitmp)
      enddo
      enddo
    enddo
    enddo

     do p=-pmax+1,pmax-1
      do n=-nmax+1,nmax-1
      Phisqr=Phisqr+conjg(Phi_kcy(p,n))*Phi_kcy(p,n)/(Omega_star**2)
      enddo
    enddo

   do n=-nmax+1,nmax-1
   do p=-pmax+1,pmax-1
     do a=-N_CY+1,N_CY-1
     do i=1,N_mu
     l=(a+N_CY-1)*(N_mu+1)+i
     Fk_cycomT(i)=Fk_cycomT(i)+F_kcya(l,p,n)
     enddo
     enddo
   enddo
   enddo
    Fk_cymuT(:)=2*PI*abs(Fk_cycomT(:))

    if(dir_dia_count==1.and.restart==0)then
      open(stdout+38,file="Phisqrcy.txt",status='replace')
      open(stdout+44,file="D_kcy.txt",status='replace')
      open(stdout+54,file="tot_Entropycy.txt",status='replace')
      open(stdout+52,file="Fk_cymuT.txt",status='replace')
      open(stdout+49,file="Chi_mucy.txt",status='replace')
      open(stdout+56,file="ne.txt",status='replace')
    elseif(dir_dia_count==1.and.(restart==1.or.restart==2))then
      open(stdout+38,file="Phisqrcy.txt",status='old')
      open(stdout+44,file="D_kcy.txt",status='old')
      open(stdout+54,file="tot_Entropycy.txt",status='old')
      open(stdout+52,file="Fk_cymuT.txt",status='old')
      open(stdout+49,file="Chi_mucy.txt",status='old')
      open(stdout+56,file="ne.txt",status='old')
    endif

    If(dir_dia_count==1.and.(restart==1.or.restart==2))then
      do i=1,nt/output_step
        read(stdout+38,*)readtemp
      enddo
      do i=1,(2*nmax-1)*(2*pmax-1)*nt/output_step
        read(stdout+44,*)readtemp
      enddo
      do i=1,nt/output_step
        read(stdout+54,*)readtemp
      enddo
      do i=1,N_mu*nt/output_step
        read(stdout+52,*)readtemp
      enddo
      do i=1,N_mu*nt/output_step
        read(stdout+49,*)readtemp
      enddo
      do i=1,(2*nmax-1)*(2*pmax-1)*nt/output_step
        read(stdout+56,*)readtempcom
      enddo      
    Endif

      write(stdout+38,102) Phisqr
      write(stdout+44,102) D_kcy
      write(stdout+54,102) tot_Entropycy
      write(stdout+49,102) Chi_mucy
      write(stdout+52,102) Fk_cymuT

        i=0
        a=0
        l=(a+N_CY-1)*(N_mu+1)+i
       do n=-nmax+1,nmax-1
       do p=-pmax+1,pmax-1
        write(stdout+56,"('(',e14.7,',',e14.7,')')")real(F_kcya(l,p,n)),aimag(F_kcya(l,p,n))
       enddo
       enddo      
  Endif
  !-----------------------------------------------
  IF(verify_NL==2)then
    if(dir_dia_count==1.and.restart==0)then
    open(stdout+58,file="GOMG.txt",status='replace')
    open(stdout+59,file="Phi_kOMG.txt",status='replace')
    endif

    If(dir_dia_count==1.and.(restart==1.or.restart==2))then
     If(Stopnt==0)then
      open(stdout+58,file="GOMG.txt",status='old')
      open(stdout+59,file="Phi_kOMG.txt",status='old')
      do i=1,(2*nmax-1)*(2*pmax-1)*(nt)/output_step
        read(stdout+58,*)readtempcom
      enddo
      do i=1,(2*nmax-1)*(2*pmax-1)*(nt)/output_step
        read(stdout+59,*)readtempcom
      enddo
     ELSE
       if(nt==Stopnt .or. nt==Stopnt+1)then
        open(stdout+58,file="GOMG.txt",status='replace')
        open(stdout+59,file="Phi_kOMG.txt",status='replace')
       else
        open(stdout+58,file="GOMG.txt",status='old')
        open(stdout+59,file="Phi_kOMG.txt",status='old')

        do i=1,(2*nmax-1)*(2*pmax-1)*(nt-Stopnt)/output_step
        read(stdout+58,*)readtempcom
        enddo
        do i=1,(2*nmax-1)*(2*pmax-1)*(nt-Stopnt)/output_step
        read(stdout+59,*)readtempcom
        enddo
       endif
     ENDIF
    Endif

      do n=-nmax+1,nmax-1
       do p=-pmax+1,pmax-1
       write(stdout+58,"('(',e14.7,',',e14.7,')')")real(GOMG(p,n)),aimag(GOMG(p,n))
       enddo
      enddo


      do n=-nmax+1,nmax-1
       do p=-pmax+1,pmax-1
      write(stdout+59,"('(',e14.7,',',e14.7,')')")real(Phi_kcy(p,n)),aimag(Phi_kcy(p,n))
       enddo
      enddo

      deallocate(GOMG)
  Endif
  !+++++++++++++++++++++++++++++++++++++++++++++++

  if(dir_dia_count==1)then   !write the directory name of the output data at the first time in a run of this 3D code, and calculate linear disperson relation Omiga
  open(stdout+990,file="data_file_name.txt",status='replace')
  open(stdout+991,file="cpu_timeCK.txt",status='old',POSITION='APPEND')
  write(stdout+990,&
"('nt',i6.6,'tsp',f7.4,'pmax',i2.2,'delta',f6.3,'muCD',f6.3,'nuDW',f6.3,'nuZF',f7.4,'lambda',f5.2,'Omega_st',f6.2,'N_CY',i1)")&!,'Int',f8.4)")&
 ntmax,tstep,pmax,delta_1,mu_HK,nu_DW,nu_ZF,lambda_0,Omega_star,N_CY
    write(stdout+991,*)"mu_point="
    do i=1,N_mu
    write(stdout+991,*)mu_point(i)
    enddo
    write(stdout+991,*)" "
    write(stdout+991,*)"w_point="
    do i=1,N_mu
    write(stdout+991,*)w_point(i)
    enddo
    write(stdout+991,*)" "
    write(stdout+991,*)"w_pointMod="
    do i=1,N_mu
    write(stdout+991,*)w_point(i)*FM(i)
    enddo
    write(stdout+991,*)" "
  FM_integ=0
  do i=1,N_mu
     FM_integ = FM_integ+w_point(i)*mu_point(i)*exp(-mu_point(i)**2)
  enddo
  write(stdout+991,*)'Test1=',FM_integ
  write(stdout+991,*)' '
  FM_integ=0
  do i=1,N_mu
     FM_integ = FM_integ+w_point(i)*mu_point(i)*exp(-mu_point(i))
  enddo
  write(stdout+991,*)'Test2=',FM_integ
  write(stdout+991,*)' '
  write(stdout+991,*)"CK"
  write(stdout+991,"('CDW=',I2)")CDW
  write(stdout+991,"('mugridtype=',I2)")mugridtype
  write(stdout+991,*)"ntmax=",ntmax
  write(stdout+991,"('tstep=',f10.7)")tstep
  write(stdout+991,"('output_step=',I5)")output_step
  write(stdout+991,"('pmax=',I2)")pmax
  write(stdout+991,"('kxmax=',f6.3)")kxmax
  write(stdout+991,"('mumax=',f4.1)")mumax
  write(stdout+991,"('N_mu=',I2)")N_mu
  write(stdout+991,"('a_Ln=',f6.3)")a_Ln
  write(stdout+991,"('a_LTi=',f6.3)")a_LTi
  write(stdout+991,"('lambda_0=',f6.3)")lambda_0
  If(CDW==1)then
  write(stdout+991,"('lambda_D=',f6.3)")lambda_D
  write(stdout+991,"('AlphaA=',f8.3)")AlphaA
  Endif
  write(stdout+991,"('delta_1=',f6.3)")delta_1
  write(stdout+991,"('mu_HK=',f6.3)")mu_HK
  write(stdout+991,"('mu_LK=',f7.4)")mu_LK
  write(stdout+991,"('nu_DW=',f6.3)")nu_DW
  write(stdout+991,"('nu_ZF=',f6.3)")nu_ZF
  write(stdout+991,"('G=',f6.3)")G
  write(stdout+991,"('Epsilon=',f6.3)")Epsilon
  write(stdout+991,"('OMG=',f5.1)")Omega_star
  write(stdout+991,"('kpar=',f7.3)")kpar
  write(stdout+991,"('uD=',f8.4)")uD
  write(stdout+991,"('CDWid=',f5.2)")CDWid
  write(stdout+991,"('gamE=',f8.4)")gamE
  write(stdout+991,"('gamIC=',f8.4)")gamIC
  write(stdout+991,"('D_IC=',f8.4)")D_IC
  write(stdout+991,"('Beta=',f6.2)")Beta
  write(stdout+991,"('N_CY=',I2)")N_CY
  write(stdout+991,"('N_FT=',I2)")N_FT
  write(stdout+991,"('CHpola=',I2)")CHpola
  close(stdout+990)
  close(stdout+991)
  endif
 Endif

 if(myid==0 .and. verify_L==1)then
  call cpu_time(DBte)
  DBt=DBt+ (DBte-DBtb)
   if(nt==ntmax)then
   write(*,*)"Dia&Back time= ",DBt
   open(stdout+103,file="cpu_timeCK.txt",status='old',POSITION='APPEND')
   write(stdout+103,*)"Dia&Back time= ",DBt
   close(stdout+103)
   endif
 endif

Endif


!=================================================================================================================================

       dir_dia_count=dir_dia_count+1
    102 format(e14.7)

end Subroutine diagnose
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine indir_diagnose(myid)

  use mpi
  use parameter
  implicit none
  integer,intent(in) :: myid
  real readtemp
  real NLconser
  real tot_Entropy
  real Phitemp
!====================================================================
!Fourier kinetics
  real tot_Entropyft
  real Phitempft
  real NLconserft
!====================================================================
!Cyclotron kinetics
  real tot_Entropycy
  real Phitempcy
  real NLconsercy

If(verify_NL==1)then
!=================================================================================================================================
!GK
IF(GK_FK_CK==0)then

    NLconser=0
    tot_Entropy=0
    Phitemp=0

 If(myid==0)then
  do p=-pmax+1,pmax-1
    do n=0,nmax-1
      do i=1,N_mu
      if(n==0)then
      NLconser=NLconser+ 2*real(w_point(i)*conjg(F_kmu(p,n,i))*NLSa(p,n,i))
      elseif (n>0)then
      NLconser=NLconser+ 2*2*real(w_point(i)*conjg(F_kmu(p,n,i))*NLSa(p,n,i))
      endif
      enddo
    enddo
  enddo


  do p=-pmax+1,pmax-1
    do n=0,nmax-1
      do i=1,N_mu
      if(n==0)then
      tot_Entropy=tot_Entropy+ real(w_point(i)*conjg(F_kmu(p,n,i))*F_kmu(p,n,i))
      elseif (n>0)then
      tot_Entropy=tot_Entropy+ 2*real(w_point(i)*conjg(F_kmu(p,n,i))*F_kmu(p,n,i))
      endif
      enddo
    enddo
  enddo


  do p=-pmax+1,pmax-1
    do n=0,nmax-1
    if(n==0)then
    Phitemp=Phitemp+Phi_k(p,n)*conjg(Phi_k(p,n))
    elseif (n>0)then
    Phitemp=Phitemp+2* Phi_k(p,n)*conjg(Phi_k(p,n))
    endif
    enddo
  enddo
  Phitemp=sqrt(Phitemp/(2*pmax-1)/(2*nmax-1))

  NLconser=NLconser/tot_Entropy/Phitemp


    if(indir_dia_count==1.and.restart==0)then
      open(stdout+31,file="NLconser.txt",status='replace')
    elseif(indir_dia_count==1.and.(restart==1.or.restart==2))then !open old output files at the beginning of the restart
      open(stdout+31,file="NLconser.txt",status='old')
    endif

    If(indir_dia_count==1.and.(restart==1.or.restart==2))then
      do i=1,nt/output_step-1
        read(stdout+31,*)readtemp
      enddo
    Endif

    write(stdout+31,102) NLconser

 Endif

Endif

!=================================================================================================================================
!CKinFH
IF(GK_FK_CK==1)then
  If(myid==0)then

  if(verify_L==1)then
  call cpu_time(DBtb)
  endif
!Fourier kinetics
    tot_Entropyft=0
    Phitempft=0
    NLconserft=0

  do p=-pmax+1,pmax-1
    do n=0,nmax-1
    if(n==0)then
      do i=1,N_mu
      do a=-N_FT+1,N_FT-1
      l=i*(2*N_FT-1)+a+(N_FT-1)
      tot_Entropyft=tot_Entropyft+ real(w_point(i)*conjg(F_kfta(p,n,l))*F_kfta(p,n,l))
      enddo
      enddo
    elseif (n>0)then
      do i=1,N_mu
      do a=-N_FT+1,N_FT-1
      l=i*(2*N_FT-1)+a+(N_FT-1)
      tot_Entropyft=tot_Entropyft+ 2*real(w_point(i)*conjg(F_kfta(p,n,l))*F_kfta(p,n,l))
      enddo
      enddo
    endif
    enddo
  enddo

  do p=-pmax+1,pmax-1
    do n=-nmax+1,nmax-1
      if(n .GE. 0)then
      do i=1,N_mu
      do a=-N_FT+1,N_FT-1
      l=i*(2*N_FT-1)+a+(N_FT-1)
      NLconserft=NLconserft+ 2*((w_point(i)*conjg(F_kfta(p,n,l))*NLSfta(p,n,l)))&
                           +conjg((w_point(i)*conjg(F_kfta(p,n,l))*NLSfta(p,n,l)))
      enddo
      enddo
      elseif (n<0)then
      do i=1,N_mu
      do a=-N_FT+1,N_FT-1
      j=i*(2*N_FT-1)-a+(N_FT-1)
      NLconserft=NLconserft+ 2*((w_point(i)*conjg(F_kfta(-p,-n,j))*NLSfta(-p,-n,j)))&
                           + conjg((w_point(i)*conjg(F_kfta(-p,-n,j))*NLSfta(-p,-n,j)))
      enddo
      enddo
      endif
    enddo
  enddo

  do p=-pmax+1,pmax-1
    do n=-nmax+1,nmax-1
    Phitempft=Phitempft+Phi_kft(p,n)*conjg(Phi_kft(p,n))
    enddo
  enddo
  Phitempft=sqrt(Phitempft/(2*pmax-1)/(2*nmax-1))


    if(indir_dia_count==1.and.restart==0)then
      open(stdout+74,file="NLconserft.txt",status='replace')
    elseif(indir_dia_count==1.and.(restart==1.or.restart==2))then !open old output files at the beginning of the restart
      open(stdout+74,file="NLconserft.txt",status='old')
    endif

    If(indir_dia_count==1.and.(restart==1.or.restart==2))then
      do i=1,nt/output_step-1
        read(stdout+74,*)readtemp
      enddo
    Endif

    write(stdout+74,102) NLconserft/tot_Entropyft/Phitempft

  if(verify_L==1)then
  call cpu_time(DBte)
  DBt=DBt+ (DBte-DBtb)
  endif

 Endif


Endif
!=================================================================================================================================
!CKinCH
IF(GK_FK_CK==2 .or. GK_FK_CK== -2)then
  If(myid==0)then

  if(verify_L==1)then
  call cpu_time(DBtb)
  endif

    NLconsercy=0
    tot_Entropycy=0
    Phitempcy=0


  IF(GK_FK_CK==2)then
   do n=-nmax+1,nmax-1
   do p=-pmax+1,pmax-1
      if(n .GE. 0)then
      do a=-N_CY+1,N_CY-1
      do i=1,N_mu
      l=(a+N_CY-1)*(N_mu+1)+i
      NLconsercy=NLconsercy+ 2*real(w_point(i)*conjg(F_kcya(l,p,n))*NLScya(l,p,n))
      enddo
      enddo
      elseif (n<0)then
      do a=-N_CY+1,N_CY-1
      do i=1,N_mu
      j=(-a+N_CY-1)*(N_mu+1)+i
      NLconsercy=NLconsercy+ 2*real(w_point(i)*conjg(F_kcya(j,-p,-n))*NLScya(j,-p,-n))
      enddo
      enddo
      endif
   enddo
   enddo
  ENDIF

  IF(GK_FK_CK== -2)then
   do n=-nmax+1,nmax-1
   do p=-pmax+1,pmax-1
      if(n .GE. 0)then
      do a=-N_CY+1,N_CY-1
      do i=1,N_mu
      l=(a+N_CY-1)*(N_mu+1)+i
      jL=(p+pmax-1)*(N_mu+1)*(2*N_CY-1)+(a+N_CY-1)*(N_mu+1)+i
      NLconsercy=NLconsercy+ 2*real(w_point(i)*conjg(F_kcya(l,p,n))*NLScyadp(jL,n))
      enddo
      enddo
      elseif (n<0)then
      do a=-N_CY+1,N_CY-1
      do i=1,N_mu
      j=(-a+N_CY-1)*(N_mu+1)+i
      jL=(-p+pmax-1)*(N_mu+1)*(2*N_CY-1)+(-a+N_CY-1)*(N_mu+1)+i
      NLconsercy=NLconsercy+ 2*real(w_point(i)*conjg(F_kcya(j,-p,-n))*NLScyadp(jL,-n))
      enddo
      enddo
      endif
   enddo
   enddo
  ENDIF

  do p=-pmax+1,pmax-1
    do n=0,nmax-1
      if(n==0)then
      do i=1,N_mu
      do a=-N_CY+1,N_CY-1
      l=(a+N_CY-1)*(N_mu+1)+i
      tot_Entropycy=tot_Entropycy+ real(w_point(i)*conjg(F_kcya(l,p,n))*F_kcya(l,p,n))
      enddo
      enddo
      elseif (n>0)then
      do i=1,N_mu
      do a=-N_CY+1,N_CY-1
      l=(a+N_CY-1)*(N_mu+1)+i
      tot_Entropycy=tot_Entropycy+ 2*real(w_point(i)*conjg(F_kcya(l,p,n))*F_kcya(l,p,n))
      enddo
      enddo
      endif
    enddo
  enddo


  do p=-pmax+1,pmax-1
    do n=-nmax+1,nmax-1
    Phitempcy=Phitempcy+Phi_kcy(p,n)*conjg(Phi_kcy(p,n))
    enddo
  enddo
  Phitempcy=sqrt(Phitempcy/(2*pmax-1)/(2*nmax-1))

  NLconsercy=NLconsercy/tot_Entropycy/Phitempcy


    if(indir_dia_count==1.and.restart==0)then
      open(stdout+74,file="NLconsercy.txt",status='replace')
    elseif(indir_dia_count==1.and.(restart==1.or.restart==2))then !open old output files at the beginning of the restart
      open(stdout+74,file="NLconsercy.txt",status='old')
    endif

    If(indir_dia_count==1.and.(restart==1.or.restart==2))then
      do i=1,nt/output_step-1
        read(stdout+74,*)readtemp
      enddo
    Endif

    write(stdout+74,102) NLconsercy


  if(verify_L==1)then
  call cpu_time(DBte)
  DBt=DBt+ (DBte-DBtb)
  endif

  Endif

Endif

!=================================================================================================================================


Endif

  indir_dia_count=indir_dia_count+1

  102 format(e14.7)



end Subroutine indir_diagnose
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gauss_legendre(x1,x2,x,w,n)
!---------------------------------------------------------
! gauss_legendre.f90
!
! PURPOSE:
!  Given the lower and upper limits of integration, x1 and
!  x2, and given n, this routine returns arrays x and w of
!  length n, containing the abscissas and weights of the
!  Gauss-Legendre n-point quadrature formula.
!
! NOTES:
!  Taken from "Numerical Recipes" routine GAULEG.
!
!  Book uses eps = 3e-14 -- this has been decreased.
!
! REVISIONS
! 06 Dec 00: jeff.candy@gat.com
!  Copied from book.
!  Recoded in f90 (with *explicit* variable declarations).
!  Added exception for n=1.
!---------------------------------------------------------

  implicit none

  real, parameter :: eps = 1e-7

  real, parameter :: pi = 3.141592653589793

  real, intent(in) :: x1
  real, intent(in) :: x2
  integer, intent(in) :: n

  real, dimension(n) :: x
  real, dimension(n) :: w

  real :: xm
  real :: xl
  real :: z
  real :: z1
  real :: p1
  real :: p2
  real :: p3
  real :: pp

  integer :: m
  integer :: j
  integer :: i

  xm = 0.5*(x2+x1)
  xl = 0.5*(x2-x1)

  ! Exception for n=1 is required:

  if (n == 1) then
     x(1) = xm
     w(1) = 2.0*xl
     return
  endif

  ! Roots are symmetric.  We need only find half of them.

  m = (n+1)/2

  ! Initialize to fail first do test

  z1 = -1.0

  ! Loop over roots.

  do i=1,m

     z = cos(pi*(i-0.25)/(n+0.5))

     do while (abs(z-z1) > eps)

        p1 = 1.0
        p2 = 0.0

        do j=1,n
           p3 = p2
           p2 = p1
           p1 = ((2*j-1)*z*p2-(j-1)*p3)/j
        enddo

        ! p1 is the Legendre polynomial.  Now compute its
        ! derivative, pp.

        pp = n*(z*p1-p2)/(z*z-1.0)
        z1 = z
        z  = z1-p1/pp

     enddo

     x(i)     = xm-xl*z
     x(n+1-i) = xm+xl*z
     w(i)     = 2.0*xl/((1.0-z*z)*pp*pp)
     w(n+1-i) = w(i)

!  write(*,*)"x",i,"=",x(i)
!  write(*,*)"w",i,"=",w(i)

  enddo

end subroutine gauss_legendre
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





     FUNCTION BESSJ (N,X)

!     This subroutine calculates the first kind modified Bessel function
!     of integer order N, for any REAL X. We use here the classical
!     recursion formula, when X > N. For X < N, the Miller's algorithm
!     is used to avoid overflows.
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.

      PARAMETER (IACC = 40,BIGNO = 1.D10, BIGNI = 1.D-10)
      REAL  X,BESSJ,BESSJ0,BESSJ1,TOX,BJM,BJ,BJP,SUM
      IF (N.EQ.0) THEN
      BESSJ = BESSJ0(X)
      RETURN
      ENDIF
      IF (N.EQ.1) THEN
      BESSJ = BESSJ1(X)
      RETURN
      ENDIF
      IF (X.EQ.0.) THEN
      BESSJ = 0.
      RETURN
      ENDIF
      TOX = 2./X
      IF (X.GT.FLOAT(N)) THEN
      BJM = BESSJ0(X)
      BJ  = BESSJ1(X)
      DO 11 J = 1,N-1
      BJP = J*TOX*BJ-BJM
      BJM = BJ
      BJ  = BJP
   11 CONTINUE
      BESSJ = BJ
      ELSE
      M = 2*((N+INT(SQRT(FLOAT(IACC*N))))/2)
      BESSJ = 0.
      JSUM = 0
      SUM = 0.
      BJP = 0.
      BJ  = 1.
      DO 12 J = M,1,-1
      BJM = J*TOX*BJ-BJP
      BJP = BJ
      BJ  = BJM
      IF (ABS(BJ).GT.BIGNO) THEN
      BJ  = BJ*BIGNI
      BJP = BJP*BIGNI
      BESSJ = BESSJ*BIGNI
      SUM = SUM*BIGNI
      ENDIF
      IF (JSUM.NE.0) SUM = SUM+BJ
      JSUM = 1-JSUM
      IF (J.EQ.N) BESSJ = BJP
   12 CONTINUE
      SUM = 2.*SUM-BJ
      BESSJ = BESSJ/SUM
      ENDIF
      RETURN
      END

      FUNCTION BESSJ0 (X)
      REAL  X,BESSJ0,AX,FR,FS,Z,FP,FQ,XX

!     This subroutine calculates the First Kind Bessel Function of
!     order 0, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.

      REAL  Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5 /1.D0,-.1098628627D-2,.2734510407D-4, &
      -.2073370639D-5,.2093887211D-6 /
      DATA Q1,Q2,Q3,Q4,Q5 /-.1562499995D-1,.1430488765D-3, &
      -.6911147651D-5,.7621095161D-6,-.9349451520D-7 /
      DATA R1,R2,R3,R4,R5,R6 /57568490574.D0,-13362590354.D0, &
      651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0 /
      DATA S1,S2,S3,S4,S5,S6 /57568490411.D0,1029532985.D0, &
      9494680.718D0,59272.64853D0,267.8532712D0,1.D0 /
      IF(X.EQ.0.D0) GO TO 1
      AX = ABS (X)
      IF (AX.LT.8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      BESSJ0 = FR/FS
      ELSE
      Z = 8./AX
      Y = Z*Z
      XX = AX-.785398164
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      BESSJ0 = SQRT(.636619772/AX)*(FP*COS(XX)-Z*FQ*SIN(XX))
      ENDIF
      RETURN
    1 BESSJ0 = 1.D0
      RETURN
      END
! ---------------------------------------------------------------------------
      FUNCTION BESSJ1 (X)
      REAL  X,BESSJ1,AX,FR,FS,Z,FP,FQ,XX
!     This subroutine calculates the First Kind Bessel Function of
!     order 1, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.
      REAL *8 Y,P1,P2,P3,P4,P5,P6,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5 /1.D0,.183105D-2,-.3516396496D-4,  &
      .2457520174D-5,-.240337019D-6 /,P6 /.636619772D0 /
      DATA Q1,Q2,Q3,Q4,Q5 /.04687499995D0,-.2002690873D-3,   &
      .8449199096D-5,-.88228987D-6,.105787412D-6 /
      DATA R1,R2,R3,R4,R5,R6 /72362614232.D0,-7895059235.D0, &
      242396853.1D0,-2972611.439D0,15704.48260D0,-30.16036606D0 /
      DATA S1,S2,S3,S4,S5,S6 /144725228442.D0,2300535178.D0, &
      18583304.74D0,99447.43394D0,376.9991397D0,1.D0 /

      AX = ABS(X)
      IF (AX.LT.8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      BESSJ1 = X*(FR/FS)
      ELSE
      Z = 8./AX
      Y = Z*Z
      XX = AX-2.35619491
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      if(X.GE.0.0)then
      BESSJ1 = SQRT(P6/AX)*(COS(XX)*FP-Z*SIN(XX)*FQ)*S6
      else
      BESSJ1 = SQRT(P6/AX)*(COS(XX)*FP-Z*SIN(XX)*FQ)*(-1.0)*S6
      endif
      ENDIF
      RETURN
      END

!End of file Tbessj.f90



! ----------------------------------------------------------------------
      FUNCTION BESSI(N,X)
!
!     This subroutine calculates the first kind modified Bessel function
!     of integer order N, for any REAL X. We use here the classical
!     recursion formula, when X > N. For X < N, the Miller's algorithm
!     is used to avoid overflows.
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
!
      PARAMETER (IACC = 40,BIGNO = 1.D10, BIGNI = 1.D-10)
      REAL *8 X,BESSI,BESSI0,BESSI1,TOX,BIM,BI,BIP
      IF (N.EQ.0) THEN
      BESSI = BESSI0(X)
      RETURN
      ENDIF
      IF (N.EQ.1) THEN
      BESSI = BESSI1(X)
      RETURN
      ENDIF
      IF(X.EQ.0.D0) THEN
      BESSI=0.D0
      RETURN
      ENDIF
      TOX = 2.D0/X
      BIP = 0.D0
      BI  = 1.D0
      BESSI = 0.D0
      M = 2*((N+INT(SQRT(FLOAT(IACC*N)))))
      DO 12 J = M,1,-1
      BIM = BIP+DFLOAT(J)*TOX*BI
      BIP = BI
      BI  = BIM
      IF (ABS(BI).GT.BIGNO) THEN
      BI  = BI*BIGNI
      BIP = BIP*BIGNI
      BESSI = BESSI*BIGNI
      ENDIF
      IF (J.EQ.N) BESSI = BIP
   12 CONTINUE
      BESSI = BESSI*BESSI0(X)/BI
      RETURN
      END
! ----------------------------------------------------------------------
! Auxiliary Bessel functions for N=0, N=1
      FUNCTION BESSI0(X)
      REAL *8 X,BESSI0,Y,P1,P2,P3,P4,P5,P6,P7,  &
      Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,1.2067492D0,  &
      0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1, &
      0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,  &
      0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF(ABS(X).LT.3.75D0) THEN
      Y=(X/3.75D0)**2
      BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
      AX=ABS(X)
      Y=3.75D0/AX
      BX=EXP(AX)/SQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI0=AX*BX
      ENDIF
      RETURN
      END
! ----------------------------------------------------------------------
      FUNCTION BESSI1(X)
      REAL *8 X,BESSI1,Y,P1,P2,P3,P4,P5,P6,P7,  &
      Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,  &
      0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1, &
      -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1, &
      -0.2895312D-1,0.1787654D-1,-0.420059D-2/
      IF(ABS(X).LT.3.75D0) THEN
      Y=(X/3.75D0)**2
      BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      AX=ABS(X)
      Y=3.75D0/AX
      BX=EXP(AX)/SQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI1=AX*BX
      ENDIF
      RETURN
      END
! end of file tbessi.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDEIGEN(Matrix, n, x,EigenVal, steps)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, steps  !n = order of matrix, steps = number of iterations
	complex, INTENT(IN), DIMENSION(n,n) :: Matrix(n,n)  !Input Matrix
	complex, INTENT(INOUT), DIMENSION(n) :: x !Eigenvector
	complex, INTENT(INOUT) :: EigenVal !Eigenvalue
	INTEGER :: i, j

	x  = (1,1) !Initialize eigen vector to any value.

	DO i = 1, steps
		CALL MULMATRIX(Matrix, x, n)       !Multiply input matrix by eigenvector
		CALL FINDLARGEST(x, n, EigenVal)   !Find eigenvalue
		IF(EigenVal == (0,0)) EXIT
		DO j = 1, n                        !Find eigenvector
			x(j) = x(j)/EigenVal
		END DO
	END DO

END SUBROUTINE FINDEIGEN

SUBROUTINE MULMATRIX(a, b, n)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n !matrix size
	complex, INTENT(IN), DIMENSION(n,n) :: a  !Matrix of order > 1
	complex, INTENT(INOUT), DIMENSION(n) :: b !1x1 matrix

	INTEGER i, j
	complex, DIMENSION(n) :: temp !temporary matrix

	temp = 0

	!These two loops to the multiplication
	DO i = 1, n
		DO j = 1, n
			temp(i) = temp(i) + a(i,j)*b(j)
		END DO
	END DO
	b = temp

END SUBROUTINE MULMATRIX

SUBROUTINE FINDLARGEST(x, n, l)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n
	complex, INTENT(IN), DIMENSION(n) :: x
	complex, INTENT(INOUT) :: l !Largest value
    real l_im

	INTEGER :: i
	!Algorithm is easy
	!Let the largest number be the first one.
	!If you find a number larger than it, store this number and then continue
	l_im= real(x(1))
    l=x(1)
	DO i = 2, n
		IF (real(x(i)) > l_im)then
        l_im = real(x(i))
        l=x(i)
        Endif
	END DO
    !write(*,*)'eigen value',x

END SUBROUTINE FINDLARGEST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION gammln(xx)
REAL gammln,xx
INTEGER j
DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
!Internal arithmetic will be done in double precision, a nicety that you can omit if ve-gure
!accuracy  is  good  enough.
SAVE cof,stp
DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,&
   24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
   -.5395239384953d-5,2.5066282746310005d0/
x=xx
y=x
tmp=x+5.5d0
tmp=(x+0.5d0)*log(tmp)-tmp
ser=1.000000000190015d0
do  j=1,6
y=y+1.d0
ser=ser+cof(j)/y
enddo
gammln=tmp+log(stp*ser/x)
return
END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE gaulag(x,w,n,alf)
INTEGER n,MAXIT
REAL alf,w(n),x(n)
DOUBLE PRECISION EPS
PARAMETER (EPS=3.D-14,MAXIT=10)   !Increase EPS if you dont have this precision.
!USES gammln
!Given alf, the parameter of the Laguerre polynomials, this routine returns arrays x(1:n)
!and w(1:n)containing the abscissas and weights of the n-point Gauss-Laguerre quadrature
!formula.  The smallest abscissa is returned in x(1), the largest in x(n).
INTEGER i,its,j
REAL ai,gammln
DOUBLE PRECISION p1,p2,p3,pp,z,z1
!High  precision is  a good idea  for this routine.
do  i=1,n !Loop over the desired roots.
if(i.eq.1)then !Initial guess for the smallest root.
z=(1.+alf)*(3.+.92*alf)/(1.+2.4*n+1.8*alf)
else if(i.eq.2)then !Initial guess for the second root.
z=z+(15.+6.25*alf)/(1.+.9*alf+2.5*n)
else! Initial guess for the other roots.
ai=i-2
z=z+((1.+2.55*ai)/(1.9*ai)+1.26*ai*alf/&
 (1.+3.5*ai))*(z-x(i-2))/(1.+.3*alf)
endif
do its=1,MAXIT !Refinement by Newtons method.
p1=1.d0
p2=0.d0
do j=1,n !Loop up the recurrence relation to get the Laguerre
p3=p2             ! polynomial evaluated at z. p3=p2
p2=p1
p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j
enddo
!p1 is now the desired Laguerre polynomial.  We next compute pp, its derivative, by
!a standard relation involving also p2, the polynomial of one lower order.
pp=(n*p1-(n+alf)*p2)/z
z1=z
z=z1-p1/pp !Newtons formula.
if(abs(z-z1).le.EPS)goto 1
enddo
write(*,*)'too many iterations in gaulag'
continue
1      x(i)=z !Store the root and the weight.4.5 Gaussian Quadratures and Orthogonal Polynomials 147
w(i)=-exp(gammln(alf+n)-gammln(float(n)))/(pp*n*p2)
enddo
return
END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
