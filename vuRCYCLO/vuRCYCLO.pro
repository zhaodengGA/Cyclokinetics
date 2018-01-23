;
;CODE NAME:
;   vuRCYCLO
;
;
; AUTHOR:
;  ZHAO Deng
;  email: zhao.deng@foxmail.com
;
;
; ADVISER: 
;   R. E. Waltz
;   
;   
; PURPOSE: 
;  To plot the data from rCYCLO code, and analyze the physics.
;   
;   
; HOW TO START:
; 1. Make sure rCYCLO output files needed for plotting and all the vuRCYCLO .pro files are in the same directory.
; 2. Compile cgplot first.
; 3. Compile vuRCYCLO
; 4. Run vuRCYCLO
; 
; 
; CODE STRUCTURE:
; Main plots:       To plot the basic information of a case, including: total Chi vs time, total_D vs time,
;                   and growth rate (and frequency) vs ky of gyrokinetics(GK), cyclokinetics in Fourier harmonic
;                   representation (CKinFH), and cyclokinetics in cyclotron harmonic representation (CKinCH) 
;                   respectively.
; 
; More details:     To plot more details e.g. the 3D plot of Chi, the distribution function in velocity space.
; Colorful plot:    To plot the colorful contour figures of growth rate and frequency, and the eddy in (x,y) space
;                   evolution in time. Pay attention to that the 'Main plot' and 'More details' applies the default.  
;                   system color, however, the color table used for colorful plot is given by the file 'color.dat'.
;                   So there may be a color display problem when operation switch between them. To solve this problem,
;                   you need to restart idl.
;                   
;                    
; Movie:            Disable, the eddy movie is plotted in 'Colorful plot'
; ColorMap:         To plot the color table for the Colorful plot
; 
;
;; INPUT FILES: 
; inputRCYCLO.txt   -- in order to read the grid parameters and physics variables.
; inputvuRCYCLO.txt -- in order to read the control parameters for IDL plotting.
; 
; 
; CONTROL VARIABLES INTRODUCE:
; The variables in inputvuRCYCLO.txt is used for control the idl plot.
; 
; avenum=4             The number of pieces for taking average of energy overplot.          
; plotstep=1           The interval of kx lines in Omega overplot.
; sattime=0.0          For setting the starting point of time average.
; Intermit_Plot=0      A switch of output intermittency and average value. 
;                      1 for output intermit&average ON, 0 for OFF 
; overplot_log=0       A switch between normal axis and logarithmic axis in 3D plot. 
;                      0 for normal, 1 for ylog, 2 for xlog and ylog
; Az=-60               To rotate a certain angle along z axis in 3D plot
; Ax=35                To rotate a certain angle along x axis in 3D plot
; ix1=0                If ix1=0 & ix2=0, take the time=0 moment data for Calculate botton.
; ix2=0                 
; Switch_sat=1         A switch of outputting the ave_time given by saturated time in colorfulplot. 
;                      1 for ON, 0 for OFF.  
; Omega_min=0.01       The minimum value of y axis. 0 for xmin=min(x), >0 for xmin=Omega_min
; Omega_max=60         The maximum value of y axis.  0 for xmax=max(x), >0 for xmax=Omega_max
; Phi_logmin=0.0001    The minimum value of z axis for 3D log plot of |Phi|
; Chi_logmin=0.000001  The minimum value of z axis for 3D log plot of |Chi|
; OMGpoint=10          The frequency region chosen to plot high frequency ion cyclotron (IC) modes in Omega vs ky plot.  
; freqPer=0.2          The pick up frequency range for choosing IC modes. 
;                      e.g. freqPer=0.2 means if the frequency of modes are 20% smaller or larger than OMGpoint, 
;                      it will be chosen as the first harmonic of IC mode. 
; begintime=0          The begin time of eddy movie show.
;
; 
; BUTTON INTRODUCE:
; *total Chi (GK) (to plot total Chi vs time) window:
; Plot:                 Start plot with the default parameters.
; t+ or t-:             To adjust the point for taking the average of a saturation period.
; aveT and aveP:        Used for take the average of two state situation.
;                       aveT is the switch between take the average of one state or two states.
;                       aveP is change the current point to be adjusted by t+ or t-.
; ZOOM in or ZOOM out:  To adjust the range of the y axis.
; Yrange:               To decide the way of adjusting the yrange. Such as: fix the top, middle or bottom part.
; log:                  A switch between normal axis and logarithmic axis.
; PDF:                  To show the plot of the particle distribution function
; Calculate:            TO calculate anything plugged in there.
; Done:                 Close this window.
;      
; *Omega (GK) (to plot  growth rate (and frequency) vs ky) window:
; Plot(re):             To plot the real part of Omega, which means frequency.
; Plot(im):             To plot the imaginary part of Omega, which means growth rate.
; MS/Init:              To switch between plotting Matrix Solver eigenvalue method result and plotting the time initial
;                       value method result.
; 3D:                   To show the 3D plot
; 
; *F_mu (to plot delta f vs mu and time) window:
; Z+ or Z-:             To rotate a certain angle along z axis.           
; X+ or X-:             To rotate a certain angle along x axis.     
; zmin*10 or zmin/10:   To change the minimum value of zrange.
; jump*2 or jump/2:     To change the interval of sampling.   
; 2D:                   To show the 2D plot
; 
; *Omega_matr (to plot the colorful contour figures of growth rate and frequency) window:
; Type:                 To switch between GK, CKinCH and CKinFH
; 
; *eddy movie (to plot the eddy in (x,y) space evolution in time) window:
; velocity:             To show the velocity vector plot.
;    
;
; 
; REVISE TIME & CONTENT:
; 2014-8-6: add annotation
;
;
;*****************************************************************************************************************
; OTHER SPECIFIED FIGURES:
; 
; *** How to plot the frequncy spectrum of Chi (|Phi|^2 is the same) new method (on 2014.7.11)
;
; 1. Apply the after '2014.7.10.f90' code, and set 'verify_NL=2'
;
; (a). You can get the 'Phi_kOMG.txt' and 'GOMG.txt' from nt=0 to nt=ntmax, by set Stopnt=0 in 'inputrCYCLO.txt'.
;
; (b). Or you can restart from any old running restart point. Then 'Phi_kOMG.txt' and 'GOMG.txt' will be from nt=Stopnt 
; to nt=ntmax. Remember to set Stopnt in 'inputrCYCLO.txt' be the restart point.
; 
; PS: If you use a different output_step value as the original one, then remeber to set fortran code text some 'output_step' 
; equal to the smaller one to prevent output/input error. 
; 
; 2. Use the 'Chi_vs_nnew.pro' code control=2, Begintime='the begin ploting point, 'jump' is jumping some time points in 
; order to reduce the calculation. (better be 1).Then IDL will output the freqency result in file 'y.txt' (output 'y.txt'
; will save time for repeating operations).  
;
; 3. Set control=2.1 or 2.2 to continue to plot Chi vs omega or Integral_Chi vs omega.
;
; 4. As an code verification, remember to recover Parseval's Theorem.
;
;
; ** A fast way to get 'y.txt' file on cluster:
; Since the file for freq spectrum plot is very large, the transfor of the file to a PC is very terrible.
; Thus, we offer a way to calculate the short 'y.txt' file on cluster. Following by:
; 1.Put a suit of idl code in to the cluster, most importantly, including 'idl2014.7.12_rCYCLO.pro'
; 
; 2.When Stopnt=0, set the 'Begintime' in the 'inputvuRCYCLO.txt' file to be the  begin plotting point, and set ntmax be the end
; plotting point. When Stopnt>0, ploting will start from Stopnt to ntmax.
; 
; 3.run 'idl2014.7.12_rCYCLO.pro', and click the 'Energy(DW, ZF, GAM)', then will generate 'y.txt' file.
;
;
;*****************************************************************************************************************
;
;
pro twoD_time_trace_see,nametemp,windownum,group=group

;  common GLOBAL 

  ;;-----------------------------------------------
  ;; Private (local) data
  ;;
;  common PRIVATE_time_trace,widget
  common plotaxis,zoomz,zoombackup0,zoombackup1,log,half,freq,i_ptype,satdelta_t,name,yrange_style
  common main_plot,plotidh,ntime,ndata,data_name,windowmark,basemark

  ;;-----------------------------------------------

  ;;------------------------------------------
  ;; Return conditions
  ;;
 ; if exists_diff eq 0 then return
 ; if xregistered('diffusion_ave_see') then return
  ;;------------------------------------------
  name=nametemp
  base = widget_base(title=name[0],$
                     /column)
   ;;set default value

  defsysv,'!aveT',0 ;; ave TYPE
  defsysv,'!aveP',0 ;; ave POINT
  defsysv,'!tminus',0
  defsysv,'!tplus',0
  defsysv,'!point1',0
  defsysv,'!point2',0
  defsysv,'!point3',0
  defsysv,'!point4',0

  i_ptype=0
  log=0
  zoomz = 1.0  ;zoomz is the parameter which used to adjust the range of the axis
  satdelta_t=0.0  ;satdelta_t is the parameter which used to adjust the saturation time when the user click the button "t+" or "t-"
  yrange_style=1  ;yrange_style is the parameter which used to choose plot the top, middle or bottom part of the picture.
                  ;for the value of yrange_style/3, 0 for top, 1 for middle, 2 for bottom;
  ;i_tp = 0
  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row1, $
                    value='Plot', $
                    uvalue=0)

  x = widget_button(row1, $
                    value='t+', $
                    uvalue=1)

  x = widget_button(row1, $
                    value='t-', $
                    uvalue=2)

  x = widget_button(row1, $
                    value='aveT', $
                    uvalue=3)

  x = widget_button(row1, $
                    value='aveP', $
                    uvalue=4)

  x = widget_button(row1, $
                    value='ZOOM in', $
                    uvalue=5)

  x = widget_button(row1, $
                    value='ZOOM out', $
                    uvalue=6)
  x = widget_button(row1, $
                    value='Yrange', $
                    uvalue=7)
  x = widget_button(row1, $
                    value='units', $
                    uvalue=8)

  x = widget_button(row1,$
                    value='log',$
                    uvalue=9)

  x = widget_button(row1,$
                    value='TYPE',$
                    /menu)

  tlevels=['Line Plot','PDF']
  for i=0,1 do begin
     x1 = widget_button(x,$
                        value=tlevels[i],$
                        uvalue=20+i)
  endfor

  x = widget_button(row1, $
                    value='Calculate', $
                    uvalue=10)

  x = widget_button(row1, $
                    value='Done', $
                    uvalue=11)

  ;;----------------------------------------------------------
  ;; DRAW WIDGET and CONTROL
  ;;----------------------------------------------------------

  draw = widget_draw(base,     $
                     xsize=700, $
                     ysize=500)

  widget_control, base, $
    ;set_uvalue=state,$
    /no_copy, $
    /realize

  ;!plotkspecid=!D.WINDOW
    (*windowmark)[windownum]=!D.window
    basemark[windownum]=base


  xmanager,'twoD_time_trace', $
    base,$
;    event='energy_time_trace_event',$
    group_leader=group


end

;*******************************************************************************
pro twoD_time_trace_event,event

  common startup,number_plot,fpath,ncolor,color_value,plotid
  common plotaxis,zoomz,zoombackup0,zoombackup1,log,half,freq,i_ptype,satdelta_t,name,yrange_style
  common main_plot,plotidh,ntime,ndata,data_name,windowmark,basemark

  widget_control, event.id, $
    get_uvalue=uvalue
 ; wset, widget

openr,lun,'inputRCYCLO.txt',/get_lun
  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Control_variables',plot_name_temp,19)
  endwhile
  thisFormat='(8x,i2,/,9x,i2,/,4x,i2)'
  readf,lun,restart,GK_FK_CK,CDW,Format=thisFormat
  thisFormat='(8x,i2,/,11x,i2,/,12x,i5,/,11x,i4,/,9x,f5.3,/,/,/,/,7x,i10)'
  readf,lun,muDtype,mugridtype,output_step,backup_num,Const_nl,Stopnt,Format=thisFormat

  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Grid_variables',plot_name_temp,16)
  endwhile
  thisFormat='(6x,i10,/,6x,f12.9,/,6x,f10.3,/,6x,f10.3,/,5x,i5,/,5x,i5,/,6x,i8.3,/,5x,i5,/,5x,i6,/,5x,i6)'
  readf,lun,ntmax,tstep,kxmax,kymax,pmax,nmax,mumax,N_mu,$
      N_FT,N_CY,Format=thisFormat

  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Physics_variables',plot_name_temp,19)
  endwhile
  thisFormat='(11x,f7.4,/,5x,f5.3,/,6x,f6.3,/,9x,f5.3,/,9x,f5.3,/,9x,f5.3,/,7x,f9.3,/,8x,f5.3,/,6x,f5.3,/,6x,f5.3,/,8x,f12.7,/,8x,f8.5,/,/,/,6x,f8.4)'
  readf,lun,Omega_star,a_Ln,a_LTi,lambda_n,lambda_0,lambda_D,AlphaA,delta_1,$
  mu_HK,mu_LK,F_k_int,Epsilon,gamIC, Format=thisFormat

free_lun,lun



openr,lun,'inputvuRCYCLO.txt',/get_lun
  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Control_variables',plot_name_temp,19)
  endwhile
  thisFormat='(7x,i5,/,9x,i5,/,8x,f10.3,/,14x,i2,/,13x,i2,/,/,/,4x,i6,/,4x,i6,/,11x,i4,/,13x,i2,/,5x,i5,/,8x,i2)'
  readf,lun,avenum,plotstep,sattime,Intermit_Plot,overplot_log,ix1,ix2,Switch_sat,quicktrigger,jump,control,Format=thisFormat

  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Range_control',plot_name_temp,15)
  endwhile
  thisFormat='(10x,f15.12,/,10x,f12.4,/,11x,f15.10,/,11x,f15.12,/,9x,i5,/,8x,f10.4,/,10x,f12.3)'
  readf,lun,Omega_min,Omega_max,Phi_logmin,Chi_logmin,OMGpoint,freqPer,begintime,Format=thisFormat

free_lun,lun

;  openr,lun,'nt_Odd.txt',/get_lun
;  readf,lun,nt_Odd
;  free_lun,lun
;  openr,lun,'nt_Even.txt',/get_lun
;  readf,lun,nt_Even
;  free_lun,lun
;  IF(ntmax GE max([nt_Odd,nt_Even]))Then Begin
;  ntmax=max([nt_Odd,nt_Even])
;  Endif
;  print,'ntmax=',ntmax

  Time=findgen(ntmax/output_step+1)*(tstep*output_step)
  kx=findgen(2*pmax-1)*kxmax/(pmax-1)-kxmax
  ky=findgen(2*nmax-1)*kymax/(nmax-1)-kymax
  kxhalf=findgen(pmax)*kxmax/(pmax-1)
  kyhalf=findgen(nmax)*kymax/(nmax-1)

  ;set the initial value of the four average points

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  ;;set the click botton action act on the window you want.
  For i=0,200-1 Do begin
  IF(basemark[i] eq event.top)then begin
    IF((*windowmark)[i] ne !D.window)then begin
    wset,(*windowmark)[i]
    Endif
    case(i) of
    0: name=['total energy (conj(F)*J0*Phi+F*J0*conj(Phi))/2 (3DGK)','energy','total_energy.txt']
    1: name=['total D (GK)','D','total_D_3G.txt']
    2: name=['total Entropy (GK)','entropy','total_entropy.txt']
    3: name=['total '+'$\chi$$\downi$'+' (GK)','$\chi$$\downi$','total_Chi.txt']
    15: name=['Energy(DW, ZF, GAM)','Energy','tot_ene_GAM.txt']
    18: name=['total '+'|$\Phi$| (GK)','|$\Phi$|','Phi_k.txt']
    29: name=['total D (CKinFH)','D','total_Dft.txt']
    30: name=['total '+'$\chi$$\downi$'+' (CKinFH)','$\chi$$\downi$','total_Chift.txt']
    31: name=['total Entropy (CKinFH)','Entropy','tot_Entropyft.txt']
    32: name=['total '+'|$\Phi$| (CKinFH)','|$\Phi$|','Phi_kft.txt']
    45: name=['total '+'|$\Phi$| (CKinCH)','|$\Phi$|','Phi_kcy.txt']
    46: name=['total entropy (CKinCH)','entropy','tot_Entropycy.txt']
    47: name=['Data processing','',' .txt']
    48: name=['Others by time','combine',' .txt']
    ; 25: name=['Phi^2 (3DGK)','|Phi| ^2','total_Phi^2.txt']
    100: name=['total '+'$\chi$$\downi$'+' (GK)','$\chi$$\downi$','total_Chi.txt']
    101: name=['total D (GK)','D','total_D_3G.txt']    
    103: name=['total '+'$\chi$$\downi$'+' (CKinFH)','$\chi$$\downi$','total_Chift.txt']
    104: name=['total D (CKinFH)','D','total_Dft.txt']
    106: name=['total '+'$\chi$$\downi$'+' (CKinCH)','$\chi$$\downi$','total_Chicy.txt']
    107: name=['total D (CKinCH)','D','total_Dcy.txt']
    109: name=['($\delta$n$\downe$/n$\down0$)$\up2$ (GK)','($\delta$n$\downe$/n$\down0$)$\up2$','sqrn_n0.txt']
    110: name=['($\delta$n$\downe$/n$\down0$)$\up2$ (CKinFH)','($\delta$n$\downe$/n$\down0$)$\up2$','sqrn_n0ft.txt']
    111: name=['($\delta$n$\downe$/n$\down0$)$\up2$ (CKinCH)','($\delta$n$\downe$/n$\down0$)$\up2$','sqrn_n0cy.txt']    
    112: name=['electrostatic Energy E (GK)','E$\upGK$','E.txt']
    113: name=['electrostatic Energy E (CKinFH)','E$\upFK$','Eft.txt']
    114: name=['electrostatic Energy E (CKinCH)','E$\upCK$','Ecy.txt']
    endcase
  Endif
  Endfor

  IF(satdelta_t EQ 0)then begin
  satdelta_t=sattime    ;;set the first left average line point to be the sattime
  Endif

  case (uvalue) of

     0:begin
     ;wset,(*ndata)[windownum]
     goto, plot_it
     end
     1: begin
        !tplus=!tplus+1
        goto, plot_it
     end

     2: begin
        !tminus=!tminus+1
        goto, plot_it
     end

     3: begin
        !aveT = (!aveT+1) mod 2
        goto, plot_it
     end

     4: begin
        !aveP = (!aveP+1) mod 4
        goto, plot_it
     end

     5: begin
        zoomz = 2.0*zoomz
        goto, plot_it
     end

     6: begin
        zoomz = zoomz/2.0
        goto, plot_it
     end

     7: begin
        yrange_style=yrange_style+1
        goto, plot_it
     end


     8: begin
        i_units = 1-i_units
        goto, plot_it
     end

     9: begin
        log = (log+1) mod 2
        goto, plot_it
     end

     10: begin
      energy_kx_kyd=fltarr(2*pmax-1,2*nmax-1,ntmax/output_step+1)
      energy_kx_kysat=fltarr(2*pmax-1,2*nmax-1)
      openr,lun,'energy_kx_ky.txt',/get_lun
      readf,lun,energy_kx_kyd
      free_lun,lun

      IF(Switch_sat eq 1)Then Begin
      ix1=max(Time)/2.0+satdelta_t
      ix2=max(Time)
      Endif
      print,'Time:  ix1=',ix1,'     ix2=',ix2
      xxx1=round(ix1/(output_step*tstep))
      xxx2=round(ix2/(output_step*tstep))
      xxx1_2=xxx2-xxx1+1
      IF((ix1 EQ 0) && (ix2 EQ 0))Then Begin
      energy_kx_kysat(*,*)=energy_kx_kyd(*,*,0)
      Endif Else Begin
      energy_kx_kysat(*,*)=total(energy_kx_kyd(*,*,xxx1:xxx2),3)/xxx1_2
      Endelse

;   ;;write out the energy data for gnuplot to plot the energy spectrum
;    openw,lun,'energy_kx_kysat.txt',/get_lun
;         For p=0,2*pmax-2 Do Begin
;               For n=0,2*nmax-2 Do Begin
;               Printf,lun,kx(p),ky(n),energy_kx_kysat(p,n)
;               Endfor
;               Printf,lun,''
;         Endfor
;    free_lun,lun
;
;    openw,lun,'energy_k_sat.txt',/get_lun
;         For p=0,2*pmax-2 Do Begin
;               For n=0,2*nmax-2 Do Begin
;               Printf,lun,energy_kx_kysat(p,n)
;               Endfor
;         Endfor
;    free_lun,lun
;
;   ;;write out the energy data for gnuplot to plot the energy growth rate spectrum of NL term
;      Enegrowthrate_of_NLd=fltarr(2*pmax-1,2*nmax-1,ntmax/output_step+1)
;      Enegrowthrate_of_NLsat=fltarr(2*pmax-1,2*nmax-1)
;      openr,lun,'Enegrowthrate_of_NL.txt',/get_lun
;      readf,lun,Enegrowthrate_of_NLd
;      free_lun,lun
;
;      IF((ix1 EQ 0) && (ix2 EQ 0))Then Begin
;      Enegrowthrate_of_NLsat(*,*)=Enegrowthrate_of_NLd(*,*,0)
;      Endif Else Begin
;      Enegrowthrate_of_NLsat(*,*)=total(Enegrowthrate_of_NLd(*,*,xxx1:xxx2),3)/xxx1_2
;      Endelse
;
;    openw,lun,'Enegrowthrate_of_NLsat.txt',/get_lun
;         For p=0,2*pmax-2 Do Begin
;               For n=0,2*nmax-2 Do Begin
;               Printf,lun,kx(p),ky(n),Enegrowthrate_of_NLsat(p,n)
;               Endfor
;               Printf,lun,''
;         Endfor
;    free_lun,lun

    energy_kk=fltarr((pmax-1)^2+(nmax-1)^2+1)
    kk_ngrid=fltarr((pmax-1)^2+(nmax-1)^2+1)
    kmax=0
    For p=0,2*pmax-2 Do Begin
       For n=0,2*nmax-2 Do Begin
       kk=(p-(pmax-1))^2+(n-(nmax-1))^2
       if(kk ne 0)then begin
       kk_ngrid(kk)=kk_ngrid(kk)+1
       energy_kk(kk)=energy_kk(kk)+energy_kx_kysat(p,n)
       endif
    Endfor
    Endfor

    For kk=1,(pmax-1)^2+(nmax-1)^2 Do begin
    if(kk_ngrid(kk) ne 0)then begin
     kmax=kmax+1
    endif
    Endfor

    k=fltarr(kmax)
    energy_k=fltarr(kmax)

    i=-1
    For kk=1,(pmax-1)^2+(nmax-1)^2 Do begin
    if(kk_ngrid(kk)ne 0)then begin
      i=i+1
      k(i)=sqrt(kk)*kxmax/pmax
      energy_kk(kk)=(energy_kk(kk)/kk_ngrid(kk))*2*3.1415926536*k(i)
      energy_k(i)=energy_kk(kk)
    endif
    Endfor

    xmax=max(k)
    xmin=min(k)
    ymax=max(energy_k)
    ymin=min(energy_k)

  !plotwindowsid=10
  window,!plotwindowsid, TITLE='Energy spectrum', xsize=650,ysize=620
  !noeras=0
  !p.color=color_value(ncolor+1)
  set_viewport,0.15,0.95,0.15,0.9
  set_xy,xmin,xmax,ymin,ymax
  !p.color=color_value(ncolor+1)

  plot,[0],[0],$
       /nodata,$
       xstyle=1,$
       xminor=0,$
       xrange=[xmin,xmax],$
;       xrange=[min(x),max(x)],$
       xtitle='k',$
       ystyle=1,$
       yminor=0,$
     yrange=yrangetemp,$
       ytitle='E(k)', $
       title='Energy(k)  (time='+strtrim(ix1,2)+'s-'+strtrim(ix2,2)+'s)',$
       /xlog,/ylog
   oplot,k,energy_k



      Phi_kd=complexarr(2*pmax-1,2*nmax-1,ntmax/output_step+1)
      phi_x=complexarr(2*pmax-1,2*nmax-1)
      Phi_xtsat=fltarr(xxx2-xxx1)
      Phi_xfsat=fltarr(xxx2-xxx1)
      openr,lun,'Phi_k.txt',/get_lun
      readf,lun,Phi_kd
      free_lun,lun
      freq=findgen(xxx2-xxx1)
      freq=freq/((xxx2-xxx1)*output_step*tstep)
      IF((ix1 EQ 0) && (ix2 EQ 0))Then Begin

      Endif Else Begin
      For i=xxx1,xxx2-1 DO begin
      phi_x=FFT(Phi_kd(*,*,i),1)
      Phi_xtsat(i-xxx1)=phi_x(pmax,nmax)^2
; a = FFT(a, -1, /OVERWRITE)
      Endfor
      Phi_xfsat=FFT(Phi_xtsat,-1)
      xmax=max(freq(0:(xxx2-xxx1)/2))
      xmin=min(freq(0:(xxx2-xxx1)/2))
      ymax=max(Phi_xfsat)
      ymin=min(Phi_xfsat)
      !plotwindowsid=11
  window,!plotwindowsid+1, TITLE='frequency spectrum of Potential  (time='+strtrim(ix1,2)+'s-'+strtrim(ix2,2)+'s)', xsize=650,ysize=600
  !noeras=0
  !p.color=color_value(ncolor+1)
  set_viewport,0.15,0.95,0.15,0.9
  set_xy,xmin,xmax,ymin,ymax
  !p.color=color_value(ncolor+1)

      plot,[0],[0],$
       /nodata,$
       xstyle=1,$
       xminor=0,$
       xrange=[xmin,xmax],$
;       xrange=[min(x),max(x)],$
       xtitle='frequency',$
       ystyle=1,$
       yminor=0,$
     yrange=yrangetemp,$
       ytitle='|Phi| ^2', $
       title='frequency spectrum of Potential',$
       /ylog;,/xlog
      oplot,freq(0:(xxx2-xxx1)/2),Phi_xfsat(0:(xxx2-xxx1)/2)

   i=pmax-1
   j=nmax+2
   x=Time
   xmin=min(x)
   xmax=max(x)
   ymin=min(REAL_PART(Phi_kd(i,j,*)))
   ymax=max(REAL_PART(Phi_kd(i,j,*)))
  window,!plotwindowsid+2, TITLE='Phi_Re(kx='+strtrim(i*kxmax/(pmax-1)-kxmax,2)+'ky='+strtrim(j*kymax/(nmax-1)-kymax,2)+')', xsize=650,ysize=550
  !noeras=0
  !p.color=color_value(ncolor+1)
  set_viewport,0.2,0.95,0.15,0.9
  set_xy,xmin,xmax,ymin,ymax
  !p.color=color_value(ncolor+1)

      plot,[0],[0],$
       /nodata,$
       xstyle=1,$
       xminor=0,$
       xrange=[xmin,xmax],$
       xtitle='time',$
       ystyle=1,$
       yminor=0,$
     yrange=[ymin,ymax],$
       ytitle='Phi_Re', $
       title='Phi_Re(kx='+strtrim(i*kxmax/(pmax-1)-kxmax,2)+'ky='+strtrim(j*kymax/(nmax-1)-kymax,2)+')';,$
       ;/ylog;,/xlog
      oplot,x,REAL_PART(Phi_kd(i,j,*))

   ymin=min(abs(Phi_kd(i,j,*)))
   ymax=max(abs(Phi_kd(i,j,*)))
  window,!plotwindowsid+3, TITLE='|Phi|(kx='+strtrim(i*kxmax/(pmax-1)-kxmax,2)+'ky='+strtrim(j*kymax/(nmax-1)-kymax,2)+')', xsize=650,ysize=550
  !noeras=0
  !p.color=color_value(ncolor+1)
  set_viewport,0.2,0.95,0.15,0.9
  set_xy,xmin,xmax,ymin,ymax
  !p.color=color_value(ncolor+1)

      plot,[0],[0],$
       /nodata,$
       xstyle=1,$
       xminor=0,$
       xrange=[xmin,xmax],$
       xtitle='time',$
       ystyle=1,$
       yminor=0,$
     yrange=[ymin,ymax],$
       ytitle='|Phi|', $
       title='|Phi|(kx='+strtrim(i*kxmax/(pmax-1)-kxmax,2)+'ky='+strtrim(j*kymax/(nmax-1)-kymax,2)+')',$
       /ylog;,/xlog
      oplot,x,abs(Phi_kd(i,j,*))

      xxx1=round(10/(output_step*tstep))
      xxx2=round(50/(output_step*tstep))
      print,'gamma=',(alog(abs(Phi_kd(i,j,xxx2)))-alog(abs(Phi_kd(i,j,xxx1))))/(50-10)
      Endelse

     end

     11: begin
     IF(!plotwindowsid eq 10)then begin
     Wdelete,!plotwindowsid
     !plotwindowsid=0
     Endif
     widget_control, event.top, /destroy
     end

     20: begin
        i_ptype = 0
        goto, plot_it
     end

     21: begin
        i_ptype = 1
        goto, plot_it
     end

  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

IF(strcmp(name[0],'Energy(DW, ZF, GAM)',19)) THEN BEGIN

      energy_kx_kyd=fltarr(2*pmax-1,2*nmax-1,ntmax/output_step+1)
      openr,lun,'energy_kx_ky.txt',/get_lun
      readf,lun,energy_kx_kyd
      free_lun,lun

      ydave=fltarr(ntmax/output_step+1,3)  ;;ydave(*,0)-DW, ydave(*,1)-ZF, ydave(*,2)-GAM
      ydave(*,1)=total(energy_kx_kyd(*,nmax-1,*),1)
      ydave(*,0)=total(total(energy_kx_kyd(*,*,*),1),1)

      ;openw,lun,'total_energyIDL.txt',/get_lun
      ;printf,lun,ydave(*,0)
      ;free_lun,lun

      ydave2=fltarr(ntmax/output_step+1)
      ydave(*,0)=ydave(*,0)-ydave(*,1)
      openr,lun,name[2],/get_lun
      readf,lun,ydave2
      free_lun,lun
      ydave(*,2)=ydave2

      x=time
      !mtitle=name[0]

       xmax=max(x)
       xmin=min(x)
       ymax=max(ydave)
       ymin=min(ydave)
      For i=0,2 Do Begin
        y=ydave(*,i)
        overplot,x,y,'time (a/cs)','',xmax,xmin,ymax,ymin,0,log,0,i,0,zoomz,yrange_style,1
      Endfor
      lines = indgen(3)                 ; for line styles
      lines=0
      col=fltarr(3)
      col=[color_value(ncolor/8+1),color_value(ncolor/4+1),color_value(ncolor*11/16+1)]
      items =['DW','ZF','GAM']        ; annotations
      sym = [0]
      legend,items,linestyle=lines,colors=col,charsize=1.5              ; vertical legend---upper left

Endif Else IF(strcmp(name[0],'F_GAM(kx=1)',11)) THEN BEGIN
      F_GAM=fltarr(3,ntmax/output_step+1)
      openr,lun,'F_GAM_1.txt',/get_lun
      readf,lun,F_GAM
      free_lun,lun

      x=time
      !mtitle=name[0]

       xmax=max(x)
       xmin=min(x)
       ymax=max(F_GAM)
       ymin=min(F_GAM)
      For i=0,2 Do Begin
        y=F_GAM(i,*)
        overplot,x,y,'time (a/cs)','',xmax,xmin,ymax,ymin,0,log,0,i,0,zoomz,yrange_style,1
      Endfor
      lines = indgen(3)                 ; for line styles
      lines=0
      col=fltarr(3)
      col=[color_value(ncolor/8+1),color_value(ncolor/4+1),color_value(ncolor*11/16+1)]
      items =['Re(F_GAM)','Im(F_GAM)','| F_GAM |']        ; annotations
      sym = [0]
      legend,items,linestyle=lines,colors=col,charsize=1.5              ; vertical legend---upper left

Endif Else IF(strcmp(name[0],'total '+'|$\Phi$| ',14) && strcmp(name[2],'Phi_k',5)) THEN BEGIN
      jump=10
      IF(ntmax/output_step GE 1000)then begin
      jump=100
      Endif Else IF(ntmax/output_step GE 10000)then begin
      jump=1000
      Endif Else IF(ntmax/output_step GE 100000)then begin
      jump=10000
      Endif Else IF(ntmax/output_step GE 1000000)then begin
      jump=100000
      Endif

      y=fltarr(ntmax/(output_step*jump)+1)
      Phi_kftnew=complexarr(2*pmax-1,2*nmax-1,ntmax/(output_step*jump)+1)
      Phi_ktmp=complexarr(2*pmax-1,2*nmax-1)

      openr,lun,name[2],/get_lun
      readf,lun,Phi_ktmp
      POINT_LUN, -lun, pos
      Print,'Pos=',pos
      For i=0LL,ntmax/output_step,jump Do begin
      POINT_LUN, lun, pos*i
      readf,lun,Phi_ktmp
      Phi_kftnew(*,*,i/jump)=Phi_ktmp
      Endfor
      free_lun,lun

      ;y(*)=total(total(abs(Phi_kftnew(*,*,*)),1),1)/Omega_star
      y(*)=sqrt(total(total(conj(Phi_kftnew(*,*,*))*Phi_kftnew(*,*,*),1),1)/(Omega_star^2))
      x=findgen(ntmax/(output_step*jump)+1)*(tstep*output_step*jump)
      !mtitle=name[0]
   plot1line,x,y,'time (a/c'+'$\downs$'+')',name[1],0,(output_step*jump),ntmax,tstep,satdelta_t,i_ptype,log,zoomz,yrange_style ;plot1line used to plot total energy and total D

   
Endif Else IF(strcmp(name[0],'Data processing',15) && strcmp(name[2],' .txt',5)) THEN BEGIN

;; quicktrigger=1; for Chi=conjg(Phi)*G freq spectrum. and output "y.txt"  [Important]!!!
;; quicktrigger=2; quick output E=0.5(k|phi|)^2 vs time in the file 'Ecy.txt'.
;; quicktrigger=3; quick output E=0.5(k|phi|)^2 vs omega in the file 'ykphisqr.txt'.
;; quicktrigger=4; quick output 'yRatio.txt' for omegaNL/Omegai vs time Finally for it vs omega_star 

;;-------------------------------------------------------------
  if(quicktrigger eq 1)then begin
   jump=1  ;; keep the jump=1!  If jump>1 the freq axies maybe expand!!!  
   ibegin=0LL
   IF(Stopnt EQ 0)then begin
   Nsat=(ntmax-begintime/tstep)/output_step     
   ibegin=round(begintime/tstep/output_step)
   ENDIF   
   IF(Stopnt GT 0)then begin
   Nsat=(ntmax-Stopnt)/output_step  
   ibegin=0  
   ENDIF
   
    print,'begintime=',begintime
    print,'endtime=',ntmax*tstep   
   
    y=fltarr(Nsat/jump)
    x=FINDGEN(Nsat/jump)
    x=(x-Nsat/jump/2+1)/(Nsat*tstep*output_step)*2*!PI
    Phik=complexarr(2*pmax-1,2*nmax-1,Nsat/jump)
    PhiOMG=complexarr(2*pmax-1,2*nmax-1,Nsat/jump)
    GOMG=complexarr(2*pmax-1,2*nmax-1,Nsat/jump)
    ChiOMGcom=complexarr(Nsat/jump) 
    ChiOMG=fltarr(Nsat/jump) 
    Phi_ktmp=complexarr(2*pmax-1,2*nmax-1)   
    
    print,''  
    print,'Start working: (quicktrigger=1) ...'
    print,'For outputing Chi=conjg(Phi)*G in freq spectrum.'
    
      pos=0LL
      openr,lun,'Phi_kOMG.txt',/get_lun
      readf,lun,Phi_ktmp
      POINT_LUN, -lun, pos
      Print,'Pos=',pos
      For i=0LL,Nsat-1,jump Do begin
      POINT_LUN, lun, pos*(i+ibegin)  ;;Pay attention that i and ibegin must be interger.
      readf,lun,Phi_ktmp
      Phik(*,*,i/jump)=Phi_ktmp
      Endfor
      free_lun,lun
      
      openr,lun,'GOMG.txt',/get_lun
      readf,lun,Phi_ktmp
      POINT_LUN, -lun, pos
      Print,'Pos=',pos
      For i=0LL,Nsat-1,jump Do begin
      POINT_LUN, lun, pos*(i+ibegin)
      readf,lun,Phi_ktmp
      GOMG(*,*,i/jump)=Phi_ktmp
      Endfor
      free_lun,lun          
          
;      openr,lun,'PhiOMG.txt',/get_lun
;      readf,lun,PhiOMG
;      free_lun,lun
      
;      openr,lun,'GOMG.txt',/get_lun
;      readf,lun,GOMG
;      free_lun,lun       
      
         
       For n=0,2*nmax-2 Do Begin
        For p=0,2*pmax-2 Do Begin
        Phik(p,n,*)=FFT(Phik(p,n,*), DIMENSION=3)
        PhiOMG(p,n,*)=conj(Phik(p,n,*))*complex(0,1)*ky(n)/a_LTi
        GOMG(p,n,*)=FFT(GOMG(p,n,*), DIMENSION=3)
        Endfor
       Endfor
      
     ChiOMGcom(*)=total(total(PhiOMG(*,*,*)*GOMG(*,*,*),1),1)
     ChiOMG(*)=REAL_PART(ChiOMGcom(*))
          
          
      y(0:Nsat/jump/2-2)=ChiOMG(Nsat/jump/2+1:Nsat/jump-1)
      y(Nsat/jump/2-1:Nsat/jump-1)=ChiOMG(0:Nsat/jump/2)    

      openw,lun,'y.txt',/get_lun
      For j=0L,Nsat/jump-1 Do begin
      printf,lun,y(j)
      Endfor
      free_lun,lun      
      print,'Finish outputing y.txt !'
      print,' '
  endif
;;---------------------------------------------------------
  if(quicktrigger eq 2)then begin

jump=1  ;;jump can only be 1 !!
      
;      ntmax=1400000
;      output_step=20
;      pmax=11
;      nmax=11
;      Omega_star=10
;      kxmax=1.5
;      kymax=1.5
;jump=100
      
      kx=findgen(2*pmax-1)*kxmax/(pmax-1)-kxmax
      ky=findgen(2*nmax-1)*kymax/(nmax-1)-kymax
  
      y=fltarr(ntmax/(output_step*jump)+1)
      Phi_kftnew=complexarr(2*pmax-1,2*nmax-1,ntmax/(output_step)+1)
      ;Phi_ktmp=complexarr(2*pmax-1,2*nmax-1)
      
      print,''  
      print,'Start working: (quicktrigger=2) ...'
      print,'For outputing electrostatic energy E vs time.'
      
      openr,lun,'Phi_kcy.txt',/get_lun
      readf,lun,Phi_kftnew
      free_lun,lun
        
      y(*)=0
      For i=0LL,ntmax/output_step, jump Do begin
        For p=0,2*pmax-2  Do Begin
        For n=0,2*nmax-2  Do Begin
        y(i/jump)=y(i/jump)+(kx(p)^2+ky(n)^2)*conj(Phi_kftnew(p,n,i))*Phi_kftnew(p,n,i)/(Omega_star^2)/2.0
        Endfor
        Endfor
      Endfor
      
      openw,lun,'Ecy.txt',/get_lun
      
      For i=0LL,ntmax/(output_step*jump) Do begin
      printf,lun,y(i)
      Endfor
      free_lun,lun
      print,'Finish outputing: Ecy.txt'
      print,''  
  endif
;;---------------------------------------------------------
  if(quicktrigger eq 3)then begin
  
 ;; ZF=0 for output the sum over all k result to 'ykphisqr.txt'
 ;; ZF=1 for output ky=0 result to 'yphiZF.txt', ky/=0 result to 'yphiNZF.txt'
  
  ZF=0
  jump=1  ;;jump can only be 1 !!
  
  
   ibegin=0LL
   IF(Stopnt EQ 0)then begin
   Nsat=(ntmax-begintime/tstep)/output_step     
   ibegin=round(begintime/tstep/output_step)
   ENDIF   
   IF(Stopnt GT 0)then begin
   Nsat=(ntmax-Stopnt)/output_step  
   ibegin=0  
   ENDIF
   
    print,'begintime=',begintime
    print,'endtime=',ntmax*tstep   
         
    y=fltarr(Nsat/jump)
    x=FINDGEN(Nsat/jump)
    x=(x-Nsat/jump/2+1)/(Nsat*tstep*output_step)*2*!PI
    Phi_k=complexarr(2*pmax-1,2*nmax-1,Nsat/jump)
    PhiOMG=fltarr(Nsat/jump)   
    Phi_ktmp=complexarr(2*pmax-1,2*nmax-1)   
      
    kx=findgen(2*pmax-1)*kxmax/(pmax-1)-kxmax
    ky=findgen(2*nmax-1)*kymax/(nmax-1)-kymax
    
    print,''        
    print,'Start working: (quicktrigger=3) ...'  
    print,'For outputing electrostatic energy E in frequency spectrum.'
          
      pos=0LL
      openr,lun,'Phi_kcy.txt',/get_lun
      readf,lun,Phi_ktmp
      POINT_LUN, -lun, pos
      Print,'Pos=',pos
      For i=0LL,Nsat-1,jump Do begin
      POINT_LUN, lun, pos*(i+ibegin) ;;Pay attention that i and ibegin must be interger.
      readf,lun,Phi_ktmp
      Phi_k(*,*,i/jump)=Phi_ktmp
      Endfor
      free_lun,lun
      
       For n=0,2*nmax-2 Do Begin
        For p=0,2*pmax-2 Do Begin
        Phi_k(p,n,*)=FFT(Phi_k(p,n,*), DIMENSION=3)
        Endfor
      Endfor
      
      
   ;----------------------------------- 
    IF(ZF eq 0)then begin
      PhiOMG(*)=0
      For i=0LL,Nsat-1,jump Do begin
        For p=0,2*pmax-2  Do Begin
        For n=0,2*nmax-2  Do Begin
        PhiOMG(i)=PhiOMG(i)+(kx(p)^2+ky(n)^2)*conj(Phi_k(p,n,i))*Phi_k(p,n,i)/(Omega_star^2)/2.0
        Endfor
        Endfor
      Endfor
      
      y(0:Nsat/jump/2-2)=PhiOMG(Nsat/jump/2+1:Nsat/jump-1)
      y(Nsat/jump/2-1:Nsat/jump-1)=PhiOMG(0:Nsat/jump/2)    

      openw,lun,'ykphisqr.txt',/get_lun
      For j=0L,Nsat/jump-1 Do begin
      printf,lun,y(j)
      Endfor
      free_lun,lun
      print,'Finish outputing: ykphisqr.txt'
      print,''  
    Endif
   ;------------------------------------
    IF(ZF eq 1)then begin
      PhiZF=fltarr(Nsat/jump) 
      PhiNZF=fltarr(Nsat/jump) 
    
      PhiZF(*)=0
      PhiNZF(*)=0
      print, 'ky=0:', ky(nmax-1)
      For i=0LL,Nsat-1,jump Do begin
      For n=0,2*nmax-2  Do Begin
      For p=0,2*pmax-2  Do Begin
        IF(n eq nmax-1)then begin
        PhiZF(i)=PhiZF(i)+(kx(p)^2+ky(n)^2)*conj(Phi_k(p,n,i))*Phi_k(p,n,i)/(Omega_star^2)/2.0
        Endif Else Begin
        PhiNZF(i)=PhiNZF(i)+(kx(p)^2+ky(n)^2)*conj(Phi_k(p,n,i))*Phi_k(p,n,i)/(Omega_star^2)/2.0
        Endelse
      Endfor  
      Endfor
      Endfor
      
      
      y(0:Nsat/jump/2-2)=PhiZF(Nsat/jump/2+1:Nsat/jump-1)
      y(Nsat/jump/2-1:Nsat/jump-1)=PhiZF(0:Nsat/jump/2) 
      openw,lun,'yphiZF.txt',/get_lun
      For j=0L,Nsat/jump-1 Do begin
      printf,lun,y(j)
      Endfor
      free_lun,lun
      print,'Finish outputing: yphiZF.txt'
      
      y(*)=0         
      y(0:Nsat/jump/2-2)=PhiNZF(Nsat/jump/2+1:Nsat/jump-1)
      y(Nsat/jump/2-1:Nsat/jump-1)=PhiNZF(0:Nsat/jump/2) 
      openw,lun,'yphiNZF.txt',/get_lun
      For j=0L,Nsat/jump-1 Do begin
      printf,lun,y(j)
      Endfor
      free_lun,lun
      print,'Finish outputing: yphiNZF.txt'      
      print,''  
    Endif    
    ;------------------------------------ 
      
  endif  
;;--------------------------------------------------------------------
  if(quicktrigger eq 4)then begin
      
      ;jump=5
      
      kx=findgen(2*pmax-1)*kxmax/(pmax-1)-kxmax
      ky=findgen(2*nmax-1)*kymax/(nmax-1)-kymax
      k1x=findgen(2*pmax-1)*kxmax/(pmax-1)-kxmax
      k1y=findgen(2*nmax-1)*kymax/(nmax-1)-kymax
        
      y=fltarr(ntmax/(output_step*jump)+1)
      Phi_k=complexarr(2*pmax-1,2*nmax-1,ntmax/output_step+1)
      ;Phi_ktmp=complexarr(2*pmax-1,2*nmax-1)
      
      print,''  
      print,'Start working: (quicktrigger=4) ...'
      print,'For outputing frequency ratio omegaNL/Omegai vs time.'
       
      openr,lun,'Phi_k.txt',/get_lun
      readf,lun,Phi_k
      free_lun,lun
        
     For i=0,ntmax/output_step, jump Do begin   
      For p=0,2*pmax-2  Do Begin  
      For n=0,2*nmax-2  Do Begin  
        For p1=0,2*pmax-2  Do Begin  
        For n1=0,2*nmax-2  Do Begin  
              
        y(i/jump)=y(i/jump)+abs(REAL_PART(complex(0,-1)*(kx(p)*k1y(n1)-k1x(p1)*ky(n))*Phi_k(p1,n1,i)/(2*pmax-1)/(2*nmax-1)))
        
        ;(k1x(p1)*k2y(n2)-k2x(p2)*k1y(n1))
         Endfor
         Endfor
       Endfor
       Endfor
      Endfor
                   
      openw,lun,'yRatio.txt',/get_lun
      
      For i=0LL,ntmax/(output_step*jump) Do begin
      printf,lun,y(i)
      Endfor
      free_lun,lun
      print,'Finish outputing: yRatio.txt'  
      print,''      
      
  endif  
  
;;--------------------------------------------------------------------  
    if(quicktrigger eq 5)then begin
      
      
      kx=findgen(2*pmax-1)*kxmax/(pmax-1)-kxmax
      ky=findgen(2*nmax-1)*kymax/(nmax-1)-kymax
      k1x=findgen(2*pmax-1)*kxmax/(pmax-1)-kxmax
      k1y=findgen(2*nmax-1)*kymax/(nmax-1)-kymax
        
      
      print,''  
      print,'jump=',jump  
      print,'Start working: (quicktrigger=5) ...'
      print,'For outputing frequency ratio omegaNL/Omegai vs time.'
       
        
      Phi_knew=complexarr(2*pmax-1,2*nmax-1,ntmax/(output_step*jump)+1)
      Phi_ktmp=complexarr(2*pmax-1,2*nmax-1)
      Phisqrsat=fltarr(2*pmax-1,2*nmax-1)
      freqRatio=fltarr(2*pmax-1,2*nmax-1)
      

      openr,lun,'Phi_k.txt',/get_lun
      readf,lun,Phi_ktmp
      POINT_LUN, -lun, pos
      Print,'Pos=',pos
      For i=0L,ntmax/output_step,jump Do begin
      POINT_LUN, lun, pos*i
      readf,lun,Phi_ktmp
      Phi_knew(*,*,i/jump)=Phi_ktmp
      Endfor
      free_lun,lun

      Phisqrsat(*,*)=complex(0,0)
      For j=round(sattime/(output_step*jump*tstep)),round(ntmax/(output_step*jump)) Do Begin
      Phisqrsat(*,*)=Phisqrsat(*,*)+Phi_knew(*,*,j)*conj(Phi_knew(*,*,j))
      Endfor
      Phisqrsat(*,*)=Phisqrsat(*,*)*(1.0/(ntmax/(output_step*jump)-sattime/(output_step*jump*tstep)))


      freqRatio(*,*)=0 
      For p=0,2*pmax-2  Do Begin  
      For n=0,2*nmax-2  Do Begin  
        For p1=0,2*pmax-2  Do Begin  
        For n1=0,2*nmax-2  Do Begin  
          For j=round(sattime/(output_step*jump*tstep)),round(ntmax/(output_step*jump)) Do Begin
          freqRatio(p,n)=freqRatio(p,n)+abs(REAL_PART(complex(0,-1)*(kx(p)*k1y(n1)-k1x(p1)*ky(n))*Phi_knew(p1,n1,j)))
           Endfor
         Endfor
         Endfor
       Endfor
       Endfor
       freqRatio(*,*)=freqRatio(*,*)*(1.0/(ntmax/(output_step*jump)-sattime/(output_step*jump*tstep)))
     
      FinalFR=0
      For p=0,2*pmax-2  Do Begin  
      For n=0,2*nmax-2  Do Begin      
      FinalFR=FinalFR+freqRatio(p,n)*Phisqrsat(p,n)
      Endfor
      Endfor
      FinalFR=FinalFR/total(Phisqrsat(*,*))
      
    print,'' 
    print,'Final frequency ratio=' 
    print,FinalFR
    
        
                    
      openw,lun,'Phisqrsat.txt',/get_lun
      printf,lun,Phisqrsat
      free_lun,lun
      
      openw,lun,'freqRatio.txt',/get_lun
      printf,lun,freqRatio
      free_lun,lun      
      
      print,'Finish outputing: Phisqrsat.txt and freqRatio.txt'  
      print,''      
      
  endif  
;;--------------------------------------------------------------------
Endif Else IF(strcmp(name[0],'Others by time',14) && strcmp(name[2],' .txt',5)) THEN BEGIN

  
   ;; quicktrigger=0  For exploring what variable conserves during the CK transport jump
   ;; quicktrigger=2  For using 'Ecy.txt' to plot electrostatic energy E (CKinCH) vs time
   ;; quicktrigger=4  For using 'yRatio.txt' to plot omega_NL/OMG vs time
   
   
   
   if (quicktrigger eq 0)then begin
   ;;------------------------------------  
;      ni_n0=fltarr(ntmax/output_step+1)
;      openr,lun,'ni_n0.txt',/get_lun
;      readf,lun,ni_n0
;      free_lun,lun
; 
;      ne_n0=fltarr(ntmax/output_step+1)
;      openr,lun,'sqrn_n0cy.txt',/get_lun
;      readf,lun,ne_n0
;      free_lun,lun
      
      Phisks=fltarr(ntmax/output_step+1)
      openr,lun,'Phisks.txt',/get_lun
      readf,lun,Phisks
      free_lun,lun
                  
;      Phisqrcy=fltarr(ntmax/output_step+1)
;      openr,lun,'Phisqrcy.txt',/get_lun
;      readf,lun,Phisqrcy
;      free_lun,lun
                      
      y=fltarr(ntmax/output_step+1)
      
;      y=(ni_n0)
;      name[0]='($\delta$n$\downi$/n$\down0$)$\up2$ (CKinCH)'      

      y=(Phisks)
      name[0]='k$\up2$$\delta$$\phi$$\up2$ (CKinCH)' 
            
;      y=(ni_n0+Phisks)/2
;      name[0]='(ni_n0+Phisks)/2 (CKinCH)'
      
;      y=(ni_n0+ne_n0)/2
;      name[0]='(ni_n0+ne_n0)/2 (CKinCH)'
      
;      y=(ni_n0+ne_n0+Phisks)/2
;      name[0]='(ni_n0+ne_n0+Phisks)/2 (CKinCH)
      
;      y=(ne_n0+Phisqrcy)/2
;      name[0]='(ne_n0+Phisqrcy)/2 (CKinCH)
                      
      x=time
      !mtitle=name[0]
  print,'plot1line yrange_style=',yrange_style
  plot1line,x,y,'time (a/c'+'$\downs$'+')',name[1],0,output_step,ntmax,tstep,satdelta_t,i_ptype,log,zoomz,yrange_style ;plot1line used to plot total energy and total D
  

 ;;------------------------------------  
   endif else if(quicktrigger eq 2)then begin
      jump=1
     
      y=fltarr(ntmax/(output_step*jump)+1)
    
      openr,lun,'Ecy.txt',/get_lun
      readf,lun,y
      free_lun,lun 
      
      x=findgen(ntmax/(output_step*jump)+1)*(tstep*output_step*jump)
      name=['electrostatic Energy E (CKinCH)','E$\upCK$','Phi_kcy.txt']
      
      !mtitle=name[0]
      plot1line,x,y,'time (a/c'+'$\downs$'+')',name[1],0,(output_step*jump),ntmax,tstep,satdelta_t,i_ptype,log,zoomz,yrange_style    
      
      
      
   endif else if(quicktrigger eq 4)then begin  
     
      y=fltarr(ntmax/(output_step*jump)+1)
    
      openr,lun,'yRatio.txt',/get_lun
      readf,lun,y
      free_lun,lun 
      
      x=findgen(ntmax/(output_step*jump)+1)*(tstep*output_step*jump)
      name=['$\omega$$\downNL$/$\Omega$$\downi$','$\omega$$\downNL$/$\Omega$$\downi$','Phi_k.txt']
      
      !mtitle=name[0]
      plot1line,x,y,'time (a/c'+'$\downs$'+')',name[1],0,(output_step*jump),ntmax,tstep,satdelta_t,i_ptype,log,zoomz,yrange_style       
   
        
   endif
;;---------------------------------------------------------------  

Endif Else Begin

      y=fltarr(ntmax/output_step+1)
      openr,lun,name[2],/get_lun
      readf,lun,y
      free_lun,lun
     ; y(0)=y(1)     ;ignore the first point of nt=0, because it's really a bad point, so I set the y(0)=y(1)
      x=time
      !mtitle=name[0]
 print,'plot1line yrange_style=',yrange_style
   plot1line,x,y,'time (a/c'+'$\downs$'+')',name[1],0,output_step,ntmax,tstep,satdelta_t,i_ptype,log,zoomz,yrange_style ;plot1line used to plot total energy and total D


;the following part is inserted to analys the energy2(mu_count)
;total_energy_zhu=fltarr(mumax,ntmax/output_step+1)
; openr,lun,'total_energy_zhu.txt',/get_lun
; readf,lun,total_energy_zhu
; free_lun,lun

;For mu_count=0,mumax-1  Do Begin
;For nt_count=ntmax/output_step-5,ntmax/output_step  Do Begin
;print,'total_energy_zhu(',mu_count+1,')',total_energy_zhu(mu_count,nt_count)
;EndFor
;EndFor


Endelse

 return
end

;*******************************************************************************


pro kspectrum1D_see,nametemp,windownum,group=group

;  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data
  ;;
  common plotaxis,zoomz,zoombackup0,zoombackup1,log,half,freq,i_ptype,satdelta_t,name,yrange_style
  common main_plot,plotidh,ntime,ndata,data_name,windowmark,basemark

  ;;-----------------------------------------------

  ;;------------------------------------------
  ;; Return conditions
  ;;
 ; if exists_diff eq 0 then return
 ; if xregistered('diffusion_ave_see') then return
  ;;------------------------------------------
  name=nametemp
  base = widget_base(title=name[0],$
                     /column)
   ;;set default value
  i_ptype=0
  log=0
  half=0
  freq=0
  zoomz = 1.0  ;zoomz is the parameter which used to adjust the range of the axis
  zoombackup0=1.0  ;zoombackup0 and zoombackup1 are used to set zoomz=1, when changes the plot between plot(re) and plot(im).
  zoombackup1=1.0
  satdelta_t=0.0  ;satdelta_t is the parameter which used to adjust the saturation time when the user click the button "t+" or "t-"
  yrange_style=1  ;yrange_style is the parameter which used to choose plot the top, middle or bottom part of the picture.
                  ;for the value of yrange_style/3, 0 for top, 1 for middle, 2 for bottom;
  defsysv,'!Omega_p',0  ;;!Omega_p--to mark the Im or Re plot of complex number plot
  defsysv,'!MS_Init',0
  defsysv,'!threeD',0
  ;i_tp = 0
  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row,$
                     /frame)

 IF(strcmp(name[0],'Omega',5) ) THEN BEGIN
  x = widget_button(row1, $
                    value='Plot(re)', $
                    uvalue=0)

  x = widget_button(row1, $
                    value='Plot(im)', $
                    uvalue=1)
 Endif Else Begin
  x = widget_button(row1, $
                    value='Plot', $
                    uvalue=0)
 Endelse

 IF(strcmp(name[0],'(CKinFH)|Phi|^2 freq spectrum',25)) THEN BEGIN
 x = widget_button(row1, $
                    value='Half', $
                    uvalue=2)
 Endif
;  x = widget_button(row1, $
;                    value='t+', $
;                    uvalue=1)
;
;  x = widget_button(row1, $
;                    value='t-', $
;                    uvalue=2)
;
  x = widget_button(row1, $
                    value='ZOOM in', $
                    uvalue=3)

  x = widget_button(row1, $
                    value='ZOOM out', $
                    uvalue=4)
  x = widget_button(row1, $
                    value='Yrange', $
                    uvalue=5)
 IF((strcmp(name[0],'Omegaft(',8) || strcmp(name[0],'Omegacy(',8)) && (!MS_Init eq 0))THEN BEGIN
  x = widget_button(row1, $
                    value='freq', $
                    uvalue=6)
 Endif

 IF(strcmp(name[0],'Omega',5)) THEN BEGIN
   x = widget_button(row1,$
                    value='MS/Init',$
                    uvalue=7)
 Endif
 
  x = widget_button(row1,$
                    value='log',$
                    uvalue=8)
;
;  x = widget_button(row1,$
;                    value='TYPE',$
;                    /menu)
;
;  tlevels=['Line Plot','PDF']
;  for i=0,1 do begin
;     x1 = widget_button(x,$
;                        value=tlevels[i],$
;                        uvalue=20+i)
;  endfor
;
;  x = widget_button(row1, $
;                    value='Cal E(k)sat', $
;                    uvalue=9)

  x = widget_button(row1, $
                    value='3D', $
                    uvalue=9)

  x = widget_button(row1, $
                    value='Done', $
                    uvalue=10)
  ;;----------------------------------------------------------
  ;; DRAW WIDGET and CONTROL
  ;;----------------------------------------------------------

  draw = widget_draw(base,     $
                     xsize=750, $
                     ysize=650)

  widget_control, base, $
    ;set_uvalue=state,$
    /no_copy, $
    /realize

  ;!plotkspecid=!D.WINDOW
    (*windowmark)[windownum]=!D.window
    basemark[windownum]=base

  xmanager,'kspectrum1D', $
    base,$
;    event='energy_time_trace_event',$
    group_leader=group


end


;*******************************************************************************
pro kspectrum1D_event,event


  common startup,number_plot,fpath,ncolor,color_value,plotid
  common plotaxis,zoomz,zoombackup0,zoombackup1,log,half,freq,i_ptype,satdelta_t,name,yrange_style
  common main_plot,plotidh,ntime,ndata,data_name,windowmark,basemark

  widget_control, event.id, $
    get_uvalue=uvalue
 ; wset, widget

openr,lun,'inputRCYCLO.txt',/get_lun
  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Control_variables',plot_name_temp,19)
  endwhile
  thisFormat='(8x,i2,/,9x,i2,/,4x,i2)'
  readf,lun,restart,GK_FK_CK,CDW,Format=thisFormat
  thisFormat='(8x,i2,/,11x,i2,/,12x,i5,/,11x,i4,/,9x,f5.3,/,/,/,/,7x,i10)'
  readf,lun,muDtype,mugridtype,output_step,backup_num,Const_nl,Stopnt,Format=thisFormat

  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Grid_variables',plot_name_temp,16)
  endwhile
  thisFormat='(6x,i10,/,6x,f12.9,/,6x,f10.3,/,6x,f10.3,/,5x,i5,/,5x,i5,/,6x,i8.3,/,5x,i5,/,5x,i6,/,5x,i6)'
  readf,lun,ntmax,tstep,kxmax,kymax,pmax,nmax,mumax,N_mu,$
      N_FT,N_CY,Format=thisFormat

  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Physics_variables',plot_name_temp,19)
  endwhile
  thisFormat='(11x,f7.4,/,5x,f5.3,/,6x,f6.3,/,9x,f5.3,/,9x,f5.3,/,9x,f5.3,/,7x,f9.3,/,8x,f5.3,/,6x,f5.3,/,6x,f5.3,/,8x,f12.7,/,8x,f8.5,/,/,/,6x,f8.4)'
  readf,lun,Omega_star,a_Ln,a_LTi,lambda_n,lambda_0,lambda_D,AlphaA,delta_1,$
  mu_HK,mu_LK,F_k_int,Epsilon,gamIC, Format=thisFormat

free_lun,lun



openr,lun,'inputvuRCYCLO.txt',/get_lun
  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Control_variables',plot_name_temp,19)
  endwhile
  thisFormat='(7x,i5,/,9x,i5,/,8x,f10.3,/,14x,i2,/,13x,i2,/,/,/,4x,i6,/,4x,i6,/,11x,i4)'
  readf,lun,avenum,plotstep,sattime,Intermit_Plot,overplot_log,ix1,ix2,Switch_sat,Format=thisFormat

  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Range_control',plot_name_temp,15)
  endwhile
  thisFormat='(10x,f15.12,/,10x,f12.4,/,11x,f15.10,/,11x,f15.12,/,9x,i5,/,8x,f10.4)'
  readf,lun,Omega_min,Omega_max,Phi_logmin,Chi_logmin,OMGpoint,freqPer,Format=thisFormat

free_lun,lun



;  openr,lun,'nt_Odd.txt',/get_lun
;  readf,lun,nt_Odd
;  free_lun,lun
;  openr,lun,'nt_Even.txt',/get_lun
;  readf,lun,nt_Even
;  free_lun,lun
;  IF(ntmax GE max([nt_Odd,nt_Even]))Then Begin
;  ntmax=max([nt_Odd,nt_Even])
;  Endif
;  print,'ntmax=',ntmax

  Time=findgen(ntmax/output_step+1)*(tstep*output_step)
  kx=findgen(2*pmax-1)*kxmax/(pmax-1)-kxmax
  ky=findgen(2*nmax-1)*kymax/(nmax-1)-kymax
  kxhalf=findgen(pmax)*kxmax/(pmax-1)
  kyhalf=findgen(nmax)*kymax/(nmax-1)

  For i=0,200-1 Do begin
  IF(basemark[i] eq event.top)then begin
    IF((*windowmark)[i] ne !D.window)then begin
    wset,(*windowmark)[i]
    Endif
    case(i) of
    4: name=['Energy(kx) (3DGK)','energy_kx_ky.txt','kx']
    5: name=['Energy(ky) (3DGK)','energy_kx_ky.txt','ky']
    6: name=['Energy (kx)int&sat (3DGK)','energy_kx_ky.txt','kx']
    7: name=['Energy (ky)int&sat (3DGK)','energy_kx_ky.txt','ky']
    8: name=['D(kx) (3DGK)','D_3G.txt','kx']
    9: name=['D(ky) (3DGK)','D_3G.txt','ky']
    10: name=['Chi(kx) (3DGK)','Chi.txt','kx']
    11: name=['Chi(ky) (3DGK)','Chi.txt','ky']
    12: name=['Chi (kx)int&sat (3DGK)','Chi.txt','kx']
    13: name=['Chi (ky)int&sat (3DGK)','Chi.txt','ky']
    14: name=['Energy(k) (3DGK)','energy_k.txt','k.txt','k']
    16: name=['('+cgGreek('omega')+') GK [Eigen]','ky','Phi_k.txt']
    17: name=['('+cgGreek('omega')+') CKinFH [Eigen]','ky','Phi_kcy.txt']
    19: name=['D (kx)int&sat (3DGK)','D_3G.txt','kx']
    20: name=['D (ky)int&sat (3DGK)','D_3G.txt','ky']
    21: name=['|Phi|^2 (kx)int&sat (3DGK)','Phi^2.txt','kx']
    22: name=['|Phi|^2 (ky)int&sat (3DGK)','Phi^2.txt','ky']
    36: name=['(CKinFH)|Phi|^2 freq spectrum','Phi_kft.txt','ky']
    37: name=['(GK)|Phi|^2 freq spectrum','Phi_k.txt','ky']
    38: name=['('+cgGreek('omega')+') CKinCH [Eigen]','ky','Phi_kcy.txt']
    102: name=['('+cgGreek('omega')+') GK [Eigen]','ky','Phi_k.txt']
    105: name=['('+cgGreek('omega')+') CKinFH [Eigen]','ky','Phi_kcy.txt']    
    108: name=['('+cgGreek('omega')+') CKinCH [Eigen]','ky','Phi_kcy.txt']    
    endcase
  Endif
  Endfor

  case (uvalue) of

;  IF(strcmp(name[0],'Omega_',6)) THEN BEGIN
     0:begin

     !Omega_p=0
     zoomz=zoombackup0
     goto, plot_it
     end

     1:begin
     !Omega_p=1
     zoomz=zoombackup1
     goto, plot_it
     end

     2: begin
        half=(half+1) mod 2
        goto, plot_it
     end

     3: begin
        zoomz = 2.0*zoomz
        IF(!Omega_p eq 0) THEN BEGIN
        zoombackup0=zoomz
        Endif Else If(!Omega_p eq 1)THEN Begin
        zoombackup1=zoomz
        Endif
        goto, plot_it
     end

     4: begin
        zoomz = zoomz/2.0
        IF(!Omega_p eq 0) THEN BEGIN
        zoombackup0=zoomz
        Endif Else If(!Omega_p eq 1)THEN Begin
        zoombackup1=zoomz
        Endif
        goto, plot_it
     end

     5: begin
        yrange_style=yrange_style+1
        goto, plot_it
     end

     6: begin
        freq=(freq+1) mod 2   ;;here using 'freq' to control the switch between Low freq and High freq
        goto, plot_it
     end

     7: begin
        !MS_Init=(!MS_Init+1) mod 2   ;;here using '!MS_Init' to control the switch between Matrix Solver result and Initial growing result
        freq=0
        goto, plot_it
     end

     8: begin
        log = (log+1) mod 3
        goto, plot_it
     end

;  Endif Else Begin
;     0:begin
;     ;wset,(*ndata)[windownum]
;     goto, plot_it
;     end
;
;  Endelse


     9: begin
        !threeD = (!threeD + 1) mod 2
        goto, plot_it
     end
     
     10: begin ;shut down this window
;     IF(!plotwindowsid eq 10)then begin
;     Wdelete,!plotwindowsid
;     !plotwindowsid=0
;     Endif
     widget_control, event.top, /destroy
     end


  endcase

  return

  plot_it:

  tempReCon=8
  
  ntnum_of_bin=round(ntmax/(output_step*avenum))
  !p.background=color_value(ncolor/2)
;=================================================================================================================
    IF(strcmp(name[0],'Energy(k)',9)) THEN BEGIN
      openr,lun,name[2],/get_lun
      readf,lun,kmax
      x=fltarr(kmax)
      readf,lun,x
      free_lun,lun

      yd=fltarr(kmax,ntmax/output_step+1)
      openr,lun,name[1],/get_lun
      readf,lun,yd
      free_lun,lun
      ydave=fltarr(kmax,avenum)
      For i=0,avenum-1 Do Begin
        ydave(*,i)=total(yd(*,ntnum_of_bin*i:ntnum_of_bin*(i+1)-1),2)
        ydave(*,i)=ydave(*,i)/ntnum_of_bin
      Endfor
      !mtitle=name[0]
       xmax=max(x)
       xmin=min(x)
       ymax=max(ydave)
       ymin=min(ydave)
      For i=0,avenum-1 Do Begin
        y=ydave(*,i)
        overplot,x,y,name[3],'',xmax,xmin,ymax,ymin,0,log,i,i,0,zoomz,yrange_style,0
      Endfor
      lines = indgen(avenum)                 ; for line styles
      linetemp1 = findgen(avenum)*(ntmax*tstep/avenum)
      linetemp2 = findgen(avenum)*(1.0*ntmax*tstep/avenum)$
        +(1.0*ntmax*tstep/avenum)
      ;line1=Uint(linetemp1)
      ;line2=Uint(linetemp2)
      line1=linetemp1
      line2=linetemp2
      items = 'time='+strtrim(line1,2)+'s-'+strtrim(line2,2)+'s'         ; annotations
      sym = [0]
      legend,items,linestyle=lines,charsize=1.5              ; vertical legend---upper left
    EndIf

;=================================================================================================================

    IF(strcmp(name[0],'Energy(kx)',10) || strcmp(name[0],'D(kx)',5) || strcmp(name[0],'Chi(kx)',7)) THEN BEGIN
      yd=fltarr(2*pmax-1,2*nmax-1,ntmax/output_step+1)
      openr,lun,name[1],/get_lun
      readf,lun,yd
      free_lun,lun
      ydave=fltarr(pmax,2*nmax-1,avenum)
      For i=0,avenum-1 Do Begin
        ydave(*,*,i)=total(yd(pmax-1:2*pmax-2,*,ntnum_of_bin*i:ntnum_of_bin*(i+1)-1),3)
        ydave(*,*,i)=ydave(*,*,i)/ntnum_of_bin
      Endfor
      x=kxhalf
      !mtitle=name[0]
      sumky_yave=fltarr(pmax,avenum)
      For i=0,avenum-1 Do Begin
      sumky_yave(*,i)=TOTAL(ydave(*,*,i),2)/(2*nmax-1)
      Endfor
       xmax=max(x)
       xmin=min(x)
       ymax=max(ydave)
       ymin=min(ydave)
      For i=0,avenum-1 Do Begin
        y=sumky_yave(*,i)
        overplot,x,y,name[2],'',xmax,xmin,ymax,ymin,0,log,i,i,0,zoomz,yrange_style,0;,xmax,xmin,ymax,ymin,0,overplot_log,i,i,0,1.0,1,0,
      Endfor
      lines = indgen(avenum)                 ; for line styles
      linetemp1 = findgen(avenum)*(ntmax*tstep/avenum)
      linetemp2 = findgen(avenum)*(1.0*ntmax*tstep/avenum)$
        +(1.0*ntmax*tstep/avenum)
      ;line1=Uint(linetemp1)
      ;line2=Uint(linetemp2)
      line1=linetemp1
      line2=linetemp2
      items = 'time='+strtrim(line1,2)+'s-'+strtrim(line2,2)+'s'         ; annotations
      sym = [0]
      legend,items,linestyle=lines,charsize=1.5              ; vertical legend---upper left
    EndIf

;=================================================================================================================

    IF(strcmp(name[0],'Energy(ky)',10) || strcmp(name[0],'D(ky)',5) || strcmp(name[0],'Chi(ky)',7)) THEN BEGIN
      yd=fltarr(2*pmax-1,2*nmax-1,ntmax/output_step+1)
      openr,lun,name[1],/get_lun
      readf,lun,yd
      free_lun,lun
      ydave=fltarr(2*pmax-1,nmax,avenum)
      For i=0,avenum-1 Do Begin
        ydave(*,*,i)=total(yd(*,nmax-1:2*nmax-2,ntnum_of_bin*i:ntnum_of_bin*(i+1)-1),3)
        ydave(*,*,i)=ydave(*,*,i)/ntnum_of_bin
      Endfor
      x=kyhalf
      !mtitle=name[0]
      sumkx_yave=fltarr(nmax,avenum)
      For i=0,avenum-1 Do Begin
      sumkx_yave(*,i)=TOTAL(ydave(*,*,i),1)/(2*pmax-1)
      Endfor
       xmax=max(x)
       xmin=min(x)
       ymax=max(ydave)
       ymin=min(ydave)
      For i=0,avenum-1 Do Begin
        y=sumkx_yave(*,i)
        overplot,x,y,name[2],'',xmax,xmin,ymax,ymin,0,log,i,i,0,zoomz,yrange_style,0;,xmax,xmin,ymax,ymin,0,overplot_log,i,i,0,1.0,1,0
      Endfor
      lines = indgen(avenum)                 ; for line styles
      linetemp1 = findgen(avenum)*(ntmax*tstep/avenum)
      linetemp2 = findgen(avenum)*(1.0*ntmax*tstep/avenum)$
        +(1.0*ntmax*tstep/avenum)
      ;line1=Uint(linetemp1)
      ;line2=Uint(linetemp2)
      line1=linetemp1
      line2=linetemp2
      items = 'time='+strtrim(line1,2)+'s-'+strtrim(line2,2)+'s'         ; annotations
      sym = [0]
      legend,items,linestyle=lines,charsize=1.5              ; vertical legend---upper left
    Endif

;=================================================================================================================

    IF(strcmp(strmid(name[0],17,18,/reverse_offset),'(kx)int&sat (3DGK)',18)) THEN BEGIN
     ; window, plotidh+2, TITLE='main_plot', xsize=600,ysize=600
      yd=fltarr(2*pmax-1,2*nmax-1,ntmax/output_step+1)
      openr,lun,name[1],/get_lun
      readf,lun,yd
      free_lun,lun
      ydave=fltarr(pmax,2*nmax-1,2)
      ydave(*,*,0)=yd(pmax-1:2*pmax-2,*,0)  ;the initial energy

      ydave(*,*,1)=total(yd(pmax-1:2*pmax-2,*,round(sattime/(output_step*tstep)+1):round(ntmax/output_step)),3);the average energy of steady state
      ydave(*,*,1)=ydave(*,*,1)/(ntmax/output_step-sattime/(output_step*tstep))

      x=kxhalf
      !mtitle=name[0]
      sumky_yave=fltarr(pmax,2)
      For i=0,1 Do Begin
      sumky_yave(*,i)=TOTAL(ydave(*,*,i),2);/(2*nmax-1)
      Endfor
      xmax=max(x)
      xmin=min(x)
      ymax=max(sumky_yave)
      ymin=min(sumky_yave)
      IF(strcmp(name[0],'Energy (kx)int&sat',18)) THEN BEGIN
      namey='Energy'
      Endif
      IF(strcmp(name[0],'Chi (kx)int&sat',15)) THEN BEGIN
      namey='Chi'
      Endif
      IF(strcmp(name[0],'|Phi|^2 (kx)int&sat',19)) THEN BEGIN
      namey='|Phi| ^2'
      Endif
      IF(strcmp(name[0],'D (kx)int&sat',13)) THEN BEGIN
      namey='D'
      Endif
      For i=0,1 Do Begin
        y=sumky_yave(*,i)
        overplot,x,y,name[2],namey,xmax,xmin,ymax,ymin,0,log,i+1,i+1,1,zoomz,yrange_style,0
      Endfor
      lines = [1,2]                 ; for line styles
      items =['time=0','time='+strtrim(round(sattime),2)+'Ln/Cs-'+strtrim(round(ntmax*tstep),2)+'Ln/Cs']
      sym = [0]
      legend,items,linestyle=lines,charsize=1.5              ; vertical legend---upper left
    EndIf
;=================================================================================================================
    IF(strcmp(strmid(name[0],17,18,/reverse_offset),'(ky)int&sat (3DGK)',18)) THEN BEGIN
     ; window, plotidh+3, TITLE='main_plot', xsize=600,ysize=600
      yd=fltarr(2*pmax-1,2*nmax-1,ntmax/output_step+1)
      openr,lun,name[1],/get_lun
      readf,lun,yd
      free_lun,lun
      ydave=fltarr(2*pmax-1,nmax,2)
      ydave(*,*,0)= yd(*,nmax-1:2*nmax-2,0)  ;the initial energy

      ydave(*,*,1)=total(yd(*,nmax-1:2*nmax-2,round(sattime/(output_step*tstep)+1):round(ntmax/output_step)),3);the average energy of steady state
      ydave(*,*,1)=ydave(*,*,1)/(ntmax/output_step-sattime/(output_step*tstep))

      x=kxhalf
      !mtitle=name[0]
      sumkx_yave=fltarr(nmax,2)
      For i=0,1 Do Begin
      sumkx_yave(*,i)=TOTAL(ydave(*,*,i),1);/(2*pmax-1)
      Endfor
      xmax=max(x)
      xmin=min(x)
      ymax=max(sumkx_yave)
      ymin=min(sumkx_yave)
      IF(strcmp(name[0],'Energy (ky)int&sat',18)) THEN BEGIN
      namey='Energy'
      Endif
      IF(strcmp(name[0],'Chi (ky)int&sat',15)) THEN BEGIN
      namey='Chi'
      Endif
      IF(strcmp(name[0],'|Phi|^2 (ky)int&sat',19)) THEN BEGIN
      namey='|Phi| ^2'
      Endif
      IF(strcmp(name[0],'D (ky)int&sat',13)) THEN BEGIN
      namey='D'
      Endif
      For i=0,1 Do Begin
        y=sumkx_yave(*,i)
        overplot,x,y,name[2],namey,xmax,xmin,ymax,ymin,0,log,i+1,i+1,1,zoomz,yrange_style,0
      Endfor
      lines = [1,2]                 ; for line styles
      items =['time=0','time='+strtrim(round(sattime),2)+'Ln/Cs-'+strtrim(round(ntmax*tstep),2)+'Ln/Cs']
      sym = [0]
      legend,items,linestyle=lines,charsize=1.5              ; vertical legend---upper left
    EndIf

;=================================================================================================================

IF(strcmp(name[0],'('+cgGreek('omega')+') GK',10)) THEN BEGIN
tempp=name[0]
print,'tempp=',tempp
  x=kyhalf
  y=fltarr(nmax)
  Omega_matr=complexarr(pmax,nmax)
  
 IF(!MS_Init eq 0)THEN BEGIN
  Omega_matOrg=complexarr(2*pmax-1,nmax,N_mu)
  openr,lun,'Omega_matr.txt',/get_lun
  readf,lun,Omega_matOrg
  n=eof(lun)
  if n ne 1 then print,'error with file load!!!!!!'
  free_lun,lun
  For p=0,pmax-1 Do Begin
    For n=0,nmax-1 Do Begin
      Omega_matr(p,n)=Omega_matOrg(p+pmax-1,n,1)
      For i=1,N_mu-1 Do Begin
      IF(imaginary(Omega_matOrg(p+pmax-1,n,i)) GT imaginary(Omega_matr(p,n)) && (abs(imaginary(Omega_matOrg(p+pmax-1,n,i))) ne 0))then begin
      Omega_matr(p,n)=Omega_matOrg(p+pmax-1,n,i)
      Endif
      Endfor
    Endfor
  Endfor
  
  
  For p=0,pmax-1 Do Begin
    For n=1,nmax-1 Do Begin
      For i=1,N_mu-1 Do Begin
      IF(imaginary(Omega_matOrg(p+pmax-1,n,i)) eq imaginary(Omega_matr(p,n)))then begin
      print,'kx=',p*kxmax/(pmax-1),'ky=',n*kymax/(nmax-1),'i=',i 
      Endif
      Endfor
    Endfor
  Endfor  


 Endif Else IF(!MS_Init eq 1)THEN BEGIN
  Omega=complexarr(pmax,nmax)
  Omega(*,*)=complex(0,0)
      Phi_kd=complexarr(2*pmax-1,2*nmax-1,ntmax/output_step+1)
      openr,lun,'Phi_k.txt',/get_lun
      readf,lun,Phi_kd
      free_lun,lun
      
      ;+-+-
      tmax=max(Time)
      xx1=tmax*0.5
      xx2=tmax-tmax/20.0
      xxx1=round(xx1/(output_step*tstep))
      xxx2=round(xx2/(output_step*tstep))
      lxx1=xx1-tmax/20.0
      rxx1=xx1+tmax/20.0
      lxx2=xx2-tmax/20.0
      rxx2=xx2+tmax/20.0
      lxxx1=round(lxx1/(output_step*tstep))
      rxxx1=round(rxx1/(output_step*tstep))
      lxxx2=round(lxx2/(output_step*tstep))
      rxxx2=round(rxx2/(output_step*tstep))
   For p=0,pmax-1 Do Begin
    For n=0,nmax-1 Do Begin
     For i=1,ntmax/output_step Do Begin
     If((phi_kd(p+pmax-1,n+nmax-1,i)+phi_kd(p+pmax-1,n+nmax-1,i-1)) eq complex(0,0))then begin
     Omega(p,n)=Omega(p,n)+complex(0,1)*(phi_kd(p+pmax-1,n+nmax-1,i)-phi_kd(p+pmax-1,n+nmax-1,i-1))/(output_step*tstep)
     Endif Else Begin
     Omega(p,n)=Omega(p,n)+2*complex(0,1)*(phi_kd(p+pmax-1,n+nmax-1,i)-phi_kd(p+pmax-1,n+nmax-1,i-1))/(output_step*tstep*(phi_kd(p+pmax-1,n+nmax-1,i)+phi_kd(p+pmax-1,n+nmax-1,i-1)))
     Endelse
     Endfor
     Omega(p,n)=Omega(p,n)/(ntmax/output_step)
     
      y1=total(alog(abs(Phi_kd(p+pmax-1,n+nmax-1,lxxx1:rxxx1))))/(rxxx1-lxxx1+1)
      y2=total(alog(abs(Phi_kd(p+pmax-1,n+nmax-1,lxxx2:rxxx2))))/(rxxx2-lxxx2+1)
           
      Omega(p,n)=complex(REAL_PART(Omega(p,n)),(y2-y1)/(xx2-xx1))
    Endfor
   Endfor
   Omega(0,0)=complex(0,0)
      
   Omega_matr(*,*)=Omega(*,*)   
   name[0]='('+cgGreek('omega')+') GK [Init]'  
    
 Endif
 
;;-----------------------------------------------------------------
;;public code below:
  If(!Omega_p eq 0)then begin
       name[0]='Re'+name[0]
      !mtitle=name[0]
       xmax=max(x)
       xmin=min(x)
       ymax=max(REAL_PART(Omega_matr(*,*)))
       ymin=min(REAL_PART(Omega_matr(*,*)))
      For p=0,pmax-1,plotstep Do Begin
      i=p/plotstep
      y(*)=REAL_PART(Omega_matr(p,*))
      overplot,x,y,name[1],'',xmax,xmin,ymax,ymin,0,log,i,i,0,zoomz,yrange_style,0
      Endfor

   Endif Else IF(!Omega_p eq 1) THEN BEGIN
      name[0]='Im'+name[0]
      !mtitle=name[0]
       xmax=max(x)
       xmin=min(x)
       ymax=max(imaginary(Omega_matr(*,*)))
       If(min(imaginary(Omega_matr(*,*))) ge -0.5*max(imaginary(Omega_matr(*,*))))Then begin
       ymin=min(imaginary(Omega_matr(*,*)))
       Endif Else begin
       ymin=-0.5*max(imaginary(Omega_matr(*,*)))
       Endelse
       ymax=ymax+(ymax-ymin)/5.0
      For p=0,pmax-1,plotstep Do Begin
      i=p/plotstep
      y(*)=imaginary(Omega_matr(p,*))
      overplot,x,y,name[1],'',xmax,xmin,ymax,ymin,0,log,i,i,0,zoomz,yrange_style,0
      Endfor
     gammad=fltarr(pmax,nmax)
     gammad=imaginary(Omega_matr(*,*))
      temp=gammad(0,0)
      For i=0,pmax-1 Do Begin
        For j=0,nmax-1 Do Begin
          IF(gammad(i,j) NE 0) THEN BEGIN
            IF(temp LE gammad(i,j)) THEN BEGIN
              temp=gammad(i,j)
            Endif
          Endif
        Endfor
      Endfor
      For i=0,pmax-1 Do Begin
        For j=0,nmax-1 Do Begin
          IF(gammad(i,j) NE 0) THEN BEGIN
            IF(gammad(i,j) EQ temp) THEN BEGIN
            gamma_max=gammad(i,j)
            ky_max=j*kymax/(nmax-1)
              print,'kx=',i*kxmax/(pmax-1),'  ','ky=',$
                j*kymax/(nmax-1),'  ','max growth rate mode = ',Omega_matr(i,j)
            Endif
          Endif
        Endfor
      Endfor
   !p.color=color_value(ncolor*11/16)  ;;red
    xyouts,kymax/5.0,gamma_max*0.9/zoomz,'Max='+strtrim(gamma_max,2),charthick=3.5,size=4.7
   !p.color=color_value(ncolor+1)
         
      usnum=0
      For i=0,pmax-1 Do Begin
        For j=0,nmax-1 Do Begin
        IF(gammad(i,j) gt 0)Then Begin
        usnum=usnum+1
        Endif
        Endfor
      Endfor
       
   !p.color=color_value(ncolor*1/8) ;green
  ; xyouts,xmin,gamma_max*1.05/zoomz,'Unstable modes: '+strtrim(usnum,2)+'/'+strtrim(round(pmax*nmax),2),charthick=2,size=2
   !p.color=color_value(ncolor+1)
            
      stabnum=0
      For i=0,pmax-1 Do Begin
        For j=0,nmax-1 Do Begin
        IF(gammad(i,j) lt 0)Then Begin
        stabnum=stabnum+1
        Endif
        Endfor
      Endfor     
      print,''
      print,'Stable modes: '+strtrim(stabnum,2)+'/'+strtrim(round(pmax*nmax),2)   
      print,''       
      
      print,'Sum(Gamma(kx>=0,ky>=0))=',total(gammad(*,*))
      print,''
      
    ;  print,'Sum(OmegaIM(kx,ky>=0,mu))=',total(imaginary(Omega_matOrg(*,*,*)))
      print,''      
      
   Endif

      lines = indgen(5)                    ; for line styles
      line = findgen(5)*kxmax/(nmax-1)*plotstep
      items = 'kx='+strtrim(line,2)           ; annotations
      legend,items,linestyle=lines,charsize=1.5;,/bottom      ; at lower left

;      lines = indgen(nmax/plotstep)                    ; for line styles
;      line = findgen(nmax/plotstep)*kxmax/(nmax-1)*plotstep
;      items = 'kx='+strtrim(line,2)           ; annotations
;      legend,items,linestyle=lines,charsize=1.5,/right;bottom      ; at lower left
 
     
     print,"line=",line(*)
 
 IF(!threeD eq 1)THEN BEGIN
  IF(!Omega_p eq 0) THEN BEGIN
   !p.background=color_value(ncolor/2)
   window, plotidh+2, TITLE='main_plot', xsize=500,ysize=500
      x=kxhalf
      y=kyhalf
      z=abs(REAL_PART(Omega_matr(*,*)))
      !mtitle=name[0]
      surface3D,z,x,y,'kx','ky',0,1,0,Az,Ax
  Endif
  IF(!Omega_p eq 1) THEN BEGIN
   !p.background=color_value(ncolor/2)
   window, plotidh+2, TITLE='main_plot', xsize=500,ysize=500
      x=kxhalf
      y=kyhalf
      z=gammad(*,*)
      !mtitle=name[0]
      surface3D,z,x,y,'kx','ky',0,1,0,Az,Ax
  Endif
 ENDIF 
  
EndIF

;=================================================================================================================
IF(strcmp(name[0],'('+cgGreek('omega')+') CKinFH',13)) THEN BEGIN
  
  tempReconst=0.7
  x=kyhalf
  y=fltarr(nmax)
  Omegaft=complexarr(pmax,nmax)
  
 IF(!MS_Init eq 0)THEN BEGIN  
  Omega_matrftLow=complexarr(pmax,nmax) 
  Omega_matrftHigh=complexarr(pmax,nmax)
  Omega_matOrgft=complexarr(2*pmax-1,nmax,N_mu*(2*N_FT-1))
  Omega_matOrgft(*,*,*)=complex(0,0)
  openr,lun,'Omega_matrft.txt',/get_lun   ;;0 Gyro-kinetic
  readf,lun,Omega_matOrgft
  n=eof(lun)
  if n ne 1 then print,'error with file load!!!!!!'
  free_lun,lun
  
  sumOMG=0.0
 For p=0,2*pmax-2 Do Begin
  For n=0,nmax-1 Do Begin
   For i=0, N_mu*(2*N_FT-1)-1 Do Begin  
  sumOMG=sumOMG+imaginary(Omega_matOrgft(p,n,i))*1.0/(2*pmax-1)/nmax/N_mu/(2*N_FT-1)
   Endfor
  Endfor
 Endfor
  
      print,''      
      print,'Sum(OmegaFK_IM(kx,ky>=0,mu,n))=',sumOMG    ;(total(double(0.00000000001*imaginary(Omega_matOrgft(*,*,*)))))
      print,''    
;-----------------------------------------------------------------------------
  Omega_matrcyHigh=complexarr(pmax,nmax)
  Omega_matOrgcy=complexarr(2*pmax-1,nmax,N_mu*(2*N_CY-1))
  Omega_matOrgcy(*,*,*)=complex(0,0)
  openr,lun,'Omega_matrcy.txt',/get_lun
  readf,lun,Omega_matOrgcy
  n=eof(lun)
  if n ne 1 then print,'error with file load!!!!!!'
  free_lun,lun

 For p=0,pmax-1 Do Begin
  For n=0,nmax-1 Do Begin
  tempRe=abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,1)))
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
      IF(tempRe LE abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j)))) THEN BEGIN
           tempRe= abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j)))
       Endif
    Endfor
    IF(N_CY GE 2)then begin
    tempRe=tempRe/Omega_star*tempReconst/(N_CY-1) ;!!!!!!!!!!!!!!!!!! 1. 2.
    Endif

    tempRe=tempReCon ;Omega_star/2.0
    
    jtemp=-1
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j))) GE tempRe) && (abs(imaginary(Omega_matOrgcy(p+pmax-1,n,j))) ne 0)) THEN BEGIN
          jtemp=j
          growthIm=imaginary(Omega_matOrgcy(p+pmax-1,n,j))
          Endif
    Endfor
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j))) GE tempRe) && (abs(imaginary(Omega_matOrgcy(p+pmax-1,n,j))) ne 0)) THEN BEGIN
            IF(imaginary(Omega_matOrgcy(p+pmax-1,n,j)) GT growthIm) THEN BEGIN
            jtemp=j
            growthIm=imaginary(Omega_matOrgcy(p+pmax-1,n,j))
            Endif
          Endif
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrcyHigh(p,n)=complex(0,0)
      IF(p EQ pmax-1 && n Eq 0)then begin
      Omega_matrcyHigh(0,0)=Omega_matOrgcy(pmax-1,0,0)
      ENdif
    Endif else begin
    Omega_matrcyHigh(p,n)=Omega_matOrgcy(p+pmax-1,n,jtemp)
    Endelse

  Endfor
 Endfor
  ;-----------------------------------------------------------------------
 For p=0,pmax-1 Do Begin
  For n=0,nmax-1 Do Begin
  tempRe=abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,0)))
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
      IF(tempRe LE abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j)))) THEN BEGIN
           tempRe= abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j)))
       Endif
    Endfor
    IF(N_FT GE 2)then begin
    tempRe=tempRe/Omega_star*tempReconst/(N_FT-1);!!!!!!!!!!!! 3. 4. !!used as the limiter of the low freq and high freq
    Endif
 
    tempRe=tempReCon ;Omega_star/2.0
   ; print,'tempRe=',tempRe
    
    jtemp=-1
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j))) LE tempRe) && (abs(imaginary(Omega_matOrgft(p+pmax-1,n,j))) ne 0) )THEN BEGIN
          jtemp=j
          growthIm=imaginary(Omega_matOrgft(p+pmax-1,n,j))
          EndIf
    Endfor
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j))) LE tempRe) && (abs(imaginary(Omega_matOrgft(p+pmax-1,n,j))) ne 0))THEN BEGIN
           IF(growthIm LT (imaginary(Omega_matOrgft(p+pmax-1,n,j)))) THEN BEGIN
           growthIm=imaginary(Omega_matOrgft(p+pmax-1,n,j))
           jtemp=j
           Endif
          Endif
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrftLow(p,n)=complex(0,0)
    Endif else begin
    Omega_matrftLow(p,n)=Omega_matOrgft(p+pmax-1,n,jtemp)
    Endelse
        
    
    if((p gt pmax/2) and (n gt nmax/2))then begin
      For j=0,N_mu*(2*N_FT-1)-1 Do Begin
      if(Omega_matOrgft(p+pmax-1,n,j) eq Omega_matrftLow(p,n) )then begin
      print,'kx=',p*kxmax/(pmax-1),'ky=',n*kymax/(nmax-1),'i=',round(j)/(2*round(N_FT)-1) 
      endif
      Endfor
    endif

    jtemp=-1
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF(abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j))) GE tempRe) THEN BEGIN
            If(abs(abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j)))-abs(REAL_PART(Omega_matrcyHigh(p,n)))) LT freqPer*abs(REAL_PART(Omega_matrcyHigh(p,n)))) THEN BEGIN ;;**
            ;If(abs(abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j)))-OMGpoint) LT freqPer*Omega_star) THEN BEGIN
             growthIm=imaginary(Omega_matOrgft(p+pmax-1,n,j))
             jtemp=j
            Endif
          EndIf
    Endfor
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF(abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j))) GE tempRe) THEN BEGIN
             If(abs(abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j)))-abs(REAL_PART(Omega_matrcyHigh(p,n)))) LT freqPer*abs(REAL_PART(Omega_matrcyHigh(p,n)))) THEN BEGIN ;;**
             ;If(abs(abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j)))-OMGpoint) LT freqPer*Omega_star) THEN BEGIN
             IF((imaginary(Omega_matOrgft(p+pmax-1,n,j))) GT growthIm) THEN BEGIN
             growthIm=imaginary(Omega_matOrgft(p+pmax-1,n,j))
             jtemp=j
             Endif
             Endif
          EndIf
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrftHigh(p,n)=complex(0,0)
    Endif else begin
    Omega_matrftHigh(p,n)=Omega_matOrgft(p+pmax-1,n,jtemp)
    Endelse
  Endfor
 Endfor

  Omega_matrftHigh(0,0)=Omega_matOrgft(pmax-1,0,0)
 ;---------------------------------------------------------------------------------
  If(freq eq 0)then begin
  Omegaft(*,*)=Omega_matrftLow(*,*)
  nameap=' Low freq'
  Endif Else IF(freq eq 1) THEN BEGIN
  Omegaft(*,*)=Omega_matrftHigh(*,*)
  nameap=' High freq'
  Endif
  
 ;---------------------------------------------------------------------------------
 Endif Else IF(!MS_Init eq 1)THEN BEGIN
  Omegaft(*,*)=complex(0,0) 
      Phi_kftd=complexarr(2*pmax-1,2*nmax-1,ntmax/output_step+1)
      openr,lun,'Phi_kft.txt',/get_lun
      readf,lun,Phi_kftd
      free_lun,lun
        
      tmax=max(Time)
      xx1=tmax*0.5
      xx2=tmax-tmax/20.0
      xxx1=round(xx1/(output_step*tstep))
      xxx2=round(xx2/(output_step*tstep))
      lxx1=xx1-tmax/20.0
      rxx1=xx1+tmax/20.0
      lxx2=xx2-tmax/20.0
      rxx2=xx2+tmax/20.0
      lxxx1=round(lxx1/(output_step*tstep))
      rxxx1=round(rxx1/(output_step*tstep))
      lxxx2=round(lxx2/(output_step*tstep))
      rxxx2=round(rxx2/(output_step*tstep))


   For p=0,pmax-1 Do Begin
    For n=0,nmax-1 Do Begin
     For i=1,ntmax/output_step Do Begin
     If((phi_kftd(p+pmax-1,n+nmax-1,i)+phi_kftd(p+pmax-1,n+nmax-1,i-1)) eq complex(0,0))then begin
     Omegaft(p,n)=Omegaft(p,n)+complex(0,1)*(phi_kftd(p+pmax-1,n+nmax-1,i)-phi_kftd(p+pmax-1,n+nmax-1,i-1))/(output_step*tstep)
     Endif Else Begin
     Omegaft(p,n)=Omegaft(p,n)+2*complex(0,1)*(phi_kftd(p+pmax-1,n+nmax-1,i)-phi_kftd(p+pmax-1,n+nmax-1,i-1))/(output_step*tstep*(phi_kftd(p+pmax-1,n+nmax-1,i)+phi_kftd(p+pmax-1,n+nmax-1,i-1)))
     Endelse
     Endfor
     Omegaft(p,n)=Omegaft(p,n)/(ntmax/output_step)
     
      y1ft=total(alog(abs(Phi_kftd(p+pmax-1,n+nmax-1,lxxx1:rxxx1))))/(rxxx1-lxxx1+1)
      y2ft=total(alog(abs(Phi_kftd(p+pmax-1,n+nmax-1,lxxx2:rxxx2))))/(rxxx2-lxxx2+1)
    
      Omegaft(p,n)=complex(REAL_PART(Omegaft(p,n)),(y2ft-y1ft)/(xx2-xx1))
          
    Endfor
   Endfor
       
    Omegaft(0,0)=complex(0,0)
    name[0]='('+cgGreek('omega')+') CKinFH [Init]' 
    nameap=' '
 Endif
 
;;-----------------------------------------------------------------
;;public code below:
  If(!Omega_p eq 0)then begin
       name[0]='Re'+name[0]+nameap
      !mtitle=name[0]
       xmax=max(x)
       xmin=min(x)
       ymax=max(REAL_PART(Omegaft(*,*)))
       ymin=min(REAL_PART(Omegaft(*,*)))
      For p=0,pmax-1,plotstep Do Begin
      i=p/plotstep
      y(*)=REAL_PART(Omegaft(p,*))
      overplot,x,y,name[1],'',xmax,xmin,ymax,ymin,0,log,i,i,0,zoomz,yrange_style,0
      Endfor

   Endif Else IF(!Omega_p eq 1) THEN BEGIN
      name[0]='Im'+name[0]+nameap
      !mtitle=name[0]
       xmax=max(x)
       xmin=min(x)
       ymax=max(imaginary(Omegaft(*,*)))
       
       ymin=min(imaginary(Omegaft(*,*)))
       ;If(min(imaginary(Omegaft(*,*))) ge -0.5*max(imaginary(Omegaft(*,*))))Then begin
       ;ymin=min(imaginary(Omegaft(*,*)))
       ;Endif Else begin
       ;ymin=-0.5*max(imaginary(Omegaft(*,*)))
       ;Endelse
       ymax=ymax+(ymax-ymin)/5.0       
      For p=0,pmax-1,plotstep Do Begin
      i=p/plotstep
      y(*)=imaginary(Omegaft(p,*))
      overplot,x,y,name[1],'',xmax,xmin,ymax,ymin,0,log,i,i,0,zoomz,yrange_style,0
      Endfor
     gammad=fltarr(pmax,nmax)
     gammad=imaginary(Omegaft(*,*))
      temp=gammad(0,1)  ;;do not set temp=0
      For i=0,pmax-1 Do Begin
        For j=0,nmax-1 Do Begin
          IF(gammad(i,j) NE 0) THEN BEGIN
            IF(temp LE gammad(i,j)) THEN BEGIN
              temp=gammad(i,j)
            Endif
          Endif
        Endfor
      Endfor
      For i=0,pmax-1 Do Begin
        For j=0,nmax-1 Do Begin
          IF(gammad(i,j) NE 0) THEN BEGIN
            IF(gammad(i,j) EQ temp) THEN BEGIN
            gamma_max=gammad(i,j)
            ky_max=j*kymax/(nmax-1)
              print,'kx=',i*kxmax/(pmax-1),'  ','ky=',$
                j*kymax/(nmax-1),'  ','max growth rate mode = ',Omegaft(i,j)
            Endif
          Endif
        Endfor
      Endfor
   !p.color=color_value(ncolor*11/16)  ;;red
   xyouts,kymax/5.0,gamma_max*0.9/zoomz,'Max='+strtrim(gamma_max,2),charthick=3.5,size=4.7
   !p.color=color_value(ncolor+1)
   Endif

      lines = indgen(5)                    ; for line styles
      line = findgen(5)*kxmax/(nmax-1)*plotstep
      items = 'kx='+strtrim(line,2)           ; annotations
      legend,items,linestyle=lines,charsize=1.5;,/bottom      ; at lower left

;      lines = indgen(nmax/plotstep)                    ; for line styles
;      line = findgen(nmax/plotstep)*kxmax/(nmax-1)*plotstep
;      items = 'kx='+strtrim(line,2)           ; annotations
;      legend,items,linestyle=lines,charsize=1.5;,/bottom      ; at lower left

     
      usnum=0
      For i=0,pmax-1 Do Begin
        For j=0,nmax-1 Do Begin
        IF(gammad(i,j) gt 0)Then Begin
        usnum=usnum+1
        Endif
        Endfor
      Endfor
       
   !p.color=color_value(ncolor*1/8) ;green
   ;xyouts,xmin,gamma_max*1.05/zoomz,'Unstable modes: '+strtrim(usnum,2)+'/'+strtrim(round(pmax*nmax),2),charthick=2,size=2
   !p.color=color_value(ncolor+1)    
            
 IF(!threeD eq 1)THEN BEGIN            
  IF(!Omega_p eq 0) THEN BEGIN
   !p.background=color_value(ncolor/2)
   window, plotidh+2, TITLE='main_plot', xsize=500,ysize=500
      x=kxhalf
      y=kyhalf
      z=abs(REAL_PART(Omegaft(*,*)))
      !mtitle=name[0]+nameap
      surface3D,z,x,y,'kx','ky',0,1,0,Az,Ax
  Endif
  IF(!Omega_p eq 1) THEN BEGIN
   !p.background=color_value(ncolor/2)
   window, plotidh+2, TITLE='main_plot', xsize=500,ysize=500
      x=kxhalf
      y=kyhalf
      z=gammad(*,*)
      !mtitle=name[0]+nameap
      surface3D,z,x,y,'kx','ky',0,1,0,Az,Ax
  Endif
  IF(freq eq 1) THEN BEGIN
  If(!Omega_p eq 0)then begin
   !p.background=color_value(ncolor/2)
   window, plotidh+3, TITLE='main_plot', xsize=600,ysize=500
   ymax=max(REAL_PART(Omega_matrftHigh(*,0)))
   ymin=min(REAL_PART(Omega_matrftHigh(*,0)))
  !mtitle='Omegaft'+'_Re High freq(ky=0)'
  plot,kxhalf,REAL_PART(Omega_matrftHigh(*,0)),xtitle='kx',yrange=[ymin,ymax]
  Endif Else IF(!Omega_p eq 1) THEN BEGIN
   !p.background=color_value(ncolor/2)
   window, plotidh+3, TITLE='main_plot', xsize=600,ysize=500
   ymax=max(imaginary(Omega_matrftHigh(*,0)))
   ymin=min(imaginary(Omega_matrftHigh(*,0)))
  !mtitle='Omegaft'+'_Im High freq(ky=0)'
  plot,kxhalf,imaginary(Omega_matrftHigh(*,0)),xtitle='kx',yrange=[ymin,ymax]
  EndIF
  Endif
 ENDIF



EndIF

;=================================================================================================================
IF(strcmp(name[0],'('+cgGreek('omega')+') CKinCH',13)) THEN BEGIN
  tempReconst=3
  x=kyhalf
  y=fltarr(nmax)
  Omegacy=complexarr(pmax,nmax)
;---------------------------------------------------------------------------------
 IF(!MS_Init eq 0)THEN BEGIN 
  Omega_matrcyLow=complexarr(pmax,nmax)
  Omega_matrcyHigh=complexarr(pmax,nmax)
  Omega_matOrgcy=complexarr(2*pmax-1,nmax,N_mu*(2*N_CY-1))
  Omega_matOrgcy(*,*,*)=complex(0,0)
  openr,lun,'Omega_matrcy.txt',/get_lun
  readf,lun,Omega_matOrgcy
  n=eof(lun)
  if n ne 1 then print,'error with file load!!!!!!'
  free_lun,lun

 For p=0,pmax-1 Do Begin
  For n=0,nmax-1 Do Begin
  tempRe=abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,0)))
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
      IF(tempRe LE abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j)))) THEN BEGIN
           tempRe= abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j)))
       Endif
    Endfor
    IF(N_CY GE 2)then begin
    tempRe=tempRe/Omega_star*tempReconst/(N_CY-1) ;!!!!!!!!!!!!!!!!!! 1. 2.
    Endif

   ; tempRe=tempReconst
   ; tempRe=9.8
    print,'tempRe=',tempRe
      
    jtemp=-1
    IF((mu_LK EQ 0) && (mu_HK EQ 0))then begin
      For j=0,N_mu*(2*N_CY-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j))) LE tempRe) ) THEN BEGIN
          ;IF((abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j))) LE tempRe) && (abs(imaginary(Omega_matOrgcy(p+pmax-1,n,j))) ne 0)) THEN BEGIN
          jtemp=j
          growthIm=imaginary(Omega_matOrgcy(p+pmax-1,n,j))
          Endif
      Endfor
      For j=0,N_mu*(2*N_CY-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j))) LE tempRe) ) THEN BEGIN
          ;IF((abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j))) LE tempRe) && (abs(imaginary(Omega_matOrgcy(p+pmax-1,n,j))) ne 0)) THEN BEGIN
           IF(growthIm LT imaginary(Omega_matOrgcy(p+pmax-1,n,j))) THEN BEGIN
           growthIm= imaginary(Omega_matOrgcy(p+pmax-1,n,j))
           jtemp=j
           Endif
          Endif
      Endfor
    Endif else begin
      For j=0,N_mu*(2*N_CY-1)-1 Do Begin
          ;IF((abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j))) LE tempRe) && (abs(imaginary(Omega_matOrgcy(p+pmax-1,n,j))) ne 0)) THEN BEGIN
          IF((abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j))) LE tempRe) ) THEN BEGIN
          jtemp=j
          growthIm=imaginary(Omega_matOrgcy(p+pmax-1,n,j))
          Endif
      Endfor
      For j=0,N_mu*(2*N_CY-1)-1 Do Begin
          ;IF((abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j))) LE tempRe) && (abs(imaginary(Omega_matOrgcy(p+pmax-1,n,j))) ne 0)) THEN BEGIN
          IF((abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j))) LE tempRe) ) THEN BEGIN
           IF(growthIm LT imaginary(Omega_matOrgcy(p+pmax-1,n,j))) THEN BEGIN
           growthIm= imaginary(Omega_matOrgcy(p+pmax-1,n,j))
           jtemp=j
           Endif
          Endif
      Endfor      
    Endelse    
    IF(jtemp EQ -1)then begin
    Omega_matrcyLow(p,n)=complex(0,0)
    Endif else begin
    Omega_matrcyLow(p,n)=Omega_matOrgcy(p+pmax-1,n,jtemp)
    Endelse

    jtemp=-1
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
      If(abs((REAL_PART(Omega_matOrgcy(p+pmax-1,n,j)))-OMGpoint*(1)) LT freqPer*Omega_star) THEN BEGIN
          IF((abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j))) GE tempRe) && (abs(imaginary(Omega_matOrgcy(p+pmax-1,n,j))) ne 0)) THEN BEGIN ;&& (abs(imaginary(Omega_matOrgcy(p+pmax-1,n,j))) ne 0)
          jtemp=j
          growthIm=imaginary(Omega_matOrgcy(p+pmax-1,n,j))
          Endif
      Endif
    Endfor
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
      If(abs((REAL_PART(Omega_matOrgcy(p+pmax-1,n,j)))-OMGpoint*(1)) LT freqPer*Omega_star) THEN BEGIN
          IF((abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j))) GE tempRe) && (abs(imaginary(Omega_matOrgcy(p+pmax-1,n,j))) ne 0)) THEN BEGIN ;&& (abs(imaginary(Omega_matOrgcy(p+pmax-1,n,j))) ne 0)
            IF(imaginary(Omega_matOrgcy(p+pmax-1,n,j)) GT growthIm) THEN BEGIN
            jtemp=j
            growthIm=imaginary(Omega_matOrgcy(p+pmax-1,n,j))
            Endif
          Endif
      Endif    
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrcyHigh(p,n)=complex(0,0)
      IF(p EQ pmax-1 && n Eq 0)then begin
      Omega_matrcyHigh(0,0)=Omega_matOrgcy(pmax-1,0,0)
      ENdif
    Endif else begin
    Omega_matrcyHigh(p,n)=Omega_matOrgcy(p+pmax-1,n,jtemp)
    Endelse

  Endfor
 Endfor

 ;-----------------------------------------------------------
  If(freq eq 0)then begin
  Omegacy(*,*)=Omega_matrcyLow(*,*)
  nameap=' Low freq'
  Endif Else IF(freq eq 1) THEN BEGIN
  Omegacy(*,*)=Omega_matrcyHigh(*,*) 
  nameap=' High freq'
  Endif
 
 ;---------------------------------------------------------------------------------
 Endif Else IF(!MS_Init eq 1)THEN BEGIN 
  Omegacy(*,*)=complex(0,0)
      Phi_kcyd=complexarr(2*pmax-1,2*nmax-1,ntmax/output_step+1)
      openr,lun,'Phi_kcy.txt',/get_lun
      readf,lun,Phi_kcyd
      free_lun,lun  
      
      tmax=max(Time)
      xx1=tmax*0.5
      xx2=tmax-tmax/20.0
      xxx1=round(xx1/(output_step*tstep))
      xxx2=round(xx2/(output_step*tstep))
      lxx1=xx1-tmax/20.0
      rxx1=xx1+tmax/20.0
      lxx2=xx2-tmax/20.0
      rxx2=xx2+tmax/20.0
      lxxx1=round(lxx1/(output_step*tstep))
      rxxx1=round(rxx1/(output_step*tstep))
      lxxx2=round(lxx2/(output_step*tstep))
      rxxx2=round(rxx2/(output_step*tstep))  
      
   For p=0,pmax-1 Do Begin
    For n=0,nmax-1 Do Begin      
     For i=1,ntmax/output_step Do Begin
     If((phi_kcyd(p+pmax-1,n+nmax-1,i)+phi_kcyd(p+pmax-1,n+nmax-1,i-1)) eq complex(0,0))then begin
     Omegacy(p,n)=Omegacy(p,n)+complex(0,1)*(phi_kcyd(p+pmax-1,n+nmax-1,i)-phi_kcyd(p+pmax-1,n+nmax-1,i-1))/(output_step*tstep)
     Endif Else Begin
     Omegacy(p,n)=Omegacy(p,n)+2*complex(0,1)*(phi_kcyd(p+pmax-1,n+nmax-1,i)-phi_kcyd(p+pmax-1,n+nmax-1,i-1))/(output_step*tstep*(phi_kcyd(p+pmax-1,n+nmax-1,i)+phi_kcyd(p+pmax-1,n+nmax-1,i-1)))
     Endelse
     Endfor
     Omegacy(p,n)=Omegacy(p,n)/(ntmax/output_step)  
     
      y1cy=total(alog(abs(Phi_kcyd(p+pmax-1,n+nmax-1,lxxx1:rxxx1))))/(rxxx1-lxxx1+1)
      y2cy=total(alog(abs(Phi_kcyd(p+pmax-1,n+nmax-1,lxxx2:rxxx2))))/(rxxx2-lxxx2+1)
      
      Omegacy(p,n)=complex(REAL_PART(Omegacy(p,n)),(y2cy-y1cy)/(xx2-xx1))           
    Endfor
   Endfor     
    Omegacy(0,0)=complex(0,0)
    name[0]='('+cgGreek('omega')+') CKinCH [Init]' 
    nameap=' '      
 Endif
 
 
  If(!Omega_p eq 0)then begin
       name[0]='Re'+name[0]+nameap
      !mtitle=name[0]
       xmax=max(x)
       xmin=min(x)
       ymax=max(REAL_PART(Omegacy(*,*)))
       ymin=min(REAL_PART(Omegacy(*,*)))
      For p=0,pmax-1,plotstep Do Begin
      i=p/plotstep
      y(*)=REAL_PART(Omegacy(p,*))
      overplot,x,y,name[1],'',xmax,xmin,ymax,ymin,0,log,i,i,0,zoomz,yrange_style,0
      Endfor

   Endif Else IF(!Omega_p eq 1) THEN BEGIN
      name[0]='Im'+name[0]+nameap
      !mtitle=name[0]
       xmax=max(x)
       xmin=min(x)
       ymax=max(imaginary(Omegacy(*,*)))
       
       ymin=min(imaginary(Omegacy(*,*)))
      ; If(min(imaginary(Omegacy(*,*))) ge -0.5*max(imaginary(Omegacy(*,*))))Then begin
      ; ymin=min(imaginary(Omegacy(*,*)))
      ; Endif Else begin
      ; ymin=-0.5*max(imaginary(Omegacy(*,*)))
      ; Endelse
       ymax=ymax+(ymax-ymin)/5.0       
      For p=0,pmax-1,plotstep Do Begin
      i=p/plotstep
      y(*)=imaginary(Omegacy(p,*))
      overplot,x,y,name[1],'',xmax,xmin,ymax,ymin,0,log,i,i,0,zoomz,yrange_style,0
      Endfor
     gammad=fltarr(pmax,nmax)
     gammad=imaginary(Omegacy(*,*))
      temp=gammad(0,1)  ;;do not set temp=0
      For i=0,pmax-1 Do Begin
        For j=0,nmax-1 Do Begin
          IF(gammad(i,j) NE 0) THEN BEGIN
            IF(temp LE gammad(i,j)) THEN BEGIN
              temp=gammad(i,j)
            Endif
          Endif
        Endfor
      Endfor
      gamma_max=0
      For i=0,pmax-1 Do Begin
        For j=0,nmax-1 Do Begin
          IF(gammad(i,j) NE 0) THEN BEGIN
            IF(gammad(i,j) EQ temp) THEN BEGIN
            gamma_max=gammad(i,j)
            ky_max=j*kymax/(nmax-1)
              print,'kx=',i*kxmax/(pmax-1),'  ','ky=',$
                j*kymax/(nmax-1),'  ','max growth rate mode = ',Omegacy(i,j)
            Endif
          Endif
        Endfor
      Endfor
   !p.color=color_value(ncolor*11/16)  ;;red
    xyouts,kymax/5.0,gamma_max*0.9/zoomz,'Max='+strtrim(gamma_max,2),charthick=3.5,size=4.7
   !p.color=color_value(ncolor+1)

     usnum=0
     For i=0,pmax-1 Do Begin
        For j=0,nmax-1 Do Begin
        IF(gammad(i,j) gt 0)Then Begin
        usnum=usnum+1
        Endif
        Endfor
     Endfor
     !p.color=color_value(ncolor*1/8) ;green
     ;xyouts,xmin,gamma_max*1.05/zoomz,'Unstable modes: '+strtrim(usnum,2)+'/'+strtrim(round(pmax*nmax),2),charthick=2,size=2
     !p.color=color_value(ncolor+1)
      
   Endif

      lines = indgen(5)                    ; for line styles
      line = findgen(5)*kxmax/(nmax-1)*plotstep
      items = 'kx='+strtrim(line,2)           ; annotations
      legend,items,linestyle=lines,charsize=1.5;,/bottom      ; at lower lecy

;      lines = indgen(nmax/plotstep)                    ; for line styles
;      line = findgen(nmax/plotstep)*kxmax/(nmax-1)*plotstep
;      items = 'kx='+strtrim(line,2)           ; annotations
;      legend,items,linestyle=lines,charsize=1.5;,/bottom      ; at lower lecy
           

 IF(!threeD eq 1)THEN BEGIN  
  IF(!Omega_p eq 0) THEN BEGIN
   !p.background=color_value(ncolor/2)
   window, plotidh+2, TITLE='main_plot', xsize=500,ysize=500
      x=kxhalf
      y=kyhalf
      z=abs(REAL_PART(Omegacy(*,*)))
      !mtitle=name[0]+nameap
      surface3D,z,x,y,'kx','ky',0,1,0,Az,Ax
  Endif
  IF(!Omega_p eq 1) THEN BEGIN
   !p.background=color_value(ncolor/2)
   window, plotidh+2, TITLE='main_plot', xsize=500,ysize=500
      x=kxhalf
      y=kyhalf
      z=gammad(*,*)
      !mtitle=name[0]+nameap
      surface3D,z,x,y,'kx','ky',0,1,0,Az,Ax
  Endif
 ENDIF

EndIF
;=================================================================================================================
IF(strcmp(name[0],'(CKinFH)|Phi|^2 freq spectrum',25) || strcmp(name[0],'(GK)|Phi|^2 freq spectrum',25)) THEN BEGIN
nametemp=['(GK)|Phi|^2 freq spectrum','Phi_k.txt','ky']

  ii=0
  ptemp=pmax-1+ii

  Nsat=round(ntmax/output_step - sattime/(output_step*tstep))
  Phi_kftsat=complexarr(2*pmax-1,2*nmax-1,Nsat)
  Phi_ktmp=complexarr(2*pmax-1,2*nmax-1)
  Phi_sqsat=fltarr(nmax,Nsat)
  Phi_sqsatO=fltarr(nmax,Nsat)

  openr,lun,name[1],/get_lun
  readf,lun,Phi_ktmp
  POINT_LUN, -lun, pos
  Print,'Pos=',pos

  POINT_LUN, lun, pos*round(sattime/(output_step*tstep))
  readf,lun,Phi_kftsat
  free_lun,lun

  For n=0,nmax-1,plotstep Do Begin
  Phi_sqsat(n,*)=Phi_kftsat(ptemp,n,*)*conj(Phi_kftsat(ptemp,n,*))
  Phi_sqsatO(n,*)=FFT(Phi_sqsat(n,*), DIMENSION=2)
  Endfor
  name[0]=name[0]+'(kx='+strtrim(ii*kxmax/(pmax-1),2)+';'+strtrim(round(sattime),2)+'-'+strtrim(round(ntmax*tstep),2)+'Ln/Cs)'
  !mtitle=name[0]

 If(half eq 0)then begin
  y=fltarr(Nsat)
  x=FINDGEN(Nsat)
  x=(x-Nsat/2+1)/(Nsat*tstep*output_step)*2*!PI
  IF(Omega_max EQ 0)then begin
  xmin=min(x)
  xmax=max(x)
  Endif Else IF(Omega_max GT 0) THEN BEGIN
  xmin=-Omega_max
  xmax=Omega_max
  Endif
  ymin=min(Phi_sqsatO)
  ymax=max(Phi_sqsatO)

  For n=0,nmax-1,plotstep Do Begin
  y(*)=0
  i=n/plotstep
  y(0:Nsat/2-2)=Phi_sqsatO(n,Nsat/2+1:Nsat-1)
  y(Nsat/2-1:Nsat-1)=Phi_sqsatO(n,0:Nsat/2)
  overplot,x,y,'Omega','|Phi|^2',xmax,xmin,ymax,ymin,0,overplot_log,i,i,0,zoomz,yrange_style,0
  Endfor

 Endif Else IF(half eq 1) THEN BEGIN
  y=fltarr(Nsat/2+1)
  x=FINDGEN(Nsat/2+1)
  x=(x)/(Nsat*tstep*output_step)*2*!PI
  IF((Omega_max EQ 0) && (Omega_min EQ 0))then begin
  xmin=min(x)
  xmax=max(x)
  Endif Else IF(Omega_max GT 0) THEN BEGIN
  xmin=Omega_min
  xmax=Omega_max
  Endif
;  ymin=min(Phi_sqsatO)
;  ymax=max(Phi_sqsatO)
  x_mintmp=round(Omega_min*(Nsat*tstep*output_step)/(2*!PI))
  x_maxtmp=round(Omega_max*(Nsat*tstep*output_step)/(2*!PI))
  ymin=min(Phi_sqsatO(*,x_mintmp:x_maxtmp))
  ymax=max(Phi_sqsatO(*,x_mintmp:x_maxtmp))
  print,'Omegamin',x(x_mintmp)
  print,'Omegamax',x(x_maxtmp)

  For n=0,nmax-1,plotstep Do Begin
  y(*)=0
  i=n/plotstep
  y(0:Nsat/2)=Phi_sqsatO(n,0:Nsat/2)
  overplot,x,y,'Omega','|Phi|^2',xmax,xmin,ymax,ymin,0,overplot_log,i,i,0,zoomz,yrange_style,0
  Endfor

 Endif

      lines = indgen(6)                    ; for line styles
      line = findgen(6)*kxmax/(nmax-1)*plotstep
      items = 'ky='+strtrim(line,2)           ; annotations
      legend,items,linestyle=lines,charsize=1.5;,/bottom      ; at lower left

;      lines = indgen(nmax/plotstep)                    ; for line styles
;      line = findgen(nmax/plotstep)*kxmax/(nmax-1)*plotstep
;      items = 'ky='+strtrim(line,2)           ; annotations
;      legend,items,linestyle=lines,charsize=1.5;,/bottom      ; at lower left


EndIF
;=================================================================================================================


end

;*******************************************************************************
pro Main_Plots_event,event

openr,lun,'inputRCYCLO.txt',/get_lun
  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Control_variables',plot_name_temp,19)
  endwhile
  thisFormat='(8x,i2,/,9x,i2,/,4x,i2)'
  readf,lun,restart,GK_FK_CK,CDW,Format=thisFormat
  thisFormat='(8x,i2,/,11x,i2,/,12x,i5,/,11x,i4,/,9x,f5.3,/,/,/,/,7x,i10)'
  readf,lun,muDtype,mugridtype,output_step,backup_num,Const_nl,Stopnt,Format=thisFormat

  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Grid_variables',plot_name_temp,16)
  endwhile
  thisFormat='(6x,i10,/,6x,f12.9,/,6x,f10.3,/,6x,f10.3,/,5x,i5,/,5x,i5,/,6x,i8.3,/,5x,i5,/,5x,i6,/,5x,i6)'
  readf,lun,ntmax,tstep,kxmax,kymax,pmax,nmax,mumax,N_mu,$
      N_FT,N_CY,Format=thisFormat

  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Physics_variables',plot_name_temp,19)
  endwhile
  thisFormat='(11x,f7.4,/,5x,f5.3,/,6x,f6.3,/,9x,f5.3,/,9x,f5.3,/,9x,f5.3,/,7x,f9.3,/,8x,f5.3,/,6x,f5.3,/,6x,f5.3,/,8x,f12.7,/,8x,f8.5,/,/,/,6x,f8.4)'
  readf,lun,Omega_star,a_Ln,a_LTi,lambda_n,lambda_0,lambda_D,AlphaA,delta_1,$
  mu_HK,mu_LK,F_k_int,Epsilon,gamIC, Format=thisFormat

free_lun,lun



openr,lun,'inputvuRCYCLO.txt',/get_lun
  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Control_variables',plot_name_temp,19)
  endwhile
  thisFormat='(7x,i5,/,9x,i5,/,8x,f10.3,/,14x,i2,/,13x,i2,/,/,/,4x,i6,/,4x,i6,/,11x,i4,/,13x,i2,/,5x,i5,/,8x,i2)'
  readf,lun,avenum,plotstep,sattime,Intermit_Plot,overplot_log,ix1,ix2,Switch_sat,quicktrigger,jump,control,Format=thisFormat

  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Range_control',plot_name_temp,15)
  endwhile
  thisFormat='(10x,f15.12,/,10x,f12.4,/,11x,f15.10,/,11x,f15.12,/,9x,i5,/,8x,f10.4,/,10x,f12.3)'
  readf,lun,Omega_min,Omega_max,Phi_logmin,Chi_logmin,OMGpoint,freqPer,begintime,Format=thisFormat

free_lun,lun



  Time=findgen(ntmax/output_step+1)*(tstep*output_step)
  kx=findgen(2*pmax-1)*kxmax/(pmax-1)-kxmax
  ky=findgen(2*nmax-1)*kymax/(nmax-1)-kymax
  kxhalf=findgen(pmax)*kxmax/(pmax-1)
  kyhalf=findgen(nmax)*kymax/(nmax-1)
  mu=fltarr(N_mu)


  common startup,number_plot,fpath,ncolor,color_value,plotid
  common main_plot,plotidh,ntime,ndata,data_name,windowmark,basemark

  widget_control,event.id,get_uvalue=choice
    !p.thick=2

  case choice of

    0: begin ;total Chi (GK)
    nametemp=['total Chi (GK)','Chi','total_Chi.txt']
    twoD_time_trace_see,nametemp,100
    end

    1: begin ;total D (GK)
    nametemp=['total D (GK)','D','total_D_3G.txt']
    twoD_time_trace_see,nametemp,101
    end
    
    2: begin 
    nametemp=['(ne/n0)^2 (GK)','(ne/n0)^2','absn_n0.txt']
    twoD_time_trace_see,nametemp,109
    end    
    
    3: begin 
    nametemp=['electrostatic Energy E (GK)','E$\upGK$','E.txt']
    twoD_time_trace_see,nametemp,112
    end 
        
    4: begin ;Omega (GK)
    nametemp=['Omega(ky)','ky','Phi_k.txt'];,'Omega_the.txt','Omega_matr.txt']
    ;;nametemp include: title, ytitle, and input data file name
     kspectrum1D_see,nametemp,102
    end

    5: begin ;total_Chi (CKinFH)
    nametemp=['total Chi (CKinFH)','Chi','total_Chift.txt']
    twoD_time_trace_see,nametemp,103
    end
    
    6: begin ;total_D (CKinFH)
    nametemp=['total D (CKinFH)','D','total_Dft.txt']
    twoD_time_trace_see,nametemp,104
    end

    7: begin
    nametemp=['(ne/n0)^2 (CKinFH)','(ne/n0)^2','absn_n0ft.txt']
    twoD_time_trace_see,nametemp,110
    end 
    
    8: begin 
    nametemp=['electrostatic Energy E (CKinFH)','E$\upFK$','Eft.txt']
    twoD_time_trace_see,nametemp,113
    end            
                
    9: begin ;Omega (CKinFH)
    nametemp=['Omegaft(ky)','ky','Phi_kcy.txt']
    ;;nametemp include: title, ytitle, and input data file name
     kspectrum1D_see,nametemp,105
    end   
           
    10: begin ;total_Chi (CKinCH)
    nametemp=['total Chi (CKinCH)','Chi','total_Chicy.txt']
    twoD_time_trace_see,nametemp,106
    end
    
    11: begin ;total_D (CKinCH)
    nametemp=['total D (CKinCH)','D','total_Dcy.txt']
    twoD_time_trace_see,nametemp,107
    end

    12: begin 
    nametemp=['(ne/n0)^2 (CKinCH)','(ne/n0)^2','absn_n0cy.txt']
    twoD_time_trace_see,nametemp,111
    end 

    13: begin 
    nametemp=['electrostatic Energy E (CKinCH)','E$\upCK$','Ecy.txt']
    twoD_time_trace_see,nametemp,114
    end 
                            
    14: begin ;Omega (CKinCH)
    nametemp=['Omegacy(ky)','ky','Phi_kcy.txt']
    ;;nametemp include: title, ytitle, and input data file name
     kspectrum1D_see,nametemp,108
    end
    

            

    
  endcase
end
;*******************************************************************************
pro More_details_event,event 

openr,lun,'inputRCYCLO.txt',/get_lun
  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Control_variables',plot_name_temp,19)
  endwhile
  thisFormat='(8x,i2,/,9x,i2,/,4x,i2)'
  readf,lun,restart,GK_FK_CK,CDW,Format=thisFormat
  thisFormat='(8x,i2,/,11x,i2,/,12x,i5,/,11x,i4,/,9x,f5.3,/,/,/,/,7x,i10)'
  readf,lun,muDtype,mugridtype,output_step,backup_num,Const_nl,Stopnt,Format=thisFormat

  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Grid_variables',plot_name_temp,16)
  endwhile
  thisFormat='(6x,i10,/,6x,f12.9,/,6x,f10.3,/,6x,f10.3,/,5x,i5,/,5x,i5,/,6x,i8.3,/,5x,i5,/,5x,i6,/,5x,i6)'
  readf,lun,ntmax,tstep,kxmax,kymax,pmax,nmax,mumax,N_mu,$
      N_FT,N_CY,Format=thisFormat

  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Physics_variables',plot_name_temp,19)
  endwhile
  thisFormat='(11x,f7.4,/,5x,f5.3,/,6x,f6.3,/,9x,f5.3,/,9x,f5.3,/,9x,f5.3,/,7x,f9.3,/,8x,f5.3,/,6x,f5.3,/,6x,f5.3,/,8x,f12.7,/,8x,f8.5,/,/,/,6x,f8.4)'
  readf,lun,Omega_star,a_Ln,a_LTi,lambda_n,lambda_0,lambda_D,AlphaA,delta_1,$
  mu_HK,mu_LK,F_k_int,Epsilon,gamIC, Format=thisFormat

free_lun,lun



openr,lun,'inputvuRCYCLO.txt',/get_lun
  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Control_variables',plot_name_temp,19)
  endwhile
  thisFormat='(7x,i5,/,9x,i5,/,8x,f10.3,/,14x,i2,/,13x,i2,/,/,/,4x,i6,/,4x,i6,/,11x,i4)'
  readf,lun,avenum,plotstep,sattime,Intermit_Plot,overplot_log,ix1,ix2,Switch_sat,Format=thisFormat

  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Range_control',plot_name_temp,15)
  endwhile
  thisFormat='(10x,f15.12,/,10x,f12.4,/,11x,f15.10,/,11x,f15.12,/,9x,i5,/,8x,f10.4,/,10x,f12.4/,3x,i4)'
  readf,lun,Omega_min,Omega_max,Phi_logmin,Chi_logmin,OMGpoint,freqPer,begintime,Nk,Format=thisFormat

free_lun,lun


;  openr,lun,'nt_Odd.txt',/get_lun
;  readf,lun,nt_Odd
;  free_lun,lun
;  openr,lun,'nt_Even.txt',/get_lun
;  readf,lun,nt_Even
;;  free_lun,lun
;  IF(ntmax GE max([nt_Odd,nt_Even]))Then Begin
;  ntmax=max([nt_Odd,nt_Even])
;  Endif
;  print,'ntmax=',ntmax

  Time=findgen(ntmax/output_step+1)*(tstep*output_step)
  kx=findgen(2*pmax-1)*kxmax/(pmax-1)-kxmax
  ky=findgen(2*nmax-1)*kymax/(nmax-1)-kymax
  kxhalf=findgen(pmax)*kxmax/(pmax-1)
  kyhalf=findgen(nmax)*kymax/(nmax-1)
  mu=fltarr(N_mu)
;  openr,lun,'cpu_timeGK.txt',/get_lun
;  plot_name_temp="aaa"
;  result1=0
;  While (result1 eq 0) Do Begin
;  readf,lun,plot_name_temp
;  result1=strcmp(' mu_point=',plot_name_temp,10)
;  endwhile
;  readf,lun,mu
;  free_lun,lun



  common startup,number_plot,fpath,ncolor,color_value,plotid
  common main_plot,plotidh,ntime,ndata,data_name,windowmark,basemark

  widget_control,event.id,get_uvalue=choice
    !p.thick=2

  case choice of

;    0: begin
;      ;plotidh=1
;      ;wdelete,plotidh        ; Closes plot windows
;     ; plotidh=!D.window
;      while(!D.window NE -1) Do begin
;    ;  IF( plotidh GT 0)then Begin
;      Wdelete,!D.window  ;;delete all the windows.
;     ; plotidh=plotidh-1
;      Endwhile
;      widget_control, event.top,/destroy
;    end

    0: begin ;total energy
    ;;nametemp include: title, ytitle, and input data file name

    nametemp=['total energy (conj(F)*J0*Phi+F*J0*conj(Phi))/2 (3DGK)','energy','total_energy.txt']
    twoD_time_trace_see,nametemp,0
    end


    1: begin ;total D
    nametemp=['total D (GK)','D','total_D_3G.txt']
    twoD_time_trace_see,nametemp,1
    end

    2: begin
    nametemp=['total Entropy (GK)','entropy','total_entropy.txt']
    twoD_time_trace_see,nametemp,2
    end

    3: begin ;total Chi
    nametemp=['total Chi (GK)','Chi','total_Chi.txt']
    twoD_time_trace_see,nametemp,3
    end

    4: begin ;Energy(kx), energy sum ky versus kx
    nametemp=['Energy(kx)','energy_kx_ky.txt','kx']
    kspectrum1D_see,nametemp,4
    end

    5: begin ;Energy(ky), energy sum kx versus ky
    nametemp=['Energy(ky)','energy_kx_ky.txt','ky']
    kspectrum1D_see,nametemp,5
    end

    6: begin ;Energy(kx)int&sat: plot Energy(kx) of time=0 and time=steady time(from sattime to the end)
    nametemp=['Energy (kx)int&sat','energy_kx_ky.txt','kx']
    kspectrum1D_see,nametemp,6
    end

    7: begin ;Energy (ky)int&sat:plot Energy(ky) of time=0 and time=steady time
    nametemp=['Energy (ky)int&sat','energy_kx_ky.txt','ky']
    kspectrum1D_see,nametemp,7
    end

    8: begin ;D(kx), D sum ky versus kx !
    nametemp=['D(kx)','D_3G.txt','kx']
    kspectrum1D_see,nametemp,8
    end

    9: begin ;D(ky),D sum kx versus ky
    nametemp=['D(ky)','D_3G.txt','ky']
    kspectrum1D_see,nametemp,9
    end

    10: begin ;Chi(kx), Chi sum ky versus kx !
    nametemp=['Chi(kx)','Chi.txt','kx']
    kspectrum1D_see,nametemp,10
    end

    11: begin ;Chi(ky),Chi sum kx versus ky
    nametemp=['Chi(ky)','Chi.txt','ky']
    kspectrum1D_see,nametemp,11
    end

    12: begin ;Chi(kx)int&sat: plot Chi(kx) of time=0 and time=steady time(from sattime to the end)
    nametemp=['Chi (kx)int&sat','Chi.txt','kx']
    kspectrum1D_see,nametemp,12
    end

    13: begin ;Chi (ky)int&sat:plot Chi(ky) of time=0 and time=steady time
    nametemp=['Chi (ky)int&sat','Chi.txt','ky']
    kspectrum1D_see,nametemp,13
    end

    14: begin
    nametemp=['Energy(k)','energy_k.txt','k.txt','k']
    kspectrum1D_see,nametemp,14
    end

    15: begin
    nametemp=['Energy(DW, ZF, GAM)','Energy','tot_ene_GAM.txt']
    ;;nametemp include: title, ytitle, and input data file name
     twoD_time_trace_see,nametemp,15
    end

    16: begin
    nametemp=['Omega(ky)','ky','Phi_k.txt'];,'Omega_the.txt','Omega_matr.txt']
    ;;nametemp include: title, ytitle, and input data file name
     kspectrum1D_see,nametemp,16
    end

    17: begin
    nametemp=['Omegaft(ky)','ky','Phi_kcy.txt']
    ;;nametemp include: title, ytitle, and input data file name
     kspectrum1D_see,nametemp,17
    end

    18: begin ;plot total Phi^2 vs time
    nametemp=['total|Phi| (GK)','|Phi|','Phi_k.txt']
    ;nametemp=['Phi^2 (3DGK)','|Phi| ^2','total_Phi^2.txt']
    twoD_time_trace_see,nametemp,18
    end

    19: begin
    nametemp=['D (kx)int&sat','D_3G.txt','kx']
    ;;nametemp include: title, ytitle, and input data file name
     kspectrum1D_see,nametemp,19
    end

    20: begin
    nametemp=['D (ky)int&sat','D_3G.txt','ky']
    ;;nametemp include: title, ytitle, and input data file name
     kspectrum1D_see,nametemp,20
    end

    21: begin
    nametemp=['|Phi|^2 (kx)int&sat','Phi^2.txt','kx']
    ;;nametemp include: title, ytitle, and input data file name
     kspectrum1D_see,nametemp,21
    end

    22: begin
    nametemp=['|Phi|^2 (ky)int&sat','Phi^2.txt','ky']
    ;;nametemp include: title, ytitle, and input data file name
     kspectrum1D_see,nametemp,22
    end

    23: begin ;plot time average |Phi| at saturated state --3D plot
         ; default window

      plotidh=number_plot
      !p.background=color_value(ncolor/2)
      window, plotidh, TITLE='main_plot', xsize=600,ysize=600
      Phiabssat=fltarr(2*pmax-1,nmax)

      jump=1
      IF(ntmax/output_step GE 1000)then begin
      jump=10
      Endif Else IF(ntmax/output_step GE 10000)then begin
      jump=100
      Endif Else IF(ntmax/output_step GE 100000)then begin
      jump=1000
      Endif Else IF(ntmax/output_step GE 1000000)then begin
      jump=10000
      Endif

      Phi_knew=complexarr(2*pmax-1,2*nmax-1,ntmax/(output_step*jump)+1)
      Phi_ktmp=complexarr(2*pmax-1,2*nmax-1)

      openr,lun,'Phi_k.txt',/get_lun
      readf,lun,Phi_ktmp
      POINT_LUN, -lun, pos
      Print,'Pos=',pos
      For i=0L,ntmax/output_step,jump Do begin
      POINT_LUN, lun, pos*i
      readf,lun,Phi_ktmp
      Phi_knew(*,*,i/jump)=Phi_ktmp
      Endfor
      free_lun,lun

For j=round(sattime/(output_step*jump*tstep)),round(ntmax/(output_step*jump)) Do Begin
Phiabssat(*,*)=Phiabssat(*,*)+abs(Phi_knew(*,nmax-1:2*nmax-2,j))
Endfor
Phiabssat(*,*)=Phiabssat(*,*)*(1.0/(ntmax/(output_step*jump)-sattime/(output_step*jump*tstep)))
      x=kx
      y=kyhalf
      !mtitle='|'+cgGreek('Phi')+'|'+' (GK) '+strtrim(round(sattime),2)+'-'+strtrim(round(ntmax*tstep),2)+'Ln/Cs'
      surface3D,Phiabssat,x,y,'kx','ky',0,1,overplot_log,Az,Ax,Phi_logmin

    end

    24: begin ;plot time average Chi at saturated state --3D plot

      plotidh=number_plot+1
      !p.background=color_value(ncolor/2)
      window, plotidh, TITLE='main_plot', xsize=600,ysize=600

      Chisat=fltarr(2*pmax-1,nmax)
      jump=1
      IF(ntmax/output_step GE 1000)then begin
      jump=10
      Endif Else IF(ntmax/output_step GE 10000)then begin
      jump=100
      Endif Else IF(ntmax/output_step GE 100000)then begin
      jump=1000
      Endif Else IF(ntmax/output_step GE 1000000)then begin
      jump=10000
      Endif

      Chi_knew=fltarr(2*pmax-1,2*nmax-1,ntmax/(output_step*jump)+1)
      Chi_ktmp=fltarr(2*pmax-1,2*nmax-1)

      openr,lun,'Chi.txt',/get_lun
      readf,lun,Chi_ktmp
      POINT_LUN, -lun, pos
      Print,'Pos=',pos
      For i=0L,ntmax/output_step,jump Do begin
      POINT_LUN, lun, pos*i
      readf,lun,Chi_ktmp
      Chi_knew(*,*,i/jump)=Chi_ktmp
      Endfor
      free_lun,lun

For j=round(sattime/(output_step*jump*tstep)),round(ntmax/(output_step*jump)) Do Begin
Chisat(*,*)=Chisat(*,*)+abs(Chi_knew(*,nmax-1:2*nmax-2,j))
Endfor
Chisat(*,*)=Chisat(*,*)*(1.0/(ntmax/(output_step*jump)-sattime/(output_step*jump*tstep)))
      x=kx
      y=kyhalf
      !mtitle=cgGreek('chi')+'i (GK) '+strtrim(round(sattime),2)+'-'+strtrim(round(ntmax*tstep),2)+'Ln/Cs'
      surface3D,Chisat,x,y,'kx','ky',0,1,overplot_log,Az,Ax,Chi_logmin
    end


    25: begin
  log=0.0
  zoomz = 1.0  ;zoomz is the parameter which used to adjust the range of the axis
  yrange_style=2  ;yrange_style is the parameter which used to choose plot the top, middle or bottom part of the picture.
                  ;for the value of yrange_style/3, 0 for top, 1 for middle, 2 for bottom;

  plotidh=number_plot
   !p.background=color_value(ncolor/2)
  window, plotidh, TITLE='main_plot', xsize=600,ysize=600
      ;name=['0Phi^2.txt','1Phi^2.txt','2Phi^2.txt']   ;for ['2D Fluid','3DGK & G=0.005','3DGK & G=1.0'] respectively
      name=['0D.txt','1D_3G.txt','2D_3G.txt']   ;for ['2D Fluid','3DGK & G=0.005','3DGK & G=1.0'] respectively
      yd0=fltarr(2*pmax-1,2*nmax-1,ntmax/output_step+1) ;'2D Fluid'
      yd1=fltarr(2*pmax-1,2*nmax-1,ntmax/output_step+1) ;
      yd2=fltarr(2*pmax-1,2*nmax-1,ntmax/output_step+1) ;
      openr,lun,name[0],/get_lun
      readf,lun,yd0
      free_lun,lun

      openr,lun,name[1],/get_lun
      readf,lun,yd1
      free_lun,lun

      openr,lun,name[2],/get_lun
      readf,lun,yd2
      free_lun,lun

      IF(strcmp(name[0],'0Phi^2.txt',10) && strcmp(name[1],'1Phi^2.txt',10)) THEN BEGIN
      yd0(*,nmax-1:2*nmax-2,*)=sqrt(yd0(*,nmax-1:2*nmax-2,*))
      yd1(*,nmax-1:2*nmax-2,*)=sqrt(yd1(*,nmax-1:2*nmax-2,*))
      yd2(*,nmax-1:2*nmax-2,*)=sqrt(yd2(*,nmax-1:2*nmax-2,*))
      Endif

      ydave=fltarr(2*pmax-1,nmax,3)

      ydave(*,*,0)=total(yd0(*,nmax-1:2*nmax-2,round(sattime/(output_step*tstep)+1):round(ntmax/output_step)),3);the average 0 of steady state
      ydave(*,*,0)=ydave(*,*,0)/(ntmax/output_step-sattime/(output_step*tstep))

      ydave(*,*,1)=total(yd1(*,nmax-1:2*nmax-2,round(sattime/(output_step*tstep)+1):round(ntmax/output_step)),3);the average 1 of steady state
      ydave(*,*,1)=ydave(*,*,1)/(ntmax/output_step-sattime/(output_step*tstep))

      ydave(*,*,2)=total(yd2(*,nmax-1:2*nmax-2,round(sattime/(output_step*tstep)+1):round(ntmax/output_step)),3);the average 2 of steady state
      ydave(*,*,2)=ydave(*,*,2)/(ntmax/output_step-sattime/(output_step*tstep))

      x=kxhalf
      IF(strcmp(name[0],'0Phi^2.txt',10) && strcmp(name[1],'1Phi^2.txt',10)) THEN BEGIN
      !mtitle='|Phi|rms'
      Endif else begin
      !mtitle='D'
      Endelse
      sumkx_yave=fltarr(nmax,3)
      For i=0,2 Do Begin
      sumkx_yave(*,i)=TOTAL(ydave(*,*,i),1);/(2*pmax-1)
      Endfor
      xmax=max(x)
      xmin=min(x)
      ymax=max(sumkx_yave)
      ;ymin=0.00005
      ymin=min(sumkx_yave)
      For i=0,2 Do Begin
        y=sumkx_yave(*,i)
        overplot,x,y,'ky',' ',xmax,xmin,ymax,ymin,0,log,i,i,0,zoomz,yrange_style,0
      Endfor
      lines = [0,1,2]                 ; for line styles
      items =['2D Fluid','3DGK & G=0.005','3DGK & G=1.0']
      sym = [0]
      legend,items,linestyle=lines,charsize=1.5,charthick=1.5              ; vertical legend---upper left

    end

    26: begin
  log=1.0
  zoomtemp = 1.0  ;zoomz is the parameter which used to adjust the range of the axis
  yrange_style=2  ;yrange_style is the parameter which used to choose plot the top, middle or bottom part of the picture.
                  ;for the value of yrange_style/3, 0 for top, 1 for middle, 2 for bottom;


      ;name=['0Phi^2.txt','1Phi^2.txt','2Phi^2.txt']   ;for [] respectively
      ;namenumsat=['0D.txt','1D.txt']   ;for [1D.txt is base case] respectively
      namenumsat=['0Chi.txt','1Chi.txt']

;  plotidh=number_plot
;   !p.background=color_value(ncolor/2)
;  window, plotidh, TITLE='main_plot', xsize=600,ysize=600
;
;      pmax=13
;      nmax=13
;      kxmax=2.0
;      kymax=2.0
;      yd0=fltarr(2*pmax-1,2*nmax-1,ntmax/output_step+1)
;      openr,lun,namenumsat[0],/get_lun
;      readf,lun,yd0
;      free_lun,lun
;      IF(strcmp(namenumsat[0],'0Phi^2.txt',10) && strcmp(namenumsat[1],'1Phi^2.txt',10)) THEN BEGIN
;      yd0(*,nmax-1:2*nmax-2,*)=sqrt(yd0(*,nmax-1:2*nmax-2,*))
;      Endif
;      ydave0=fltarr(2*pmax-1,nmax)
;      ydave0(*,*)=total(yd0(*,nmax-1:2*nmax-2,round(sattime/(output_step*tstep)+1):round(ntmax/output_step)),3);the average 0 of steady state
;      ydave0(*,*)=ydave0(*,*)/(ntmax/output_step-sattime/(output_step*tstep))
;      sumkx_yave0=fltarr(nmax)
;      sumkx_yave0(*)=TOTAL(ydave0(*,*),1)/(2*pmax-1)
;      x0=findgen(nmax)*kymax/(nmax-1)
;
;      pmax=13
;      nmax=13
;      kxmax=2.0
;      kymax=2.0
;      yd1=fltarr(2*pmax-1,2*nmax-1,ntmax/output_step+1) ;
;      openr,lun,namenumsat[1],/get_lun
;      readf,lun,yd1
;      free_lun,lun
;      IF(strcmp(namenumsat[0],'0Phi^2.txt',10) && strcmp(namenumsat[1],'1Phi^2.txt',10)) THEN BEGIN
;      yd1(*,nmax-1:2*nmax-2,*)=sqrt(yd1(*,nmax-1:2*nmax-2,*))
;      Endif
;      ydave1=fltarr(2*pmax-1,nmax)
;      ydave1(*,*)=total(yd1(*,nmax-1:2*nmax-2,round(sattime/(output_step*tstep)+1):round(ntmax/output_step)),3);the average 0 of steady state
;      ydave1(*,*)=ydave1(*,*)/(ntmax/output_step-sattime/(output_step*tstep))
;      sumkx_yave1=fltarr(nmax)
;      sumkx_yave1(*)=TOTAL(ydave1(*,*),1)/(2*pmax-1)
;      x1=findgen(nmax)*kymax/(nmax-1)
;
;  xmax=max([max(x0),max(x1)])
;  xmin=min([min(x0),min(x1)])
;  ymax=max([max(sumkx_yave0),max(sumkx_yave1)])
;  ymin=min([min(sumkx_yave0),min(sumkx_yave1)])
;  !noeras=0
;  !p.color=color_value(ncolor+1)
;  set_viewport,0.15,0.95,0.15,0.9
;  set_xy,xmin,xmax,ymin,ymax
;  !p.color=color_value(ncolor+1)
;   yrangetemp=Dblarr(2)
;   if yrange_style ge 3  then begin
;   yrange_style=yrange_style-3
;   endif
;   if yrange_style eq 1  then begin  ;;plot the middle part of the data
;   yrangetemp=[(ymin+ymax)/2-(ymax/2-ymin/3)/zoomtemp,(ymin+ymax)/2+(ymax/2-ymin/3)/zoomtemp]
;   endif
;   if yrange_style eq 2 then begin  ;;plot from the min one of the data
;   yrangetemp=[ymin,(ymin+ymax)/2+(ymax/2-ymin/2)/zoomtemp]
;   endif
;   if yrange_style eq 0 then begin  ;;plot upto the max one of the data
;   yrangetemp=[(ymin+ymax)/2-(ymax/2-ymin/2)/zoomtemp,ymax]
;   endif
;     ; !mtitle='|Phi|rms'
;     ;  !mtitle='D'
;      !mtitle='Chi'
;      plot,[0],[0],$
;       /nodata,$
;       xstyle=1,$
;       xminor=0,$
;       xrange=[xmin,xmax],$
;       xtitle='ky',$
;       ystyle=1,$
;       yminor=0,$
;     ;yrange=[0.000000001,1],$
;      yrange=yrangetemp,$
;       ytitle=' ';,$
;      ; /ylog
;
;     oplot,x0,sumkx_yave0,thick=1.0,linestyle=0;,Psym=-5
;     oplot,x1,sumkx_yave1,thick=1.0,linestyle=1;,Psym=-5
;
;      lines = [0,1]                 ; for line styles
;      items =['Gyro-kinetic','Drift kinetic']
;      sym = [0]
;      legend,items,linestyle=lines,charsize=1.5,charthick=1.5              ; vertical legend---upper left
;
;


 ;; plot two lines of gamma vs ky (kx=0)
  plotidh=number_plot
   !p.background=color_value(ncolor/2)
  window, plotidh+1, TITLE='main_plot', xsize=600,ysize=600

  Omega_matOrg=complexarr(2*pmax-1,2*nmax-1,N_mu)
  Omega_matr=complexarr(2*nmax-1)
  openr,lun,'Omega_matr0.txt',/get_lun   ;;0 Gyro-kinetic
  readf,lun,Omega_matOrg
  n=eof(lun)
  if n ne 1 then print,'error with file load!!!!!!'
  free_lun,lun
    For n=0,2*nmax-2 Do Begin
      Omega_matr(n)=Omega_matOrg(pmax-1,n,0)
      For i=1,N_mu-1 Do Begin
      IF(imaginary(Omega_matOrg(pmax-1,n,i)) GT imaginary(Omega_matr(n)))then begin
      Omega_matr(n)=Omega_matOrg(pmax-1,n,i)
      Endif
      Endfor
    Endfor
     gammad0=fltarr(nmax)
     gammad0=imaginary(Omega_matr(nmax-1:2*nmax-2))

  openr,lun,'Omega_matr1.txt',/get_lun   ;;1 Drift-kinetic
  readf,lun,Omega_matOrg
  n=eof(lun)
  if n ne 1 then print,'error with file load!!!!!!'
  free_lun,lun
    For n=0,2*nmax-2 Do Begin
      Omega_matr(n)=Omega_matOrg(pmax-1,n,0)
      For i=1,N_mu-1 Do Begin
      IF(imaginary(Omega_matOrg(pmax-1,n,i)) GT imaginary(Omega_matr(n)))then begin
      Omega_matr(n)=Omega_matOrg(pmax-1,n,i)
      Endif
      Endfor
    Endfor
     gammad1=fltarr(nmax)
     gammad1=imaginary(Omega_matr(nmax-1:2*nmax-2))

     x0=findgen(nmax)*kymax/(nmax-1)
     x1=findgen(nmax)*kymax/(nmax-1)

  xmax=max([max(x0),max(x1)])
  xmin=min([min(x0),min(x1)])
  ymax=max([max(gammad0),max(gammad1)])
  ymax=ymax*1.3
  ymin=min([min(gammad0),min(gammad1)])
  !noeras=0
  !p.color=color_value(ncolor+1)
  set_viewport,0.15,0.95,0.15,0.9
  set_xy,xmin,xmax,ymin,ymax
  !p.color=color_value(ncolor+1)
   yrangetemp=Dblarr(2)
   if yrange_style ge 3  then begin
   yrange_style=yrange_style-3
   endif
   if yrange_style eq 1  then begin  ;;plot the middle part of the data
   yrangetemp=[(ymin+ymax)/2-(ymax/2-ymin/3)/zoomtemp,(ymin+ymax)/2+(ymax/2-ymin/3)/zoomtemp]
   endif
   if yrange_style eq 2 then begin  ;;plot from the min one of the data
   yrangetemp=[ymin,(ymin+ymax)/2+(ymax/2-ymin/2)/zoomtemp]
   endif
   if yrange_style eq 0 then begin  ;;plot upto the max one of the data
   yrangetemp=[(ymin+ymax)/2-(ymax/2-ymin/2)/zoomtemp,ymax]
   endif
     ; !mtitle='|Phi|rms'
     ;  !mtitle='D'
      !mtitle='gamma vs ky (kx=0)'
      plot,[0],[0],$
       /nodata,$
       xstyle=1,$
       xminor=0,$
       xrange=[xmin,xmax],$
       xtitle='ky',$
       ystyle=1,$
       yminor=0,$
     ;yrange=[0.000000001,1],$
      yrange=yrangetemp,$
       ytitle=' ';,$
      ; /ylog

     oplot,x0,gammad0,thick=1.0,linestyle=0;,Psym=-5
     oplot,x1,gammad1,thick=1.0,linestyle=1;,Psym=-5

      lines = [0,1]                 ; for line styles
      items =['Gyro-kinetic','Drift-kinetic']
      sym = [0]
      legend,items,linestyle=lines,charsize=1.5,charthick=1.5              ; vertical legend---upper left

  end



;================================================================================================================
    27: begin   ;;Cyclo vs Fourier (ky)
  ;;We should set N_CY=2!!! Because Omegacyhigh is used as the reference standard for judge the first harmonic of Fourier kinetic high
  ;;frequency motion. Thus we should only keep the first harmonic of Cyclo-Kinetics(n=-1,0,1).
  ;; This part compare the linear dispersion relation of 4D Cyclotron with 4D Fourier method.
  ;; The 3D linear gyro-kinietic dispersion relation is chosen as the standard, and only low
  ;; frequency modes are compared in first and second plots, only high frequency modes are
  ;; compared in the third and fourth plots. Thus we need to take out the low frequency mode
  ;; first in the code, then ...
  tempReconst=0.8
  ptemp=pmax-1   ;choose kx. Start from ptemp=pmax-1 for kx=0

  Omega_matOrg=complexarr(2*pmax-1,nmax,N_mu)
  Omega_matrFinal=complexarr(nmax)
  openr,lun,'Omega_matr.txt',/get_lun   ;;0 Gyro-kinetic
  readf,lun,Omega_matOrg
  n=eof(lun)
  if n ne 1 then print,'error with file load!!!!!!'
  free_lun,lun
   For n=0,nmax-1 Do Begin
     jtemp=0
     growthIm=imaginary(Omega_matOrg(ptemp,n,0))
     For j=0,N_mu-1 Do Begin
       IF(imaginary(Omega_matOrg(ptemp,n,j)) GT growthIm)then begin
       jtemp=j
       growthIm=imaginary(Omega_matOrg(ptemp,n,j))
       Endif
     Endfor
     IF(jtemp EQ -1)then begin
     Omega_matrFinal(n)=complex(0,0)
     Endif else begin
     Omega_matrFinal(n)=Omega_matOrg(ptemp,n,jtemp)
     Endelse

   Endfor

  ;;For cyclotron harmonics:
  ;;Low frequency (<tempRe ) and high growth rate (>tempIm) is the condition of unstable low frequency mode, in which the most unstable one is what we need.
  ;; If the growth rate of all the low frequency modes <tempIm, we set the Omega to be (0,0).
  ;;High frequency (>tempRe) and high growth rate (>tempIm) is the condition of unstable high frequency mode, however we choose the lowest frequency one as
  ;; we only need the lowest ion cyclotron harmonic. To be the same, if the growth rate of all the high frequency modes <tempIm, we set the Omega to be (0,0).
  Omega_matrcyLow=complexarr(nmax)
  Omega_matrcyHigh=complexarr(nmax)
  Omega_matOrgcy=complexarr(2*pmax-1,nmax,N_mu*(2*N_CY-1))
  Omega_matOrgcy(*,*,*)=complex(0,0)

  openr,lun,'Omega_matrcy.txt',/get_lun
  readf,lun,Omega_matOrgcy
  n=eof(lun)
  if n ne 1 then print,'error with file load!!!!!!'
  free_lun,lun

  For n=0,nmax-1 Do Begin
  tempRe=abs(REAL_PART(Omega_matOrgcy(ptemp,n,0)))
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
      IF(tempRe LE abs(REAL_PART(Omega_matOrgcy(ptemp,n,j)))) THEN BEGIN
           tempRe= abs(REAL_PART(Omega_matOrgcy(ptemp,n,j)))
       Endif
    Endfor
    IF(N_CY GE 2)then begin
    tempRe=tempRe/Omega_star*tempReconst ;!!!!!!!!!!!!!!!!!! 1. 2.
    Endif

    jtemp=-1
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgcy(ptemp,n,j))) LE tempRe))THEN BEGIN
          jtemp=j
          growthIm=imaginary(Omega_matOrgcy(ptemp,n,j))
          Endif
    Endfor
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgcy(ptemp,n,j))) LE tempRe))THEN BEGIN
           IF(growthIm LT imaginary(Omega_matOrgcy(ptemp,n,j))) THEN BEGIN
           growthIm= imaginary(Omega_matOrgcy(ptemp,n,j))
           jtemp=j
           Endif
          Endif
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrcyLow(n)=complex(0,0)
    Endif else begin
    Omega_matrcyLow(n)=Omega_matOrgcy(ptemp,n,jtemp)
    Endelse

    jtemp=-1
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgcy(ptemp,n,j))) GE tempRe)) THEN BEGIN
          jtemp=j
          growthIm=imaginary(Omega_matOrgcy(ptemp,n,j))
          Endif
    Endfor
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgcy(ptemp,n,j))) GE tempRe)) THEN BEGIN
            IF(imaginary(Omega_matOrgcy(ptemp,n,j)) GT growthIm) THEN BEGIN
            jtemp=j
            growthIm=imaginary(Omega_matOrgcy(ptemp,n,j))
            Endif
          Endif
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrcyHigh(n)=complex(0,0)
    Endif else begin
    Omega_matrcyHigh(n)=Omega_matOrgcy(ptemp,n,jtemp)
    Endelse

  Endfor

  ;;For Fourier harmonics:
  ;;The same treating of the low frequency modes to Cyclotron harmonics.
  ;;To high frequency modes we need to choose the one with the closest frequency to ion cyclotron modes comes from Cyclotron kinetics. Gives the condition as below:
  ;;High frequency (>tempRe) and high growth rate (>tempIm) is the condition of unstable high frequency mode, however we choose the one who possess the closest frequency
  ;;to ion cyclotron modes comes from Cyclotron kinetics. To be the same, if the growth rate of all the high frequency modes <tempIm, we set the Omega to be (0,0).
  Omega_matrftLow=complexarr(nmax)
  Omega_matrftHigh=complexarr(nmax)
  Omega_matrftHigh2=complexarr(nmax)
  Omega_matrftHigh3=complexarr(nmax)
  Omega_matOrgft=complexarr(2*pmax-1,nmax,N_mu*(2*N_FT-1))
  Omega_matOrgft(*,*,*)=complex(0,0)

  openr,lun,'Omega_matrft.txt',/get_lun   ;;0 Gyro-kinetic
  readf,lun,Omega_matOrgft
  n=eof(lun)
  if n ne 1 then print,'error with file load!!!!!!'
  free_lun,lun

  For n=0,nmax-1 Do Begin
  tempRe=abs(REAL_PART(Omega_matOrgft(ptemp,n,0)))
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
      IF(tempRe LE abs(REAL_PART(Omega_matOrgft(ptemp,n,j)))) THEN BEGIN
           tempRe= abs(REAL_PART(Omega_matOrgft(ptemp,n,j)))
       Endif
    Endfor
    IF(N_FT GE 2)then begin
    tempRe=tempRe/Omega_star*tempReconst;!!!!!!!!!!!! 3. 4. !!used as the limiter of the low freq and high freq
    Endif


    jtemp=-1
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgft(ptemp,n,j))) LE tempRe))THEN BEGIN
          jtemp=j
          growthIm=(imaginary(Omega_matOrgft(ptemp,n,j)))
          EndIf
    Endfor
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgft(ptemp,n,j))) LE tempRe))THEN BEGIN
           IF(growthIm LT (imaginary(Omega_matOrgft(ptemp,n,j)))) THEN BEGIN
           growthIm=imaginary(Omega_matOrgft(ptemp,n,j))
           jtemp=j
           Endif
          Endif
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrftLow(n)=complex(0,0)
    Endif else begin
    Omega_matrftLow(n)=Omega_matOrgft(ptemp,n,jtemp)
    Endelse


    jtemp=-1
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF(abs(REAL_PART(Omega_matOrgft(ptemp,n,j))) GE tempRe) THEN BEGIN
            If(abs(REAL_PART(Omega_matOrgft(ptemp,n,j))-REAL_PART(Omega_matrcyHigh(n))) LT freqPer*abs(REAL_PART(Omega_matrcyHigh(n)))) THEN BEGIN
             growthIm=imaginary(Omega_matOrgft(ptemp,n,j))
             jtemp=j
            Endif
          EndIf
    Endfor
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF(abs(REAL_PART(Omega_matOrgft(ptemp,n,j))) GE tempRe) THEN BEGIN
            If(abs(REAL_PART(Omega_matOrgft(ptemp,n,j))-REAL_PART(Omega_matrcyHigh(n))) LT freqPer*abs(REAL_PART(Omega_matrcyHigh(n)))) THEN BEGIN
             IF((imaginary(Omega_matOrgft(ptemp,n,j))) GT growthIm) THEN BEGIN
             growthIm=imaginary(Omega_matOrgft(ptemp,n,j))
             jtemp=j
             Endif
            Endif
          EndIf
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrftHigh(n)=complex(0,0)
    Endif else begin
    Omega_matrftHigh(n)=Omega_matOrgft(ptemp,n,jtemp)
    Endelse


    jtemp=-1
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF(abs(REAL_PART(Omega_matOrgft(ptemp,n,j))) GE tempRe) THEN BEGIN
            If(abs(REAL_PART(Omega_matOrgft(ptemp,n,j))-2*REAL_PART(Omega_matrcyHigh(n))) LT freqPer*abs(REAL_PART(Omega_matrcyHigh(n)))) THEN BEGIN
             growthIm=imaginary(Omega_matOrgft(ptemp,n,j))
             jtemp=j
            Endif
          EndIf
    Endfor
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF(abs(REAL_PART(Omega_matOrgft(ptemp,n,j))) GE tempRe) THEN BEGIN
            If(abs(REAL_PART(Omega_matOrgft(ptemp,n,j))-2*REAL_PART(Omega_matrcyHigh(n))) LT freqPer*abs(REAL_PART(Omega_matrcyHigh(n)))) THEN BEGIN
             IF((imaginary(Omega_matOrgft(ptemp,n,j))) GT growthIm) THEN BEGIN
             growthIm=imaginary(Omega_matOrgft(ptemp,n,j))
             jtemp=j
             Endif
            Endif
          EndIf
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrftHigh2(n)=complex(0,0)
    Endif else begin
    Omega_matrftHigh2(n)=Omega_matOrgft(ptemp,n,jtemp)
    Endelse

    jtemp=-1
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF(abs(REAL_PART(Omega_matOrgft(ptemp,n,j))) GE tempRe) THEN BEGIN
            If(abs(REAL_PART(Omega_matOrgft(ptemp,n,j))-3*REAL_PART(Omega_matrcyHigh(n))) LT freqPer*abs(REAL_PART(Omega_matrcyHigh(n)))) THEN BEGIN
             growthIm=imaginary(Omega_matOrgft(ptemp,n,j))
             jtemp=j
            Endif
          EndIf
    Endfor
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF(abs(REAL_PART(Omega_matOrgft(ptemp,n,j))) GE tempRe) THEN BEGIN
            If(abs(REAL_PART(Omega_matOrgft(ptemp,n,j))-3*REAL_PART(Omega_matrcyHigh(n))) LT freqPer*abs(REAL_PART(Omega_matrcyHigh(n)))) THEN BEGIN
             IF((imaginary(Omega_matOrgft(ptemp,n,j))) GT growthIm) THEN BEGIN
             growthIm=imaginary(Omega_matOrgft(ptemp,n,j))
             jtemp=j
             Endif
            Endif
          EndIf
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrftHigh3(n)=complex(0,0)
    Endif else begin
    Omega_matrftHigh3(n)=Omega_matOrgft(ptemp,n,jtemp)
    Endelse

  Endfor


  plotidh=number_plot
  !p.background=color_value(ncolor/2)
  window, plotidh+1, TITLE='main_plot', xsize=600,ysize=650   ;;1. Low frequency mode, Real part
  xmax=max(kyhalf)
  xmin=min(kyhalf)
  ymax=max([max(REAL_PART(Omega_matrcyLow)),max(REAL_PART(Omega_matrftLow)),max(REAL_PART(Omega_matrFinal))])
  ymin=min([min(REAL_PART(Omega_matrcyLow)),min(REAL_PART(Omega_matrftLow)),min(REAL_PART(Omega_matrFinal))])
  ymin=ymin-(ymax-ymin)/20
  ymax=ymax+(ymax-ymin)/20
  !noeras=0
  !p.color=color_value(ncolor+1)
  set_viewport,0.15,0.95,0.15,0.9
  set_xy,xmin,xmax,ymin,ymax
  !p.color=color_value(ncolor+1)
  ;!mtitle='Omega(Re) vs ky (kx=0, n=0,+1,-1)'
    !mtitle=strjoin(['Omega(Re) Lowfreq(kx=0,N_CY=',strtrim(string(fix(N_CY)),2),',N_FT=',strtrim(string(fix(N_FT)),2),')'])
      plot,[0],[0],$
       /nodata,$
       Charsize=1.7,$
       xstyle=1,$
       xminor=0,$
       xrange=[xmin,xmax],$
       xtitle='ky',$
       ystyle=1,$
       yminor=0,$
       ;yrange=[-1,0],$
       yrange=[ymin,ymax],$
       ytitle=' ';,$
       ;/ylog


     oplot,kyhalf,REAL_PART(Omega_matrFinal(*)),thick=1.5,linestyle=0,Psym=-6,symsize=2
     !p.color=color_value(ncolor*11/16)
     oplot,kyhalf,REAL_PART(Omega_matrcyLow(*)),thick=1.5,linestyle=0,Psym=-2,symsize=1.5
      !p.color=color_value(ncolor*1/8)
     oplot,kyhalf,REAL_PART(Omega_matrftLow(*)),thick=1.5,linestyle=0,Psym=-5,symsize=1.5
     !p.color=color_value(ncolor+1)
      lines = [0,0,0]                 ; for line styles
      col=fltarr(3)
      col=[color_value(ncolor+1),color_value(ncolor*11/16),color_value(ncolor*1/8) ]
      items =['Gyro-kinetic','Cyclotron','Fourier']
      sym = [0]
      legend,items,linestyle=lines,colors=col,charsize=1.5,charthick=1.5,psym=[6,2,5]              ; vertical legend---upper left

  !p.background=color_value(ncolor/2)
  window, plotidh+2, TITLE='main_plot', xsize=600,ysize=650   ;;2. Low frequency mode, imaginary part
  xmax=max(kyhalf)
  xmin=min(kyhalf)
  ymax=max([max(imaginary(Omega_matrcyLow)),max(imaginary(Omega_matrftLow)),max(imaginary(Omega_matrFinal))])
  ymin=min([min(imaginary(Omega_matrcyLow)),min(imaginary(Omega_matrftLow)),min(imaginary(Omega_matrFinal))])
  ymin=ymin-(ymax-ymin)/20
  ymax=ymax+(ymax-ymin)/20
  !noeras=0
  !p.color=color_value(ncolor+1)
  set_viewport,0.15,0.95,0.15,0.9
  set_xy,xmin,xmax,ymin,ymax
  !p.color=color_value(ncolor+1)
    ;!mtitle='Omega(Im) vs ky (kx=0, n=0,+1,-1)'
    !mtitle=strjoin(['Omega(Im) Lowfreq(kx=0,N_CY=',strtrim(string(fix(N_CY)),2),',N_FT=',strtrim(string(fix(N_FT)),2),')'])
      plot,[0],[0],$
       /nodata,$
       Charsize=1.7,$
       xstyle=1,$
       xminor=0,$
       xrange=[xmin,xmax],$
       xtitle='ky',$
       ystyle=1,$
       yminor=0,$
     ;yrange=[0,0.4],$
      yrange=[ymin,ymax],$
       ytitle=' ';,$
       ;/ylog

     oplot,kyhalf,imaginary(Omega_matrFinal(*)),thick=1.5,linestyle=0,Psym=-6,symsize=2
     !p.color=color_value(ncolor*11/16)
     oplot,kyhalf,imaginary(Omega_matrcyLow(*)),thick=1.5,linestyle=0,Psym=-2,symsize=1.5
     !p.color=color_value(ncolor*1/8)
     oplot,kyhalf,imaginary(Omega_matrftLow(*)),thick=1.5,linestyle=0,Psym=-5,symsize=1.5
     !p.color=color_value(ncolor+1)
      lines = [0,0,0]                 ; for line styles
      col=fltarr(3)
      col=[color_value(ncolor+1),color_value(ncolor*11/16),color_value(ncolor*1/8) ]
      items =['Gyro-kinetic','Cyclotron','Fourier']
      sym = [0]
      legend,items,linestyle=lines,colors=col,charsize=1.5,charthick=1.5,psym=[6,2,5]      ; vertical legend---upper left

   i=pmax-1
   jf=nmax
print,'Omega(kx='+strtrim(String(i*kxmax/(pmax-1)-kxmax,FORMAT='(f0.2)'),2)+'ky='+strtrim(String(jf*kymax/(nmax-1)-kymax,FORMAT='(f0.2)'),2)+')=',imaginary(Omega_matrFinal(jf-nmax+1))
print,'OmegacyLow(kx='+strtrim(String(i*kxmax/(pmax-1)-kxmax,FORMAT='(f0.2)'),2)+'ky='+strtrim(String(jf*kymax/(nmax-1)-kymax,FORMAT='(f0.2)'),2)+')=',imaginary(Omega_matrcyLow(jf-nmax+1))
print,'OmegaftLow(kx='+strtrim(String(i*kxmax/(pmax-1)-kxmax,FORMAT='(f0.2)'),2)+'ky='+strtrim(String(jf*kymax/(nmax-1)-kymax,FORMAT='(f0.2)'),2)+')=',imaginary(Omega_matrftLow(jf-nmax+1))

  !p.background=color_value(ncolor/2)
  window, plotidh+3, TITLE='main_plot', xsize=600,ysize=650   ;;3. High frequency mode, Real part
  xmax=max(kyhalf)
  xmin=min(kyhalf)
  ymax=max([max(REAL_PART(Omega_matrcyHigh)),max(REAL_PART(Omega_matrftHigh))])
  ymin=min([min(REAL_PART(Omega_matrcyHigh)),min(REAL_PART(Omega_matrftHigh))])
  ymin=ymin-(ymax-ymin)/20
  ymax=ymax+(ymax-ymin)/20
  !noeras=0
  !p.color=color_value(ncolor+1)
  set_viewport,0.15,0.95,0.15,0.9
  set_xy,xmin,xmax,ymin,ymax
  !p.color=color_value(ncolor+1)
  ;!mtitle='Omega(Re) vs ky (kx=0, n=0,+1,-1)'
    !mtitle=strjoin(['Omega(Re) Highfreq(kx=0,N_CY=',strtrim(string(fix(N_CY)),2),',N_FT=',strtrim(string(fix(N_FT)),2),')'])
      plot,[0],[0],$
       /nodata,$
       xstyle=1,$
       Charsize=1.7,$
       xminor=0,$
       xrange=[xmin,xmax],$
       xtitle='ky',$
       ystyle=1,$
       yminor=0,$
     ;yrange=[-800,800],$
     yrange=[ymin,ymax],$
       ytitle=' ';,$
       ;/ylog

     !p.color=color_value(ncolor*11/16)
     oplot,kyhalf,REAL_PART(Omega_matrcyHigh(*)),thick=1.5,linestyle=0,Psym=-2,symsize=1.5
      !p.color=color_value(ncolor*1/8)
     oplot,kyhalf,REAL_PART(Omega_matrftHigh(*)),thick=1.5,linestyle=0,Psym=-5,symsize=1.5
     !p.color=color_value(ncolor+1)
      lines = [0,0]                 ; for line styles
      col=fltarr(2)
      col=[color_value(ncolor*11/16),color_value(ncolor*1/8) ]
      items =['Cyclotron','Fourier']
      sym = [0]
      legend,items,linestyle=lines,colors=col,charsize=1.5,charthick=1.5,psym=[2,5]              ; vertical legend---upper left

  !p.background=color_value(ncolor/2)
  window, plotidh+4, TITLE='main_plot', xsize=600,ysize=650   ;;4. High frequency mode, imaginary part
  xmax=max(kyhalf)
  xmin=min(kyhalf)
  ymax=max([max(imaginary(Omega_matrcyHigh)),max(imaginary(Omega_matrftHigh))])
  ymin=min([min(imaginary(Omega_matrcyHigh)),min(imaginary(Omega_matrftHigh))])
  ymin=ymin-(ymax-ymin)/20
  ymax=ymax+(ymax-ymin)/20
  !noeras=0
  !p.color=color_value(ncolor+1)
  set_viewport,0.15,0.95,0.15,0.9
  set_xy,xmin,xmax,ymin,ymax
  !p.color=color_value(ncolor+1)
    ;!mtitle='Omega(Im) vs ky (kx=0, n=0,+1,-1)'
    !mtitle=strjoin(['Omega(Im) Highfreq(kx=0,N_CY=',strtrim(string(fix(N_CY)),2),',N_FT=',strtrim(string(fix(N_FT)),2),')'])
      plot,[0],[0],$
       /nodata,$
       Charsize=1.6,$
       xstyle=1,$
       xminor=0,$
       xrange=[xmin,xmax],$
       xtitle='ky',$
       ystyle=1,$
       yminor=0,$
      ;yrange=[0,0.00015],$
      yrange=[ymin,ymax],$
       ytitle=' ';,$
       ;/ylog

     !p.color=color_value(ncolor*11/16)
     oplot,kyhalf,imaginary(Omega_matrcyHigh(*)),thick=1.5,linestyle=0,Psym=-2,symsize=1.5
     !p.color=color_value(ncolor*1/8)
     oplot,kyhalf,imaginary(Omega_matrftHigh(*)),thick=1.5,linestyle=0,Psym=-5,symsize=1.5
     !p.color=color_value(ncolor+1)
      lines = [0,0]                 ; for line styles
      col=fltarr(2)
      col=[color_value(ncolor*11/16),color_value(ncolor*1/8) ]
      items =['Cyclotron','Fourier']
      sym = [0]
      legend,items,linestyle=lines,colors=col,charsize=1.5,charthick=1.5,psym=[2,5]      ; vertical legend---upper left

   i=pmax-1
print,'OmegacyHigh(kx='+strtrim(String(i*kxmax/(pmax-1)-kxmax,FORMAT='(f0.2)'),2)+'ky='+strtrim(String(jf*kymax/(nmax-1)-kymax,FORMAT='(f0.2)'),2)+')=',imaginary(Omega_matrcyHigh(jf-nmax+1))
print,'OmegaftHigh(kx='+strtrim(String(i*kxmax/(pmax-1)-kxmax,FORMAT='(f0.2)'),2)+'ky='+strtrim(String(jf*kymax/(nmax-1)-kymax,FORMAT='(f0.2)'),2)+')=',imaginary(Omega_matrftHigh(jf-nmax+1))

   ;;
  !p.background=color_value(ncolor/2)
  window, plotidh+8, TITLE='main_plot', xsize=600,ysize=650   ;;3. High frequency mode, Real part
  xmax=max(kyhalf)
  xmin=min(kyhalf)
  ymax=max([max(REAL_PART(Omega_matrcyHigh)),max(REAL_PART(Omega_matrftHigh)),max(REAL_PART(Omega_matrftHigh2)),max(REAL_PART(Omega_matrftHigh3))])
  ymin=min([min(REAL_PART(Omega_matrcyHigh)),min(REAL_PART(Omega_matrftHigh)),min(REAL_PART(Omega_matrftHigh2)),min(REAL_PART(Omega_matrftHigh3))])
  ymin=ymin-(ymax-ymin)/20
  ymax=ymax+(ymax-ymin)/20
  !noeras=0
  !p.color=color_value(ncolor+1)
  set_viewport,0.15,0.95,0.15,0.9
  set_xy,xmin,xmax,ymin,ymax
  !p.color=color_value(ncolor+1)
  ;!mtitle='Omega(Re) vs ky (kx=0, n=0,+1,-1)'
    !mtitle=strjoin(['Omega(Re) 3*Hf(kx=0,N_CY=',strtrim(string(fix(N_CY)),2),',N_FT=',strtrim(string(fix(N_FT)),2),')'])
      plot,[0],[0],$
       /nodata,$
       xstyle=1,$
       Charsize=1.7,$
       xminor=0,$
       xrange=[xmin,xmax],$
       xtitle='ky',$
       ystyle=1,$
       yminor=0,$
     ;yrange=[-800,800],$
     yrange=[ymin,ymax],$
       ytitle=' ';,$
       ;/ylog
     !p.color=color_value(ncolor*11/16)
     oplot,kyhalf,REAL_PART(Omega_matrcyHigh(*)),thick=1.5,linestyle=0,Psym=-2,symsize=2
      !p.color=color_value(ncolor*1/8)
     oplot,kyhalf,REAL_PART(Omega_matrftHigh(*)),thick=1.5,linestyle=0,Psym=-5,symsize=1.5
     oplot,kyhalf,REAL_PART(Omega_matrftHigh2(*)),thick=1.5,linestyle=1,Psym=-5,symsize=1.5
     oplot,kyhalf,REAL_PART(Omega_matrftHigh3(*)),thick=1.5,linestyle=2,Psym=-5,symsize=1.5
     !p.color=color_value(ncolor+1)
      lines = [0,0,1,2]                 ; for line styles
      col=fltarr(4)
      col=[color_value(ncolor*11/16),color_value(ncolor*1/8),color_value(ncolor*1/8),color_value(ncolor*1/8) ]
      items =['CY 1st','FT 1st','FT 2nd','FT 3rd']
      sym = [0]
      legend,items,linestyle=lines,colors=col,charsize=1.5,charthick=1.5,psym=[-2,-5,-5,-5]                ; vertical legend---upper left

  !p.background=color_value(ncolor/2)
  window, plotidh+9, TITLE='main_plot', xsize=600,ysize=650   ;;4. High frequency mode, imaginary part
  xmax=max(kyhalf)
  xmin=min(kyhalf)
  ymax=max([max(imaginary(Omega_matrcyHigh)),max(imaginary(Omega_matrftHigh)),max(imaginary(Omega_matrftHigh2)),max(imaginary(Omega_matrftHigh3))])
  ymin=min([min(imaginary(Omega_matrcyHigh)),min(imaginary(Omega_matrftHigh)),min(imaginary(Omega_matrftHigh2)),min(imaginary(Omega_matrftHigh3))])
  ymin=ymin-(ymax-ymin)/20
  ymax=ymax+(ymax-ymin)/20
  !noeras=0
  !p.color=color_value(ncolor+1)
  set_viewport,0.15,0.95,0.15,0.9
  set_xy,xmin,xmax,ymin,ymax
  !p.color=color_value(ncolor+1)
    ;!mtitle='Omega(Im) vs ky (kx=0, n=0,+1,-1)'
    !mtitle=strjoin(['Omega(Im) 3*Hf(kx=0,N_CY=',strtrim(string(fix(N_CY)),2),',N_FT=',strtrim(string(fix(N_FT)),2),')'])
      plot,[0],[0],$
       /nodata,$
       Charsize=1.6,$
       xstyle=1,$
       xminor=0,$
       xrange=[xmin,xmax],$
       xtitle='ky',$
       ystyle=1,$
       yminor=0,$
      ;yrange=[0,0.00015],$
      yrange=[ymin,ymax],$
       ytitle=' ';,$
       ;/ylog
     !p.color=color_value(ncolor*11/16)
     oplot,kyhalf,imaginary(Omega_matrcyHigh(*)),thick=1.5,linestyle=0,Psym=-2,symsize=1.5
      !p.color=color_value(ncolor*1/8)
     oplot,kyhalf,imaginary(Omega_matrftHigh(*)),thick=1.5,linestyle=0,Psym=-5,symsize=1.5
     oplot,kyhalf,imaginary(Omega_matrftHigh2(*)),thick=1.5,linestyle=1,Psym=-5,symsize=1.5
     oplot,kyhalf,imaginary(Omega_matrftHigh3(*)),thick=1.5,linestyle=2,Psym=-5,symsize=1.5
     !p.color=color_value(ncolor+1)
      lines = [0,0,1,2]                 ; for line styles
      col=fltarr(4)
      col=[color_value(ncolor*11/16),color_value(ncolor*1/8),color_value(ncolor*1/8),color_value(ncolor*1/8) ]
      items =['CY 1st','FT 1st','FT 2nd','FT 3rd']
      sym = [0]
      legend,items,linestyle=lines,colors=col,charsize=1.5,charthick=1.5,psym=[-2,-5,-5,-5]                ; vertical legend---upper left
    end


;;====================================================================================================================================

    28: begin
   ;;Cyclo vs Fourier (Mult FT), here Mult means run three different N_FT case to put in the same plot.
   ;;We should set N_CY=2!!! Because Omegacyhigh is used as the reference standard for judge the first harmonic of Fourier kinetic high
   ;;frequency motion. Thus we should only keep the first harmonic of Cyclo-Kinetics(n=-1,0,1).
  N_FTmult=intarr(3)
  N_FTmult(*)=[2,6,8]
  tempReconst=0.8
  ptemp=pmax-1   ;choose kx. Start from ptemp=pmax-1 for kx=0

  Omega_matOrg=complexarr(2*pmax-1,nmax,N_mu)
  Omega_matrFinal=complexarr(nmax)
  openr,lun,'Omega_matr.txt',/get_lun   ;;0 Gyro-kinetic
  readf,lun,Omega_matOrg
  n=eof(lun)
  if n ne 1 then print,'error with file load!!!!!!'
  free_lun,lun
   For n=0,nmax-1 Do Begin
     jtemp=0
     growthIm=imaginary(Omega_matOrg(ptemp,n,0))
     For j=0,N_mu-1 Do Begin
       IF(imaginary(Omega_matOrg(ptemp,n,j)) GT growthIm)then begin
       jtemp=j
       growthIm=imaginary(Omega_matOrg(ptemp,n,j))
       Endif
     Endfor
     IF(jtemp EQ -1)then begin
     Omega_matrFinal(n)=complex(0,0)
     Endif else begin
     Omega_matrFinal(n)=Omega_matOrg(ptemp,n,jtemp)
     Endelse

   Endfor

  ;;For cyclotron harmonics:
  ;;Low frequency (<tempRe ) and high growth rate (>tempIm) is the condition of unstable low frequency mode, in which the most unstable one is what we need.
  ;; If the growth rate of all the low frequency modes <tempIm, we set the Omega to be (0,0).
  ;;High frequency (>tempRe) and high growth rate (>tempIm) is the condition of unstable high frequency mode, however we choose the lowest frequency one as
  ;; we only need the lowest ion cyclotron harmonic. To be the same, if the growth rate of all the high frequency modes <tempIm, we set the Omega to be (0,0).
  Omega_matrcyLow=complexarr(nmax)
  Omega_matrcyHigh=complexarr(nmax)
  Omega_matOrgcy=complexarr(2*pmax-1,nmax,N_mu*(2*N_CY-1))
  Omega_matOrgcy(*,*,*)=complex(0,0)


  openr,lun,'Omega_matrcy.txt',/get_lun
  readf,lun,Omega_matOrgcy
  n=eof(lun)
  if n ne 1 then print,'error with file load!!!!!!'
  free_lun,lun

  For n=0,nmax-1 Do Begin
  tempRe=abs(REAL_PART(Omega_matOrgcy(ptemp,n,0)))
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
      IF(tempRe LE abs(REAL_PART(Omega_matOrgcy(ptemp,n,j)))) THEN BEGIN
           tempRe= abs(REAL_PART(Omega_matOrgcy(ptemp,n,j)))
       Endif
    Endfor
    IF(N_CY GE 2)then begin
    tempRe=tempRe/Omega_star*tempReconst/(N_CY-1) ;!!!!!!!!!!!!!!!!!! 1. 2.
    Endif

    jtemp=-1
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgcy(ptemp,n,j))) LE tempRe))THEN BEGIN
          jtemp=j
          growthIm=imaginary(Omega_matOrgcy(ptemp,n,j))
          Endif
    Endfor
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgcy(ptemp,n,j))) LE tempRe))THEN BEGIN
           IF(growthIm LT imaginary(Omega_matOrgcy(ptemp,n,j))) THEN BEGIN
           growthIm= imaginary(Omega_matOrgcy(ptemp,n,j))
           jtemp=j
           Endif
          Endif
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrcyLow(n)=complex(0,0)
    Endif else begin
    Omega_matrcyLow(n)=Omega_matOrgcy(ptemp,n,jtemp)
    Endelse

    jtemp=-1
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgcy(ptemp,n,j))) GE tempRe)) THEN BEGIN
          jtemp=j
          growthIm=imaginary(Omega_matOrgcy(ptemp,n,j))
          Endif
    Endfor
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgcy(ptemp,n,j))) GE tempRe)) THEN BEGIN
            IF(imaginary(Omega_matOrgcy(ptemp,n,j)) GT growthIm) THEN BEGIN
            jtemp=j
            growthIm=imaginary(Omega_matOrgcy(ptemp,n,j))
            Endif
          Endif
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrcyHigh(n)=complex(0,0)
    Endif else begin
    Omega_matrcyHigh(n)=Omega_matOrgcy(ptemp,n,jtemp)
    Endelse

  Endfor

  ;;For Fourier harmonics:
  ;;The same treating of the low frequency modes to Cyclotron harmonics.
  ;;To high frequency modes we need to choose the one with the closest frequency to ion cyclotron modes comes from Cyclotron kinetics. Gives the condition as below:
  ;;High frequency (>tempRe) and high growth rate (>tempIm) is the condition of unstable high frequency mode, however we choose the one who possess the closest frequency
  ;;to ion cyclotron modes comes from Cyclotron kinetics. To be the same, if the growth rate of all the high frequency modes <tempIm, we set the Omega to be (0,0).
  Omega_matrftLowmult=complexarr(nmax,3)
  Omega_matrftHighmult=complexarr(nmax,3)
  Omega_matOrgftmult=complexarr(2*pmax-1,nmax,N_mu*(2*N_FTmult(2)-1),3)
  Omega_matOrgftmult_0=complexarr(2*pmax-1,nmax,N_mu*(2*N_FTmult(0)-1))
  Omega_matOrgftmult_1=complexarr(2*pmax-1,nmax,N_mu*(2*N_FTmult(1)-1))
  Omega_matOrgftmult_2=complexarr(2*pmax-1,nmax,N_mu*(2*N_FTmult(2)-1))

  Omegaftname=string(fix(N_FTmult(0)))
  Omegaftname=strjoin(['Omega_matrftN_FT', strtrim(Omegaftname,2), '.txt'])
  openr,lun,Omegaftname,/get_lun   ;;0 Gyro-kinetic
  readf,lun,Omega_matOrgftmult_0
  n=eof(lun)
  if n ne 1 then print,'error with file load!!!!!!'
  free_lun,lun

  Omegaftname=string(fix(N_FTmult(1)))
  Omegaftname=strjoin(['Omega_matrftN_FT', strtrim(Omegaftname,2), '.txt'])
  openr,lun,Omegaftname,/get_lun   ;;0 Gyro-kinetic
  readf,lun,Omega_matOrgftmult_1
  n=eof(lun)
  if n ne 1 then print,'error with file load!!!!!!'
  free_lun,lun

  Omegaftname=string(fix(N_FTmult(2)))
  Omegaftname=strjoin(['Omega_matrftN_FT', strtrim(Omegaftname,2), '.txt'])
  openr,lun,Omegaftname,/get_lun   ;;0 Gyro-kinetic
  readf,lun,Omega_matOrgftmult_2
  n=eof(lun)
  if n ne 1 then print,'error with file load!!!!!!'
  free_lun,lun

  Omega_matOrgftmult(*,*,*,*)=complex(0,0)
  Omega_matOrgftmult(*,*,0:N_mu*(2*N_FTmult(0)-1)-1,0)=Omega_matOrgftmult_0
  Omega_matOrgftmult(*,*,0:N_mu*(2*N_FTmult(1)-1)-1,1)=Omega_matOrgftmult_1
  Omega_matOrgftmult(*,*,0:N_mu*(2*N_FTmult(2)-1)-1,2)=Omega_matOrgftmult_2

For count=0,2 Do Begin
  For n=0,nmax-1 Do Begin
  tempRe=abs(REAL_PART(Omega_matOrgftmult(ptemp,n,0,count)))
    For j=0,N_mu*(2*N_FTmult(count)-1)-1 Do Begin
      IF(tempRe LE abs(REAL_PART(Omega_matOrgftmult(ptemp,n,j,count)))) THEN BEGIN
           tempRe= abs(REAL_PART(Omega_matOrgftmult(ptemp,n,j,count)))
       Endif
    Endfor
    IF(N_FTmult(count) GE 2)then begin
    tempRe=tempRe/Omega_star*tempReconst;!!!!!!!!!!!! 3. 4.
    Endif

    jtemp=-1
    For j=0,N_mu*(2*N_FTmult(count)-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgftmult(ptemp,n,j,count))) LE tempRe))THEN BEGIN
          jtemp=j
          growthIm=(imaginary(Omega_matOrgftmult(ptemp,n,j,count)))
          EndIf
    Endfor
    For j=0,N_mu*(2*N_FTmult(count)-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgftmult(ptemp,n,j,count))) LE tempRe))THEN BEGIN
           IF(growthIm LT (imaginary(Omega_matOrgftmult(ptemp,n,j,count)))) THEN BEGIN
           growthIm=imaginary(Omega_matOrgftmult(ptemp,n,j,count))
           jtemp=j
           Endif
          Endif
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrftLowmult(n,count)=complex(0,0)
    Endif else begin
    Omega_matrftLowmult(n,count)=Omega_matOrgftmult(ptemp,n,jtemp,count)
    Endelse

    jtemp=-1
    For j=0,N_mu*(2*N_FTmult(count)-1)-1 Do Begin
          IF(abs(REAL_PART(Omega_matOrgftmult(ptemp,n,j,count))) GE tempRe) THEN BEGIN
            If(abs(REAL_PART(Omega_matOrgftmult(ptemp,n,j,count))-REAL_PART(Omega_matrcyHigh(n))) LT freqPer*abs(REAL_PART(Omega_matrcyHigh(n)))) THEN BEGIN
             growthIm=imaginary(Omega_matOrgftmult(ptemp,n,j,count))
             jtemp=j
            Endif
          EndIf
    Endfor
    For j=0,N_mu*(2*N_FTmult(count)-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgftmult(ptemp,n,j,count))) GE tempRe)) THEN BEGIN
            If(abs(REAL_PART(Omega_matOrgftmult(ptemp,n,j,count))-REAL_PART(Omega_matrcyHigh(n))) LT freqPer*abs(REAL_PART(Omega_matrcyHigh(n)))) THEN BEGIN
             IF(imaginary(Omega_matOrgftmult(ptemp,n,j,count)) GT growthIm) THEN BEGIN
             growthIm=imaginary(Omega_matOrgftmult(ptemp,n,j,count))
             jtemp=j
             Endif
            Endif
          EndIf
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrftHighmult(n,count)=complex(0,0)
    Endif else begin
    Omega_matrftHighmult(n,count)=Omega_matOrgftmult(ptemp,n,jtemp,count)
    Endelse

  Endfor
Endfor


  plotidh=number_plot
  !p.background=color_value(ncolor/2)
  window, plotidh+1, TITLE='main_plot', xsize=600,ysize=650   ;;1. Low frequency mode, Real part
  xmax=max(kyhalf)
  xmin=min(kyhalf)
  ymax=max([max(REAL_PART(Omega_matrcyLow)),max(REAL_PART(Omega_matrftLowmult)),max(REAL_PART(Omega_matrFinal))])
  ymin=min([min(REAL_PART(Omega_matrcyLow)),min(REAL_PART(Omega_matrftLowmult)),min(REAL_PART(Omega_matrFinal))])
  ymin=ymin-(ymax-ymin)/20
  ymax=ymax+(ymax-ymin)/20
  !noeras=0
  !p.color=color_value(ncolor+1)
  set_viewport,0.15,0.95,0.15,0.9
  set_xy,xmin,xmax,ymin,ymax
  !p.color=color_value(ncolor+1)
  ;!mtitle='Omega(Re) vs ky (kx=0, n=0,+1,-1)'
    !mtitle=strjoin(['Omega(Re) Low freq(kx=0,N_CY=',strtrim(string(fix(N_CY)),2),')'])
      plot,[0],[0],$
       /nodata,$
       Charsize=1.7,$
       xstyle=1,$
       xminor=0,$
       xrange=[xmin,xmax],$
       xtitle='ky',$
       ystyle=1,$
       yminor=0,$
       ;yrange=[-1,0],$
       yrange=[ymin,ymax],$
       ytitle=' ';,$
       ;/ylog


     oplot,kyhalf,REAL_PART(Omega_matrFinal(*)),thick=1.5,linestyle=0,Psym=-6,symsize=2.1
     !p.color=color_value(ncolor*11/16)
     oplot,kyhalf,REAL_PART(Omega_matrcyLow(*)),thick=1.5,linestyle=0,Psym=-2,symsize=2
      !p.color=color_value(ncolor*1/8)
     oplot,kyhalf,REAL_PART(Omega_matrftLowmult(*,0)),thick=1.5,linestyle=1,Psym=-5,symsize=1.5
     oplot,kyhalf,REAL_PART(Omega_matrftLowmult(*,1)),thick=1.5,linestyle=2,Psym=-5,symsize=1.5
     oplot,kyhalf,REAL_PART(Omega_matrftLowmult(*,2)),thick=1.5,linestyle=0,Psym=-5,symsize=1.5
     !p.color=color_value(ncolor+1)
      lines = [0,0,1,2,0]                 ; for line styles
      col=fltarr(5)
      col=[color_value(ncolor+1),color_value(ncolor*11/16),color_value(ncolor*1/8),color_value(ncolor*1/8),color_value(ncolor*1/8) ]
      items =['Gyro-kinetic',strjoin(['Cyclo N_CY=',strtrim(string(fix(N_CY)),2)]),strjoin(['Fourier N_FT=',strtrim(string(fix(N_FTmult(0))),2)]),$
      strjoin(['Fourier N_FT=',strtrim(string(fix(N_FTmult(1))),2)]),strjoin(['Fourier N_FT=',strtrim(string(fix(N_FTmult(2))),2)])]
      sym = [0]
      legend,items,linestyle=lines,colors=col,charsize=1.5,charthick=1.5,psym=[-6,-2,-5,-5,-5]              ; vertical legend---upper left

  !p.background=color_value(ncolor/2)
  window, plotidh+2, TITLE='main_plot', xsize=600,ysize=650   ;;2. Low frequency mode, imaginary part
  xmax=max(kyhalf)
  xmin=min(kyhalf)
  ymax=max([max(imaginary(Omega_matrcyLow)),max(imaginary(Omega_matrftLowmult)),max(imaginary(Omega_matrFinal))])
  ymin=min([min(imaginary(Omega_matrcyLow)),min(imaginary(Omega_matrftLowmult)),min(imaginary(Omega_matrFinal))])
  ;ymin=ymin-(ymax-ymin)/5
  ;ymax=ymax+(ymax-ymin)/1
  !noeras=0
  !p.color=color_value(ncolor+1)
  set_viewport,0.15,0.95,0.15,0.9
  set_xy,xmin,xmax,ymin,ymax
  !p.color=color_value(ncolor+1)
    ;!mtitle='Omega(Im) vs ky (kx=0, n=0,+1,-1)'
    !mtitle=strjoin(['Omega(Im) Low freq(kx=0,N_CY=',strtrim(string(fix(N_CY)),2),')'])
      plot,[0],[0],$
       /nodata,$
       Charsize=1.7,$
       xstyle=1,$
       xminor=0,$
       xrange=[xmin,xmax],$
       xtitle='ky',$
       ystyle=1,$
       yminor=0,$
     ;yrange=[0,0.4],$
      yrange=[ymin,ymax],$
       ytitle=' ';,$
       ;/ylog

     oplot,kyhalf,imaginary(Omega_matrFinal(*)),thick=1.5,linestyle=0,Psym=-6,symsize=2.1
     !p.color=color_value(ncolor*11/16)
     oplot,kyhalf,imaginary(Omega_matrcyLow(*)),thick=1.5,linestyle=0,Psym=-2,symsize=2
     !p.color=color_value(ncolor*1/8)
     oplot,kyhalf,imaginary(Omega_matrftLowmult(*,0)),thick=1.5,linestyle=1,Psym=-5,symsize=1.5
     oplot,kyhalf,imaginary(Omega_matrftLowmult(*,1)),thick=1.5,linestyle=2,Psym=-5,symsize=1.5
     oplot,kyhalf,imaginary(Omega_matrftLowmult(*,2)),thick=1.5,linestyle=0,Psym=-5,symsize=1.5
     !p.color=color_value(ncolor+1)
      lines = [0,0,1,2,0]                 ; for line styles
      col=fltarr(5)
      col=[color_value(ncolor+1),color_value(ncolor*11/16),color_value(ncolor*1/8),color_value(ncolor*1/8),color_value(ncolor*1/8) ]
      items =['Gyro-kinetic',strjoin(['Cyclo N_CY=',strtrim(string(fix(N_CY)),2)]),strjoin(['Fourier N_FT=',strtrim(string(fix(N_FTmult(0))),2)])$
      ,strjoin(['Fourier N_FT=',strtrim(string(fix(N_FTmult(1))),2)]),strjoin(['Fourier N_FT=',strtrim(string(fix(N_FTmult(2))),2)])]
      sym = [0]
      legend,items,linestyle=lines,colors=col,charsize=1.5,charthick=1.5,psym=[-6,-2,-5,-5,-5]      ; vertical legend---upper left

  !p.background=color_value(ncolor/2)
  window, plotidh+3, TITLE='main_plot', xsize=600,ysize=650   ;;3. High frequency mode, Real part
  xmax=max(kyhalf)
  xmin=min(kyhalf)
  ymax=max([max(REAL_PART(Omega_matrcyHigh)),max(REAL_PART(Omega_matrftHighmult))])
  ymin=min([min(REAL_PART(Omega_matrcyHigh)),min(REAL_PART(Omega_matrftHighmult))])
  ymin=ymin-(ymax-ymin)/20
  ymax=ymax+(ymax-ymin)/20
  !noeras=0
  !p.color=color_value(ncolor+1)
  set_viewport,0.15,0.95,0.15,0.9
  set_xy,xmin,xmax,ymin,ymax
  !p.color=color_value(ncolor+1)
  ;!mtitle='Omega(Re) vs ky (kx=0, n=0,+1,-1)'
    !mtitle=strjoin(['Omega(Re) High freq(kx=0,N_CY=',strtrim(string(fix(N_CY)),2),')'])
      plot,[0],[0],$
       /nodata,$
       xstyle=1,$
       Charsize=1.7,$
       xminor=0,$
       xrange=[xmin,xmax],$
       xtitle='ky',$
       ystyle=1,$
       yminor=0,$
     ;yrange=[-150,150],$
     yrange=[ymin,ymax],$
       ytitle=' ';,$
       ;/ylog

     !p.color=color_value(ncolor*11/16)
     oplot,kyhalf,REAL_PART(Omega_matrcyHigh(*)),thick=1.5,linestyle=0,Psym=-2,symsize=2
      !p.color=color_value(ncolor*1/8)
     oplot,kyhalf,REAL_PART(Omega_matrftHighmult(*,0)),thick=1.5,linestyle=1,Psym=-5,symsize=1.5
     oplot,kyhalf,REAL_PART(Omega_matrftHighmult(*,1)),thick=1.5,linestyle=2,Psym=-5,symsize=1.5
     oplot,kyhalf,REAL_PART(Omega_matrftHighmult(*,2)),thick=1.5,linestyle=0,Psym=-5,symsize=1.5
     !p.color=color_value(ncolor+1)
      lines = [0,1,2,0]                 ; for line styles
      col=fltarr(4)
      col=[color_value(ncolor*11/16),color_value(ncolor*1/8),color_value(ncolor*1/8),color_value(ncolor*1/8) ]
 items =[strjoin(['Cyclo N_CY=',strtrim(string(fix(N_CY)),2)]),strjoin(['Fourier N_FT=',strtrim(string(fix(N_FTmult(0))),2)]),$
 strjoin(['Fourier N_FT=',strtrim(string(fix(N_FTmult(1))),2)]),strjoin(['Fourier N_FT=',strtrim(string(fix(N_FTmult(2))),2)])]
      sym = [0]
      legend,items,linestyle=lines,colors=col,charsize=1.5,charthick=1.5,psym=[-2,-5,-5,-5]              ; vertical legend---upper left

  !p.background=color_value(ncolor/2)
  window, plotidh+4, TITLE='main_plot', xsize=600,ysize=650   ;;4. High frequency mode, imaginary part
  xmax=max(kyhalf)
  xmin=min(kyhalf)
  ymax=max([max(imaginary(Omega_matrcyHigh)),max(imaginary(Omega_matrftHighmult))])
  ymin=min([min(imaginary(Omega_matrcyHigh)),min(imaginary(Omega_matrftHighmult))])
  ;ymin=ymin-(ymax-ymin)/5
  ;ymax=ymax+(ymax-ymin)/1
  !noeras=0
  !p.color=color_value(ncolor+1)
  set_viewport,0.15,0.95,0.15,0.9
  set_xy,xmin,xmax,ymin,ymax
  !p.color=color_value(ncolor+1)
    ;!mtitle='Omega(Im) vs ky (kx=0, n=0,+1,-1)'
    !mtitle=strjoin(['Omega(Im) High freq(kx=0,N_CY=',strtrim(string(fix(N_CY)),2),')'])
      plot,[0],[0],$
       /nodata,$
       Charsize=1.6,$
       xstyle=1,$
       xminor=0,$
       xrange=[xmin,xmax],$
       xtitle='ky',$
       ystyle=1,$
       yminor=0,$
      yrange=[ymin,ymax],$
      ;yrange=[ymin,ymax],$
       ytitle=' ';,$
       ;/ylog

     !p.color=color_value(ncolor*11/16)
     oplot,kyhalf,imaginary(Omega_matrcyHigh(*)),thick=1.5,linestyle=0,Psym=-2,symsize=2
     !p.color=color_value(ncolor*1/8)
     oplot,kyhalf,imaginary(Omega_matrftHighmult(*,0)),thick=1.5,linestyle=1,Psym=-5,symsize=1.5
     oplot,kyhalf,imaginary(Omega_matrftHighmult(*,1)),thick=1.5,linestyle=2,Psym=-5,symsize=1.5
     oplot,kyhalf,imaginary(Omega_matrftHighmult(*,2)),thick=1.5,linestyle=0,Psym=-5,symsize=1.5
     !p.color=color_value(ncolor+1)
      lines = [0,1,2,0]                 ; for line styles
      col=fltarr(4)
      col=[color_value(ncolor*11/16),color_value(ncolor*1/8),color_value(ncolor*1/8),color_value(ncolor*1/8) ]
items =[strjoin(['Cyclo N_CY=',strtrim(string(fix(N_CY)),2)]),strjoin(['Fourier N_FT=',strtrim(string(fix(N_FTmult(0))),2)]),$
strjoin(['Fourier N_FT=',strtrim(string(fix(N_FTmult(1))),2)]),strjoin(['Fourier N_FT=',strtrim(string(fix(N_FTmult(2))),2)])]
      sym = [0]
      legend,items,linestyle=lines,colors=col,charsize=1.5,charthick=1.5,psym=[-2,-5,-5,-5]      ; vertical legend---upper left
    end

    29: begin ;total_Dft
    nametemp=['total D (CKinFH)','D','total_Dft.txt']
    twoD_time_trace_see,nametemp,29
    end

    30: begin ;total_Chift
    nametemp=['total Chi (CKinFH)','Chi','total_Chift.txt']
    twoD_time_trace_see,nametemp,30
    end

    31: begin ;total_Entropyft
    nametemp=['total Entropy (CKinFH)','Entropy','tot_Entropyft.txt']
    twoD_time_trace_see,nametemp,31
    end

    32: begin ;plot total Phi^2 vs time
    nametemp=['total|Phi| (CKinFH)','|Phi|','Phi_kft.txt']
    ;nametemp=['Phi^2 (3DGK)','|Phi| ^2','total_Phi^2.txt']
    twoD_time_trace_see,nametemp,32
    end


   33: begin ;plot time average D at saturated state --3D plot of D(CKinFH)
         ; default window
  plotidh=number_plot+1
   !p.background=color_value(ncolor/2)
  window, plotidh, TITLE='main_plot', xsize=600,ysize=600
    Dftd=fltarr(2*pmax-1,2*nmax-1,ntmax/output_step+1)
Dftdsat=fltarr(2*pmax-1,nmax)
openr,lun,'D_kft.txt',/get_lun
readf,lun,Dftd
free_lun,lun
For j=sattime/(output_step*tstep),ntmax/output_step Do Begin
Dftdsat(*,*)=Dftdsat(*,*)+Dftd(*,nmax-1:2*nmax-2,j)
Endfor
Dftdsat(*,*)=Dftdsat(*,*)*(1.0/(ntmax/output_step-sattime/(output_step*tstep)))
      x=kx
      y=kyhalf
      !mtitle='D(CKinFH) '+strtrim(round(sattime),2)+'-'+strtrim(round(ntmax*tstep),2)+'Ln/Cs'
      surface3D,Dftdsat,x,y,'kx','ky',0,1,overplot_log,Az,Ax,0.00000001
    end

    34: begin ;plot time average D at saturated state --3D plot of Chi(CKinFH)
         ; default window
      plotidh=number_plot+2
      !p.background=color_value(ncolor/2)
      window, plotidh, TITLE='main_plot', xsize=600,ysize=600

      Chiftsat=fltarr(2*pmax-1,nmax)
      jump=1
      IF(ntmax/output_step GE 1000)then begin
      jump=10
      Endif Else IF(ntmax/output_step GE 10000)then begin
      jump=100
      Endif Else IF(ntmax/output_step GE 100000)then begin
      jump=1000
      Endif Else IF(ntmax/output_step GE 1000000)then begin
      jump=10000
      Endif
      jump=1
      Chi_kftnew=fltarr(2*pmax-1,2*nmax-1,ntmax/(output_step*jump)+1)
      Chi_ktmp=fltarr(2*pmax-1,2*nmax-1)

      openr,lun,'Chi_kft.txt',/get_lun
      readf,lun,Chi_ktmp
      POINT_LUN, -lun, pos
      Print,'Pos=',pos
      For i=0L,ntmax/output_step,jump Do begin
      POINT_LUN, lun, pos*i
      readf,lun,Chi_ktmp
      Chi_kftnew(*,*,i/jump)=Chi_ktmp
      Endfor
      free_lun,lun

For j=round(sattime/(output_step*jump*tstep)),round(ntmax/(output_step*jump))-1 Do Begin
Chiftsat(*,*)=Chiftsat(*,*)+abs(Chi_kftnew(*,nmax-1:2*nmax-2,j))
Endfor
Chiftsat(*,*)=Chiftsat(*,*)*(1.0/(ntmax/(output_step*jump)-sattime/(output_step*jump*tstep)))
      x=kx
      y=kyhalf
      !mtitle=cgGreek('chi')+'i (CKinFH) '+strtrim(round(sattime),2)+'-'+strtrim(round(ntmax*tstep),2)+'Ln/Cs'
      surface3D,Chiftsat,x,y,'kx','ky',0,1,overplot_log,Az,Ax,Chi_logmin
    end

    35: begin ;plot time average |Phi| at saturated state --3D plot of |Phi|(CKinFH)
         ; default window
      plotidh=number_plot+3
      !p.background=color_value(ncolor/2)
      window, plotidh, TITLE='main_plot', xsize=600,ysize=600
      Phiabssat=fltarr(2*pmax-1,nmax)

      jump=1
      IF(ntmax/output_step GE 1000)then begin
      jump=10
      Endif Else IF(ntmax/output_step GE 10000)then begin
      jump=100
      Endif Else IF(ntmax/output_step GE 100000)then begin
      jump=1000
      Endif Else IF(ntmax/output_step GE 1000000)then begin
      jump=10000
      Endif

      Phi_kftnew=complexarr(2*pmax-1,2*nmax-1,ntmax/(output_step*jump)+1)
      Phi_ktmp=complexarr(2*pmax-1,2*nmax-1)

      openr,lun,'Phi_kft.txt',/get_lun
      readf,lun,Phi_ktmp
      POINT_LUN, -lun, pos
      Print,'Pos=',pos
      For i=0,ntmax/(output_step*jump) Do begin
      POINT_LUN, lun, pos*i*jump
      readf,lun,Phi_ktmp
      Phi_kftnew(*,*,i)=Phi_ktmp
      Endfor
      free_lun,lun

For j=round(sattime/(output_step*jump*tstep)),round(ntmax/(output_step*jump)) Do Begin
Phiabssat(*,*)=Phiabssat(*,*)+abs(Phi_kftnew(*,nmax-1:2*nmax-2,j))
Endfor
Phiabssat(*,*)=Phiabssat(*,*)*(1.0/(ntmax/(output_step*jump)-sattime/(output_step*jump*tstep)))
      x=kx
      y=kyhalf
      !mtitle='|'+cgGreek('Phi')+'|'+' (CKinFH) '+strtrim(round(sattime),2)+'-'+strtrim(round(ntmax*tstep),2)+'Ln/Cs'
      surface3D,Phiabssat,x,y,'kx','ky',0,1,overplot_log,Az,Ax,Phi_logmin
    end

    36: begin ;(CKinFH)|Phi|^2 freq spectrum
    nametemp=['(CKinFH)|Phi|^2 freq spectrum','Phi_kft.txt','ky']
    ;;nametemp include: title, ytitle, and input data file name
     kspectrum1D_see,nametemp,36
    end

    37: begin ;(CKinFH)|Phi|^2 freq spectrum
    nametemp=['(GK)|Phi|^2 freq spectrum','Phi_k.txt','ky']
    ;;nametemp include: title, ytitle, and input data file name
     kspectrum1D_see,nametemp,37
    end

    38: begin
    nametemp=['Omegacy(ky)','ky','Phi_kcy.txt']
    ;;nametemp include: title, ytitle, and input data file name
     kspectrum1D_see,nametemp,38
    end

    39: begin  ;; plot Chi_mu (GK)
    nametemp=['Chi(GK)','Chi_mu.txt']
    plot3D_mu_time_see,nametemp,39
    end

    40: begin  ;; plot Chi_mu (CKinFH)
    nametemp=['Chi(CKinFH)','Chi_muft.txt']
    plot3D_mu_time_see,nametemp,40
    end

    41: begin  ;; plot F_mu (GK)
    nametemp=['F (GK)','Fk_mu.txt']
    plot3D_mu_time_see,nametemp,41
    end

    42: begin  ;; plot F_muft (CKinFH)
    nametemp=['F (CKinFH)','Fk_ftmu.txt']
    plot3D_mu_time_see,nametemp,42
    end

    43: begin  ;; plot F_tot_muft (CKinFH)
    nametemp=['F_tot (CKinFH)','Fk_ftmuT.txt']
    plot3D_mu_time_see,nametemp,43
    end
    
    44: begin  ;;
    Omega_matrft=complexarr(2*pmax-1,nmax,N_mu*(2*N_FT-1))
    Ratave=fltarr(nmax-1,ntmax)
    delRat=fltarr(nmax-1,ntmax)
    Time2=findgen(ntmax/output_step+1)*(tstep*output_step)
    
    openr,lun,'Omega_matrft.txt',/get_lun   ;;0 Gyro-kinetic
    readf,lun,Omega_matrft
    n=eof(lun)
    if n ne 1 then print,'error with file load!!!!!!'
    free_lun,lun    
    openr,lun,'Ratave.txt',/get_lun   ;;0 Gyro-kinetic
    readf,lun,Ratave
    n=eof(lun)
    if n ne 1 then print,'error with file load!!!!!!'
    free_lun,lun 
    openr,lun,'delRat.txt',/get_lun   ;;0 Gyro-kinetic
    readf,lun,delRat
    n=eof(lun)
    if n ne 1 then print,'error with file load!!!!!!'
    free_lun,lun       
    
      plotidh=number_plot
      !p.background=color_value(ncolor/2)
      window, plotidh, TITLE='main_plot', xsize=700,ysize=500
      Phiabssat=fltarr(2*pmax-1,nmax) 
      
      ;rangeval= (abs(min(Ratave(Nk,*))/2+max(Ratave(Nk,*)))/2)/6 
      ;ymax=rangeval
      ;ymin=-rangeval
      ymin=0.4
      ymax=0.8
      
   cgplot,[0],[0],$
       /nodata,$
       Charsize=2,$
       xrange=[0,max(Time2)],$
       yrange=[ymin,ymax],$ 
       xtitle='time (a/c'+'$\downs$)',$
       ytitle='$\gamma$',$
       title='$\gamma$ [CKinCH] kx=0 ky='+strtrim((Nk+1)*kymax/(nmax-1),2)
       
   !p.color=color_value(ncolor*1/16)  ;;light blue
   For i=0,ntmax-1 do begin
   oplot,[time2(i),time2(i)],[Ratave(Nk,i)-abs(delRat(Nk,i)),Ratave(Nk,i)+abs(delRat(Nk,i))],psym=-3,symsize=4,thick=0.7
   Endfor
   
   !p.color=color_value(ncolor+1)  ;black 
   oplot,time2,Ratave(Nk,*),thick=2.4
   xyouts,max(Time2)*3/5,(ymax+ymin)*1.0/2,'Init gam='+strtrim(Ratave(Nk,ntmax-1),2)
   
   !p.color=color_value(ncolor*11/16)  ;;red
   oplot,[min(Time2),max(Time2)],[max(imaginary(Omega_matrft(pmax-1,Nk+1,*))),max(imaginary(Omega_matrft(pmax-1,Nk+1,*)))],thick=1.2
   xyouts,max(Time2)*3/5,(ymax+ymin)*0.9/2,'Eigen gam='+strtrim(max(imaginary(Omega_matrft(pmax-1,Nk+1,*))),2)
   
   !p.color=color_value(ncolor+1)  ;black  
    end
        
    45: begin  ;; plot total |Phi| (CK) by time
    nametemp=['Phi_k (CKinCH)','|Phi|','Phi_kcy.txt']
    twoD_time_trace_see,nametemp,45
    end
    
        
    46: begin  ;; plot total |Phi| (CK) by time
    nametemp=['total entropy (CKinCH)','entropy','tot_Entropycy.txt']
    twoD_time_trace_see,nametemp,46
    end    
  
  
    47: begin  ;; combine
    nametemp=['Data processing','',' .txt']
    twoD_time_trace_see,nametemp,47
    end
    
      
    48: begin  
    nametemp=['Others by time','combine',' .txt']
    twoD_time_trace_see,nametemp,48
    end    
    
    49: begin ;; 3D plot of |ne|
      plotidh=number_plot
      !p.background=color_value(ncolor/2)
      window, plotidh, TITLE='main_plot', xsize=600,ysize=600
      nesat=fltarr(2*pmax-1,nmax)

      jump=1
      IF(ntmax/output_step GE 1000)then begin
      jump=10
      Endif Else IF(ntmax/output_step GE 10000)then begin
      jump=100
      Endif Else IF(ntmax/output_step GE 100000)then begin
      jump=1000
      Endif Else IF(ntmax/output_step GE 1000000)then begin
      jump=10000
      Endif

      nsqr_k=complexarr(2*pmax-1,2*nmax-1,ntmax/(output_step*jump)+1)
      nsqrtmp=complexarr(2*pmax-1,2*nmax-1)

      openr,lun,'ne_k.txt',/get_lun
      readf,lun,nsqrtmp
      POINT_LUN, -lun, pos
      Print,'Pos=',pos
      For i=0L,ntmax/output_step,jump Do begin
      POINT_LUN, lun, pos*i
      readf,lun,nsqrtmp
      nsqr_k(*,*,i/jump)=nsqrtmp
      Endfor
      free_lun,lun

For j=round(sattime/(output_step*jump*tstep)),round(ntmax/(output_step*jump)) Do Begin
nesat(*,*)=nesat(*,*)+abs(nsqr_k(*,nmax-1:2*nmax-2,j))
Endfor
nesat(*,*)=nesat(*,*)*(1.0/(ntmax/(output_step*jump)-sattime/(output_step*jump*tstep)))
      x=kx
      y=kyhalf
      !mtitle='|ne|'+' (GK) '+strtrim(round(sattime),2)+'-'+strtrim(round(ntmax*tstep),2)+'Ln/Cs'
      surface3D,nesat,x,y,'kx','ky',0,1,overplot_log,Az,Ax,Phi_logmin
            
    end
    
    50: begin ;plot time average (Phi)^2 of GK at saturated state --3D plot
         ; default window

      plotidh=number_plot
      !p.background=color_value(ncolor/2)
      window, plotidh, TITLE='main_plot', xsize=600,ysize=600
      Phiabssat=fltarr(2*pmax-1,nmax)

      jump=1
      IF(ntmax/output_step GE 1000)then begin
      jump=10
      Endif Else IF(ntmax/output_step GE 10000)then begin
      jump=100
      Endif Else IF(ntmax/output_step GE 100000)then begin
      jump=1000
      Endif Else IF(ntmax/output_step GE 1000000)then begin
      jump=10000
      Endif

      Phi_knew=complexarr(2*pmax-1,2*nmax-1,ntmax/(output_step*jump)+1)
      Phi_ktmp=complexarr(2*pmax-1,2*nmax-1)

      openr,lun,'Phi_k.txt',/get_lun
      readf,lun,Phi_ktmp
      POINT_LUN, -lun, pos
      Print,'Pos=',pos
      For i=0L,ntmax/output_step,jump Do begin
      POINT_LUN, lun, pos*i
      readf,lun,Phi_ktmp
      Phi_knew(*,*,i/jump)=Phi_ktmp
      Endfor
      free_lun,lun

For j=round(sattime/(output_step*jump*tstep)),round(ntmax/(output_step*jump)) Do Begin
;Phiabssat(*,*)=Phiabssat(*,*)+abs(Phi_knew(*,nmax-1:2*nmax-2,j))^2
Phiabssat(*,*)=Phiabssat(*,*)+(Phi_knew(*,nmax-1:2*nmax-2,j))*conj(Phi_knew(*,nmax-1:2*nmax-2,j))
Endfor
Phiabssat(*,*)=Phiabssat(*,*)*(1.0/(ntmax/(output_step*jump)-sattime/(output_step*jump*tstep)))
      x=kx
      y=kyhalf
      !mtitle=''
      ;!mtitle=cgGreek('Phi')+'^2'+' (GK) '+strtrim(round(sattime),2)+'-'+strtrim(round(ntmax*tstep),2)+'Ln/Cs'
      surface3D,Phiabssat,x,y,'kx','ky',0,1,overplot_log,Az,Ax,Phi_logmin


    end    
      
      
    51: begin ;plot frequencyNL/OMGi GK at saturated state --3D plot
         ; default window

      plotidh=number_plot
      !p.background=color_value(ncolor/2)
      window, plotidh, TITLE='main_plot', xsize=600,ysize=600
      freqRatio=fltarr(2*pmax-1,nmax)


      inputvab=fltarr(2*pmax-1,2*nmax-1)

      openr,lun,'freqRatio.txt',/get_lun
      readf,lun,inputvab
      free_lun,lun
      
      freqRatio(*,*)=inputvab(*,nmax-1:2*nmax-2)
      

      x=kx
      y=kyhalf
      !mtitle=''
      ;!mtitle=cgGreek('Phi')+'^2'+' (GK) '+strtrim(round(sattime),2)+'-'+strtrim(round(ntmax*tstep),2)+'Ln/Cs'
      surface3D,freqRatio/Omega_star,x,y,'kx','ky',0,1,overplot_log,Az,Ax,Phi_logmin

    Print,'a=',complex(1,1)
    Print,'a^2=',complex(1,1)^2
    Print,'abs(a)^2=',abs(complex(1,1))^2
  
    end    
                
  endcase
end

;*******************************************************************************
pro eplot,x,y,sigyup,sigylo,_extra=_extra,barlinestyle=barlinestyle, $
          color=color,linestyle=linestyle,thick=thick,noclip=noclip, $
          t3d=t3d
;+
; NAME:
;
;       EPLOT
;
; PURPOSE:
;
;       Plot x vs y, with vertical error bars on y.
;
; CALLING SEQUENCE:
;
;       EPLOT,Y,SIGY
;       EPLOT,X,Y,SIGY
;       EPLOT,Y,SIGY_UP,SIGY_DOWN
;       EPLOT,X,Y,SIGY_UP,SIGY_DOWN
;
; INPUTS:
;
;       X, Y -  1-D arrays
;
;       SIGY - Uncertainty in Y, i.e. Y+/-SIGY
;
;       SIGY_UP, SIGY_DOWN - +/- uncertainties in Y, i.e.,
;                           Y +SIGY_UP -SIGY_DOWN
;
; KEYWORD PARAMETERS:
;
;       BARLINESTYLE = Linestyle for error bars.
;
;               plus all valid IDL plot keywords.  Only the COLOR,
;               THICK, NOCLIP, and T3D keywords apply to the error
;               bars.
;
; MODIFICATION HISTORY:
;
;      D. L. Windt, Bell Laboratories, November 1989
;      Replaced specific plot/oplot keywords with _EXTRA,
;      April, 1997
;
;      windt@bell-labs.com
;-
on_error,2

if n_params() lt 3 then message,'Usage: EPLOT,X,Y,SIGY'

if n_elements(color) eq 0 then color=!p.color
if n_elements(linestyle) eq 0 then linestyle=!p.linestyle
if n_elements(thick) eq 0 then thick=!p.thick
if n_elements(noclip) eq 0 then noclip=!p.noclip
if n_elements(t3d) eq 0 then t3d=!p.t3d
if n_elements(barlinestyle) eq 0 then barlinestyle=linestyle

plot,x,y,_extra=_extra, $
  color=color,linestyle=linestyle,thick=thick,noclip=noclip,t3d=t3d
psym=4 ;!p.psym
!p.psym=0
xt=fltarr(2)
xb=fltarr(2)
yt=xt
if n_params() eq 3 then sigylo=sigyup
for i=0,n_elements(x)-1 do begin
    xt(0)=x(i)
    xt(1)=x(i)
    yt(0)=y(i)
    yt(1)=y(i)+sigyup(i)
    oplot,xt,yt, $
      color=color,linestyle=barlinestyle,thick=thick,noclip=noclip,t3d=t3d
;    xb(0)=x(i);-(x(n_elements(x)-1)-x(0))/60
;    xb(1)=x(i)+(x(n_elements(x)-1)-x(0))/60
;    oplot,xb,yt(1), $
;      color=color,linestyle=barlinestyle,thick=thick,noclip=noclip,t3d=t3d
    yt(1)=y(i)-sigylo(i)
    oplot,xt,yt, $
      color=color,linestyle=barlinestyle,thick=thick,noclip=noclip,t3d=t3d
;    oplot,xb,yt(1), $
;      color=color,linestyle=barlinestyle,thick=thick,noclip=noclip,t3d=t3d
endfor
;   plot,x(3),y(3),psym=4,symsize=2.5
!p.psym=psym
return
end
;*******************************************************************************
;*******************************************************************************
;
;pro surface3D1,z,x,y,x_title,y_title,nbegin,logtezoomzemp
;
;  common startup,number_plot,fpath,ncolor,color_value,plotid
;
;  zmax=max(z)
;  zmin=min(z)
;  xmax=max(x)
;  xmin=min(x)
;  ymax=max(y)
;  ymin=min(y)
;;  if !D.name eq 'X' then wset,plotid
;  !noeras=0
;  !p.color=color_value(ncolor+1)
;;  set_viewport,0.15,0.95,0.15,0.9
;;  set_xy,xmin,xmax,ymin,ymax
;;
;  if logtemp eq 1 then surface,z,x,y,$
;       Skirt=0,xtitle=x_title,$
;       xrange=[xmin,xmax],$
;       yrange=[ymin,ymax],$
;       zrange=[zmin,zmax]/zoomtemp,$
;       ytitle=y_title,charsize=2.0,/zlog
;
;  if logtemp eq 0 then surface,z,x,y,$
;       Skirt=0,xtitle=x_title,$
;       xrange=[xmin,xmax],$
;       yrange=[ymin,ymax],$
;       zrange=[zmin,zmax]/zoomtemp,$
;       ytitle=y_title,charsize=2.0
;
;  return
;end
;********************************************************************************

pro surface3D,z,x,y,x_title,y_title,nbegin,ntime,logtemp,Aztemp,Axtemp,zminlog

  common startup,number_plot,fpath,ncolor,color_value,plotid

  zmax=max(z)
  zmin=min(z)
  xmax=max(x)
  xmin=min(x)
  ymax=max(y)
  ymin=min(y)
;
;  if !D.name eq 'X' then wset,plotid
  !noeras=0
  !p.color=color_value(ncolor+1)
;  set_viewport,0.15,0.95,0.15,0.9
;  set_xy,xmin,xmax,ymin,ymax
;
  if logtemp eq 1 then surface,z,x,y,Skirt=0,$
     xrange=[xmin,xmax],$
     yrange=[ymin,ymax],$
     zrange=[zminlog,zmax],$
     Az=Aztemp,Ax=Axtemp,$
     xtitle=x_title,ytitle=y_title,charsize=3.0,thick=1.0,/zlog
  if logtemp eq 0 then surface,z,x,y,Skirt=0,$
     xrange=[xmin,xmax],yrange=[ymin,ymax],$
     zrange=[zmin,zmax],$
     Az=Aztemp,Ax=Axtemp,$
     xtitle=x_title,ytitle=y_title,charsize=3.0,thick=1.0

  return
end

;*******************************************************************************

pro plot1line,x,y,x_title,ytitle,nbegin,output_step,ntmax,tstep,sattime,i_ptype,logtemp,zoomtemp,yrange_style

  common startup,number_plot,fpath,ncolor,color_value,plotid
  common plot_variables,y_bin,pdf,nbin,d_y


;print,'zoomtemp=',zoomtemp
;print,'ntmax=',ntmax
;print,'tstep=',tstep
;print,'x_title=',x_title
;print,'ytitle=',ytitle
;print,'nbegin=',nbegin
;print,'output_step=',output_step
;print,'yrange_style=',yrange_style

IF(!aveT eq 0)then begin
  xmax=max(x)
  xmin=min(x)
  ;ymin=0
  ymin=min(y)

  ix1=sattime + !tplus*(ntmax/40)*tstep - !tminus*(ntmax/60)*tstep
  if (ix1 ge (ntmax)*tstep or ix1 le 0) then begin
  !tplus=0
  !tminus=0
  ix1=sattime
  endif
  ix2=max(x)
  xxx1=round(ix1/(output_step*tstep))
  xxx2=round(ix2/(output_step*tstep))
  ymax=max(y(xxx1:xxx2))+max(y(xxx1:xxx2))/10.0


 !point1=sattime
 !point4=max(x)
 !point2=(sattime+max(x))/2.0
 !point3=(sattime+max(x))/2.0


 ; plotid=number_plot
 ; print,'plotid========',plotid
 ;if !D.name eq 'X' then wset,plotid
  !noeras=0
  !p.color=color_value(ncolor+1)
  set_viewport,0.15,0.95,0.2,0.9
  set_xy,xmin,xmax,ymin,ymax
  !p.color=color_value(ncolor+1)

   if i_ptype eq 0 then begin
     ;; Diffusion TIME trace
   yrangetemp=Dblarr(2)
   if yrange_style ge 3  then begin
   yrange_style=yrange_style-3
   endif
   if yrange_style eq 1  then begin  ;;plot the middle part of the data
   yrangetemp=[0, 1.1*ymax/zoomtemp]
   endif
   if yrange_style eq 2 then begin  ;;plot from the min one of the data
   yrangetemp=[0.3*ymin/zoomtemp,1.1*ymax/zoomtemp]
   endif
   if yrange_style eq 0 then begin  ;;plot upto the max one of the data
   yrangetemp=[0, 1.1*ymax/zoomtemp/zoomtemp*0.8]
   endif

  if logtemp eq 1 then  cgplot,[0],[0],$
       /nodata,$
       xstyle=1,$
      ; xminor=0,$
       xrange=[xmin,xmax],$
;       xrange=[min(x),max(x)],$
       xtitle=x_title,$
       ystyle=1,$
       yminor=0,$
     yrange=yrangetemp,$
     title=!mtitle,$
       ytitle=ytitle, $
       /ylog
;       color=line

  if logtemp eq 0 then  cgplot,[0],[0],$
       /nodata,$
       xstyle=1,$
       ;xminor=0,$
       xrange=[min(x),max(x)],$
       xtitle=x_title,$
       ystyle=1,$
       yminor=0,$
       yrange=yrangetemp,$
       title=!mtitle,$
       ytitle=ytitle

    ;; Diffusion trace
;     !p.color=color_value(ncolor+1)
     cgplot,x,y,thick=1.2,/overplot;,color=color_vec[0]


     ;; Average line
     ; ix1=delta_t
      ix2=max(x)
      xxx1=round(ix1/(output_step*tstep))
      xxx2=round(ix2/(output_step*tstep))
      xxx1_2=xxx2-xxx1+1
      sat_ave=total(y(xxx1:xxx2))/xxx1_2
      ysat_ave=fltarr(2)
      ysat_ave[*]=sat_ave

delta_Sat=sqrt(total((y(xxx1:xxx2)-sat_ave)^2)/xxx1_2 )
Intermittency=delta_Sat/sat_ave
print, ' from',ix1,'s to',ix2,'s       sat_ave & delta_Sat=',sat_ave,delta_Sat
print, ' from',ix1,'s to',ix2,'s Intermittency=', Intermittency
print,' '

   !p.color=color_value(ncolor*11/16) ;red
     cgplot,[ix1,ix2],ysat_ave,thick=2.5,Color='red',/overplot
   xyouts,ix1,ysat_ave+(yrangetemp[1]-yrangetemp[0])/60,ysat_ave,charthick=3.8,size=2.6  ;+'+/-'+strtrim(delta_Sat,2)

  !p.color=color_value(ncolor+1)
   xyouts,xmin,yrangetemp[1]-(yrangetemp[1]-yrangetemp[0])/15,'Intermittency='+strtrim(Intermittency,2)
     ;; RMS deviation bars

      ixmid=ix1+(ix2-ix1)/2.0
      xxxmid=round(ixmid/(output_step*tstep))
      xxx1_mid=xxxmid-xxx1+1
      xxxmid_2=xxx2-xxxmid+1
      sat_ave_front=total(y(xxx1:xxxmid))/xxx1_mid
      sat_ave_back=total(y(xxxmid:xxx2))/xxxmid_2
      ysat_ave_front=fltarr(2)
      ysat_ave_front[*]=sat_ave_front
      ysat_ave_back=fltarr(2)
      ysat_ave_back[*]=sat_ave_back
     !p.color=color_value(ncolor*1/8) ;green
     cgplot,[ix1,ixmid],ysat_ave_front,thick=2.0,Color='green',/overplot
   xyouts,ix1,ysat_ave_front-(yrangetemp[1]-yrangetemp[0])/20,strtrim(ysat_ave_front,2),charthick=3.8,size=2.6
      !p.color=color_value(ncolor*1/8)
     cgplot,[ixmid,ix2],ysat_ave_back ,thick=2.0,Color='green',/overplot
   xyouts,ixmid,ysat_ave_back+(yrangetemp[1]-yrangetemp[0])/70,ysat_ave_back,charthick=3.8,size=2.6
     !p.color=color_value(ncolor+1)  ;black



      endif else begin
  ;; Diffusion PDF histogram

      nbin=15  ;set default value
      it1=max(x)/2.0+ix1
      it2=max(x)

     pdf_statistics,x,y,it1,it2,tstep,output_step;,y_bin,pdf

     y_bin_plot = fltarr(2*nbin)
     pdf_plot   = fltarr(2*nbin)

     i  = indgen(nbin)
     dy = y_bin[1]-y_bin[0]

     y_bin_plot[2*i]   = y_bin[i]-0.5*dy
     y_bin_plot[2*i+1] = y_bin[i]+0.5*dy
     pdf_plot[2*i]     = pdf[i]
     pdf_plot[2*i+1]   = pdf[i]

     xmin = min(y_bin_plot)
     xmax = max(y_bin_plot)
     ymin = min(pdf_plot)
     ymax = max(pdf_plot)

     yave = total(y_bin[*]*pdf[*])
  if logtemp eq 1 then  plot,[0],[0],$
       /nodata,$
;       title=title,$
       xstyle=1,$
       xminor=0,$
       xrange=[xmin,xmax],$
       xtitle=ytitle,$
       ystyle=1,$
       yminor=0,$
       yrange=[ymin,ymax],$
       ytitle='!3Probability Density Function (saturated state)', $
       color=line,$
       /ylog

   if logtemp eq 0 then  plot,[0],[0],$
       /nodata,$
;       title=title,$
       xstyle=1,$
       xminor=0,$
       xrange=[xmin,xmax],$
       xtitle=ytitle,$
       ystyle=1,$
       yminor=0,$
       yrange=[ymin,ymax],$
       ytitle='!3Probability Density Function (saturated state)', $
       color=line


     oplot,y_bin_plot,pdf_plot;,color=color_vec[0]
       !p.color=color_value(ncolor*13/16)
     oplot,yave*[1,1],100*[-1,1],linestyle=1,thick=1.5
     xyouts,yave-(xmax-xmin)/4.0,(ymin+3*ymax)/4.0,'average value of saturate state',size=1.5
  endelse
Endif

;==================================================================================


IF(!aveT eq 1)then begin
  xmax=max(x)
  xmin=min(x)
  ymin=min(y)

 print,'!tplus=',!tplus
 print,'!tminus=',!tminus

 ;get the four points
 print,' '
 print,'Point',!aveP+1
 print,' '
 IF(!aveP eq 0)then begin
  !point1=sattime + !tplus*(ntmax/40)*tstep - !tminus*(ntmax/60)*tstep
  if (!point1 ge !point2 or !point1 le 0) then begin
  !tplus=0
  !tminus=0
  !point1=min(sattime,!point2)
  endif
 Endif
 IF(!aveP eq 1)then begin
  !point2=(sattime+max(x))/2.0 + !tplus*(ntmax/40)*tstep - !tminus*(ntmax/60)*tstep
  if (!point2 ge !point3 or !point2 le !point1) then begin
  !tplus=0
  !tminus=0
  !point2 =max(!point1,min((sattime+max(x))/2.0,!point3))
  endif
 Endif
 IF(!aveP eq 2)then begin
  !point3=(sattime+max(x))/2.0 + !tplus*(ntmax/40)*tstep - !tminus*(ntmax/60)*tstep
  if (!point3 ge !point4 or !point3 le !point2) then begin
  !tplus=0
  !tminus=0
  !point3 =max(!point2,min((sattime+max(x))/2.0,!point4))
  endif
 Endif
 IF(!aveP eq 3)then begin
  !point4=max(x) + !tplus*(ntmax/40)*tstep - !tminus*(ntmax/60)*tstep
  if (!point4 ge max(x) or !point4 le !point3) then begin
  !tplus=0
  !tminus=0
  !point4 =max(x)
  endif
 Endif


  point1=round(!point1/(output_step*tstep))
  point4=round(!point4/(output_step*tstep))
  ymax=max(y(point1:point4))*1.1


 ; plotid=number_plot
 ; print,'plotid========',plotid
 ;if !D.name eq 'X' then wset,plotid
  !noeras=0
  !p.color=color_value(ncolor+1)
  set_viewport,0.15,0.95,0.2,0.9
  set_xy,xmin,xmax,ymin,ymax
  !p.color=color_value(ncolor+1)


   if i_ptype eq 0 then begin
     ;; Diffusion TIME trace
   yrangetemp=Dblarr(2)
   if yrange_style ge 3  then begin
   yrange_style=yrange_style-3
   endif
   if yrange_style eq 1  then begin  ;;plot the middle part of the data
   yrangetemp=[0, 1.1*ymax/zoomtemp]
   endif
   if yrange_style eq 2 then begin  ;;plot from the min one of the data
   yrangetemp=[0.1*ymin/zoomtemp,1.1*ymax/zoomtemp]
   endif
   if yrange_style eq 0 then begin  ;;plot upto the max one of the data
   yrangetemp=[0, 1.1*ymax/zoomtemp*0.8]
   endif

  if logtemp eq 1 then  cgplot,[0],[0],$
       /nodata,$
       xstyle=1,$
       xminor=0,$
       xrange=[xmin,xmax],$
;       xrange=[min(x),max(x)],$
       xtitle=x_title,$
       ystyle=1,$
       yminor=0,$
     yrange=yrangetemp,$
      title=!mtitle,$
       ytitle=ytitle, $
       /ylog
;       color=line

  if logtemp eq 0 then  cgplot,[0],[0],$
       /nodata,$
       xstyle=1,$
       xminor=0,$
       xrange=[min(x),max(x)],$
       xtitle=x_title,$
       ystyle=1,$
       yminor=0,$
     yrange=yrangetemp,$
      title=!mtitle,$
      ytitle=ytitle

    ;; Diffusion trace
;     !p.color=color_value(ncolor+1)
     cgplot,x,y,thick=1.2,/overplot;,color=color_vec[0]



      point1=round(!point1/(output_step*tstep))
      point2=round(!point2/(output_step*tstep))
      num1_2=point2-point1+1
      ave1_2=total(y(point1:point2))/num1_2
      yave1_2=fltarr(2)
      yave1_2[*]=ave1_2
   !p.color=color_value(ncolor*11/16) ;red
   cgplot,[!point1,!point2],yave1_2,thick=1.5,Color='red',/overplot
   xyouts,!point1+(!point2-!point1)/3,ave1_2+(yrangetemp[1]-yrangetemp[0])/60,strtrim(ave1_2,2),charthick=1.5,size=1.6

      mid=!point1+(!point2-!point1)/2.0
      xmid=round(mid/(output_step*tstep))
      num1_mid=xmid-point1+1
      nummid_2=point2-xmid+1
      ave_front=total(y(point1:xmid))/num1_mid
      ave_back=total(y(xmid:point2))/nummid_2
      yave_front=fltarr(2)
      yave_front[*]=ave_front
      yave_back=fltarr(2)
      yave_back[*]=ave_back
      if ave_front ge ave_back then sign=1
      if ave_front le ave_back then sign=-1
     !p.color=color_value(ncolor*1/8) ;green
     cgplot,[!point1,mid],yave_front,thick=1.0,Color='green',/overplot
   xyouts,!point1,yave_front+sign*(yrangetemp[1]-yrangetemp[0])/15,strtrim(ave_front,2),charthick=1.5,size=1.5
      !p.color=color_value(ncolor*1/8)
     cgplot,[mid,!point2],yave_back ,thick=1.0,Color='green',/overplot
   xyouts,mid,yave_back-sign*(yrangetemp[1]-yrangetemp[0])/15,strtrim(ave_back,2),charthick=1.5,size=1.5
     !p.color=color_value(ncolor+1)  ;black


      point3=round(!point3/(output_step*tstep))
      point4=round(!point4/(output_step*tstep))
      num3_4=point4-point3+1
      ave3_4=total(y(point3:point4))/num3_4
      yave3_4=fltarr(2)
      yave3_4[*]=ave3_4
   !p.color=color_value(ncolor*11/16) ;red
   cgplot,[!point3,!point4],yave3_4,thick=1.5,Color='red',/overplot
   xyouts,!point3+(!point4-!point3)/3,ave3_4+(yrangetemp[1]-yrangetemp[0])/60,strtrim(ave3_4,2),charthick=1.5,size=1.6

      mid=!point3+(!point4-!point3)/2.0
      xmid=round(mid/(output_step*tstep))
      num3_mid=xmid-point3+1
      nummid_4=point4-xmid+1
      ave_front=total(y(point3:xmid))/num3_mid
      ave_back=total(y(xmid:point4))/nummid_4
      yave_front=fltarr(2)
      yave_front[*]=ave_front
      yave_back=fltarr(2)
      yave_back[*]=ave_back
      if ave_front ge ave_back then sign=1
      if ave_front le ave_back then sign=-1
     !p.color=color_value(ncolor*1/8) ;green
     cgplot,[!point3,mid],yave_front,thick=1.0,Color='green',/overplot
   xyouts,!point3,yave_front+sign*(yrangetemp[1]-yrangetemp[0])/15,strtrim(ave_front,2),charthick=1.5,size=1.5
      !p.color=color_value(ncolor*1/8)
     cgplot,[mid,!point4],yave_back ,thick=1.0,Color='green',/overplot
   xyouts,mid,yave_back-sign*(yrangetemp[1]-yrangetemp[0])/15,strtrim(ave_back,2),charthick=1.5,size=1.5
     !p.color=color_value(ncolor+1)  ;black

  ;;calculate the intermetancy
  delta_Sat1=sqrt(total((y(point1:point2)-ave1_2)^2)/num1_2 )
  Intermittency1=delta_Sat1/ave1_2

  print, ' from',!point1,'s to',!point2,'s   ave1_2=',ave1_2
  print, ' from',!point1,'s to',!point2,'s Intermittency1=', Intermittency1
  print,' '

  !p.color=color_value(ncolor+1) ;black
  xyouts,xmin+xmax/100,ymax-(yrangetemp[1]-yrangetemp[0])/15,'Intermittency1='+strtrim(Intermittency1,2),charthick=1.2,size=1.2


  delta_Sat2=sqrt(total((y(point3:point4)-ave3_4)^2)/num3_4 )
  Intermittency2=delta_Sat2/ave3_4

  print, ' from',!point3,'s to',!point4,'s   ave3_4=',ave3_4
  print, ' from',!point3,'s to',!point4,'s Intermittency2=', Intermittency2
  print,' '

  !p.color=color_value(ncolor+1) ;black
  xyouts,xmin+xmax/100,ymax-(yrangetemp[1]-yrangetemp[0])/10,'Intermittency2='+strtrim(Intermittency2,2),charthick=1.2,size=1.2

;;----------------------------------------------------------------------------

      endif else begin
  ;; Diffusion PDF histogram

      nbin=15  ;set default value
      it1=max(x)/2.0+ix1
      it2=max(x)

     pdf_statistics,x,y,it1,it2,tstep,output_step;,y_bin,pdf

     y_bin_plot = fltarr(2*nbin)
     pdf_plot   = fltarr(2*nbin)

     i  = indgen(nbin)
     dy = y_bin[1]-y_bin[0]

     y_bin_plot[2*i]   = y_bin[i]-0.5*dy
     y_bin_plot[2*i+1] = y_bin[i]+0.5*dy
     pdf_plot[2*i]     = pdf[i]
     pdf_plot[2*i+1]   = pdf[i]

     xmin = min(y_bin_plot)
     xmax = max(y_bin_plot)
     ymin = min(pdf_plot)
     ymax = max(pdf_plot)

     yave = total(y_bin[*]*pdf[*])
  if logtemp eq 1 then  plot,[0],[0],$
       /nodata,$
;       title=title,$
       xstyle=1,$
       xminor=0,$
       xrange=[xmin,xmax],$
       xtitle=ytitle,$
       ystyle=1,$
       yminor=0,$
       yrange=[ymin,ymax],$
       ytitle='!3Probability Density Function (saturated state)', $
       color=line,$
       /ylog

   if logtemp eq 0 then  plot,[0],[0],$
       /nodata,$
;       title=title,$
       xstyle=1,$
       xminor=0,$
       xrange=[xmin,xmax],$
       xtitle=ytitle,$
       ystyle=1,$
       yminor=0,$
       yrange=[ymin,ymax],$
       ytitle='!3Probability Density Function (saturated state)', $
       color=line


     oplot,y_bin_plot,pdf_plot;,color=color_vec[0]
       !p.color=color_value(ncolor*13/16)
     oplot,yave*[1,1],100*[-1,1],linestyle=1,thick=1.5
     xyouts,yave-(xmax-xmin)/4.0,(ymin+3*ymax)/4.0,'average value of saturate state',size=1.5
  endelse
Endif


  return
end
;*******************************************************************************

pro overplot,x,y,x_title,y_title,xmax,xmin,ymax,ymin,nbegin,logtemp,linestylemum,i,i_begin,zoomtemp,yrange_style,col_sty

  common startup,number_plot,fpath,ncolor,color_value,plotid

;  xmax=max(x)
;  xmin=min(x)
;  ymax=max(y)
;  ymin=min(y)
  IF (i EQ i_begin) THEN BEGIN
;  if !D.name eq 'X' then wset,plotid
  !noeras=0
  !p.color=color_value(ncolor+1)
  set_viewport,0.2,0.95,0.2,0.9
  IF(ymax eq ymin) THEN BEGIN
  ymin=ymin/2
  ymax=ymax*2
  Endif
  set_xy,xmin,xmax,ymin,ymax

   if yrange_style ge 3  then begin
   yrange_style=yrange_style-3
   endif
   if yrange_style eq 1  then begin  ;;plot the middle part of the data
   yrangetemp=[0,ymax*1.2/zoomtemp]   
   endif
   if yrange_style eq 2 then begin  ;;plot from the min one of the data
   yrangetemp=[-(ymax-ymin)/2/zoomtemp,(ymax-ymin)/2/zoomtemp]
   endif
   if yrange_style eq 0 then begin  ;;plot upto the max one of the data
   yrangetemp=[-(ymax-ymin)/2/zoomtemp,ymax*1.05]   
   endif

   ;yrangetemp=[0,0.6]
   
    if logtemp eq 0 then plot,[0],[0],$
    /nodata,$
    xrange=[xmin,xmax],$
    yrange=yrangetemp,$
    ;yrange=[0,1.0],$    
    xtitle=x_title,$
    ytitle=y_title,$
    linestyle=linestylemum,$
    charsize=1.8

    if logtemp eq 1 then plot,[0],[0],$
    /nodata,$
    xrange=[xmin,xmax],$
    yrange=yrangetemp,$
    xtitle=x_title,$
    ytitle=y_title,$
    linestyle=linestylemum,$
    charsize=1.8,$
    /ylog

    if logtemp eq 2 then plot,[0],[0],$
    /nodata,$
    xrange=[xmin,xmax],$
    yrange=yrangetemp,$
    xtitle=x_title,$
    ytitle=y_title,$
    linestyle=linestylemum,$
    charsize=1.8,$
    /xlog,$
    /ylog
  Endif
;  if logtemp eq 1 then oplot,x,y,linestyle=linestylemum
;  if logtemp eq 0 then oplot,x,y,linestyle=linestylemum

IF(col_sty eq 0)then begin   ;;col_sty=0 -- no color
IF(linestylemum LE 5)THEN BEGIN
  oplot,x,y,linestyle=linestylemum,thick=2.5
ENDif ELSE IF ((linestylemum GT 5) && (linestylemum LE 11))THEN BEGIN
linestylemum=linestylemum-6
  oplot,x,y,linestyle=linestylemum,Psym=-2,thick=2.5
ENDif ELSE IF ((linestylemum GT 11) && (linestylemum LE 17))THEN BEGIN
linestylemum=linestylemum-12
  oplot,x,y,linestyle=linestylemum,Psym=-5,thick=2.5
Endif
EndIF

oplot,[xmin,xmax],[0,0],thick=1

IF(col_sty eq 1)then begin    ;;col_sty=1 -- plot three color lines for Energy(DW,ZF,GAM)
  !p.color=color_value(ncolor*(i+1)/8+1)
  IF(i eq 2)then begin
  !p.color=color_value(ncolor*11/16+1)
  endif
  oplot,x,y,linestyle=linestylemum,thick=2.5
  !p.color=color_value(ncolor+1)
Endif Else IF(col_sty eq 2)then begin     ;;col_sty=2 -- plot three color lines for Omega_sim,_the,_matr
  IF(i eq 0)then begin
  !p.color=color_value(ncolor+1)  ;;black
  Endif Else If(i eq 1)then begin
  !p.color=color_value(ncolor*13/16)  ;;red line
  EndIF
  oplot,x,y,linestyle=linestylemum,thick=2.5
  !p.color=color_value(ncolor+1)
Endif Else IF(col_sty eq 3)then begin     ;;col_sty=3 -- for 'Omega_Phi_k','Omega_Cyclor','Omega_3GK'
  IF(i eq 0)then begin
  !p.color=color_value(ncolor+1)  ;;black
  Endif Else If(i eq 1)then begin
  !p.color=color_value(ncolor*1/8)  ;;green line
  Endif Else If(i eq 2)then begin
  !p.color=color_value(ncolor*13/16)  ;;red line
  EndIF
  oplot,x,y,linestyle=linestylemum,thick=2.5
  !p.color=color_value(ncolor+1)
Endif

  return
end
;********************************************************************************
pro movie_event,event

  common startup,number_plot,fpath,ncolor,color_value
  common movie,plotidm,nvelocity,ngrid,ntime,delv,delx,delt,dataf,wtime,nt1,nt2,nv1,nv2,nf

  widget_control,event.id,get_uvalue=choice

  case choice of
    0: begin
      wdelete, plotidm        ; Closes plot windows
      widget_control, event.top,/destroy
    end

    1: begin
      print,'input 1 for delta_f/f_0, 2 for delta_f, 3 for log(f)'
      read,nf
    end

    2: begin
      print,'current waiting time=',wtime
      print,'input new waiting time in second'
      read,wtime
    end

    3: begin
      print,'current frame number, start=',nt1,'end=',nt2
      print,'input start and end frame number [0,',ntime-1,']'
      read,nt1,nt2
    end

    4: begin
      print,'current velocity grid number, lower=',nv1,'upper=',nv2
      print,'input lower and upper velocity grid number [0,',nvelocity-1,']'
      read,nv1,nv2
    end

    5: begin

      number_plot=number_plot+1
      plotidm=plotidm+1
      !p.background=color_value(ncolor/2)
      !noeras=0
      window, plotidm, TITLE='movie', xsize=600,ysize=600

      ; setup grids
      time=delt*(nt1+indgen(nt2-nt1+1))
      xgrid=delx*indgen(ngrid)
      vgrid=delv*(0.5+nv1+indgen(nv2-nv1+1))
      vmaxwel=exp(-0.5*vgrid*vgrid)

      xmax=max(vgrid)
      xmin=min(vgrid)
      ymax=max(xgrid)
      ymin=min(xgrid)
      set_viewport,0.05,0.95,0.05,0.95
      set_xy,xmin,xmax,ymin,ymax

      ; particle data = delta_f/f_0
      field=fltarr(ngrid,nt2-nt1+1)
      particle=fltarr(nv2-nv1+1,nt2-nt1+1)
      for i=0,nt2-nt1 do begin
        field(*,i)=dataf(0,*,i+nt1)
        for j=0,nv2-nv1 do begin
          ; plot delta_f/f_0
          if nf eq 1 then particle(j,i)=total(dataf(j+nv1+1,*,i+nt1))/(ngrid)
          ; plot delta_f
          if nf eq 2 then particle(j,i)=vmaxwel(j)*total(dataf(j+nv1+1,*,i+nt1))/(ngrid)
          ; plot log(f(v))
          if nf eq 3 then particle(j,i)=alog(vmaxwel(j)*(1.0e-5+1.0+total(dataf(j+nv1+1,*,i+nt1))/ngrid))
        endfor
      endfor

      ; normalize PDF
      pmax=max(particle)
      pmin=min(particle)
      particle=ymin+0.4*(ymax-ymin)*(1.25+particle/max([pmax,abs(pmin)]))
      print,'range of delta_f(v)=',pmin,pmax

      ; normalize field
      phimax=max(field)
      phimin=min(field)
      field=xmin+0.4*(xmax-xmin)*(1.25+field/max([phimax,abs(phimin)]))
      print,'range of field(x)=',phimin,phimax

      p1dv=fltarr(nv2-nv1+1)
      for i=0,nt2-nt1 do begin

        ; plot zero line
        !mtitle='time='+string(time(i))
        !p.color=color_value(ncolor+1)
        yzero=0.5*(ymax-ymin)+ymin
        plot,[xmin,xmax],[yzero,yzero],/xstyle,/ystyle
        xzero=0.5*(xmax-xmin)+xmin
        oplot,[xzero,xzero],[ymin,ymax]

        ; plot distribution function
        !p.color=color_value(ncolor*3/4+1)
        oplot,vgrid,particle(*,i)

        ; plot potential
        !p.color=color_value(ncolor/4+1)
        oplot,field(*,i),xgrid

        if i eq 0 then cursor,xc,yc,4
        wait,wtime
      endfor

      ; plot initial distribution function
      !p.color=color_value(ncolor+1)
      oplot,vgrid,particle(*,0)

    end

    6: begin

      ; open window
      number_plot=number_plot+1
      plotidm=plotidm+1
       !p.background=color_value(ncolor/2)
      window, plotidm, TITLE='movie', xsize=600,ysize=600

      ; color table
      set_viewport,0.2,0.8,0.85,0.95
      !p.thick=16
      !noeras=0
      !mtitle='color table'
      xmin=-float(ncolor)/2
      xmax=float(ncolor)/2
      ymin=0.0
      ymax=1.0
      set_xy,xmin,xmax,ymin,ymax
      !p.color=color_value(ncolor/2)
      plot,[0.0],[1,1]
      for i=0,ncolor-1 do begin
        xa=xmin+float(i)
        x=[xa,xa]
        y=[0.0,1.0]
        !p.color=color_value(i)
        oplot,x,y
      endfor
      !p.color=color_value(ncolor+1)
      xyouts,0.1,0.9,charsize=2.0,'low',/normal
      xyouts,0.85,0.9,charsize=2.0,'high',/normal
      !p.thick=1

      ; setup grids
      time=delt*(nt1+indgen(nt2-nt1+1))
      xgrid=delx*indgen(ngrid)
      vgrid=delv*(0.5+nv1+indgen(nv2-nv1+1))
      vmaxwel=exp(-0.5*vgrid*vgrid)/sqrt(2.0*3.14159265)

      ; (x,v) range
      xmax=max(vgrid)
      xmin=min(vgrid)
      ymax=max(xgrid)
      ymin=min(xgrid)
      !noeras=1
      set_viewport,0.1,0.9,0.1,0.8
      set_xy,xmin,xmax,ymin,ymax

      ; particle data = delta_f/f_0
      field=fltarr(ngrid,nt2-nt1+1)
      particle=fltarr(nv2-nv1+1,ngrid,nt2-nt1+1)
      for i=0,nt2-nt1 do begin
        for j=0,ngrid-1	do begin
          field(j,i)=dataf(0,j,i+nt1)
          ; plot delta_f/_0
          if nf eq 1 then particle(*,j,i)=dataf(nv1+1:nv2+1,j,i+nt1)
          ; plot delta_f
          if nf eq 2 then particle(*,j,i)=vmaxwel*dataf(nv1+1:nv2+1,j,i+nt1)
          ; plot log(f(x,v))
          if nf eq 3 then begin
            for k=0,nv2-nv1 do begin
              particle(k,j,i)=alog(max([vmaxwel(nv2-nv1),vmaxwel(k)*(1.0+dataf(k+nv1+1,j,i+nt1))]))
            endfor
          endif
        endfor
      endfor

      ; f range
      zmax=max(particle)
      zmin=min(particle)
      print,'range of deltaf(x,v)=',zmin,zmax
      zmax=max([zmax,abs(zmin)])
      zmin=-zmax

      ; plot initial PDF
      !mtitle='time='+string(time(0))
      contour,particle(*,*,0),vgrid,xgrid,nlevels=ncolor,c_colors=color_value(indgen(ncolor)),max_value=zmax,min_value=zmin,/xstyle,/ystyle,/fill
      xyouts,0.04,0.5,charsize=4.0,'x',/normal
      xyouts,0.5,0.04,charsize=4.0,'v',/normal

      ; plot initial potential
      phimax=max(abs(field))
      field=xmin+0.4*(xmax-xmin)*(1.25+field/phimax)
      !p.color=color_value(ncolor+1)
      oplot,field(*,0),xgrid
      xyouts,0.6,0.5,charsize=2.0,'Phi(x)',/normal

      cursor,xc,yc,4
      set_viewport,0.05,0.95,0.05,0.95

      for i=0,nt2-nt1 do begin

        !mtitle='time='+string(time(i))
        !noeras=0
        contour,particle(*,*,i),vgrid,xgrid,nlevels=ncolor,c_colors=color_value(indgen(ncolor)),max_value=zmax,min_value=zmin,/xstyle,/ystyle,/fill

        ; plot potential
        oplot,field(*,i),xgrid

        wait,wtime
      endfor
    end

  endcase
end
;******************************************************************************

pro movie

  common startup,number_plot,fpath,ncolor,color_value
  common movie,plotidm,nvelocity,ngrid,ntime,delv,delx,delt,dataf,wtime,nt1,nt2,nv1,nv2,nf

  ; default window
  plotidm=number_plot

  openr, plotidm, 'movie.out'

  ; # of velocity, spatial grids and time steps
  nvelocity=1
  ngrid=1
  ntime=1
  readf,plotidm,nvelocity,ngrid,ntime,delv,delx,delt

  ; read data
  dataf=fltarr(nvelocity+1,ngrid,ntime)
  readf,plotidm,dataf
  close,plotidm

  wtime=0.0
  nt1=0
  nt2=ntime-1
  nv1=0
  nv2=nvelocity-1
  nf=3

  pname=strarr(7)
  pname=["Exit movie","function","Speed","frame number","velocity range","phi(x,t),f(v,t)","f(x,v,t)"]
  xmenu,pname,BASE=pbase,SPACE=10,TITLE='movie',xpad=20,ypad=20
  widget_control,pbase,/realize
  xmanager,"movie",pbase,/NO_BLOCK

end

;********************************************************************************
pro spectrum,x,yy,nbegin,ntime,plotidh

  common startup,number_plot,fpath,ncolor,color_value

  if !D.name eq 'X' then wset,plotidh
  !p.color=color_value(ncolor+1)

  ; panel 1: mode history of real and imaginary components
  yr=yy(*,0)
  yi=yy(*,1)
  xmax=max(x)
  xmin=min(x)
  ymax=max([yr,yi])
  ymin=min([yr,yi])

  !noeras=0
  !linetype=0
  set_viewport,0.14,0.54,0.55,0.95
  set_xy,xmin,xmax,ymin,ymax
  plot,x,yr
  xyouts,0.1,0.9,charsize=2.0,'real',/normal

  !p.color=color_value(ncolor*3/4)
  oplot,x,yi
  xyouts,0.1,0.86,charsize=2.0,'imaginary',/normal

  ; panel 2: mode amplitude history
  ya=sqrt(yr*yr+yi*yi)
  ymax=max(ya)
  ymin=min(ya)

  !noeras=1
  !p.color=color_value(ncolor+1)
  !mtitle="mode amplitude .vs. t"
  set_viewport,0.59,0.99,0.55,0.95
  set_xy,xmin,xmax,ymin,ymax
  plot,x,ya,/ylog

  xpeak=fltarr(11)
  ypeak=fltarr(11)
  npeak=0

  for i=nbegin+1,ntime-2 do begin
    if ya(i) gt ya(i-1) and ya(i) gt ya(i+1) then begin
      xpeak(npeak)=i
      ypeak(npeak)=ya(i)
      npeak=min([10,npeak+1])
    end
  end
  npeak=min([9,npeak])
  x1=x(xpeak(0))
  x2=x(xpeak(npeak-1))
  y1=ypeak(0)
  y2=ypeak(npeak-1)
  !p.color=color_value(ncolor*3/4)
  oplot,[x1,x2],[y1,y2]

  gamma=alog(y2/y1)/(x2-x1)
  print,'real frequency=',(npeak-1)*3.14159265/(x2-x1),"    growth rate=",gamma

  ;dispersion relation of plasma oscillation from Brunner
  ;k*lambda_D  omega_r   omega_i
  ; print,'theoretical (k=0.2)=1.06398447, -0.00005523'
  ; print,'theoretical (k=0.3)=1.15984650, -0.01262042'
  print,'theoretical (k=0.4)=1.28505698, -0.06612800'
  ; print,'theoretical (k=0.5)=1.41566189, -0.15335949'
  ; print,'theoretical (k=0.6)=1.54570677, -0.26411036'
  ; print,'theoretical (k=0.7)=1.67386598, -0.39240143'
  ; print,'theoretical (k=0.8)=1.79989932, -0.53455237'
  ; print,'theoretical (k=0.9)=1.92386517, -0.68810933'
  ; print,'theoretical (k=1.0)=2.04590486, -0.85133046'

  ; panel 3: mode amplitude normalized by growth rate
  ; gamma=0.0
  yr=yr/exp(gamma*x)
  yi=yi/exp(gamma*x)
  ysize=size(yr)
  mean=total(yr)/ysize(1)
  yr=yr-mean
  ymin=min([yr,yi])
  ymax=max([yr,yi])

  !noeras=1
  !mtitle="mode"
  !p.color=color_value(ncolor+1)
  set_viewport,0.14,0.54,0.05,0.45
  set_xy,xmin,xmax,ymin,ymax
  plot,x,yr
  !p.color=color_value(ncolor*3/4)
  oplot,x,yi

  ; panel 4: power spectrum
  np=(ntime-nbegin)/16
  power=complex(yr,yi)
  power=fft(power,-1)
  ypow=abs(power)
  yp=fltarr(2*np)
  xp=fltarr(2*np)
  for i=0,np-2 do begin
    yp(i)=ypow(i+ntime-nbegin-np+1)
    xp(i)=(i-np+1)*6.283185/(x(ntime)-x(nbegin))
  end
  for i=0, np do begin
    yp(np-1+i)=ypow(i)
    xp(np-1+i)=i*6.283185/(x(ntime)-x(nbegin))
  end
  xmax=max(xp)
  xmin=min(xp)
  ymax=max(yp)
  ymin=min(yp)

  !noeras=1
  !p.color=color_value(ncolor+1)
  !mtitle="frequency spectrum"
  set_viewport,0.59,0.99,0.05,0.45
  set_xy,xmin,xmax,ymin,ymax
  plot,xp,yp,/xstyle

  ; plasma frequency for k=0.5
  k=0.4
  freqr=sqrt(1.0+3.0*k*k)
  freqi=sqrt(3.14159265/2)*freqr/(k*k*k)*exp(-freqr*freqr/(2.0*k*k))

  ; set display back to terminal after printing to PS file
  if !D.name eq 'PS' then begin
    device,/close
    set_plot,'X'
    !p.thick=1
  endif

  return
end
;********************************************************************************
pro plot3D_mu_time_see,nametemp,windownum,group=group

;  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data
  ;;
  common plotaxis,zoomz,zoombackup0,zoombackup1,log,half,freq,i_ptype,satdelta_t,name,yrange_style
  common main_plot,plotidh,ntime,ndata,data_name,windowmark,basemark

  ;;-----------------------------------------------

  ;;------------------------------------------
  ;; Return conditions
  ;;
 ; if exists_diff eq 0 then return
 ; if xregistered('diffusion_ave_see') then return
  ;;------------------------------------------
  name=nametemp
  base = widget_base(title=name[0],$
                     /column)
   ;;set default value
  i_ptype=0
  log=0
  zoomz = 1.0  ;zoomz is the parameter which used to adjust the range of the axis
  satdelta_t=0.0  ;satdelta_t is the parameter which used to adjust the saturation time when the user click the button "t+" or "t-"
  yrange_style=1  ;yrange_style is the parameter which used to choose plot the top, middle or bottom part of the picture.
                  ;for the value of yrange_style/3, 0 for top, 1 for middle, 2 for bottom;
  defsysv,'!Azdelta',0 ;;to rotate the z axis
  defsysv,'!Axdelta',0 ;;to rotate the x axis
  defsysv,'!jumpmul',1.0 ;;to modify the jump  ;; Must get float initial value, like 1.0 !
  defsysv,'!zmintem',1.0 ;;to modify the jump  !zmintem  ;; Must get float initial value, like 1.0 !
  defsysv,'!twoD',0  ;;decide plot 2D or not
  
  ;i_tp = 0
  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row1, $
                    value='Plot', $
                    uvalue=0)

  x = widget_button(row1, $
                    value='t+', $
                    uvalue=1)

  x = widget_button(row1, $
                    value='t-', $
                    uvalue=2)

  x = widget_button(row1, $
                    value='Z+', $
                    uvalue=3)

  x = widget_button(row1, $
                    value='Z-', $
                    uvalue=4)

  x = widget_button(row1, $
                    value='X+', $
                    uvalue=5)

  x = widget_button(row1, $
                    value='X-', $
                    uvalue=6)

  x = widget_button(row1, $
                    value='zmin*10', $
                    uvalue=7)

  x = widget_button(row1, $
                    value='zmin/10', $
                    uvalue=8)

  x = widget_button(row1, $
                    value='jump*2', $
                    uvalue=9)

  x = widget_button(row1, $
                    value='jump/2', $
                    uvalue=10)

  x = widget_button(row1,$
                    value='log',$
                    uvalue=11)

  x = widget_button(row1, $
                    value='2D', $
                    uvalue=12)
                    
  x = widget_button(row1, $
                    value='Done', $
                    uvalue=13)
  ;;----------------------------------------------------------
  ;; DRAW WIDGET and CONTROL
  ;;----------------------------------------------------------

  draw = widget_draw(base,     $
                     xsize=600, $
                     ysize=650)

  widget_control, base, $
    ;set_uvalue=state,$
    /no_copy, $
    /realize

  ;!plotkspecid=!D.WINDOW
    (*windowmark)[windownum]=!D.window
    basemark[windownum]=base

  xmanager,'plot3D_mu_time', $
    base,$
;    event='energy_time_trace_event',$
    group_leader=group


end

;*******************************************************************************
pro plot3D_mu_time_event,event

  common startup,number_plot,fpath,ncolor,color_value,plotid
  common plotaxis,zoomz,zoombackup0,zoombackup1,log,half,freq,i_ptype,satdelta_t,name,yrange_style
  common main_plot,plotidh,ntime,ndata,data_name,windowmark,basemark

  widget_control, event.id, $
    get_uvalue=uvalue
 ; wset, widget

openr,lun,'inputRCYCLO.txt',/get_lun
  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Control_variables',plot_name_temp,19)
  endwhile
  thisFormat='(8x,i2,/,9x,i2,/,4x,i2)'
  readf,lun,restart,GK_FK_CK,CDW,Format=thisFormat
  thisFormat='(8x,i2,/,11x,i2,/,12x,i5,/,11x,i4,/,9x,f5.3,/,/,/,/,7x,i10)'
  readf,lun,muDtype,mugridtype,output_step,backup_num,Const_nl,Stopnt,Format=thisFormat

  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Grid_variables',plot_name_temp,16)
  endwhile
  thisFormat='(6x,i10,/,6x,f12.9,/,6x,f10.3,/,6x,f10.3,/,5x,i5,/,5x,i5,/,6x,i8.3,/,5x,i5,/,5x,i6,/,5x,i6)'
  readf,lun,ntmax,tstep,kxmax,kymax,pmax,nmax,mumax,N_mu,$
      N_FT,N_CY,Format=thisFormat

  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Physics_variables',plot_name_temp,19)
  endwhile
  thisFormat='(11x,f7.4,/,5x,f5.3,/,6x,f6.3,/,9x,f5.3,/,9x,f5.3,/,9x,f5.3,/,7x,f9.3,/,8x,f5.3,/,6x,f5.3,/,6x,f5.3,/,8x,f12.7,/,8x,f8.5,/,/,/,6x,f8.4)'
  readf,lun,Omega_star,a_Ln,a_LTi,lambda_n,lambda_0,lambda_D,AlphaA,delta_1,$
  mu_HK,mu_LK,F_k_int,Epsilon,gamIC, Format=thisFormat

free_lun,lun



openr,lun,'inputvuRCYCLO.txt',/get_lun
  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Control_variables',plot_name_temp,19)
  endwhile
  thisFormat='(7x,i5,/,9x,i5,/,8x,f10.3,/,14x,i2,/,13x,i2,/,3x,i3,/,3x,i3,/,4x,i6,/,4x,i6,/,11x,i4)'
  readf,lun,avenum,plotstep,sattime,Intermit_Plot,overplot_log,Az,Ax,ix1,ix2,Switch_sat,Format=thisFormat

  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Range_control',plot_name_temp,15)
  endwhile
  thisFormat='(10x,f15.12,/,10x,f12.4,/,11x,f15.10,/,11x,f15.12,/,9x,i5,/,8x,f10.4)'
  readf,lun,Omega_min,Omega_max,Phi_logmin,Chi_logmin,OMGpoint,freqPer,Format=thisFormat

free_lun,lun



;  openr,lun,'nt_Odd.txt',/get_lun
;  readf,lun,nt_Odd
;  free_lun,lun
;  openr,lun,'nt_Even.txt',/get_lun
;  readf,lun,nt_Even
;  free_lun,lun
;  IF(ntmax GE max([nt_Odd,nt_Even]))Then Begin
;  ntmax=max([nt_Odd,nt_Even])
;  Endif
;  print,'ntmax=',ntmax

  Time=findgen(ntmax/output_step+1)*(tstep*output_step)
  kx=findgen(2*pmax-1)*kxmax/(pmax-1)-kxmax
  ky=findgen(2*nmax-1)*kymax/(nmax-1)-kymax
  kxhalf=findgen(pmax)*kxmax/(pmax-1)
  kyhalf=findgen(nmax)*kymax/(nmax-1)
  mu=fltarr(N_mu)


  For i=0,ndata-1 Do begin
  IF(basemark[i] eq event.top)then begin
    IF((*windowmark)[i] ne !D.window)then begin
    wset,(*windowmark)[i]
    Endif
    case(i) of
    39: name=[cgGreek('chi')+'i (GK)','Chi_mu.txt','cpu_timeGK.txt']
    40: name=[cgGreek('chi')+'i (CKinFH)','Chi_muft.txt','cpu_timeFK.txt']
    41: name=['F (GK)','Fk_mu.txt','cpu_timeGK.txt']
    42: name=['F (CKinFH)','Fk_ftmu.txt','cpu_timeFK.txt']
    43: name=['F_tot (CKinFH)','Fk_ftmuT.txt','cpu_timeFK.txt']
    endcase
  Endif
  Endfor

  openr,lun,name[2],/get_lun
  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp(' mu_point=',plot_name_temp,10)
  endwhile
  readf,lun,mu
  free_lun,lun
  
  case (uvalue) of

     0:begin
     goto, plot_it
     end

     1:begin
        satdelta_t = satdelta_t+(ntmax/40)*tstep  ;satdelta_t is the variable which used to adjust the saturation time when the user click the button "t+" or "t-"
        if (satdelta_t ge (ntmax)*tstep-sattime) then satdelta_t = 0
        goto, plot_it
     goto, plot_it
     end

     2: begin
        satdelta_t = satdelta_t-(ntmax/40)*tstep
        if (satdelta_t le -sattime) then satdelta_t = 0
        goto, plot_it
     end

     3: begin
        !Azdelta=!Azdelta+10
        goto, plot_it
     end

     4: begin
        !Azdelta=!Azdelta-10
        goto, plot_it
     end

     5: begin
        !Axdelta=!Axdelta+5
        goto, plot_it
     end

     6: begin
        !Axdelta=!Axdelta-5
        goto, plot_it
     end

     7: begin
        !zmintem=!zmintem*10.0
        goto, plot_it
     end

     8: begin
        !zmintem=!zmintem/10.0
        goto, plot_it
     end

     9: begin
        !jumpmul=!jumpmul*2.0
        goto, plot_it
     end

     10: begin
        !jumpmul=!jumpmul/2.0
        goto, plot_it
     end

     11: begin
        log = (log+1) mod 2
        goto, plot_it
     end

     12: begin ;shut down this window
      !twoD=(!twoD+1) mod 2
      goto, plot_it
     end

     13: begin ;shut down this window
;     IF(!plotwindowsid eq 10)then begin
;     Wdelete,!plotwindowsid
;     !plotwindowsid=0
;     Endif
     widget_control, event.top, /destroy
     end
     
  endcase

  return

  plot_it:


  ntnum_of_bin=round(ntmax/(output_step*avenum))
  !p.background=color_value(ncolor/2)
;=================================================================================================================
;      plotidh=number_plot+1
;      !p.background=color_value(ncolor/2)
;      window, plotidh, TITLE='main_plot', xsize=700,ysize=600

      jump=1*!jumpmul
      IF(ntmax/output_step GE 1000)then begin
      jump=100*!jumpmul
      Endif Else IF(ntmax/output_step GE 10000)then begin
      jump=1000*!jumpmul
      Endif Else IF(ntmax/output_step GE 100000)then begin
      jump=10000*!jumpmul
      Endif Else IF(ntmax/output_step GE 1000000)then begin
      jump=100000*!jumpmul
      Endif

      Z=fltarr(N_mu,ntmax/(output_step*jump)+1)
      temp=fltarr(N_mu)

      openr,lun,name[1],/get_lun
      readf,lun,temp
      POINT_LUN, -lun, pos
      Print,'Pos=',pos
      For i=0L,ntmax/output_step,jump Do begin
      POINT_LUN, lun, pos*i
      readf,lun,temp
      Z(*,i/jump)=temp
      Endfor
      free_lun,lun

      sattime=sattime+satdelta_t
      Az=Az+!Azdelta
      Ax=Ax+!Axdelta
      Chi_logmin=Chi_logmin*!zmintem

      x=mu
      y=findgen(ntmax/(output_step*jump)+1)*(tstep*output_step*jump)
      !mtitle=name[0]+' '+strtrim(round(sattime),2)+'-'+strtrim(round(ntmax*tstep),2)+'Ln/Cs'
      surface3D,Z(*,round(sattime/(output_step*jump*tstep)):round(ntmax/(output_step*jump))),x,$
      y(round(sattime/(output_step*jump*tstep)):round(ntmax/(output_step*jump))),cgGreek('mu'),'time',0,1,log,Az,Ax,Chi_logmin

      if(!twoD eq 1)then begin

     !p.background=color_value(ncolor/2)
      window, plotidh+2, TITLE='main_plot', xsize=500,ysize=550
      chimu=fltarr(N_mu)
      chimu=total(Z(*,round(sattime/(output_step*jump*tstep)):round(ntmax/(output_step*jump))),2)/round((ntmax-sattime+1)/(output_step*jump))
      plot,x,chimu,$
      xtitle=cgGreek('mu'),$
      xrange=[min(x),max(x)],$
      yrange=[min(chimu),1.2*max(chimu)],$
       ytitle=name[0], $
       title=!mtitle
      endif
;=================================================================================================================

end
;********************************************************************************
pro Main_plots_name

  common startup,number_plot,fpath,ncolor,color_value,plotid
  common main_plot,plotidh,ntime,ndata,data_name,windowmark,basemark

  defsysv,'!plotkspecid',4   ;;define systerm variable to count the plot times of E(kx) and D(kx)...
  defsysv,'!plotwindowsid',0   ;;defined to mark the windows number

;  ; default window
  plotidh=number_plot
   !p.background=color_value(ncolor/2)
;  window, plotidh, TITLE='main_plot', xsize=600,ysize=600

  openr, plotidh, 'inputvuRCYCLO.txt'
  ; # of time steps and # of data
  ntime=1
  ndata=1
  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,plotidh,plot_name_temp
  result1=strcmp('##Main_plots_name',plot_name_temp,17)
  endwhile
  readf,plotidh,ntime,ndata

  ; names of data
  data_name=strarr(ndata)
  readf,plotidh,data_name
  close, plotidh


  ; widget panel
  hname=strarr(ndata)
  hname(0:ndata-1)=data_name

  windowmark=ptr_new(intarr(200)) ;;create a point array to record the ID of windows.
  basemark=intarr(200)

  xmenu,hname,BASE=hbase,SPACE=10,TITLE='inputvuRCYCLO',column=3,xpad=20,ypad=20
  widget_control,hbase,/realize
  xmanager,"Main_plots",hbase,/NO_BLOCK


end
;********************************************************************************
pro More_details_name

  common startup,number_plot,fpath,ncolor,color_value,plotid
  common main_plot,plotidh,ntime,ndata,data_name,windowmark,basemark

  defsysv,'!plotkspecid',4   ;;define systerm variable to count the plot times of E(kx) and D(kx)...
  defsysv,'!plotwindowsid',0   ;;defined to mark the windows number

;  ; default window
  plotidh=number_plot
   !p.background=color_value(ncolor/2)
;  window, plotidh, TITLE='main_plot', xsize=600,ysize=600

  openr, plotidh, 'inputvuRCYCLO.txt'
  ; # of time steps and # of data
  ntime=1
  ndata=1
  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,plotidh,plot_name_temp
  result1=strcmp('##More_details_name',plot_name_temp,19)
  endwhile
  readf,plotidh,ntime,ndata

  ; names of data
  data_name=strarr(ndata)
  readf,plotidh,data_name
  close, plotidh


  ; widget panel
  hname=strarr(ndata)
  hname(0:ndata-1)=data_name

  windowmark=ptr_new(intarr(200)) ;;create a point array to record the ID of windows.
  basemark=intarr(200)

  xmenu,hname,BASE=hbase,SPACE=10,TITLE='More details',column=3,xpad=20,ypad=20
  widget_control,hbase,/realize
  xmanager,"More_details",hbase,/NO_BLOCK


end
;*******************************************************************************
pro colorful_plot_name

  common startup,number_plot,fpath,ncolor,color_value,plotid
  common colorful_plot,Colndata,Colwindowmark,Colbasemark

;  ; default window
  lun=number_plot
   !p.background=color_value(ncolor/2)

  openr, lun, 'inputvuRCYCLO.txt',/get_lun
  ; # of time steps and # of data
  Colndata=1
  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##colorful_plot_name',plot_name_temp,20)
  endwhile
  readf,lun,Colndata

  ; names of data
  data_name=strarr(Colndata)
  readf,lun,data_name
  free_lun,lun

  ; widget panel
  hname=strarr(Colndata)
  hname(0:Colndata-1)=data_name

  Colwindowmark=ptr_new(intarr(Colndata)) ;;create a point array to record the ID of windows.
  Colbasemark=intarr(Colndata)

  xmenu,hname,BASE=hbase,SPACE=10,TITLE='colorful_plot',column=2,xpad=20,ypad=20
  widget_control,hbase,/realize
  xmanager,"colorful_plot",hbase,/NO_BLOCK


end

;*******************************************************************************

pro colorful_Plot_event,event

  common startup,number_plot,fpath,ncolor,color_value,plotid
  common colorful_plot,Colndata,Colwindowmark,Colbasemark

  widget_control,event.id,get_uvalue=choice

  !p.thick=2

  case choice of

    0: begin
      ;plotidh=1
      ;wdelete,plotidh        ; Closes plot windows
     ; plotidh=!D.window
      while(!D.window NE -1) Do begin
    ;  IF( plotidh GT 0)then Begin
      Wdelete,!D.window  ;;delete all the windows.
     ; plotidh=plotidh-1
      Endwhile
      widget_control, event.top,/destroy
    end

    1: begin ;total energy
    nametemp=['Energy(k)','energy','energy_k_sat.txt']
    ;;nametemp include: title, ytitle, and input data file name
     Colorfulplot_see,nametemp,1
    end

    2: begin ;total energy
    nametemp=['Omega_sim(k)','Omega','Phi_k.txt']
    Colorfulplot_see,nametemp,2
    end

    3: begin ;total energy
    nametemp=['Omega_the(k)','Omega','Omega_the.txt']
    Colorfulplot_see,nametemp,3
    end

    4: begin ;total energy
    nametemp=['Omega_matr(k)','Omega','Omega_matr.txt']
    Colorfulplot_see,nametemp,4
    end

    5: begin ;eddy: phi vs (x,y)
    nametemp=['|ne|','|ne|','ne_k.txt']
    Colorfulplot_see,nametemp,5
    end
    
    6: begin ;eddy: phi vs (x,y)
    nametemp=['|Phi|^2','|Phi|^2','Phi_kOMG.txt']
    Colorfulplot_see,nametemp,6
    end
    
  endcase
end

;*******************************************************************************

pro Colorfulplot_see,nametemp,windownum,group=group

;  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data
  ;;
;  common PRIVATE_time_trace,widget
  common plotaxis,zoomz,zoombackup0,zoombackup1,log,half,freq,i_ptype,satdelta_t,name,yrange_style
  ;common main_plot,plotidh,ntime,ndata,data_name,windowmark,basemark
  common colorful_plot,Colndata,Colwindowmark,Colbasemark

  ;;-----------------------------------------------

  ;;------------------------------------------
  ;; Return conditions
  ;;
 ; if exists_diff eq 0 then return
 ; if xregistered('diffusion_ave_see') then return
  ;;------------------------------------------
  name=nametemp
  base = widget_base(title=name[0],$
                     /column)
   ;;set default value
  defsysv,'!Fignum',1
  defsysv,'!Fignum2',1
  defsysv,'!Omegap',0
  defsysv,'!velocity',0  ;0 for no velocity vector plot, >0 for ploting velocity vector at a certain moment.
  i_ptype=0
  log=0
  zoomz = 1.0  ;zoomz is the parameter which used to adjust the range of the axis
  satdelta_t=0.0  ;satdelta_t is the parameter which used to adjust the saturation time when the user click the button "t+" or "t-"
  yrange_style=1  ;yrange_style is the parameter which used to choose plot the top, middle or bottom part of the picture.
                  ;for the value of yrange_style/3, 0 for top, 1 for middle, 2 for bottom;
  ;i_tp = 0
  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row,$
                     /frame)

  IF(strcmp(name[0],'Omega_',6)) THEN BEGIN
  x = widget_button(row1, $
                    value='Plot(re)', $
                    uvalue=0)

  x = widget_button(row1, $
                    value='Plot(im)', $
                    uvalue=1)
  x = widget_button(row1,$
                    value='Type',$
                    uvalue=2)
  Endif Else Begin
  x = widget_button(row1, $
                    value='Plot', $
                    uvalue=0)
                    
  x = widget_button(row1,$
                    value='log',$
                    uvalue=7)    

  x = widget_button(row1,$
                    value='velocity',$
                    uvalue=3)                                       
  Endelse
;  x = widget_button(row1, $
;                    value='t+', $
;                    uvalue=1)
;
;  x = widget_button(row1, $
;                    value='t-', $
;                    uvalue=2)
;
;  x = widget_button(row1, $
;                    value='ZOOM in', $
;                    uvalue=3)
;
;  x = widget_button(row1, $
;                    value='ZOOM out', $
;                    uvalue=4)
;  x = widget_button(row1, $
;                    value='Yrange', $
;                    uvalue=5)
;  x = widget_button(row1, $
;                    value='units', $
;                    uvalue=6)
;
;  x = widget_button(row1,$
;                    value='log',$
;                    uvalue=7)
;
;  x = widget_button(row1,$
;                    value='TYPE',$
;                    /menu)
;
;  tlevels=['Line Plot','PDF']
;  for i=0,1 do begin
;     x1 = widget_button(x,$
;                        value=tlevels[i],$
;                        uvalue=20+i)
;  endfor
;
;  x = widget_button(row1, $
;                    value='Cal E(k)sat', $
;                    uvalue=9)

  x = widget_button(row1, $
                    value='Done', $
                    uvalue=8)

  ;;----------------------------------------------------------
  ;; DRAW WIDGET and CONTROL
  ;;----------------------------------------------------------

  draw = widget_draw(base,     $
                     xsize=600, $
                     ysize=500)

  widget_control, base, $
    ;set_uvalue=state,$
    /no_copy, $
    /realize

  ;!plotkspecid=!D.WINDOW
    (*colwindowmark)[windownum]=!D.window
    colbasemark[windownum]=base

  xmanager,'Colorfulplot', $
    base,$
;    event='energy_time_trace_event',$
    group_leader=group


end

;*******************************************************************************

pro Colorfulplot_event,event

  common startup,number_plot,fpath,ncolor,color_value,plotid
  common plotaxis,zoomz,zoombackup0,zoombackup1,log,half,freq,i_ptype,satdelta_t,name,yrange_style
  common colorful_plot,Colndata,Colwindowmark,Colbasemark

  widget_control, event.id, $
    get_uvalue=uvalue
 ; wset, widget

openr,lun,'inputRCYCLO.txt',/get_lun
  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Control_variables',plot_name_temp,19)
  endwhile
  thisFormat='(8x,i2,/,9x,i2,/,4x,i2)'
  readf,lun,restart,GK_FK_CK,CDW,Format=thisFormat
  thisFormat='(8x,i2,/,11x,i2,/,12x,i5,/,11x,i4,/,9x,f5.3,/,/,/,/,7x,i10)'
  readf,lun,muDtype,mugridtype,output_step,backup_num,Const_nl,Stopnt,Format=thisFormat

  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Grid_variables',plot_name_temp,16)
  endwhile
  thisFormat='(6x,i10,/,6x,f12.9,/,6x,f10.3,/,6x,f10.3,/,5x,i5,/,5x,i5,/,6x,i8.3,/,5x,i5,/,5x,i6,/,5x,i6)'
  readf,lun,ntmax,tstep,kxmax,kymax,pmax,nmax,mumax,N_mu,$
      N_FT,N_CY,Format=thisFormat

  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Physics_variables',plot_name_temp,19)
  endwhile
  thisFormat='(11x,f7.4,/,5x,f5.3,/,6x,f6.3,/,9x,f5.3,/,9x,f5.3,/,9x,f5.3,/,7x,f9.3,/,8x,f5.3,/,6x,f5.3,/,6x,f5.3,/,8x,f12.7,/,8x,f8.5,/,/,/,6x,f8.4)'
  readf,lun,Omega_star,a_Ln,a_LTi,lambda_n,lambda_0,lambda_D,AlphaA,delta_1,$
  mu_HK,mu_LK,F_k_int,Epsilon,gamIC, Format=thisFormat

free_lun,lun

N_mu=N_mu+1

openr,lun,'inputvuRCYCLO.txt',/get_lun
  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Control_variables',plot_name_temp,19)
  endwhile
  thisFormat='(7x,i5,/,9x,i5,/,8x,f10.3,/,14x,i2,/,13x,i2,/,/,/,4x,i6,/,4x,i6,/,11x,i4)'
  readf,lun,avenum,plotstep,sattime,Intermit_Plot,overplot_log,ix1,ix2,Switch_sat,Format=thisFormat

  plot_name_temp="aaa"
  result1=0
  While (result1 eq 0) Do Begin
  readf,lun,plot_name_temp
  result1=strcmp('##Range_control',plot_name_temp,15)
  endwhile
  thisFormat='(10x,f15.12,/,10x,f12.4,/,11x,f15.10,/,11x,f15.12,/,9x,i5,/,8x,f10.4)'
  readf,lun,Omega_min,Omega_max,Phi_logmin,Chi_logmin,OMGpoint,freqPer,Format=thisFormat

free_lun,lun



;  openr,lun,'nt_Odd.txt',/get_lun
;  readf,lun,nt_Odd
;  free_lun,lun
;  openr,lun,'nt_Even.txt',/get_lun
;  readf,lun,nt_Even
;  free_lun,lun
;  IF(ntmax GE max([nt_Odd,nt_Even]))Then Begin
;  ntmax=max([nt_Odd,nt_Even])
;  Endif
;  print,'ntmax=',ntmax

  Time=findgen(ntmax/output_step+1)*(tstep*output_step)
  kx=findgen(2*pmax-1)*kxmax/(pmax-1)-kxmax
  ky=findgen(2*nmax-1)*kymax/(nmax-1)-kymax
  kxhalf=findgen(pmax)*kxmax/(pmax-1)
  kyhalf=findgen(nmax)*kymax/(nmax-1)

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  ;;set the click botton action act on the window you want.
  For i=0,colndata-1 Do begin
  IF(colbasemark[i] eq event.top)then begin
    IF((*colwindowmark)[i] ne !D.window)then begin
    wset,(*colwindowmark)[i]
    Endif
    case(i) of
    1: name=['Energy(k)','energy','energy_k_sat.txt']
    2: name=['Omega_sim(k)','Omega','Phi_k.txt']
    3: name=['Omega_the(k)','Omega','Omega_the.txt']
    4: name=['Omega_matr(k)','Omega','Omega_matr.txt']
    5: name=['|ne|','|ne|','ne_k.txt']
    6: name=['|Phi|^2','|Phi|^2','Phi_kOMG.txt'] ;;eddy movie
    endcase
  Endif
  Endfor


  case (uvalue) of

;  IF(strcmp(name[0],'Omega_',6)) THEN BEGIN
     0:begin
     !Omegap=0
     goto, plot_it
     end

     1:begin
     !Omegap=1
     goto, plot_it
     end

     2:begin
     !Fignum=!Fignum+1
     if !Fignum gt 3  then begin
     !Fignum=!Fignum-3
     endif
     !Fignum2=!Fignum2+1
     if !Fignum2 gt 7  then begin
     !Fignum2=!Fignum2-7
     endif
     goto, plot_it
     end

     3:begin
     !velocity=!velocity+1
     goto, plot_it
     end
          
     7: begin
        log = (log+1) mod 2
        goto, plot_it
     end
          
     8: begin ;shut down this window
;     IF(!plotwindowsid eq 10)then begin
;     Wdelete,!plotwindowsid
;     !plotwindowsid=0
;     Endif
     widget_control, event.top, /destroy
     end


  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------
  neps=0  ;;0 output on screen; 1 for .eps file.
  rms=2   ;?张画不同的图的不同载入方式
  pbar=100  ;画colorbar分多少块来画。
  levs=256  ;假彩色的颜色数量
  contr=1.5 ;画线的粗细度
  contrs=1.5
 device,decomposed=0
 ;;Set this keyword to 0 to cause color values to be interpreted as indices into a color lookup table.
 ;;Set this keyword to 1 to cause color values to be interpreted as 24-bit color specifications.
  loadct,39
  !p.thick=contr
  !x.thick=contr
  !y.thick=contr
  !p.charsize=contrs
  !p.charthick=contr
  !p.background=255
  !p.color=0
  ; posi1=[0.25,0.06,0.85,0.08]
   posi1=[0.9,0.3,0.92,0.9]
   posi2=[0.13,0.3,0.73,0.9]  ;;左下角与右上角两点。

   if rms eq 1 then begin
   loadct,39
     ccol=indgen(levs)/double(levs)*(256)
     ;ccol(0:55)=25    ;用同一种颜色来显示大于（小于）某个值的区域
     ;ccol(200:255)=230 ;用同一种颜色来显示大于（小于）某个值的区域
   endif else begin
     ccol=255-indgen(levs)/double(levs)*(256-25)
     ;ccol(0:25)=255  ;用同一种颜色来显示大于（小于）某个值的区域
   endelse

  zd=fltarr(2*pmax-1,2*nmax-1)

  IF(strcmp(name[0],'Omega_sim(k)',12)) THEN BEGIN
  Omega=complexarr(pmax,nmax)
  Omegacy=complexarr(pmax,nmax)
  Omegaft=complexarr(pmax,nmax)
      Phi_kd=complexarr(2*pmax-1,2*nmax-1,ntmax/output_step+1)
      openr,lun,'Phi_k.txt',/get_lun
      readf,lun,Phi_kd
      free_lun,lun
      Phi_kftd=complexarr(2*pmax-1,2*nmax-1,ntmax/output_step+1)
      openr,lun,'Phi_kft.txt',/get_lun
      readf,lun,Phi_kftd
      free_lun,lun
      Phi_kcyd=complexarr(2*pmax-1,2*nmax-1,ntmax/output_step+1)
      openr,lun,'Phi_kcy.txt',/get_lun
      readf,lun,Phi_kcyd
      free_lun,lun
  Omega(*,*)=complex(0,0)
  Omegacy(*,*)=complex(0,0)
  Omegaft(*,*)=complex(0,0)

      tmax=max(Time)
      xx1=tmax*0.5
      xx2=tmax-tmax/40.0
      xxx1=round(xx1/(output_step*tstep))
      xxx2=round(xx2/(output_step*tstep))
      lxx1=xx1-tmax/40.0
      rxx1=xx1+tmax/40.0
      lxx2=xx2-tmax/40.0
      rxx2=xx2+tmax/40.0
      lxxx1=round(lxx1/(output_step*tstep))
      rxxx1=round(rxx1/(output_step*tstep))
      lxxx2=round(lxx2/(output_step*tstep))
      rxxx2=round(rxx2/(output_step*tstep))
   For p=0,pmax-1 Do Begin
    For n=0,nmax-1 Do Begin
     For i=1,ntmax/output_step Do Begin
     If((phi_kd(p+pmax-1,n+nmax-1,i)+phi_kd(p+pmax-1,n+nmax-1,i-1)) eq complex(0,0))then begin
     Omega(p,n)=Omega(p,n)+complex(0,1)*(phi_kd(p+pmax-1,n+nmax-1,i)-phi_kd(p+pmax-1,n+nmax-1,i-1))/(output_step*tstep)
     Endif Else Begin
     Omega(p,n)=Omega(p,n)+2*complex(0,1)*(phi_kd(p+pmax-1,n+nmax-1,i)-phi_kd(p+pmax-1,n+nmax-1,i-1))/(output_step*tstep*(phi_kd(p+pmax-1,n+nmax-1,i)+phi_kd(p+pmax-1,n+nmax-1,i-1)))
     Endelse
     Endfor
     Omega(p,n)=Omega(p,n)/(ntmax/output_step)

     For i=1,ntmax/output_step Do Begin
     If((phi_kcyd(p+pmax-1,n+nmax-1,i)+phi_kcyd(p+pmax-1,n+nmax-1,i-1)) eq complex(0,0))then begin
     Omegacy(p,n)=Omegacy(p,n)+complex(0,1)*(phi_kcyd(p+pmax-1,n+nmax-1,i)-phi_kcyd(p+pmax-1,n+nmax-1,i-1))/(output_step*tstep)
     Endif Else Begin
     Omegacy(p,n)=Omegacy(p,n)+2*complex(0,1)*(phi_kcyd(p+pmax-1,n+nmax-1,i)-phi_kcyd(p+pmax-1,n+nmax-1,i-1))/(output_step*tstep*(phi_kcyd(p+pmax-1,n+nmax-1,i)+phi_kcyd(p+pmax-1,n+nmax-1,i-1)))
     Endelse
     Endfor
     Omegacy(p,n)=Omegacy(p,n)/(ntmax/output_step)

     For i=1,ntmax/output_step Do Begin
     If((phi_kftd(p+pmax-1,n+nmax-1,i)+phi_kftd(p+pmax-1,n+nmax-1,i-1)) eq complex(0,0))then begin
     Omegaft(p,n)=Omegaft(p,n)+complex(0,1)*(phi_kftd(p+pmax-1,n+nmax-1,i)-phi_kftd(p+pmax-1,n+nmax-1,i-1))/(output_step*tstep)
     Endif Else Begin
     Omegaft(p,n)=Omegaft(p,n)+2*complex(0,1)*(phi_kftd(p+pmax-1,n+nmax-1,i)-phi_kftd(p+pmax-1,n+nmax-1,i-1))/(output_step*tstep*(phi_kftd(p+pmax-1,n+nmax-1,i)+phi_kftd(p+pmax-1,n+nmax-1,i-1)))
     Endelse
     Endfor
     Omegaft(p,n)=Omegaft(p,n)/(ntmax/output_step)

      y1=total(alog(abs(Phi_kd(p+pmax-1,n+nmax-1,lxxx1:rxxx1))))/(rxxx1-lxxx1+1)
      y2=total(alog(abs(Phi_kd(p+pmax-1,n+nmax-1,lxxx2:rxxx2))))/(rxxx2-lxxx2+1)
      y1cy=total(alog(abs(Phi_kcyd(p+pmax-1,n+nmax-1,lxxx1:rxxx1))))/(rxxx1-lxxx1+1)
      y2cy=total(alog(abs(Phi_kcyd(p+pmax-1,n+nmax-1,lxxx2:rxxx2))))/(rxxx2-lxxx2+1)
      y1ft=total(alog(abs(Phi_kftd(p+pmax-1,n+nmax-1,lxxx1:rxxx1))))/(rxxx1-lxxx1+1)
      y2ft=total(alog(abs(Phi_kftd(p+pmax-1,n+nmax-1,lxxx2:rxxx2))))/(rxxx2-lxxx2+1)

      Omega(p,n)=complex(REAL_PART(Omega(p,n)),(y2-y1)/(xx2-xx1))
      Omegacy(p,n)=complex(REAL_PART(Omegacy(p,n)),(y2cy-y1cy)/(xx2-xx1))
      Omegaft(p,n)=complex(REAL_PART(Omegaft(p,n)),(y2ft-y1ft)/(xx2-xx1))
      Omega(0,0)=complex(0,0)
      Omegacy(0,0)=complex(0,0)
      Omegaft(0,0)=complex(0,0)
    Endfor
   Endfor
   
;;----------------------------------------------------------------------
Endif Else IF( strcmp(name[0],'Omega_matr(k)',13) ) THEN BEGIN
  tempReconst=1  ;!!!!!!!!!!!!!!!!
  zdh=fltarr(pmax,nmax)
  Omega_matOrg=complexarr(2*pmax-1,nmax,N_mu)
  Omega_matrFinal=complexarr(pmax,nmax)
  openr,lun,'Omega_matr.txt',/get_lun   ;;0 Gyro-kinetic
  readf,lun,Omega_matOrg
  n=eof(lun)
  if n ne 1 then print,'error with file load!!!!!!'
  free_lun,lun
  tempRe=100
  For p=pmax-1,2*pmax-2 Do Begin
   For n=0,nmax-1 Do Begin
    jtemp=-1
    For j=0,N_mu-1 Do Begin
          IF(abs(REAL_PART(Omega_matOrg(p,n,j))) LE tempRe)THEN BEGIN
          jtemp=j
          growthIm=imaginary(Omega_matOrg(p,n,j))
          Endif
    Endfor

     For j=0,N_mu-1 Do Begin
       IF(abs(REAL_PART(Omega_matOrg(p,n,j))) LE tempRe)THEN BEGIN
       IF(imaginary(Omega_matOrg(p,n,j)) GT growthIm)then begin
       jtemp=j
       growthIm=imaginary(Omega_matOrg(p,n,j))
       Endif
       Endif
     Endfor
     IF(jtemp EQ -1)then begin
     Omega_matrFinal(p-pmax+1,n)=complex(0,0)
     Endif else begin
     Omega_matrFinal(p-pmax+1,n)=Omega_matOrg(p,n,jtemp)
     Endelse
   Endfor
  Endfor

  Omega_matrcyLow=complexarr(pmax,nmax)
  Omega_matrcyHigh=complexarr(pmax,nmax)
  Omega_matOrgcy=complexarr(2*pmax-1,nmax,N_mu*(2*N_CY-1))
  Omega_matOrgcy(*,*,*)=complex(0,0)
  openr,lun,'Omega_matrcy.txt',/get_lun
  readf,lun,Omega_matOrgcy
  n=eof(lun)
  if n ne 1 then print,'error with file load!!!!!!'
  free_lun,lun
  For p=0,pmax-1 Do Begin
   For n=0,nmax-1 Do Begin
  tempRe=abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,0)))
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
      IF(tempRe LE abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j)))) THEN BEGIN
           tempRe= abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j)))
       Endif
    Endfor
    IF(N_CY GE 2)then begin
    tempRe=tempRe/Omega_star*tempReconst/(N_CY-1) ;!!!!!!!!!!!!!!!!!! 1. 2.
    ;tempRe=1
    Endif

    jtemp=-1
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
          IF(abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j))) LE tempRe)THEN BEGIN
          jtemp=j
          growthIm=imaginary(Omega_matOrgcy(p+pmax-1,n,j))
          Endif
    Endfor
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
          IF(abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j))) LE tempRe)THEN BEGIN
           IF(growthIm LT imaginary(Omega_matOrgcy(p+pmax-1,n,j))) THEN BEGIN
           growthIm= imaginary(Omega_matOrgcy(p+pmax-1,n,j))
           jtemp=j
           Endif
          Endif
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrcyLow(p,n)=complex(0,0)
    Endif else begin
    Omega_matrcyLow(p,n)=Omega_matOrgcy(p+pmax-1,n,jtemp)
    Endelse

    jtemp=-1
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j))) GE tempRe)) THEN BEGIN
          jtemp=j
          growthIm=imaginary(Omega_matOrgcy(p+pmax-1,n,j))
          Endif
    Endfor
    For j=0,N_mu*(2*N_CY-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgcy(p+pmax-1,n,j))) GE tempRe)) THEN BEGIN
            IF(imaginary(Omega_matOrgcy(p+pmax-1,n,j)) GT growthIm) THEN BEGIN
            jtemp=j
            growthIm=imaginary(Omega_matOrgcy(p+pmax-1,n,j))
            Endif
          Endif
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrcyHigh(p,n)=complex(0,0)
    Endif else begin
    Omega_matrcyHigh(p,n)=Omega_matOrgcy(p+pmax-1,n,jtemp)
    Endelse
  Endfor
 Endfor

  Omega_matrftLow=complexarr(pmax,nmax)
  Omega_matrftHigh=complexarr(pmax,nmax)
  Omega_matrftHigh2=complexarr(pmax,nmax)
  Omega_matrftHigh3=complexarr(pmax,nmax)
  Omega_matOrgft=complexarr(2*pmax-1,nmax,N_mu*(2*N_FT-1))
  Omega_matOrgft(*,*,*)=complex(0,0)
  openr,lun,'Omega_matrft.txt',/get_lun   ;;0 Gyro-kinetic
  readf,lun,Omega_matOrgft
  n=eof(lun)
  if n ne 1 then print,'error with file load!!!!!!'
  free_lun,lun
 For p=0,pmax-1 Do Begin
  For n=0,nmax-1 Do Begin
  tempRe=abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,0)))
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
      IF(tempRe LE abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j)))) THEN BEGIN
           tempRe= abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j)))
       Endif
    Endfor
    IF(N_FT GE 2)then begin
    tempRe=tempRe/Omega_star*tempReconst/(N_FT-1);!!!!!!!!!!!! 3. 4. !!used as the limiter of the low freq and high freq
    Endif

    jtemp=-1
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF((abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j))) LE tempRe))THEN BEGIN
          jtemp=j
          growthIm=(imaginary(Omega_matOrgft(p+pmax-1,n,j)))
          EndIf
    Endfor
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF(abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j))) LE tempRe)THEN BEGIN
           IF(growthIm LT (imaginary(Omega_matOrgft(p+pmax-1,n,j)))) THEN BEGIN
           growthIm=imaginary(Omega_matOrgft(p+pmax-1,n,j))
           jtemp=j
           Endif
          Endif
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrftLow(p,n)=complex(0,0)
    Endif else begin
    Omega_matrftLow(p,n)=Omega_matOrgft(p+pmax-1,n,jtemp)
    Endelse

    jtemp=-1
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF(abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j))) GE tempRe) THEN BEGIN
            If(abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j))-REAL_PART(Omega_matrcyHigh(p,n))) LT freqPer*abs(REAL_PART(Omega_matrcyHigh(p,n)))) THEN BEGIN
             growthIm=imaginary(Omega_matOrgft(p+pmax-1,n,j))
             jtemp=j
            Endif
          EndIf
    Endfor
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF(abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j))) GE tempRe) THEN BEGIN
            If(abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j))-REAL_PART(Omega_matrcyHigh(p,n))) LT freqPer*abs(REAL_PART(Omega_matrcyHigh(p,n)))) THEN BEGIN
             IF((imaginary(Omega_matOrgft(p+pmax-1,n,j))) GT growthIm) THEN BEGIN
             growthIm=imaginary(Omega_matOrgft(p+pmax-1,n,j))
             jtemp=j
             Endif
            Endif
          EndIf
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrftHigh(p,n)=complex(0,0)
    Endif else begin
    Omega_matrftHigh(p,n)=Omega_matOrgft(p+pmax-1,n,jtemp)
    Endelse


    jtemp=-1
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF(abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j))) GE tempRe) THEN BEGIN
            If(abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j))-2*REAL_PART(Omega_matrcyHigh(p,n))) LT freqPer*abs(REAL_PART(Omega_matrcyHigh(p,n)))) THEN BEGIN
             growthIm=imaginary(Omega_matOrgft(p+pmax-1,n,j))
             jtemp=j
            Endif
          EndIf
    Endfor
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF(abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j))) GE tempRe) THEN BEGIN
            If(abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j))-2*REAL_PART(Omega_matrcyHigh(p,n))) LT freqPer*abs(REAL_PART(Omega_matrcyHigh(p,n)))) THEN BEGIN
             IF((imaginary(Omega_matOrgft(p+pmax-1,n,j))) GT growthIm) THEN BEGIN
             growthIm=imaginary(Omega_matOrgft(p+pmax-1,n,j))
             jtemp=j
             Endif
            Endif
          EndIf
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrftHigh2(p,n)=complex(0,0)
    Endif else begin
    Omega_matrftHigh2(p,n)=Omega_matOrgft(p+pmax-1,n,jtemp)
    Endelse

    jtemp=-1
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF(abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j))) GE tempRe) THEN BEGIN
            If(abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j))-3*REAL_PART(Omega_matrcyHigh(p,n))) LT freqPer*abs(REAL_PART(Omega_matrcyHigh(p,n)))) THEN BEGIN
             growthIm=imaginary(Omega_matOrgft(p+pmax-1,n,j))
             jtemp=j
            Endif
          EndIf
    Endfor
    For j=0,N_mu*(2*N_FT-1)-1 Do Begin
          IF(abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j))) GE tempRe) THEN BEGIN
            If(abs(REAL_PART(Omega_matOrgft(p+pmax-1,n,j))-3*REAL_PART(Omega_matrcyHigh(p,n))) LT freqPer*abs(REAL_PART(Omega_matrcyHigh(p,n)))) THEN BEGIN
             IF((imaginary(Omega_matOrgft(p+pmax-1,n,j))) GT growthIm) THEN BEGIN
             growthIm=imaginary(Omega_matOrgft(p+pmax-1,n,j))
             jtemp=j
             Endif
            Endif
          EndIf
    Endfor
    IF(jtemp EQ -1)then begin
    Omega_matrftHigh3(p,n)=complex(0,0)
    Endif else begin
    Omega_matrftHigh3(p,n)=Omega_matOrgft(p+pmax-1,n,jtemp)
    Endelse

  Endfor
 Endfor

   xbar=indgen(pbar)/float(pbar)
   ybar=indgen(pbar)/float(pbar)
   bars=dblarr(pbar,pbar)
;~~~color bar generation
  Omegacy=complex(0,0)
  Omegaft=complex(0,0)
  IF(!Fignum2 eq 1) Then Begin
   If(!Omegap eq 0)then begin
   zdh=REAL_PART(Omega_matrFinal)
   name[0]='Re('+cgGreek('omega')+') (GK)'
   Endif Else IF(!Omegap eq 1) THEN BEGIN
   zdh=imaginary(Omega_matrFinal)
   name[0]='Im('+cgGreek('omega')+') (GK)'
   Endif
  Endif Else IF(!Fignum2 eq 2) THEN BEGIN
   If(!Omegap eq 0)then begin
   zdh=REAL_PART(Omega_matrcyLow)
   name[0]='Re('+cgGreek('omega')+') (CKinCH) low freq'
   Endif Else IF(!Omegap eq 1) THEN BEGIN
   zdh=imaginary(Omega_matrcyLow)
   name[0]='Im('+cgGreek('omega')+') (CKinCH) low freq'
   Endif
  Endif Else IF(!Fignum2 eq 3) THEN BEGIN
   If(!Omegap eq 0)then begin
   zdh=REAL_PART(Omega_matrcyHigh)
   name[0]='Re('+cgGreek('omega')+') (CKinCH) high freq'
   Endif Else IF(!Omegap eq 1) THEN BEGIN
   zdh=imaginary(Omega_matrcyHigh)
   name[0]='Im('+cgGreek('omega')+') (CKinCH) high freq'
   Endif
  Endif Else IF(!Fignum2 eq 4) THEN BEGIN
   If(!Omegap eq 0)then begin
   zdh=REAL_PART(Omega_matrftLow)
   name[0]='Re('+cgGreek('omega')+') (CKinFH) low freq' 
   Endif Else IF(!Omegap eq 1) THEN BEGIN
   zdh=imaginary(Omega_matrftLow)
   name[0]='Im('+cgGreek('omega')+') (CKinFH) low freq'  
   Endif
  Endif Else IF(!Fignum2 eq 5) THEN BEGIN
   If(!Omegap eq 0)then begin
   zdh=REAL_PART(Omega_matrftHigh)
   name[0]='Re('+cgGreek('omega')+') (CKinFH) high freq'  
   Endif Else IF(!Omegap eq 1) THEN BEGIN
   zdh=imaginary(Omega_matrftHigh)
   name[0]='Im('+cgGreek('omega')+') (CKinFH) high freq' 
   Endif
  Endif Else IF(!Fignum2 eq 6) THEN BEGIN
   If(!Omegap eq 0)then begin
   zdh=REAL_PART(Omega_matrftHigh2)
   name[0]='Re('+cgGreek('omega')+') (CKinFH) high freq 2' 
   Endif Else IF(!Omegap eq 1) THEN BEGIN
   zdh=imaginary(Omega_matrftHigh2)
   name[0]='Im('+cgGreek('omega')+') (CKinFH) high freq 2'
   Endif
  Endif Else IF(!Fignum2 eq 7) THEN BEGIN
   If(!Omegap eq 0)then begin
   zdh=REAL_PART(Omega_matrftHigh3)
   name[0]='Re('+cgGreek('omega')+') (CKinFH) high freq 3'  
   Endif Else IF(!Omegap eq 1) THEN BEGIN
   zdh=imaginary(Omega_matrftHigh3)
   name[0]='Im('+cgGreek('omega')+') (CKinFH) high freq 3' 
   Endif
  Endif
   zmax=max(zdh)+(max(zdh)-min(zdh))/10
   zmin=min(zdh)-(max(zdh)-min(zdh))/10
   ;zdh=zdh+0.005
   If((zmax eq 0) and (zmin eq 0))then begin
   zmax=1
   zmin=-1
   endif
   vb1=zmax
   vb2=zmin
   for i=0,pbar-1 do begin
       for j=0,pbar-1 do begin
           bars(i,j)=j/double(pbar)*(vb1-vb2)+vb2
       endfor
   endfor
   ybar=ybar*(vb1-vb2)+vb2
   contour,bars,xbar,ybar,$
   title=' ',$
   xticks=2,$                ;;让x轴上不标字
   xtickname=[' ',' ',' '],$ ;;让x轴上不标字
   /fill,$
   nlevels=levs,$
   c_color=ccol,$
   position=posi1,$
   zstyle=1,xstyle=1,ystyle=1,$  ;;值1让图像充满坐标轴
   yrange=[zmin,zmax]

   contour,zdh(*,*),kxhalf,kyhalf,$
   title=name[0],xtitle='kx',ytitle='ky',$
   /fill,nlevels=levs,$
   c_color=ccol,$
   position=posi2,$
   zstyle=1,ystyle=1,$
   xstyle=1,$
   zrange=[zmin,zmax],$
   /noerase

   goto, Matrixtheend
   


;;----------------------------------------------------------------------   

Endif Else IF( strcmp(name[0],'|ne|',4) ) THEN BEGIN
;      zd=fltarr(2*pmax-1,nmax)
;      nsqr_k=fltarr(2*pmax-1,nmax,ntmax/output_step+1)
;      openr,lun,'nsqr_k.txt',/get_lun
;      readf,lun,nsqr_k
;      free_lun,lun
;      
;      zd=total(nsqr_k(*,*,round(sattime/(output_step*tstep)+1):round(ntmax/output_step)),3);the average energy of steady state
;      zd=zd/(ntmax/output_step-sattime/(output_step*tstep))
  
      zd=fltarr(2*pmax-1,nmax)
      
      nesat=fltarr(2*pmax-1,nmax)
      jump=1
      IF(ntmax/output_step GE 1000)then begin
      jump=10
      Endif Else IF(ntmax/output_step GE 10000)then begin
      jump=100
      Endif Else IF(ntmax/output_step GE 100000)then begin
      jump=1000
      Endif Else IF(ntmax/output_step GE 1000000)then begin
      jump=10000
      Endif

      nsqr_k=complexarr(2*pmax-1,2*nmax-1,ntmax/(output_step*jump)+1)
      nsqrtmp=complexarr(2*pmax-1,2*nmax-1)

      openr,lun,'ne_k.txt',/get_lun
      readf,lun,nsqrtmp
      POINT_LUN, -lun, pos
      Print,'Pos=',pos
      For i=0L,ntmax/output_step,jump Do begin
      POINT_LUN, lun, pos*i
      readf,lun,nsqrtmp
      nsqr_k(*,*,i/jump)=nsqrtmp
      Endfor
      free_lun,lun

For j=round(sattime/(output_step*jump*tstep)),round(ntmax/(output_step*jump)) Do Begin
;nesat(*,*)=nesat(*,*)+abs(nsqr_k(*,nmax-1:2*nmax-2,j))
nesat(*,*)=nesat(*,*)+sqrt(conj(nsqr_k(*,nmax-1:2*nmax-2,j))*nsqr_k(*,nmax-1:2*nmax-2,j))
Endfor
nesat(*,*)=nesat(*,*)*(1.0/(ntmax/(output_step*jump)-sattime/(output_step*jump*tstep)))
   
      zd=nesat
  
            
   xbar=indgen(pbar)/float(pbar)
   ybar=indgen(pbar)/float(pbar)
   bars=dblarr(pbar,pbar)
      
   zmax=max(zd)+(max(zd)-min(zd))/10
   zmin=min(zd)-(max(zd)-min(zd))/10
   If((zmax eq 0) and (zmin eq 0))then begin
   zmax=1
   zmin=-1
   endif
   vb1=zmax
   vb2=zmin
   for i=0,pbar-1 do begin
       for j=0,pbar-1 do begin
           bars(i,j)=j/double(pbar)*(vb1-vb2)+vb2
       endfor
   endfor
   
   ybarp=indgen(pbar)/float(pbar)
   ybarp=ybarp*(vb1-vb2)+vb2


   contour,bars,xbar,ybarp,$
   title=' ',$
   xticks=2,$                ;;让x轴上不标字
   xtickname=[' ',' ',' '],$ ;;让x轴上不标字
   /fill,$
   nlevels=levs,$
   c_color=ccol,$
   position=posi1,$
   zstyle=1,xstyle=1,ystyle=1,$  ;;值1让图像充满坐标轴
   yrange=[zmin,zmax]

   contour,zd(*,*),kx,kyhalf,$
   title=name[0],xtitle='kx',ytitle='ky',$
   /fill,nlevels=levs,$
   c_color=ccol,$
   position=posi2,$
   zstyle=1,ystyle=1,$
   xstyle=1,$
   zrange=[zmin,zmax],$
   /noerase


Endif

;;----------------------------------------------------------------------
IF( strcmp(name[0],'|Phi|^2',7) ) THEN BEGIN
  tempReconst=1  ;!!!!!!!!!!!!!!!!
  
  zdh=fltarr(2*pmax-1,2*nmax-1)

   xbar=indgen(pbar)/float(pbar)
   ybar=indgen(pbar)/float(pbar)
   bars=dblarr(pbar,pbar)
;~~~color bar generation
  
 ;;--------------------------------------------------------- 

;    Print,"input jump:"
;    read,jump
;    Print,"input imax_deno"
;    read,imax_deno
    
    jump=100
    imax_deno=1
    
    
   ibegin=0LL
   IF(Stopnt EQ 0)then begin
   Nsat=round(ntmax-begintime/tstep)/output_step     
   ibegin=round(begintime/tstep/output_step)
   ENDIF   
   IF(Stopnt GT 0)then begin
   Nsat=(ntmax-Stopnt)/output_step  
   ibegin=0  
   ENDIF
       

    imax=Nsat/imax_deno
    Time=findgen(imax/jump)*(jump*tstep*output_step) + Stopnt*tstep + begintime    
    Phi_ktmp=complexarr(2*pmax-1,2*nmax-1)  
    z=complexarr(2*pmax-1,2*nmax-1) 
    zphi=fltarr(2*pmax-1,2*nmax-1,imax/jump) 
    Phi_k=complexarr(2*pmax-1,2*nmax-1,imax/jump)  
    x=FINDGEN(2*pmax-1)
    y=FINDGEN(2*pmax-1)
    x=(x-pmax+1)/(kxmax/(pmax-1))/(2*pmax-2)*2*!PI
    y=(y-nmax+1)/(kymax/(nmax-1))/(2*nmax-2)*2*!PI
    
      pos=0LL
      openr,lun,'Phi_kOMG.txt',/get_lun
      readf,lun,Phi_ktmp
      POINT_LUN, -lun, pos
      Print,'Pos=',pos
      For i=0LL,imax-1,jump Do begin
      POINT_LUN, lun, pos*(i+ibegin) 
      readf,lun,Phi_ktmp
      Phi_k(*,*,i/jump)=Phi_ktmp
      Endfor
      free_lun,lun       
      
   For i=0LL,imax-1,jump Do Begin
   z=Phi_k(*,*,i/jump)
   z=shift(z,-pmax+1,-nmax+1)
   z = FFT(z,/INVERSE)
   If(log eq 0)then begin
   zphi(*,*,i/jump) = (REAL_PART(z))
   Endif Else Begin
     For p=0,2*pmax-2 Do Begin
     For n=0,2*nmax-2 Do Begin
     IF(REAL_PART(z(p,n)) GT 0)then begin
     zphi(*,*,i/jump) = Alog10(abs(REAL_PART(z(p,n))))
     ENDIF
     IF(REAL_PART(z(p,n)) LT 0)then begin
     zphi(*,*,i/jump) = -Alog10(abs(REAL_PART(z(p,n))))
     ENDIF
     Endfor 
     Endfor 
   Endelse
   
   ;zphi(*,*,i/jump)=zphi(*,*,i/jump)-total(total(zphi(*,*,i/jump),1),1)/(2*pmax-1)/(2*nmax-1)
   Endfor 
        
   print,'average=',total(total(total(zphi(*,*,*),1),1),1)/(2*pmax-1)/(2*nmax-1)/(imax/jump)        
   ;zphi=zphi-total(total(total(zphi(*,*,*),1),1),1)/(2*pmax-1)/(2*nmax-1)/(imax/jump)
   
    zdh=zphi(*,*,0)
        
;   top=max(zphi)
;   bottom=min(zphi)
;   zmax=min([abs(top),abs(bottom)])*1
;   zmin=-min([abs(top),abs(bottom)])*1
   
   zmax=20
   zmin=-20
   
   print,'zmax=',max(zphi)
   print,'zmin=',min(zphi)

   
   If((zmax eq 0) and (zmin eq 0))then begin
   zmax=1
   zmin=-1
   endif
   vb1=zmax
   vb2=zmin
   for i=0,pbar-1 do begin
       for j=0,pbar-1 do begin
           bars(i,j)=j/double(pbar)*(vb1-vb2)+vb2
       endfor
   endfor
   ybar=ybar*(vb1-vb2)+vb2
   contour,bars,xbar,ybar,$
   title='',$
   xticks=2,$                ;;让x轴上不标字
   xtickname=[' ',' ',' '],$ ;;让x轴上不标字
   /fill,$
   nlevels=levs,$
   c_color=ccol,$
   position=posi1,$
   zstyle=1,xstyle=1,ystyle=1,$  ;;值1让图像充满坐标轴
   yrange=[zmin,zmax]

   contour,zdh(*,*),x,y,$
   title='',$
   xtitle='x',ytitle='y',$
   /fill,nlevels=levs,$
   c_color=ccol,$
   position=posi2,$
   zstyle=1,ystyle=1,$
   xstyle=1,$
   /noerase,$
   zrange=[zmin,zmax]
;;------------------------------------------------     
;Test the 2 dimentional Fourier Transform
;
; WINDOW,2, XSIZE=600, YSIZE=500   
;   zmax=max(Phi_k(*,*,imax/jump/10))
;   zmin=min(Phi_k(*,*,imax/jump/10))
;   If((zmax eq 0) and (zmin eq 0))then begin
;   zmax=1
;   zmin=-1
;   endif
;   vb1=zmax
;   vb2=zmin
;   for i=0,pbar-1 do begin
;       for j=0,pbar-1 do begin
;           bars(i,j)=j/double(pbar)*(vb1-vb2)+vb2
;       endfor
;   endfor
;   ybar=ybar*(vb1-vb2)+vb2
;   
;   contour,bars,xbar,ybar,$
;   title=' ',$
;   xticks=2,$                ;;让x轴上不标字
;   xtickname=[' ',' ',' '],$ ;;让x轴上不标字
;   /fill,$
;   nlevels=levs,$
;   ;c_color=ccol,$
;   position=posi1,$
;   zstyle=1,xstyle=1,ystyle=1,$  ;;值1让图像充满坐标轴
;   yrange=[zmin,zmax] 
;   
;    contour,Phi_k(*,*,imax/jump/10),kx,ky,$
;   title='1',$
;   xtitle='kx',ytitle='ky',$
;   /fill,nlevels=levs,$
;   ;c_color=ccol,$
;   position=posi2,$
;   zstyle=1,ystyle=1,$
;   xstyle=1,$
;   /noerase,$
;   zrange=[zmin,zmax]
;;------------------------------------------------     
; WINDOW,3, XSIZE=600, YSIZE=500  
;   z3=zdh(*,*)
;   z3 = FFT(z3)    
;   z3=shift(z3,pmax-1,nmax-1)
;
;  
;   zmax=max(z3)
;   zmin=min(z3)
;   If((zmax eq 0) and (zmin eq 0))then begin
;   zmax=1
;   zmin=-1
;   endif
;   vb1=zmax
;   vb2=zmin
;   for i=0,pbar-1 do begin
;       for j=0,pbar-1 do begin
;           bars(i,j)=j/double(pbar)*(vb1-vb2)+vb2
;       endfor
;   endfor
;   ybar=ybar*(vb1-vb2)+vb2
;   
;   contour,bars,xbar,ybar,$
;   title=' ',$
;   xticks=2,$                ;;让x轴上不标字
;   xtickname=[' ',' ',' '],$ ;;让x轴上不标字
;   /fill,$
;   nlevels=levs,$
;   ;c_color=ccol,$
;   position=posi1,$
;   zstyle=1,xstyle=1,ystyle=1,$  ;;值1让图像充满坐标轴
;   yrange=[zmin,zmax] 
;   
;    contour,z3,kx,ky,$
;   title='3',$
;   xtitle='kx',ytitle='ky',$
;   /fill,nlevels=levs,$
;   ;c_color=ccol,$
;   position=posi2,$
;   zstyle=1,ystyle=1,$
;   xstyle=1,$
;   /noerase,$
;   zrange=[zmin,zmax] 
;   
;;------------------------------------------------   
; WINDOW,4, XSIZE=600, YSIZE=500  
;   z4=z3
;   z4=shift(z4,-pmax+1,-nmax+1)
;   z4 = FFT(z4,/INVERSE)    
;   
;   zmax=max(z4)
;   zmin=min(z4)
;   If((zmax eq 0) and (zmin eq 0))then begin
;   zmax=1
;   zmin=-1
;   endif
;   vb1=zmax
;   vb2=zmin
;   for i=0,pbar-1 do begin
;       for j=0,pbar-1 do begin
;           bars(i,j)=j/double(pbar)*(vb1-vb2)+vb2
;       endfor
;   endfor
;   ybar=ybar*(vb1-vb2)+vb2
;   
;   contour,bars,xbar,ybar,$
;   title=' ',$
;   xticks=2,$                ;;让x轴上不标字
;   xtickname=[' ',' ',' '],$ ;;让x轴上不标字
;   /fill,$
;   nlevels=levs,$
;   ;c_color=ccol,$
;   position=posi1,$
;   zstyle=1,xstyle=1,ystyle=1,$  ;;值1让图像充满坐标轴
;   yrange=[zmin,zmax] 
;   
;    contour,z4,kx,ky,$
;   title='4',$
;   xtitle='x',ytitle='y',$
;   /fill,nlevels=levs,$
;   ;c_color=ccol,$
;   position=posi2,$
;   zstyle=1,ystyle=1,$
;   xstyle=1,$
;   /noerase,$
;   zrange=[zmin,zmax] 
 ;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 
 IF(!velocity EQ 0)then Begin
 ;------------------------------------------------     
  For i=0LL,imax-1,jump Do Begin  ;;the main loop of the movie
   zdh=zphi(*,*,i/jump)
   contour,zdh(*,*),x,y,$
   ;title='phi:  '+strtrim(Time(i/jump),2)+'a/cs',$
   xtitle='x',ytitle='y',$
   /fill,nlevels=levs,$
   c_color=ccol,$
   position=posi2,$
   zstyle=1,ystyle=1,$
   xstyle=1,$
   zrange=[zmin,zmax],$
   /noerase,$
   /OVERPLOT
   
   xyouts,-10,17,'phi:  '+strtrim(Time(i/jump),2)+'   a/cs'
   wait,0.05
  Endfor
 ;------------------------------------------------
  yline=total(total((zphi),1),1)/(2*pmax-1)/(2*nmax-1)
  WINDOW,2, XSIZE=700, YSIZE=700   
   plot,[0],[0],$
       /nodata,$
       xstyle=1,$
       ystyle=1,$
       yminor=0,$
       xrange=[min(Time),max(Time)],$
       yrange=[min(yline)*1,max(yline)*1.02],$       
       xtitle='Time',$
        title="Total "+cgGreek('phi')+' (CK)'$
        ,psym=-5 ,symsize=1.5;,/xlog;,thich=8.0
        
    oplot,Time,yline,thick=0.1    
 ;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
 Endif else if(!velocity GT 0)then Begin
   print,'!velocity=',!velocity
   i=!velocity*jump
   zdh=zphi(*,*,i/jump)
   
   contour,zdh(*,*),x,y,$
   ;title='phi:  '+strtrim(Time(i/jump),2)+'a/cs',$
   xtitle='x',ytitle='y',$
   /fill,nlevels=levs,$
   c_color=ccol,$
   position=posi2,$
   zstyle=1,ystyle=1,$
   xstyle=1,$
   zrange=[zmin,zmax],$
   /noerase,$
   /OVERPLOT   
   xyouts,-10,23,'phi:  '+strtrim(Time(i/jump),2)+'   a/cs'
   
   deltax=1/(kxmax/(pmax-1))/(2*pmax-2)*2*!PI   
   deltay=1/(kymax/(nmax-1))/(2*nmax-2)*2*!PI
   const=1.8
      
      
   For p=0,2*pmax-3 Do begin
     For n=0,2*nmax-3 Do begin
     oplot,[x(p),x(p)+(zdh(p,n+1)-zdh(p,n))/deltay/const],$
     [y(n),y(n)+(-zdh(p+1,n)+zdh(p,n))/deltax/const],$
     thick=1.0
     
     !p.color=color_value(ncolor*8/16)  ;;white 
     oplot,[x(p)+(zdh(p,n+1)-zdh(p,n))/deltay/const,x(p)+(zdh(p,n+1)-zdh(p,n))/deltay/const],$
     [y(n)+(-zdh(p+1,n)+zdh(p,n))/deltax/const,y(n)+(-zdh(p+1,n)+zdh(p,n))/deltax/const],$
     thick=1,psym=4 ,symsize=0.2
     
      !p.color=color_value(ncolor+1)  ;black  
     Endfor
   Endfor

;--------------------------------------------
; plot velocity vector
   !P.BACKGROUND=color_value(ncolor/2)
    WINDOW,2, XSIZE=500, YSIZE=500  
   ; !P.COLOR=color_value(ncolor/2)
    
     v_eddy=fltarr(2*pmax-1,2*nmax-1)  
   For p=0,2*pmax-2 Do begin
     For n=0,2*nmax-2 Do begin
       if(n eq 2*nmax-2)then begin ;;cyclic boundary condition
       n_p=0
       endif else begin
       n_p=n+1
       endelse
       if(p eq 2*pmax-2)then begin ;;cyclic boundary condition
       p_p=0
       endif else begin
       p_p=p+1
       endelse   
            
       IF(zdh(p,n) GE 0)then begin    
       v_eddy(p,n)=sqrt(((zdh(p,n_p)-zdh(p,n))/deltay)^2+((-zdh(p_p,n)+zdh(p,n))/deltax)^2)
       Endif Else begin
       v_eddy(p,n)=-sqrt(((zdh(p,n_p)-zdh(p,n))/deltay)^2+((-zdh(p_p,n)+zdh(p,n))/deltax)^2)
       Endelse
     Endfor
   Endfor
   

   zmin=-max([abs(min(v_eddy)),abs(max(v_eddy))])
   zmax=max([abs(min(v_eddy)),abs(max(v_eddy))])
;   zmin=-10
;   zmax=10
   print,'vmax=',zmin
   print,'vmin=',zmax
   
   If((zmax eq 0) and (zmin eq 0))then begin
   zmax=1
   zmin=-1
   endif
   vb1=zmax
   vb2=zmin
   for i=0,pbar-1 do begin
       for j=0,pbar-1 do begin
           bars(i,j)=j/double(pbar)*(vb1-vb2)+vb2
       endfor
   endfor
   ybar=ybar*(vb1-vb2)+vb2
   contour,bars,xbar,ybar,$
   title='',$
   xticks=2,$                ;;让x轴上不标字
   xtickname=[' ',' ',' '],$ ;;让x轴上不标字
   /fill,$
   nlevels=levs,$
   c_color=ccol,$
   position=posi1,$
   zstyle=1,xstyle=1,ystyle=1,$  ;;值1让图像充满坐标轴
   yrange=[zmin,zmax]
         
     
   contour,v_eddy(*,*),x,y,$
   xtitle='x',ytitle='y',$
   /fill,nlevels=levs,$
   c_color=ccol,$
   position=posi2,$
   zstyle=1,ystyle=1,$
   xstyle=1,$
   zrange=[zmin,zmax],$
   /noerase;,$
   ;/OVERPLOT   
   
   i=!velocity*jump
   xyouts,-10,23,'|v| :  '+strtrim(Time(i/jump),2)+'   a/cs'
   
 Endif
   
Endif else begin
;;--------------------------------------------------------------------------------------

   if neps eq 0 then begin
   ;   window,3,xsize=700,ysize=700;,/pixmap
   endif else begin
      snapname=name[0]+string(rms,format='(I1)')+'bin'+'.eps'
      set_plot,'PS'  ;;Change the IDL graphics device to PostScript（为了 输出为eps图像）
      device,color=1,/encapsul,file=snapname,xsize=20,ysize=20,bits_per_pixel=8
      ;;provides device-dependent control over the current graphics device (as set by the SET_PLOT routine).
     ;;/encapsul :Set this keyword to create an encapsulated PostScript file, suitable for importing into another document (e.g., a LaTeX or FrameMaker document).
     ;;pixel像素的意思
   endelse

   xbar=indgen(pbar)/float(pbar)
   ybar=indgen(pbar)/float(pbar)
   bars=dblarr(pbar,pbar)

;~~~color bar generation
  IF(!Fignum eq 1) Then Begin
   If(!Omegap eq 0)then begin
   zd=REAL_PART(Omega)
   name[0]='Re('+cgGreek('omega')+') (GK) sim' 
   Endif Else IF(!Omegap eq 1) THEN BEGIN
   zd=imaginary(Omega)
   name[0]='Im('+cgGreek('omega')+') (GK) sim'
   Endif
  Endif Else IF(!Fignum eq 2) THEN BEGIN
   If(!Omegap eq 0)then begin
   zd=REAL_PART(Omegacy)
   name[0]='Re('+cgGreek('omega')+') (CKinCH) sim'  
   Endif Else IF(!Omegap eq 1) THEN BEGIN
   zd=imaginary(Omegacy)
   name[0]='Im('+cgGreek('omega')+') (CKinCH) sim'  
   Endif
  Endif Else IF(!Fignum eq 3) THEN BEGIN
   If(!Omegap eq 0)then begin
   zd=REAL_PART(Omegaft)
   name[0]='Re('+cgGreek('omega')+') (CKinFH) sim' 
   Endif Else IF(!Omegap eq 1) THEN BEGIN
   zd=imaginary(Omegaft)
   name[0]='Im('+cgGreek('omega')+') (CKinFH) sim'  
   Endif
  Endif
 CKinCHzmax=max(zd)+(max(zd)-min(zd))/10
   zmin=min(zd)-(max(zd)-min(zd))/10
   If((zmax eq 0) and (zmin eq 0))then begin
   zmax=1
   zmin=-1
   endif
   vb1=zmax
   vb2=zmin
   for i=0,pbar-1 do begin
       for j=0,pbar-1 do begin
           bars(i,j)=j/double(pbar)*(vb1-vb2)+vb2
       endfor
   endfor
   
   ybarp=indgen(pbar)/float(pbar)
   ybarp=ybarp*(vb1-vb2)+vb2


   contour,bars,xbar,ybarp,$
   title=' ',$
   xticks=2,$                ;;让x轴上不标字
   xtickname=[' ',' ',' '],$ ;;让x轴上不标字
   /fill,$
   nlevels=levs,$
   c_color=ccol,$
   position=posi1,$
   zstyle=1,xstyle=1,ystyle=1,$  ;;值1让图像充满坐标轴
   yrange=[zmin,zmax]

   contour,zd(*,*),kxhalf,kyhalf,$
   title=name[0],xtitle='kx',ytitle='ky',$
   /fill,nlevels=levs,$
   c_color=ccol,$
   position=posi2,$
   zstyle=1,ystyle=1,$
   xstyle=1,$
   zrange=[zmin,zmax],$
   /noerase

   Matrixtheend:

   if neps ne 0 then begin
      device,/close_file
      set_plot,!device_org
   End

help,loadct

;; open window
;   number_plot=number_plot+1
;   ;plotidm=plotidm+1
;   !p.background=color_value(ncolor/2)
;   window, 1, TITLE='movie', xsize=600,ysize=600
;; color table
;   set_viewport,0.2,0.8,0.85,0.95
;   !p.thick=2
;   !noeras=0
;   !mtitle='color table'
;   xmin=-float(256)/2
;   xmax=float(256)/2
;   ymin=0.0
;   ymax=1.0
;   set_xy,xmin,xmax,ymin,ymax
;   !p.color=color_value(ncolor/2)
;   plot,[0.0],[1,1]
;   for i=0,256-1 do begin
;  xa=xmin+float(i)
;  x=[xa,xa]
;  y=[0.0,1.0]
;  !p.color=color_value(i)
;  oplot,x,y
;    endfor

Endelse
 return
end

;*******************************************************************************


pro turbulence3D_event,event,group_leader=group

  common startup,number_plot,fpath,ncolor,color_value,plotid

  widget_control,event.id,get_uvalue=choice

  case choice of

    0: begin            ;exit idl
      widget_control, event.top,/destroy
      exit
    end

    1: begin            ;hplt
      number_plot=number_plot+1
      Main_plots_name
    end

    2: begin            
      number_plot=number_plot+1
      More_details_name
    end
    
    3: begin            ;hplt
      number_plot=number_plot+1
      Colorful_plot_name
    end

    4: begin            ;movie for 2d data
      number_plot=number_plot+1
      Movie
    end

    5: begin            ;colormap

      number_plot=number_plot+1
      window, number_plot, TITLE='snap', xsize=600,ysize=600
      !p.thick=16

      xmin=-float(ncolor+2)/2
      xmax=float(ncolor+2)/2
      ymin=0.0
      ymax=1.0
      set_xy,xmin,xmax,ymin,ymax

      for i=0,ncolor+1 do begin
        xa=xmin+float(i)
        x=[xa,xa]
        y=[0.0,1.0]
        !p.color=color_value(i)
        oplot,x,y
      endfor
      !p.thick=1

    end

    5: begin            ;exit idl
      widget_control, event.top,/destroy
      exit
    end

  endcase
end

    ; plot program for 3D turbulence code  ;;!!Main process
common startup,number_plot,fpath,ncolor,color_value,plotid

; read color map
 openr,1,'color.dat'
 ncolor=1
 readf,1,ncolor
 red=intarr(ncolor+2)
 green=intarr(ncolor+2)
 blue=intarr(ncolor+2)
 readf,1,red,green,blue
 close,1

 color_value=indgen(ncolor+2)

nbyte=3
; print, 'input 1 for 8-bit color display, 3 for true color''
; read, nbyte

;true color
 if nbyte eq 3 then color_value=red+256l*(green+256l*blue)

; load 8-bits color table
 if nbyte eq 1 then tvlct,red,green,blue

; default setting
; set_plot,'X'    ;================= We must remove this line on my computer
 defsysv,'!device_org',!D.NAME
 c=1
 !linetype=0
 !p.thick=2
 !p.charsize=2
 !p.charthick=2


Number_plot=0
fpath='.'

plotid=Number_plot

pname=strarr(5)
pname=["Exit IDL","Main plots","More details","Colorful Plot","Movie","ColorMap"]
xmenu,pname,BASE=pbase,SPACE=10,TITLE='3D turbulence',xpad=20,ypad=20
widget_control,pbase,/realize
xmanager,"turbulence3D",pbase,group_leader=group,/NO_BLOCK

end


;***********************************************************************


