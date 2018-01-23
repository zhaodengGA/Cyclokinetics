pro pdf_statistics,t,y,i1,i2,tstep,output_step

  common plot_variables,y_bin,pdf,nbin,d_y

  pdf = fltarr(nbin)
  y_bin = fltarr(nbin)

  ymin = min(y[i1/(output_step*tstep):i2/(output_step*tstep)])
  ymax = max(y[i1/(output_step*tstep):i2/(output_step*tstep)])

  d_y = (ymax-ymin)/(nbin-1)

  y_bin = ymin+findgen(nbin)*d_y

  for ibin=0,nbin-1 do begin
     pdf[ibin] = 0.0 
     for i=i1,i2-1 do begin

        dt = t[(i+1)/(output_step*tstep)]-t[i/(output_step*tstep)]
        f  = 0.5*(y[i/(output_step*tstep)+1]+y[i/(output_step*tstep)])
        
        f_min = y_bin[ibin]-0.5*d_y
        f_max = y_bin[ibin]+0.5*d_y 
        
        if (f ge f_min and f lt f_max) then begin
           pdf[ibin] = pdf[ibin]+dt
        endif
        
     endfor
  endfor

  pdf[*] = pdf[*]/total(pdf)

  return

end

   

   
