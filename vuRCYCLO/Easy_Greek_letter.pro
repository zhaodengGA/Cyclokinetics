pro Easy_Greek_letter

!PATH = Expand_Path('+C:\idlfiles\coyote\') + ';' + !PATH

;!P.Font=1

print,"!P.Font=",!P.Font
;greekLetter = '!4' + String("166B) + '!X'

plot,[1,2],[3,4], XTitle='This title contains '+ 'Re('+cgGreek('omega')+')'

  
;cgplot,[1,1],[2,3],/nodata,XTitle='$\Omega$$\exp\\lambda$',$

;cgplot,[1,1],[2,3],/nodata,XTitle='time (a/c'+'$\downs$'+')',$
;Charsize=2.0,xrange=[0,4],yrange=[0,4]
;cgplot,[1,2],[2,4],/overplot
;cgplot,[1,2],[3,4],/overplot;,output='eps'

;cgPlot,[1,2],[3,4], XTitle='Length ($\mu$M)', YTitle='Distance ($\Angstrom$$\up2$)', $
;       Output='embedsymbols_1.png', Aspect = 0.66
  
;cgPlot, [1,2],[3,4], XTitle='$\Omega$$\exp\\lambda$', Charsize=2.0
   
;cgplot,[1,2],[3,4],xtitle='Distance ($\Angstrom$$\up2$)',ytitle='Length ($\mu$M)';,Font=-1

;p=plot([1,2],[3,4],xtitle='($\AA$)')


;t = TEXT(0.5, 0.5, 'Hello', 'r')


end