; docformat = 'rst'
;
;About how to write superscript and subscript see:
;<http://www.idlcoyote.com/cg_tips/embedsymbols.php>
; NAME:
;   cgPlot
;
; PURPOSE:
;   The purpose of cgPlot is to create a wrapper for the traditional IDL graphics
;   command, Plot. The primary purpose of this is to create plot commands that work
;   and look identically both on the display and in PostScript files.
;
;******************************************************************************************;
;                                                                                          ;
;  Copyright (c) 2010, by Fanning Software Consulting, Inc. All rights reserved.           ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;
;+
; The purpose of cgPlot is to create a wrapper for the traditional IDL graphics
; command, Plot. The primary purpose of this is to create plot commands that work
; and look identically both on the display and in PostScript files.
; 
; Program default colors will depend on the IDL graphics window. If no windows are currently
; open when the program is called, cgDisplay is used to create a window.
; 
; The program requires the `Coyote Library <http://www.idlcoyote.com/documents/programs.php>`
; to be installed on your machine.
;
; :Categories:
;    Graphics
;    
; :Params:
;    x: in, required, type=any
;         If X is provided without Y, a vector representing the dependent values to be 
;         plotted If both X and Y are provided, X is the independent parameter and 
;         Y is the dependent parameter to be plotted.
;    y: in, optional, type=any
;         A vector representing the dependent values to be plotted.
;       
; :Keywords:
;     addcmd: in, optional, type=boolean, default=0
;        Set this keyword to add the command to the resizeable graphics window cgWindow.
;     aspect: in, optional, type=float, default=none
;        Set this keyword to a floating point ratio that represents the aspect ratio 
;        (ysize/xsize) of the resulting plot. The plot position may change as a result
;        of setting this keyword. Note that `Aspect` cannot be used when plotting with
;        !P.MULTI.
;     axiscolor: in, optional, type=string/integer, default='opposite'
;        If this keyword is a string, the name of the axis color. 
;        Otherwise, the keyword is assumed to be a color index into the current color table.
;     axescolor: in, optional, type=string/integer
;        Provisions for bad spellers.
;     background: in, optional, type=string/integer, default='background'
;        If this keyword is a string, the name of the background color. 
;        Otherwise, the keyword is assumed to be a color index into the current color table.
;     charsize: in, optional, type=float, default=cgDefCharSize()
;        The character size for axes annotations. Uses cgDefCharSize to select default
;        character size, unless !P.Charsize is set, in which case !P.Charsize is always used.
;     color: in, optional, type=string/integer, default='black'
;        If this keyword is a string, the name of the data color. By default, 'black'.
;        Color names are those used with cgColor. Otherwise, the keyword is assumed 
;        to be a color index into the current color table.
;     font: in, optional, type=integer, default=!P.Font
;        The type of font desired for axis annotation.
;     isotropic: in, optional, type=boolean, default=0
;        Maintain the same scale on both axes.
;     label: in, optional, type=string
;        A label is similar to a plot title, but it is aligned to the left edge
;        of the plot and is written in hardware fonts. Use of the label keyword
;        will suppress the plot title.
;     layout: in, optional, type=intarr(3)
;        This keyword specifies a grid with a graphics window and determines where the
;        graphic should appear. The syntax of LAYOUT is three numbers: [ncolumns, nrows, location].
;        The grid is determined by the number of columns (ncolumns) by the number of 
;        rows (nrows). The location of the graphic is determined by the third number. The
;        grid numbering starts in the upper left (1) and goes sequentually by column and then
;        by row.
;     legends: in, optional, type=object
;        One or more cgLegendItem objects that are to be drawn on the plot.
;     nodata: in, optional, type=boolean, default=0
;        Set this keyword to draw axes, but no data.
;     noerase: in, optional, type=boolean, default=0
;        Set this keyword to draw the plot without erasing the display first.
;     oplots: in, optional, type=object
;        A single cgOverPlot object, or an array of cgOverPlot objects that will be
;        overplot on the axes set up by the original data. The user will be responsible
;        for destroying the objects. The cgPlot program will simply draw the objects.
;     outfilename: in, optional, type=string
;        If the `Output` keyword is set, the user will be asked to supply an output
;        filename, unless this keyword is set to a non-null string. In that case, the
;        value of this keyword will be used as the filename and there will be no dialog
;        presented to the user.
;     output: in, optional, type=string, default=""
;        Set this keyword to the type of output desired. Possible values are these::
;            
;            'PS'   - PostScript file
;            'EPS'  - Encapsulated PostScript file
;            'PDF'  - PDF file
;            'BMP'  - BMP raster file
;            'GIF'  - GIF raster file
;            'JPEG' - JPEG raster file
;            'PNG'  - PNG raster file
;            'TIFF' - TIFF raster file
;            
;        Or, you can simply set this keyword to the name of the output file, and the type of
;        file desired will be determined by the file extension. If you use this option, the
;        user will not be prompted to supply the name of the output file.  
;            
;        All raster file output is created through PostScript intermediate files (the
;        PostScript files will be deleted), so ImageMagick and Ghostview MUST be installed 
;        to produce anything other than PostScript output. (See cgPS2PDF and PS_END for 
;        details.) And also note that you should NOT use this keyword when doing multiple 
;        plots. The keyword is to be used as a convenient way to get PostScript or raster 
;        output for a single graphics command. Output parameters can be set with cgWindow_SetDefs.
;     overplot: in, optional, type=boolean, default=0
;        Set this keyword if you wish to overplot data on an already exisiting set of
;        axes. It is like calling the IDL OPLOT command.
;     position: in, optional, type=vector
;        The usual four-element position vector for the Plot comamnd. Only monitored and
;        possibly set if the `Aspect` keyword is used.
;     psym: in, optional, type=integer
;        Any normal IDL PSYM values, plus any value supported by the Coyote Library
;        routine cgSYMCAT. An integer between 0 and 46. This may also be set to the
;        "name" of a symbol, such as returned from Print, cgSymCat(/Names).
;     symcolor: in, optional, type=string/integer, default='black'
;        If this keyword is a string, the name of the symbol color. By default, 'black'.
;        Otherwise, the keyword is assumed to be a color index into the current color table.
;     symsize: in, optional, type=float, default=1.0
;        The symbol size.
;     title: in, optional, type=string
;         The title of the plot.
;     traditional: in, optional, type=boolean, default=0
;        If this keyword is set, the traditional color scheme of a black background for
;        graphics windows on the display is used and PostScript files always use a white background.
;     window: in, optional, type=boolean, default=0
;        Set this keyword to replace all the commands in a current cgWindow or to
;        create a new cgWindow for displaying this command.
;     xrange: in, optional
;         Set this keyword to a two-element vector setting the X axis range for the plot.
;         If this keyword is used, and the `XStyle` keyword is NOT used, then XSTYLE is set to 1.
;     xstyle: in, optional, type=integer
;         This keyword is a bit map that allows a variety of axis options, depending upon which bit
;         is set. Bits are set by adding the following values together when setting the value of
;         the keyword::
;            Value    Description
;              0      Allow axis autoscaling.
;              1      Turn axis autoscaling off, force exact axis range.
;              2      Extend axis range.
;              4      Suppress entire axis.
;              8      Suppress box style axis. Draw only main axis.
;         To suppress box axis style and force exact axis range, for example, set the keyword to 8+1=9::
;             cgPlot, cgDemoData(1), XRange=[15,78], XStyle=9
;     xtitle: in, optional, type=string
;         The X title of the plot.
;     yrange: in, optional
;         Set this keyword to a two-element vector setting the Y axis range for the plot.
;         If this keyword is used, and the `YStyle` keyword is NOT used, then YSTYLE is set to 1.
;     ystyle: in, optional, type=integer
;         This keyword is a bit map that allows a variety of axis options, depending upon which bit
;         is set. Bits are set by adding the following values together when setting the value of
;         the keyword::
;            Value    Description
;              0      Allow axis autoscaling.
;              1      Turn axis autoscaling off, force exact axis range.
;              2      Extend axis range.
;              4      Suppress entire axis.
;              8      Suppress box style axis. Draw only main axis.
;             16      Inhibt setting the Y axis minimum value to 0.
;         To suppress box axis style and force exact axis range, for example, set the keyword to 8+1=9::
;             cgPlot, cgDemoData(1), YRange=[15,28], YStyle=9
;     ytitle: in, optional, type=string
;         The Y title of the plot.
;     _ref_extra: in, optional, type=any
;        Any `IDL Plot keyword <http://www.exelisvis.com/docs/PLOT_Procedure.html>` 
;        not defined here is allowed in the program.
;
; :Examples:
;    Use as you would use the IDL PLOT command::
;       cgPlot, Findgen(11)
;       cgPlot, Findgen(11), Aspect=1.0
;       cgPlot, Findgen(11), Color='olive', AxisColor='red', Thick=2
;       cgPlot, Findgen(11), Color='blue', SymColor='red', PSym=-16
;       
; :Author:
;    FANNING SOFTWARE CONSULTING::
;        David W. Fanning 
;        1645 Sheely Drive
;        Fort Collins, CO 80526 USA
;        Phone: 970-221-0438
;        E-mail: david@idlcoyote.com
;        Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Written, 12 November 2010. DWF.
;        Added SYMCOLOR keyword, and allow all 46 symbols from cgSYMCAT. 15 November 2010. DWF.
;        Added NODATA keyword. 15 November 2010. DWF.
;        Now setting decomposition state by calling SetDecomposedState. 16 November 2010. DWF.
;        Final color table restoration skipped in Z-graphics buffer. 17 November 2010. DWF.
;        Fixed a problem with overplotting with symbols. 17 November 2010. DWF.
;        Background keyword now applies in PostScript file as well. 17 November 2010. DWF.
;        Many changes after BACKGROUND changes to get !P.MULTI working again! 18 November 2010. DWF.
;        Fixed a small problem with the OVERPLOT keyword. 18 Nov 2010. DWF.
;        Changes so that color inputs don't change type. 23 Nov 2010. DWF.
;        Added WINDOW keyword to allow graphic to be displayed in a resizable graphics window. 8 Dec 2010. DWF
;        Modifications to allow cgPlot to be drop-in replacement for old PLOT commands in 
;            indexed color mode. 24 Dec 2010. DWF.
;        Previous changes introduced problems with OVERPLOT that have now been fixed. 28 Dec 2010. DWF.
;        Set NOERASE keyword from !P.NoErase system variable when appropriate. 28 Dec 2010. DWF.
;        Additional problems with NOERASE discovered and solved. 29 Dec 2010. DWF.
;        In some cases, I was turning BYTE values to strings without converting to 
;            INTEGERS first. 30 Dec 2010. DWF.  
;         Selecting character size now with cgDefCharSize. 11 Jan 2011. DWF.   
;         Moved setting to decomposed color before color selection process to avoid PostScript
;             background problems when passed 24-bit color integers. 12 Jan 2011. DWF. 
;         Changed _EXTRA to _REF_EXTRA on procedure definition statement to be able to return
;             plot keywords such as XGET_TICKS. 13 Jan 2011. DWF.  
;         Added SYMSIZE keyword. 16 Jan 2011. DWF.
;         Fixed a problem in which I assumed the background color was a string. 18 Jan 2011. DWF.  
;         Added ADDCMD keyword. 26 Jan 2011. DWF.
;         Added LAYOUT keyword. 28 Jan 2011. DWF.
;         Made a modification that allows THICK and COLOR keywords apply to symbols, too. 24 Feb 2011. DWF.
;         Modified error handler to restore the entry decomposition state if there is an error. 17 March 2011. DWF
;         Somehow I had gotten independent and dependent data reversed in the code. Put right. 16 May 2011. DWF.
;         Allowed ASPECT (and /ISOTROPIC) to take into account input POSITION. 15 June 2011. Jeremy Bailin.
;         Updated the BACKGROUND color selection from lessons learned in 27 Oct 2011 cgContour 
;             corrections. 27 Oct 2011. DWF.
;         Added the ability to send the output directly to a file via the OUTPUT keyword. 9 Dec 2011, DWF.
;         PostScript, PDF, and Imagemagick parameters can now be tailored with cgWindow_SetDefs. 14 Dec 2011. DWF.
;         Modified to use cgDefaultColor for default color selection. 24 Dec 2011. DWF.
;         Over-zealous use of _STRICT_EXTRA when overplotting resulted in errors. Now use _EXTRA. 1 Jan 2012. DWF.
;         Changes to allow better default colors, based on changes to cgColor and cgDefaultColor. 1 Feb 2012. DWF.
;         Now allowing the user to draw in the "background" color, if the COLOR or AXISCOLOR is "BACKGROUND". 19 March 2012. DWF.
;         Scalar input parameters are changed to 1-element vectors to avoid annoying error messages from PLOT. 6 April 2012. DWF.
;         Added a LABEL keyword. 12 July 2012. DWF.
;         Yikes! Bad choice of variable names in LABEL work yesterday has severe consequences. Changed names. 13 July 2012. DWF.
;         Added OPLOTS keyword to allow cgOverplot objects. 18 July 2012. DWF.
;         Added the ability to specify a symbol name with the PSYM keyword. 19 Juyl 2012. DWF.
;         Added the ability to use escape characters in plot titles to specify cgSymbol symbols. 27 July 2012. DWF.
;         Fixed an interaction with the LABEL keyword that prevented a Title from appearing. 2 Oct 2012. DWF.
;         Modified the way default colors are selected when the background color is "white". 4 Dec 2012. DWF.
;         Still trying to accommodate users who incorrectly specify LONG integers while using INDEXED color. 26 Dec 2012. DWF.
;         Modified code that checks to see if COLOR and AXISCOLOR keywords are the same as BACKGROUND and changes them.
;              This precludes drawing in background color on non-white backgrounds. Now only change the
;              colors if it is possible to draw a background color. 12 Feb 2013. DWF.
;         Problem using symbol names (e.g., 'opencircle') in cgWindows is fixed. 10 May 2013. DWF.
;         Changed the meaning of ISOTROPIC to its true meaning of keeping the same scale on both axes. 21 June 2013. DWF.
;         Added XRANGE, XSTYLE, YRANGE, and YSTYLE keywords. This allows exact axis scaling if the XRANGE or YRANGE
;             keywords are used without setting the XSTYLE or YSTYLE keywords, which is more intuitive. 15 July 2013. DWF.
;         
; :Copyright:
;     Copyright (c) 2010-2013, Fanning Software Consulting, Inc.
;-
PRO cgPlot, x, y, $
    ADDCMD=addcmd, $
    ASPECT=aspect, $
    AXISCOLOR=saxiscolor, $
    AXESCOLOR=saxescolor, $
    BACKGROUND=sbackground, $
    CHARSIZE=charsize, $
    COLOR=scolor, $
    FONT=font, $
    ISOTROPIC=isotropic, $
    LABEL=label, $
    LAYOUT=layout, $
    LEGENDS=legends, $
    NODATA=nodata, $
    NOERASE=noerase, $
    OPLOTS=oplots, $
    OUTFILENAME=outfilename, $
    OUTPUT=output, $
    OVERPLOT=overplot, $
    POSITION=position, $
    PSYM=psymIn, $
    SYMCOLOR=ssymcolor, $
    SYMSIZE=symsize, $
    TITLE=title, $
    TRADITIONAL=traditional, $
    WINDOW=window, $
    XRANGE=xrange, $
    XSTYLE=xstyle, $
    XTITLE=xtitle, $
    YRANGE=yrange, $
    YSTYLE=ystyle, $
    YTITLE=ytitle, $
    _REF_EXTRA=extra
    
!PATH = Expand_Path('+C:\idlfiles\coyote\') + ';' + !PATH

    
    Compile_Opt idl2
    
;!PATH = Expand_Path('+C:\idlfiles\coyote\') + ';' + !PATH

    Catch, theError
   ; theError=0
    
    IF theError NE 0 THEN BEGIN
        Catch, /CANCEL
        void = Error_Message()
        IF N_Elements(thisMulti) NE 0 THEN !P.Multi = thisMulti
        IF N_Elements(currentState) NE 0 THEN SetDecomposedState, currentState
        IF (N_Elements(output) NE 0) THEN PS_END, /NOFIX
        RETURN
    ENDIF
    
    ; Set up PostScript device for working with colors.
    IF !D.Name EQ 'PS' THEN Device, COLOR=1, BITS_PER_PIXEL=8
    
    ; Check parameters.
    IF N_Params() EQ 0 THEN BEGIN
        Print, 'USE SYNTAX: cgPlot, x, y'
        RETURN
    ENDIF
    
    ; Pay attention to !P.Noerase in setting the NOERASE kewyord. This must be
    ; done BEFORE checking the LAYOUT properties.
    IF !P.NoErase NE 0 THEN noerase = !P.NoErase ELSE noerase = Keyword_Set(noerase)
    
    ; Do they want this plot in a resizeable graphics window?
    IF Keyword_Set(addcmd) THEN window = 1
    IF Keyword_Set(window) AND ((!D.Flags AND 256) NE 0) THEN BEGIN
    
        ; If you are using a layout, you can't ever erase.
        IF N_Elements(layout) NE 0 THEN noerase = 1
        
        ; Special treatment for overplotting or adding a command.
        IF Keyword_Set(overplot) OR Keyword_Set(addcmd) THEN BEGIN
            cgWindow, 'cgPlot', x, y, $
                ASPECT=aspect, $
                AXISCOLOR=saxiscolor, $
                AXESCOLOR=saxescolor, $
                BACKGROUND=sbackground, $
                CHARSIZE=charsize, $
                COLOR=scolor, $
                FONT=font, $
                ISOTROPIC=isotropic, $
                LABEL=label, $
                LAYOUT=layout, $
                LEGENDS=legends, $
                NODATA=nodata, $
                NOERASE=noerase, $
                OPLOTS=oplots, $
                OVERPLOT=overplot, $
                POSITION=position, $
                PSYM=psymIn, $
                SYMCOLOR=ssymcolor, $
                SYMSIZE=symsize, $
                TITLE=title, $
                TRADITIONAL=traditional, $
                XTITLE=xtitle, $
                XRANGE=xrange, $
                XSTYLE=xstyle, $
                YRANGE=yrange, $
                YSTYLE=ystyle, $
                YTITLE=ytitle, $
                ADDCMD=1, $
                _Extra=extra
             RETURN
       ENDIF
        
        ; Open a new window or replace the current commands, as required.
        currentWindow = cgQuery(/CURRENT, COUNT=wincnt)
        IF wincnt EQ 0 THEN replaceCmd = 0 ELSE replaceCmd=1
        cgWindow, 'cgPlot', x, y, $
            ASPECT=aspect, $
            AXISCOLOR=saxiscolor, $
            AXESCOLOR=saxescolor, $
            BACKGROUND=sbackground, $
            CHARSIZE=charsize, $
            COLOR=scolor, $
            FONT=font, $
            ISOTROPIC=isotropic, $
            LABEL=label, $
            LAYOUT=layout, $
            LEGENDS=legends, $
            NODATA=nodata, $
            NOERASE=noerase, $
            OPLOTS=oplots, $
            OVERPLOT=overplot, $
            POSITION=position, $
            PSYM=psymIn, $
            SYMCOLOR=ssymcolor, $
            SYMSIZE=symsize, $
            TITLE=title, $
            TRADITIONAL=traditional, $
            XTITLE=xtitle, $
            XRANGE=xrange, $
            XSTYLE=xstyle, $
            YRANGE=yrange, $
            YSTYLE=ystyle, $
            YTITLE=ytitle, $
            REPLACECMD=replaceCmd, $
            _Extra=extra
            
         RETURN
    ENDIF
    
    ; Sort out which is the dependent and which is independent data.
    CASE N_Params() OF
      
       1: BEGIN
       dep = x
       indep = Findgen(N_Elements(dep))
       ENDCASE
    
       2: BEGIN
       dep = y
       indep = x
       ENDCASE
    
    ENDCASE
    
    ; If either of these input vectors are scalars, make them vectors.
    IF N_Elements(dep) EQ 1 THEN dep = [dep]
    IF N_Elements(indep) EQ 1 THEN indep = [indep]
    
    
    ; Check to see if psymIn is a string. If so, covert it here.
    IF N_Elements(psymIn) NE 0 THEN BEGIN
        IF Size(psymIn, /TNAME) EQ 'STRING' THEN BEGIN
            names = cgSymCat(/Names)
            index = Where(STRUPCASE(StrCompress(names, /REMOVE_ALL)) EQ STRUPCASE(StrCompress(psymIN, /REMOVE_ALL)), count)
            IF count GT 0 THEN psym = index[0] ELSE Message, 'Cannot resolve the PSYM value: ' + psymIn
        ENDIF ELSE psym = psymIn
    ENDIF
    
    ; Are we doing some kind of output?
    IF (N_Elements(output) NE 0) && (output NE "") THEN BEGIN
    
       ; If the output string has a dot character, then this must be a
       ; filename, and we will determine the type of file from the filename extension.
       IF StrPos(output, '.') NE -1 THEN BEGIN
             root_name = cgRootName(output, DIRECTORY=theDir, EXTENSION=ext)
             IF theDir EQ "" THEN CD, CURRENT=theDir
             outfilename = output
             outputSelection = StrUpCase(ext)
       ENDIF
    
       IF N_Elements(outputSelection) EQ 0 THEN outputSelection = StrUpCase(output)
       typeOfOutput = ['PS','EPS','PDF','BMP','GIF','JPEG','JPG','PNG','TIFF', 'TIF']
       void = Where(typeOfOutput EQ outputSelection, count)
       IF count EQ 0 THEN Message, 'Cannot find ' + outputSelection + ' in allowed output types.'
       
       ; Set things up.
       CASE outputSelection OF
          'PS': BEGIN
              ext = '.ps'
              delete_ps = 0
              END    
          'EPS': BEGIN
              ext = '.eps'
              encapsulated = 1
              delete_ps = 0
              END
          'PDF': BEGIN
              ext = '.pdf'
              pdf_flag = 1
              delete_ps = 1
              END     
          'BMP': BEGIN
              ext = '.bmp'
              bmp_flag = 1
              delete_ps = 1
              END      
          'GIF': BEGIN
              ext = '.gif'
              gif_flag = 1
              delete_ps = 1
              END
          'JPEG': BEGIN
              ext = '.jpg'
              jpeg_flag = 1
              delete_ps = 1
              END      
          'JPG': BEGIN
              ext = '.jpg'
              jpeg_flag = 1
              delete_ps = 1
              END
          'PNG': BEGIN
              ext = '.png'
              png_flag = 1
              delete_ps = 1
              END      
          'TIFF': BEGIN
              ext = '.tif'
              tiff_flag = 1
              delete_ps = 1
              END
          'TIF': BEGIN
              ext = '.tif'
              tiff_flag = 1
              delete_ps = 1
              END    
       ENDCASE
       
       ; Do you need a filename?
       IF ( (N_Elements(outfilename) EQ 0) || (outfilename EQ "") ) THEN BEGIN 
            filename = 'cgplot' + ext
            outfilename = cgPickfile(FILE=filename, TITLE='Select Output File Name...', $
                FILTER=ext, /WRITE)
            IF outfilename EQ "" THEN RETURN
       ENDIF
       
       ; We need to know the root name of the file, because we have to make a PostScript
       ; file of the same name. At least we do if the type is not PS or EPS.
       IF (outputSelection NE 'PS') && (outputSelection NE 'EPS') THEN BEGIN
           root_name = cgRootName(outfilename, DIRECTORY=theDir)
           IF theDir EQ "" THEN CD, CURRENT=theDir
           ps_filename = Filepath(ROOT_DIR=theDir, root_name + '.ps')
       ENDIF ELSE ps_filename = outfilename
       
       ; Get the output default values.
       cgWindow_GetDefs, $
         PS_Charsize = ps_charsize, $          ; The PostScript character size.
         PS_FONT = ps_font, $                  ; Select the font for PostScript output.
         PS_Decomposed = ps_decomposed, $      ; Sets the PostScript color mode.
         PS_Delete = ps_delete, $              ; Delete PS file when making IM raster.
         PS_Metric = ps_metric, $              ; Select metric measurements in PostScript output.
         PS_Scale_factor = ps_scale_factor, $  ; Select the scale factor for PostScript output.
         PS_TT_Font = ps_tt_font               ; Select the true-type font to use for PostScript output.   
       
       ; Set up the PostScript device.
       PS_Start, $
          CHARSIZE=ps_charsize, $
          DECOMPOSED=ps_decomposed, $
          FILENAME=ps_filename, $
          FONT=ps_font , $
          ENCAPSULATED=encapsulated, $
          METRIC=ps_metric, $
          SCALE_FACTOR=ps_scale_factor, $
          TT_FONT=ps_tt_font, $
          QUIET=1
    
    
    ENDIF
   
    ; Get the current color table vectors.
    TVLCT, rr, gg, bb, /GET
    
    ; Going to do this in decomposed color, if possible.
    SetDecomposedState, 1, CURRENTSTATE=currentState
    
    ; If current state is "indexed color" and colors are represented as long integers then "fix" them.
    IF (currentState EQ 0) THEN BEGIN
      IF Size(sbackground, /TNAME) EQ 'LONG' THEN sbackground = Fix(sbackground)
      IF Size(saxiscolor, /TNAME) EQ 'LONG' THEN saxiscolor = Fix(saxiscolor)
      IF Size(saxescolor, /TNAME) EQ 'LONG' THEN saxescolor = Fix(saxescolor)
      IF Size(scolor, /TNAME) EQ 'LONG' THEN scolor = Fix(scolor)
      IF Size(ssymcolor, /TNAME) EQ 'LONG' THEN ssymcolor = Fix(ssymcolor)
    ENDIF
    
    ; Set up the layout, if necessary.
    IF N_Elements(layout) NE 0 THEN BEGIN
       thisMulti = !P.Multi
       totalPlots = layout[0]*layout[1]
       !P.Multi = [0,layout[0], layout[1], 0, 0]
       IF layout[2] EQ 1 THEN BEGIN
            noerase = 1
            !P.Multi[0] = 0
       ENDIF ELSE BEGIN
            !P.Multi[0] = totalPlots - layout[2] + 1
       ENDELSE
    ENDIF

    ; Check the color keywords.
    IF N_Elements(title) EQ 0 THEN title = ""
    IF N_Elements(xtitle) EQ 0 THEN xtitle = ""
    IF N_Elements(ytitle) EQ 0 THEN ytitle = ""
    IF N_Elements(xrange) NE 0 THEN BEGIN
       IF N_Elements(xstyle) EQ 0 THEN xstyle = 1 
    ENDIF
    IF N_Elements(yrange) NE 0 THEN BEGIN
        IF N_Elements(ystyle) EQ 0 THEN ystyle = 1
    ENDIF
    title = cgCheckForSymbols(title)
    xtitle = cgCheckForSymbols(xtitle)
    ytitle = cgCheckForSymbols(ytitle)
    IF (N_Elements(label) NE 0) && (label NE "") THEN title = ""
    traditional = Keyword_Set(traditional)
    background = cgDefaultColor(sbackground, /BACKGROUND, TRADITIONAL=traditional)
    IF Size(background, /TNAME) EQ 'STRING' && (StrUpCase(background[0]) EQ 'WHITE') THEN BEGIN
       IF (N_Elements(saxisColor) EQ 0) && (N_Elements(saxesColor) NE 0) THEN saxisColor = saxesColor
       axisColor = cgDefaultColor(saxisColor, DEFAULT='black', TRADITIONAL=traditional)
       color = cgDefaultColor(sColor, DEFAULT='black', TRADITIONAL=traditional)
    ENDIF ELSE BEGIN
       IF (N_Elements(saxisColor) EQ 0) && (N_Elements(saxesColor) NE 0) THEN saxisColor = saxesColor
       axisColor = cgDefaultColor(saxisColor, TRADITIONAL=traditional)
       color = cgDefaultColor(sColor, DEFAULT=axisColor, TRADITIONAL=traditional)
    ENDELSE

    ; If color is the same as background, do something. Since this precludes drawing the 
    ; background color (perhaps you want to "erase" something), I offer an exception. If the
    ; COLOR is "Background", I am going to assume you know what you are doing!
    IF ColorsAreIdentical(background, color) THEN BEGIN
        IF ((!D.Flags AND 256) NE 0) THEN BEGIN
           IF (!P.Multi[0] EQ 0) && (~Keyword_Set(overplot) && ~noerase) THEN cgErase, background
        ENDIF
        IF (Size(color, /TNAME) EQ 'STRING') THEN BEGIN
            IF (StrUpCase(color) NE 'BACKGROUND') THEN BEGIN
                IF ~noerase && ~Keyword_Set(overplot) THEN color = 'OPPOSITE'
            ENDIF
        ENDIF ELSE BEGIN
            IF ~noerase && ~Keyword_Set(overplot) THEN color = 'OPPOSITE'
        ENDELSE
    ENDIF
    IF ColorsAreIdentical(background, axiscolor) THEN BEGIN
        IF ((!D.Flags AND 256) NE 0) THEN BEGIN
           IF (!P.Multi[0] EQ 0) && (~Keyword_Set(overplot) && ~noerase) THEN cgErase, background
        ENDIF
        IF (Size(axiscolor, /TNAME) EQ 'STRING') THEN BEGIN
           IF (StrUpCase(axiscolor) NE 'BACKGROUND') THEN BEGIN
               IF ~noerase && ~Keyword_Set(overplot) THEN axiscolor = 'OPPOSITE'
           ENDIF
        ENDIF ELSE BEGIN
             IF ~noerase && ~Keyword_Set(overplot) THEN axiscolor = 'OPPOSITE'
        ENDELSE
    ENDIF
    symcolor = cgDefaultColor(ssymcolor, DEFAULT=color, TRADITIONAL=traditional)
    
    ; Character size has to be determined *after* the layout has been decided.
    IF N_Elements(font) EQ 0 THEN font = !P.Font
    IF N_Elements(charsize) EQ 0 THEN charsize = cgDefCharSize(FONT=font)
    
    ; Other keywords.
    IF N_Elements(symsize) EQ 0 THEN symsize = 1.0
    IF Keyword_Set(isotropic) THEN BEGIN
        yscale = Max(dep)-Min(dep)
        xscale = Max(indep)-Min(indep)
        aspect = Float(yscale)/xscale
        xstyle=1
        ystyle=1
    ENDIF
    IF N_Elements(psym) EQ 0 THEN psym = 0
    IF (N_Elements(aspect) NE 0) AND (Total(!P.MULTI) EQ 0) THEN BEGIN
    
        ; If position is set, then fit the plot into those bounds.
        IF (N_Elements(position) GT 0) THEN BEGIN
          trial_position = Aspect(aspect, margin=0.)
          trial_width = trial_position[2]-trial_position[0]
          trial_height = trial_position[3]-trial_position[1]
          pos_width = position[2]-position[0]
          pos_height = position[3]-position[1]

          ; Same logic as cgImage: try to fit image width, then if you can't get the right aspect
          ; ratio, fit the image height instead.
          fit_ratio = pos_width / trial_width
          IF trial_height * fit_ratio GT pos_height THEN $
             fit_ratio = pos_height / trial_height

          ; new width and height
          trial_width *= fit_ratio
          trial_height *= fit_ratio

          ; calculate position vector based on trial_width and trial_height
          position[0] += 0.5*(pos_width - trial_width)
          position[2] -= 0.5*(pos_width - trial_width)
          position[1] += 0.5*(pos_height - trial_height)
          position[3] -= 0.5*(pos_height - trial_height)
        ENDIF ELSE position=Aspect(aspect)   ; if position isn't set, just use output of Aspect
        
    ENDIF
    
    ; If you get here with no position defined, and no layout, and no !P.Multi and no nothing,
    ; then for God's sake, define a reasonable position in the window!
    IF (N_Elements(position) EQ 0) && (Total(!P.Position) EQ 0) && (N_Elements(layout) EQ 0) && (Total(!P.Multi) LE 0) THEN BEGIN
        position = [0.125, 0.125, 0.925, 0.9] 
    ENDIF
           
    ; Do you need a PostScript background color? Lot's of problems here!
    ; Basically, I MUST draw a plot to advance !P.MULTI. But, drawing a
    ; plot of any sort erases the background color. So, I have to draw a 
    ; plot, store the new system variables, then draw my background, etc.
    ; I have tried LOTS of options. This is the only one that worked.
    IF !D.Name EQ 'PS' THEN BEGIN
         IF ~noerase THEN BEGIN
       
           ; I only have to do this, if this is the first plot.
           IF !P.MULTI[0] EQ 0 THEN BEGIN
           
                IF Keyword_Set(overplot) NE 1 THEN BEGIN
                
                    ; Save the current system variables. Will need to restore later.
                    bangx = !X
                    bangy = !Y
                    bangp = !P
                    
                    ; Draw the plot that doesn't draw anything.
                    Plot, indep, dep, POSITION=position, CHARSIZE=charsize, /NODATA, $
                        FONT=font, XRANGE=xrange, XSTYLE=xstyle, YRANGE=yrange, YSTYLE=ystyle, $
                        _STRICT_EXTRA=extra  
                    
                    ; Save the "after plot" system variables. Will use later. 
                    afterx = !X
                    aftery = !Y
                    afterp = !P     
                    
                    ; Draw the background color and set the variables you will need later.
                    PS_Background, background
                    psnodraw = 1
                    tempNoErase = 1
                    
                    ; Restore the original system variables so that it is as if you didn't
                    ; draw the invisible plot.
                    !X = bangx
                    !Y = bangy
                    !P = bangp
                
                ENDIF
            ENDIF ELSE tempNoErase = noerase
        ENDIF ELSE tempNoErase = noerase
     ENDIF ELSE tempNoErase = noerase
 
    
    ; Load the drawing colors. If needed create a window first, so the drawing
    ; colors are correct for the window you want to draw into.
    IF ((!D.Flags AND 256) NE 0) && (!D.Window LT 0) THEN cgDisplay
    IF Size(axiscolor, /TNAME) EQ 'STRING' THEN axiscolor = cgColor(axiscolor)
    IF Size(color, /TNAME) EQ 'STRING' THEN color = cgColor(color)
    IF Size(background, /TNAME) EQ 'STRING' THEN background = cgColor(background)
    IF Size(symcolor, /TNAME) EQ 'STRING' THEN symcolor = cgColor(symcolor)
    
    ; Draw the plot.
    IF Keyword_Set(overplot) THEN BEGIN
       IF psym LE 0 THEN OPlot, indep, dep, COLOR=color, _EXTRA=extra
    ENDIF ELSE BEGIN
      Plot, indep, dep, BACKGROUND=background, COLOR=axiscolor, CHARSIZE=charsize, $
            POSITION=position, /NODATA, NOERASE=tempNoErase, FONT=font, TITLE=title, $
            XTITLE=xtitle, YTITLE=ytitle, XRANGE=xrange, YRANGE=yrange, $
            XSTYLE=xstyle, YSTYLE=ystyle, _STRICT_EXTRA=extra
        IF psym LE 0 THEN BEGIN
           IF ~Keyword_Set(nodata) THEN OPlot, indep, dep, COLOR=color, _EXTRA=extra  
        ENDIF  
    ENDELSE
    IF Abs(psym) GT 0 THEN BEGIN
        asymbol = cgSymCat(Abs(psym), COLOR=symcolor, _Extra=extra)
        IF ~Keyword_Set(nodata) THEN OPlot, indep, dep, COLOR=symcolor, $
            PSYM=asymbol, SYMSIZE=symsize, _EXTRA=extra
    ENDIF 
    
    ; Do you have overplot objects to plot?
    IF N_Elements(oplots) NE 0 THEN BEGIN
        FOR j=0,N_Elements(oplots)-1 DO BEGIN
           thisObject = oplots[j]
           IF Obj_Valid(thisObject) THEN thisObject -> Draw
        ENDFOR
    ENDIF
    
    ; Do you have legend objects to draw?
    IF N_Elements(legends) NE 0 THEN BEGIN
        FOR j=0,N_Elements(legends)-1 DO BEGIN
           thisObject = legends[j]
           IF Obj_Valid(thisObject) THEN thisObject -> Draw
        ENDFOR
    ENDIF

    ; Need a label on the plot?
    IF N_Elements(label) NE 0 THEN BEGIN
        xx = !X.Window[0]
        yy = !Y.Window[1] + 0.015
        labelfont = (!D.Name EQ 'PS') ? 1 : 0
        cgText, xx, yy, /NORMAL, label, FONT=labelfont, COLOR=axiscolor
    ENDIF
         
    ; If this is the first plot in PS, then we have to make it appear that we have
    ; drawn a plot, even though we haven't.
    IF N_Elements(psnodraw) EQ 1 THEN BEGIN
        !X = afterX
        !Y = afterY
        !P = afterP
    ENDIF
    
    ; Are we producing output? If so, we need to clean up here.
    IF (N_Elements(output) NE 0) && (output NE "") THEN BEGIN
    
       ; Get the output default values.
       cgWindow_GetDefs, $
           IM_Density = im_density, $                      ; Sets the density parameter on ImageMagick convert command.
           IM_Options = im_options, $                      ; Sets extra ImageMagick options on the ImageMagick convert command.
           IM_Resize = im_resize, $                        ; Sets the resize parameter on ImageMagick convert command.
           IM_Transparent = im_transparent, $              ; Sets the "alpha" keyword on ImageMagick convert command.
           IM_Width = im_width, $                          ; Sets the width of raster file output created with ImageMagick.
           PDF_Unix_Convert_Cmd = pdf_unix_convert_cmd, $  ; Command to convert PS to PDF.
           PDF_Path = pdf_path                             ; The path to the Ghostscript conversion command.
    
        ; Close the PostScript file and create whatever output is needed.
        PS_END, DELETE_PS=delete_ps, $
             ALLOW_TRANSPARENT=im_transparent, $
             BMP=bmp_flag, $
             DENSITY=im_density, $
             GIF=gif_flag, $
             GS_PATH=pdf_path, $
             IM_OPTIONS=im_options, $
             JPEG=jpeg_flag, $
             PDF=pdf_flag, $
             PNG=png_flag, $
             RESIZE=im_resize, $
             TIFF=tiff_flag, $
             UNIX_CONVERT_CMD=pdf_unix_convert_cmd, $
             WIDTH=im_width

         basename = File_Basename(outfilename)
         dirname = File_Dirname(outfilename)
         IF dirname EQ "." THEN CD, CURRENT=dirname
         Print, 'Output File: ' + Filepath(ROOT_DIR=dirname, basename)
    ENDIF
    
    ; Restore the decomposed color state if you can.
    SetDecomposedState, currentState
    
    ; Restore the color table. Can't do this for the Z-buffer or
    ; the snap shot will be incorrect.
    IF (!D.Name NE 'Z') THEN TVLCT, rr, gg, bb
    
    ; Clean up if you are using a layout.
    IF N_Elements(layout) NE 0 THEN !P.Multi = thisMulti
    
END
    
;;==============================================================================
;;+
; NAME:
;  cgGREEK
;
; PURPOSE:
;
;   This function provides a device-independent way to ask for a Greek letter as
;   a string that can be included, for example, in a plot title. It uses the Greek
;   simplex font (!4) when used with Hershey fonts, and the Symbol font (!9) when
;   used with PostScript or True-Type fonts. Selects the type of Greek character to 
;   return based on value of !P.FONT. Updated now to return the UNICODE values for 
;   Greek characters for those fonts that support them (Macintosh?).
;   
; AUTHOR:
;
;   FANNING SOFTWARE CONSULTING
;   David Fanning, Ph.D.
;   1645 Sheely Drive
;   Fort Collins, CO 80526 USA
;   Phone: 970-221-0438
;   E-mail: david@idlcoyote.com
;   Coyote's Guide to IDL Programming: http://www.idlcoyote.com/
;
; CATEGORY:
;
;   Graphics, Utilities
;
; CALLING SEQUENCE:
;
;   greekString = cgGreek(greekLetter)
;
; RETURN VALUE:
;
;   greekString    A string that represents the Greek letter.
;
; ARGUMENTS:
;
;  greekLetter:    The name of the Greek letter desired. A string. Default: 'alpha'.
;                  Valid string names are the 24 characters of the Greek alphabet.
;                     alpha        nu
;                     beta         xi
;                     gamma        omicron
;                     delta        pi
;                     epsilon      rho
;                     zeta         sigma
;                     eta          tau
;                     theta        upsilon
;                     iota         phi
;                     kappa        chi
;                     lambda       psi
;                     mu           omega
;                    
;                   Note that if the first letter of the name is capitalized, this is
;                   the equivalent of setting the CAPITAL keyword. 
;
; KEYWORDRS:
;
;  CAPTIAL:        If this keyword is set, the captial Greek letter is returned rather 
;                  than the lowercase Greek letter. An alternative way of capitalizing
;                  the letter is to make the first letter of the name an uppercase letter.
;                  
;  EXAMPLE:        If this keyword is set, the names of the Greek characters and their
;                  symbols are written out in the current graphics window.
;                  
;  PS:             Normally, the PostScript version of the greek letter is returned if
;                  the current device is PostScript and !P.Font is 0 or 1. But, the 
;                  PostScript version of the greek letter can be obtained at any time
;                  and in any device, by setting this keyword.
;                                    
;  UNICODE:        If this keyword is set, the function returns the Unicode for the Greek
;                  letter.
;
; EXAMPLE:
;
;  Lowercase PSI:
;  
;     IDL> Plot, findgen(11), XTitle='This title contains ' + $
;           cgGreek('psi') + ' as a Greek letter' 
;
;  Uppercase PSI:
;  
;     IDL> Plot, findgen(11), XTitle='This title contains ' + $
;           cgGreek('Psi') + ' as a Greek letter' 
; NOTES:
; 
;  See the following article for additional information: 
;  
;       http://www.idlcoyote.com/ps_tips/greeksym.html
;       
; RESTRICTIONS:
; 
;  For this program to work correctly on your graphics display, you should be using
;  Hershey fonts (!P.Font=-1). It will work correctly in PostScript with either 
;  hardware fonts (!P.Font=0) or True-Type fonts (!P.Font=1).
;  
; MODIFICATION HISTORY:
;
;  Written by: David W. Fanning, 9 January 2010.
;  An alternative way to get an uppercase letter is to make the first letter of
;     the Greek name uppercase. (Maarten Sneep's suggestion!) 11 Jan 2010. DWF
;  Had the wrong value for the PostScript version of Phi. 26 January 2010. DWF
;  Added UNICODE keyword and values for Greek characters. 11 June 2010. DWF.
;  Changed the branching from !D.NAME EQ 'PS' to !P.FONT NE -1. (This is actually
;      what the documentation says, and what I intended.) 13 Dec 2010. DWF.
;  I don't think the last change did quite want I wanted. More tweaking to make
;      this more responsive to being in a PostScript file. 31 July 2011. DWF.
;  Added PS keyword so the return value is the PostScript file. This is for
;      convenience only, as the return value will be the PostScript value if
;      the current graphics device is PS and FONT is not equal to -1. 30 Aug 2011. DWF.
;  Retired Greek and replaced by cgGreek. 23 Dec 2012. DWF.
;-
;******************************************************************************************;
;  Copyright (c) 2010, by Fanning Software Consulting, Inc.                                ;
;  All rights reserved.                                                                    ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
Forward_Function cgGreek

PRO cgGreek_Example, UNICODE=unicode, PS=ps

    Compile_Opt hidden
    
    ; The 24 Greek letters.
    letter = [  'alpha', 'beta', 'gamma', 'delta', 'epsilon',  'zeta', $
                'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu', $
                'nu', 'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau', $
                'upsilon', 'phi', 'chi', 'psi', 'omega' ]
    
    ; Output positions.
    x = [0.25, 0.6]
    y = Reverse((Indgen(12) + 1) * (1.0 / 13))
    
    ; Need PostScript output?
    IF Keyword_Set(ps) || (!D.Name EQ 'PS') THEN PS_Start
    
    ; Create a window.
    cgDisplay, 600, 500
    
    ; Output the letters.
    FOR j=0,11 DO BEGIN
        cgText, x[0], y[j], letter[j] + ': ' + $
            cgGreek(letter[j], UNICODE=unicode) + cgGreek(letter[j], /CAPITAL, UNICODE=unicode), $
            /NORMAL, CHARSIZE=1.5
        cgText, x[1], y[j], letter[j+12] + ': ' + $
            cgGreek(letter[j+12], UNICODE=unicode) + cgGreek(letter[j+12], /CAPITAL, UNICODE=unicode), $
            /NORMAL, CHARSIZE=1.5
    ENDFOR
    
    ; Clean up PostScript, if needed.
    IF Keyword_Set(ps) || (!D.Name EQ 'PS') THEN PS_End

    ; Restore the users window.
    IF N_Elements(thisWindow) NE 0 THEN BEGIN
       IF thisWindow GE 0 THEN WSet, thisWindow
    ENDIF
    
END ; --------------------------------------------------------------------------------------


FUNCTION cgGreek, letter, CAPITAL=capital, EXAMPLE=example, PS=PS, UNICODE=unicode

    Compile_Opt idl2
    
    ; Return to caller on error.
    ON_Error, 2
    
    ; Do you wish to see an example?
    IF Keyword_Set(example) THEN BEGIN
        cgGreek_Example, UNICODE=unicode, PS=ps
        RETURN, ""
    ENDIF

    ; A symbol name is required.
    IF N_Elements(letter) EQ 0 THEN BEGIN
       Print, 'Syntax: Greek("theLetter")
       RETURN, ""
    ENDIF
    
    ; Return quietly if the Greek letter is a null string.
    IF letter  EQ "" THEN RETURN, ""
    
     ; Set up PostScript device for working with Greek letters.
    IF !D.Name EQ 'PS' THEN Device, ISOLATIN1=1
    
    ; Check keywords.
    IF N_Elements(letter) EQ 0 THEN letter = 'alpha'
    capital = Keyword_Set(capital)
    
    ; If the first letter of the "letter" variable is uppercase
    ; the user wants an uppercase Greek letter.
    firstLetter = StrMid(letter, 0, 1)
    IF firstLetter EQ StrUpCase(firstLetter) THEN capital = 1
    
    IF Keyword_Set(unicode) THEN BEGIN
    
        CASE StrLowCase(letter) OF
            'alpha':   greekLetter = (capital) ? '!Z(0391)' : '!Z(03B1)'
            'beta':    greekLetter = (capital) ? '!Z(0392)' : '!Z(03B2)'
            'gamma':   greekLetter = (capital) ? '!Z(0393)' : '!Z(03B3)'
            'delta':   greekLetter = (capital) ? '!Z(0394)' : '!Z(03B4)'
            'epsilon': greekLetter = (capital) ? '!Z(0395)' : '!Z(03B5)'
            'zeta':    greekLetter = (capital) ? '!Z(0396)' : '!Z(03B6)'
            'eta':     greekLetter = (capital) ? '!Z(0397)' : '!Z(03B7)'
            'theta':   greekLetter = (capital) ? '!Z(0398)' : '!Z(03B8)'
            'iota':    greekLetter = (capital) ? '!Z(0399)' : '!Z(03B9)'
            'kappa':   greekLetter = (capital) ? '!Z(039A)' : '!Z(03BA)'
            'lambda':  greekLetter = (capital) ? '!Z(039B)' : '!Z(03BB)'
            'mu':      greekLetter = (capital) ? '!Z(039C)' : '!Z(03BC)'
            'nu':      greekLetter = (capital) ? '!Z(039D)' : '!Z(03BD)'
            'xi':      greekLetter = (capital) ? '!Z(039E)' : '!Z(03BE)'
            'omicron': greekLetter = (capital) ? '!Z(039F)' : '!Z(03BF)'
            'pi':      greekLetter = (capital) ? '!Z(03A0)' : '!Z(03C0)'
            'rho':     greekLetter = (capital) ? '!Z(03A1)' : '!Z(03C1)'
            'sigma':   greekLetter = (capital) ? '!Z(03A3)' : '!Z(03C3)'
            'tau':     greekLetter = (capital) ? '!Z(03A4)' : '!Z(03C4)'
            'upsilon': greekLetter = (capital) ? '!Z(03A5)' : '!Z(03C5)'
            'phi':     greekLetter = (capital) ? '!Z(03A6)' : '!Z(03C6)'
            'chi':     greekLetter = (capital) ? '!Z(03A7)' : '!Z(03C7)'
            'psi':     greekLetter = (capital) ? '!Z(03A8)' : '!Z(03C8)'
            'omega':   greekLetter = (capital) ? '!Z(03A9)' : '!Z(03C9)'
            ELSE: Message, 'The greek letter ' + letter + ' is unrecognized.'
       ENDCASE
           
       ; Return the UNICODE greek letter.
       RETURN, greekLetter
       
    ENDIF
    
    IF ((!D.Name EQ 'PS') && (!P.FONT NE -1)) || Keyword_Set(PS)  THEN BEGIN
    
        ; Make sure ISOLATIN1 encoding is turned on.
        IF (!D.Name EQ 'PS') THEN DEVICE, /ISOLATIN1
        
        CASE StrLowCase(letter) OF
            'alpha':   greekLetter = (capital) ? '!9' + String("101B) + '!X' : '!9' + String("141B) + '!X'
            'beta':    greekLetter = (capital) ? '!9' + String("102B) + '!X' : '!9' + String("142B) + '!X'
            'gamma':   greekLetter = (capital) ? '!9' + String("107B) + '!X' : '!9' + String("147B) + '!X'
            'delta':   greekLetter = (capital) ? '!9' + String("104B) + '!X' : '!9' + String("144B) + '!X'
            'epsilon': greekLetter = (capital) ? '!9' + String("105B) + '!X' : '!9' + String("145B) + '!X'
            'zeta':    greekLetter = (capital) ? '!9' + String("132B) + '!X' : '!9' + String("172B) + '!X'
            'eta':     greekLetter = (capital) ? '!9' + String("110B) + '!X' : '!9' + String("150B) + '!X'
            'theta':   greekLetter = (capital) ? '!9' + String("121B) + '!X' : '!9' + String("161B) + '!X'
            'iota':    greekLetter = (capital) ? '!9' + String("111B) + '!X' : '!9' + String("151B) + '!X'
            'kappa':   greekLetter = (capital) ? '!9' + String("113B) + '!X' : '!9' + String("153B) + '!X'
            'lambda':  greekLetter = (capital) ? '!9' + String("114B) + '!X' : '!9' + String("154B) + '!X'
            'mu':      greekLetter = (capital) ? '!9' + String("115B) + '!X' : '!9' + String("155B) + '!X'
            'nu':      greekLetter = (capital) ? '!9' + String("116B) + '!X' : '!9' + String("156B) + '!X'
            'xi':      greekLetter = (capital) ? '!9' + String("130B) + '!X' : '!9' + String("170B) + '!X'
            'omicron': greekLetter = (capital) ? '!9' + String("117B) + '!X' : '!9' + String("157B) + '!X'
            'pi':      greekLetter = (capital) ? '!9' + String("120B) + '!X' : '!9' + String("160B) + '!X'
            'rho':     greekLetter = (capital) ? '!9' + String("122B) + '!X' : '!9' + String("162B) + '!X'
            'sigma':   greekLetter = (capital) ? '!9' + String("123B) + '!X' : '!9' + String("163B) + '!X'
            'tau':     greekLetter = (capital) ? '!9' + String("124B) + '!X' : '!9' + String("164B) + '!X'
            'upsilon': greekLetter = (capital) ? '!9' + String("125B) + '!X' : '!9' + String("165B) + '!X'
            'phi':     greekLetter = (capital) ? '!9' + String("106B) + '!X' : '!9' + String("152B) + '!X'
            'chi':     greekLetter = (capital) ? '!9' + String("103B) + '!X' : '!9' + String("143B) + '!X'
            'psi':     greekLetter = (capital) ? '!9' + String("131B) + '!X' : '!9' + String("171B) + '!X'
            'omega':   greekLetter = (capital) ? '!9' + String("127B) + '!X' : '!9' + String("167B) + '!X'
            ELSE: Message, 'The greek letter ' + letter + ' is unrecognized.'
        ENDCASE
    
    ENDIF ELSE BEGIN
    
        CASE StrLowCase(letter) OF
            'alpha':   greekLetter = (capital) ? '!4' + String("101B) + '!X' : '!4' + String("141B) + '!X'
            'beta':    greekLetter = (capital) ? '!4' + String("102B) + '!X' : '!4' + String("142B) + '!X'
            'gamma':   greekLetter = (capital) ? '!4' + String("103B) + '!X' : '!4' + String("143B) + '!X'
            'delta':   greekLetter = (capital) ? '!4' + String("104B) + '!X' : '!4' + String("144B) + '!X'
            'epsilon': greekLetter = (capital) ? '!4' + String("105B) + '!X' : '!4' + String("145B) + '!X'
            'zeta':    greekLetter = (capital) ? '!4' + String("106B) + '!X' : '!4' + String("146B) + '!X'
            'eta':     greekLetter = (capital) ? '!4' + String("107B) + '!X' : '!4' + String("147B) + '!X'
            'theta':   greekLetter = (capital) ? '!4' + String("110B) + '!X' : '!4' + String("150B) + '!X'
            'iota':    greekLetter = (capital) ? '!4' + String("111B) + '!X' : '!4' + String("151B) + '!X'
            'kappa':   greekLetter = (capital) ? '!4' + String("112B) + '!X' : '!4' + String("152B) + '!X'
            'lambda':  greekLetter = (capital) ? '!4' + String("113B) + '!X' : '!4' + String("153B) + '!X'
            'mu':      greekLetter = (capital) ? '!4' + String("114B) + '!X' : '!4' + String("154B) + '!X'
            'nu':      greekLetter = (capital) ? '!4' + String("115B) + '!X' : '!4' + String("155B) + '!X'
            'xi':      greekLetter = (capital) ? '!4' + String("116B) + '!X' : '!4' + String("156B) + '!X'
            'omicron': greekLetter = (capital) ? '!4' + String("117B) + '!X' : '!4' + String("157B) + '!X'
            'pi':      greekLetter = (capital) ? '!4' + String("120B) + '!X' : '!4' + String("160B) + '!X'
            'rho':     greekLetter = (capital) ? '!4' + String("121B) + '!X' : '!4' + String("161B) + '!X'
            'sigma':   greekLetter = (capital) ? '!4' + String("122B) + '!X' : '!4' + String("162B) + '!X'
            'tau':     greekLetter = (capital) ? '!4' + String("123B) + '!X' : '!4' + String("163B) + '!X'
            'upsilon': greekLetter = (capital) ? '!4' + String("124B) + '!X' : '!4' + String("164B) + '!X'
            'phi':     greekLetter = (capital) ? '!4' + String("125B) + '!X' : '!4' + String("165B) + '!X'
            'chi':     greekLetter = (capital) ? '!4' + String("126B) + '!X' : '!4' + String("166B) + '!X'
            'psi':     greekLetter = (capital) ? '!4' + String("127B) + '!X' : '!4' + String("167B) + '!X'
            'omega':   greekLetter = (capital) ? '!4' + String("130B) + '!X' : '!4' + String("170B) + '!X'
            ELSE: Message, 'The greek letter ' + letter + ' is unrecognized.'      
        ENDCASE
        
      ENDELSE

    RETURN, greekLetter
    
END    

;;==============================================================================
; docformat = 'rst'
;
; NAME:
;   cgSymbol
;
; PURPOSE:
; This function provides a device-independent way to ask for commonly-used
; symbols (e.g., less than or equal to, Angstrom, degree symbol, etc.),
; including all 24 Greek symbols in upper and lower case, to display with text.
; 
;******************************************************************************************;
;                                                                                          ;
;  Copyright (c) 2011, by Fanning Software Consulting, Inc. All rights reserved.           ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;
;+
; This function provides a device-independent way to ask for commonly-used
; symbols (e.g., less than or equal to, Angstrom, degree symbol, subscripts, superscripts, etc.),
; including all 24 Greek symbols in upper and lower case, to display with text.
;
; For this program to work correctly on your graphics display, you should be using
; Hershey fonts (!P.Font=-1). It will work correctly in PostScript with either 
; hardware fonts (!P.Font=0) or True-Type fonts (!P.Font=1).
; 
; `Greek Symbols <http://www.idlcoyote.com/ps_tips/greeksym.html>` are created by
; calling the Coyote Library routine `cgGreek <http://www.idlcoyote.com/programs/cggreek.pro>' 
; from this program.
; 
; .. image:: cgsymbol.png 
; 
; Normally, rather than calling cgSymbol, the symbols are embedded in Coyote Graphics
; text that are used for axis annotation and so forth. See `Embedding Symbols in Coyote
; Graphics Output <http://www.idlcoyote.com/cg_tips/embedsymbols.php>`. Embedded subscripts
; and superscripts are implemented like this::
;     $\upTEXT$ : Raise the text in TEXT to superscript level and reduce character size by 62%.
;     $\expTEXT$ : Raise the text in TEXT to superscript level and reduce character size by 44%.
;     $\downTEXT$ : Lower the text in TEXT to subscript level and reduce character size by 62%.
;     $\subTEXT$ : Lower the text in TEXT to subscript level and reduce character size by 44%.
;
; :Categories:
;    Graphics
;    
; :Returns:
;    A string variable that represents the requested symbol and can be used
;    in a textual context.
;    
; :Examples:
;     To create a lowercase Greek psi symbol::
;  
;        IDL> cgPlot, findgen(11), XTitle='This title contains ' + $
;                 cgSymbol('psi') + ' as a Greek letter' 
;        IDL> cgPlot, findgen(11), XTitle='This title contains $\psi$ as a Greek letter'
;
;     To create an Angstrom sign::
;  
;        IDL> cgPlot, findgen(11), XTitle='This title contains (' + $
;                cgSymbol('Angstrom') +  ') an Angstrom sign.' 
;        IDL> cgPlot, findgen(11), XTitle='This title contains ($\Angstrom$) an Angstrom sign.'
;       
;     To create subscripts or superscripts::
;  
;        IDL> cgPlot, findgen(11), XTitle='E=mc$\up2$' 
;        IDL> cgPlot, findgen(11), XTitle='H$\sub2$O'
;        IDL> cgPlot, findgen(11), XTitle='H$\upSuper$MT $\Omega$$\subSubscript$', Charsize=2.0
;        
;     It is possible to use Greek characters as superscripts and subscripts. Do so by
;     prepending the Greek character with "\\" inside the normal superscript or subscript
;     notation. For example, to use lambda as an exponent to the Greek character Omega, you
;     can type this::
;
;        IDL> cgPlot, cgDemoData(1), XTitle='$\Omega$$\exp\\lambda$', Charsize=2.0
;
;     To use lambda as a subscript, type this:
;
;         IDL> cgPlot, cgDemoData(1), XTitle='$\Omega$$\sub\\lambda$', Charsize=2.0
;
; :Author:
;       FANNING SOFTWARE CONSULTING::
;           David W. Fanning 
;           1645 Sheely Drive
;           Fort Collins, CO 80526 USA
;           Phone: 970-221-0438
;           E-mail: david@idlcoyote.com
;           Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Written by: David W. Fanning, 2 September 2011. 
;        Added plus-minus symbol. 2 Nov 2011. DWF.
;        Added "up", "down", "exp" "sub" and "n" symbols for subscripting and superscripting. 9 Nov 2012. DWF.
;        Added "division" and "times" signs. 24 Nov 2012. DWF.
;        Updated UNICODE values to display capital letters correctly. 23 Dec 2012. DWF.
;
; :Copyright:
;     Copyright (c) 2011, Fanning Software Consulting, Inc.
;-
;+
; Displays the symbols and their names in a graphics window.
; 
; :Keywords:
;    ps: in, optional, type=boolean, default=0                                
;        Normally, the PostScript version of the symbol is returned automatically if
;        the current device is PostScript and !P.Font is 0 or 1. But, the PostScript version 
;        of the symbol can be obtained at any time and in any device, by setting this keyword.
;    unicode: in, optional, type=boolean, default=0                                
;        If this keyword is set, the function returns the Unicode value for the symbol.
; 
;-
PRO cgSymbol_Example, PS=ps, UNICODE=unicode

    Forward_Function cgSymbol

    Compile_Opt hidden
    
    ; The symbols..
    symbol = [  'alpha', 'beta', 'gamma', 'delta', 'epsilon',  'zeta', $
                'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu', $
                'nu', 'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau', $
                'upsilon', 'phi', 'chi', 'psi', 'omega', 'leq', 'geq', $
                'neq', 'deg', '+-', 'equiv', 'prime', 'angstrom', 'sun', $
                'varphi', 'infinity', 'copyright' ]
    
    ; Output positions.
    x = [0.15, 0.4, 0.65]
    y = Reverse((Indgen(13) + 1) * (1.0 / 14))
    
    ; Need PostScript output?
    IF Keyword_Set(ps) || (!D.Name EQ 'PS') THEN PS_Start
    
    ; Create a window.
    cgDisplay, 800, 500
    
    ; Output the symbols.
    FOR j=0,11 DO BEGIN
        cgText, x[0], y[j], symbol[j] + ': ' + $
            cgSymbol(symbol[j], UNICODE=unicode) + cgSymbol(symbol[j], /CAPITAL, UNICODE=unicode), $
            /NORMAL, CHARSIZE=1.5
        cgText, x[1], y[j], symbol[j+12] + ': ' + $
            cgSymbol(symbol[j+12], UNICODE=unicode) + cgSymbol(symbol[j+12], /CAPITAL, UNICODE=unicode), $
            /NORMAL, CHARSIZE=1.5
        IF j NE 12 THEN cgText, x[2], y[j], symbol[j+24] + ': ' + $
            cgSymbol(symbol[j+24], UNICODE=unicode), $
            /NORMAL, CHARSIZE=1.5
    ENDFOR
    
    ; Clean up PostScript, if needed.
    IF Keyword_Set(ps) || (!D.Name EQ 'PS') THEN PS_End
    
    ; Restore the users window.
    IF N_Elements(thisWindow) NE 0 THEN BEGIN
       IF thisWindow GE 0 THEN WSet, thisWindow
    ENDIF
    
END ; --------------------------------------------------------------------------------------


;+
; This function provides a device-independent way to ask for commonly-used
; symbols (e.g., less than or equal to, Angstrom, degree symbol, etc.),
; including all 24 Greek symbols in upper and lower case, to display with text.
;
; :Params:
;    symbol: in, required, type=string, default="alpha"
;       The name of the symbol desired as a string. Valid string names are the 24 
;       characters of the Greek alphabet, plus commonly used graphical symbols::
;       
;          alpha        nu        leq
;          beta         xi        geq
;          gamma        omicron   neg
;          delta        pi        deg
;          epsilon      rho       equiv
;          zeta         sigma     prime
;          eta          tau       angstrom
;          theta        upsilon   sun
;          iota         phi       varphi
;          kappa        chi       infinity
;          lambda       psi       copyright
;          mu           omega
;                    
;       Note that if the first letter of the name is capitalized, this is
;       the equivalent of setting the `Capital` keyword. Only Greek characters
;       use this method of selecting symbols.
;       
; :Keywords:
;    capital: in, optional, type=boolean, default=0
;        If this keyword is set, the captial Greek letter is returned rather than the lowercase 
;        Greek letter. An alternative way of capitalizing the letter is to make the first letter 
;        of the name an uppercase letter.
;    example: in, optional, type=boolean, default=0             
;        If this keyword is set, the names of symbols and the symbols themselves are written 
;        out in the current graphics window.
;    ps: in, optional, type=boolean, default=0                                
;        Normally, the PostScript version of the symbol is returned automatically if
;        the current device is PostScript and !P.Font is 0 or 1. But, the PostScript version 
;        of the symbol can be obtained at any time and in any device, by setting this keyword.
;    unicode: in, optional, type=boolean, default=0                                
;        If this keyword is set, the function returns the Unicode value for the symbol.
;          
;-
FUNCTION cgSymbol, symbol, CAPITAL=capital, EXAMPLE=example, PS=PS, UNICODE=unicode

    Compile_Opt idl2
    
    ; Return to caller on error.
    ON_Error, 2
    
    ; A common block to communicate with PS_Start.
    COMMON _$FSC_PS_START_, ps_struct
    
    ; Do you wish to see an example?
    IF Keyword_Set(example) THEN BEGIN
        cgSymbol_Example, UNICODE=unicode, PS=ps
        RETURN, ""
    ENDIF

    ; A symbol name is required.
    IF N_Elements(symbol) EQ 0 THEN BEGIN
       Print, 'Syntax: cgSymbol("theSymbol")
       RETURN, ""
    ENDIF
    
    ; Return quietly if the symbol is a null string.
    IF symbol EQ "" THEN RETURN, ""
    
    ; Set up PostScript device for working with Greek letters.
    IF !D.Name EQ 'PS' THEN Device, ISOLATIN1=1
    
    ; Check keywords.
    capital = Keyword_Set(capital)
    
    ; If the first letter of the "symbol" variable is uppercase
    ; the user wants an uppercase symbol, if available.
    firstLetter = StrMid(symbol, 0, 1)
    IF firstLetter EQ StrUpCase(firstLetter) THEN capital = 1
    
    IF Keyword_Set(unicode) THEN BEGIN

        CASE StrLowCase(symbol) OF
            'alpha':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'beta':    retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'gamma':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'delta':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'epsilon': retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'zeta':    retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'eta':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'theta':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'iota':    retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'kappa':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'lambda':  retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'mu':      retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'nu':      retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'xi':      retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'omicron': retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'pi':      retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'rho':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'sigma':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'tau':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'upsilon': retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'phi':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'chi':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'psi':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'omega':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'leq':     retSymbol = '!Z(2266)'
            'geq':     retSymbol = '!Z(2267)'
            'neq':     retSymbol = '!Z(2260)'
            'deg':     retSymbol = '!Z(00B0)'
            '+-':      retSymbol = '!Z(00B1)'
            'equiv':   retSymbol = '!Z(2261)'
            'prime':   retSymbol = '!Z(2232)'
            'angstrom':retSymbol = '!Z(00C5)'
            'sun':     retSymbol = '!Z(2609)'
            'varphi':  retSymbol = '!Z(03D5)'
            'infinity':retSymbol = '!Z(221E)'
            'copyright': retSymbol = '!Z(00A9)'
            'up':      retSymbol = '!U'
            'down':    retSymbol = '!D'
            'exp':     retSymbol = '!E'
            'sub':     retSymbol = '!I'
            'n':       retSymbol = '!N'
            'div':     retSymbol = '!Z(00F7)'
            'times':   retSymbol = '!Z(00D7)'
            ELSE: Message, 'The symbol ' + symbol + ' is unrecognized.'
       ENDCASE
           
       ; Return the UNICODE letter.
       RETURN, retSymbol
       
    ENDIF
    
    IF ((!D.Name EQ 'PS') && (!P.FONT NE -1)) || Keyword_Set(PS)  THEN BEGIN
    
        ; Make sure ISOLATIN1 encoding is turned on.
        IF (!D.Name EQ 'PS') THEN DEVICE, /ISOLATIN1
        
        CASE StrLowCase(symbol) OF
            'alpha':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'beta':    retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'gamma':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'delta':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'epsilon': retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'zeta':    retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'eta':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'theta':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'iota':    retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'kappa':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'lambda':  retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'mu':      retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'nu':      retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'xi':      retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'omicron': retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'pi':      retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'rho':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'sigma':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'tau':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'upsilon': retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'phi':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'chi':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'psi':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'omega':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'leq':     retSymbol = '!9' + String("243B) + '!X'
            'geq':     retSymbol = '!9' + String("263B) + '!X'
            'neq':     retSymbol = '!9' + String("271B) + '!X' 
            'deg':     retSymbol = '!9' + String("260B) + '!X'
            '+-':      retSymbol = '!9' + String("261B) + '!X'
            'equiv':   retSymbol = '!9' + String("272B) + '!X'
            'prime':   retSymbol = '!9' + String("242B) + '!X'
            'angstrom':retSymbol = '!3' + String("305B) + '!X'
            'sun':     BEGIN
                       IF (ps_struct.font EQ 1) && (StrUpCase(ps_struct.tt_font) EQ 'DEJAVUSANS') THEN BEGIN
                          retSymbol = '!Z(2609)'
                       ENDIF ELSE BEGIN
                         thisDevice = !D.Name
                         Set_Plot, 'PS'
                         Device, /AVANTGARDE, ISOLATIN1=0, /BOOK, FONT_INDEX = 20
                         retSymbol = '!20!S!DO!R!I ' + string(183b) + '!X!N'
                         Set_Plot, thisDevice
                       ENDELSE
                       END
            'varphi':  retSymbol = '!9' + String("152B) + '!X'
            'infinity':retSymbol = '!9' + String("245B) + '!X'
            'copyright':retSymbol = '!9' + String("323B) + '!X'
            'up':      retSymbol = '!U'
            'down':    retSymbol = '!D'
            'exp':     retSymbol = '!E'
            'sub':     retSymbol = '!I'
            'n':       retSymbol = '!N'
            'div':     retSymbol = '!9' + String("270B) + '!X'
            'times':   retSymbol = '!9' + String("264B) + '!X'
           ELSE: Message, 'The symbol ' + symbol + ' is unrecognized.'
        ENDCASE
    
    ENDIF ELSE BEGIN
    
        CASE StrLowCase(symbol) OF
            'alpha':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'beta':    retSymbol = cgGreek(symbol, CAPITAL=capital)
            'gamma':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'delta':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'epsilon': retSymbol = cgGreek(symbol, CAPITAL=capital)
            'zeta':    retSymbol = cgGreek(symbol, CAPITAL=capital)
            'eta':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'theta':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'iota':    retSymbol = cgGreek(symbol, CAPITAL=capital)
            'kappa':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'lambda':  retSymbol = cgGreek(symbol, CAPITAL=capital)
            'mu':      retSymbol = cgGreek(symbol, CAPITAL=capital)
            'nu':      retSymbol = cgGreek(symbol, CAPITAL=capital)
            'xi':      retSymbol = cgGreek(symbol, CAPITAL=capital)
            'omicron': retSymbol = cgGreek(symbol, CAPITAL=capital)
            'pi':      retSymbol = cgGreek(symbol, CAPITAL=capital)
            'rho':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'sigma':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'tau':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'upsilon': retSymbol = cgGreek(symbol, CAPITAL=capital)
            'phi':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'chi':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'psi':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'omega':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'leq':     retSymbol = '!9' + String("154B) + '!X'
            'geq':     retSymbol = '!9' + String("142B) + '!X'
            'neq':     retSymbol = '!9' + String("75B) + '!X' 
            'deg':     retSymbol = '!9' + String("45B) + '!X'
            '+-':      retSymbol = '!9' + String("53B) + '!X'
            'equiv':   retSymbol = '!9' + String("72B) + '!X'
            'prime':   retSymbol = '!9' + String("140B) + '!X'
            'angstrom':retSymbol = '!3' + String("305B) + '!X'
            'sun':     retSymbol = '!D!9n!N!X'
            'varphi':  retSymbol = '!9' + String("120B) + '!X'
            'infinity':retSymbol = '!9' + String("44B) + '!X'
            'copyright':retSymbol = '!3' + String("251B) + '!X'
            'up':      retSymbol = '!U'
            'down':    retSymbol = '!D'
            'exp':     retSymbol = '!E'
            'sub':     retSymbol = '!I'
            'n':       retSymbol = '!N'
            'div':     retSymbol = '!9' + String("57B) + '!X'
            'times':   retSymbol = '!9' + String("130B) + '!X'
            ELSE: Message, 'The symbol ' + symbol + ' is unrecognized.'      
        ENDCASE
        
      ENDELSE

    RETURN, retSymbol
    
END

;;==============================================================================
; docformat = 'rst'
;
; NAME:
;   cgCheckForSymbols
;
; PURPOSE:
;   Checks a string for symbols that should be revolved with cgSymbol.
;
;******************************************************************************************;
;                                                                                          ;
;  Copyright (c) 2012, by Fanning Software Consulting, Inc. All rights reserved.           ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;
;+
; The purpose of this funciton is to check a string for symbols, encased in escape
; characters, that should be revolved with cgSymbol. The cgSymbol name will appear
; with the characters "$\" prepended to the name, and the character "$" appended. All
; Greek characters and other symbols supported by cgSymbol are allowed. Also,
; subscripts and superscripts are allowed::
;     $\upTEXT$ : Raise the text in TEXT to superscript level and reduce character size by 62%.
;     $\extTEXT$ : Raise the text in TEXT to superscript level and reduce character size by 44%.
;     $\downTEXT$ : Lower the text in TEXT to subscript level and reduce character size by 62%.
;     $\subTEXT$ : Lower the text in TEXT to subscript level and reduce character size by 44%.
;     
;
; :Categories:
;    Utilities
;    
; :Returns:
;    The modified string with the cgSymbol code embedded in place of the
;    escaped symbol name.
;    
; :Params:
;    astring: in, required, type=string
;       The string that should be searched for cgSymbol values.
;          
; :Examples:
;    To create a plot that uses the Greek mu character on the X axis and
;    the Angstrom squared symbol on the Y axis::
;    
;       cgPlot, cgDemoData(1), XTitle='Length ($\mu$M)', YTitle='Distance ($\Angstrom$$\up2$)'
;       
;    It is possible to use Greek characters as superscripts and subscripts. Do so by
;    prepending the Greek character with "\\" inside the normal superscript or subscript
;    notation. For example, to use lambda as an exponent to the Greek character Omega, you
;    can type this::
;    
;       cgPlot, cgDemoData(1), XTitle='$\Omega$$\exp\\lambda$', Charsize=2.0
;       
;    To use lambda as a subscript, type this:
;    
;        cgPlot, cgDemoData(1), XTitle='$\Omega$$\sub\\lambda$', Charsize=2.0
;
;    The program has been modified to accept TexToIDL tokens. They must be preceed by
;    a "\tex" prefix. For example, to draw a right arrow between 5 and 3, you would
;    construct the embedded string like this::
;
;         aString = '5 $\tex\rightarrow$ 3'
;         cgText, 0.5, 0.5, /Normal, Align=0.5, Charsize=3.0, aString
;       
; :Author:
;    FANNING SOFTWARE CONSULTING::
;       David W. Fanning 
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: david@idlcoyote.com
;       Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Written by David W. Fanning, 27 July 2012.
;        Modified to check for superscript and subscript codes. 9 November 2012. DWF.
;        Modified to allow the user to use the TexToIDL program from embedded codes.
;            To use a right arrow, for example, aString = '5 $\tex\rightarrow$ 3'
;        Added the ability to use Greek letters as subscripts and superscripts. See
;            the examples for details. 21 April 2013. DWF.
;
; :Copyright:
;     Copyright (c) 2012, Fanning Software Consulting, Inc.
;-
FUNCTION cgCheckForSymbols, aString

    Compile_Opt idl2
    
    ; Error handling.
    Catch, theError
    IF theError NE 0 THEN BEGIN
        Catch, /CANCEL
        void = Error_Message()
        IF N_Elements(aString) EQ 0 THEN RETURN, "" ELSE RETURN, aString
    ENDIF
    
    ; Must have a parameter.
    IF N_Elements(aString) EQ 0 THEN Message, 'Must pass a string.'
    
    ; What kind of thing is the parameter?
    type = Size(aString, /TNAME)
    
    
    ; If this is a string, then do your thing.
    IF type EQ 'STRING' THEN BEGIN
    
        FOR j=0,N_Elements(aString)-1 DO BEGIN
            thisString = aString[j]
            
            ; Can you find an escape sequence (i.e., "{\") in this string?
            locstart = StrPos(thisString, '$\')
            IF locStart NE -1 THEN BEGIN
            
               finalLoc = StrPos(StrMid(thisString, locstart+2), '$')
               IF finalLoc NE -1 THEN BEGIN
                  token = StrMid(thisString, locStart+2, finalLoc)
                  strToReplace = StrMid(thisString, locStart, finalLoc+3)
                  
                  ; Special handling for superscripts and subscripts. Check to see if there are embedded
                  ; special strings within the superscripts and subscripts. These are identified by "\\".
                  CASE 1 OF
                    
                      StrUpCase(StrMid(strToReplace, 2, 2)) EQ 'UP': BEGIN
                         check = StrPos(thisString, '\\')
                         IF check NE -1 THEN BEGIN
                            replaceStr = '!U' + cgCheckForSymbols('$\' + StrMid(strToReplace,6,StrLen(strToReplace)-5)) + '!N'
                         ENDIF ELSE BEGIN
                            replaceStr = '!U' + StrMid(strToReplace,4,StrLen(strToReplace)-5) + '!N'
                         ENDELSE
                         newString = StrMid(thisString, 0, locstart) + replaceStr + StrMid(thisString, locstart+3+StrLen(token))
                         END
                         
                      StrUpCase(StrMid(strToReplace, 2, 4)) EQ 'DOWN': BEGIN
                         check = StrPos(thisString, '\\')
                         IF check NE -1 THEN BEGIN
                            replaceStr = '!D' + cgCheckForSymbols('$\' + StrMid(strToReplace,8,StrLen(strToReplace)-7)) + '!N'
                         ENDIF ELSE BEGIN
                            replaceStr = '!D' + StrMid(strToReplace,6,StrLen(strToReplace)-7) + '!N'
                         ENDELSE
                         newString = StrMid(thisString, 0, locstart) + replaceStr + StrMid(thisString, locstart+3+StrLen(token))
                         END
                         
                      StrUpCase(StrMid(strToReplace, 2, 3)) EQ 'EXP': BEGIN
                         check = StrPos(thisString, '\\')
                         IF check NE -1 THEN BEGIN
                            replaceStr = '!E' + cgCheckForSymbols('$\' + StrMid(strToReplace,7,StrLen(strToReplace)-6)) + '!N'
                         ENDIF ELSE BEGIN
                            replaceStr = '!E' + StrMid(strToReplace,5,StrLen(strToReplace)-6) + '!N'
                         ENDELSE
                         newString = StrMid(thisString, 0, locstart) + replaceStr + StrMid(thisString, locstart+3+StrLen(token))
                         END
                         
                      StrUpCase(StrMid(strToReplace, 2, 3)) EQ 'SUB': BEGIN
                         check = StrPos(thisString, '\\')
                         IF check NE -1 THEN BEGIN
                            replaceStr = '!I' + cgCheckForSymbols('$\' + StrMid(strToReplace,7,StrLen(strToReplace)-6)) + '!N'
                         ENDIF ELSE BEGIN
                            replaceStr = '!I' + StrMid(strToReplace,5,StrLen(strToReplace)-6) + '!N'
                         ENDELSE
                         newString = StrMid(thisString, 0, locstart) + replaceStr + StrMid(thisString, locstart+3+StrLen(token))
                         END

                      ; Provide a mechanism whereby people can use TexToIDL.
                      StrUpCase(StrMid(strToReplace, 2, 3)) EQ 'TEX': BEGIN
                         replaceStr = Call_Function('TexToIDL', StrMid(token,3), FONT=!P.FONT)
                         newString = StrMid(thisString, 0, locstart) + replaceStr + StrMid(thisString, locstart+3+StrLen(token))
                         END
                     ELSE: newString = StrMid(thisString, 0, locstart) + cgSymbol(token) + StrMid(thisString, locstart+3+StrLen(token))

                  ENDCASE
                  
                  thisString = cgCheckForSymbols(newString)
               ENDIF
               
            ENDIF
            
            aString[j] = thisString
         
        ENDFOR
        
        RETURN, aString
        
    ENDIF

END ;----------------------------------------------------------------------------------------------------------------
;;==============================================================================
; docformat = 'rst'
;
; NAME:
;   cgColor
;
; PURPOSE:
; The purpose of this function is to obtain drawing colors
; by name and in a device/decomposition independent way.
; The color names and values may be read in as a file, or 192 color
; names and values are supplied with the program. These colors were
; obtained from the file rgb.txt, found on most X-Window distributions,
; and from colors in the Brewer color tables (http://colorbrewer2.org/).
; Representative colors were chosen from across the color spectrum. 
; If the color names '0', '1', '2', ..., '255' are used, they will
; correspond to the colors in the current color table in effect at
; the time the cgColor program is called.
;
;******************************************************************************************;
;                                                                                          ;
;  Copyright (c) 2011, by Fanning Software Consulting, Inc. All rights reserved.           ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;
;+
; The purpose of this function is to obtain drawing colors
; by name and in a device/decomposition independent way.
; The color names and values may be read in as a file, or 192 color
; names and values are supplied with the program. These colors were
; obtained from the file rgb.txt, found on most X-Window distributions,
; and from colors in the `Brewer color tables <http://colorbrewer2.org/>`.
; Representative colors were chosen from across the color spectrum. 
; If the color names '0', '1', '2', ..., '255' are used, they will
; correspond to the colors in the current color table in effect at
; the time the `cgColor` program is called.
; 
; Please note that all Coyote Graphics routines use cgColor internally to specify
; their colors in a color-model independent way. It is not necessary to call
; cgColor yourself unless you are using it with a traditional IDL command (e.g., Plot).
; For example::
;  
;     Plot, data, Color=cgColor('dodger blue')
;     
; But, it is not needed with Coyote Graphics commands::
; 
;     cgPlot, data, Color='dodger blue'
; 
; The program requires the `Coyote Library <http://www.idlcoyote.com/documents/programs.php>`
; to be installed on your machine.
;
; :Categories:
;    Graphics
;    
; :Examples:
;    To get drawing colors in a device-decomposed independent way::
;
;        axisColor = cgColor("Green", !D.Table_Size-2)
;        backColor = cgColor("Charcoal", !D.Table_Size-3)
;        dataColor = cgColor("Yellow", !D.Table_Size-4)
;        Plot, Findgen(11), Color=axisColor, Background=backColor, /NoData
;        OPlot, Findgen(11), Color=dataColor
;
;    To set the viewport color in object graphics::
;
;        theView = Obj_New('IDLgrView', Color=cgColor('Charcoal', /Triple))
;
;    To change the viewport color later::
;
;        theView->SetProperty, Color=cgColor('Antique White', /Triple)
;
;    To load the drawing colors "red", "green", and "yellow" at indices 100-102, type this::
;
;        IDL> TVLCT, cgColor(["red", "green", "yellow"], /Triple), 100
;           
;    To interactively choose a color, set the SELECTCOLOR keyword::
;    
;        IDL> color = cgColor(/SelectColor)
;        
;    The cgPickColorName program is a good way to learn the names of the colors available::
;    
;        IDL> color = cgPickColorName()
;
; .. image:: cgpickcolorname.png
;
; :Author:
;    FANNING SOFTWARE CONSULTING::
;       David W. Fanning 
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: david@idlcoyote.com
;       Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Written by: David W. Fanning
;        Modified FSC_COLOR to create cgColor 9 February 2011. DWF.
;        Modified to allow a three-element color triple to be used in place of the color
;           name parameter. This allows one user-defined color to be used. 4 Dec 2011. DWF.
;        Modified to allow byte and 16-bit integer values to be used to specify colors
;           in the current color table. 5 Dec 2011. DWF.
;        Modified to allow the "opposite" pixel to be determined in the Z-graphics buffer. 24 Dec 2011. DWF.
;        Modified the code to handle long integers depending on the current color mode and the
;            number of values passed in. 10 January 2012. DWF.
;        Made sure the return values are BYTES not INTEGERS, in cases where this is expected. 10 Jan 2012. DWF.
;        Added "Background" as a color name. The opposite of "Opposite". 1 Feb 2012. DWF.
;        When returning a vector of color values, now making sure to return a byte array if 
;             in indexed color mode. 27 Feb 2012. DWF.
;        Added Compile Opt id2 to all file modules. 22 July 2012. DWF.
;        Added "opposite" and "background" colors to Brewer colors. 14 August 2012. DWF.
;        Some versions of IDL report the size of widget windows incorrectly, so instead of
;              sampling the very top-right pixel, I now back off a little. 1 Nov 2012. DWF.
;        For numerical values less than 256, in indexed color state, I now return the values
;              directly to the user. This should significantly speed up many Coyote Graphics
;              processes. 14 December 2012. DWF.
;        Removed cgColor_Color24 module in favor of using Coyote Library routine cgColor24. 5 Jan 2013. DWF.
;        The keyword ROW was being ignored if multiple colors were specified with TRIPLE keyword. Fixed. 10 July 2013. DWF.
;        
; :Copyright:
;     Copyright (c) 2009-2013, Fanning Software Consulting, Inc.
;-
;
;+
; The purpose of this function is to obtain drawing colors
; by name and in a device and color model independent way.
; 
; :Returns:
;     The return value depends on the color mode in effect at the time
;     the program is called and which keyword is used with the program.
;     In normal operation, if the graphics device is using indexed color
;     mode, the program will load a color at a unique (or specified)
;     index and return that index number. If the graphics device is using
;     decomposed color mode, the program will create a 24-bit color value
;     that can be used to specify the particular color desired. In this
;     case, no color is loaded in the color table. This is the preferred
;     mode for working with colors in IDL.
;     
; :Params:
;    theColour: required, optional, type=varies
;        Normally the name of the color desired. However, this can also be
;        a string index number (e.g., '215') or a byte or short integer
;        value (e.g, 215B or 215S). If this is the case, the color
;        in the current color table at this index number is used for the 
;        color that is returned. The value may also be a vector of color names. 
;        The color may also be a three-element byte or integer array specifying a 
;        user-defined color triple. Only one color triple is allowed.
;
;        To see a list of the color names available set the NAMES keyword. Colors available are these::
;
;           Active            Almond     Antique White        Aquamarine             Beige            Bisque
;           Black               Blue       Blue Violet             Brown         Burlywood        Cadet Blue
;           Charcoal       Chartreuse         Chocolate             Coral   Cornflower Blue          Cornsilk
;           Crimson              Cyan    Dark Goldenrod         Dark Gray        Dark Green        Dark Khaki
;           Dark Orchid      Dark Red       Dark Salmon   Dark Slate Blue         Deep Pink       Dodger Blue
;           Edge                 Face         Firebrick      Forest Green             Frame              Gold
;           Goldenrod            Gray             Green      Green Yellow         Highlight          Honeydew
;           Hot Pink       Indian Red             Ivory             Khaki          Lavender        Lawn Green
;           Light Coral    Light Cyan        Light Gray      Light Salmon   Light Sea Green      Light Yellow
;           Lime Green          Linen           Magenta            Maroon       Medium Gray     Medium Orchid
;           Moccasin             Navy             Olive        Olive Drab            Orange        Orange Red
;           Orchid     Pale Goldenrod        Pale Green            Papaya              Peru              Pink
;           Plum          Powder Blue            Purple               Red              Rose        Rosy Brown
;           Royal Blue   Saddle Brown            Salmon       Sandy Brown         Sea Green          Seashell
;           Selected           Shadow            Sienna          Sky Blue        Slate Blue        Slate Gray
;           Snow         Spring Green        Steel Blue               Tan              Teal              Text
;           Thistle            Tomato         Turquoise            Violet        Violet Red             Wheat
;           White              Yellow
;
;        Here are the Brewer color names::
;
;           WT1        WT2       WT3       WT4       WT5       WT6       WT7       WT8
;           TAN1      TAN2      TAN3      TAN4      TAN5      TAN6      TAN7      TAN8
;           BLK1      BLK2      BLK3      BLK4      BLK5      BLK6      BLK7      BLK8
;           GRN1      GRN2      GRN3      GRN4      GRN5      GRN6      GRN7      GRN8
;           BLU1      BLU2      BLU3      BLU4      BLU5      BLU6      BLU7      BLU8
;           ORG1      ORG2      ORG3      ORG4      ORG5      ORG6      ORG7      ORG8
;           RED1      RED2      RED3      RED4      RED5      RED6      RED7      RED8
;           PUR1      PUR2      PUR3      PUR4      PUR5      PUR6      PUR7      PUR8
;           PBG1      PBG2      PBG3      PBG4      PBG5      PBG6      PBG7      PBG8
;           YGB1      YGB2      YGB3      YGB4      YGB5      YGB6      YGB7      YGB8
;           RYB1      RYB2      RYB3      RYB4      RYB5      RYB6      RYB7      RYB8
;           TG1        TG2       TG3       TG4       TG5       TG6       TG7       TG8
;            
;        The color name "OPPOSITE" is also available. It chooses a color "opposite" to the 
;        color of the pixel in the upper-right corner of the display, if a window is open.
;        Otherwise, this color is "black" in PostScript and "white" everywhere else.
;        The color OPPOSITE is used if this parameter is absent or a color name is mis-spelled.
;        
;         The color name "BACKGROUND" can similarly be used to select the color of the pixel
;         in the upper-right corner of the display, if a window is open.
;           
;    colorindex: in, optional, type=byte
;        The color table index where the color should be loaded. Colors are
;        loaded into the color table only if using indexed color mode in the
;        current graphics device. If this parameter is missing, the color will
;        be loaded at a unique color index number, if necessary.
;        
; :Keywords:
;     allcolors: in, optional, type=boolean, default=0
;        Set this keyword to return indices, or 24-bit values, or color
;        triples, for all the known colors, instead of for a single color.
;     brewer: in, optional, type=boolean, default=0
;        An obsolete keyword. If used, only Brewer colors are loaded into the color
;        vectors internally.
;     cancel: out, optional, type=boolean, default=0
;        This keyword is always set to 0, unless that SELECTCOLOR keyword is used.
;        Then it will correspond to the value of the CANCEL output keyword in cgPickColorName.
;     check_connection: in, optional, type=boolean, default=0
;         An obsolete keyword now completely ignored.
;     colorstructure: out, optional, type=structure
;        This output keyword (if set to a named variable) will return a
;        structure in which the fields will be the known color names (without spaces)
;        and the values of the fields will be either color table index numbers or
;        24-bit color values. If you have specified a vector of color names, then
;        this will be a structure containing just those color names as fields.
;     decomposed: in, optional, type=boolean
;        Set this keyword to 0 or 1 to force the return value to be
;        a color table index or a 24-bit color value, respectively. This
;        keyword is normally set by the color state of the current graphics device.
;     filename: in, optional, type=string
;        The  name of an ASCII file that can be opened to read in color values and color 
;        names. There should be one color per row in the file. Please be sure there are 
;        no blank lines in the file. The format of each row should be::
;
;           redValue  greenValue  blueValue  colorName
;
;        Color values should be between 0 and 255. Any kind of white-space
;        separation (blank characters, commas, or tabs) are allowed. The color
;        name should be a string, but it should NOT be in quotes. A typical
;        entry into the file would look like this::
;
;           255   255   0   Yellow
;     names: in, optional, type=boolian, default=0
;        If this keyword is set, the return value of the function is a string array 
;        containing the names of the colors. These names would be appropriate, for example, 
;        in building a list widget with the names of the colors. If the NAMES
;        keyword is set, the COLOR and INDEX parameters are ignored.
;     ncolors: out, optional, type=integer
;        Returns the number of colors that cgColor "knows" about. Currently ncolors=193.
;     nodisplay: in, optional, type=boolean, default=0
;        An obsolete keyword, now totally ignored.
;     row: in, optional, type=boolean
;        If this keyword is set, the return value of the function when the TRIPLE
;        keyword is set is returned as a row vector, rather than as the default
;        column vector. This is required, for example, when you are trying to
;        use the return value to set the color for object graphics objects. This
;        keyword is completely ignored, except when used in combination with the
;        TRIPLE keyword.
;     selectcolor: in, optional, type=boolean
;       Set this keyword if you would like to select the color name with
;       the cgPickColorName program. Selecting this keyword automaticallys sets
;       the INDEX positional parameter. If this keyword is used, any keywords
;       appropriate for cgPickColorName can also be used. If this keyword is used,
;       the first positional parameter can be a color name that will appear in
;       the SelectColor box.
;     triple: in, optional, type=boolean
;        Setting this keyword will force the return value of the function to
;        always be a color triple, regardless of color decomposition state or
;        visual depth of the machine. The value will be a three-element column
;        vector unless the ROW keyword is also set.
;     _ref_extra: in, optional
;        Any keyword parameter appropriate for cgPickColorName can be used.
;       These include BOTTOM, COLUMNS, GROUP_LEADER, INDEX, and TITLE.
;
;-
FUNCTION cgColor, theColour, colorIndex, $
   ALLCOLORS=allcolors, $
   BREWER=brewer, $ ; This keyword is no longer used.
   CHECK_CONNECTION=check_connection, $ ; This keyword is completely ignored.
   COLORSTRUCTURE=colorStructure, $
   CANCEL=cancelled, $
   DECOMPOSED=decomposedState, $
   FILENAME=filename, $
   NAMES=names, $
   NCOLORS=ncolors, $
   NODISPLAY=nodisplay, $ ; This keyword is completely ignored.
   ROW=row, $
   SELECTCOLOR=selectcolor, $
   TRIPLE=triple, $
  _REF_EXTRA=extra
  
    Compile_Opt idl2
   
    ; Return to caller as the default error behavior.
    On_Error, 2
        
    ; Error handling for the rest of the program.
    Catch, theError
    IF theError NE 0 THEN BEGIN
       Catch, /Cancel
       ok = Error_Message()
       cancelled = 1
       RETURN, !P.Color
    ENDIF
    
    ; Get the current color state. This will help you determine what to 
    ; do with the input color.
    IF N_Elements(decomposedState) NE 0 THEN colorState = Keyword_Set(decomposedState) ELSE colorState = GetDecomposedState()
    
    ; Set up PostScript device for working with colors.
    IF !D.Name EQ 'PS' THEN Device, COLOR=1, BITS_PER_PIXEL=8
    
    ; I don't want to change the original variable.
    IF N_Elements(theColour) NE 0 THEN theColor = theColour ELSE $
        theColor = 'OPPOSITE'
        
     ; Allow the color values to be something other than a string.
     ; There will be some ambiguity between a color triple and a number
     ; array of three elements, but I am willing to live with this.
     IF Size(theColor, /TNAME) NE 'STRING' THEN BEGIN
     
        ; Make sure this is not a 1x3 array, which we want to treat as a color triple.
        IF (N_Elements(theColor) EQ 3) && (Size(theColor, /N_DIMENSIONS) EQ 2) THEN BEGIN
            theColor = Reform(theColor)
        ENDIF
     
        ; Allow the color to be a three-element array of byte values.
        ; If it is, we will define the USERDEF color with these values.
        ; Otherwise the USERDEF color will be unused.
        IF (Size(theColor, /N_DIMENSIONS) EQ 1) && $
           (N_Elements(theColor) EQ 3) && $
           (Max(theColor) LE 255) && $
           (Min(theColor) GE 0) THEN BEGIN
           usercolor = theColor
           theColor = 'USERDEF'
        ENDIF
        
        ; If the input didn't qualify as a color triple, then see if you 
        ; can somehow make sense of the number values.
        IF Size(theColor, /TNAME) NE 'STRING' THEN BEGIN
        
          ; We can assume that any number that is a byte or short integer must
          ; be an index into the color table. Return that value directly, if
          ; you are currently in an indexed color state.
          IF (Size(theColor, /TYPE) LE 2) THEN BEGIN
             IF (colorState EQ 1) THEN theColor = StrTrim(Fix(theColor),2) ELSE RETURN, theColor
          ENDIF 
          
          ; Long integers are problematic. If the current color mode is INDEXED, then
          ; we will treat long integers as color indices. If it is DECOMPOSED, then if
          ; there is just one value, we can handle this as a color triple.
          IF (Size(theColor, /TYPE) EQ 3) THEN BEGIN
             
               IF (colorState EQ 1) THEN BEGIN
                   IF N_Elements(theColor) EQ 1 THEN BEGIN
                      usercolor = [theColor MOD 2L^8, (theColor MOD 2L^16)/2L^8, theColor/2L^16]
                      theColor = 'USERDEF'
                   ENDIF ELSE Message, 'Do not know how to handle a vector of LONG integers!
               ENDIF ELSE BEGIN
                   IF N_Elements(theColor) EQ 1 THEN BEGIN
                      IF theColor LE 255 THEN BEGIN
                          RETURN, theColor
                      ENDIF ELSE Message, 'Long integer ' + StrTrim(theColor,2) + ' is out of indexed color range.'
                   ENDIF ELSE Message, 'Do not know how to handle a vector of LONG integers!
               ENDELSE
               
          ENDIF
          
          ; Anything that is not an BYTE, INTEGER, LONG, or STRING causes problems.
          IF (Size(theColor, /TYPE) GT 4) && (Size(theColor, /TNAME) NE 'STRING') THEN BEGIN
             IF (colorstate EQ 0) AND (theColor LE 255) THEN BEGIN
                RETURN, theColor
             ENDIF ELSE BEGIN
                Message, 'Use BYTE, INTEGER, or STRING data to specify a color.'
             ENDELSE
          ENDIF
        ENDIF
     ENDIF
        
    ; Make sure the color parameter is a string.
    varName = Size(theColor, /TNAME)
    IF varName NE 'STRING' THEN $
       Message, 'The color name parameter must be a string.', /NoName
       
    ; We don't want underscores in color names. Turn all underscores
    ; to spaces.
    FOR j=0,N_Elements(theColor)-1 DO BEGIN
        theColor[j] = StrJoin( StrSplit(theColor[j], '_', /Extract, $
           /Preserve_Null), ' ')
    ENDFOR
    
    ; Make sure the color is compressed and uppercase.   
    theColor = StrUpCase(StrCompress(StrTrim(theColor,2), /Remove_All))
    
    ; Get the pixel value of the "opposite" color. This is the pixel color
    ; opposite the pixel color in the upper right corner of the display.
    ; Because Windows versions of IDL (through at least IDL 8.2.1) report
    ; the size of draw widget windows inaccurately until you make them the
    ; current graphics window, I back off from the very corner pixel by two
    ; pixels to read Windows windows.
    IF ((!D.Window GE 0) && ((!D.Flags AND 256) NE 0)) || (!D.Name EQ 'Z') THEN BEGIN
       IF StrUpCase(!Version.OS_Family) EQ 'WINDOWS' THEN BEGIN
          opixel = cgSnapshot(!D.X_Size-3, !D.Y_Size-3, 1, 1)
       ENDIF ELSE BEGIN
          opixel = cgSnapshot(!D.X_Size-1, !D.Y_Size-1, 1, 1)
       ENDELSE
       IF N_Elements(opixel) NE 3 THEN BEGIN
            IF (!D.Name NE 'NULL') THEN TVLCT, rrr, ggg, bbb, /Get
            opixel = [rrr[opixel], ggg[opixel], bbb[opixel]]
       ENDIF
    ENDIF ELSE BEGIN
       IF (!D.Name EQ 'PS') THEN opixel = [255,255,255] ELSE opixel = [0,0,0]
    ENDELSE
    IF N_Elements(opixel) EQ 0 THEN opixel = [0,0,0]
    bgcolor = opixel
    opixel = 255 - opixel
    
    ; Read the first color as bytes. If none of the bytes are less than 48
    ; or greater than 57, then this is a "number" string and you should
    ; assume the current color table is being used.
    bytecheck = Byte(theColor[0])
    i = Where(bytecheck LT 48, lessthan)
    i = Where(bytecheck GT 57, greaterthan)
    IF (lessthan + greaterthan) EQ 0 THEN useCurrentColors = 1 ELSE useCurrentColors = 0
    
    ; Get the decomposed state of the IDL session right now.
    IF N_Elements(decomposedState) EQ 0 THEN BEGIN
       IF Float(!Version.Release) GE 5.2 THEN BEGIN
          IF (!D.Name EQ 'X' OR !D.Name EQ 'WIN' OR !D.Name EQ 'MAC') THEN BEGIN
             Device, Get_Decomposed=decomposedState
          ENDIF ELSE decomposedState = 0
       ENDIF ELSE decomposedState = 0
       IF (Float(!Version.Release) GE 6.4) AND (!D.NAME EQ 'Z') THEN BEGIN
          Device, Get_Decomposed=decomposedState, Get_Pixel_Depth=theDepth
          IF theDepth LT 24 THEN decomposedState = 0
       ENDIF
    ENDIF ELSE decomposedState = Keyword_Set(decomposedState)
    
    ; Get depth of visual display (and decomposed state for PostScript devices).
    IF (!D.Flags AND 256) NE 0 THEN Device, Get_Visual_Depth=theDepth ELSE theDepth = 8
    IF (Float(!Version.Release) GE 6.4) AND (!D.NAME EQ 'Z') THEN Device, Get_Pixel_Depth=theDepth
    IF (!D.NAME EQ 'PS') AND (Float(!Version.Release) GE 7.1) THEN BEGIN
       decomposedState = DecomposedColor(DEPTH=theDepth)
    ENDIF

    ; Need brewer colors?
    brewer = Keyword_Set(brewer)
    
    ; Load the colors.
    IF N_Elements(filename) NE 0 THEN BEGIN
    
       ; Count the number of rows in the file.
       ncolors = File_Lines(filename)
    
       ; Read the data.
       OpenR, lun, filename, /Get_Lun
       rvalue = BytArr(NCOLORS)
       gvalue = BytArr(NCOLORS)
       bvalue = BytArr(NCOLORS)
       colors = StrArr(NCOLORS)
       redvalue = 0B
       greenvalue = 0B
       bluevalue = 0B
       colorvalue = ""
       FOR j=0L, NCOLORS-1 DO BEGIN
          ReadF, lun, redvalue, greenvalue, bluevalue, colorvalue
          rvalue[j] = redvalue
          gvalue[j] = greenvalue
          bvalue[j] = bluevalue
          colors[j] = colorvalue
       ENDFOR
       Free_Lun, lun
    
       ; Trim the colors array of blank characters.
       colors = StrTrim(colors, 2)
    
    ENDIF ELSE BEGIN
    
       ; Set up the color vectors.
       IF Keyword_Set(Brewer) THEN BEGIN
       
           ; Set up the color vectors.
           colors = [ 'WT1', 'WT2', 'WT3', 'WT4', 'WT5', 'WT6', 'WT7', 'WT8']
           rvalue = [  255,   255,   255,   255,   255,   245,   255,   250 ]
           gvalue = [  255,   250,   255,   255,   248,   245,   245,   240 ]
           bvalue = [  255,   250,   240,   224,   220,   220,   238,   230 ]
           colors = [ colors, 'TAN1', 'TAN2', 'TAN3', 'TAN4', 'TAN5', 'TAN6', 'TAN7', 'TAN8']
           rvalue = [ rvalue,   250,   255,    255,    255,    255,    245,    222,    210 ]
           gvalue = [ gvalue,   235,   239,    235,    228,    228,    222,    184,    180 ]
           bvalue = [ bvalue,   215,   213,    205,    196,    181,    179,    135,    140 ]
           colors = [ colors, 'BLK1', 'BLK2', 'BLK3', 'BLK4', 'BLK5', 'BLK6', 'BLK7', 'BLK8']
           rvalue = [ rvalue,   250,   230,    210,    190,    128,     110,    70,       0 ]
           gvalue = [ gvalue,   250,   230,    210,    190,    128,     110,    70,       0 ]
           bvalue = [ bvalue,   250,   230,    210,    190,    128,     110,    70,       0 ]
           colors = [ colors, 'GRN1', 'GRN2', 'GRN3', 'GRN4', 'GRN5', 'GRN6', 'GRN7', 'GRN8']
           rvalue = [ rvalue,   250,   223,    173,    109,     53,     35,      0,       0 ]
           gvalue = [ gvalue,   253,   242,    221,    193,    156,     132,    97,      69 ]
           bvalue = [ bvalue,   202,   167,    142,    115,     83,      67,    52,      41 ]
           colors = [ colors, 'BLU1', 'BLU2', 'BLU3', 'BLU4', 'BLU5', 'BLU6', 'BLU7', 'BLU8']
           rvalue = [ rvalue,   232,   202,    158,     99,     53,     33,      8,       8 ]
           gvalue = [ gvalue,   241,   222,    202,    168,    133,    113,     75,      48 ]
           bvalue = [ bvalue,   250,   240,    225,    211,    191,    181,    147,     107 ]
           colors = [ colors, 'ORG1', 'ORG2', 'ORG3', 'ORG4', 'ORG5', 'ORG6', 'ORG7', 'ORG8']
           rvalue = [ rvalue,   254,    253,    253,    250,    231,    217,    159,    127 ]
           gvalue = [ gvalue,   236,    212,    174,    134,     92,     72,     51,     39 ]
           bvalue = [ bvalue,   217,    171,    107,     52,     12,      1,      3,      4 ]
           colors = [ colors, 'RED1', 'RED2', 'RED3', 'RED4', 'RED5', 'RED6', 'RED7', 'RED8']
           rvalue = [ rvalue,   254,    252,    252,    248,    225,    203,    154,    103 ]
           gvalue = [ gvalue,   232,    194,    146,     97,     45,     24,     12,      0 ]
           bvalue = [ bvalue,   222,    171,    114,     68,     38,     29,     19,     13 ]
           colors = [ colors, 'PUR1', 'PUR2', 'PUR3', 'PUR4', 'PUR5', 'PUR6', 'PUR7', 'PUR8']
           rvalue = [ rvalue,   244,    222,    188,    152,    119,    106,     80,     63 ]
           gvalue = [ gvalue,   242,    221,    189,    148,    108,     82,     32,      0 ]
           bvalue = [ bvalue,   248,    237,    220,    197,    177,    163,    139,    125 ]
           colors = [ colors, 'PBG1', 'PBG2', 'PBG3', 'PBG4', 'PBG5', 'PBG6', 'PBG7', 'PBG8']
           rvalue = [ rvalue,   243,    213,    166,     94,     34,      3,      1,      1 ]
           gvalue = [ gvalue,   234,    212,    189,    164,    138,    129,    101,     70 ]
           bvalue = [ bvalue,   244,    232,    219,    204,    171,    139,     82,     54 ]
           colors = [ colors, 'YGB1', 'YGB2', 'YGB3', 'YGB4', 'YGB5', 'YGB6', 'YGB7', 'YGB8']
           rvalue = [ rvalue,   244,    206,    127,     58,     30,     33,     32,      8 ]
           gvalue = [ gvalue,   250,    236,    205,    175,    125,     95,     48,     29 ]
           bvalue = [ bvalue,   193,    179,    186,    195,    182,    168,    137,     88 ]
           colors = [ colors, 'RYB1', 'RYB2', 'RYB3', 'RYB4', 'RYB5', 'RYB6', 'RYB7', 'RYB8']
           rvalue = [ rvalue,   201,    245,    253,    251,    228,    193,    114,     59 ]
           gvalue = [ gvalue,    35,    121,    206,    253,    244,    228,    171,     85 ]
           bvalue = [ bvalue,    38,    72,     127,    197,    239,    239,    207,    164 ]
           colors = [ colors, 'TG1', 'TG2', 'TG3', 'TG4', 'TG5', 'TG6', 'TG7', 'TG8']
           rvalue = [ rvalue,  84,    163,   197,   220,   105,    51,    13,     0 ]
           gvalue = [ gvalue,  48,    103,   141,   188,   188,   149,   113,    81 ]
           bvalue = [ bvalue,   5,     26,    60,   118,   177,   141,   105,    71 ]
           colors = [ colors, 'OPPOSITE', 'BACKGROUND']
           rvalue = [ rvalue,  opixel[0],  bgcolor[0]]
           gvalue = [ gvalue,  opixel[1],  bgcolor[1]]
           bvalue = [ bvalue,  opixel[2],  bgcolor[2]]
       
       ENDIF ELSE BEGIN
       
           ; Set up the color vectors. Both original and Brewer colors.
           colors= ['White']
           rvalue = [ 255]
           gvalue = [ 255]
           bvalue = [ 255]
           colors = [ colors,   'Snow',     'Ivory','Light Yellow', 'Cornsilk',     'Beige',  'Seashell' ]
           rvalue = [ rvalue,     255,         255,       255,          255,          245,        255 ]
           gvalue = [ gvalue,     250,         255,       255,          248,          245,        245 ]
           bvalue = [ bvalue,     250,         240,       224,          220,          220,        238 ]
           colors = [ colors,   'Linen','Antique White','Papaya',     'Almond',     'Bisque',  'Moccasin' ]
           rvalue = [ rvalue,     250,        250,        255,          255,          255,          255 ]
           gvalue = [ gvalue,     240,        235,        239,          235,          228,          228 ]
           bvalue = [ bvalue,     230,        215,        213,          205,          196,          181 ]
           colors = [ colors,   'Wheat',  'Burlywood',    'Tan', 'Light Gray',   'Lavender','Medium Gray' ]
           rvalue = [ rvalue,     245,        222,          210,      230,          230,         210 ]
           gvalue = [ gvalue,     222,        184,          180,      230,          230,         210 ]
           bvalue = [ bvalue,     179,        135,          140,      230,          250,         210 ]
           colors = [ colors,  'Gray', 'Slate Gray',  'Dark Gray',  'Charcoal',   'Black',  'Honeydew', 'Light Cyan' ]
           rvalue = [ rvalue,      190,      112,          110,          70,         0,         240,          224 ]
           gvalue = [ gvalue,      190,      128,          110,          70,         0,         255,          255 ]
           bvalue = [ bvalue,      190,      144,          110,          70,         0,         255,          240 ]
           colors = [ colors,'Powder Blue',  'Sky Blue', 'Cornflower Blue', 'Cadet Blue', 'Steel Blue','Dodger Blue', 'Royal Blue',  'Blue' ]
           rvalue = [ rvalue,     176,          135,          100,              95,            70,           30,           65,            0 ]
           gvalue = [ gvalue,     224,          206,          149,             158,           130,          144,          105,            0 ]
           bvalue = [ bvalue,     230,          235,          237,             160,           180,          255,          225,          255 ]
           colors = [ colors,  'Navy', 'Pale Green','Aquamarine','Spring Green',  'Cyan' ]
           rvalue = [ rvalue,        0,     152,          127,          0,            0 ]
           gvalue = [ gvalue,        0,     251,          255,        250,          255 ]
           bvalue = [ bvalue,      128,     152,          212,        154,          255 ]
           colors = [ colors, 'Turquoise', 'Light Sea Green', 'Sea Green','Forest Green',  'Teal','Green Yellow','Chartreuse', 'Lawn Green' ]
           rvalue = [ rvalue,      64,          143,               46,          34,             0,      173,           127,         124 ]
           gvalue = [ gvalue,     224,          188,              139,         139,           128,      255,           255,         252 ]
           bvalue = [ bvalue,     208,          143,               87,          34,           128,       47,             0,           0 ]
           colors = [ colors, 'Green', 'Lime Green', 'Olive Drab',  'Olive','Dark Green','Pale Goldenrod']
           rvalue = [ rvalue,      0,        50,          107,        85,            0,          238 ]
           gvalue = [ gvalue,    255,       205,          142,       107,          100,          232 ]
           bvalue = [ bvalue,      0,        50,           35,        47,            0,          170 ]
           colors = [ colors,     'Khaki', 'Dark Khaki', 'Yellow',  'Gold', 'Goldenrod','Dark Goldenrod']
           rvalue = [ rvalue,        240,       189,        255,      255,      218,          184 ]
           gvalue = [ gvalue,        230,       183,        255,      215,      165,          134 ]
           bvalue = [ bvalue,        140,       107,          0,        0,       32,           11 ]
           colors = [ colors,'Saddle Brown',  'Rose',   'Pink', 'Rosy Brown','Sandy Brown', 'Peru']
           rvalue = [ rvalue,     139,          255,      255,        188,        244,        205 ]
           gvalue = [ gvalue,      69,          228,      192,        143,        164,        133 ]
           bvalue = [ bvalue,      19,          225,      203,        143,         96,         63 ]
           colors = [ colors,'Indian Red',  'Chocolate',  'Sienna','Dark Salmon',   'Salmon','Light Salmon' ]
           rvalue = [ rvalue,    205,          210,          160,        233,          250,       255 ]
           gvalue = [ gvalue,     92,          105,           82,        150,          128,       160 ]
           bvalue = [ bvalue,     92,           30,           45,        122,          114,       122 ]
           colors = [ colors,  'Orange',      'Coral', 'Light Coral',  'Firebrick', 'Dark Red', 'Brown',  'Hot Pink' ]
           rvalue = [ rvalue,       255,         255,        240,          178,        139,       165,        255 ]
           gvalue = [ gvalue,       165,         127,        128,           34,          0,        42,        105 ]
           bvalue = [ bvalue,         0,          80,        128,           34,          0,        42,        180 ]
           colors = [ colors, 'Deep Pink',    'Magenta',   'Tomato', 'Orange Red',   'Red', 'Crimson', 'Violet Red' ]
           rvalue = [ rvalue,      255,          255,        255,        255,          255,      220,        208 ]
           gvalue = [ gvalue,       20,            0,         99,         69,            0,       20,         32 ]
           bvalue = [ bvalue,      147,          255,         71,          0,            0,       60,        144 ]
           colors = [ colors,    'Maroon',    'Thistle',       'Plum',     'Violet',    'Orchid','Medium Orchid']
           rvalue = [ rvalue,       176,          216,          221,          238,         218,        186 ]
           gvalue = [ gvalue,        48,          191,          160,          130,         112,         85 ]
           bvalue = [ bvalue,        96,          216,          221,          238,         214,        211 ]
           colors = [ colors,'Dark Orchid','Blue Violet',  'Purple']
           rvalue = [ rvalue,      153,          138,       160]
           gvalue = [ gvalue,       50,           43,        32]
           bvalue = [ bvalue,      204,          226,       240]
           colors = [ colors, 'Slate Blue',  'Dark Slate Blue']
           rvalue = [ rvalue,      106,            72]
           gvalue = [ gvalue,       90,            61]
           bvalue = [ bvalue,      205,           139]
           colors = [ colors, 'WT1', 'WT2', 'WT3', 'WT4', 'WT5', 'WT6', 'WT7', 'WT8']
           rvalue = [ rvalue,  255,   255,   255,   255,   255,   245,   255,   250 ]
           gvalue = [ gvalue,  255,   250,   255,   255,   248,   245,   245,   240 ]
           bvalue = [ bvalue,  255,   250,   240,   224,   220,   220,   238,   230 ]
           colors = [ colors, 'TAN1', 'TAN2', 'TAN3', 'TAN4', 'TAN5', 'TAN6', 'TAN7', 'TAN8']
           rvalue = [ rvalue,   250,   255,    255,    255,    255,    245,    222,    210 ]
           gvalue = [ gvalue,   235,   239,    235,    228,    228,    222,    184,    180 ]
           bvalue = [ bvalue,   215,   213,    205,    196,    181,    179,    135,    140 ]
           colors = [ colors, 'BLK1', 'BLK2', 'BLK3', 'BLK4', 'BLK5', 'BLK6', 'BLK7', 'BLK8']
           rvalue = [ rvalue,   250,   230,    210,    190,    128,     110,    70,       0 ]
           gvalue = [ gvalue,   250,   230,    210,    190,    128,     110,    70,       0 ]
           bvalue = [ bvalue,   250,   230,    210,    190,    128,     110,    70,       0 ]
           colors = [ colors, 'GRN1', 'GRN2', 'GRN3', 'GRN4', 'GRN5', 'GRN6', 'GRN7', 'GRN8']
           rvalue = [ rvalue,   250,   223,    173,    109,     53,     35,      0,       0 ]
           gvalue = [ gvalue,   253,   242,    221,    193,    156,     132,    97,      69 ]
           bvalue = [ bvalue,   202,   167,    142,    115,     83,      67,    52,      41 ]
           colors = [ colors, 'BLU1', 'BLU2', 'BLU3', 'BLU4', 'BLU5', 'BLU6', 'BLU7', 'BLU8']
           rvalue = [ rvalue,   232,   202,    158,     99,     53,     33,      8,       8 ]
           gvalue = [ gvalue,   241,   222,    202,    168,    133,    113,     75,      48 ]
           bvalue = [ bvalue,   250,   240,    225,    211,    191,    181,    147,     107 ]
           colors = [ colors, 'ORG1', 'ORG2', 'ORG3', 'ORG4', 'ORG5', 'ORG6', 'ORG7', 'ORG8']
           rvalue = [ rvalue,   254,    253,    253,    250,    231,    217,    159,    127 ]
           gvalue = [ gvalue,   236,    212,    174,    134,     92,     72,     51,     39 ]
           bvalue = [ bvalue,   217,    171,    107,     52,     12,      1,      3,      4 ]
           colors = [ colors, 'RED1', 'RED2', 'RED3', 'RED4', 'RED5', 'RED6', 'RED7', 'RED8']
           rvalue = [ rvalue,   254,    252,    252,    248,    225,    203,    154,    103 ]
           gvalue = [ gvalue,   232,    194,    146,     97,     45,     24,     12,      0 ]
           bvalue = [ bvalue,   222,    171,    114,     68,     38,     29,     19,     13 ]
           colors = [ colors, 'PUR1', 'PUR2', 'PUR3', 'PUR4', 'PUR5', 'PUR6', 'PUR7', 'PUR8']
           rvalue = [ rvalue,   244,    222,    188,    152,    119,    106,     80,     63 ]
           gvalue = [ gvalue,   242,    221,    189,    148,    108,     82,     32,      0 ]
           bvalue = [ bvalue,   248,    237,    220,    197,    177,    163,    139,    125 ]
           colors = [ colors, 'PBG1', 'PBG2', 'PBG3', 'PBG4', 'PBG5', 'PBG6', 'PBG7', 'PBG8']
           rvalue = [ rvalue,   243,    213,    166,     94,     34,      3,      1,      1 ]
           gvalue = [ gvalue,   234,    212,    189,    164,    138,    129,    101,     70 ]
           bvalue = [ bvalue,   244,    232,    219,    204,    171,    139,     82,     54 ]
           colors = [ colors, 'YGB1', 'YGB2', 'YGB3', 'YGB4', 'YGB5', 'YGB6', 'YGB7', 'YGB8']
           rvalue = [ rvalue,   244,    206,    127,     58,     30,     33,     32,      8 ]
           gvalue = [ gvalue,   250,    236,    205,    175,    125,     95,     48,     29 ]
           bvalue = [ bvalue,   193,    179,    186,    195,    182,    168,    137,     88 ]
           colors = [ colors, 'RYB1', 'RYB2', 'RYB3', 'RYB4', 'RYB5', 'RYB6', 'RYB7', 'RYB8']
           rvalue = [ rvalue,   201,    245,    253,    251,    228,    193,    114,     59 ]
           gvalue = [ gvalue,    35,    121,    206,    253,    244,    228,    171,     85 ]
           bvalue = [ bvalue,    38,    72,     127,    197,    239,    239,    207,    164 ]
           colors = [ colors, 'TG1', 'TG2', 'TG3', 'TG4', 'TG5', 'TG6', 'TG7', 'TG8']
           rvalue = [ rvalue,  84,    163,   197,   220,   105,    51,    13,     0 ]
           gvalue = [ gvalue,  48,    103,   141,   188,   188,   149,   113,    81 ]
           bvalue = [ bvalue,   5,     26,    60,   118,   177,   141,   105,    71 ]
           colors = [ colors, 'OPPOSITE', 'BACKGROUND']
           rvalue = [ rvalue,  opixel[0],  bgcolor[0]]
           gvalue = [ gvalue,  opixel[1],  bgcolor[1]]
           bvalue = [ bvalue,  opixel[2],  bgcolor[2]]
         ENDELSE
   ENDELSE
   
    ; If you have a USERDEF color (from a color triple) then load it here.
    IF N_Elements(usercolor) NE 0 THEN BEGIN
       colors = [colors, 'USERDEF']
       rvalue = [rvalue, usercolor[0]]
       gvalue = [gvalue, usercolor[1]]
       bvalue = [bvalue, usercolor[2]]
    ENDIF
       
    ; Load the colors from the current color table, if you need them.
    IF useCurrentColors THEN BEGIN
        IF (!D.Name NE 'NULL') THEN TVLCT, rrr, ggg, bbb, /GET
        IF decomposedState EQ 0 THEN BEGIN
            colors = SIndgen(256)
            rvalue = rrr
            gvalue = ggg
            bvalue = bbb           
        ENDIF ELSE BEGIN
            colors = [colors, SIndgen(256)]
            rvalue = [rvalue, rrr]
            gvalue = [gvalue, ggg]
            bvalue = [bvalue, bbb]
        ENDELSE
    ENDIF
    
    ; Make sure we are looking at compressed, uppercase names.
    colors = StrUpCase(StrCompress(StrTrim(colors,2), /Remove_All))

    ; Check synonyms of color names.
    FOR j=0, N_Elements(theColor)-1 DO BEGIN
       IF StrUpCase(theColor[j]) EQ 'GREY' THEN theColor[j] = 'GRAY'
       IF StrUpCase(theColor[j]) EQ 'LIGHTGREY' THEN theColor[j] = 'LIGHTGRAY'
       IF StrUpCase(theColor[j]) EQ 'MEDIUMGREY' THEN theColor[j] = 'MEDIUMGRAY'
       IF StrUpCase(theColor[j]) EQ 'SLATEGREY' THEN theColor[j] = 'SLATEGRAY'
       IF StrUpCase(theColor[j]) EQ 'DARKGREY' THEN theColor[j] = 'DARKGRAY'
       IF StrUpCase(theColor[j]) EQ 'AQUA' THEN theColor[j] = 'AQUAMARINE'
       IF StrUpCase(theColor[j]) EQ 'SKY' THEN theColor[j] = 'SKYBLUE'
       IF StrUpCase(theColor[j]) EQ 'NAVYBLUE' THEN theColor[j] = 'NAVY'
       IF StrUpCase(theColor[j]) EQ 'CORNFLOWER' THEN theColor[j] = 'CORNFLOWERBLUE'
       IF StrUpCase(theColor[j]) EQ 'BROWN' THEN theColor[j] = 'SIENNA'
    ENDFOR
    
    ; How many colors do we have?
    ncolors = N_Elements(colors)
    
    ; Check for offset.
    IF (theDepth EQ 8) OR (decomposedState EQ 0) THEN offset = !D.Table_Size - ncolors - 2 ELSE offset = 0
    IF (useCurrentColors) AND (decomposedState EQ 0) THEN offset = 0
        
    ; Did the user want to select a color name? If so, we set
    ; the color name and color index, unless the user provided
    ; them. In the case of a single positional parameter, we treat
    ; this as the color index number as long as it is not a string.
    cancelled = 0.0
    IF Keyword_Set(selectcolor) THEN BEGIN
    
       CASE N_Params() OF
          0: BEGIN
             theColor = cgPickColorName(Filename=filename, _Strict_Extra=extra, Cancel=cancelled, BREWER=brewer)
             IF cancelled THEN RETURN, !P.Color
             IF theDepth GT 8 AND (decomposedState EQ 1) THEN BEGIN
                   colorIndex = Fix(!P.Color < (!D.Table_Size - 1))
             ENDIF ELSE BEGIN
                   colorIndex = Where(StrUpCase(colors) EQ StrUpCase(StrCompress(theColor, /Remove_All)), count) + offset
                   colorIndex = Fix(colorIndex[0])
                   IF count EQ 0 THEN Message, 'Cannot find color: ' + StrUpCase(theColor), /NoName
             ENDELSE
    
             END
          1: BEGIN
             IF Size(theColor, /TName) NE 'STRING' THEN BEGIN
                colorIndex = Fix(theColor)
                theColor = brewer ? 'WT1' : 'White'
             ENDIF ELSE colorIndex = Fix(!P.Color < 255)
             theColor = cgPickColorName(theColor, Filename=filename, _Strict_Extra=extra, Cancel=cancelled, BREWER=brewer)
             IF cancelled THEN RETURN, !P.Color
             END
          2: BEGIN
             theColor = cgPickColorName(theColor, Filename=filename, _Strict_Extra=extra, Cancel=cancelled, BREWER=brewer)
             IF cancelled THEN RETURN, !P.Color
             END
       ENDCASE
    ENDIF
    
    ; Make sure you have a color name and color index.
    CASE N_Elements(theColor) OF
       0: BEGIN
             theColor = brewer ? 'WT1' : 'White'
             IF N_Elements(colorIndex) EQ 0 THEN BEGIN
                IF theDepth GT 8 THEN BEGIN
                   colorIndex = Fix(!P.Color < (!D.Table_Size - 1))
                ENDIF ELSE BEGIN
                   colorIndex = Where(colors EQ theColor, count) + offset
                   colorIndex = Fix(colorIndex[0])
                   IF count EQ 0 THEN Message, 'Cannot find color: ' + theColor, /NoName
                ENDELSE
             ENDIF ELSE colorIndex = 0S > colorIndex < Fix((!D.Table_Size - 1))
          ENDCASE
    
       1: BEGIN
             type = Size(theColor, /TNAME)
             IF type NE 'STRING' THEN Message, 'The color must be expressed as a color name.'
             theColor = theColor[0] ; Make it a scalar or you run into a WHERE function "feature". :-(
             IF N_Elements(colorIndex) EQ 0 THEN BEGIN
                IF (theDepth GT 8) AND (decomposedState EQ 1) THEN BEGIN
                   colorIndex = Fix(!P.Color < (!D.Table_Size - 1))
                ENDIF ELSE BEGIN
                   colorIndex = Where(colors EQ theColor, count) + offset
                   colorIndex = Fix(colorIndex[0])
                   IF count EQ 0 THEN Message, 'Cannot find color: ' + theColor, /NoName
                ENDELSE
             ENDIF ELSE colorIndex = 0S > colorIndex < Fix(!D.Table_Size - 1)
             ENDCASE
    
       ELSE: BEGIN
             type = Size(theColor, /TNAME)
             IF type NE 'STRING' THEN Message, 'The colors must be expressed as color names.'
             ncolors = N_Elements(theColor)
             CASE N_Elements(colorIndex) OF
                0: colorIndex = Fix(Indgen(ncolors) + (!D.Table_Size - (ncolors + 1)))
                1: colorIndex = Fix(Indgen(ncolors) + colorIndex)
                ELSE: IF N_Elements(colorIndex) NE ncolors THEN $
                   Message, 'Index vector must be the same length as color name vector.'
             ENDCASE
    
                ; Did the user want color triples?
    
             IF Keyword_Set(triple) THEN BEGIN
                colors = BytArr(ncolors, 3)
                FOR j=0,ncolors-1 DO colors[j,*] = cgColor(theColor[j], colorIndex[j], Filename=filename, $
                   Decomposed=decomposedState, /Triple, BREWER=brewer)
                IF Keyword_Set(row) THEN RETURN, Transpose(Byte(colors)) ELSE RETURN, Byte(colors)
             ENDIF ELSE BEGIN
                colors = LonArr(ncolors)
                FOR j=0,ncolors-1 DO colors[j] = cgColor(theColor[j], colorIndex[j], Filename=filename, $
                   Decomposed=decomposedState, BREWER=brewer)
                IF decomposedState THEN RETURN, colors ELSE RETURN, Byte(colors)
            ENDELSE
          END
    ENDCASE
    
    ; Did the user ask for the color names? If so, return them now.
    IF Keyword_Set(names) THEN RETURN, Reform(colors, 1, ncolors)
    
    ; Process the color names.
    theNames = StrUpCase( StrCompress(colors, /Remove_All ) )
    
    ; Find the asked-for color in the color names array.
    theIndex = Where(theNames EQ StrUpCase(StrCompress(theColor, /Remove_All)), foundIt)
    theIndex = theIndex[0]
    
    ; If the color can't be found, report it and continue with the color set to "OPPOSITE."
    IF foundIt EQ 0 THEN BEGIN
       Message, "Can't find color " + theColor + ". Substituting 'OPPOSITE'.", /Informational
       theColor = 'OPPOSITE'
       theIndex = Where(StrUpCase(colors) EQ 'OPPOSITE')
    ENDIF
    
    ; Get the color triple for this color.
    r = rvalue[theIndex]
    g = gvalue[theIndex]
    b = bvalue[theIndex]
    
    ; Did the user want a color triple? If so, return it now.
    IF Keyword_Set(triple) THEN BEGIN
       IF Keyword_Set(allcolors) THEN BEGIN
          IF Keyword_Set(row) $
             THEN RETURN, Byte(Transpose([[rvalue], [gvalue], [bvalue]])) $
             ELSE RETURN, Byte([[rvalue], [gvalue], [bvalue]])
       ENDIF ELSE BEGIN
          IF Keyword_Set(row) THEN RETURN, Byte([r, g, b]) ELSE RETURN, Byte([[r], [g], [b]])
       ENDELSE
    ENDIF
    
    ; Otherwise, we are going to return either an index
    ; number where the color has been loaded, or a 24-bit
    ; value that can be decomposed into the proper color.
    IF decomposedState THEN BEGIN
    
       ; Need a color structure?
       IF Arg_Present(colorStructure) THEN BEGIN
          theColors = cgColor24([[rvalue], [gvalue], [bvalue]])
          colorStructure = Create_Struct(theNames[0], theColors[0])
          FOR j=1, ncolors-1 DO colorStructure = Create_Struct(colorStructure, theNames[j], theColors[j])
       ENDIF
    
       IF Keyword_Set(allcolors) THEN BEGIN
          RETURN, cgColor24([[rvalue], [gvalue], [bvalue]])
       ENDIF ELSE BEGIN
          RETURN, cgColor24([r, g, b])
       ENDELSE
    
    ENDIF ELSE BEGIN
    
       IF Keyword_Set(allcolors) THEN BEGIN
    
          ; Need a color structure?
          IF Arg_Present(colorStructure) THEN BEGIN
             allcolorIndex = !D.Table_Size - ncolors - 2
             IF allcolorIndex LT 0 THEN $
                Message, 'Number of colors exceeds available color table values. Returning.', /NoName
             IF (allcolorIndex + ncolors) GT 255 THEN $
                Message, 'Number of colors exceeds available color table indices. Returning.', /NoName
             theColors = IndGen(ncolors) + allcolorIndex
             colorStructure = Create_Struct(theNames[0],  theColors[0])
             FOR j=1, ncolors-1 DO colorStructure = Create_Struct(colorStructure, theNames[j], theColors[j])
          ENDIF
    
          IF N_Elements(colorIndex) EQ 0 THEN colorIndex = Fix(!D.Table_Size - ncolors - 2)
          IF colorIndex LT 0 THEN $
             Message, 'Number of colors exceeds available color table values. Returning.', /NoName
          IF (colorIndex + ncolors) GT 255 THEN BEGIN
             colorIndex = Fix(!D.Table_Size - ncolors - 2)
          ENDIF
          IF (!D.Name NE 'PRINTER') AND (!D.Name NE 'NULL') THEN TVLCT, rvalue, gvalue, bvalue, colorIndex
          RETURN, BIndGen(ncolors) + colorIndex
       ENDIF ELSE BEGIN
    
          ; Need a color structure?
          IF Arg_Present(colorStructure) THEN BEGIN
             colorStructure = Create_Struct(theColor,  colorIndex)
          ENDIF
    
          IF (!D.Name NE 'PRINTER') AND (!D.Name NE 'NULL') THEN $
              TVLCT, rvalue[theIndex], gvalue[theIndex], bvalue[theIndex], colorIndex
          RETURN, Byte(colorIndex)
       ENDELSE
    
    
    ENDELSE

END ;-------------------------------------------------------------------------------------------------------
;;==============================================================================
; docformat = 'rst'
;
; NAME:
;   cgColor24
;
; PURPOSE:
;   The purpose of this function is to convert a RGB color triple
;   into the equivalent 24-bit long integer. The 24-bit integer
;   can be decomposed into the appropriate color by interpreting
;   the lowest 8 bits as red, the middle 8 bits as green, and the
;   highest 8 bits as blue.
;
;******************************************************************************************;
;                                                                                          ;
;  Copyright (c) 2012, by Fanning Software Consulting, Inc. All rights reserved.           ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;
;+
; The purpose of this function is to convert a RGB color triple
; into the equivalent 24-bit long integer. The 24-bit integer
; can be decomposed into the appropriate color by interpreting
; the lowest 8 bits as red, the middle 8 bits as green, and the
; highest 8 bits as blue. This routine was written to be used with 
; device-independent color programs like `cgColor`.
; 
; :Categories:
;    Graphics, Utilities
;    
; :Returns:
;    A 24-bit long integer that can be decomposed into a color triple value.
;    
; :Params:
;    color: in, required
;       A three-element column or row array representing a color triple. Or an 
;       N-by-three element array of color triples. The values of the elements 
;       must be between 0 and 255.
;
; :Examples:
;    To convert the color triple for the color YELLOW, (255, 255, 0), to the 
;    hexadecimal value '00FFFF'x or the decimal number 65535, type::
;    
;       color = COLOR24([255, 255, 0])
;
; :Author:
;    FANNING SOFTWARE CONSULTING::
;       David W. Fanning 
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: david@idlcoyote.com
;       Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Written by:  David Fanning, 3 February 96.
;        Completely revised the algorithm to accept color arrays. 19 October 2000. DWF.
;            
; :Copyright:
;     Copyright (c) 1996-2012, Fanning Software Consulting, Inc.
;-
FUNCTION cgCOLOR24, color

    ON_ERROR, 2
    
    s = Size(color)
    
    IF s[0] EQ 1 THEN BEGIN
       IF s[1] NE 3 THEN Message, 'Input color parameter must be a 3-element vector.'
       RETURN, color[0] + (color[1] * 2L^8) + (color[2] * 2L^16)
    ENDIF ELSE BEGIN
       IF s[2] GT 3 THEN Message, 'Input color parameter must be an N-by-3 array.'
       RETURN, color[*,0] + (color[*,1] * 2L^8) + (color[*,2] * 2L^16)
    ENDELSE

END ;--------------------------------------------------------------------------------------------
;;==============================================================================
; docformat = 'rst'
;
; NAME:
;   cgDefaultColor
;
; PURPOSE:
;   The purpose of this function is to choose a default color for Coyote Graphics routines.
;
;******************************************************************************************;
;                                                                                          ;
;  Copyright (c) 2011, by Fanning Software Consulting, Inc. All rights reserved.           ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;
;+
; The purpose of this function is to choose a default color for Coyote Graphics routines.
; 
; :Categories:
;    Graphics
;    
; :Returns:
;    Returns a scalar or vector (depends on the type of the input color) of color 
;    names (strings) to be used as the "color" in Coyote Graphics routines. If the
;    MODE is 1 and the inputColor is of type LONG, then the return value is an array
;    or scalar of type LONG.
; 
; :Params:
;    inputcolour: in, optional
;        The input color. May be undefined, a byte, integer, long, or string. If the
;        device is in indexed color mode at the time the request is made, then all,
;        byte, integer, and long values will be treated as indices into the current
;        color table. The variable may be a vector.
;        
; :Keywords:
;     background: in, optional, type=boolean
;        If this keyword is set, the color is treated as a background color. Otherwise,
;        it is treated as a drawing color.
;     default: in, optional
;         A color of any type allowed as the `inputColour`. Used if the `inputColour` is undefined.
;     mode: in, optional, type=boolean
;         The color mode. A 0 mean indexed color mode. A 1 means decomposed color mode.
;         If not supplied in the call, the color mode is determined at run-time with `GetDecomposedState`.
;         The color mode is *always* determined at run time if the current graphics device 
;         is the PostScript device.
;     traditional: in, optional, type=boolean, default=0
;         Set this keyword if you are using the traditional color scheme of white foreground
;         and black background. If this keyword is set, and the current graphics device is
;         the PostScript device, the colors will be reversed, in the traditional IDL graphics
;         way.
;         
; :Examples:
;    Use as a device independent way to get a color::
;       background = cgDefaultColor(bColor, /Background)
;       color = cgDefaultColor(bColor)
;       cgPlot, cgDemoData, Background=background, Color=color
;       
;       
; :Author:
;    FANNING SOFTWARE CONSULTING::
;        David W. Fanning 
;        1645 Sheely Drive
;        Fort Collins, CO 80526 USA
;        Phone: 970-221-0438
;        E-mail: david@idlcoyote.com
;        Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Written, 24 December 2011. David W. Fanning.
;        Modified to make sure a LONG integer in indexed color mode is in the range 0-255. 10 Jan 2012. DWF.
;        Modified to make sure MODE is always determined at run-time for PostScript device. 14 Jan 2012. DWF.
;        Allow other data types to be treated as color table index numbers, as long as they are in the
;           range 0 to 255, and the MODE indicates indexed color. 7 March 2012. DWF.
;        Modified so that the variable MODE will not change in the calling program program. 8 March 2012. DWF.
;        Made FOR loop counter a LONG integer. 3 July 2012. DWF.
;        More work on getting the MODE correct to handle LONG integers properly. Now, in PostScript, we get
;            the current mode, except if we have been told what mode to be in. 21 July 2012. DWF.
;        For numerical values less than 256, in indexed color state, I now return the values
;              directly to the user. This should significantly speed up many Coyote Graphics
;              processes. 14 December 2012. DWF.
;        Modified to return byte and integer values as LONG decomposed integers, if decomposed
;              mode is currently in effect. 14 December 2012. DWF.
;        Now setting all background colors to WHITE, not BACKGROUND. 16 Dec 2012. DWF.
;        Whoops! Forgot one of the background colors in the Traditional section. 11 Jan 2013. DWF.
;        Made the counters in the loops long integers to accommodate large color vectors. 29 Apr 2013. Joe Sapp.
;        Better error handling if/when long integers are used in indexed color mode. 17 Sept 2013. DWF.
;        
; :Copyright:
;     Copyright (c) 2011, Fanning Software Consulting, Inc.
;-
FUNCTION cgDefaultColor, inputColour, $
    BACKGROUND=background, $
    DEFAULT=default, $
    MODE=mode, $
    TRADITIONAL=traditional
    
    ; Return to the caller on error.
    On_Error, 2
    
    ; Default values and variables needed for the program.
    background = Keyword_Set(background)
    IF N_Elements(inputColour) NE 0 THEN inputColor = inputColour
    
    ; Getting the mode correct is critical for handling LONG integers correctly.
    ; If we are doing this in PostScript, then we have to get the mode while we
    ; are in the PostScript file (since we try to draw in DECOMPOSED mode always).
    ; But, if we are TOLD what mode to use, we have to honor that, even in PostScript.
    IF ((N_Elements(mode) EQ 0) || (!D.Name EQ 'PS')) THEN BEGIN
       thisMode = GetDecomposedState()
       IF (N_Elements(mode) NE 0) THEN thisMode = Keyword_Set(mode)
    ENDIF ELSE thisMode = Keyword_Set(mode)
    traditional = Keyword_Set(traditional)
    thisDevice = !D.Name
    
    ; Is the color undefined? If so, is there a default color to assign?
    IF N_Elements(inputColor) EQ 0 THEN BEGIN
        IF N_Elements(default) NE 0 THEN inputColor = default
    ENDIF
    
    ; Is the input color still undefined?
    IF (N_Elements(inputColor) EQ 0) THEN BEGIN
       IF background THEN BEGIN
       
          IF traditional THEN BEGIN
             CASE thisDevice OF
                 'PS': inputColor = 'WHITE'
                 'Z': inputColor = 'BLACK'
                 ELSE: inputColor = 'BLACK'
             ENDCASE
          ENDIF ELSE BEGIN
              CASE thisDevice OF
                 'PS': inputColor = 'WHITE'
                 'Z': inputColor = 'WHITE'
                 ELSE: inputColor = 'WHITE'
             ENDCASE
          ENDELSE
       
       ENDIF ELSE BEGIN
       
          IF traditional THEN BEGIN
              CASE thisDevice OF
                 'PS': inputColor = 'BLACK'
                 'Z': inputColor = 'WHITE'
                 ELSE: inputColor = 'WHITE'
             ENDCASE
          ENDIF ELSE BEGIN
              CASE thisDevice OF
                 'PS': inputColor = 'BLACK'
                 'Z': inputColor = 'BLACK'
                 ELSE: inputColor = 'BLACK'
             ENDCASE
          ENDELSE

       ENDELSE
       
    ENDIF
    
    ; If we get here, the input color is defined as *something*.
    ; This is the crux of the problem. What does this color mean?
    ; If it is a byte or integer value, we assume this is an index
    ; into the current color table. We do the same thing with a long
    ; integer if the MODE is 0, or indexed color, and the value is 
    ; between 0 and 255. If it is a string, we can return it directly.
    thisType = Size(inputcolor, /TNAME)
    
    TVLCT, r, g, b, /GET
    IF (thisType EQ 'BYTE') && (thisMode EQ 0) THEN BEGIN
        RETURN, inputColor 
    ENDIF 
    IF (thisType EQ 'BYTE') && (thisMode EQ 1) THEN BEGIN
        ncolors = N_Elements(inputColor)
        colors = LonArr(ncolors)
        FOR j=0L,ncolors-1 DO colors[j] = cgColor24([r[inputColor[j]], g[inputColor[j]], b[inputColor[j]]]) 
        IF N_Elements(colors) EQ 1 THEN RETURN, colors[0] ELSE RETURN, colors
    ENDIF 
    IF (thisType EQ 'INT') && (thisMode EQ 0) THEN BEGIN
        RETURN, inputColor 
    ENDIF 
    IF (thisType EQ 'INT') && (thisMode EQ 1) THEN BEGIN
        ncolors = N_Elements(inputColor)
        colors = LonArr(ncolors)
        
        ; Make sure the input color is in the range 0 to 255
        index = Where((inputcolor GT 255) OR (inputcolor LT 0), count)
        IF count GT 0 THEN BEGIN
            Message, 'Improper input color. It is possible 24-bit colors (LONGs) are being used in indexed color mode to specify colors.'
        ENDIF
        FOR j=0L,ncolors-1 DO colors[j] = cgColor24([r[inputColor[j]], g[inputColor[j]], b[inputColor[j]]]) 
        IF N_Elements(colors) EQ 1 THEN RETURN, colors[0] ELSE RETURN, colors
    ENDIF 
    IF thisType EQ 'LONG' && (thisMode EQ 0) THEN BEGIN
        RETURN, inputColor 
    ENDIF
    IF thisType EQ 'LONG' && (thisMode EQ 1) THEN BEGIN
        theColors = LonArr(N_Elements(inputColor))    
    ENDIF 
    IF N_Elements(theColors) EQ 0 THEN theColors = StrArr(N_Elements(inputColor))

    ; Fill the color return array.
    FOR j=0L,N_Elements(theColors)-1 DO BEGIN
        thisColor = inputColor[j]
        CASE thisType OF
        
            'STRING': theColors[j] = StrTrim(thisColor,2)
            'LONG': BEGIN
                IF (thisMode EQ 0) && ( (thisColor GE 0) && (thisColor LT 256) ) THEN BEGIN
                    theColors[j] = StrTrim(thisColor, 2)
                ENDIF ELSE BEGIN
                   IF (thisMode EQ 0) && ( (thisColor LT 0) || (thisColor GT 255) ) THEN BEGIN
                       Message, 'Value of LONG integer ' + StrTrim(thisColor,2) + ' is out of indexed color range.'
                   ENDIF ELSE BEGIN
                      theColors[j] = thisColor
                   ENDELSE
                ENDELSE
                END
            ELSE: BEGIN
              ; I feel like an enabler of bad programming practices. Sigh...
              IF (thisColor GE 0) && (thisColor LE 255) && (thisMode EQ 0) THEN BEGIN
                  theColors[j] = StrTrim(Fix(thisColor),2)
              ENDIF ELSE BEGIN
                  Message, 'Cannot determine a color from a value of ' + $
                            StrTrim(thisColor,2) + ' of data type ' + thisType + '.'
              ENDELSE
              END
        ENDCASE
        
    ENDFOR
    
    ; Return the colors after making sure to return a scalar if there is only one element.
    IF N_Elements(theColors) EQ 1 THEN theColors = theColors[0]
    RETURN, theColors
    
END
;;==============================================================================
; docformat = 'rst'
;
; NAME:
;   cgDisplay
;
; PURPOSE:
;   The purpose of cgDisplay is to open a graphics window on the display, or in the
;   PostScript device, or in the Z-graphics buffer, depending upon the current graphics
;   device. In PostScript a window of the proper aspect ratio is created with PSWindow.
;   Using cgDisplay to open "windows" will allow you to more easily write device-independent
;   IDL programs.
;
;******************************************************************************************;
;                                                                                          ;
;  Copyright (c) 2010, by Fanning Software Consulting, Inc. All rights reserved.           ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;
;+
;   The purpose of cgDisplay is to open a graphics window on the display, or in the
;   PostScript device, or in the Z-graphics buffer, depending upon the current graphics
;   device. In PostScript a window of the proper aspect ratio is created with PSWindow.
;   Using cgDisplay to open "windows" will allow you to more easily write device-independent
;   IDL programs.
;
; :Categories:
;    Graphics
;    
; :Params:
;    pxsize: in, optional, type=integer, default=640
;         The X size of the graphics window created. By default, 640.
;    pysize: in, optional, type=integer, default=512
;         The Y size of the graphics window created. By default, 512.
;         
; :Keywords:
;    aspect, in, optional, type=float
;        Set this keyword to create a window with this aspect ratio (ysize/xsize).
;        If aspect is greater than 1, then the ysize will be used in the aspect
;        ratio calculation. If the aspect is less than or equal to 1, then the
;        xsize will be used in the aspect ratio calculation of the final window size.
;        If the input to the ASPECT keyword is an image, then the aspect ratio will
;        be calculated from the image itself.
;    color: in, optional, type=string/integer, default='white'
;        If this keyword is a string, the name of the data color. By default, 'white'.
;        Color names are those used with cgColor. Otherwise, the keyword is assumed 
;        to be a color index into the current color table. The color is not used if
;        the "window" is opened in PostScript on the Z-graphics buffer.
;    force: in, optional, type=boolean, default=0
;         Because of the way cgDisplay is designed to work in many devices and in resizeable
;         graphics windows, it is sometimes the case that it won't create a window for you.
;         If you set this keyword, a graphics window will be created while in any device that 
;         supports graphics windows.
;    free: in, optional, type=boolean, default=0
;         Set this keyword to open a window with a free or unused window index number.
;         This keyword applied only to graphics windows created on the computer display.
;    location: in, optional, type=integer
;         Set this keyword to a two-element integer array indicated the pixel position of
;         the upper-left corner of the graphics window from the upper-left corner of the display.
;    match: in, optional, type=boolean, default=0
;         If this keyword is set, the new display window will match the size of the current
;         display window, if there is one.
;    pixmap: in, optional, type=boolean, default=0
;         Set this keyword to create a pixmap window (a window in memory only).
;    retain: in, optional, type=integer
;         Set this keyword to the values 0, 1, or 2, to indicate no backing store, server
;         provided backing store, or IDL provided backing store, respectively. By default,
;         set to 1 for Windows users and to 2 for UNIX users.
;    title: in, optional, type=string
;         Set this keyword to a string that is used as the window title.
;    wid: in, optional, type=integer, default=0
;         The window index number of the IDL graphics window to create.
;    window: in, optional, type=integer, default=0
;         Because I want to use cgDisplay everywhere, including in resizeable graphics
;         windows, and I don't want it opening windows then, it first checks to be sure
;         there are no resizeable graphics windows on the display before it creates a window.
;         Setting this keyword will overrule this check and create a normal IDL graphics window
;         on the display. This will allow you to open a normal graphics window at the same
;         time a resizeable graphics window exists on the display.
;    xpos: in, optional, type=integer
;         The X position of the window, specified in device coordinates. On Motif platforms, 
;         XPOS specifies the X position of the lower left corner and is measured from the 
;         lower left corner of the screen. On Windows platforms, XPOS specifies the X position 
;         of the upper left corner and is measured from the upper left corner of the screen.
;         This value can also be specified as the first element in the `Location` keyword.
;    xsize: in, optional, type=integer, default=640
;         The X size of the graphics window created. By default, 640. The PXSIZE parameter 
;         is used in preference to the XSIZE keyword value.
;    ypos: in, optional, type=integer
;         The Y position of the window, specified in device coordinates. On Motif platforms, 
;         YPOS specifies the Y position of the lower left corner and is measured from the 
;         lower left corner of the screen. On Windows platforms, YPOS specifies the Y position 
;         of the upper left corner and is measured from the upper left corner of the screen.
;         This value can also be specified as the second element in the `Location` keyword.
;    ysize: in, optional, type=integer, default=512
;         The Y size of the graphics window created. By default, 512. The PYSIZE parameter 
;         is used in preference to the YSIZE keyword value.
;    _extra: in, optional, type=any
;         Any keywords supported by the WINDOW command are allowed.
;         
; :Examples:
;    Use like the IDL WINDOW command::
;       IDL> cgDisplay, XSIZE=500 YSIZE=400
;       IDL> cgDisplay, 500, 500, WID=1, COLOR='gray'
;       
; :Author:
;       FANNING SOFTWARE CONSULTING::
;           David W. Fanning 
;           1645 Sheely Drive
;           Fort Collins, CO 80526 USA
;           Phone: 970-221-0438
;           E-mail: david@idlcoyote.com
;           Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Written, 15 November 2010. DWF.
;        Changes so that color variables don't change type. 23 Nov 2010. DWF.
;        Moved the window index argument to the WID keyword. 9 Dec 2010. DWF.
;        Modified to produce a window in PostScript and the Z-buffer, too. 15 Dec 2010. DWF.
;        Added the FREE keyword. 3 January 2011. DWF.
;        I made a change that allows you to call cgDisplay inside a program that is
;           going to be added to a cgWindow. The program will not open a graphics window
;           if the current graphics window ID is found in a list of cgWindow window IDs.
;           It is now possible to use cgDisplay in any graphics program, even those that
;           will be run in cgWindow. 17 Nov 2011. DWF.
;        Added ASPECT keyword. 18 Nov 2011. DWF.
;        Allowed the window ASPECT to be set with an image argument. 25 Nov 2011. DWF.
;        Now use Scope_Level to always create a display when cgDisplay is called from
;           the main IDL level. 7 Feb 2012. DWF.
;        Added FORCE and MATCH keywords. 16 Feb 2012. DWF.
;        Added PIXMAP, RETAIN, TITLE, XPOS, YPOS, and LOCATION keywords. 4 Sept 2012. DWF.
;
; :Copyright:
;     Copyright (c) 2010-2012, Fanning Software Consulting, Inc.
;-
PRO cgDisplay, pxsize, pysize, $
    ASPECT=aspect, $
    COLOR=scolor, $
    FREE=free, $
    FORCE=force, $
    LOCATION=location, $
    MATCH=match, $
    PIXMAP=pixmap, $
    RETAIN=retain, $
    TITLE=title, $
    WID=windowIndex, $
    WINDOW=window, $
    XPOS=xpos, $
    XSIZE=xsize, $
    YPOS=ypos, $
    YSIZE=ysize, $
    _EXTRA=extra

    Compile_Opt idl2

    ; Error handling.
    Catch, theError
    IF theError NE 0 THEN BEGIN
        Catch, /CANCEL
        void = Error_Message()
        RETURN
    ENDIF
    
    ; Set up PostScript device for working with colors.
    IF !D.Name EQ 'PS' THEN Device, COLOR=1, BITS_PER_PIXEL=8
    
    ; Check parameters and keywords.
    free = Keyword_Set(free)
    IF Keyword_Set(match) THEN BEGIN
       IF !D.Window GE 0 THEN BEGIN
          xsize = !D.X_Size
          ysize = !D.Y_Size
       ENDIF
    ENDIF
    IF N_Elements(location) NE 0 THEN BEGIN
       IF N_Elements(xpos) EQ 0 THEN xpos = location[0]
       IF N_Elements(ypos) EQ 0 THEN ypos = location[1]
    ENDIF
    IF N_Elements(retain) EQ 0 THEN retain = (StrUpCase(!Version.OS_Family) EQ 'UNIX') ? 2 : 1
    IF N_Elements(scolor) EQ 0 THEN color = 'white' ELSE color = scolor
    IF N_Elements(windowIndex) EQ 0 THEN windowIndex = 0
    IF N_Elements(xsize) EQ 0 THEN xsize = 640
    IF N_Elements(ysize) EQ 0 THEN ysize = 512
    IF N_Elements(pxsize) EQ 0 THEN pxsize = xsize
    IF N_Elements(pysize) EQ 0 THEN pysize = ysize
    
    ; Do you need a window with a particular aspect ratio? Can't do this
    ; if you are matching a window.
    IF (N_Elements(aspect) NE 0) && (Keyword_Set(match) EQ 0) THEN BEGIN
    
       ; If aspect is not a scalar, but an image. Use the aspect
       ; ratio of the image to determine the aspect ratio of the
       ; display.
       ndims = Size(aspect, /N_DIMENSIONS)
       IF  (ndims GE 2) && (ndims LE 4) THEN BEGIN 
           void = Image_Dimensions(aspect, XSIZE=xsize, YSIZE=ysize)
           waspect = Float(ysize) / xsize
       ENDIF ELSE waspect = aspect
       
       IF waspect GT 1.0 THEN BEGIN
          pxsize = pysize / waspect
       ENDIF ELSE BEGIN
          pysize = pxsize * waspect
       ENDELSE
       
    ENDIF
    
    ; If you are on a machine that supports windows, you can create a window
    ; if the current graphics window ID cannot be found in the list of cgWindow IDs.
    ; This will allow you to create a window in a program that can still run in
    ; a resizeable graphics window. If you absolutely MUST have a graphics window,
    ; set the window keyword to force a normal IDL graphics window.
    IF (!D.Flags AND 256) NE 0 THEN BEGIN
    
        ; Assume you can create a window.
        createWindow = 1
        
        IF ~Keyword_Set(window) THEN BEGIN
            ; Get the window ids of all cgWindows.
            windowIDs = cgQuery(COUNT=windowCnt)
            IF windowCnt NE 0 THEN BEGIN
                index = Where(windowIDs EQ !D.Window, foundit)
                IF foundit && (Scope_Level() NE 2) THEN createWindow = 0
            ENDIF 
        ENDIF
        
        ; If you are forcing it, then always create a window.
        IF Keyword_Set(force) THEN createWindow = 1
        
        ; If you are not running this program in a cgWindow, feel
        ; free to create a window!
        IF createWindow THEN BEGIN
            Window, windowIndex, XSIZE=pxsize, YSIZE=pysize, PIXMAP=pixmap, $
                 FREE=free, TITLE=title, XPOS=xpos, YPOS=ypos, $
                 RETAIN=retain, _STRICT_EXTRA=extra
            
            ; cgErase will take care of sorting out what kind of 
            ; "color" indicator we are using (string, long, etc.)
            ; so we don't have to worry about that here.
            cgErase, color   
            
        ENDIF
    ENDIF ELSE BEGIN
        CASE !D.Name OF
        
            ; There can be some strange interactions with PS_START if PS_START
            ; is called with no current windows open, and cgDisplay is called with
            ; an aspect ratio that results in a PORTRAIT mode display. This checks
            ; for that problem.
            'PS': BEGIN
                COMMON _$FSC_PS_START_, ps_struct
                keywords = PSWindow(AspectRatio=Float(pysize)/pxsize)
                Device, _Extra=keywords
                IF N_Elements(ps_struct) NE 0 THEN ps_struct.landscape = keywords.landscape
                END
            'Z': Device, Set_Resolution=[pxsize,pysize]
            ELSE:
        ENDCASE
    ENDELSE
END
;;==============================================================================
; docformat = 'rst'
;
; NAME:
;   cgErase
;
; PURPOSE:
;   Provides a device-independent and color-model-independent way to erase a graphics
;   window with a particular color.
;
;******************************************************************************************;
;                                                                                          ;
;  Copyright (c) 2010, by Fanning Software Consulting, Inc. All rights reserved.           ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;
;+
;   Provides a device-independent and color-model-independent way to erase a graphics
;   window with a particular color.
;
; :Categories:
;    Graphics
;    
; :Params:
;    background_color: in, optional, type=string/integer/long, default='white'
;         The color to use in erasing the graphics window. Default is "white."
;         Color names are those used with cgColor.
;       
; :Keywords:
;     color: in, optional, type=string/integer/long, default='white'
;         An alternative way to specify the color to use in erasing the graphics window.
;         Color names are those used with cgColor. This parameter is used in
;         preference to the background_color parameter.
;     layout: in, optional, type=intarr(3)
;         This keyword specifies a grid with a graphics window and determines where the
;         graphic should appear. The syntax of LAYOUT is three numbers: [ncolumns, nrows, location].
;         The grid is determined by the number of columns (ncolumns) by the number of 
;         rows (nrows). The location of the graphic is determined by the third number. The
;         grid numbering starts in the upper left (1) and goes sequentually by column and then
;         by row. If this keyword is used, only this portion of the window is erased.
;     window: in, optional, type=boolean, default=0
;         Set this keyword to erase the current cgWindow application. "Erasing" in
;         this case means removing all the current commands.
;          
; :Examples:
;    Used to "erase" various things::
;       IDL> cgErase
;       IDL> cgErase, 'gray'
;       IDL> cgErase, COLOR='charcoal'
;       
;       IDL> cgPlot, cgDemoData(1), /Window
;       IDL> cgErase, /Window
;       
;       IDL> cgPlot, cgDemoData(17), Layout=[2,2,1]
;       IDL> cgPlot, cgDemoData(17), Layout=[2,2,4]
;       IDL> cgErase, Layout=[2,2,1]
;       
; :Author:
;       FANNING SOFTWARE CONSULTING::
;           David W. Fanning 
;           1645 Sheely Drive
;           Fort Collins, CO 80526 USA
;           Phone: 970-221-0438
;           E-mail: david@idlcoyote.com
;           Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Written, 12 November 2010. DWF.
;        Modified so that input variables are not changed. 18 Nov 2010. DWF.
;        Got my color selection algorithm right. COLOR keyword takes precedence
;          over the parameter. 18 Nov 2010. DWF.
;        Modified to erase in decomposed color, if possible.
;        In some cases, I was turning BYTE values to strings without converting to 
;            INTEGERS first. 30 Dec 2010. DWF.   
;        Added WINDOW keyword. 26 Jan 2011. DWF. 
;        Added LAYOUT keyword. 1 Feb 2011. DWF.    
;        Modified error handler to restore the entry decomposition state if there is an error. 17 March 2011. DWF
;
; :Copyright:
;     Copyright (c) 2010, Fanning Software Consulting, Inc.
;-
PRO cgErase, background_color, COLOR=color, LAYOUT=layout, WINDOW=window

    Catch, theError
    IF theError NE 0 THEN BEGIN
       Catch, /Cancel
       IF N_Elements(currentState) NE 0 THEN SetDecomposedState, currentState
       RETURN
    ENDIF

    ; Are we erasing an cgWindow application?
    IF Keyword_Set(window) THEN BEGIN
    
        currentID = cgQuery(COUNT=wincnt, OBJECTREF=currentObj, /Current)
        IF wincnt GT 0 THEN BEGIN
            currentObj -> DeleteCommand, /All
            currentObj -> ExecuteCommands
        ENDIF
        RETURN
        
    ENDIF

    ; Set up PostScript device for working with colors.
    IF !D.Name EQ 'PS' THEN Device, COLOR=1, BITS_PER_PIXEL=8
    
    ; Get a color for erasing.
    IF N_Elements(background_color) EQ 0 THEN thisColor = 'white' ELSE thisColor = background_color
    IF N_Elements(color) NE 0 THEN thisColor = color 
    IF Size(thisColor, /TYPE) EQ 3 THEN IF GetDecomposedState() EQ 0 THEN thisColor = Byte(thisColor)
    IF Size(thisColor, /TYPE) LE 2 THEN thisColor = StrTrim(Fix(thisColor),2)

    ; Get the current color vectors.
    TVLCT, rr, gg, bb, /Get
    
    ; Do this in decomposed color, if possible.
    SetDecomposedState, 1, CURRENTSTATE=currentState
    
    IF Size(thisColor, /TNAME) EQ 'STRING' THEN thisColor = cgColor(thisColor)
    
    ; Set up the layout, if necessary.
    IF N_Elements(layout) NE 0 THEN BEGIN
       thisMulti = !P.Multi
       totalPlots = layout[0]*layout[1]
       !P.Multi = [0,layout[0], layout[1], 0, 0]
       IF layout[2] EQ 1 THEN BEGIN
            noerase = 1
            !P.Multi[0] = 0
       ENDIF ELSE BEGIN
            !P.Multi[0] = totalPlots - layout[2] + 1
       ENDELSE
       
       ; Draw an invisible plot in this space.
       Plot, [0,1], XStyle=4, YStyle=4, /NoErase
       
       ; Fill the plot area with the color.
       x = !X.Region
       y = !Y.Region
       PolyFill, [x[0], x[0], x[1], x[1]], $
                 [y[0], y[1], y[1], y[0]], /FIll, $
                 Color=thisColor, /Normal
       !P.Multi = thisMulti
    ENDIF ELSE Erase, thisColor

    
    ; Clean up.
    SetDecomposedState, currentState
    IF !D.Name NE 'Z' THEN TVLCT, rr, gg, bb
   
END
;;==============================================================================
;+
; NAME:
;  cgGREEK
;
; PURPOSE:
;
;   This function provides a device-independent way to ask for a Greek letter as
;   a string that can be included, for example, in a plot title. It uses the Greek
;   simplex font (!4) when used with Hershey fonts, and the Symbol font (!9) when
;   used with PostScript or True-Type fonts. Selects the type of Greek character to 
;   return based on value of !P.FONT. Updated now to return the UNICODE values for 
;   Greek characters for those fonts that support them (Macintosh?).
;   
; AUTHOR:
;
;   FANNING SOFTWARE CONSULTING
;   David Fanning, Ph.D.
;   1645 Sheely Drive
;   Fort Collins, CO 80526 USA
;   Phone: 970-221-0438
;   E-mail: david@idlcoyote.com
;   Coyote's Guide to IDL Programming: http://www.idlcoyote.com/
;
; CATEGORY:
;
;   Graphics, Utilities
;
; CALLING SEQUENCE:
;
;   greekString = cgGreek(greekLetter)
;
; RETURN VALUE:
;
;   greekString    A string that represents the Greek letter.
;
; ARGUMENTS:
;
;  greekLetter:    The name of the Greek letter desired. A string. Default: 'alpha'.
;                  Valid string names are the 24 characters of the Greek alphabet.
;                     alpha        nu
;                     beta         xi
;                     gamma        omicron
;                     delta        pi
;                     epsilon      rho
;                     zeta         sigma
;                     eta          tau
;                     theta        upsilon
;                     iota         phi
;                     kappa        chi
;                     lambda       psi
;                     mu           omega
;                    
;                   Note that if the first letter of the name is capitalized, this is
;                   the equivalent of setting the CAPITAL keyword. 
;
; KEYWORDRS:
;
;  CAPTIAL:        If this keyword is set, the captial Greek letter is returned rather 
;                  than the lowercase Greek letter. An alternative way of capitalizing
;                  the letter is to make the first letter of the name an uppercase letter.
;                  
;  EXAMPLE:        If this keyword is set, the names of the Greek characters and their
;                  symbols are written out in the current graphics window.
;                  
;  PS:             Normally, the PostScript version of the greek letter is returned if
;                  the current device is PostScript and !P.Font is 0 or 1. But, the 
;                  PostScript version of the greek letter can be obtained at any time
;                  and in any device, by setting this keyword.
;                                    
;  UNICODE:        If this keyword is set, the function returns the Unicode for the Greek
;                  letter.
;
; EXAMPLE:
;
;  Lowercase PSI:
;  
;     IDL> Plot, findgen(11), XTitle='This title contains ' + $
;           cgGreek('psi') + ' as a Greek letter' 
;
;  Uppercase PSI:
;  
;     IDL> Plot, findgen(11), XTitle='This title contains ' + $
;           cgGreek('Psi') + ' as a Greek letter' 
; NOTES:
; 
;  See the following article for additional information: 
;  
;       http://www.idlcoyote.com/ps_tips/greeksym.html
;       
; RESTRICTIONS:
; 
;  For this program to work correctly on your graphics display, you should be using
;  Hershey fonts (!P.Font=-1). It will work correctly in PostScript with either 
;  hardware fonts (!P.Font=0) or True-Type fonts (!P.Font=1).
;  
; MODIFICATION HISTORY:
;
;  Written by: David W. Fanning, 9 January 2010.
;  An alternative way to get an uppercase letter is to make the first letter of
;     the Greek name uppercase. (Maarten Sneep's suggestion!) 11 Jan 2010. DWF
;  Had the wrong value for the PostScript version of Phi. 26 January 2010. DWF
;  Added UNICODE keyword and values for Greek characters. 11 June 2010. DWF.
;  Changed the branching from !D.NAME EQ 'PS' to !P.FONT NE -1. (This is actually
;      what the documentation says, and what I intended.) 13 Dec 2010. DWF.
;  I don't think the last change did quite want I wanted. More tweaking to make
;      this more responsive to being in a PostScript file. 31 July 2011. DWF.
;  Added PS keyword so the return value is the PostScript file. This is for
;      convenience only, as the return value will be the PostScript value if
;      the current graphics device is PS and FONT is not equal to -1. 30 Aug 2011. DWF.
;  Retired Greek and replaced by cgGreek. 23 Dec 2012. DWF.
;-
;******************************************************************************************;
;  Copyright (c) 2010, by Fanning Software Consulting, Inc.                                ;
;  All rights reserved.                                                                    ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
Forward_Function cgGreek

PRO cgGreek_Example, UNICODE=unicode, PS=ps

    Compile_Opt hidden
    
    ; The 24 Greek letters.
    letter = [  'alpha', 'beta', 'gamma', 'delta', 'epsilon',  'zeta', $
                'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu', $
                'nu', 'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau', $
                'upsilon', 'phi', 'chi', 'psi', 'omega' ]
    
    ; Output positions.
    x = [0.25, 0.6]
    y = Reverse((Indgen(12) + 1) * (1.0 / 13))
    
    ; Need PostScript output?
    IF Keyword_Set(ps) || (!D.Name EQ 'PS') THEN PS_Start
    
    ; Create a window.
    cgDisplay, 600, 500
    
    ; Output the letters.
    FOR j=0,11 DO BEGIN
        cgText, x[0], y[j], letter[j] + ': ' + $
            cgGreek(letter[j], UNICODE=unicode) + cgGreek(letter[j], /CAPITAL, UNICODE=unicode), $
            /NORMAL, CHARSIZE=1.5
        cgText, x[1], y[j], letter[j+12] + ': ' + $
            cgGreek(letter[j+12], UNICODE=unicode) + cgGreek(letter[j+12], /CAPITAL, UNICODE=unicode), $
            /NORMAL, CHARSIZE=1.5
    ENDFOR
    
    ; Clean up PostScript, if needed.
    IF Keyword_Set(ps) || (!D.Name EQ 'PS') THEN PS_End

    ; Restore the users window.
    IF N_Elements(thisWindow) NE 0 THEN BEGIN
       IF thisWindow GE 0 THEN WSet, thisWindow
    ENDIF
    
END ; --------------------------------------------------------------------------------------


FUNCTION cgGreek, letter, CAPITAL=capital, EXAMPLE=example, PS=PS, UNICODE=unicode

    Compile_Opt idl2
    
    ; Return to caller on error.
    ON_Error, 2
    
    ; Do you wish to see an example?
    IF Keyword_Set(example) THEN BEGIN
        cgGreek_Example, UNICODE=unicode, PS=ps
        RETURN, ""
    ENDIF

    ; A symbol name is required.
    IF N_Elements(letter) EQ 0 THEN BEGIN
       Print, 'Syntax: Greek("theLetter")
       RETURN, ""
    ENDIF
    
    ; Return quietly if the Greek letter is a null string.
    IF letter  EQ "" THEN RETURN, ""
    
     ; Set up PostScript device for working with Greek letters.
    IF !D.Name EQ 'PS' THEN Device, ISOLATIN1=1
    
    ; Check keywords.
    IF N_Elements(letter) EQ 0 THEN letter = 'alpha'
    capital = Keyword_Set(capital)
    
    ; If the first letter of the "letter" variable is uppercase
    ; the user wants an uppercase Greek letter.
    firstLetter = StrMid(letter, 0, 1)
    IF firstLetter EQ StrUpCase(firstLetter) THEN capital = 1
    
    IF Keyword_Set(unicode) THEN BEGIN
    
        CASE StrLowCase(letter) OF
            'alpha':   greekLetter = (capital) ? '!Z(0391)' : '!Z(03B1)'
            'beta':    greekLetter = (capital) ? '!Z(0392)' : '!Z(03B2)'
            'gamma':   greekLetter = (capital) ? '!Z(0393)' : '!Z(03B3)'
            'delta':   greekLetter = (capital) ? '!Z(0394)' : '!Z(03B4)'
            'epsilon': greekLetter = (capital) ? '!Z(0395)' : '!Z(03B5)'
            'zeta':    greekLetter = (capital) ? '!Z(0396)' : '!Z(03B6)'
            'eta':     greekLetter = (capital) ? '!Z(0397)' : '!Z(03B7)'
            'theta':   greekLetter = (capital) ? '!Z(0398)' : '!Z(03B8)'
            'iota':    greekLetter = (capital) ? '!Z(0399)' : '!Z(03B9)'
            'kappa':   greekLetter = (capital) ? '!Z(039A)' : '!Z(03BA)'
            'lambda':  greekLetter = (capital) ? '!Z(039B)' : '!Z(03BB)'
            'mu':      greekLetter = (capital) ? '!Z(039C)' : '!Z(03BC)'
            'nu':      greekLetter = (capital) ? '!Z(039D)' : '!Z(03BD)'
            'xi':      greekLetter = (capital) ? '!Z(039E)' : '!Z(03BE)'
            'omicron': greekLetter = (capital) ? '!Z(039F)' : '!Z(03BF)'
            'pi':      greekLetter = (capital) ? '!Z(03A0)' : '!Z(03C0)'
            'rho':     greekLetter = (capital) ? '!Z(03A1)' : '!Z(03C1)'
            'sigma':   greekLetter = (capital) ? '!Z(03A3)' : '!Z(03C3)'
            'tau':     greekLetter = (capital) ? '!Z(03A4)' : '!Z(03C4)'
            'upsilon': greekLetter = (capital) ? '!Z(03A5)' : '!Z(03C5)'
            'phi':     greekLetter = (capital) ? '!Z(03A6)' : '!Z(03C6)'
            'chi':     greekLetter = (capital) ? '!Z(03A7)' : '!Z(03C7)'
            'psi':     greekLetter = (capital) ? '!Z(03A8)' : '!Z(03C8)'
            'omega':   greekLetter = (capital) ? '!Z(03A9)' : '!Z(03C9)'
            ELSE: Message, 'The greek letter ' + letter + ' is unrecognized.'
       ENDCASE
           
       ; Return the UNICODE greek letter.
       RETURN, greekLetter
       
    ENDIF
    
    IF ((!D.Name EQ 'PS') && (!P.FONT NE -1)) || Keyword_Set(PS)  THEN BEGIN
    
        ; Make sure ISOLATIN1 encoding is turned on.
        IF (!D.Name EQ 'PS') THEN DEVICE, /ISOLATIN1
        
        CASE StrLowCase(letter) OF
            'alpha':   greekLetter = (capital) ? '!9' + String("101B) + '!X' : '!9' + String("141B) + '!X'
            'beta':    greekLetter = (capital) ? '!9' + String("102B) + '!X' : '!9' + String("142B) + '!X'
            'gamma':   greekLetter = (capital) ? '!9' + String("107B) + '!X' : '!9' + String("147B) + '!X'
            'delta':   greekLetter = (capital) ? '!9' + String("104B) + '!X' : '!9' + String("144B) + '!X'
            'epsilon': greekLetter = (capital) ? '!9' + String("105B) + '!X' : '!9' + String("145B) + '!X'
            'zeta':    greekLetter = (capital) ? '!9' + String("132B) + '!X' : '!9' + String("172B) + '!X'
            'eta':     greekLetter = (capital) ? '!9' + String("110B) + '!X' : '!9' + String("150B) + '!X'
            'theta':   greekLetter = (capital) ? '!9' + String("121B) + '!X' : '!9' + String("161B) + '!X'
            'iota':    greekLetter = (capital) ? '!9' + String("111B) + '!X' : '!9' + String("151B) + '!X'
            'kappa':   greekLetter = (capital) ? '!9' + String("113B) + '!X' : '!9' + String("153B) + '!X'
            'lambda':  greekLetter = (capital) ? '!9' + String("114B) + '!X' : '!9' + String("154B) + '!X'
            'mu':      greekLetter = (capital) ? '!9' + String("115B) + '!X' : '!9' + String("155B) + '!X'
            'nu':      greekLetter = (capital) ? '!9' + String("116B) + '!X' : '!9' + String("156B) + '!X'
            'xi':      greekLetter = (capital) ? '!9' + String("130B) + '!X' : '!9' + String("170B) + '!X'
            'omicron': greekLetter = (capital) ? '!9' + String("117B) + '!X' : '!9' + String("157B) + '!X'
            'pi':      greekLetter = (capital) ? '!9' + String("120B) + '!X' : '!9' + String("160B) + '!X'
            'rho':     greekLetter = (capital) ? '!9' + String("122B) + '!X' : '!9' + String("162B) + '!X'
            'sigma':   greekLetter = (capital) ? '!9' + String("123B) + '!X' : '!9' + String("163B) + '!X'
            'tau':     greekLetter = (capital) ? '!9' + String("124B) + '!X' : '!9' + String("164B) + '!X'
            'upsilon': greekLetter = (capital) ? '!9' + String("125B) + '!X' : '!9' + String("165B) + '!X'
            'phi':     greekLetter = (capital) ? '!9' + String("106B) + '!X' : '!9' + String("152B) + '!X'
            'chi':     greekLetter = (capital) ? '!9' + String("103B) + '!X' : '!9' + String("143B) + '!X'
            'psi':     greekLetter = (capital) ? '!9' + String("131B) + '!X' : '!9' + String("171B) + '!X'
            'omega':   greekLetter = (capital) ? '!9' + String("127B) + '!X' : '!9' + String("167B) + '!X'
            ELSE: Message, 'The greek letter ' + letter + ' is unrecognized.'
        ENDCASE
    
    ENDIF ELSE BEGIN
    
        CASE StrLowCase(letter) OF
            'alpha':   greekLetter = (capital) ? '!4' + String("101B) + '!X' : '!4' + String("141B) + '!X'
            'beta':    greekLetter = (capital) ? '!4' + String("102B) + '!X' : '!4' + String("142B) + '!X'
            'gamma':   greekLetter = (capital) ? '!4' + String("103B) + '!X' : '!4' + String("143B) + '!X'
            'delta':   greekLetter = (capital) ? '!4' + String("104B) + '!X' : '!4' + String("144B) + '!X'
            'epsilon': greekLetter = (capital) ? '!4' + String("105B) + '!X' : '!4' + String("145B) + '!X'
            'zeta':    greekLetter = (capital) ? '!4' + String("106B) + '!X' : '!4' + String("146B) + '!X'
            'eta':     greekLetter = (capital) ? '!4' + String("107B) + '!X' : '!4' + String("147B) + '!X'
            'theta':   greekLetter = (capital) ? '!4' + String("110B) + '!X' : '!4' + String("150B) + '!X'
            'iota':    greekLetter = (capital) ? '!4' + String("111B) + '!X' : '!4' + String("151B) + '!X'
            'kappa':   greekLetter = (capital) ? '!4' + String("112B) + '!X' : '!4' + String("152B) + '!X'
            'lambda':  greekLetter = (capital) ? '!4' + String("113B) + '!X' : '!4' + String("153B) + '!X'
            'mu':      greekLetter = (capital) ? '!4' + String("114B) + '!X' : '!4' + String("154B) + '!X'
            'nu':      greekLetter = (capital) ? '!4' + String("115B) + '!X' : '!4' + String("155B) + '!X'
            'xi':      greekLetter = (capital) ? '!4' + String("116B) + '!X' : '!4' + String("156B) + '!X'
            'omicron': greekLetter = (capital) ? '!4' + String("117B) + '!X' : '!4' + String("157B) + '!X'
            'pi':      greekLetter = (capital) ? '!4' + String("120B) + '!X' : '!4' + String("160B) + '!X'
            'rho':     greekLetter = (capital) ? '!4' + String("121B) + '!X' : '!4' + String("161B) + '!X'
            'sigma':   greekLetter = (capital) ? '!4' + String("122B) + '!X' : '!4' + String("162B) + '!X'
            'tau':     greekLetter = (capital) ? '!4' + String("123B) + '!X' : '!4' + String("163B) + '!X'
            'upsilon': greekLetter = (capital) ? '!4' + String("124B) + '!X' : '!4' + String("164B) + '!X'
            'phi':     greekLetter = (capital) ? '!4' + String("125B) + '!X' : '!4' + String("165B) + '!X'
            'chi':     greekLetter = (capital) ? '!4' + String("126B) + '!X' : '!4' + String("166B) + '!X'
            'psi':     greekLetter = (capital) ? '!4' + String("127B) + '!X' : '!4' + String("167B) + '!X'
            'omega':   greekLetter = (capital) ? '!4' + String("130B) + '!X' : '!4' + String("170B) + '!X'
            ELSE: Message, 'The greek letter ' + letter + ' is unrecognized.'      
        ENDCASE
        
      ENDELSE

    RETURN, greekLetter
    
END
;;==============================================================================
;+
; NAME:
;  cgGREEK
;
; PURPOSE:
;
;   This function provides a device-independent way to ask for a Greek letter as
;   a string that can be included, for example, in a plot title. It uses the Greek
;   simplex font (!4) when used with Hershey fonts, and the Symbol font (!9) when
;   used with PostScript or True-Type fonts. Selects the type of Greek character to 
;   return based on value of !P.FONT. Updated now to return the UNICODE values for 
;   Greek characters for those fonts that support them (Macintosh?).
;   
; AUTHOR:
;
;   FANNING SOFTWARE CONSULTING
;   David Fanning, Ph.D.
;   1645 Sheely Drive
;   Fort Collins, CO 80526 USA
;   Phone: 970-221-0438
;   E-mail: david@idlcoyote.com
;   Coyote's Guide to IDL Programming: http://www.idlcoyote.com/
;
; CATEGORY:
;
;   Graphics, Utilities
;
; CALLING SEQUENCE:
;
;   greekString = cgGreek(greekLetter)
;
; RETURN VALUE:
;
;   greekString    A string that represents the Greek letter.
;
; ARGUMENTS:
;
;  greekLetter:    The name of the Greek letter desired. A string. Default: 'alpha'.
;                  Valid string names are the 24 characters of the Greek alphabet.
;                     alpha        nu
;                     beta         xi
;                     gamma        omicron
;                     delta        pi
;                     epsilon      rho
;                     zeta         sigma
;                     eta          tau
;                     theta        upsilon
;                     iota         phi
;                     kappa        chi
;                     lambda       psi
;                     mu           omega
;                    
;                   Note that if the first letter of the name is capitalized, this is
;                   the equivalent of setting the CAPITAL keyword. 
;
; KEYWORDRS:
;
;  CAPTIAL:        If this keyword is set, the captial Greek letter is returned rather 
;                  than the lowercase Greek letter. An alternative way of capitalizing
;                  the letter is to make the first letter of the name an uppercase letter.
;                  
;  EXAMPLE:        If this keyword is set, the names of the Greek characters and their
;                  symbols are written out in the current graphics window.
;                  
;  PS:             Normally, the PostScript version of the greek letter is returned if
;                  the current device is PostScript and !P.Font is 0 or 1. But, the 
;                  PostScript version of the greek letter can be obtained at any time
;                  and in any device, by setting this keyword.
;                                    
;  UNICODE:        If this keyword is set, the function returns the Unicode for the Greek
;                  letter.
;
; EXAMPLE:
;
;  Lowercase PSI:
;  
;     IDL> Plot, findgen(11), XTitle='This title contains ' + $
;           cgGreek('psi') + ' as a Greek letter' 
;
;  Uppercase PSI:
;  
;     IDL> Plot, findgen(11), XTitle='This title contains ' + $
;           cgGreek('Psi') + ' as a Greek letter' 
; NOTES:
; 
;  See the following article for additional information: 
;  
;       http://www.idlcoyote.com/ps_tips/greeksym.html
;       
; RESTRICTIONS:
; 
;  For this program to work correctly on your graphics display, you should be using
;  Hershey fonts (!P.Font=-1). It will work correctly in PostScript with either 
;  hardware fonts (!P.Font=0) or True-Type fonts (!P.Font=1).
;  
; MODIFICATION HISTORY:
;
;  Written by: David W. Fanning, 9 January 2010.
;  An alternative way to get an uppercase letter is to make the first letter of
;     the Greek name uppercase. (Maarten Sneep's suggestion!) 11 Jan 2010. DWF
;  Had the wrong value for the PostScript version of Phi. 26 January 2010. DWF
;  Added UNICODE keyword and values for Greek characters. 11 June 2010. DWF.
;  Changed the branching from !D.NAME EQ 'PS' to !P.FONT NE -1. (This is actually
;      what the documentation says, and what I intended.) 13 Dec 2010. DWF.
;  I don't think the last change did quite want I wanted. More tweaking to make
;      this more responsive to being in a PostScript file. 31 July 2011. DWF.
;  Added PS keyword so the return value is the PostScript file. This is for
;      convenience only, as the return value will be the PostScript value if
;      the current graphics device is PS and FONT is not equal to -1. 30 Aug 2011. DWF.
;  Retired Greek and replaced by cgGreek. 23 Dec 2012. DWF.
;-
;******************************************************************************************;
;  Copyright (c) 2010, by Fanning Software Consulting, Inc.                                ;
;  All rights reserved.                                                                    ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
Forward_Function cgGreek

PRO cgGreek_Example, UNICODE=unicode, PS=ps

    Compile_Opt hidden
    
    ; The 24 Greek letters.
    letter = [  'alpha', 'beta', 'gamma', 'delta', 'epsilon',  'zeta', $
                'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu', $
                'nu', 'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau', $
                'upsilon', 'phi', 'chi', 'psi', 'omega' ]
    
    ; Output positions.
    x = [0.25, 0.6]
    y = Reverse((Indgen(12) + 1) * (1.0 / 13))
    
    ; Need PostScript output?
    IF Keyword_Set(ps) || (!D.Name EQ 'PS') THEN PS_Start
    
    ; Create a window.
    cgDisplay, 600, 500
    
    ; Output the letters.
    FOR j=0,11 DO BEGIN
        cgText, x[0], y[j], letter[j] + ': ' + $
            cgGreek(letter[j], UNICODE=unicode) + cgGreek(letter[j], /CAPITAL, UNICODE=unicode), $
            /NORMAL, CHARSIZE=1.5
        cgText, x[1], y[j], letter[j+12] + ': ' + $
            cgGreek(letter[j+12], UNICODE=unicode) + cgGreek(letter[j+12], /CAPITAL, UNICODE=unicode), $
            /NORMAL, CHARSIZE=1.5
    ENDFOR
    
    ; Clean up PostScript, if needed.
    IF Keyword_Set(ps) || (!D.Name EQ 'PS') THEN PS_End

    ; Restore the users window.
    IF N_Elements(thisWindow) NE 0 THEN BEGIN
       IF thisWindow GE 0 THEN WSet, thisWindow
    ENDIF
    
END ; --------------------------------------------------------------------------------------


FUNCTION cgGreek, letter, CAPITAL=capital, EXAMPLE=example, PS=PS, UNICODE=unicode

    Compile_Opt idl2
    
    ; Return to caller on error.
    ON_Error, 2
    
    ; Do you wish to see an example?
    IF Keyword_Set(example) THEN BEGIN
        cgGreek_Example, UNICODE=unicode, PS=ps
        RETURN, ""
    ENDIF

    ; A symbol name is required.
    IF N_Elements(letter) EQ 0 THEN BEGIN
       Print, 'Syntax: Greek("theLetter")
       RETURN, ""
    ENDIF
    
    ; Return quietly if the Greek letter is a null string.
    IF letter  EQ "" THEN RETURN, ""
    
     ; Set up PostScript device for working with Greek letters.
    IF !D.Name EQ 'PS' THEN Device, ISOLATIN1=1
    
    ; Check keywords.
    IF N_Elements(letter) EQ 0 THEN letter = 'alpha'
    capital = Keyword_Set(capital)
    
    ; If the first letter of the "letter" variable is uppercase
    ; the user wants an uppercase Greek letter.
    firstLetter = StrMid(letter, 0, 1)
    IF firstLetter EQ StrUpCase(firstLetter) THEN capital = 1
    
    IF Keyword_Set(unicode) THEN BEGIN
    
        CASE StrLowCase(letter) OF
            'alpha':   greekLetter = (capital) ? '!Z(0391)' : '!Z(03B1)'
            'beta':    greekLetter = (capital) ? '!Z(0392)' : '!Z(03B2)'
            'gamma':   greekLetter = (capital) ? '!Z(0393)' : '!Z(03B3)'
            'delta':   greekLetter = (capital) ? '!Z(0394)' : '!Z(03B4)'
            'epsilon': greekLetter = (capital) ? '!Z(0395)' : '!Z(03B5)'
            'zeta':    greekLetter = (capital) ? '!Z(0396)' : '!Z(03B6)'
            'eta':     greekLetter = (capital) ? '!Z(0397)' : '!Z(03B7)'
            'theta':   greekLetter = (capital) ? '!Z(0398)' : '!Z(03B8)'
            'iota':    greekLetter = (capital) ? '!Z(0399)' : '!Z(03B9)'
            'kappa':   greekLetter = (capital) ? '!Z(039A)' : '!Z(03BA)'
            'lambda':  greekLetter = (capital) ? '!Z(039B)' : '!Z(03BB)'
            'mu':      greekLetter = (capital) ? '!Z(039C)' : '!Z(03BC)'
            'nu':      greekLetter = (capital) ? '!Z(039D)' : '!Z(03BD)'
            'xi':      greekLetter = (capital) ? '!Z(039E)' : '!Z(03BE)'
            'omicron': greekLetter = (capital) ? '!Z(039F)' : '!Z(03BF)'
            'pi':      greekLetter = (capital) ? '!Z(03A0)' : '!Z(03C0)'
            'rho':     greekLetter = (capital) ? '!Z(03A1)' : '!Z(03C1)'
            'sigma':   greekLetter = (capital) ? '!Z(03A3)' : '!Z(03C3)'
            'tau':     greekLetter = (capital) ? '!Z(03A4)' : '!Z(03C4)'
            'upsilon': greekLetter = (capital) ? '!Z(03A5)' : '!Z(03C5)'
            'phi':     greekLetter = (capital) ? '!Z(03A6)' : '!Z(03C6)'
            'chi':     greekLetter = (capital) ? '!Z(03A7)' : '!Z(03C7)'
            'psi':     greekLetter = (capital) ? '!Z(03A8)' : '!Z(03C8)'
            'omega':   greekLetter = (capital) ? '!Z(03A9)' : '!Z(03C9)'
            ELSE: Message, 'The greek letter ' + letter + ' is unrecognized.'
       ENDCASE
           
       ; Return the UNICODE greek letter.
       RETURN, greekLetter
       
    ENDIF
    
    IF ((!D.Name EQ 'PS') && (!P.FONT NE -1)) || Keyword_Set(PS)  THEN BEGIN
    
        ; Make sure ISOLATIN1 encoding is turned on.
        IF (!D.Name EQ 'PS') THEN DEVICE, /ISOLATIN1
        
        CASE StrLowCase(letter) OF
            'alpha':   greekLetter = (capital) ? '!9' + String("101B) + '!X' : '!9' + String("141B) + '!X'
            'beta':    greekLetter = (capital) ? '!9' + String("102B) + '!X' : '!9' + String("142B) + '!X'
            'gamma':   greekLetter = (capital) ? '!9' + String("107B) + '!X' : '!9' + String("147B) + '!X'
            'delta':   greekLetter = (capital) ? '!9' + String("104B) + '!X' : '!9' + String("144B) + '!X'
            'epsilon': greekLetter = (capital) ? '!9' + String("105B) + '!X' : '!9' + String("145B) + '!X'
            'zeta':    greekLetter = (capital) ? '!9' + String("132B) + '!X' : '!9' + String("172B) + '!X'
            'eta':     greekLetter = (capital) ? '!9' + String("110B) + '!X' : '!9' + String("150B) + '!X'
            'theta':   greekLetter = (capital) ? '!9' + String("121B) + '!X' : '!9' + String("161B) + '!X'
            'iota':    greekLetter = (capital) ? '!9' + String("111B) + '!X' : '!9' + String("151B) + '!X'
            'kappa':   greekLetter = (capital) ? '!9' + String("113B) + '!X' : '!9' + String("153B) + '!X'
            'lambda':  greekLetter = (capital) ? '!9' + String("114B) + '!X' : '!9' + String("154B) + '!X'
            'mu':      greekLetter = (capital) ? '!9' + String("115B) + '!X' : '!9' + String("155B) + '!X'
            'nu':      greekLetter = (capital) ? '!9' + String("116B) + '!X' : '!9' + String("156B) + '!X'
            'xi':      greekLetter = (capital) ? '!9' + String("130B) + '!X' : '!9' + String("170B) + '!X'
            'omicron': greekLetter = (capital) ? '!9' + String("117B) + '!X' : '!9' + String("157B) + '!X'
            'pi':      greekLetter = (capital) ? '!9' + String("120B) + '!X' : '!9' + String("160B) + '!X'
            'rho':     greekLetter = (capital) ? '!9' + String("122B) + '!X' : '!9' + String("162B) + '!X'
            'sigma':   greekLetter = (capital) ? '!9' + String("123B) + '!X' : '!9' + String("163B) + '!X'
            'tau':     greekLetter = (capital) ? '!9' + String("124B) + '!X' : '!9' + String("164B) + '!X'
            'upsilon': greekLetter = (capital) ? '!9' + String("125B) + '!X' : '!9' + String("165B) + '!X'
            'phi':     greekLetter = (capital) ? '!9' + String("106B) + '!X' : '!9' + String("152B) + '!X'
            'chi':     greekLetter = (capital) ? '!9' + String("103B) + '!X' : '!9' + String("143B) + '!X'
            'psi':     greekLetter = (capital) ? '!9' + String("131B) + '!X' : '!9' + String("171B) + '!X'
            'omega':   greekLetter = (capital) ? '!9' + String("127B) + '!X' : '!9' + String("167B) + '!X'
            ELSE: Message, 'The greek letter ' + letter + ' is unrecognized.'
        ENDCASE
    
    ENDIF ELSE BEGIN
    
        CASE StrLowCase(letter) OF
            'alpha':   greekLetter = (capital) ? '!4' + String("101B) + '!X' : '!4' + String("141B) + '!X'
            'beta':    greekLetter = (capital) ? '!4' + String("102B) + '!X' : '!4' + String("142B) + '!X'
            'gamma':   greekLetter = (capital) ? '!4' + String("103B) + '!X' : '!4' + String("143B) + '!X'
            'delta':   greekLetter = (capital) ? '!4' + String("104B) + '!X' : '!4' + String("144B) + '!X'
            'epsilon': greekLetter = (capital) ? '!4' + String("105B) + '!X' : '!4' + String("145B) + '!X'
            'zeta':    greekLetter = (capital) ? '!4' + String("106B) + '!X' : '!4' + String("146B) + '!X'
            'eta':     greekLetter = (capital) ? '!4' + String("107B) + '!X' : '!4' + String("147B) + '!X'
            'theta':   greekLetter = (capital) ? '!4' + String("110B) + '!X' : '!4' + String("150B) + '!X'
            'iota':    greekLetter = (capital) ? '!4' + String("111B) + '!X' : '!4' + String("151B) + '!X'
            'kappa':   greekLetter = (capital) ? '!4' + String("112B) + '!X' : '!4' + String("152B) + '!X'
            'lambda':  greekLetter = (capital) ? '!4' + String("113B) + '!X' : '!4' + String("153B) + '!X'
            'mu':      greekLetter = (capital) ? '!4' + String("114B) + '!X' : '!4' + String("154B) + '!X'
            'nu':      greekLetter = (capital) ? '!4' + String("115B) + '!X' : '!4' + String("155B) + '!X'
            'xi':      greekLetter = (capital) ? '!4' + String("116B) + '!X' : '!4' + String("156B) + '!X'
            'omicron': greekLetter = (capital) ? '!4' + String("117B) + '!X' : '!4' + String("157B) + '!X'
            'pi':      greekLetter = (capital) ? '!4' + String("120B) + '!X' : '!4' + String("160B) + '!X'
            'rho':     greekLetter = (capital) ? '!4' + String("121B) + '!X' : '!4' + String("161B) + '!X'
            'sigma':   greekLetter = (capital) ? '!4' + String("122B) + '!X' : '!4' + String("162B) + '!X'
            'tau':     greekLetter = (capital) ? '!4' + String("123B) + '!X' : '!4' + String("163B) + '!X'
            'upsilon': greekLetter = (capital) ? '!4' + String("124B) + '!X' : '!4' + String("164B) + '!X'
            'phi':     greekLetter = (capital) ? '!4' + String("125B) + '!X' : '!4' + String("165B) + '!X'
            'chi':     greekLetter = (capital) ? '!4' + String("126B) + '!X' : '!4' + String("166B) + '!X'
            'psi':     greekLetter = (capital) ? '!4' + String("127B) + '!X' : '!4' + String("167B) + '!X'
            'omega':   greekLetter = (capital) ? '!4' + String("130B) + '!X' : '!4' + String("170B) + '!X'
            ELSE: Message, 'The greek letter ' + letter + ' is unrecognized.'      
        ENDCASE
        
      ENDELSE

    RETURN, greekLetter
    
END
;;==============================================================================
; docformat = 'rst'
;
; NAME:
;   cgQuery
;
; PURPOSE:
;   Provides information about any cgWindow applications currently on the display. Returns
;   the window index numbers of any cgWindow applications current on the display.
;
;******************************************************************************************;
;                                                                                          ;
;  Copyright (c) 2011, by Fanning Software Consulting, Inc. All rights reserved.           ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;
;+
;   Provides information about any cgWindow applications currently on the display. Returns
;   the window index numbers of any cgWindow applications current on the display.
;
; :Categories:
;    Graphics
;    
; :Keywords:
;     count: out, optional, type=long
;         The number of cgWindow applications currently on the display.
;     current: in, optional, type=boolean
;         If set, the current cgWindow application information is returned in the result
;         of the function and in the information keywords.
;     dimensions: out, optional, type=integer
;         The dimensions of the ctWindow application, [xdim, ydim, n].
;     objectref: out, optional, type=object
;         A vector of FSC_CMDWINDOW object references for each cgWindow application currently 
;         on the display.
;     title: out, optional, type=string
;         A vector of window titles for each cgWindow application currently on the display.
;     widgetID: out, optional, type=long
;         A vector of widget identifiers of the top-level base widget for each cgWindow
;         application currently on the display.
;          
; :Returns:
;      windowIndexID: out, type=long
;          An array of window index numbers for each cgWindow application currently on the display.
;          
; :Examples:
;    Used as a query routine::
;       IDL> wids = cgQuery(TITLE=titles, COUNT=count)
;       IDL> index = Where(StrUpCase(titles) EQ 'PLOT WINDOW', tcnt)
;       IDL> IF tcnt GT 0 THEN cgSet, wids[index]
;       IDL> cgWindow, 'Oplot', thisData, /AddCmd
;       
; :Author:
;       FANNING SOFTWARE CONSULTING::
;           David W. Fanning 
;           1645 Sheely Drive
;           Fort Collins, CO 80526 USA
;           Phone: 970-221-0438
;           E-mail: david@idlcoyote.com
;           Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Written, 23 January 2011. DWF.
;        Added DIMENSIONS keyword to return current dimensions of cgWindows. 24 Feb 2011. DWF.
;        Made sure this program only returns information on devices that support windows. 20 July 2011. DWF.
;        
; :Copyright:
;     Copyright (c) 2011, Fanning Software Consulting, Inc.
;-
FUNCTION cgQuery, $
    COUNT=count, $
    CURRENT=current, $
    DIMENSIONS=dimensions, $
    OBJECTREF=objectRef, $
    TITLE=title, $
    WIDGETID=widgetID
    
    ; This can only be done in devices that support windows.
    IF ~((!D.Flags AND 256) NE 0) THEN BEGIN
        count = 0
        RETURN, -1
    ENDIF

    ; Are there cgWindow applications around?
    DefSysV, '!FSC_WINDOW_LIST', EXISTS=exists
    
    ; If it doesn't exist, or it is invalid, leave.
    IF ~exists THEN BEGIN
        count = 0
        RETURN, -1
    ENDIF ELSE BEGIN
        IF ~Obj_Valid(!FSC_WINDOW_LIST) THEN BEGIN
            count = 0
            RETURN, -1
        ENDIF   
    ENDELSE
    
    ; Get the window list and find out how many windows there are.
    list = !FSC_WINDOW_LIST
    count = list -> Get_Count()
    
    ; Make arrays to hold the values.
    widgetID = LonArr(count)
    objectRef = ObjArr(count)
    title = StrArr(count)
    windowIndex = IntArr(count)
    dimensions = IntArr(2, count)
    
    ; Fill them up.
    thisWindow = !D.Window
    FOR j=0,count-1 DO BEGIN
        thisItem = list -> Get_Item(j, /DEREFERENCE)
        
        ; Make sure this is a valid widget.
        IF Widget_Info(thisItem.tlb, /VALID_ID) THEN BEGIN
           widgetID[j] = thisItem.tlb
        ENDIF ELSE BEGIN
           widgetID[j] = -1
           CONTINUE
        ENDELSE
        objectRef[j] = thisItem.windowobj
        title[j] = thisItem.title
        windowIndex[j] = thisItem.wid
        dimensions[*,j] = [!D.X_Size, !D.Y_Size]
    ENDFOR
    
    ; IF you have bad indices, clean things up.
    badIndices = Where(widgetID EQ -1, badCount, COMPLEMENT=goodIndices, NCOMPLEMENT=count)
    IF badCount GT 0 THEN BEGIN
       widgetID = widgetID[goodIndices]
       objectRef[j] = objectRef[goodIndices]
       title[j] =  title[goodIndices]
       windowIndex[j] =  windowIndex[goodIndices]
       dimensions[*,j] =  dimensions[*,goodIndices]
       FOR k=0,badCount-1 DO list -> Destroy, badIndices[k], /Destroy
    ENDIF
    IF (thisWindow GE 0) && WindowAvailable(thisWindow) THEN WSet, thisWindow ELSE WSet, -1
    
    ; Return just the current values if the CURRENT keyword is set.
    IF Keyword_Set(current) THEN BEGIN
        IF (count-1) GE 0 THEN BEGIN
            windowIndex = windowIndex[count-1]
            widgetID = widgetID[count-1]
            objectRef = objectRef[count-1]
            title = title[count-1] 
            dimensions = dimensions[*,count-1]
        ENDIF ELSE BEGIN
            windowIndex = -1
            widgetID = -1
            objectRef = Obj_New()
            title = ""
            dimensions = [0,0]
        ENDELSE
    ENDIF

    ; Make sure scalar is returned if just one element.
    IF count EQ 1 THEN BEGIN
        widgetID = widgetID[0]
        objectRef = objectRef[0]
        title = title[0]
        windowIndex = windowIndex[0]
        dimensions = dimensions[*,0]
    ENDIF
    
    
    RETURN, windowIndex
    
END
;;==============================================================================
; docformat = 'rst'
;
; NAME:
;   cgSnapshot
;
; PURPOSE:
;   To get accurate screen dumps with the IDL command TVRD on 24-bit
;   PC and Macintosh computers, you have to be sure to set color
;   decomposition on. This program adds that capability automatically.
;   In addition, the program will optionally write BMP, GIF, JPEG,
;   PICT, PNG, and TIFF color image files of the screen dump.
;
;******************************************************************************************;
;                                                                                          ;
;  Copyright (c) 2011, by Fanning Software Consulting, Inc. All rights reserved.           ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;
;+
; To get accurate screen dumps with the IDL command TVRD on 24-bit
; PC and Macintosh computers, you have to be sure to set color
; decomposition on. This program adds that capability automatically.
; In addition, the program will optionally write BMP, GIF, JPEG,
; PICT, PNG, and TIFF color image files of the screen dump.
;
; :Categories:
;    Graphics
;    
; :Returns:
;    The returned image will be a 2D image on 8-bit systems and a 24-bit pixel 
;    interleaved true-color image on 24-bit systems. A -1 will be returned if a 
;    file output keyword is used (e.g., JPEG, TIFF, etc.).
;    
; :Params:
;    xstart: in, optional, type=integer, default=0
;       The starting column index of the rectantular area that is to be copied.
;    ystart: in, optional, type=integer, default=0
;       The starting row index of the rectantular area that is to be copied.
;    ncols: in, optional, type=integer
;       The number of columns to read in the rectantular area that is to be 
;       copied. By default, !D.X_Size - xstart.
;    nrows: in, optional, type=integer
;       The number of rows to read in the rectantular area that is to be 
;       copied. By default, !D.Y_Size - ystart.
;
; :Keywords:
;    bmp: in, optional, type=boolean, default=0
;        Set this keyword to write the screen dump as a color BMP file.
;    cancel: out, optional, type=boolean, default=0
;        An output keyword set to 1 if the user cancels out of a filename dialog. 
;        Set to 0 otherwise.
;    colors: in, optional, type=integer, default=256
;        If a 24-bit image has to be quantized, this will set the number of colors in 
;        the output image. Applies to BMP, GIF, PICT, and PNG formats written from 
;        24-bit displays.(See the COLOR_QUAN documentation for details.)
;    cube: in, optional, type=integer
;        If this keyword is set to a value between 2 and 6 the color quantization will 
;        use a cubic method of quantization. Applies to BMP, GIF, PICT, and PNG formats 
;        written from 24-bit displays.(See the COLOR_QUAN documentation for details.)
;    dither: in, optional, type=boolean, default=0 
;        If this keyword is set the quantized image will be dithered. Applies to BMP, 
;        GIF, PICT, and PNG formats written from 24-bit displays.(See the COLOR_QUAN 
;        documentation for details.)
;    filename: in, optional, type=string
;        The name of the output file. If you specify a name with a file extension of the
;        type of file you want to create (e.g., *.jpg, *.png, etc), then you do not have
;        to use the file type keywords (e.g., JPEG, PNG, etc.). Otherwise, you can specify
;        the name of the the file without an extension, use the file keywords, and a file
;        extension will be added to the filename automatically, depending upon the type of
;        output file selected.
;    gif: in, optional, type=boolean, default=0
;        Set this keyword to write the screen dump as a color GIF file.
;    jpeg: in, optional, type=boolean, default=0
;        Set this keyword to write the screen dump as a color JPEG file.
;    nodialog: in, optional, type=boolean, default=0        
;        Set this keyword if you wish to avoid the DIALOG_PICKFILE dialog that asks you 
;        to name the output file. This keyword should be set, for example, if you are 
;        processing screens in batch mode.
;    order: in, optional, type=boolean, default=0
;        Set this keyword to determine the image order for reading the display. Corresponds to 
;        !Order and set to such as the default.
;    overwrite_prompt: in, optional, type=boolean, default=0       
;        Set this keyword if you would like to get a prompt if you are overwriting a file. 
;        This applies only to operations involving DIALOG_PICKFILE.
;    pict: in, optional, type=boolean, default=0
;        Set this keyword to write the screen dump as a color PICT file.
;    png: in, optional, type=boolean, default=0
;        Set this keyword to write the screen dump as a color PNG file.
;    position: in, optional, type=float
;        An alternative way of setting the `xstart`, `ystart`, `ncols` and `nrows` parameters
;        by specifying a four-element normalized array, [x0,y0,x1,y1].
;    tiff: in, optional, type=boolean, default=0
;        Set this keyword to write the screen dump as a color TIFF file.
;    true: in, optional, type=integer, default=1
;        Set this keyword to the type of interleaving you want. 1 = Pixel interleaved, 
;        2 = row interleaved, 3 = band interleaved.
;    type: in, optional, type=string
;        Set this keyword to the type of file to write. Use this instead of
;        setting BMP, GIF, JPEG, PICT, PNG, or TIFF keywords: TYPE='JPEG'. The
;        primary purpose of this is to make widget event handlers easier to write.
;    quality: in, optional, type=integer, default=75
;        This keyword sets the amount of compression for JPEG images. It should be set to a 
;        value between 0 and 100. (See the WRITE_JPEG documentation for details.)
;    wid: in, optional, type=integer
;        The index number of the window to read from. The current graphics window
;        (!D.Window) is selected by default. An error is issued if no windows are
;         currently open on a device that supports windows.
;    _ref_extra: in, optional
;        Any keywords that are appropriate for the WRITE_*** routines are also accepted via 
;        keyword inheritance.
;          
; :Examples:
;    To obtain an image of the current graphics window::
;    
;       IDL> image = cgSnapshot()
;       
;    To create a PNG file, named "test.png", of the current graphics window::
;    
;       IDL> void = cgSnapshot(FILENAME='test.png')
;       
;    To obtain the lower quadrant of a 512-by-512 graphics window as a
;    band interleaved image::
;    
;       IDL> image = cgSnapshot(0, 0, 256, 256, TRUE=3)
;       
; :Author:
;       FANNING SOFTWARE CONSULTING::
;           David W. Fanning 
;           1645 Sheely Drive
;           Fort Collins, CO 80526 USA
;           Phone: 970-221-0438
;           E-mail: david@idlcoyote.com
;           Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Renamed TVRead to cgSnapshot and retired TVRead. 20 February 2011. DWF.
;        Added the ability to get the file type from the file name extension. 26 Dec 2011. DWF.
;        Added a POSITION keyword to select a position inside the window for capture. 20 October 2012. DWF.
;        Fixed a problem with not setting back to incoming decomposed state on an error. 20 Nov 2012. DWF.
; :Copyright:
;     Copyright (c) 2011, Fanning Software Consulting, Inc.
;-
FUNCTION cgSnapshot, xstart, ystart, ncols, nrows, $
   BMP=bmp, $
   Cancel=cancel, $
   Colors=colors, $
   Cube=cube, $
   Dither=dither, $
   Filename=filename, $
   GIF=gif, $
   JPEG=jpeg, $
   NoDialog=nodialog, $
   Order=order, $
   Overwrite_Prompt=overwrite_prompt, $
   PICT=pict, $
   PNG=png, $
   POSITION=position, $
   TIFF=tiff, $
   True=true, $
   Type=type, $
   Quality=quality, $
   WID=wid, $
   _Ref_Extra=extra
   

   ; Error handling.
   Catch, theError
   IF theError NE 0 THEN BEGIN
       Catch, /Cancel
       ok = Error_Message()
       IF N_Elements(thisWindow) EQ 0 THEN RETURN, -1
       IF thisWindow GE 0 THEN WSet, thisWindow
       
       ; Need to set color decomposition back?
       IF (N_Elements(theDecomposedState) NE 0) && (theDepth GT 0) THEN BEGIN
           Device, Decomposed=theDecomposedState
       ENDIF
       RETURN, -1
    ENDIF
    
    cancel = 0
    
    ; Check for availability of GIF files.
    thisVersion = Float(!Version.Release)
    IF (thisVersion LT 5.3) OR (thisVersion GE 6.1) THEN haveGif = 1 ELSE haveGIF = 0
    
    ; Go to correct window.
    IF N_Elements(wid) EQ 0 THEN wid =!D.Window
    thisWindow = !D.Window
    IF (!D.Flags AND 256) NE 0 THEN WSet, wid
    
    ; Did the user specify a normalized position in the window?
    IF N_Elements(position) NE 0 THEN BEGIN
       xstart = position[0] * !D.X_VSize
       ystart = position[1] * !D.Y_VSize
       ncols = (position[2]*!D.X_VSize) - xstart
       nrows = (position[3]*!D.Y_VSize) - ystart
    ENDIF
    
    ; Check keywords and parameters. Define values if necessary.
    IF N_Elements(xstart) EQ 0 THEN xstart = 0
    IF N_Elements(ystart) EQ 0 THEN ystart = 0
    IF N_Elements(ncols) EQ 0 THEN ncols = !D.X_VSize - xstart
    IF N_Elements(nrows) EQ 0 THEN nrows = !D.Y_VSize - ystart
    IF N_Elements(order) EQ 0 THEN order = !Order
    IF N_Elements(true) EQ 0 THEN true = 1
    dialog = 1 - Keyword_Set(nodialog)
    
    ; Is the FILENAME keyword being used? If so, get the type of the
    ; file from the filename extension.
    IF N_Elements(filename) NE 0 THEN BEGIN
       root_name = cgRootName(filename, DIRECTORY=theDir, EXTENSION=ext)
       IF ext NE "" THEN BEGIN
           type = StrUpCase(ext)
           typeFromExtension = 1
       ENDIF ELSE typeFromExtension = 0
    ENDIF ELSE typeFromExtension = 0
    
    ; Do you want to write an image file instead of capturing an image?
    IF N_Elements(type) NE 0 THEN BEGIN
       CASE StrUpCase(type) OF
          'BMP': bmp = 1
          'GIF': gif = 1
          'JPEG': jpeg = 1
          'JPG': jpeg = 1
          'PICT': pict = 1
          'PNG': png = 1
          'TIFF': tiff = 1
          'TIF': tif = 1
          ELSE: Message, 'Cannot write a file of type: ' + StrUpCase(type) + '.'
       ENDCASE
    ENDIF
    writeImage = 0
    fileType = ""
    extention = ""
    IF Keyword_Set(bmp)THEN BEGIN
       writeImage = 1
       fileType = 'BMP'
       extension = 'bmp'
    ENDIF
    IF Keyword_Set(gif) THEN BEGIN
       IF havegif THEN BEGIN
          writeImage = 1
          fileType = 'GIF'
          extension = 'gif'
        ENDIF ELSE BEGIN
           ok = Dialog_Message('GIF files not supported in this IDL version. Replacing with JPEG.')
           writeImage = 1
          fileType = 'JPEG'
          extension = 'jpg'
       ENDELSE
    ENDIF
    IF Keyword_Set(jpeg) THEN BEGIN
       writeImage = 1
       fileType = 'JPEG'
       extension = 'jpg'
    ENDIF
    IF Keyword_Set(PICT) THEN BEGIN
       writeImage = 1
       fileType = 'PICT'
       extension = 'pict'
    ENDIF
    IF Keyword_Set(png) THEN BEGIN
       writeImage = 1
       fileType = 'PNG'
       extension = 'png'
    ENDIF
    IF Keyword_Set(tiff) THEN BEGIN
       writeImage = 1
       fileType = 'TIFF'
       extension = 'tif'
    ENDIF
    
    IF N_Elements(colors) EQ 0 THEN colors = 256
    IF N_Elements(quality) EQ 0 THEN quality = 75
    dither = Keyword_Set(dither)
    
    ; On 24-bit displays, make sure color decomposition is ON.
    IF (!D.Flags AND 256) NE 0 THEN BEGIN
       Device, Get_Decomposed=theDecomposedState, Get_Visual_Depth=theDepth
       IF theDepth GT 8 THEN BEGIN
          Device, Decomposed=1
          IF theDepth EQ 24 THEN truecolor = true ELSE truecolor = 0
       ENDIF ELSE truecolor = 0
       IF wid LT 0 THEN $
          Message, 'No currently open windows. Returning.', /NoName
    ENDIF ELSE BEGIN
       truecolor = 0
       theDepth = 8
    ENDELSE
    
    ; Fix for 24-bit Z-buffer.
    IF (Float(!Version.Release) GE 6.4) AND (!D.NAME EQ 'Z') THEN BEGIN
       Device, Get_Decomposed=theDecomposedState, Get_Pixel_Depth=theDepth
       IF theDepth EQ 24 THEN truecolor = true ELSE truecolor = 0
    ENDIF
    
   ; Get the screen dump. 2D image on 8-bit displays. 3D image on 24-bit displays.
    image = TVRD(xstart, ystart, ncols, nrows, True=truecolor, Order=order)
    
    ; Need to set color decomposition back?
    IF theDepth GT 8 THEN Device, Decomposed=theDecomposedState
    
    ; If we need to write an image, do it here.
    IF writeImage THEN BEGIN
    
       ; Get the name of the output file.
       IF N_Elements(filename) EQ 0 THEN BEGIN
          filename = 'idl.' + StrLowCase(extension)
       ENDIF ELSE BEGIN
          IF typeFromExtension EQ 0 THEN filename = filename + "." + StrLowCase(extension)
       ENDELSE
       IF dialog THEN filename = Dialog_Pickfile(/Write, File=filename, OVERWRITE_PROMPT=Keyword_Set(overwrite_prompt))
    
       IF filename EQ "" THEN BEGIN
          cancel = 1
          RETURN, image
       ENDIF
    
       ; Write the file.
       CASE fileType OF
    
          'BMP': BEGIN
             IF truecolor THEN BEGIN
                ; BMP files assume blue, green, red planes.
                temp = image[0,*,*]
                image[0,*,*] = image[2,*,*]
                image[2,*,*] = temp
                Write_BMP, filename, image, _Extra=extra
             ENDIF ELSE BEGIN
                TVLCT, r, g, b, /Get
                Write_BMP, filename, image, r, g, b, _Extra=extra
             ENDELSE
             END
    
          'GIF': BEGIN
             IF truecolor THEN BEGIN
                CASE Keyword_Set(cube) OF
                   0: image2D = Color_Quan(image, 1, r, g, b, Colors=colors, Dither=dither)
                   1: image2D = Color_Quan(image, 1, r, g, b, Cube=2 > cube < 6)
                ENDCASE
             ENDIF ELSE BEGIN
                TVLCT, r, g, b, /Get
                image2D = image
             ENDELSE
             Write_GIF, filename, image2D, r, g, b, _Extra=extra
             END
    
          'JPEG': BEGIN
             IF truecolor THEN BEGIN
                image3D = image
             ENDIF ELSE BEGIN
                s = Size(image, /Dimensions)
                image3D = BytArr(3, s[0], s[1])
                TVLCT, r, g, b, /Get
                image3D[0,*,*] = r[image]
                image3D[1,*,*] = g[image]
                image3D[2,*,*] = b[image]
             ENDELSE
             Write_JPEG, filename, image3D, True=1, Quality=quality, _Extra=extra
             END
    
          'PICT': BEGIN
             IF truecolor THEN BEGIN
                CASE Keyword_Set(cube) OF
                   0: image2D = Color_Quan(image, 1, r, g, b, Colors=colors, Dither=dither)
                   1: image2D = Color_Quan(image, 1, r, g, b, Cube=2 > cube < 6)
                ENDCASE
             ENDIF ELSE BEGIN
                TVLCT, r, g, b, /Get
                image2D = image
             ENDELSE
             Write_PICT, filename, image2D, r, g, b
             END
    
          'PNG': BEGIN
             IF truecolor THEN BEGIN
                Write_PNG, filename, image, _Extra=extra
             ENDIF ELSE BEGIN
                TVLCT, r, g, b, /Get
                image2D = image
                Write_PNG, filename, image2D, r, g, b, _Extra=extra
             ENDELSE
             END
    
          'TIFF': BEGIN
             IF truecolor THEN BEGIN
                image3D = Reverse(image,3)
             ENDIF ELSE BEGIN
                s = Size(image, /Dimensions)
                image3D = BytArr(3, s[0], s[1])
                TVLCT, r, g, b, /Get
                image3D[0,*,*] = r[image]
                image3D[1,*,*] = g[image]
                image3D[2,*,*] = b[image]
                image3D = Reverse(Temporary(image3D), 3)
             ENDELSE
             Write_TIFF, filename, image3D, 1, _Extra=extra
             END
       ENDCASE
       RETURN, -1
    ENDIF
    
    ; Return the screen dump image.
    RETURN, image
    
END ;-------------------------------------------------------------------------------
;;==============================================================================
; docformat = 'rst'
;
; NAME:
;   cgSymbol
;
; PURPOSE:
; This function provides a device-independent way to ask for commonly-used
; symbols (e.g., less than or equal to, Angstrom, degree symbol, etc.),
; including all 24 Greek symbols in upper and lower case, to display with text.
; 
;******************************************************************************************;
;                                                                                          ;
;  Copyright (c) 2011, by Fanning Software Consulting, Inc. All rights reserved.           ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;
;+
; This function provides a device-independent way to ask for commonly-used
; symbols (e.g., less than or equal to, Angstrom, degree symbol, subscripts, superscripts, etc.),
; including all 24 Greek symbols in upper and lower case, to display with text.
;
; For this program to work correctly on your graphics display, you should be using
; Hershey fonts (!P.Font=-1). It will work correctly in PostScript with either 
; hardware fonts (!P.Font=0) or True-Type fonts (!P.Font=1).
; 
; `Greek Symbols <http://www.idlcoyote.com/ps_tips/greeksym.html>` are created by
; calling the Coyote Library routine `cgGreek <http://www.idlcoyote.com/programs/cggreek.pro>' 
; from this program.
; 
; .. image:: cgsymbol.png 
; 
; Normally, rather than calling cgSymbol, the symbols are embedded in Coyote Graphics
; text that are used for axis annotation and so forth. See `Embedding Symbols in Coyote
; Graphics Output <http://www.idlcoyote.com/cg_tips/embedsymbols.php>`. Embedded subscripts
; and superscripts are implemented like this::
;     $\upTEXT$ : Raise the text in TEXT to superscript level and reduce character size by 62%.
;     $\expTEXT$ : Raise the text in TEXT to superscript level and reduce character size by 44%.
;     $\downTEXT$ : Lower the text in TEXT to subscript level and reduce character size by 62%.
;     $\subTEXT$ : Lower the text in TEXT to subscript level and reduce character size by 44%.
;
; :Categories:
;    Graphics
;    
; :Returns:
;    A string variable that represents the requested symbol and can be used
;    in a textual context.
;    
; :Examples:
;     To create a lowercase Greek psi symbol::
;  
;        IDL> cgPlot, findgen(11), XTitle='This title contains ' + $
;                 cgSymbol('psi') + ' as a Greek letter' 
;        IDL> cgPlot, findgen(11), XTitle='This title contains $\psi$ as a Greek letter'
;
;     To create an Angstrom sign::
;  
;        IDL> cgPlot, findgen(11), XTitle='This title contains (' + $
;                cgSymbol('Angstrom') +  ') an Angstrom sign.' 
;        IDL> cgPlot, findgen(11), XTitle='This title contains ($\Angstrom$) an Angstrom sign.'
;       
;     To create subscripts or superscripts::
;  
;        IDL> cgPlot, findgen(11), XTitle='E=mc$\up2$' 
;        IDL> cgPlot, findgen(11), XTitle='H$\sub2$O'
;        IDL> cgPlot, findgen(11), XTitle='H$\upSuper$MT $\Omega$$\subSubscript$', Charsize=2.0
;        
;     It is possible to use Greek characters as superscripts and subscripts. Do so by
;     prepending the Greek character with "\\" inside the normal superscript or subscript
;     notation. For example, to use lambda as an exponent to the Greek character Omega, you
;     can type this::
;
;        IDL> cgPlot, cgDemoData(1), XTitle='$\Omega$$\exp\\lambda$', Charsize=2.0
;
;     To use lambda as a subscript, type this:
;
;         IDL> cgPlot, cgDemoData(1), XTitle='$\Omega$$\sub\\lambda$', Charsize=2.0
;
; :Author:
;       FANNING SOFTWARE CONSULTING::
;           David W. Fanning 
;           1645 Sheely Drive
;           Fort Collins, CO 80526 USA
;           Phone: 970-221-0438
;           E-mail: david@idlcoyote.com
;           Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Written by: David W. Fanning, 2 September 2011. 
;        Added plus-minus symbol. 2 Nov 2011. DWF.
;        Added "up", "down", "exp" "sub" and "n" symbols for subscripting and superscripting. 9 Nov 2012. DWF.
;        Added "division" and "times" signs. 24 Nov 2012. DWF.
;        Updated UNICODE values to display capital letters correctly. 23 Dec 2012. DWF.
;
; :Copyright:
;     Copyright (c) 2011, Fanning Software Consulting, Inc.
;-
;+
; Displays the symbols and their names in a graphics window.
; 
; :Keywords:
;    ps: in, optional, type=boolean, default=0                                
;        Normally, the PostScript version of the symbol is returned automatically if
;        the current device is PostScript and !P.Font is 0 or 1. But, the PostScript version 
;        of the symbol can be obtained at any time and in any device, by setting this keyword.
;    unicode: in, optional, type=boolean, default=0                                
;        If this keyword is set, the function returns the Unicode value for the symbol.
; 
;-
PRO cgSymbol_Example, PS=ps, UNICODE=unicode

    Forward_Function cgSymbol

    Compile_Opt hidden
    
    ; The symbols..
    symbol = [  'alpha', 'beta', 'gamma', 'delta', 'epsilon',  'zeta', $
                'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu', $
                'nu', 'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau', $
                'upsilon', 'phi', 'chi', 'psi', 'omega', 'leq', 'geq', $
                'neq', 'deg', '+-', 'equiv', 'prime', 'angstrom', 'sun', $
                'varphi', 'infinity', 'copyright' ]
    
    ; Output positions.
    x = [0.15, 0.4, 0.65]
    y = Reverse((Indgen(13) + 1) * (1.0 / 14))
    
    ; Need PostScript output?
    IF Keyword_Set(ps) || (!D.Name EQ 'PS') THEN PS_Start
    
    ; Create a window.
    cgDisplay, 800, 500
    
    ; Output the symbols.
    FOR j=0,11 DO BEGIN
        cgText, x[0], y[j], symbol[j] + ': ' + $
            cgSymbol(symbol[j], UNICODE=unicode) + cgSymbol(symbol[j], /CAPITAL, UNICODE=unicode), $
            /NORMAL, CHARSIZE=1.5
        cgText, x[1], y[j], symbol[j+12] + ': ' + $
            cgSymbol(symbol[j+12], UNICODE=unicode) + cgSymbol(symbol[j+12], /CAPITAL, UNICODE=unicode), $
            /NORMAL, CHARSIZE=1.5
        IF j NE 12 THEN cgText, x[2], y[j], symbol[j+24] + ': ' + $
            cgSymbol(symbol[j+24], UNICODE=unicode), $
            /NORMAL, CHARSIZE=1.5
    ENDFOR
    
    ; Clean up PostScript, if needed.
    IF Keyword_Set(ps) || (!D.Name EQ 'PS') THEN PS_End
    
    ; Restore the users window.
    IF N_Elements(thisWindow) NE 0 THEN BEGIN
       IF thisWindow GE 0 THEN WSet, thisWindow
    ENDIF
    
END ; --------------------------------------------------------------------------------------


;+
; This function provides a device-independent way to ask for commonly-used
; symbols (e.g., less than or equal to, Angstrom, degree symbol, etc.),
; including all 24 Greek symbols in upper and lower case, to display with text.
;
; :Params:
;    symbol: in, required, type=string, default="alpha"
;       The name of the symbol desired as a string. Valid string names are the 24 
;       characters of the Greek alphabet, plus commonly used graphical symbols::
;       
;          alpha        nu        leq
;          beta         xi        geq
;          gamma        omicron   neg
;          delta        pi        deg
;          epsilon      rho       equiv
;          zeta         sigma     prime
;          eta          tau       angstrom
;          theta        upsilon   sun
;          iota         phi       varphi
;          kappa        chi       infinity
;          lambda       psi       copyright
;          mu           omega
;                    
;       Note that if the first letter of the name is capitalized, this is
;       the equivalent of setting the `Capital` keyword. Only Greek characters
;       use this method of selecting symbols.
;       
; :Keywords:
;    capital: in, optional, type=boolean, default=0
;        If this keyword is set, the captial Greek letter is returned rather than the lowercase 
;        Greek letter. An alternative way of capitalizing the letter is to make the first letter 
;        of the name an uppercase letter.
;    example: in, optional, type=boolean, default=0             
;        If this keyword is set, the names of symbols and the symbols themselves are written 
;        out in the current graphics window.
;    ps: in, optional, type=boolean, default=0                                
;        Normally, the PostScript version of the symbol is returned automatically if
;        the current device is PostScript and !P.Font is 0 or 1. But, the PostScript version 
;        of the symbol can be obtained at any time and in any device, by setting this keyword.
;    unicode: in, optional, type=boolean, default=0                                
;        If this keyword is set, the function returns the Unicode value for the symbol.
;          
;-
FUNCTION cgSymbol, symbol, CAPITAL=capital, EXAMPLE=example, PS=PS, UNICODE=unicode

    Compile_Opt idl2
    
    ; Return to caller on error.
    ON_Error, 2
    
    ; A common block to communicate with PS_Start.
    COMMON _$FSC_PS_START_, ps_struct
    
    ; Do you wish to see an example?
    IF Keyword_Set(example) THEN BEGIN
        cgSymbol_Example, UNICODE=unicode, PS=ps
        RETURN, ""
    ENDIF

    ; A symbol name is required.
    IF N_Elements(symbol) EQ 0 THEN BEGIN
       Print, 'Syntax: cgSymbol("theSymbol")
       RETURN, ""
    ENDIF
    
    ; Return quietly if the symbol is a null string.
    IF symbol EQ "" THEN RETURN, ""
    
    ; Set up PostScript device for working with Greek letters.
    IF !D.Name EQ 'PS' THEN Device, ISOLATIN1=1
    
    ; Check keywords.
    capital = Keyword_Set(capital)
    
    ; If the first letter of the "symbol" variable is uppercase
    ; the user wants an uppercase symbol, if available.
    firstLetter = StrMid(symbol, 0, 1)
    IF firstLetter EQ StrUpCase(firstLetter) THEN capital = 1
    
    IF Keyword_Set(unicode) THEN BEGIN

        CASE StrLowCase(symbol) OF
            'alpha':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'beta':    retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'gamma':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'delta':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'epsilon': retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'zeta':    retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'eta':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'theta':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'iota':    retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'kappa':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'lambda':  retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'mu':      retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'nu':      retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'xi':      retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'omicron': retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'pi':      retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'rho':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'sigma':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'tau':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'upsilon': retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'phi':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'chi':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'psi':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'omega':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'leq':     retSymbol = '!Z(2266)'
            'geq':     retSymbol = '!Z(2267)'
            'neq':     retSymbol = '!Z(2260)'
            'deg':     retSymbol = '!Z(00B0)'
            '+-':      retSymbol = '!Z(00B1)'
            'equiv':   retSymbol = '!Z(2261)'
            'prime':   retSymbol = '!Z(2232)'
            'angstrom':retSymbol = '!Z(00C5)'
            'sun':     retSymbol = '!Z(2609)'
            'varphi':  retSymbol = '!Z(03D5)'
            'infinity':retSymbol = '!Z(221E)'
            'copyright': retSymbol = '!Z(00A9)'
            'up':      retSymbol = '!U'
            'down':    retSymbol = '!D'
            'exp':     retSymbol = '!E'
            'sub':     retSymbol = '!I'
            'n':       retSymbol = '!N'
            'div':     retSymbol = '!Z(00F7)'
            'times':   retSymbol = '!Z(00D7)'
            ELSE: Message, 'The symbol ' + symbol + ' is unrecognized.'
       ENDCASE
           
       ; Return the UNICODE letter.
       RETURN, retSymbol
       
    ENDIF
    
    IF ((!D.Name EQ 'PS') && (!P.FONT NE -1)) || Keyword_Set(PS)  THEN BEGIN
    
        ; Make sure ISOLATIN1 encoding is turned on.
        IF (!D.Name EQ 'PS') THEN DEVICE, /ISOLATIN1
        
        CASE StrLowCase(symbol) OF
            'alpha':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'beta':    retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'gamma':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'delta':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'epsilon': retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'zeta':    retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'eta':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'theta':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'iota':    retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'kappa':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'lambda':  retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'mu':      retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'nu':      retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'xi':      retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'omicron': retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'pi':      retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'rho':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'sigma':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'tau':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'upsilon': retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'phi':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'chi':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'psi':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'omega':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'leq':     retSymbol = '!9' + String("243B) + '!X'
            'geq':     retSymbol = '!9' + String("263B) + '!X'
            'neq':     retSymbol = '!9' + String("271B) + '!X' 
            'deg':     retSymbol = '!9' + String("260B) + '!X'
            '+-':      retSymbol = '!9' + String("261B) + '!X'
            'equiv':   retSymbol = '!9' + String("272B) + '!X'
            'prime':   retSymbol = '!9' + String("242B) + '!X'
            'angstrom':retSymbol = '!3' + String("305B) + '!X'
            'sun':     BEGIN
                       IF (ps_struct.font EQ 1) && (StrUpCase(ps_struct.tt_font) EQ 'DEJAVUSANS') THEN BEGIN
                          retSymbol = '!Z(2609)'
                       ENDIF ELSE BEGIN
                         thisDevice = !D.Name
                         Set_Plot, 'PS'
                         Device, /AVANTGARDE, ISOLATIN1=0, /BOOK, FONT_INDEX = 20
                         retSymbol = '!20!S!DO!R!I ' + string(183b) + '!X!N'
                         Set_Plot, thisDevice
                       ENDELSE
                       END
            'varphi':  retSymbol = '!9' + String("152B) + '!X'
            'infinity':retSymbol = '!9' + String("245B) + '!X'
            'copyright':retSymbol = '!9' + String("323B) + '!X'
            'up':      retSymbol = '!U'
            'down':    retSymbol = '!D'
            'exp':     retSymbol = '!E'
            'sub':     retSymbol = '!I'
            'n':       retSymbol = '!N'
            'div':     retSymbol = '!9' + String("270B) + '!X'
            'times':   retSymbol = '!9' + String("264B) + '!X'
           ELSE: Message, 'The symbol ' + symbol + ' is unrecognized.'
        ENDCASE
    
    ENDIF ELSE BEGIN
    
        CASE StrLowCase(symbol) OF
            'alpha':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'beta':    retSymbol = cgGreek(symbol, CAPITAL=capital)
            'gamma':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'delta':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'epsilon': retSymbol = cgGreek(symbol, CAPITAL=capital)
            'zeta':    retSymbol = cgGreek(symbol, CAPITAL=capital)
            'eta':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'theta':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'iota':    retSymbol = cgGreek(symbol, CAPITAL=capital)
            'kappa':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'lambda':  retSymbol = cgGreek(symbol, CAPITAL=capital)
            'mu':      retSymbol = cgGreek(symbol, CAPITAL=capital)
            'nu':      retSymbol = cgGreek(symbol, CAPITAL=capital)
            'xi':      retSymbol = cgGreek(symbol, CAPITAL=capital)
            'omicron': retSymbol = cgGreek(symbol, CAPITAL=capital)
            'pi':      retSymbol = cgGreek(symbol, CAPITAL=capital)
            'rho':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'sigma':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'tau':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'upsilon': retSymbol = cgGreek(symbol, CAPITAL=capital)
            'phi':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'chi':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'psi':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'omega':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'leq':     retSymbol = '!9' + String("154B) + '!X'
            'geq':     retSymbol = '!9' + String("142B) + '!X'
            'neq':     retSymbol = '!9' + String("75B) + '!X' 
            'deg':     retSymbol = '!9' + String("45B) + '!X'
            '+-':      retSymbol = '!9' + String("53B) + '!X'
            'equiv':   retSymbol = '!9' + String("72B) + '!X'
            'prime':   retSymbol = '!9' + String("140B) + '!X'
            'angstrom':retSymbol = '!3' + String("305B) + '!X'
            'sun':     retSymbol = '!D!9n!N!X'
            'varphi':  retSymbol = '!9' + String("120B) + '!X'
            'infinity':retSymbol = '!9' + String("44B) + '!X'
            'copyright':retSymbol = '!3' + String("251B) + '!X'
            'up':      retSymbol = '!U'
            'down':    retSymbol = '!D'
            'exp':     retSymbol = '!E'
            'sub':     retSymbol = '!I'
            'n':       retSymbol = '!N'
            'div':     retSymbol = '!9' + String("57B) + '!X'
            'times':   retSymbol = '!9' + String("130B) + '!X'
            ELSE: Message, 'The symbol ' + symbol + ' is unrecognized.'      
        ENDCASE
        
      ENDELSE

    RETURN, retSymbol
    
END
;;==============================================================================
; docformat = 'rst'
;
; NAME:
;   cgSymbol
;
; PURPOSE:
; This function provides a device-independent way to ask for commonly-used
; symbols (e.g., less than or equal to, Angstrom, degree symbol, etc.),
; including all 24 Greek symbols in upper and lower case, to display with text.
; 
;******************************************************************************************;
;                                                                                          ;
;  Copyright (c) 2011, by Fanning Software Consulting, Inc. All rights reserved.           ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;
;+
; This function provides a device-independent way to ask for commonly-used
; symbols (e.g., less than or equal to, Angstrom, degree symbol, subscripts, superscripts, etc.),
; including all 24 Greek symbols in upper and lower case, to display with text.
;
; For this program to work correctly on your graphics display, you should be using
; Hershey fonts (!P.Font=-1). It will work correctly in PostScript with either 
; hardware fonts (!P.Font=0) or True-Type fonts (!P.Font=1).
; 
; `Greek Symbols <http://www.idlcoyote.com/ps_tips/greeksym.html>` are created by
; calling the Coyote Library routine `cgGreek <http://www.idlcoyote.com/programs/cggreek.pro>' 
; from this program.
; 
; .. image:: cgsymbol.png 
; 
; Normally, rather than calling cgSymbol, the symbols are embedded in Coyote Graphics
; text that are used for axis annotation and so forth. See `Embedding Symbols in Coyote
; Graphics Output <http://www.idlcoyote.com/cg_tips/embedsymbols.php>`. Embedded subscripts
; and superscripts are implemented like this::
;     $\upTEXT$ : Raise the text in TEXT to superscript level and reduce character size by 62%.
;     $\expTEXT$ : Raise the text in TEXT to superscript level and reduce character size by 44%.
;     $\downTEXT$ : Lower the text in TEXT to subscript level and reduce character size by 62%.
;     $\subTEXT$ : Lower the text in TEXT to subscript level and reduce character size by 44%.
;
; :Categories:
;    Graphics
;    
; :Returns:
;    A string variable that represents the requested symbol and can be used
;    in a textual context.
;    
; :Examples:
;     To create a lowercase Greek psi symbol::
;  
;        IDL> cgPlot, findgen(11), XTitle='This title contains ' + $
;                 cgSymbol('psi') + ' as a Greek letter' 
;        IDL> cgPlot, findgen(11), XTitle='This title contains $\psi$ as a Greek letter'
;
;     To create an Angstrom sign::
;  
;        IDL> cgPlot, findgen(11), XTitle='This title contains (' + $
;                cgSymbol('Angstrom') +  ') an Angstrom sign.' 
;        IDL> cgPlot, findgen(11), XTitle='This title contains ($\Angstrom$) an Angstrom sign.'
;       
;     To create subscripts or superscripts::
;  
;        IDL> cgPlot, findgen(11), XTitle='E=mc$\up2$' 
;        IDL> cgPlot, findgen(11), XTitle='H$\sub2$O'
;        IDL> cgPlot, findgen(11), XTitle='H$\upSuper$MT $\Omega$$\subSubscript$', Charsize=2.0
;        
;     It is possible to use Greek characters as superscripts and subscripts. Do so by
;     prepending the Greek character with "\\" inside the normal superscript or subscript
;     notation. For example, to use lambda as an exponent to the Greek character Omega, you
;     can type this::
;
;        IDL> cgPlot, cgDemoData(1), XTitle='$\Omega$$\exp\\lambda$', Charsize=2.0
;
;     To use lambda as a subscript, type this:
;
;         IDL> cgPlot, cgDemoData(1), XTitle='$\Omega$$\sub\\lambda$', Charsize=2.0
;
; :Author:
;       FANNING SOFTWARE CONSULTING::
;           David W. Fanning 
;           1645 Sheely Drive
;           Fort Collins, CO 80526 USA
;           Phone: 970-221-0438
;           E-mail: david@idlcoyote.com
;           Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Written by: David W. Fanning, 2 September 2011. 
;        Added plus-minus symbol. 2 Nov 2011. DWF.
;        Added "up", "down", "exp" "sub" and "n" symbols for subscripting and superscripting. 9 Nov 2012. DWF.
;        Added "division" and "times" signs. 24 Nov 2012. DWF.
;        Updated UNICODE values to display capital letters correctly. 23 Dec 2012. DWF.
;
; :Copyright:
;     Copyright (c) 2011, Fanning Software Consulting, Inc.
;-
;+
; Displays the symbols and their names in a graphics window.
; 
; :Keywords:
;    ps: in, optional, type=boolean, default=0                                
;        Normally, the PostScript version of the symbol is returned automatically if
;        the current device is PostScript and !P.Font is 0 or 1. But, the PostScript version 
;        of the symbol can be obtained at any time and in any device, by setting this keyword.
;    unicode: in, optional, type=boolean, default=0                                
;        If this keyword is set, the function returns the Unicode value for the symbol.
; 
;-
PRO cgSymbol_Example, PS=ps, UNICODE=unicode

    Forward_Function cgSymbol

    Compile_Opt hidden
    
    ; The symbols..
    symbol = [  'alpha', 'beta', 'gamma', 'delta', 'epsilon',  'zeta', $
                'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu', $
                'nu', 'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau', $
                'upsilon', 'phi', 'chi', 'psi', 'omega', 'leq', 'geq', $
                'neq', 'deg', '+-', 'equiv', 'prime', 'angstrom', 'sun', $
                'varphi', 'infinity', 'copyright' ]
    
    ; Output positions.
    x = [0.15, 0.4, 0.65]
    y = Reverse((Indgen(13) + 1) * (1.0 / 14))
    
    ; Need PostScript output?
    IF Keyword_Set(ps) || (!D.Name EQ 'PS') THEN PS_Start
    
    ; Create a window.
    cgDisplay, 800, 500
    
    ; Output the symbols.
    FOR j=0,11 DO BEGIN
        cgText, x[0], y[j], symbol[j] + ': ' + $
            cgSymbol(symbol[j], UNICODE=unicode) + cgSymbol(symbol[j], /CAPITAL, UNICODE=unicode), $
            /NORMAL, CHARSIZE=1.5
        cgText, x[1], y[j], symbol[j+12] + ': ' + $
            cgSymbol(symbol[j+12], UNICODE=unicode) + cgSymbol(symbol[j+12], /CAPITAL, UNICODE=unicode), $
            /NORMAL, CHARSIZE=1.5
        IF j NE 12 THEN cgText, x[2], y[j], symbol[j+24] + ': ' + $
            cgSymbol(symbol[j+24], UNICODE=unicode), $
            /NORMAL, CHARSIZE=1.5
    ENDFOR
    
    ; Clean up PostScript, if needed.
    IF Keyword_Set(ps) || (!D.Name EQ 'PS') THEN PS_End
    
    ; Restore the users window.
    IF N_Elements(thisWindow) NE 0 THEN BEGIN
       IF thisWindow GE 0 THEN WSet, thisWindow
    ENDIF
    
END ; --------------------------------------------------------------------------------------


;+
; This function provides a device-independent way to ask for commonly-used
; symbols (e.g., less than or equal to, Angstrom, degree symbol, etc.),
; including all 24 Greek symbols in upper and lower case, to display with text.
;
; :Params:
;    symbol: in, required, type=string, default="alpha"
;       The name of the symbol desired as a string. Valid string names are the 24 
;       characters of the Greek alphabet, plus commonly used graphical symbols::
;       
;          alpha        nu        leq
;          beta         xi        geq
;          gamma        omicron   neg
;          delta        pi        deg
;          epsilon      rho       equiv
;          zeta         sigma     prime
;          eta          tau       angstrom
;          theta        upsilon   sun
;          iota         phi       varphi
;          kappa        chi       infinity
;          lambda       psi       copyright
;          mu           omega
;                    
;       Note that if the first letter of the name is capitalized, this is
;       the equivalent of setting the `Capital` keyword. Only Greek characters
;       use this method of selecting symbols.
;       
; :Keywords:
;    capital: in, optional, type=boolean, default=0
;        If this keyword is set, the captial Greek letter is returned rather than the lowercase 
;        Greek letter. An alternative way of capitalizing the letter is to make the first letter 
;        of the name an uppercase letter.
;    example: in, optional, type=boolean, default=0             
;        If this keyword is set, the names of symbols and the symbols themselves are written 
;        out in the current graphics window.
;    ps: in, optional, type=boolean, default=0                                
;        Normally, the PostScript version of the symbol is returned automatically if
;        the current device is PostScript and !P.Font is 0 or 1. But, the PostScript version 
;        of the symbol can be obtained at any time and in any device, by setting this keyword.
;    unicode: in, optional, type=boolean, default=0                                
;        If this keyword is set, the function returns the Unicode value for the symbol.
;          
;-
FUNCTION cgSymbol, symbol, CAPITAL=capital, EXAMPLE=example, PS=PS, UNICODE=unicode

    Compile_Opt idl2
    
    ; Return to caller on error.
    ON_Error, 2
    
    ; A common block to communicate with PS_Start.
    COMMON _$FSC_PS_START_, ps_struct
    
    ; Do you wish to see an example?
    IF Keyword_Set(example) THEN BEGIN
        cgSymbol_Example, UNICODE=unicode, PS=ps
        RETURN, ""
    ENDIF

    ; A symbol name is required.
    IF N_Elements(symbol) EQ 0 THEN BEGIN
       Print, 'Syntax: cgSymbol("theSymbol")
       RETURN, ""
    ENDIF
    
    ; Return quietly if the symbol is a null string.
    IF symbol EQ "" THEN RETURN, ""
    
    ; Set up PostScript device for working with Greek letters.
    IF !D.Name EQ 'PS' THEN Device, ISOLATIN1=1
    
    ; Check keywords.
    capital = Keyword_Set(capital)
    
    ; If the first letter of the "symbol" variable is uppercase
    ; the user wants an uppercase symbol, if available.
    firstLetter = StrMid(symbol, 0, 1)
    IF firstLetter EQ StrUpCase(firstLetter) THEN capital = 1
    
    IF Keyword_Set(unicode) THEN BEGIN

        CASE StrLowCase(symbol) OF
            'alpha':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'beta':    retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'gamma':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'delta':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'epsilon': retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'zeta':    retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'eta':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'theta':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'iota':    retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'kappa':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'lambda':  retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'mu':      retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'nu':      retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'xi':      retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'omicron': retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'pi':      retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'rho':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'sigma':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'tau':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'upsilon': retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'phi':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'chi':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'psi':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'omega':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'leq':     retSymbol = '!Z(2266)'
            'geq':     retSymbol = '!Z(2267)'
            'neq':     retSymbol = '!Z(2260)'
            'deg':     retSymbol = '!Z(00B0)'
            '+-':      retSymbol = '!Z(00B1)'
            'equiv':   retSymbol = '!Z(2261)'
            'prime':   retSymbol = '!Z(2232)'
            'angstrom':retSymbol = '!Z(00C5)'
            'sun':     retSymbol = '!Z(2609)'
            'varphi':  retSymbol = '!Z(03D5)'
            'infinity':retSymbol = '!Z(221E)'
            'copyright': retSymbol = '!Z(00A9)'
            'up':      retSymbol = '!U'
            'down':    retSymbol = '!D'
            'exp':     retSymbol = '!E'
            'sub':     retSymbol = '!I'
            'n':       retSymbol = '!N'
            'div':     retSymbol = '!Z(00F7)'
            'times':   retSymbol = '!Z(00D7)'
            ELSE: Message, 'The symbol ' + symbol + ' is unrecognized.'
       ENDCASE
           
       ; Return the UNICODE letter.
       RETURN, retSymbol
       
    ENDIF
    
    IF ((!D.Name EQ 'PS') && (!P.FONT NE -1)) || Keyword_Set(PS)  THEN BEGIN
    
        ; Make sure ISOLATIN1 encoding is turned on.
        IF (!D.Name EQ 'PS') THEN DEVICE, /ISOLATIN1
        
        CASE StrLowCase(symbol) OF
            'alpha':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'beta':    retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'gamma':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'delta':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'epsilon': retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'zeta':    retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'eta':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'theta':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'iota':    retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'kappa':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'lambda':  retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'mu':      retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'nu':      retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'xi':      retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'omicron': retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'pi':      retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'rho':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'sigma':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'tau':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'upsilon': retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'phi':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'chi':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'psi':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'omega':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'leq':     retSymbol = '!9' + String("243B) + '!X'
            'geq':     retSymbol = '!9' + String("263B) + '!X'
            'neq':     retSymbol = '!9' + String("271B) + '!X' 
            'deg':     retSymbol = '!9' + String("260B) + '!X'
            '+-':      retSymbol = '!9' + String("261B) + '!X'
            'equiv':   retSymbol = '!9' + String("272B) + '!X'
            'prime':   retSymbol = '!9' + String("242B) + '!X'
            'angstrom':retSymbol = '!3' + String("305B) + '!X'
            'sun':     BEGIN
                       IF (ps_struct.font EQ 1) && (StrUpCase(ps_struct.tt_font) EQ 'DEJAVUSANS') THEN BEGIN
                          retSymbol = '!Z(2609)'
                       ENDIF ELSE BEGIN
                         thisDevice = !D.Name
                         Set_Plot, 'PS'
                         Device, /AVANTGARDE, ISOLATIN1=0, /BOOK, FONT_INDEX = 20
                         retSymbol = '!20!S!DO!R!I ' + string(183b) + '!X!N'
                         Set_Plot, thisDevice
                       ENDELSE
                       END
            'varphi':  retSymbol = '!9' + String("152B) + '!X'
            'infinity':retSymbol = '!9' + String("245B) + '!X'
            'copyright':retSymbol = '!9' + String("323B) + '!X'
            'up':      retSymbol = '!U'
            'down':    retSymbol = '!D'
            'exp':     retSymbol = '!E'
            'sub':     retSymbol = '!I'
            'n':       retSymbol = '!N'
            'div':     retSymbol = '!9' + String("270B) + '!X'
            'times':   retSymbol = '!9' + String("264B) + '!X'
           ELSE: Message, 'The symbol ' + symbol + ' is unrecognized.'
        ENDCASE
    
    ENDIF ELSE BEGIN
    
        CASE StrLowCase(symbol) OF
            'alpha':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'beta':    retSymbol = cgGreek(symbol, CAPITAL=capital)
            'gamma':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'delta':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'epsilon': retSymbol = cgGreek(symbol, CAPITAL=capital)
            'zeta':    retSymbol = cgGreek(symbol, CAPITAL=capital)
            'eta':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'theta':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'iota':    retSymbol = cgGreek(symbol, CAPITAL=capital)
            'kappa':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'lambda':  retSymbol = cgGreek(symbol, CAPITAL=capital)
            'mu':      retSymbol = cgGreek(symbol, CAPITAL=capital)
            'nu':      retSymbol = cgGreek(symbol, CAPITAL=capital)
            'xi':      retSymbol = cgGreek(symbol, CAPITAL=capital)
            'omicron': retSymbol = cgGreek(symbol, CAPITAL=capital)
            'pi':      retSymbol = cgGreek(symbol, CAPITAL=capital)
            'rho':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'sigma':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'tau':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'upsilon': retSymbol = cgGreek(symbol, CAPITAL=capital)
            'phi':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'chi':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'psi':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'omega':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'leq':     retSymbol = '!9' + String("154B) + '!X'
            'geq':     retSymbol = '!9' + String("142B) + '!X'
            'neq':     retSymbol = '!9' + String("75B) + '!X' 
            'deg':     retSymbol = '!9' + String("45B) + '!X'
            '+-':      retSymbol = '!9' + String("53B) + '!X'
            'equiv':   retSymbol = '!9' + String("72B) + '!X'
            'prime':   retSymbol = '!9' + String("140B) + '!X'
            'angstrom':retSymbol = '!3' + String("305B) + '!X'
            'sun':     retSymbol = '!D!9n!N!X'
            'varphi':  retSymbol = '!9' + String("120B) + '!X'
            'infinity':retSymbol = '!9' + String("44B) + '!X'
            'copyright':retSymbol = '!3' + String("251B) + '!X'
            'up':      retSymbol = '!U'
            'down':    retSymbol = '!D'
            'exp':     retSymbol = '!E'
            'sub':     retSymbol = '!I'
            'n':       retSymbol = '!N'
            'div':     retSymbol = '!9' + String("57B) + '!X'
            'times':   retSymbol = '!9' + String("130B) + '!X'
            ELSE: Message, 'The symbol ' + symbol + ' is unrecognized.'      
        ENDCASE
        
      ENDELSE

    RETURN, retSymbol
    
END
;;==============================================================================
; docformat = 'rst'
;
; NAME:
;   cgSymbol
;
; PURPOSE:
; This function provides a device-independent way to ask for commonly-used
; symbols (e.g., less than or equal to, Angstrom, degree symbol, etc.),
; including all 24 Greek symbols in upper and lower case, to display with text.
; 
;******************************************************************************************;
;                                                                                          ;
;  Copyright (c) 2011, by Fanning Software Consulting, Inc. All rights reserved.           ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;
;+
; This function provides a device-independent way to ask for commonly-used
; symbols (e.g., less than or equal to, Angstrom, degree symbol, subscripts, superscripts, etc.),
; including all 24 Greek symbols in upper and lower case, to display with text.
;
; For this program to work correctly on your graphics display, you should be using
; Hershey fonts (!P.Font=-1). It will work correctly in PostScript with either 
; hardware fonts (!P.Font=0) or True-Type fonts (!P.Font=1).
; 
; `Greek Symbols <http://www.idlcoyote.com/ps_tips/greeksym.html>` are created by
; calling the Coyote Library routine `cgGreek <http://www.idlcoyote.com/programs/cggreek.pro>' 
; from this program.
; 
; .. image:: cgsymbol.png 
; 
; Normally, rather than calling cgSymbol, the symbols are embedded in Coyote Graphics
; text that are used for axis annotation and so forth. See `Embedding Symbols in Coyote
; Graphics Output <http://www.idlcoyote.com/cg_tips/embedsymbols.php>`. Embedded subscripts
; and superscripts are implemented like this::
;     $\upTEXT$ : Raise the text in TEXT to superscript level and reduce character size by 62%.
;     $\expTEXT$ : Raise the text in TEXT to superscript level and reduce character size by 44%.
;     $\downTEXT$ : Lower the text in TEXT to subscript level and reduce character size by 62%.
;     $\subTEXT$ : Lower the text in TEXT to subscript level and reduce character size by 44%.
;
; :Categories:
;    Graphics
;    
; :Returns:
;    A string variable that represents the requested symbol and can be used
;    in a textual context.
;    
; :Examples:
;     To create a lowercase Greek psi symbol::
;  
;        IDL> cgPlot, findgen(11), XTitle='This title contains ' + $
;                 cgSymbol('psi') + ' as a Greek letter' 
;        IDL> cgPlot, findgen(11), XTitle='This title contains $\psi$ as a Greek letter'
;
;     To create an Angstrom sign::
;  
;        IDL> cgPlot, findgen(11), XTitle='This title contains (' + $
;                cgSymbol('Angstrom') +  ') an Angstrom sign.' 
;        IDL> cgPlot, findgen(11), XTitle='This title contains ($\Angstrom$) an Angstrom sign.'
;       
;     To create subscripts or superscripts::
;  
;        IDL> cgPlot, findgen(11), XTitle='E=mc$\up2$' 
;        IDL> cgPlot, findgen(11), XTitle='H$\sub2$O'
;        IDL> cgPlot, findgen(11), XTitle='H$\upSuper$MT $\Omega$$\subSubscript$', Charsize=2.0
;        
;     It is possible to use Greek characters as superscripts and subscripts. Do so by
;     prepending the Greek character with "\\" inside the normal superscript or subscript
;     notation. For example, to use lambda as an exponent to the Greek character Omega, you
;     can type this::
;
;        IDL> cgPlot, cgDemoData(1), XTitle='$\Omega$$\exp\\lambda$', Charsize=2.0
;
;     To use lambda as a subscript, type this:
;
;         IDL> cgPlot, cgDemoData(1), XTitle='$\Omega$$\sub\\lambda$', Charsize=2.0
;
; :Author:
;       FANNING SOFTWARE CONSULTING::
;           David W. Fanning 
;           1645 Sheely Drive
;           Fort Collins, CO 80526 USA
;           Phone: 970-221-0438
;           E-mail: david@idlcoyote.com
;           Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Written by: David W. Fanning, 2 September 2011. 
;        Added plus-minus symbol. 2 Nov 2011. DWF.
;        Added "up", "down", "exp" "sub" and "n" symbols for subscripting and superscripting. 9 Nov 2012. DWF.
;        Added "division" and "times" signs. 24 Nov 2012. DWF.
;        Updated UNICODE values to display capital letters correctly. 23 Dec 2012. DWF.
;
; :Copyright:
;     Copyright (c) 2011, Fanning Software Consulting, Inc.
;-
;+
; Displays the symbols and their names in a graphics window.
; 
; :Keywords:
;    ps: in, optional, type=boolean, default=0                                
;        Normally, the PostScript version of the symbol is returned automatically if
;        the current device is PostScript and !P.Font is 0 or 1. But, the PostScript version 
;        of the symbol can be obtained at any time and in any device, by setting this keyword.
;    unicode: in, optional, type=boolean, default=0                                
;        If this keyword is set, the function returns the Unicode value for the symbol.
; 
;-
PRO cgSymbol_Example, PS=ps, UNICODE=unicode

    Forward_Function cgSymbol

    Compile_Opt hidden
    
    ; The symbols..
    symbol = [  'alpha', 'beta', 'gamma', 'delta', 'epsilon',  'zeta', $
                'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu', $
                'nu', 'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau', $
                'upsilon', 'phi', 'chi', 'psi', 'omega', 'leq', 'geq', $
                'neq', 'deg', '+-', 'equiv', 'prime', 'angstrom', 'sun', $
                'varphi', 'infinity', 'copyright' ]
    
    ; Output positions.
    x = [0.15, 0.4, 0.65]
    y = Reverse((Indgen(13) + 1) * (1.0 / 14))
    
    ; Need PostScript output?
    IF Keyword_Set(ps) || (!D.Name EQ 'PS') THEN PS_Start
    
    ; Create a window.
    cgDisplay, 800, 500
    
    ; Output the symbols.
    FOR j=0,11 DO BEGIN
        cgText, x[0], y[j], symbol[j] + ': ' + $
            cgSymbol(symbol[j], UNICODE=unicode) + cgSymbol(symbol[j], /CAPITAL, UNICODE=unicode), $
            /NORMAL, CHARSIZE=1.5
        cgText, x[1], y[j], symbol[j+12] + ': ' + $
            cgSymbol(symbol[j+12], UNICODE=unicode) + cgSymbol(symbol[j+12], /CAPITAL, UNICODE=unicode), $
            /NORMAL, CHARSIZE=1.5
        IF j NE 12 THEN cgText, x[2], y[j], symbol[j+24] + ': ' + $
            cgSymbol(symbol[j+24], UNICODE=unicode), $
            /NORMAL, CHARSIZE=1.5
    ENDFOR
    
    ; Clean up PostScript, if needed.
    IF Keyword_Set(ps) || (!D.Name EQ 'PS') THEN PS_End
    
    ; Restore the users window.
    IF N_Elements(thisWindow) NE 0 THEN BEGIN
       IF thisWindow GE 0 THEN WSet, thisWindow
    ENDIF
    
END ; --------------------------------------------------------------------------------------


;+
; This function provides a device-independent way to ask for commonly-used
; symbols (e.g., less than or equal to, Angstrom, degree symbol, etc.),
; including all 24 Greek symbols in upper and lower case, to display with text.
;
; :Params:
;    symbol: in, required, type=string, default="alpha"
;       The name of the symbol desired as a string. Valid string names are the 24 
;       characters of the Greek alphabet, plus commonly used graphical symbols::
;       
;          alpha        nu        leq
;          beta         xi        geq
;          gamma        omicron   neg
;          delta        pi        deg
;          epsilon      rho       equiv
;          zeta         sigma     prime
;          eta          tau       angstrom
;          theta        upsilon   sun
;          iota         phi       varphi
;          kappa        chi       infinity
;          lambda       psi       copyright
;          mu           omega
;                    
;       Note that if the first letter of the name is capitalized, this is
;       the equivalent of setting the `Capital` keyword. Only Greek characters
;       use this method of selecting symbols.
;       
; :Keywords:
;    capital: in, optional, type=boolean, default=0
;        If this keyword is set, the captial Greek letter is returned rather than the lowercase 
;        Greek letter. An alternative way of capitalizing the letter is to make the first letter 
;        of the name an uppercase letter.
;    example: in, optional, type=boolean, default=0             
;        If this keyword is set, the names of symbols and the symbols themselves are written 
;        out in the current graphics window.
;    ps: in, optional, type=boolean, default=0                                
;        Normally, the PostScript version of the symbol is returned automatically if
;        the current device is PostScript and !P.Font is 0 or 1. But, the PostScript version 
;        of the symbol can be obtained at any time and in any device, by setting this keyword.
;    unicode: in, optional, type=boolean, default=0                                
;        If this keyword is set, the function returns the Unicode value for the symbol.
;          
;-
FUNCTION cgSymbol, symbol, CAPITAL=capital, EXAMPLE=example, PS=PS, UNICODE=unicode

    Compile_Opt idl2
    
    ; Return to caller on error.
    ON_Error, 2
    
    ; A common block to communicate with PS_Start.
    COMMON _$FSC_PS_START_, ps_struct
    
    ; Do you wish to see an example?
    IF Keyword_Set(example) THEN BEGIN
        cgSymbol_Example, UNICODE=unicode, PS=ps
        RETURN, ""
    ENDIF

    ; A symbol name is required.
    IF N_Elements(symbol) EQ 0 THEN BEGIN
       Print, 'Syntax: cgSymbol("theSymbol")
       RETURN, ""
    ENDIF
    
    ; Return quietly if the symbol is a null string.
    IF symbol EQ "" THEN RETURN, ""
    
    ; Set up PostScript device for working with Greek letters.
    IF !D.Name EQ 'PS' THEN Device, ISOLATIN1=1
    
    ; Check keywords.
    capital = Keyword_Set(capital)
    
    ; If the first letter of the "symbol" variable is uppercase
    ; the user wants an uppercase symbol, if available.
    firstLetter = StrMid(symbol, 0, 1)
    IF firstLetter EQ StrUpCase(firstLetter) THEN capital = 1
    
    IF Keyword_Set(unicode) THEN BEGIN

        CASE StrLowCase(symbol) OF
            'alpha':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'beta':    retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'gamma':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'delta':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'epsilon': retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'zeta':    retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'eta':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'theta':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'iota':    retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'kappa':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'lambda':  retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'mu':      retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'nu':      retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'xi':      retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'omicron': retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'pi':      retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'rho':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'sigma':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'tau':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'upsilon': retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'phi':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'chi':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'psi':     retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'omega':   retSymbol = cgGreek(symbol, /UNICODE, CAPITAL=capital)
            'leq':     retSymbol = '!Z(2266)'
            'geq':     retSymbol = '!Z(2267)'
            'neq':     retSymbol = '!Z(2260)'
            'deg':     retSymbol = '!Z(00B0)'
            '+-':      retSymbol = '!Z(00B1)'
            'equiv':   retSymbol = '!Z(2261)'
            'prime':   retSymbol = '!Z(2232)'
            'angstrom':retSymbol = '!Z(00C5)'
            'sun':     retSymbol = '!Z(2609)'
            'varphi':  retSymbol = '!Z(03D5)'
            'infinity':retSymbol = '!Z(221E)'
            'copyright': retSymbol = '!Z(00A9)'
            'up':      retSymbol = '!U'
            'down':    retSymbol = '!D'
            'exp':     retSymbol = '!E'
            'sub':     retSymbol = '!I'
            'n':       retSymbol = '!N'
            'div':     retSymbol = '!Z(00F7)'
            'times':   retSymbol = '!Z(00D7)'
            ELSE: Message, 'The symbol ' + symbol + ' is unrecognized.'
       ENDCASE
           
       ; Return the UNICODE letter.
       RETURN, retSymbol
       
    ENDIF
    
    IF ((!D.Name EQ 'PS') && (!P.FONT NE -1)) || Keyword_Set(PS)  THEN BEGIN
    
        ; Make sure ISOLATIN1 encoding is turned on.
        IF (!D.Name EQ 'PS') THEN DEVICE, /ISOLATIN1
        
        CASE StrLowCase(symbol) OF
            'alpha':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'beta':    retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'gamma':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'delta':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'epsilon': retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'zeta':    retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'eta':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'theta':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'iota':    retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'kappa':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'lambda':  retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'mu':      retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'nu':      retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'xi':      retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'omicron': retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'pi':      retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'rho':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'sigma':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'tau':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'upsilon': retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'phi':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'chi':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'psi':     retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'omega':   retSymbol = cgGreek(symbol, CAPITAL=capital, /PS)
            'leq':     retSymbol = '!9' + String("243B) + '!X'
            'geq':     retSymbol = '!9' + String("263B) + '!X'
            'neq':     retSymbol = '!9' + String("271B) + '!X' 
            'deg':     retSymbol = '!9' + String("260B) + '!X'
            '+-':      retSymbol = '!9' + String("261B) + '!X'
            'equiv':   retSymbol = '!9' + String("272B) + '!X'
            'prime':   retSymbol = '!9' + String("242B) + '!X'
            'angstrom':retSymbol = '!3' + String("305B) + '!X'
            'sun':     BEGIN
                       IF (ps_struct.font EQ 1) && (StrUpCase(ps_struct.tt_font) EQ 'DEJAVUSANS') THEN BEGIN
                          retSymbol = '!Z(2609)'
                       ENDIF ELSE BEGIN
                         thisDevice = !D.Name
                         Set_Plot, 'PS'
                         Device, /AVANTGARDE, ISOLATIN1=0, /BOOK, FONT_INDEX = 20
                         retSymbol = '!20!S!DO!R!I ' + string(183b) + '!X!N'
                         Set_Plot, thisDevice
                       ENDELSE
                       END
            'varphi':  retSymbol = '!9' + String("152B) + '!X'
            'infinity':retSymbol = '!9' + String("245B) + '!X'
            'copyright':retSymbol = '!9' + String("323B) + '!X'
            'up':      retSymbol = '!U'
            'down':    retSymbol = '!D'
            'exp':     retSymbol = '!E'
            'sub':     retSymbol = '!I'
            'n':       retSymbol = '!N'
            'div':     retSymbol = '!9' + String("270B) + '!X'
            'times':   retSymbol = '!9' + String("264B) + '!X'
           ELSE: Message, 'The symbol ' + symbol + ' is unrecognized.'
        ENDCASE
    
    ENDIF ELSE BEGIN
    
        CASE StrLowCase(symbol) OF
            'alpha':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'beta':    retSymbol = cgGreek(symbol, CAPITAL=capital)
            'gamma':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'delta':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'epsilon': retSymbol = cgGreek(symbol, CAPITAL=capital)
            'zeta':    retSymbol = cgGreek(symbol, CAPITAL=capital)
            'eta':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'theta':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'iota':    retSymbol = cgGreek(symbol, CAPITAL=capital)
            'kappa':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'lambda':  retSymbol = cgGreek(symbol, CAPITAL=capital)
            'mu':      retSymbol = cgGreek(symbol, CAPITAL=capital)
            'nu':      retSymbol = cgGreek(symbol, CAPITAL=capital)
            'xi':      retSymbol = cgGreek(symbol, CAPITAL=capital)
            'omicron': retSymbol = cgGreek(symbol, CAPITAL=capital)
            'pi':      retSymbol = cgGreek(symbol, CAPITAL=capital)
            'rho':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'sigma':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'tau':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'upsilon': retSymbol = cgGreek(symbol, CAPITAL=capital)
            'phi':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'chi':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'psi':     retSymbol = cgGreek(symbol, CAPITAL=capital)
            'omega':   retSymbol = cgGreek(symbol, CAPITAL=capital)
            'leq':     retSymbol = '!9' + String("154B) + '!X'
            'geq':     retSymbol = '!9' + String("142B) + '!X'
            'neq':     retSymbol = '!9' + String("75B) + '!X' 
            'deg':     retSymbol = '!9' + String("45B) + '!X'
            '+-':      retSymbol = '!9' + String("53B) + '!X'
            'equiv':   retSymbol = '!9' + String("72B) + '!X'
            'prime':   retSymbol = '!9' + String("140B) + '!X'
            'angstrom':retSymbol = '!3' + String("305B) + '!X'
            'sun':     retSymbol = '!D!9n!N!X'
            'varphi':  retSymbol = '!9' + String("120B) + '!X'
            'infinity':retSymbol = '!9' + String("44B) + '!X'
            'copyright':retSymbol = '!3' + String("251B) + '!X'
            'up':      retSymbol = '!U'
            'down':    retSymbol = '!D'
            'exp':     retSymbol = '!E'
            'sub':     retSymbol = '!I'
            'n':       retSymbol = '!N'
            'div':     retSymbol = '!9' + String("57B) + '!X'
            'times':   retSymbol = '!9' + String("130B) + '!X'
            ELSE: Message, 'The symbol ' + symbol + ' is unrecognized.'      
        ENDCASE
        
      ENDELSE

    RETURN, retSymbol
    
END
;;==============================================================================
; docformat = 'rst'
;
; NAME:
;   ColorsAreIdentical
;
; PURPOSE:
;   Returns a 1 if the two input colors refer to the same color, otherwise returns a 0.
;
;******************************************************************************************;
;                                                                                          ;
;  Copyright (c) 2010, by Fanning Software Consulting, Inc. All rights reserved.           ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;
;+
; :Description:
;   Returns a 1 if the two input colors refer to the same color, otherwise returns a 0.
;
; :Categories:
;    Graphics Utility
;    
; :Params:
;    color_1: in, required, type=string/integer/long
;         The first color to compare for "equality".
;    color_2: in, required, type=string/integer/long
;         The second color to compare for "equality".
;       
; :Keywords:
;     None.
;          
; :Examples:
;    Used to compare if two different colors are the same color::
;       IDL> Print, ColorsAreIdentical('white', cgColor('white'))
;       IDL> Print, ColorsAreIdentical(252, !P.Color)
;       IDL> Print, ColorsAreIdentical('white', '255')
;       
; :Author:
;       FANNING SOFTWARE CONSULTING::
;           David W. Fanning 
;           1645 Sheely Drive
;           Fort Collins, CO 80526 USA
;           Phone: 970-221-0438
;           E-mail: david@idlcoyote.com
;           Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Written, 24 December 2010. DWF.
;        Fixed a typo when first color is INTEGER and second color is STRING. 3 Jan 2011. DWF.
;        Added error handling for out of bounds color values. 25 May 2011. DWF.
;
; :Copyright:
;     Copyright (c) 2010, Fanning Software Consulting, Inc.
;-
FUNCTION ColorsAreIdentical, color_1, color_2

    On_Error, 2

    answer = 0
    IF N_Params() NE 2 THEN Message, 'Two positional parameters (colors) are required.'
    IF N_Elements(color_1) EQ 0 THEN Message, 'Color 1 is undefined.'
    IF N_Elements(color_2) EQ 0 THEN Message, 'Color 2 is undefined.'
    
    ; If the color is a "number" string, turn it back into a number.
    ; If all the bytes are between 48 and 57, then this string is a "number".
    IF Size(color_1, /TNAME) EQ 'STRING' THEN BEGIN
        bytecheck = Byte(StrUpCase(color_1))
        i = Where(bytecheck LT 48, lessthan)
        i = Where(bytecheck GT 57, greaterthan)
        IF (lessthan + greaterthan) EQ 0 THEN c1 = Fix(color_1)
    ENDIF
    IF N_Elements(c1) EQ 0 THEN c1 = color_1
    
    ; Colors have to be between 0 and 255.
    IF Size(c1, /TYPE) LE 2 THEN BEGIN
        IF (c1 LT 0) || (c1 GT 255) THEN BEGIN
            msg = 'Color value of ' + StrTrim(c1,2) + ' is outside expected color range of 0 to 255.'
            Message, msg
        ENDIF
    ENDIF
    
    IF Size(color_2, /TNAME) EQ 'STRING' THEN BEGIN
        bytecheck = Byte(StrUpCase(color_2))
        i = Where(bytecheck LT 48, lessthan)
        i = Where(bytecheck GT 57, greaterthan)
        IF (lessthan + greaterthan) EQ 0 THEN c2 = Fix(color_2)
    ENDIF
    IF N_Elements(c2) EQ 0 THEN c2 = color_2
    
    ; Colors have to be between 0 and 255.
    IF Size(c2, /TYPE) LE 2 THEN BEGIN
        IF (c2 LT 0) || (c2 GT 255) THEN BEGIN
            msg = 'Color value of ' + StrTrim(c2,2) + ' is outside expected color range of 0 to 255.'
            Message, msg
        ENDIF
    ENDIF
    
    ; If the colors are the same type, compare them directly
    IF Size(c1, /TYPE) EQ Size(c2, /TYPE) THEN BEGIN
        IF Size(c1, /TNAME) EQ 'STRING' THEN BEGIN
            IF StrUpCase(c1) EQ StrUpCase(c2) THEN answer = 1 ELSE answer = 0
        ENDIF ELSE BEGIN
            IF c1 EQ c2 THEN answer = 1 ELSE answer = 0
        ENDELSE
     ENDIF ELSE BEGIN
    
        ; Get colors so we can restore later.
        TVLCT, r, g, b, /GET
        
        ; Check different purtubations of data type.
        CASE 1 OF
        
           ; First color a STRING, second color an INTEGER.
           (Size(c1, /TYPE) EQ 7) AND (Size(c2, /TYPE) LE 2): BEGIN
                v1 = Transpose(cgColor(c1, /TRIPLE))
                v2 = [r[c2], g[c2], b[c2]]
                answer = Array_Equal(v1, v2)
                END
 
           ; First color a STRING, second color a LONG.
           (Size(c1, /TYPE) EQ 7) AND (Size(c2, /TYPE) EQ 3): BEGIN
                answer = Array_Equal(cgColor(c1, /DECOMPOSED), c2)
                END
             
           ; First color an INTEGER, second color a STRING.    
           (Size(c1, /TYPE) LE 2) AND (Size(c2, /TYPE) EQ 7): BEGIN
                v1 = Transpose(cgColor(c2, /TRIPLE))
                v2 = [r[c1], g[c1], b[c1]]
                answer = Array_Equal(v1, v2)
                END

           ; First color a LONG, second color a STRING.
           (Size(c1, /TYPE) EQ 3) AND (Size(c2, /TYPE) EQ 7): BEGIN
                answer = Array_Equal(cgColor(c2, /DECOMPOSED), c1)
                END

            ; First color a LONG and second color an INTEGER.
           (Size(c1, /TYPE) EQ 3) AND (Size(c2, /TYPE) LE 2): BEGIN
                answer = Array_Equal(c1, cgColor24([r[c2], g[c2], b[c2]]))
                END

           ; First color an INTEGER, second color a LONG.
           (Size(c1, /TYPE) LE 2) AND (Size(c2, /TYPE) EQ 3): BEGIN
                answer = Array_Equal(c2, cgColor24([r[c1], g[c1], b[c1]]))
                END

           ELSE: Message, 'Colors do not meet type expectations. Unsure how to proceed.'
        
        ENDCASE
        
        ; Restore the color vectors, in case they were changed.
        TVLCT, r, g, b
    ENDELSE

    RETURN, answer
END
;;==============================================================================
;+
; NAME:
;  DECOMPOSEDCOLOR
;
; PURPOSE:
;
;   This function is used to determine, in a device independent way, if the 
;   current graphics device is using color decomposition. The function returns
;   a 1 if color decomposition is turned on, and a 0 if it is turned off. When
;   color decomposition is turned on, we say the device is using a true-color
;   display. If color decomposition is turned off, we say the device is using
;   an indexed color display.
;
; AUTHOR:
;
;   FANNING SOFTWARE CONSULTING
;   David Fanning, Ph.D.
;   1645 Sheely Drive
;   Fort Collins, CO 80526 USA
;   Phone: 970-221-0438
;   E-mail: david@idlcoyote.com
;   Coyote's Guide to IDL Programming: http://www.idlcoyote.com/
;
; CATEGORY:
;
;   Utilities
;
; CALLING SEQUENCE:
;
;   result = DecomposedColor()
;
; RETURN VALUE:
;
;   result:       A 1 if color decomposition is turned on. A 0 if color decomposition is turned off.
;
; ARGUMENTS:
;
;  device:        The IDL graphics device whose color decomposition state you wish to know the
;                 current value of. If undefined, the current graphics device is used.
;
; KEYWORDRS:
;
;  DEPTH:          An output keyword that returns the depth of the graphics device. Normally,
;                  either 8 for index color devices, with color decomposition turned off, or 24
;                  for true-color devices with color decomposition turned on.
;
; EXAMPLE:
;
;  IDL> Print, DecomposedColor()     ; Color decomposition state of current graphics device.
;       1
;  IDL> Print, DecomposedColor('PS') ; Color decomposition state of PostScript graphics device.
;       0
;
; MODIFICATION HISTORY:
;
;  Written by: David W. Fanning, May 24, 2009.
;  Modified the way decomposition was obtained for PostScript devices IDL 7.1 and higher. 12 Dec 2010. DWF.
;  Fixed a problem in the CASE statement with ELSE clause and added a NULL device segment. 4 Jan 2011. DWF.
;  It now appears 24-bit PostScript support was added in IDL 7.1, although the Get_Decomposed keyword
;      didn't work until IDL 7.1.1. 13 January 2011. DWF
;
;-
;******************************************************************************************;
;  Copyright (c) 2009, by Fanning Software Consulting, Inc.                                ;
;  All rights reserved.                                                                    ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
FUNCTION DecomposedColor, device, DEPTH=depth

    ; Return to caller on error.
    ON_ERROR, 2

    ; Was a graphics device passed in?
    IF N_Elements(device) EQ 0 THEN device = !D.NAME
    
    ; If the asked for graphics device is not the same as the current device,
    ; load the one the user asked for.
    IF StrUpCase(device) NE !D.NAME THEN BEGIN
        thisDevice = !D.NAME
        Set_Plot, device
    ENDIF

    ; Which graphics device are you interested in?
    CASE !D.NAME OF
    
        'PS': BEGIN ; PostScript
           CASE 1 OF
                Float(!Version.Release) EQ 7.1: BEGIN
                    Help, /DEVICE, OUTPUT=outstr
                    psinfo = outstr[4]
                    parts = StrSplit(psinfo, ':', /EXTRACT)
                    IF StrUpCase(StrCompress(parts[1], /REMOVE_ALL)) EQ 'DECOMPOSED' THEN BEGIN
                        decomposed = 1
                        depth = 24
                    ENDIF ELSE BEGIN
                        decomposed = 0
                        depth = 8
                    ENDELSE
                END
                Float(!Version.Release) GT 7.1: BEGIN
                    Device, GET_DECOMPOSED=decomposed
                    IF decomposed THEN depth = 24 ELSE depth = 8
                    END
                ELSE: BEGIN
                    decomposed = 0
                    depth = 8
                    END
            ENDCASE
           END
           
        'Z': BEGIN ; Z-graphics buffer.
            IF (Float(!Version.Release) GE 6.4) THEN BEGIN
                Device, GET_DECOMPOSED=decomposed
                Device, GET_PIXEL_DEPTH=depth
            ENDIF ELSE BEGIN
                decomposed = 0
                depth = 8
            ENDELSE
            END
            
        'X': Device, GET_DECOMPOSED=decomposed, GET_VISUAL_DEPTH=depth
        
        'WIN': Device, GET_DECOMPOSED=decomposed, GET_VISUAL_DEPTH=depth
        
        'MAC': BEGIN
            IF (Float(!Version.Release) GE 5.2) THEN BEGIN
                Device, Get_Decomposed=decomposedState, GET_VISUAL_DEPTH=depth
            ENDIF ELSE BEGIN
                decomposed = 0
                depth = 8
            ENDELSE
            END
            
         'NULL': BEGIN ; Setting up in decomposed mode will make sure
                       ; drawing colors are never loaded, which is not
                       ; allowed for the NULL device.
            decomposed = 1
            depth = 24
            END
            
        ELSE: BEGIN ; All other devices are 8-bit oldsters.
            decomposed = 0
            depth = 8
            END
    ENDCASE
 
    ; Need to clean up?
    IF N_Elements(thisDevice) NE 0 THEN Set_Plot, thisDevice
    
    ; Return the result.
    RETURN, decomposed
    
END
;;==============================================================================
; docformat = 'rst'
;
; NAME:
;   GetDecomposedState
;
; PURPOSE:
;   Provides a device-independent way to get the color decomposition state of the
;   current graphics device. 
;
;******************************************************************************************;
;                                                                                          ;
;  Copyright (c) 2010, by Fanning Software Consulting, Inc. All rights reserved.           ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;
;+
;  Provides a device-independent way to get the color decomposition state of the
;  current graphics device. 
;
; :Categories:
;    Graphics, Utilities
;    
; :Returns:
;     Returns a 1 if color decomposition is turned on and a 0 if indexed color is used.
;       
; :Keywords:
;     Depth: out, optional, type=integer
;         The depth of the color display. Typically 8 for indexed color devices
;         and 24 for true-color devices.
;          
; :Examples:
;       IDL> currentState = GetDecomposedState()
;       
; :Author:
;       FANNING SOFTWARE CONSULTING::
;           David W. Fanning 
;           1645 Sheely Drive
;           Fort Collins, CO 80526 USA
;           Phone: 970-221-0438
;           E-mail: david@idlcoyote.com
;           Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Written, 12 December 2010 as a better named wrapper for DECOMPOSEDCOLOR program. DWF.
;
; :Copyright:
;     Copyright (c) 2010, Fanning Software Consulting, Inc.
;-
FUNCTION GetDecomposedState, DEPTH=depth

    Compile_Opt idl2
    
    RETURN, DecomposedColor(Depth=depth)
    
END    
;;==============================================================================
; docformat = 'rst'
;
; NAME:
;   SetDecomposedState
;
; PURPOSE:
;   Provides a device-independent way to set the color decomposition state of the
;   current graphics device. Devices that do not have a DECOMPOSED keyword to the
;   DEVICE command are assumed to be in indexed color mode always.
;
;******************************************************************************************;
;                                                                                          ;
;  Copyright (c) 2010, by Fanning Software Consulting, Inc. All rights reserved.           ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;
;+
;   Provides a device-independent way to set the color decomposition state of the
;   current graphics device. Devices that do not have a DECOMPOSED keyword to the
;   DEVICE command are assumed to be in indexed color mode always.
;   
;   I have removed the Z-graphics buffer from being controlled by this program. I
;   do so reluctantly, but I am pragmatic enough to realize that progress forward
;   is necessarily slow and that I must recognize the reality of legacy IDL code.
;   
;   My personal view is that all graphics routines should use 24-bit decomposed
;   color. There are myriad advantages, but basically they boil down to this:
;   (1) You have 16.7 million colors available to you simultaneously, and (2) you
;   don't have to contaminate color tables to use drawing colors. Coyote Graphics
;   routines are in the habit of switching out of whatever color mode the user happens 
;   to be using, into 24-bit decomposed color mode, then switching back when finished
;   with its business. But, it is impossible to do this correctly in the Z-graphics
;   buffer.
;   
;   The reason for this is that in the Z-graphics buffer, you need to switch not only
;   the color mode, but also the pixel depth. In other words, I would prefer to set
;   the Z-graphics buffer up like this::
;   
;       Set_Plot, 'Z'
;       Device, Decomposed=1, Set_Pixel_Depth=24
;       
;   But, if I do that, then I need to set it back (for 99 people out of a 100) like this::
;   
;       Device, Decomposed=0, Set_Pixel_Depth=8
;       
;   Unfortunately, using this command will erase whatever is in the Z-graphics buffer!
;   
;   My solution to this problem is to leave it to the user to configure the Z-graphics buffer
;   the way they want it. If you configure it to use 24-bit decomposed color, all of the Coyote
;   Graphics routines will work as they normally do. If you configure it to use 8-bit indexed color,
;   which is the default configuration, then it will work "normally", but you will be in great
;   danger of contaminating the system color table. The reason for this is that Coyote Graphics
;   cannot "restore" the entry color table in the Z-buffer without obliterating what is already
;   in the graphics window. Users who work with indexed color are probably already very much
;   aware of this problem, so it shouldn't surprise them. (They might not expect this with
;   Coyote Graphics, but this is the price that has to be paid.)
;   
;   My suggestion is to put the Z-graphics configuration in your IDL startup file. Set it
;   up in 24-bit decomposed color mode, use Coyote Graphics to do all your graphical output,
;   and you will find things working perfectly. 
;   See `Configuring the Z-Graphics Buffer for Coyote Graphics <http://www.idlcoyote.com/cg_tips/configz.php>` 
;   for additional information.
;
; :Categories:
;    Graphics, Utilities
;    
; :Params:
;    state: in, required, type=integer, default=0
;         Set to 1 to set the current graphics device to decomposed color. Set to 0 
;         to set the current graphics device to indexed color. Devices lacking a 
;         DECOMPOSED keyword are assumed to be perpetually in indexed color mode.
;         The Z-graphics buffer is always unchanged after the 24 Dec 2011 update. 
;       
; :Keywords:
;     currentstate: out, optional, type=integer
;         The current decomposition state of the current graphics device when the
;         program is called. A 1 indicates decomposed color. A 0 indicates indexed 
;         color.
;     depth: out, optional, type=integer
;         The currnet pixel depth of the graphics device. 
;     zdepth: in, optional, type=integer
;         The pixel depth of the Z-graphics device. Set to 8 or 24. Applies ONLY 
;         when setting the Z-graphics device state to 0. If undefined, the current 
;         depth of the Z-graphics device is unchanged from its current state. No
;         longer used after 24 Dec 2011 update. Code still in place, however.
;          
; :Examples:
;     To set the device in color decomposition mode, then return it, do something like this::
;     
;        SetDecomposedState, 1, CurrentState=mode
;        ...
;        SetDecomposeState, mode
;       
; :Author:
;       FANNING SOFTWARE CONSULTING::
;           David W. Fanning 
;           1645 Sheely Drive
;           Fort Collins, CO 80526 USA
;           Phone: 970-221-0438
;           E-mail: david@idlcoyote.com
;           Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Written, 16 November 2010. DWF.
;        Changes to include SET_PIXEL_DEPTH in Z-graphics buffer. 19 Nov 2010. DWF.
;        Allow PostScript 7.0 to set the decomposition keyword. 12 Dec 2010. DWF.
;        Added DEPTH and ZDEPTH keywords. 31 Dec 2010. DWF.
;        Added a do-nothing NULL device to Case statement. 4 Jan 2011. DWF.
;        Removed the Z-graphics buffer from control by this program. 24 Dec 2011. DWF.
;        Added back the ability to set decomposed state for the Z-buffer, but only
;           if the depth buffer is 24-bits or higher. 25 May 2012. DWF.
;
; :Copyright:
;     Copyright (c) 2010, Fanning Software Consulting, Inc.
;-
PRO SetDecomposedState, state, CURRENTSTATE=currentState, DEPTH=depth, ZDEPTH=zdepth

    Compile_Opt idl2
    
    On_Error, 2
    
    ; Set to indexed color if you are not told differently.
    state = Keyword_Set(state)
    
    ; Get the current state and pixel depth.
    currentState = DecomposedColor(DEPTH=depth)
    
    ; Set the decomposition state, if you are able. Otherwise assume
    ; indexed color mode.
    CASE StrUpCase(!D.Name) OF
    
       'PS': IF (!Version.Release GE 7.1) THEN Device, Decomposed=state
                             
        'Z': BEGIN ; Z-Graphics Buffer ; See program notes for why the Z-buffer is no longer controlled here.
;             IF Float(!Version.Release) GE 6.4 THEN BEGIN
;                IF N_Elements(zdepth) NE 0 THEN BEGIN
;                   IF (zdepth NE 8) AND (zdepth NE 24) THEN Message, 'ZDEPTH must be set to 8 or 24.'
;                ENDIF
;                CASE state OF
;                    0: BEGIN
;                    
;                        ; Set depth if you have it, otherwise, just leave it alone.
;                        IF N_Elements(zdepth) NE 0 THEN BEGIN
;                            Device, Decomposed=state, Set_Pixel_Depth=zdepth
;                        ENDIF ELSE BEGIN
;                            Device, Decomposed=state
;                        ENDELSE
;                        END
;                    1: Device, Decomposed=state, Set_Pixel_Depth=24
;                 ENDCASE
;             ENDIF
             
             ; I am going to allow the user to set the decomposed state, but only if
             ; they are using a buffer depth that will allow that.
             IF (depth GE 24) THEN Device, Decomposed=state
             END
             
        'MAC': BEGIN
             IF (Float(!Version.Release) GE 5.2) THEN Device, Decomposed=state
             END
             
         'X': Device, Decomposed=state
         
         'WIN': Device, Decomposed=state
         
         'NULL':
         
         ELSE: Message, 'Unrecognized device. Assuming indexed color state.', /INFORMATIONAL
         
    ENDCASE
    
END    
;;==============================================================================
; $Id: //depot/Release/ENVI50_IDL82/idl/idldir/lib/strsplit.pro#1 $

; Copyright (c) 1999-2012, Exelis Visual Information Solutions, Inc. All
;       rights reserved. Unauthorized reproduction is prohibited.

;+
; NAME:
;       STRSPLIT
;
; PURPOSE:
;   Wrapper on the build in system routine STRTOK that implements exactly
;   the same interface as STRTOK, but with the STRSPLIT name.
;
;       The reason for doing this is so that if a user has their own
;   STRSPLIT in their local user library, their version will superceed
;   this one. RSI does not recommend this practice, but it is
;   allowed for backwards compatability reasons. See the
;       documentation for STRSPLIT in the IDL Reference manual
;   for details on arguments, keywords, and results.
;
;
; MODIFICATION HISTORY:
;   14 October 1999, AB, RSI.
;   AB, 5/4/2001, Switch from using _EXTRA to _STRICT_EXTRA, so that
;       incorrect keywords passed to STRTOK will issue proper
;       error messages instead of being silently ignored.
;   CT, 5/12/2010: Add support for string & pattern arrays, returns LIST.
;-

function strsplit, stringIn, pattern, $
  COUNT=count, LENGTH=length, EXTRACT=extract, $
  _REF_EXTRA=extra

    compile_opt idl2, hidden
    ON_ERROR, 2  ; return to caller
    
    ; Handle array input. Return a LIST of results.
    n = N_ELEMENTS(stringIn)
    if (n gt 1) then begin

      np = N_ELEMENTS(pattern)
      if (np gt 1 && np ne n) then $
        MESSAGE, 'PATTERN must be a scalar or have the same number of elements as STRING.'

      result = LIST()
      count = LONARR(n)
      if (ARG_PRESENT(length)) then length = LIST()

      for i=0,n-1 do begin
        result1 = (N_PARAMS() eq 1) ? $
          STRTOK(stringIn[i], COUNT=c, LENGTH=len, EXTRACT=extract, _STRICT_EXTRA=extra) : $
          STRTOK(stringIn[i], pattern[[i]], COUNT=c, LENGTH=len, EXTRACT=extract, _STRICT_EXTRA=extra)
        result.Add, result1
        count[i] = c
        if (ARG_PRESENT(length)) then length.Add, len
      endfor

      return, result
    endif


    ; Scalar input.
    RETURN, (n_params() eq 1) ? $
      STRTOK(stringIn, COUNT=count, LENGTH=length, EXTRACT=extract, _STRICT_EXTRA=extra) : $
      STRTOK(stringIn, pattern, COUNT=count, LENGTH=length, EXTRACT=extract, _STRICT_EXTRA=extra)

end
;;========================================================================
; docformat = 'rst'
;
; NAME:
;   cgDefCharSize
;
; PURPOSE:
;   Defines a default character size for Coyote Graphics routines (cgPlot, cgContour, etc.)
;   IF !P.Charsize is set, the function simply returns !P.Charsize.
;
;******************************************************************************************;
;                                                                                          ;
;  Copyright (c) 2010, by Fanning Software Consulting, Inc. All rights reserved.           ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;
;+
;   Defines a default character size for Coyote Graphics routines (cgPlot, cgContour, etc.)
;   IF !P.Charsize is set, the function simply returns !P.Charsize.
;
; :Categories:
;    Graphics, Utilities
;    
; :Keywords:
;     adjustsize: in, optional, type=boolean, default=0
;        If this keyword is set, the output character size is adjusted to
;        fit the size of the output graphics window. No adjustment is ever
;        done in PostScript. Applies only when !P.Charsize=0.
;     font: in, optional, type=integer, default=!P.Font
;        The font type: -1 = Hershey, 0 = hardware, 1 = true-type. 
;          
; :Author:
;       FANNING SOFTWARE CONSULTING::
;           David W. Fanning 
;           1645 Sheely Drive
;           Fort Collins, CO 80526 USA
;           Phone: 970-221-0438
;           E-mail: david@idlcoyote.com
;           Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Written, 11 January 2011. DWF.      
;        Added an ADJUSTSIZE keyword to allow adjustable sizing of characters
;           in resizeable graphics windows. 24 April 2011. DWF.  
;        Made sure this program only adjusts text size on devices that support 
;           windows. 20 July 2011. DWF.
;        Made PostScript default font sizes on Windows slightly larger to conform 
;           with function graphics output in IDL 8.2.3. 17 June 2013. DWF.
;
; :Copyright:
;     Copyright (c) 2011, Fanning Software Consulting, Inc.
;-
FUNCTION cgDefCharSize, ADJUSTSIZE=adjustsize, FONT=font

    Compile_Opt idl2
    
    ; Return to caller on an error.
    On_Error, 2
    
    ; Check parameters
    IF N_Elements(font) EQ 0 THEN font = !P.Font
    
    ; If the current window is a cgWindow, then the ADJUSTSIZE property of the
    ; window is used to set the AdjustSize keyword. This can only be done on
    ; devices that support windows.
    IF ~((!D.Flags AND 256) NE 0) THEN adjustsize = 0
    IF (N_Elements(adjustsize) EQ 0) THEN BEGIN
    
        ; Each instance of cgWindow will store evidence of its
        ; existance in a linked list.
        DefSysV, '!FSC_WINDOW_LIST', EXISTS=exists
        IF ~exists THEN BEGIN
            adjustsize = 0 
        ENDIF ELSE BEGIN
            IF Obj_Valid(!FSC_WINDOW_LIST) THEN BEGIN
                wid = cgQuery(/Current, COUNT=count)
                IF count GT 0 THEN BEGIN
                   IF wid EQ !D.Window THEN BEGIN
                       void = cgQuery(ObjectRef=windowObj, /Current)
                       IF Obj_Valid(windowObj) THEN windowObj -> GetProperty, AdjustSize=adjustsize
                   ENDIF ELSE adjustsize = 0
                ENDIF ELSE adjustsize = 0
            ENDIF ELSE adjustsize = 0
        ENDELSE
            
    ENDIF

    ; Calculate a default character size. We absolutely do not want to
    ; do this if !P.Charsize is not set to its default value of 0.
    IF !P.Charsize EQ 0 THEN BEGIN
        
            CASE StrUpCase(!Version.OS_Family) OF
            
                'WINDOWS': BEGIN
                    IF Total(!P.MULTI) EQ 0 THEN BEGIN
                        thisCharsize = 1.25                         
                    ENDIF ELSE BEGIN
                        totalplots = !P.Multi[1]*!P.Multi[2]*(!P.Multi[3] > 1)
                        CASE 1 OF
                            totalplots LE 4: thisCharsize = 1.25
                            totalplots GT 4: thisCharsize = 1.00
                        ENDCASE
                    ENDELSE
                    IF (font EQ 1) THEN BEGIN
                        IF Total(!P.MULTI) EQ 0 THEN BEGIN
                            thisCharsize = 1.75 
                        ENDIF ELSE BEGIN
                            totalplots = !P.Multi[1]*!P.Multi[2]*(!P.Multi[3] > 1)
                            CASE 1 OF
                                totalplots LE 4: thisCharsize = 1.75
                                totalplots GT 4: thisCharsize = 1.50
                            ENDCASE
                        ENDELSE
                    ENDIF
                    END
                    
                ELSE: BEGIN
                    IF Total(!P.MULTI) EQ 0 THEN BEGIN
                        thisCharsize = 1.50 
                    ENDIF ELSE BEGIN
                         totalplots = !P.Multi[1]*!P.Multi[2]*(!P.Multi[3] > 1)
                        CASE 1 OF
                            totalplots LE 4: thisCharsize = 1.50
                            totalplots GT 4: thisCharsize = 1.25
                        ENDCASE
                    ENDELSE 
                    IF (font EQ 1) THEN BEGIN
                        IF Total(!P.MULTI) EQ 0 THEN BEGIN
                            thisCharsize = 1.75 
                        ENDIF ELSE BEGIN
                            totalplots = !P.Multi[1]*!P.Multi[2]*(!P.Multi[3] > 1)
                            CASE 1 OF
                                totalplots LE 4: thisCharsize = 1.75
                                totalplots GT 4: thisCharsize = 1.50
                            ENDCASE
                        ENDELSE
                    ENDIF
                    END
            
            ENDCASE
             
        ; Adjust this size for the size of the window. Can't do this in PostScript
        ; for some reason, as it creates an extra page of output.
        IF !D.Name NE 'PS' THEN BEGIN
        
          ; The adjustment attempts to compensate for multiple plots in the window.
          IF Keyword_Set(adjustsize) THEN BEGIN
                thisCharsize = Str_Size('This is the text size for a normal window', $
                    0.65 * ((10-!P.multi[1])/10. > 0.5), INITSIZE=thisCharsize)
           ENDIF
        ENDIF

    ENDIF ELSE thisCharSize = !P.Charsize
        
    RETURN, thisCharSize
    
END

    