;------------------------------------------------------------------
; Script for plotting station data over Greece using land masking
; Author: Stavros Dafis - National Observatory of Athens 2016
; e-mail: sdafis@noa.gr
;------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "./shapefile_utils.ncl"

err = NhlGetErrorObjectId()
  setvalues err
    "errLevel" : "Fatal" ; only report Fatal errors
end setvalues 

begin
  
  ascii_filename = "eaa_antonis.txt"

    ncols = numAsciiCol(ascii_filename)
    nflds = numAsciiRow(ascii_filename)  ; number of available stations

 stations = asciiread(ascii_filename,(/nflds,ncols/),"float")  ;for fixed number of stations
  
  print("Number of stations = " + nflds)
  print("Number of variables = " + ncols)

  x = stations(:,2)  ; Column 3 of file contains LONGITUDE values.
  y = stations(:,1)  ; Column 2 of file contains LATITUDE values.
  z = stations(:,9)  ; Column 10 of file contains TEMPERATURE values.
;print(z)
  x = where(x .le. 33 .and. x .ge. 18,x, 33)
  y = where(y .le. 42.9 .and. y .ge. 34,y,35)
  z = where(z .le. 46 .and. z .ge. -17,z,10)  ; Check for wrong values


;print(dimsizes(z))
;print(dimsizes(x))

;print(min(x))

;exit


;---------------------------- DATES ----------------------------------------
 month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                    "Oct","Nov","Dec"/)

   year   = tointeger(stations(1,3))   ; Convert to integer for
   month  = tointeger(stations(1,4))    ; use sprinti 
   day    = tointeger(stations(1,5))
   hour   = tointeger(stations(1,6))
   minute = tointeger(stations(1,7))

date_str =  sprinti("%0.2i ", day) + month_abbr(month) + " " + sprinti("%0.4i", year)+ \
             " at " +sprinti("%0.2i", hour)+":"+sprinti("%0.2i ", minute) + "local time"
;---------------------------------------------------------------------

  numxout = 90  ; Define output grid for call to "natgrids".
  numyout = 90
  xmin    = 19.15
  ymin    = 34.60
  xmax    = 29.9
  ymax    = 41.9
  zmin    = min(z)
  zmax    = max(z)

print("Min Temp = " + zmin)
print("Max Temp = " + zmax)

  xc      = (xmax-xmin)/(numxout-1)
  yc      = (ymax-ymin)/(numyout-1)
  xo      = xmin + ispan(0,numxout-1,1)*xc
  yo      = ymin + ispan(0,numyout-1,1)*yc

  ;zo = natgrid(x, y, z, xo, yo)  ; Interpolate.
  ;zo = natgrids(x, y, z, xo, yo) ; Interpolate.
  ;zo = dspnt2(x, y, z, xo, yo)  ; Interpolate.

 ; zo = dsgrid2(x, y, z, xo, yo)  ; Interpolate.
   
   
   ;---DGRID2

zo= transpose(dsgrid2(x, y, z, xo, yo))

zo!0 = "yo"
zo!1 = "xo"
zo&yo=  yo ;(::-1)                      ; don't know why
zo&xo=  xo ;(::-1)
print("---")
printVarSummary(zo)
printMinMax(zo,0)

;wrf_smooth_2d(zo, 3)
   
   
  xo@long_name = "Lon"   ; Define some attributes.
  yo@long_name = "Lat"
  zo@long_name = "Temperature 2m in ~S~o~N~C"
 
 klat = dimsizes(yo)
klon = dimsizes(xo)

real_temp = generate_2d_array(-5,40,min(zo),max(zo),50,(/klat,klon/))

  wks_type = "png" 
  wks_type@wkWidth = 1200  
  wks_type@wkHeight = 1200
  xwks   = gsn_open_wks(wks_type,"noa_temp")
 
  
 gsn_define_colormap(xwks,"BlAqGrYeOrReVi200") ; For Temperature / Dew Point
 ;gsn_reverse_colormap(xwks) 

 ;gsn_define_colormap(xwks,"t2m_29lev")
 ;gsn_define_colormap(xwks,"BlAqGrYeOrRe ") 

;----------- Begin plot ----------------------------------------------------------

  pltres = True 

  resources                       = True
  resources@sfXArray              = xo            ; X axes data points
  resources@sfYArray              = yo            ; Y axes data points

;Map/boundary options
;---------------------
  resources@mpGeophysicalLineColor      = "Black"
  resources@mpNationalLineColor         = "Black"
  resources@mpGridLineColor             = "transparent"
  resources@mpLimbLineColor             = "transparent"
  resources@mpPerimOn                   =  True
  resources@mpPerimLineColor            = "Black"
  resources@mpGeophysicalLineThicknessF = 2.0
  resources@mpGridLineThicknessF        = 1.0
  resources@mpLimbLineThicknessF        = 1.0
  resources@mpNationalLineThicknessF    = 0.5
  resources@mpProjection                = "Stereographic"

  
 resources@mpRelativeCenterLon         = True          ; set a center lon
 
 ;resources@mpCenterLonF                = lon2d(19,22)  ; center lon
 ;resources@mpRelativeCenterLat         = True		     ; set a center lat
 ;resources@mpCenterLatF                = 40.           ; center lat 
 ;resources@mpLambertParallel1F   = 0.01
 ;resources@mpLambertParallel2F   = 89.9
 ;resources@mpLambertMeridianF    = 40.
 
  resources@mpDataSetName         = "Earth..4" 
  resources@mpOutlineBoundarySets = "AllBoundaries"
  resources@mpDataBaseVersion     = "MediumRes"     ; Important for masking (not HighRes)
  resources@mpDataResolution      = "Finest"
  
  resources@mpLimitMode                 = "LatLon"
  resources@mpMinLatF    = ymin ;34.7   ;ymin - 0.45     
  resources@mpMaxLatF    = ymax ;41.35  ;ymax + 0.25   
  resources@mpMinLonF    = xmin ;19.25  ;xmin - 0.68  
  resources@mpMaxLonF    = xmax ;29.8   ;xmax + 1.49 

;-- mask areas
mask_specs = (/"Greece_1"/)

;-- fill areas
  fill_specs = (/"Turkey","Bulgaria","Albania","Macedonia", "Serbia", "Montenegro"/)

  resources@mpFillOn 	             = True
  resources@mpLandFillColor        = "Transparent"
  resources@mpOceanFillColor       = "LightBlue"
  resources@mpInlandWaterFillColor = "LightBlue"
  resources@mpOutlineOn            = True
  resources@mpFillBoundarySets     = "Geophysical"
 ;resources@mpFillAreaSpecifiers   = "Geophysical"
  resources@mpAreaMaskingOn        =  True
  resources@mpMaskAreaSpecifiers   =  mask_specs
  resources@mpFillAreaSpecifiers   =  fill_specs
  resources@mpAreaTypes            = 3
  resources@mpSpecifiedFillColors  = (/"grey","grey","grey","grey","grey","grey"/)
  resources@mpGridMaskMode         = "MaskOcean"  ;-- don't draw grid over water
  resources@mpMaskAreaSpecifiers   = mask_specs
 ;resources@mpGeophysicalLineColor = "red"        ; color of continental outlines
  resources@mpGridLineDashPattern  = 2            ; lat/lon lines as dashed

  country_border     =  True   ;-- True: draw country borderlines
  use_black_bg       =  True   ;-- True: use grey18/black background and white foreground color

if(use_black_bg) then
     xwks@wkBackgroundColor   = "grey18"
     xwks@wkForegroundColor   = "white"
  end if

  resources@tfDoNDCOverlay          = True
  resources@tiMainString            = "Temperature data from automatic meteorological stations"  ; Main title
  resources@tiMainFont              = "helvetica"
  resources@tiMainFontHeightF       = 0.03
  resources@tiXAxisString           = ""    ; X axis label.
  resources@tiYAxisString           = ""    ; Y axis label.
  resources@pmTickMarkDisplayMode   = "Always"
  resources@tmXBMode                = "Explicit"
  resources@tmXBValues              = (/x/)
  resources@tmYBLabels              = (/y/)

  resources@tmXTOn                  = False ; Turn off tickmarks and labels 
  resources@tmYROn                  = False 
  resources@tmXBLabelFontHeightF    = 0.015
  resources@tmYLLabelFontHeightF    = 0.015
  resources@tmYLDataBottomF         ="Min Temp: "+ zmin +"~S~o~N~C/Max Temp:"+ zmax +"~S~o~N~C"
  
  resources@lbTitleString           = "Temperature 2m in ~S~o~N~C"
  resources@lbTitleOn               = False 
  resources@lbOrientation           = "Vertical"   ;"Horizontal" ; Draw it horizontally.
  resources@lbPerimOn               = False        ; Turn off perimeter.
  resources@lbJustification         = "TopLeft"
  resources@lbLabelAlignment        = "InteriorEdges"
  resources@lbLabelJust             = "CenterLeft"
  resources@lbLabelOffsetF          = 0.5
  resources@lbBoxMinorExtentF       = 1.0
  resources@lbLabelFontHeightF      = 0.015
  resources@lbBottomMarginF         = 0.0
  resources@lbLeftMarginF           = 0.0
  resources@lbRightMarginF          = 0.0
  resources@lbTopMarginF            = 0.0
  resources@pmLabelBarWidthF        = 0.03

  resources@lgAutoManage            = False
  resources@lgLabelFontHeightF      = 0.005
  resources@pmLabelBarDisplayMode   = "Always"     ; Turn on a label bar.
  resources@pmLabelBarSide          = "Bottom"     ; Change location of label bar.
  ;resources@vpYF                   = 0.9          ; Change Y location of plot.
  
  ;resources@mpShapeMode = "FreeAspect"
  ;resources@vpWidthF = 1.2
  ;resources@vpHeightF = 1.0

  resources@gsnMaximize        = True
  resources@gsnLeftString      = "Min Temp: "+ zmin +"~S~o~N~C / Max Temp: "+ zmax +"~S~o~N~C ~C~"+ \
                                 "data from "+ nflds +" stations"

 ;resources@gsnRightString     = "Date: "+ day +"/"+ mon +"/"+ year +"~C~" + \
                               ; "Time:"+ hour +":"+ minute +" (local time)"

   resources@gsnRightString    = ""+date_str+""

  resources@gsnRightStringFontHeightF =  0.015
  resources@gsnLeftStringFontHeightF  =  0.015

  resources@cnGridBoundFillPattern  = True
  resources@cnFillDrawOrder         = "Predraw"
  resources@cnFillOn                = True     ; Turn on contour fill.
  resources@cnInfoLabelOn           = False    ; Turn off info label.
  resources@cnLinesOn               = False
  resources@cnLineLabelsOn          = True    ; Turn off line labels.
  resources@cnFillOpacityF          = 0.9
  resources@cnLinesOn               = False
  resources@cnMonoLineColor         = True
  resources@cnMonoLevelFlag         = True
  resources@cnLineColor             = "Black"
  resources@cnLineLabelDensityF     = 1.
  resources@cnLineLabelInterval     = 1.
  resources@gsnContourLineThicknessesScale = 0.1
  resources@cnLineLabelPlacementMode       = "Computed"
  resources@cnExplicitLegendLabelsOn       = True
  resources@cnHighLabelsOn          = False
  resources@cnLowLabelsOn           = False
  resources@cnLineLabelPerimOn      = False   
  resources@cnLevelSelectionMode    = "ManualLevels"
  resources@cnLevelSpacingF         =  2.
  resources@cnMinLevelValF          = -8. ;zmin 
  resources@cnMaxLevelValF          = 42. ;zmax  
 ;resources@cnSmoothingOn           = True
 ;resources@cnSmoothingTensionF     = -2
  resources@ContourParameters       = (/ 0., 20., 0.5 /)
 ;resources@cnLevels                = (/ 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20 /)

  ;zo!0 = "i"  ; Name the dimensions of "zo".
  ;zo!1 = "j"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;INSERT RUN INFO
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
   ;txres                  = True         ; additional info
   ;resources@txFontHeightF = 0.001        ; size of meta data
   resources@txFont        = "helvetica"
   resources@txJust        = "CenterCenter"
   info= "(c)-www.metar.gr  ~C~" + \
               "National Observatory of Athens, Greece  ~C~"
 ;gsn_text_ndc(xwks, info, 0.02, 0.02, resources)    ; note when plotted

;*********************************************************************
 txres = True
 txres@txFontHeightF = 0.015       ; size of meta data
 txres@txFont        = "helvetica"
 txres@txJust        = "BottomLeft"
 
;txres@txString      = "(c)- data from noa.gr & metar.gr  ~C~" 

   info= "National Observatory of Athens, Greece - meteo.gr"
 
;-------------------------------------------------------------------------
;Plot
;-------------------------------------------------------------------------
 resources@gsnDraw      = False
 resources@gsnFrame     = False

contour = gsn_csm_contour_map(xwks,zo,resources) 

;plot = gsn_csm_contour_map_overlay(xwks,contour,id,resources,resources)
;contour = gsn_contour(xwks,zo(j|:,i|:),resources) ; Draw a contour plot.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;---------------------------- Europe -----------------------------------------
dir2 = "/home/dafis/SHP_Files/europe/eu_countries/"

  lnres2                   = True         
  lnres2@gsLineThicknessF  = 2.6          
  lnres2@gsLineColor       = "Black"  

shp_name2  = systemfunc("ls " + dir2 + "CNTR_RG_10M_2014.shp")

;---Attach the Europe outlines
  ;  id2 = gsn_add_shapefile_polylines(xwks,contour,shp_name2,lnres2)

dir = "/home/dafis/SHP_Files/greece/Greek_shp/"

;---Set some options for the polylines
  lnres                   = True        
  lnres@gsLineThicknessF  = 1.5          
  lnres@gsLineColor       = "Black" 

shp_name  = systemfunc("ls " + dir + "GRC_adm2.shp")

lnres@return_mask = True

;---Attach the Greek outlines
   id = gsn_add_shapefile_polylines(xwks,contour,shp_name,lnres)

;maximize_output(wks,True)

gsn_text_ndc(xwks,info, 0.08, 0.11, txres)    ; note when plotted

; Maximize plots in frame.
;psres = True                                                               
;maximize_output(xwks,psres)  ; calls draw and frame for you
;======================================================================
txres1 = True
txres1@txFontHeightF = 0.0070       ; size of meta data
txres1@txFont        = "helvetica"
;txres1@txJust        = "BottomLeft"
 
lat2 = new(dimsizes(real_temp),typeof(yo),-1)
copy_VarAtts(yo,lat2)

;lat2!0 = "south_north"
;lat2&south_north = yo
;lat2!1 = "west_east"
;lat2&west_east = xo

lon2 = new(dimsizes(real_temp),typeof(xo),-1)
copy_VarAtts(xo,lon2) 

nb = dimsizes(real_temp) 

;do nl=0,nb(0)-1,1          
;do ml=0,nb(1)-1,1

;lat2(nl,ml) = yo(nl) 
;lon2(nl,ml) = xo(ml)    

;end do 
;end do

;tl = niceLatLon2D(lat2,lon2)  ; check
;print(tl)

  opt2             = True
  opt2@debug       = True
  ;opt2@shape_var   = "NAME_1"
  ;opt2@shape_names = "VARNAME_2"

;z_new = zo(j|:,i|:)
;z_new!0 = "j"  ; Name the dimensions of "zo".
;z_new!1 = "i"
;z_new&i = xo
;z_new&j = yo

;maskv = shapefile_mask_data(z_new,shp_name,opt2)
;maskv!0 = "i" 
;maskv!1 = "j"

;printVarSummary(maskv)
;print(maskv)
;exit

;z11 = maskv ;(j|:,i|:)
;z22= decimalPlaces(z11,0,True)
;z22@missing_value = -20
;z22@_FillValue = -20

;print(z22)
;exit

;temps1 = sprintf("%3.0f", z22(::3,::3))
;temps= where(temps1.eq.-20, " ", temps1)

;printVarSummary(temps)
;printVarSummary(zo)
;printVarSummary(lon2)
;printVarSummary(lat2)
;print(temps)
;exit


;ath_t = sprintf("%3.0f ", temps)  
;gsn_text_ndc(xwks,ath_t+"~S~o~N~C", 0.55, 0.50, txres1) 

;text2 = gsn_add_text(xwks,contour,temps,lon2(::3,::3),lat2(::3,::3),txres1)

text = gsn_add_text(xwks,contour,sprintf("%3.0f",z),x,y,txres1)

;
; Retrieve the ids of the text strings that were just added so we can
; retrieve their heights, widths, and locations. This information 
; will be used to determine if any text strings are overlaying other text
; strings.

;nmsg = num(ismissing(z))

;print(z)


  nstrs=dimsizes(z)
  getvalues contour
    "pmAnnoViews" : text_ids
  end getvalues

  xndc    = new(nstrs,float)     ;  X center of box in NDC coords
  yndc    = new(nstrs,float)     ;  Y center of box in NDC coords
  rgt     = new(nstrs,float)     ; Hold right position of text box.
  lft     = new(nstrs,float)     ;  "   left     "     "   "    " 
  top     = new(nstrs,float)     ;  "   top      "     "   "    " 
  bot     = new(nstrs,float)     ;  "   bottom   "     "   "    " 
  width   = new(nstrs,float)     ;  width of each text box
  height  = new(nstrs,float)     ;  height of each text box

; Loop through and get the width and height of each box.
;
  do i=0,nstrs-1
    getvalues text_ids(i)
      "vpWidthF"  : width(i)
      "vpHeightF" : height(i)
    end getvalues
  end do

; Convert the lat/lon center of each box to NDC coordinates, since this
; is what the width and height values are in.
;
  datatondc(contour,x,y,xndc,yndc)

; Calculate the four corners of each text box in NDC coordinates.
;
  top = yndc + height/2.
  bot = yndc - height/2.
  lft = xndc -  width/2.
  rgt = xndc +  width/2.

;
; Now we're going to loop through each text string that has been added,
; and compare it to all other text strings that have been added to see if
; they overlap. If an overlap is encountered, one of the strings will be
; tagged for removal, and a box will be drawn around it (just to show,
; for debug purposes, which text strings are getting removed).
;
  removed_list = new(nstrs,integer)
  rcount = 0                         ; Number of removed text boxes.

;
; ibox1 is index of box that we're checking to see if it overlaps.
; ibox2 is the index of the box that we're checking box ibox1
; against.
;
  do ibox1 = 0,nstrs-1
    ibox2 = 0
    overlap_found = False 

    do while(ibox2.lt.nstrs.and..not.overlap_found)
;
; Check if any one of the corners represented by box "ibox1" is inside
; the box representd by "ibox2".  If so, remove it.  Make sure you are
; not checking a box against itself, or against a box that has already
; been removed.
;
      if (ibox1.ne.ibox2.and.(.not.any(ibox2.eq.removed_list))) then
        if ( (top(ibox1).le.top(ibox2).and.top(ibox1).ge.bot(ibox2).and. \
              lft(ibox1).le.rgt(ibox2).and.lft(ibox1).ge.lft(ibox2)).or. \
             (bot(ibox1).le.top(ibox2).and.bot(ibox1).ge.bot(ibox2).and. \
              lft(ibox1).le.rgt(ibox2).and.lft(ibox1).ge.lft(ibox2)).or. \
             (bot(ibox1).le.top(ibox2).and.bot(ibox1).ge.bot(ibox2).and. \
              rgt(ibox1).le.rgt(ibox2).and.rgt(ibox1).ge.lft(ibox2)).or. \
             (top(ibox1).le.top(ibox2).and.top(ibox1).ge.bot(ibox2).and. \
              rgt(ibox1).le.rgt(ibox2).and.rgt(ibox1).ge.lft(ibox2))) then

          overlap_found = True         ; Mark that an overlap has been found.

          removed_list(rcount) = ibox1    ; Update removed_list and
          rcount = rcount + 1             ; its counter.
        end if
      end if
      ibox2 = ibox2 + 1               ; Update the box counter
    end do
  end do

;
; Here's the code that removes the overlapping text strings.
;
  NhlRemoveAnnotation(contour,text_ids(removed_list(:rcount-1)))

;
; Now redraw the plot, which should have all of the overlapping text 
; strings removed.
;
  setvalues contour

  ;resources@tiMainFont           = 33 ;greek
  ;resources@tiMainString         = "Qermokras"+iacute+"a ap"+oacute+" aut"+oacute+"matouV ~C~ metewrologiko"+uacute+"V staqmo"+uacute+"V ~S~o~N~~F21~C ~C~ "  ; Main title  
  
 ;"tiMainString" :  "Overlapping text strings removed"
  end setvalues




draw(contour)
frame(xwks)

;-----------------------------------------------------------------------------------
logo = "./meteo.jpg"

       outfile = "noa_temp"
     
    system("convert -trim +repage " + outfile + ".png " + outfile + ".trimmed.png")
    system("mv " + outfile + ".trimmed.png " + outfile + ".png")

delete(xwks)

cmd = "composite -geometry 165x165+820+760 /home/dafis/LOGOS/meteo.gif noa_temp.png noa_temp.png"

system(cmd)
;************************************************************************************************
;------------------------------------------------------------------------------------------------
end