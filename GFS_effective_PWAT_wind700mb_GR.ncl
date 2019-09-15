;========================================================================================
; Plot effective PWAT and wind speed at 700 hPa from GFS grib forecast files
; Author: Stavros Dafis
; sdafis@noa.gr - 2018
;========================================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;------------------------------------------------------------------------
err = NhlGetErrorObjectId() 
  setvalues err 
    "errLevel" : "Fatal" ; only report Fatal errors 
  end setvalues 
;------------------------------------------------------------------------

begin

DataDir = "/dafis/GFS_operational/gribs"

files = systemfunc ("ls  /dafis/MODELS/GFS/OPERATIONAL/gribs_GRE_EM/gfs.20*.grb") ; read all the available gribs

b  = addfiles (files, "r")   ; note the "s" of addfile

ListSetType (b, "join")

;----------------------------------------------------------------------
numfiles = dimsizes(files)
print("Going to work on  " + numfiles + " GFS gribs")

do it = 1, numfiles-1,1
print("Working on file: " + files(it))

;========================================================================================================
; --------------------- surface variables --------------------------
lat0  = b[:]->lat_0(0,:)
lon0  = b[:]->lon_0(0,:)
slp0  = b[it]->PRMSL_P0_L101_GLL0 
;-----------------------------------------------------------------------------------
tk0 = b[it]->TMP_P0_L100_GLL0 ; Temp in Kelvin
u0  =  b[it]->UGRD_P0_L100_GLL0   ;( lv_ISBL0, lat_0, lon_0 )
v0  =  b[it]->VGRD_P0_L100_GLL0   ;( lv_ISBL0, lat_0, lon_0 )
pw  =  b[it]->PWAT_P0_L200_GLL0   ;( lat_0, lon_0 ) GFS PWAT
rh0 = b[it]->RH_P0_L100_GLL0       ; RH_P0_L100_GLL0 ( lv_ISBL0, lat_0, lon_0 ) relative hum gfs

iso_lvl0 = b[it]->lv_ISBL0 ; vertical levels mb

;============================================================================================================
strs = str_split(files(it), ".") ; take the strings for creating initial time/date strings on plots
;print("'" + strs(1) + "'")  ; eg 2016072800
;print("'" + strs(4) + "'")  ; f0**

;Initial time:
dtime_map = (/4,2,2,2/)
initial_time = str_split_by_length(strs(1), dtime_map) ; split the name and get the dates
;print("" + initial_time + "")

;Valid time:
vtime_map = (/2,2/)
valid_time1 = str_split_by_length(strs(4), vtime_map) ; split the name and get the dates
;valid_time1 = tk@forecast_time
print("=========================================================")
print("Valid time: t+" +  tk0@forecast_time + "z")
print("------------------")

valid_time_t = stringtoint(tk0@forecast_time)

;=======================================================================================
; Open a working station
;-----------------------------------------------------------------------------
;type = "x11"
 type = "png"

type@wkWidth = 1800  
type@wkHeight = 1800 

outfile = "GFS_PW_GR_t" +valid_time_t

wks = gsn_open_wks(type,outfile)

 gsn_define_colormap(wks,"GMT_haxby")
 gsn_reverse_colormap(wks)  

;========================================================================
initial_time00 = str_split(tk0@initial_time, "/ (:)")
;print("" + initial_time + "")

;-------------
; TIME
;-------------
;Initial time:

 month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                    "Oct","Nov","Dec"/)

  year_init = tointeger(initial_time(0))
 month_init = tointeger(initial_time(1))
   day_init = tointeger(initial_time(2))
  hour_init = tointeger(initial_time(3))

dow2 = day_of_week(year_init,month_init,day_init) 

if (dow2.eq.0) then
namday2 = "Sun"
end if
if (dow2.eq.1) then
namday2 = "Mon"
end if 
if (dow2.eq.2) then
namday2 = "Tue"
end if
if (dow2.eq.3) then
namday2 = "Wed"
end if
if (dow2.eq.4) then
namday2 = "Thu"
end if
if (dow2.eq.5) then
namday2 = "Fri"
end if
if (dow2.eq.6) then
namday2 = "Saturday"
end if

date_initial = sprinti("%0.2i ", day_init) + month_abbr(month_init) + " " + sprinti("%0.4i", year_init) + " " +sprinti("%0.2iz ", hour_init)
 
   print("")
   print("Init. time : "+namday2+" "+date_initial+"") 

;===============================================

mm1 = stringtoint(str_get_cols(tk0@initial_time,0,1)) 
dd1 = stringtoint(str_get_cols(tk0@initial_time,3,4)) 
yyyy1 = stringtoint(str_get_cols(tk0@initial_time,6,9)) 
hh1 = stringtoint(str_get_cols(tk0@initial_time,12,13)) 

initial_time11 = ut_inv_calendar(yyyy1, mm1,dd1, hh1, 00, 00, "hours since 1800-01-01 00:00:00",0) 
;print(initial_time11) ; hours since 1800

  utc_date = ut_calendar(initial_time11, 0)
year_val   = tointeger(utc_date(:,0))    ; Convert to integer for
month_val  = tointeger(utc_date(:,1))    ; use sprinti 
day_val    = tointeger(utc_date(:,2))
hour_val   = tointeger(utc_date(:,3))
minute_val = tointeger(utc_date(:,4))

date = ut_inv_calendar(year_val,month_val,day_val,hour_val,minute_val,00, "hours since 1800-01-01 00:00:00", 0) ; metatropi se wres gia na ginei h prosthesi
valid_time11 = date + stringtoint(tk0@forecast_time) ; add hours with time_step
valid_time11@units = initial_time11@units 
;print(valid_time11)

utc_datef = ut_calendar(valid_time11, 0)   ; hours in strings
year_valf   = tointeger(utc_datef(:,0))    ; Convert to integer for
month_valf  = tointeger(utc_datef(:,1))    ; sprinti 
day_valf    = tointeger(utc_datef(:,2))
hour_valf   = tointeger(utc_datef(:,3))
minute_valf = tointeger(utc_datef(:,4))


dow3 = day_of_week(year_valf,month_valf,day_valf) 

if (dow3.eq.0) then
namday3 = "Sun"
end if
if (dow3.eq.1) then
namday3 = "Mon"
end if 
if (dow3.eq.2) then
namday3 = "Tue"
end if
if (dow3.eq.3) then
namday3 = "Wed"
end if
if (dow3.eq.4) then
namday3 = "Thu"
end if
if (dow3.eq.5) then
namday3 = "Fri"
end if
if (dow3.eq.6) then
namday3 = "Sat"
end if

date_valid = sprinti("%0.2i ", day_valf) + month_abbr(month_valf) + " " + sprinti("%0.4i", year_valf) + " " +sprinti("%0.2iz ", hour_valf)
print("Valid time : "+namday3+" "+date_valid+"")
print("")

;====================================================================

  dlatmin = min(lat0)
  dlatmax = max(lat0)
  dlonmin = min(lon0)
  dlonmax = max(lon0)

;====================================================================
; PWAT

pw!0 = "lat"
pw!1 = "lon"
pw&lat = lat0
pw&lon = lon0


minpw = min(pw)
maxpw= max(pw)
print("Max PW: "+maxpw +" mm")


; effective PWAT (PW - Saturation Deficit)

; Saturation Deficit (difference between mixing ratio and saturated mixing ratio 
; at that level)


 tk = tk0(17:33,:,:) 
 rh = rh0(14:30,:,:) 
iso_lvl = iso_lvl0(17:33)

;VPD=(1-RH/100)×0.611 exp⁡(17.27T/(T+237.3))
q_uns    = mixhum_ptrh(conform(tk,iso_lvl*0.01,0), tk ,rh,1)*10  ; mix ratio (kg/kg)
;printVarSummary(q_uns)

rh    = where(rh .eq. 0, 1,rh )
q_sat    = (q_uns/rh)*100

SD = (q_sat - q_uns) ; kg/kg
SD_limit = SD(0:16,:,:)  ; 1000-300 hpa
;printMinMax(SD_limit,0)

SD_limit!0 = "lev"
SD_limit!1 = "lat"
SD_limit!2 = "lon"
SD_limit&lat = lat0
SD_limit&lon = lon0

sat_def = vibeta(iso_lvl(16:0)*0.01,SD_limit(lat|:,lon|:,lev|:),1,slp0*0.01,1000,300) ;kg/kg

pw2 = pw - sat_def ; [mm] - [kg/kg]
;printVarSummary(pw2)
;print(""+pw2+"")

;print(""+SD(14:30,15,15)+"   "+q_uns(14:30,15,15)+"     "+q_sat(14:30,15,15)+"  "+tk(14:30,15,15)+" "+rh0(14:30,15,15)+" "+iso_lvl1(14:30)*0.01+"")

pw2@long_name = "total column precipitable water" 
pw2@units = "kg/m2" 
;copy_VarCoords(ps,pw) 

print("Max effective PWAT: "+max(pw2) +" mm")
print("")


;------------------------------------------------------------------------
; WIND 700 hPa

u = u0(25,:,:)
v = v0(25,:,:)

u!0 = "lat"
u!1 = "lon"
u&lat = lat0
u&lon = lon0

v!0 = "lat"
v!1 = "lon"
v&lat = lat0
v&lon = lon0

wind = sqrt(u^2 + v^2)
wind!0 = "lat"
wind!1 = "lon"
wind&lat = lat0
wind&lon = lon0

minwind = min(wind)
maxwind = max(wind)

print("Min wind speed 700 hPa: "+minwind+" m/s")
print("Max wind speed 700 hPa: "+maxwind+" m/s")
print("")


u_1000 = b[it]->UGRD_P0_L100_GLL0 ;(30,:,:)   ;( lv_ISBL0, lat_0, lon_0 )
u_10_sur = u_1000(30,:,:)
v_1000 = b[it]->VGRD_P0_L100_GLL0 ;(30,:,:)   ;( lv_ISBL0, lat_0, lon_0 )
v_10_sur = v_1000(30,:,:)

u = u_10_sur(:,:)
v = v_10_sur(:,:)

u@units = "m/s"
v@units = "m/s"

;Calculate Convergence at 1000 hPa
	conv = u
	conv = uv2dv_cfd(u,v,lat0,lon0,2)
	conv_scaled = - conv/1e-4

conv_scaled!0 = "lat"
conv_scaled!1 = "lon"
conv_scaled&lat = lat0
conv_scaled&lon = lon0

;======================================================================================
; Vertically integrated moisture flux

;;--vertical integrated moisture flux
um = u0(30:14,:,:)
vm = v0(30:14,:,:)

mix_v    = mixhum_ptrh(conform(tk,iso_lvl*0.01,0), tk ,rh,1) *10 ; mix ratio (kg/kg)
;printVarSummary(mix_v)
;printVarSummary(um)

uq          =   um * mix_v(16:0,:,:)                       ;vapor flux
vq          =   vm * mix_v(16:0,:,:)

uq!0 = "lev"
uq!1 = "lat"
uq!2 = "lon"
uq&lat = lat0
uq&lon = lon0
vq!0 = "lev"
vq!1 = "lat"
vq!2 = "lon"
vq&lat = lat0
vq&lon = lon0

uqflux      =   vibeta(iso_lvl(16:0)*0.01,uq(lat|:,lon|:,lev|:),2,slp0*0.01,1010,700)
vqflux      =   vibeta(iso_lvl(16:0)*0.01,vq(lat|:,lon|:,lev|:),2,slp0*0.01,1010,700)
                                                           
qflux       =   uqflux + vqflux  ;sqrt(uqflux^2+vqflux^2)  ;kg / m/s
;printVarSummary(qflux)
printMinMax(qflux,0)

;-------------------------------------------------------------------------------------
 res                    = True
        
 ;res@mpProjection                = "CylindricalEqualArea"
  res@mpProjection                = "Mercator"
 ;res@mpProjection                = "LambertConformal"
 ;res@mpProjection                = "Orthographic" 
  res@mpCenterLonF                = 20.             ; choose center lon
  res@mpCenterLatF                = 40.             ; choose center lat
   
  res@mpGeophysicalLineColor      = "Black"
  res@mpNationalLineColor         = "Black"
  res@mpUSStateLineColor          = "Black"
  res@mpGridAndLimbOn             = False
  res@mpGridLineColor             = "Black"
  res@mpLimbLineColor             = "Black"
  res@mpPerimLineColor            = "Black"
  res@mpGridLatSpacingF           = 5.0
  res@mpGridLonSpacingF           = 5.0
  res@mpGeophysicalLineThicknessF = 5.0
  res@mpGridLineThicknessF        = 1.0
  res@mpLimbLineThicknessF        = 1.0
  res@mpNationalLineThicknessF    = 0.5
  res@mpUSStateLineThicknessF     = 0.5
  res@mpDataBaseVersion           = "HighRes"
  res@mpDataResolution            = "Medium"
  res@trGridType                  ="TriangularMesh"
  res@gsnAddCyclic                = False

   res@mpLimitMode                 = "LatLon"
   ;res@mpRelativeCenterLon         = True 
   res@mpMinLatF    =  34.
   res@mpMaxLatF    =  42.
   res@mpMinLonF    =  18. 
   res@mpMaxLonF    =  30.

 ;res@gsnRightString    = "Valid time: "+date_str(it)+""
  
res@gsnRightString  = "  Init. time: "+namday2+" "+date_initial+" ~C~Valid time: "+namday3+" "+date_valid+""
res@gsnLeftString   = "Effective Precipitable water (mm) (shaded) ~C~Winds at 700 hPa (m/s) (arrows)~C~GFS 0.25~S~o~N~ t+"+valid_time_t+"z"
;res@gsnLeftString   = "Effective Precipitable water (mm) (shaded) ~C~Integrated Water Vapor Transport (kg m~S~-1~N~ s~S~-1~N~) (arrows)~C~GFS 0.25~S~o~N~ t+"+valid_time_t+"z"

  res@gsnRightStringFontHeightF =  0.014
  res@gsnLeftStringFontHeightF  =  0.014
 
  res@sfXArray     = lon0    
  res@sfYArray     = lat0
  res@txFont       = "helvetica"
  res@tiMainFont              = "helvetica-bold"
  res@tiMainFontHeightF       = 0.020
  res@lgAutoManage            = False
  res@lgLabelFontHeightF      = 0.008
  res@pmTickMarkDisplayMode   = "Always"
  res@pmLegendHeightF         = 0.01
  res@tmXBLabelFontHeightF    = 0.014
  res@tmYLLabelFontHeightF    = 0.014 
  res@tmXTOn                  = False 
  res@tmYROn                  = False 
  res@pmLabelBarDisplayMode   = "Always"     
  res@lgLabelFontHeightF      = 0.005

  res@lbTitleOn               = False
 ;res@lbTitleString           = "gpdm"
  res@lbTitleFontHeightF      = 0.018
  res@lbLabelFontHeightF      = 0.015
  res@lbLabelAlignment        = "InteriorEdges"
  res@lbLabelJust             = "TopLeft"
  res@lbLabelOffsetF          = 0.05    
  res@pmLabelBarDisplayMode   = "Always"     
 ;res@pmLabelBarSide          = "Right"  
  
 ;res@lbAutoManage            = False
  res@lgAutoManage            = False
  res@lbJustification         = "TopLeft"
  res@lbPerimOn               = False   
  res@lbOrientation           = "vertical" 
  res@lbBoxMinorExtentF       = 0.3
  res@lbBottomMarginF         = 0.0
  res@lbLeftMarginF           = 0.0
  res@lbRightMarginF          = 0.0
  res@lbTopMarginF            = 0.0
  res@pmLabelBarOrthogonalPosF= 0.003
 ;res@pmLabelBarParallelPosF  = 0.5   ;up /down the whole bar
  res@pmLabelBarWidthF        = 0.1
  res@pmLabelBarHeightF       = 0.65
 
  res@cnLinesOn               = True
  res@cnLineLabelsOn          = True
  res@cnLineThicknessF        = 1.2
  res@cnLineLabelFontHeightF  = 0.009
  res@cnLineLabelFontThicknessF = 0.85
   res@cnLineLabelFont = "helvetica-bold"
  res@cnLineLabelFont         = "helvetica"
  res@cnLineColor             = "grey"
  res@cnLineLabelFontColor    = "Black"
  res@cnLevelSelectionMode = "ManualLevels"
  ;res@cnLevelSpacingF = 20.
  ;res@cnMinLevelValF  = 0.
  ;res@cnMaxLevelValF  = 100.
  res@cnLineLabelDensityF  = 1.5
  res@cnLineLabelInterval = 2.

;res@vpWidthF = 0.8
;res@vpHeightF = 0.7

;--------------------------------------
res@sfXArray     = lon0    
  res@sfYArray     = lat0
  res@txFont       = "helvetica"
  res@tiMainFont              = "helvetica-bold"
  res@tiMainFontHeightF       = 0.020
  res@lgAutoManage            = False
  res@lgLabelFontHeightF      = 0.008
  res@pmTickMarkDisplayMode   = "Always"
  res@pmLegendHeightF         = 0.01
  res@tmXBLabelFontHeightF    = 0.014
  res@tmYLLabelFontHeightF    = 0.014 
  res@tmXTOn                  = False 
  res@tmYROn                  = False 
  res@pmLabelBarDisplayMode   = "Always"     
  res@lgLabelFontHeightF      = 0.005

  res@lbTitleOn               = False
 ;res@lbTitleString           = "gpdm"
  res@lbTitleFontHeightF      = 0.018
  res@lbLabelFontHeightF      = 0.015
  res@lbLabelAlignment        = "InteriorEdges"
  res@lbLabelJust             = "TopLeft"
  res@lbLabelOffsetF          = 0.05    
  res@pmLabelBarDisplayMode   = "Always"     
 ;res@pmLabelBarSide          = "Right"  
  
 ;res@lbAutoManage            = False
  res@lgAutoManage            = False
  res@lbJustification         = "TopLeft"
  res@lbPerimOn               = False   
  res@lbOrientation           = "vertical" 
  res@lbBoxMinorExtentF       = 0.3
  res@lbBottomMarginF         = 0.0
  res@lbLeftMarginF           = 0.0
  res@lbRightMarginF          = 0.0
  res@lbTopMarginF            = 0.0
  res@pmLabelBarOrthogonalPosF= 0.003
 ;res@pmLabelBarParallelPosF  = 0.5   ;up /down the whole bar
  res@pmLabelBarWidthF        = 0.1
  res@pmLabelBarHeightF       = 0.65
 
  res@cnLinesOn               = False
  res@cnLineLabelsOn          = True
  res@cnLineThicknessF        = 0.5
  res@cnLineLabelFontHeightF  = 0.009
  res@cnLineLabelFontThicknessF = 0.9
  res@cnLineLabelFont         = "helvetica-bold"
  res@cnLineColor             = "grey"
  res@cnLineLabelFontColor    = "Black"
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLineLabelDensityF  = 1.5
  res@cnLineLabelInterval = 1.


  res@cnFillOn                = True
  res@cnInfoLabelOn           = False    
;res@cnLevelSelectionMode    = "ManualLevels" 
 res@cnLevelSelectionMode    = "ExplicitLevels"
  res@cnMinLevelValF          =    5.  
  res@cnMaxLevelValF          =   35.  
  res@cnLevelSpacingF         =    5.
  ;res@cnSmoothingOn           = True
  res@cnSmoothingTensionF     = -0.5
  res@cnSmoothingDistanceF    = 0.01

;res@vpWidthF = 0.8
;res@vpHeightF = 0.7

gsn_define_colormap(wks,"GMT_globe")
res@cnFillColors = (/ -1,129,125,113,97,65,49,33,28,2/)
res@cnLevels  = (/ 1,5,10,15,20,25,30,35,40 /)


;===========================


res_conv = True
 res_conv@cnFillOn                = False
  
  res_conv@cnLinesOn               = True
  res_conv@cnLineLabelsOn          = False
  res_conv@cnLineThicknessF        = 8.5
  res_conv@cnLineLabelFontHeightF  = 0.0098
  res_conv@cnLineLabelFontThicknessF = 0.095
  res_conv@cnLineLabelFont         = "helvetica-bold"
  res_conv@cnLineColor             = "Red"
  res_conv@cnLineLabelFontColor    = "Black"
  res_conv@cnLevelSelectionMode = "ManualLevels"
  res_conv@cnLineLabelDensityF  = 2.
  res_conv@cnLineLabelInterval =  1.
  res_conv@cnInfoLabelOn           = False    
  res_conv@cnLevelSelectionMode    = "ManualLevels" 
 ;res_conv@cnLevelSelectionMode    = "ExplicitLevels"
  res_conv@cnMinLevelValF          =  1.5  
  res_conv@cnMaxLevelValF          =  3.5  
  res_conv@cnLevelSpacingF         =  0.5
  ;res_conv@cnSmoothingOn           = True
  ;res_conv@cnSmoothingTensionF     = -0.5
  ;res_conv@cnSmoothingDistanceF    = 0.01

;----------------------------------------------------------------------
; Wind vectors

        res_vc                            = True
        res_vc@gsnDraw   = False
        res_vc@gsnFrame  = False
           
     ;res_vc@vcLevelSelectionMode     = "ExplicitLevels"
     ;res_vc@vcLevels     =  (/5,10,15,20,25,30,35,40,45,50/) 
     ;res_vc@vcLevels     =  (/5,10,15,20,25,30,35,40,45,50,55,60/) 

       res_vc@vcLevelSelectionMode     = "ManualLevels"
       res_vc@vcMinLevelValF           = 2.0
       res_vc@vcMaxLevelValF           = 50.0
       res_vc@vcLevelSpacingF          = 5.    

       ;res_vc@vcGlyphStyle               = "CurlyVector"
        res_vc@vcGlyphStyle               = "FillArrow"
        res_vc@vcMonoFillArrowFillColor   = True
       ;res_vc@vcMonoLineArrowColor       = True
        res_vc@vcFillArrowFillColor       = "grey30"
        res_vc@vcFillArrowsOn             = True
        res_vc@vcLineArrowThicknessF      = 1.0
        res_vc@vcMinDistanceF             = 0.022
        res_vc@vcRefLengthF               = 0.015
       ;res_vc@NumVectors                 = 35 
       ;res_vc@vcMinFracLengthF           = 0.01
        res_vc@vcMinMagnitudeF            = 2.5 ;100;

        cmap = read_colormap_file("BlAqGrYeOrReVi200")

        res_vc@vcLevelPalette      = cmap(6:193,:)

;;wind barb resources don't apply
;;      res_vc@vcGlyphStyle               = "WindBarb"
;;      res_vc@vcWindBarbLineThicknessF   = 1.5
;;      res_vc@vcWindBarbColor            = "Gray40"
        res_vc@vcFillArrowsOn             = True
     
        res_vc@vcRefAnnoOn    = False
        res_vc@vcRefMagnitudeF           =  10. ;100 ;
        res_vc@vcRefAnnoJust             = "BottomRight"
        res_vc@vcRefAnnoSide             = "Bottom"
        res_vc@vcRefAnnoString1          = "Reference speed"
        res_vc@vcRefAnnoString1          = "10"
        res_vc@vcRefAnnoSide             = "Bottom"
        res_vc@vcRefAnnoString2On        = False
        res_vc@vcRefAnnoPerimOn          = True
        res_vc@vcRefAnnoOrthogonalPosF   = -0.12
        res_vc@vcRefAnnoParallelPosF     = 0.999
       ;res_vc@vcRefAnnoBackgroundColor  = "Purple"
        res_vc@vcVectorDrawOrder         = "PostDraw"
      
        ;res_vc@vcRefAnnoOn    = True
        ;res_vc@vcRefAnnoString1          = "  IVT"
        ;res_vc@vcRefAnnoString2On        = True
        ;res_vc@vcRefAnnoPerimOn          = True
        ;res_vc@vcRefAnnoString2          = "     100~C~kg m~S~-1~N~ s~S~-1~N~"
        ;res_vc@vcRefAnnoOrthogonalPosF   = -0.28
        ;res_vc@vcRefAnnoParallelPosF     = 0.999

  res_vc@gsnLeftString     = ""
  res_vc@gsnRightString     = ""
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;INSERT INFO
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   txres               = True         
   txres@txFontHeightF = 0.014        
   txres@txJust        = "CenterRight"
   txres@txFont        = "helvetica"

   ;info= "(c) LMD, Ecole Polytechnique"

    info= "National Observatory of Athens, Greece - meteo.gr"

 gsn_text_ndc  (wks, info, 0.93, 0.08, txres)  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
   txtres4 = True
   txtres4@txFontHeightF = 0.0115
   txtres4@txJust        = "CenterLeft"
   txtres4@txFont        = "helvetica"
   txtres4@txFontColor   = "Red"

 gsn_text_ndc(wks,"Red contours: strong low-level convergence",0.16,0.08,txtres4)
 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
res@gsnMaximize  = False
res@gsnFrame     = False
res_conv@gsnMaximize  = False
res_conv@gsnFrame     = False
res_conv@gsnDraw = False

plot = gsn_csm_contour_map(wks,pw2,res)
converge = gsn_csm_contour(wks,conv_scaled,res_conv)
vector  = gsn_csm_vector(wks,u,v,res_vc)

;vector  = gsn_csm_vector(wks,uqflux,vqflux,res_vc)

overlay(plot,converge)
overlay(plot,vector)
;plot = gsn_csm_contour_map(wks,geopo,res)
;plot = gsn_csm_contour_map(wks,mixing,res)

;plot = gsn_csm_contour_map_overlay(wks,PW,wind,res,opts)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;---------------------------- Europe -----------------------------------------
dir2 = "/dafis/SHP_Files/europe/eu_countries/"

  lnres2                   = True         
  lnres2@gsLineThicknessF  = 3.0          
  lnres2@gsLineColor       = "Black"  

shp_name2  = systemfunc("ls " + dir2 + "CNTR_RG_10M_2014.shp")

;---Attach the Europe outlines
    id2 = gsn_add_shapefile_polylines(wks,plot,shp_name2,lnres2)


dir = "/dafis/SHP_Files/greece/Greek_shp/"

;---Set some options for the polylines
  lnres                   = True        
  lnres@gsLineThicknessF  = 1.5          
  lnres@gsLineColor       = "Black" 

shp_name  = systemfunc("ls " + dir + "GRC_adm2.shp")

lnres@return_mask = True

;---Attach the Greek outlines
   id = gsn_add_shapefile_polylines(wks,plot,shp_name,lnres)

 draw(plot)       

;************************************************************************************************
 
frame(wks)
     
    system("convert -trim +repage " + outfile + ".png " + outfile + ".trimmed.png")
    system("mv " + outfile + ".trimmed.png " + outfile + ".png")

delete(wks)

cmd = "composite -geometry 220x220+20+1320 /home/dafis/LOGOS/meteo.jpg " + outfile + ".png " + outfile + ".png "

system(cmd)


end do
   
end
