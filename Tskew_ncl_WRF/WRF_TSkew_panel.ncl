; T-Skew and hodograph plots for WRF-ARW outputs
;
; The script will have to be modified for post-processed WRF netcdf files. 
; Author: Stavros Dafis, NOA-GR - modified sections by previous versions provided
; by the NCL community. 
;
; It is recommended to use ncl version 6.2.0 or higher
;
; For any question:
; sdafis@noa.gr
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/skewt_func.ncl"
;
; Necessary libraries:
;
load "./Panel_skewt_func_dafis.ncl"
load "./Panel_hodograph_dafis.ncl"
load "./hodo_cartesian.ncl"

begin

print("				")

a = addfile("wrfout_d0X_XXXX-XX-XX_XX:XX:XX.nc","r") 	

;===============================================================================

locs = (/ "Location1", "Location2"   /)

lats = (/ 44.5, 45.0 /)
lons = (/ 8.5, 9.0 /) 

nlocs = dimsizes(locs)    

nlocs = dimsizes(locs)      
loc = wrf_user_ll_to_ij(a, lons, lats, True)  ; use faster functions in 6.6.2

;================================================================================
; TIME and DATES

times  = wrf_user_list_times(a)  
ntimes = dimsizes(times)  

do it = 0, ntimes-1, 1

print("Working on: " + times(it))

namday = new(dimsizes(times(it)),string)	
date = stringtochar(times(it))			
yy = date(0:3)			
year = charactertostring(yy)	
yr = stringtointeger(year)			
mm = date(5:6)
month = charactertostring(mm)
mn = stringtointeger(month)			
dd = date(8:9)
day = charactertostring(dd)
da = stringtointeger(day)				
hh = date(11:12)
hour = charactertostring(hh)
hr = stringtointeger(hour) 			
;
dow = day_of_week(yr,mn,da) 
if (dow.eq.0) then
namday = "Sunday"
end if
if (dow.eq.1) then
namday = "Monday"
end if 
if (dow.eq.2) then
namday = "Tuesday"
end if
if (dow.eq.3) then
namday = "Wednesday"
end if
if (dow.eq.4) then
namday = "Thursday"
end if
if (dow.eq.5) then
namday = "Friday"
end if
if (dow.eq.6) then
namday = "Saturday"
end if
;
if (mn.eq.1) then			
namon = "Jan"
end if
if (mn.eq. 2) then
namon = "Feb"
end if
if (mn.eq. 3) then
namon = "Mar"
end if  
if (mn.eq. 4) then
namon = "Apr"
end if
if (mn.eq. 5) then
namon = "May"
end if
if (mn.eq. 6) then
namon = "Jun"
end if
if (mn.eq. 7) then
namon = "Jul"
end if
if (mn.eq. 8) then
namon = "Aug"
end if
if (mn.eq. 9) then
namon = "Sep"
end if
if (mn.eq. 10) then
namon = "Oct"
end if
if (mn.eq. 11) then
namon = "Nov"
end if
if (mn.eq. 12) then
namon = "Dec"
end if

;================================================================

   tc = wrf_user_getvar(a, "tc", it)
   td = wrf_user_getvar(a, "td", it)
theta = wrf_user_getvar(a, "theta", it)
    p = wrf_user_getvar(a, "pressure", it)
    z = wrf_user_getvar(a, "z", it)
uvmet = wrf_user_getvar(a, "uvmet", it)
    u = uvmet(0,:,:,:) * 1.94386
    v = uvmet(1,:,:,:) * 1.94386


  do iloc = 0, nlocs-1
  print("=================================================")
  print("Plotting Skew-T and Hodograph for " + locs(iloc))

;----------------------------------------------------------------------------

  locX = loc(0, iloc) - 1
  locY = loc(1, iloc) - 1

; For hodo:
ums = u/1.94386 ; [m/s] 
vms = v/1.94386

wspd = sqrt(ums(:, locY, locX)*ums(:, locY, locX) + vms(:, locY, locX)*vms(:, locY, locX)) 
wspd01 = wspd(::-1)

wdir = wind_direction(ums(:, locY, locX),vms(:, locY, locX),1)
wdir1 = wdir(::-1)

z_h = z(:, locY, locX)
z1 = z_h(::-1)

;print(z1(38))                ; check lowest model level height.
;exit

p1  = p(::-1, locY, locX)
p11  = p(:, locY, locX)

tc1  = tc(::-1, locY, locX)
tc11  = tc(:, locY, locX)

tdc1  = td(::-1, locY, locX)
tdc11  = td(:, locY, locX)
tdk = td +273.15

thte1 = theta(:, locY, locX)

pascal = p*100

  mix  = mixhum_ptd(pascal(:, locY, locX),tdk(:, locY, locX),1)*1000    ; mixing ratio [g/kg] 
  mix1 = mix(::-1)

;============================================================================
; wks = gsn_open_wks("png","./WRF-SkewT-" + locs(iloc) + "-T" + it)
  
type = "png"
type@wkWidth = 2600  
type@wkHeight = 1800 

outfile = "WRF-TSkew-Hodo_-" + locs(iloc) + "-T"+it

wks = gsn_open_wks(type,outfile)

gsn_define_colormap(wks,"WhViBlGrYeOrReWh")
print("				")
;=======================================================================

  resSkew = True
  resSkew@DrawHeightScale = True
  resSkew@DrawHeightScaleFt = False
  resSkew@DrawStandardAtm = True
  resSkew@DrawFahrenheit = False
  resSkew@DrawColAreaFill = True
  resSkew@PrintOpts = False
  resSkew@DrawFahrenheit  = False   ; default is True
  resSkew@DrawColAreaColor  = 3      ; Light Grey for WhViBlGrYeOrReWh color table
  resSkew@lineThicknessMixRatio     = 6.0

  ;resSkew@tiMainString =    locs(iloc) + ", valid for " + namday + " " + da + " " \
   ;                      + namon + " " + yr + " " + hour + "Z"
  resSkew@tiMainFontHeightF    = 0.012
  resSkew@tiMainFont    = "helvetica-bold"

  resSkew@tmXBLabelFontHeightF    = 0.0216
  resSkew@tmYLLabelFontHeightF    = 0.026 

 resSkew@vpXF = 0.046
 resSkew@vpYF = 0.9
 resSkew@vpWidthF = 0.43
 resSkew@vpHeightF = 0.62

 
; Get the skew-T background
  skewt_bkgd = skewT_BackGround(wks, resSkew)
  draw(skewt_bkgd)
  
; Set resources for Skew-T data
  resData = True
  resData@Parcel = 1
  resData@WspdWdir = False
  resData@HspdHdir = True
  resData@PlotWindH = False
  resData@ThermoInfo = False
  resData@DrawColAreaFill = True    ; default is False
  resData@DrawMixRatio    = True 
  resData@Wthin     = 1        ; plot every n-th wind barb
  resData@colTemperature  = "red"
  resData@colDewPt        = "chartreuse3"
  resData@colCape         = "chocolate3"
  resData@colWindP        = "blue"
  resData@colWindZ        = "black"                                        
  resData@colWindH        = "black" 
 
  resData@lineThicknessDewPt = 5.5
  resData@lineThicknessTemperature = 5.5
  resData@linePatternCape    = 17
  resData@xpWind             = 45

;=========================================================================
; Calculate Thermodynamic info

z1000_1 = wrf_user_intrp3d(z,p,"h",1000,0.,False)  ; use faster functions in 6.6.2
z500_1  = wrf_user_intrp3d(z,p,"h",500,0.,False)
z700_1  = wrf_user_intrp3d(z,p,"h",700,0.,False)
z850_1  = wrf_user_intrp3d(z,p,"h",850,0.,False)
z600_1  = wrf_user_intrp3d(z,p,"h",600,0.,False)

;print(p(::-1,locY, locX))

z1000 = z1000_1(locY, locX)
z500 = z500_1(locY, locX)
z700 = z700_1(locY, locX)
z850 = z850_1(locY, locX)
z600 = z600_1(locY, locX)

t500_1  = wrf_user_intrp3d(tc,p,"h",500,0.,False)  ; use faster functions in 6.6.2
t700_1  = wrf_user_intrp3d(tc,p,"h",700,0.,False)
t850_1  = wrf_user_intrp3d(tc,p,"h",850,0.,False)
tc500 = t500_1(locY, locX)
tc700 = t700_1(locY, locX)
tc850 = t850_1(locY, locX)

;------------------------------------------------------------------------------------

 fmsg = default_fillvalue(typeof(tc))   ; get default missing 

 nlvls= dimsizes(p)
 plcl = fmsg             ; p (hPa) Lifting Condensation Lvl (lcl)
 tlcl = fmsg             ; temperature (C) of lcl

 ptlcl_skewt(p(0, locY, locX),tc(0, locY, locX),td(0, locY, locX),plcl,tlcl) 

 plcl1 = decimalPlaces(plcl,0,True)
 tlcl1 = decimalPlaces(tlcl,1,True)
 print("LCL level: "+plcl1+" hPa / "+tlcl1+" oC")
      
  shox = showal_skewt(p(:, locY, locX),tc(:, locY, locX),td(:, locY, locX))           ; Showwalter Index
 shox1 = decimalPlaces(shox,1,True)
 print("Showalter Index: "+shox1+"")

  pwat = pw_skewt(td(:, locY, locX),p(:, locY, locX))                  ; precipitable water (cm)
 pwat1 = decimalPlaces(pwat,1,True)
 print("Precipitable water: "+pwat1*10+" mm")

  iprnt = 0   ; debug only (>0)
   nlLcl= 0   ; LCL                          
   nlLfc= 0   ; LFC
 nlCross= 0   ; EL

 cape = cape_thermo(p(:, locY, locX),tc(:, locY, locX),plcl,iprnt)       ; MUCAPE (J/kg)

  tpar = cape@tparcel                       ; temp of the parcel
  nlLcl= cape@jlcl
  nlLfc= cape@jlfc
  nlCross= cape@jcross                      ; Equillibrium Level index

nlLfc = where(nlLfc.ge.0,nlLfc,1)

                                            ; Equillibrium Level
EL_p = p(nlCross, locY, locX)               ; hPa
EL_z = z(nlCross, locY, locX) /1000         ; in km
EL = decimalPlaces(EL_z,1,True)
print("EL = "+EL+" km")

LFC_p = p(nlLfc, locY, locX)                 ; LFC in hPa
LFC_z = z(nlLfc, locY, locX) /1000           ; in km
LFC = decimalPlaces(LFC_z,2,True)
print("LFC = "+LFC+" km")

 cape1 = decimalPlaces(cape,0,True)
 print("MUCAPE: "+cape1+" J/kg")

;------ 
cape2d2 = wrf_user_getvar(a,"cape_2d",it) 
cape3 = cape2d2(0,locY, locX) ; MU
  ;mcape = cinfo(0,locY, locX)

  mcin3  = cape2d2(1,locY, locX)
 ;mcin3  = where(mcin3.gt.500.or.mcin3.lt.0,1,mcin3)
  mcin3  = decimalPlaces(mcin3,0,True)

;print(mcin3)
;exit

   lcl3  = cape2d2(2,locY, locX)/1000
   lcl3 =decimalPlaces(lcl3,1,True)
   lfc3  = cape2d2(3,locY, locX)/1000
   lfc3 =decimalPlaces(lfc3,1,True)
;cape3@description = "MUCAPE"
;cape3@units = "J kg~S~-1~N~"


;print(z1(38))
;print(p1(38))
;cinfo= wrf_cape_3d(p(:, locY, locX), tc(:, locY, locX), mix*1000, z(:, locY, locX), 100, 1010, True)
;print(cinfo(0,:))
;print(cape3)
;print(mcin3)
;print(lcl3)
;print(lfc3)
;exit
;-------------
;-------------------------------------------------------------------------------------------
i0 = closest_val(10.,z_h)
u0 = ums(i0, locY, locX)  
v0 = vms(i0, locY, locX)

i6 = closest_val(6000.,z_h)
u6 = ums(i6, locY, locX)               
v6 = vms(i6, locY, locX) 

i3 = closest_val(3000.,z_h)
u3 = ums(i3, locY, locX)               
v3 = vms(i3, locY, locX)

i1 = closest_val(1000.,z_h)
u1 = ums(i1, locY, locX)               
v1 = vms(i1, locY, locX)  

u06 = u6 - u0
v06 = v6 - v0

u03 = u3 - u0
v03 = v3 - v0
shear_03 = sqrt(u03^2 + v03^2)
shear03 = decimalPlaces(shear_03,1,True)

u01 = u1 - u0
v01 = v1 - v0
shear_01 = sqrt(u01^2 + v01^2)
shear01 = decimalPlaces(shear_01,1,True)

shear_06 = sqrt(u06^2 + v06^2)
shear06 = decimalPlaces(shear_06,1,True)

BRN0 = (cape / (0.5*(shear06)^2)) 

BRN  = decimalPlaces(BRN0,1,True)               ; Bulk Richardson Number
print("BRN: "+BRN+" ") 

;------------------------------------------------------------------------
i850 = closest_val(850.,p1)
;tc850 = tc1(i850)
tdc850 = tdc1(i850)

i500 = closest_val(500.,p1)
tc500 = tc1(i500)

i700 = closest_val(700.,p1)
tc700 = tc1(i700)
tdc700 = tdc1(i700)

;K-Index = ( T850 - T500 ) + TD850 - ( T700 - TD700 )

KI0 = tc850 - tc500 + tdc850 - tc700 + tdc700
KI = decimalPlaces(KI0,0,True)           ; K-index (approximation at height)
print("K-index: "+KI+" oC") 

;TOTL= ( T850 - T500 ) + ( TD850 - T500 )

TT0 = tc850 + tdc850 - 2 * (tc500)
TT = decimalPlaces(TT0,0,True)           ; Total Totals (approximation at height)
print("TotalTotals: "+TT+" oC") 

;---------------------------------------------------------------------------------
;imix2 = closest_val(1000.,z_h) ; or check below

mix2 = mix(0:13)
mixing0 = avg(mix2)
mixing = decimalPlaces(mixing0,1,True)

print("Mixing ratio: "+mixing+" g/kg ")    ; Mean Mixing ratio 0-1000m (PBL)

;---------------------------------------------------------------------------------

;KO-Index = 0.5 * (ThetaE700 + ThetaE500) - 0.5 * (ThetaE1000 + ThetaE850)

ithe1  = wrf_user_intrp3d(theta,p,"h",980,0.,False)
ithe7  = wrf_user_intrp3d(theta,p,"h",700,0.,False)
ithe5  = wrf_user_intrp3d(theta,p,"h",500,0.,False)
ithe8  = wrf_user_intrp3d(theta,p,"h",850,0.,False)

;ithe7 = closest_val(700.,p1)
 ;the7 = thte1(ithe7)
the7 = ithe7(locY, locX) 

;ithe5 = closest_val(500.,p1)
 ;the5 = thte1(ithe5)
the5 = ithe5(locY, locX) 

;ithe8 = closest_val(850.,p1)
; the8 = thte1(ithe8)
the8 = ithe8(locY, locX) 

;ithe1 = closest_val(950.,p1)
 ;the1 = thte1(ithe1)
the1 = ithe1(locY, locX) 

ko0 = 0.5*(the7 + the5) - 0.5*(the1 + the8)
ko = decimalPlaces(ko0,1,True)

print("KO-Index: "+ko+" K")             ; KO-Index (smaller than 2 indicates storm probability)

;-------------------------------------------------------------------------------------
;SWEAT = 12*(TD850)+20*(TTI -49)+2*(WS850)+(WS500)+ 125*(sin(WD500-WD850)+ 0.2)

iws85 = closest_val(850.,p1)
ws850 = wspd01(iws85)              ; in knots
wd850 = wdir1(iws85)*3.14/180      ; in radians

iws50 = closest_val(500.,p1)
ws500 = wspd01(iws50)
wd500 = wdir1(iws50)*3.14/180      ; in radians

S0 = 125*(sin(wd500 - wd850) + 0.2)
S1 = where(wdir1(iws85).ge.130 .and. wdir1(iws85).le.250 .or. wdir1(iws50).ge.210 .and. wdir1(iws50).le.310, S0, 0)
S = where(S1.gt.0, S1, 0)

Tot = 20*(TT-49)
Tot = where(TT.lt.49,0,Tot)

Dew850 = 12*tdc850
Dew850 = where(Dew850.ge.1,Dew850,0)

sweat0 = Dew850 + Tot + 2*(ws850) + ws500 + S

sweat = decimalPlaces(sweat0,0,True)

print("SWEAT: "+sweat+" ")                  ; SWEAT index


;-------------------------------------------------------------------------------------

; Showalter index = t500 - tparcel500

g_adiab = 0.0095          ;g = DT/DZ = -9.5 C/km
   g_td = 0.002           ;g = DTd/DZ = -2.0 C/km
  g_wet = 0.0065          ;g = DT/DZ = -6.5 C/km

iz850 = closest_val(850.,p1)
iz500 = closest_val(500.,p1)
;z850 = z1(iz850)             ; in meters
;z500 = z1(iz500)

ip = closest_val(5650.,z_h)      ; approximation

tpar500 = cape@tparcel(ip)     ; temp of the parcel at 500hPa                   
;tparcel500 = tpar500

tparcel500 = - g_adiab*(z500 - z850) + tc850
;print(tparcel500) 
swi0 = tc500 - tparcel500
swi = decimalPlaces(swi0,0,True)            ; Showalter index (needs calculations)

print("Showalter index: "+swi+" ")

;-----------------------------------------------------------------------------------

iz0 = closest_val(500.,z_h)
tc0 = avg(tc11(0:iz0))       ; a parcel with avg Temp at the lowest 500m
z0 = z_h(iz0)

if (z0.le.600) then

tparcel500_li = - 0.0070*(z500 - z0) + tc0

;print(tc500)
;print(tparcel500_li)

LI0 = tc500 - tparcel500_li        ; Lifted Index (needs calculations, use with caution)
LI = decimalPlaces(LI0,1,True) 

print("Lifted index: "+LI+" oC")

else

LI = -999

print("Lifted index: -nan oC")

end if

;----------------------------------------------
;Thompson Index = KI - LI  (>35 interested)

Thompson = KI - LI
Thompson = where(Thompson.le.60,Thompson,KI)

print("Thompson index: "+Thompson+" oC")         ; Thompson Index 


;--------------------------------------------------------------------------------------
;850-600 Lapse Rate 

iz6 = closest_val(600.,p1)
tc600 = tc1(iz6)
;z600 = z1(iz6)

g0_68 = -((tc600 - tc850)/(z600 - z850))*1000   ; in km

g_68 = decimalPlaces(g0_68,1,True) 

print("850 - 600 hPa lapse rate: "+g_68+" oC/km")        ;850-600 Lapse Rate

;--------------------------------------------------------------------------------------
;TQ = (T850 + Td850) � 1.7?(T700). 

TQ0 = (tc850 + tdc850) - 1.7*(tc700)

TQ = decimalPlaces(TQ0,0,True)                         ; TQ Index for low-topped convection

print("TQ Index: "+TQ+" oC")

;-------------------------------------------------------------------------------------------
; Freezing point

frz_lvl1 = where(tc1.lt.1, 0, 1)
frz_lvl = maxind(frz_lvl1)

zfrz = z1(frz_lvl)                    ; Freezing level 
pfrz = p1(frz_lvl)

print("Freezing level: "+zfrz+" m / near "+pfrz+" hPa ")

;iwb = closest_val(0,tc1)
;twb = wetbulb(p1, tc1, tdc1)    ; ===> 9.3C
;print("Wet bulb: "+twb+" ")
;print("Wet bulb 08C: "+zwb0+" km / at "+pwb0+" hPa ")          ; Wet bulb 0C (needs more calculations.)
																; use ncl  v 6.6.2

z_surf = closest_val(6,z_h)
t_surf = tc11(z_surf)

print("Surface temp: "+t_surf+" oC")

;------------------------------------------------------------------------------------------
; Thickness 1000-500, 1000-850, 1000-700

th10_500 = (z500 - z1(38)) + 115    ; conversion due to lack of 500hPa layer
th10_850 = z850 - z1(38) + 115 ;- z1000
th10_700 = z700 - z1(38)

print("Thickness 1000-500: "+th10_500+" gpm")            ; Thickness
print("Thickness 1000-850: "+th10_850+" gpm")
print("Thickness 1000-700: "+th10_700+" gpm")

;print(z1(38))
;print(z1000)
;exit

;t850_1  = wrf_user_intrp3d(tc,p,"h",850,0.,False)
;t700_1  = wrf_user_intrp3d(tc,p,"h",700,0.,False)
;t500_1  = wrf_user_intrp3d(tc,p,"h",500,0.,False)

;t500 = t500_1(locY, locX)
;t700 = t700_1(locY, locX)
;t850 = t850_1(locY, locX)

print("T850: "+tc850 +" oC")            ; Temperature
print("T700: "+tc700 +" oC")
print("T500: "+tc500 +" oC")

;------------------------------------------------------------------------------------------
;******************************************************************************************************

; Storm Relative Helicity

; SREH0-3km & SREH0-1km

pp0 = closest_val(0.,z_h)
u1000 = ums(pp0, locY, locX)  
v1000 = vms(pp0, locY, locX)

pp1 = closest_val(400.,z_h)
u950 = ums(pp1, locY, locX)  
v950 = vms(pp1, locY, locX)

pp2 = closest_val(900.,z_h)
u900 = ums(pp2, locY, locX)  
v900 = vms(pp2, locY, locX)

pp3 = closest_val(1450.,z_h)
u850 = ums(pp3, locY, locX)  
v850 = vms(pp3, locY, locX)

pp4 = closest_val(1900.,z_h)
u800 = ums(pp4, locY, locX)  
v800 = vms(pp4, locY, locX)

pp5 = closest_val(2200.,z_h)
u750 = ums(pp5, locY, locX)  
v750 = vms(pp5, locY, locX)

pp6 = closest_val(3000.,z_h)
u700 = ums(pp6, locY, locX)  
v700 = vms(pp6, locY, locX)

pp7 = closest_val(3600.,z_h)
u650 = ums(pp7, locY, locX)  
v650 = vms(pp7, locY, locX)

pp8 = closest_val(4000.,z_h)
u600 = ums(pp8, locY, locX)  
v600 = vms(pp8, locY, locX)

pp9 = closest_val(4800.,z_h)
u550 = ums(pp9, locY, locX)  
v550 = vms(pp9, locY, locX)

pp10 = closest_val(5600.,z_h)
u500 = ums(pp10, locY, locX)  
v500 = vms(pp10, locY, locX)

pp11 = closest_val(6100.,z_h)
u450 = ums(pp11, locY, locX)  
v450 = vms(pp11, locY, locX)

pp12 = closest_val(7300.,z_h)
u400 = ums(pp12, locY, locX)  
v400 = vms(pp12, locY, locX)

umean=(u1000+u950+u900+u850+u800+u750+u700+u650+u600+u550+u500+u450+u400)/13.0
vmean=(v1000+v950+v900+v850+v800+v750+v700+v650+v600+v550+v500+v450+v400)/13.0
ushear=u500-u1000
vshear=v500-v1000
shear=sqrt(ushear*ushear+vshear*vshear)
umotion=((umean+(7.5/(shear))*vshear))
vmotion=((vmean-(7.5/(shear))*ushear))

srh1=((u950-umotion)*(v1000-vmotion)-(u1000-umotion)*(v950-vmotion))
srh2=((u900-umotion)*(v950-vmotion)-(u950-umotion)*(v900-vmotion))
srh3=((u850-umotion)*(v900-vmotion)-(u900-umotion)*(v850-vmotion))
srh4=((u800-umotion)*(v850-vmotion)-(u850-umotion)*(v800-vmotion))
srh5=((u750-umotion)*(v800-vmotion)-(u800-umotion)*(v750-vmotion))
srh6=((u700-umotion)*(v750-vmotion)-(u750-umotion)*(v700-vmotion))

srh3km = srh1+srh2+srh3+srh4+srh5+srh6
srh1km = srh1+srh2

srh3km = decimalPlaces(srh3km,1,True) 
srh1km = decimalPlaces(srh1km,1,True) 

print("SREH0-3: "+srh3km+" m2/s2")
print("SREH0-1: "+srh1km+" m2/s2")

;------------------------------------------------------------------------
;EHI

EHI0 = (cape*srh3km)/160000                      ; Energy Helicity Index
EHI  = decimalPlaces(EHI0,0,True)
print("EHI: "+EHI+" ") 
;------------------------------------------------------------------------

;Significant tornado parameter

stp = (cape/1500)*((2000-plcl)/1000)*(srh1km/150)*(shear06/20)

stp  = decimalPlaces(stp,0,True)

;print(stp)

;------------------------------------------------------------------------

;Supercell Composite Parameter

SCP = (cape / 1000) * (srh3km / 100) * (BRN0 / 40)

SCP  = decimalPlaces(SCP,0,True)

print("SCP: "+SCP+"")


;===================================================
;Plot

;------------------------------------------------------------------------------------
  skewT_data = skewT_PlotData(wks, skewt_bkgd, p(:, locY, locX), \
                              tc(:, locY, locX), \
                              td(:, locY, locX), \
                              z(:, locY, locX), \
                              u(:, locY, locX), \
                              v(:, locY, locX), \
                              resData)
  ;frame(wks) 

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

;**************************************************************************

;---------------------------------- Hodograph area ------------------------
;==========================================================================
  ;wks_type = "png" 
  ;wks_type@wkWidth = 1600  
  ;wks_type@wkHeight = 1600

;outhodo = "WRF-Hodograph_GR-6km-" + locs(iloc) + "-T"+it

;xwks  = gsn_open_wks (wks_type ,outhodo)

  xres                      = True 
  ;xres@gsnMaximize          = True    ; Maximize plot in window.
  xres@gsnDraw              = False   ; Don't draw plot
  xres@gsnFrame             = False   ; Don't advance frame

  xres@gsnLeftString        = ""
  xres@gsnRightString       = ""
  xres@xyMarkLineMode       = "MarkLines"  ; markers and lines
  xres@xyMarkers            = 1           ; filled dot
  ;xres@xyMarkerSizeF        = 0.02
  ;xres@xyMarkerThicknessF   = 0.55

 ;xres@tiMainString         = "16622 LGTS Thessaloniki at "+hour+":00 UTC - "+day+" "+mon+" "+year+""
 ; xres@tiMainString  = locs(iloc) + ", valid for " + namday + " " + da + " " \
            ;             + namon + " " + yr + " " + hour + "Z"

  xres@tiMainFontHeightF    = 0.007
  xres@tiMainFont           = "helvetica-bold"
  xres@tiYAxisString        = "V-Wind"
  xres@tiXAxisString        = "U-Wind"
  
  ;popts          = True                 ; local and default options
  ;circFr = popts@circFr
  ;SpdMax = popts@SpdMax
  ;zMax   = popts@zMax

  ;SpdMaxNumCircles = floattointeger(SpdMax/circFr)
  ;SpdMaxNum        = circFr*SpdMaxNumCircles   ; number of circles

  ;extraSpace             = max((/3.,circFr/3./))    ; Extra space beyond
  ;xres@trXMinF           = -SpdMaxNum-extraSpace  ; min X 
  ;xres@trXMaxF           =  SpdMaxNum+extraSpace  ; max X
  ;xres@trYMinF           = -SpdMaxNum-extraSpace  ; min Y
  ;xres@trYMaxF           =  SpdMaxNum+extraSpace  ; max Y

  xres@trXMinF	           = -40 ;max( wspd ) ;
  xres@trXMaxF	           =  40 ;max( wspd ) ;  
  xres@trYMinF             = -40 ;max( wspd ) ;
  xres@trYMaxF             =  40 ;max( wspd ) ; 

  xres@vpWidthF  = 0.25
  xres@vpHeightF = 0.25
  xres@vpXF     = 0.60 
  xres@vpYF     = 0.9 

  xres@hodo_labelCircleFontHeightF = 0.006 


;************************************************************************************

;================== TSkew plot ========================
 txress = True                        
 txress@txFontColor       = "Black"
 txress@txFontHeightF     = 0.010
 txress@txFont      = "helvetica"

 gsn_text_ndc(wks,"MUCAPE = "+cape1+" J/kg",0.595,.55,txress)
 gsn_text_ndc(wks,"PBL mixing = "+mixing+" g/kg ",0.60,.53,txress)
 gsn_text_ndc(wks,"LCL hgt = "+lcl3+" km ",0.596,.51,txress)
 gsn_text_ndc(wks,"LFC hgt = "+lfc3+" km",0.597,.49,txress)
 gsn_text_ndc(wks,"PWAT = "+pwat1*10+" mm",0.592,.47,txress)
 gsn_text_ndc(wks,"TQ index = "+TQ+" ~S~o~N~C",0.592,.45,txress)
 gsn_text_ndc(wks,"K-index = "+KI+" ~S~o~N~C",0.592,.43,txress)
 gsn_text_ndc(wks,"EL = "+EL+" km",0.597,.41,txress)
 gsn_text_ndc(wks,"SWI = "+shox1+" ~S~o~N~C",0.592,.39,txress)
 gsn_text_ndc(wks,"SWEAT = "+sweat+"",0.592,.37,txress)
 gsn_text_ndc(wks,"KO-index = "+ko+" K",0.592,.35,txress)
 gsn_text_ndc(wks,"Total Totals = "+TT+" ~S~o~N~C",0.592,.33,txress)
 gsn_text_ndc(wks,"Lifted index = "+LI+" ~S~o~N~C",0.592,.31,txress)
 gsn_text_ndc(wks,"LR 850-600 = "+g_68+" ~S~o~N~C/km",0.592,.29,txress)

 gsn_text_ndc(wks,"Thompson index = "+Thompson+" ~S~o~N~C",0.59,.27,txress)

 txrex = True                        
 txrex@txFontColor       = "Blue"
 txrex@txFontHeightF     = 0.010
 txrex@txFont      = "helvetica"

tc850 = decimalPlaces(tc850,1,True)
tc700 = decimalPlaces(tc700,1,True)
tc500 = decimalPlaces(tc500,1,True)
th10_500 = decimalPlaces(th10_500,0,True)
th10_850 = decimalPlaces(th10_850,0,True)
th10_700 = decimalPlaces(th10_700,0,True)
zfrz = decimalPlaces(zfrz,0,True)

 gsn_text_ndc(wks,"Th1000-500 = "+th10_500+" gpm",0.79,.55,txrex)
 gsn_text_ndc(wks,"Th1000-850 = "+th10_850+" gpm",0.79,.53,txrex)
 gsn_text_ndc(wks,"Th1000-700 = "+th10_700+" gpm",0.79,.51,txrex)
 gsn_text_ndc(wks,"Freezing level = "+zfrz+" m",0.79,.49,txrex)
 gsn_text_ndc(wks,"T850 = "+tc850+" ~S~o~N~C",0.79,.47,txrex)
 gsn_text_ndc(wks,"T700 = "+tc700+" ~S~o~N~C",0.79,.45,txrex)
 gsn_text_ndc(wks,"T500 = "+tc500+" ~S~o~N~C",0.79,.43,txrex)

 txrez = True                        
 txrez@txFontColor       = "Chocolate4"
 txrez@txFontHeightF     = 0.010
 txrez@txFont     = "helvetica"

 gsn_text_ndc(wks,"Shear 0-6km = "+shear06+" m/s",0.79,.41,txrez)
 gsn_text_ndc(wks,"Shear 0-3km = "+shear03+" m/s",0.79,.39,txrez)
 gsn_text_ndc(wks,"Shear 0-1km = "+shear01+" m/s",0.79,.37,txrez)
 gsn_text_ndc(wks,"SREH0-3km = "+srh3km+" m~S~2~N~/s~S~2~N~ ",0.79,.35,txrez)
 gsn_text_ndc(wks,"SREH0-1km = "+srh1km+" m~S~2~N~/s~S~2~N~ ",0.79,.33,txrez)

 gsn_text_ndc(wks,"EHI0-3km = "+EHI+" ",0.79,.31,txrez)
 gsn_text_ndc(wks,"STP = "+stp+" ",0.79,.29,txrez)
 gsn_text_ndc(wks,"BRN = "+BRN+" ",0.79,.27,txrez)
 gsn_text_ndc(wks,"SCP = "+SCP+" ",0.79,.25,txrez)
 

 txrer = True                        
 txrer@txFontColor       = "Black"
 txrer@txFontHeightF     = 0.010
 txrer@txFont     = "helvetica-bold"


;gsn_text_ndc(wks,"   Macedonia airport",0.90,.16,txrer)
;gsn_text_ndc(wks,"   Station latitude: 40~S~o~N~.51~S~'~N~ ",0.90,.14,txrer)
;gsn_text_ndc(wks,"Station longitude: 22~S~o~N~.96~S~'~N~ ",0.90,.12,txrer)
;gsn_text_ndc(wks,"Station elevation: 4 m ",0.89,.10,txrer)


 txrer1 = True                        
 txrer1@txFontColor       = "Black"
 txrer1@txFontHeightF     = 0.011
 txrer1@txFont     = "helvetica-bold"
gsn_text_ndc(wks," Thermodynamic parameters ",0.70,.58,txrer1)

;==================================================================

;------------------------------------------------------------------------
;Initial time
namdayi = new(dimsizes(times(0)),string)	
datei = stringtochar(times(0))			
yyi = datei(0:3)			
yeari = charactertostring(yyi)	
yri = stringtointeger(yeari)			
mmi = datei(5:6)
monthi = charactertostring(mmi)
mni = stringtointeger(monthi)			
ddi = datei(8:9)
dayi = charactertostring(ddi)
dai = stringtointeger(dayi)				
hhi = datei(11:12)
houri = charactertostring(hhi)
hri = stringtointeger(houri) 

dowi = day_of_week(yri,mni,dai) 
if (dowi.eq.0) then
namdayi = "Sun"
end if
if (dowi.eq.1) then
namdayi = "Monday"
end if 
if (dowi.eq.2) then
namdayi = "Tuesday"
end if
if (dowi.eq.3) then
namdayi = "Wednesday"
end if
if (dowi.eq.4) then
namdayi = "Thursday"
end if
if (dowi.eq.5) then
namdayi = "Friday"
end if
if (dowi.eq.6) then
namdayi = "Saturday"
end if
;
if (mni.eq.1) then			
namoni = "Jan"
end if
if (mni.eq. 2) then
namoni = "Feb"
end if
if (mni.eq. 3) then
namoni = "Mar"
end if  
if (mni.eq. 4) then
namoni = "Apr"
end if
if (mni.eq. 5) then
namoni = "May"
end if
if (mni.eq. 6) then
namoni = "Jun"
end if
if (mni.eq. 7) then
namoni = "Jul"
end if
if (mni.eq. 8) then
namoni = "Aug"
end if
if (mni.eq. 9) then
namoni = "Sep"
end if
if (mni.eq. 10) then
namoni = "Oct"
end if
if (mni.eq. 11) then
namoni = "Nov"
end if
if (mni.eq. 12) then
namoni = "Dec"
end if

txtres1 = True
txtres1@txFontHeightF = 0.009
txtres1@txJust = "CenterRight"
txtres1@txFont  = "helvetica-bold"

;gsn_text_ndc(wks,"Init.:  " + namdayi + " " + dai + " " + namoni + " " + yri + \
;              " at" + " " + houri + "Z",0.90,0.92,txtres1)

;==================================================

gsn_text_ndc(wks,"WRF-ARW output - Init. NCEP: " + namdayi + " " + dai + " " + namoni + " " + yri + \
	" at" + " " + hri + "Z",0.33,0.215,txtres1)
;gsn_text_ndc(wks,"***info*** ",0.89,0.215,txtres1)

txtres11 = True
txtres11@txFontHeightF = 0.015
txtres11@txJust = "CenterRight"
txtres11@txFont  = "helvetica-bold"
gsn_text_ndc(wks,locs(iloc) + ", valid for " + namday + " " + da + " " \
                         + namon + " " + yr + " " + hour + "Z",0.65,0.925,txtres11)


plot = hodograph(wks,wspd,wdir,z_h,xres)

;================================= Hodograph ===============================
 draw (plot)

frame(wks)

print("=========================================")
print(" ")
     
        system("convert -trim +repage " + outfile + ".png " + outfile + ".trimmed.png")   ; trim the image
        system("mv " + outfile + ".trimmed.png " + outfile + ".png")


  delete(wks)
  delete(resSkew)
  delete(resData)
  delete(skewT_data)
  delete(skewt_bkgd)
  delete(z1000)
  delete(z500)
  delete(z850)
  delete(z700)
  delete(tc500)
  delete(tc850)
  delete(tc700)


  	print("				")
end do
	print("				")
	print("				")
	print("		 Next time step		")
	print("		****************		")
	print("				")
 print("CPU time consumed: " + get_cpu_time()+" seconds")
end do

	print("				")
	print("		****************		")
	print(" Finish computations - script version: 2017")
	print("		****************		")
	print("				")

	print("Total CPU time: " + get_cpu_time()+" seconds")

end