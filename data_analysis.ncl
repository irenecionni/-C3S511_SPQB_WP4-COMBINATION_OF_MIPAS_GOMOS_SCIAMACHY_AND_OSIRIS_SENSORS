load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
;************************************************
undef("coswgt_area_average")
function coswgt_area_average(field:numeric)
;
; Arguments
;    field: numeric field.
;
; Return value
;    The area average using cosine lat weights.
;
; Description
;    Computes the area average using cosine lat weights and lon weights=1.
;
local funcname, scriptname, lat, wgt_lat, lon, lon_size, wgt_lon, ave
begin

  funcname = "coswgt_area_average"
  scriptname = "diag_scripts/shared/statistics.ncl"

  lat = field&lat
  wgt_lat = tofloat(NormCosWgtGlobe(lat))

  lon = field&lon
  lon_size = dimsizes(lon)
  wgt_lon = new((/lon_size(0)/), float)
  wgt_lon = 1.0

  ave = wgt_areaave_Wrap(field, wgt_lat, wgt_lon, 0)
  return(ave)

end
; #############################################################################
undef("coswgt_area_std")
function coswgt_area_std(field:numeric)
;
; Arguments
;    field: numeric field.
;
; Return value
;    The area average using cosine lat weights.
;
; Description
;    Computes the area average using cosine lat weights and lon weights=1.
;
begin

  funcname = "coswgt_area_std"

   lat = field&lat
   lon = field&lon
   re   = 6.37122e06
   rad  = 4.0 * atan(1.0) / 180.0
   con  = re * rad
   clat = cos(lat * rad)           ; cosine of latitude
   dlon = (lon(2) - lon(1))        ; assume dlon is constant
   dlat = (lat(2) - lat(1))        ; assume dlat is constant

   dx   = con * dlon * clat        ; dx at each latitude
   dy   = con * dlat               ; dy is constant
   dydx = dy * dx                  ; dydx(nlat)
   nlat=dimsizes(lat)
   mlon=dimsizes(lon)
   wgt  = new((/nlat, mlon/), typeof(dydx))
   wgt=conform (wgt, dydx, 0)



   ave = wgt_areaave2(field, wgt, 0)
   d2 = (field - ave) ^ 2 ;
     ; Calculate weighted standard deviation
   ss=sum(d2)
   arg = wgt_areasum2(d2, wgt, 0)

   v1 = sum(where(ismissing(field), wgt@_FillValue, wgt))
   v2 = sum(where(ismissing(field), wgt@_FillValue, wgt) ^ 2)
   dd=(v1 / (v1 ^ 2 - v2))

   out = sqrt( arg*dd)

  return(out)


end
    ; Pi [1]
    PI = 3.1415927
    PI@units = "1"
    PI@name = "Greek Pi"


    EARTHradius = 6371000.
    EARTHradius@units = "m"
    EARTHradius@name = "Earth radius"
    ; Earth surface area [m^2]
  
    EARTHarea = 4. * PI * EARTHradius ^ 2
    EARTHarea@units = "m2"
    EARTHarea@name = "Earth surface area"
; #############################################################################
undef("gridcell_area")
function gridcell_area(deltax[1]: numeric,
                       lat_lo[1]: numeric,
                       lat_hi[1]: numeric)
local funcname, scriptname, deg2rad
begin

  funcname = "gridcell_area"
  scriptname = "diag_scripts/shared/latlon.ncl"
  ;enter_msg(scriptname, funcname)

  ; Convert to radiants
  deg2rad = acos(-1.0) / 180.

  deltay = sin(lat_hi * deg2rad) - sin(lat_lo * deg2rad)
  out = abs(EARTHradius ^ 2 * deltay * deltax * deg2rad)
  out@units = EARTHarea@units

  ;leave_msg(scriptname, funcname)
  return(out)

end

; #############################################################################
; #############################################################################
undef("map_area")
function map_area(lat[*]:numeric,
                  lon[*]:numeric)
local funcname, scriptname, loclat, loclon, nlat, lati, xresol
begin

  funcname = "map_area"
  scriptname = "diag_scripts/shared/latlon.ncl"
  ;enter_msg(scriptname, funcname)

  ; Type conversion
  loclat = tofloat(lat)
  loclon = tofloat(lon)
  nlat = dimsizes(loclat)

  ; Check ordering (must be S-N)
  if (isMonotonic(loclat).ne.1) then
    ;error_msg("f", scriptname, funcname, "input latitude is not " + \
;              "monotonically increasing")
  end if

  ; Calculate longitude resolution (assumed constant)
  xresol = abs(loclon(1) - loclon(0))
  dsize = dimsizes(loclon)
  if (dsize(0).gt.2) then
    if ((xresol - abs(loclon(2) - loclon(1))).gt.1.0e-4) then
      ;error_msg("f", scriptname, funcname, "input longitude " + \
 ;               "resolution is not constant")
    end if
  end if

  ; Calculate latitude interfaces
  lati = new(dimsizes(loclat) + 1, float)
  lati(0) = max((/(3 * loclat(0) - loclat(1)) / 2., -90./))
  do ii = 1, dimsizes(lati) - 2
    lati(ii) = 0.5 * (loclat(ii - 1) + loclat(ii))
  end do
  lati(dimsizes(lati) - 1) = \
    min((/(3 * loclat(nlat - 1) - loclat(nlat - 2)) / 2., 90./))

  ; Calculate area
  area = new((/dimsizes(loclat), dimsizes(loclon)/), float)
  do ii = 0, dimsizes(loclat) - 1
    area(ii, :) = gridcell_area(xresol, lati(ii), lati(ii + 1))
  end do

  ; Assign coordinates
  area!0   = "lat"
  area&lat = lat
  area!1   = "lon"
  area&lon = lon

  ;leave_msg(scriptname, funcname)
  return(area)



end
;************************************************
function area_operations(field:numeric,
                         latmin[1]:numeric,
                         latmax[1]:numeric,
                         lonmin[1]:numeric,
                         lonmax[1]:numeric,
                         opt[1]:string,
                         l_wgt[1]:logical)

local funcname, scriptname, verbosity, rank, dimnames, wgt, conwgt, f_latmin, \
    f_latmax, f_lonmin, f_lonmax, delta, idx_latmin, idx_latmax, latidx, \
    idx_lonmin, idx_lonmax, lonidx, subfield, subwgt, dnames, didx_r, didx_l, \
    norm
begin

    funcname = "area_operations"
    scriptname = "diag_scripts/lib/ncl/latlon.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    ;enter_msg(scriptname, funcname, 6)

    ;; Check rank
    if (dimsizes(dimsizes(field)) .lt. 2) then
       ; error_msg("f", scriptname, funcname, "rank of input field <2")
    end if

    ;; Check dimensions
    rank = dimsizes(dimsizes(field))
  dimnames = getvardims(field)
    if (dimnames(rank - 1) .ne. "lon") then
        ;error_msg("f", scriptname, funcname, "the rightmost dimension " + \
              ;    "of the input field must be lon")
    end if
    if (dimnames(rank - 2) .ne. "lat") then
        ;error_msg("f", scriptname, funcname, "the second rightmost " + \
               ;   "dimension of the input field must be lat")
    end if
    delete(dimnames)

    ;; Check that grid is regular
    if (dimsizes(dimsizes(field&lat)) .ne. 1 .or. \
        dimsizes(dimsizes(field&lon)) .ne. 1) then
        ;error_msg("f", scriptname, funcname, "irregular grids not supported")
    end if

    ;; Calculate weights
    wgt = map_area(field&lat, field&lon)
    if (.not.l_wgt) then
        wgt = 1.
    end if

    ;; Conform weights and filter for missing values in the input field
    if (rank.eq.4) then
        conwgt = conform(field, wgt, (/2, 3/))
    end if
    if (rank.eq.3) then
        conwgt = conform(field, wgt, (/1, 2/))
    end if
    if (rank.eq.2) then
        conwgt = wgt
    end if
    conwgt = where(ismissing(field), conwgt@_FillValue, conwgt)
    copy_VarCoords(field, conwgt)
    delete(wgt)

    ;; Copy arguments to local (swap lat limits if necessary)
    f_latmin = todouble(min((/latmin, latmax/)))
    f_latmax = todouble(max((/latmin, latmax/)))
    f_lonmin = todouble(where(lonmin .lt. 0., lonmin + 360., lonmin))
    f_lonmax = todouble(where(lonmax .lt. 0., lonmax + 360., lonmax))
; If single gridbox lon = 360, replace with 0 to account for periodicity
    if (f_lonmin .eq. 360 .and. f_lonmax .eq. 360) then
        f_lonmin = 0.d0
        f_lonmax = 0.d0
    end if

    ;; Find nearest points to latitude boundaries
    ;; In the (unlikely) case that the selected boundary point is at the
    ;; midpoint between two gridboxes, ind will return two points.
    ;; The minimum index (= minumum value) is chosen in such case,
    ;; consistently with NCL coordinate subscripting
    ;; http://www.ncl.ucar.edu/Document/Language/subscript_coord.shtml
    delta = abs(field&lat - f_latmin)
    idx_latmin = min(ind(delta .eq. min(delta)))
    delta = abs(field&lat - f_latmax)
    idx_latmax = min(ind(delta .eq. min(delta)))
    delete(delta)
    if (ismissing(idx_latmin) .or. ismissing(idx_latmax)) then
        ;error_msg("f", scriptname, funcname, "no values found in the " + \
               ;   "given latitude range")
 end if
    if (latmin .eq. latmax .and. idx_latmin .ne. idx_latmax) then
        ;error_msg("f", scriptname, funcname, "multiple indexes found " + \
               ;   "for single gridbox latitude selection")
    end if

    ;; Define latitude indexes
    latidx = ispan(idx_latmin, idx_latmax, 1)  ; array of indexes
    delete(idx_latmin)
    delete(idx_latmax)

    ;; Find neartest points to longitude boundaries
    ;; In the (unlikely) case that the selected boundary point is at the
    ;; midpoint between two gridboxes, ind will return two points.
    ;; The minimum index (= minumum value) is chosen in such case,
    ;; consistently with NCL coordinate subscripting
    ;; http://www.ncl.ucar.edu/Document/Language/subscript_coord.shtml
    if (f_lonmin.gt.max(field&lon)) then
        idx_lonmin = 0
    else
        delta = abs(field&lon - f_lonmin)
        idx_lonmin = min(ind(delta .eq. min(delta)))
        delete(delta)
    end if
    if (f_lonmax.gt.max(field&lon)) then
        idx_lonmax = dimsizes(field&lon) - 1
    else
        delta = abs(field&lon - f_lonmax)
        idx_lonmax = min(ind(delta .eq. min(delta)))
        delete(delta)
    end if
    if (ismissing(idx_lonmin) .or. ismissing(idx_lonmax)) then
        ;error_msg("f", scriptname, funcname, "no values found in the " + \
                  ;"given longitude range")
    end if
    ;; Special case: single gridbox selection beteen max(field&lon) and 0
    if (lonmin .eq. lonmax .and. lonmin.gt.max(field&lon)) then
        idx_lonmin = where(abs(360. - lonmin) .lt. \
                           abs(max(field&lon) - lonmin), \
                           0, dimsizes(field&lon) - 1)
        idx_lonmax = idx_lonmin
    end if
    if (lonmin .eq. lonmax .and. idx_lonmin .ne. idx_lonmax) then
        ;error_msg("f", scriptname, funcname, "multiple indexes found " + \
                  ;"for single gridbox longitude selection")
    end if

    ;; Define longitude indexes
    if (idx_lonmin .le. idx_lonmax) then  ; Greenwich outside
        lonidx = ispan(idx_lonmin, idx_lonmax, 1)
    else  ; Greenwich inside
        lonidx = array_append_record( \
            ispan(idx_lonmin, dimsizes(field&lon) - 1, 1), \
            ispan(0, idx_lonmax, 1), 0)
    end if
    delete(idx_lonmin)
    delete(idx_lonmax)
;; Extract region
    if (rank .eq. 4) then
        subfield = field(:, :, latidx, lonidx)
        subwgt = conwgt(:, :, latidx, lonidx)
    end if
    if (rank .eq. 3) then
        subfield = field(:, latidx, lonidx)
        subwgt = conwgt(:, latidx, lonidx)
    end if
    if (rank .eq. 2) then
        subfield = field(latidx, lonidx)
        subwgt = conwgt(latidx, lonidx)
    end if
    delete(conwgt)

    if (opt .eq. "extract") then
        ;;leave _msg(scriptname, funcname, 6)
        return(subfield)
    end if

    ;; Find dimensions to be averaged or summed
    dnames = getVarDimNames(subfield)
    didx_r = ind(dnames.eq."lat" .or. dnames.eq."lon")
    didx_l = ind(dnames.ne."lat" .and. dnames.ne."lon")
    if (all(ismissing(didx_r))) then
        ;error_msg("f", scriptname, funcname, "cannot perform " + \
        ;          "average/sum on single gridcell, use the extract option")
    end if

    ;; Normalize weights
    if (opt .eq. "average") then
        if (dimsizes(didx_r).eq.dimsizes(dnames)) then
            subwgt = subwgt / sum(subwgt)
        else
            norm = dim_sum_n_Wrap(subwgt, didx_r)
            subwgt = subwgt / conform(subwgt, norm, didx_l)
        end if
    end if

    ;; Calculate average/sum. The same code applies, as the only difference
    ;; is in the normalization of the weights, which is calculated above
    ;; only for opt="average"
    if (opt .eq. "average" .or. opt .eq. "sum") then
        subfield = subfield * subwgt
        if (dimsizes(didx_r).eq.dimsizes(dnames)) then
            out = dim_sum_Wrap(dim_sum_Wrap(subfield))
        else
            out = dim_sum_n_Wrap(subfield, didx_r)
        end if
       ; ;leave _msg(scriptname, funcname, 6)
        return(out)
    end if
  ;error_msg("f", scriptname, funcname, "unrecognized option " + opt)

end
begin

;************************************************
; read in netCDF file
;************************************************
 dir="/Users/irenecionni/"
 a = addfile(dir+"/obs_c3s/Tier3/CDS-OZONE-L3-LP/OBS_sat_ESACCI-OZONE-LMZ-MERGED_L3_T3M_vmro3_1985-2017.nc" ,"r")
;************************************************
; read in zonal winds
;************************************************
  vmro3=a->vmro3(:,:,:,0)      

  vmro31=a->vmro3                          
  vmro3=(/vmro3*10^6/)
  vmro31=(/vmro31*10^6/)
  plev=a->plev
  lat=a->lat
  plev@units="hPa"
  plev@long_name="Pressure (hPa)" 
  vmro3&plev=plev
  vmro31&plev=plev
  dim=dimsizes(vmro3)
  yyyymm=  cd_calendar(vmro3&time, 1) 
  time=vmro3&time
  date= cd_calendar(time, 0)
  yy=toint(date(:,0)) 
  ts_avg=coswgt_area_average(vmro31)
  plev=vmro3&plev
  ts1_avg=new((/dimsizes(ts_avg&time),3/),typeof(ts_avg))
  ts1_avg(:,0)=ts_avg(:,{2.})
  ts1_avg(:,1)=ts_avg(:,{10.})
  ts1_avg(:,2)=ts_avg(:,{70.})   
  ts1_lat=new((/dimsizes(vmro3&time),dimsizes(lat),3/),typeof(ts_avg))
  ao=vmro3(:,{2.},:)
  ts1_lat(:,:,0)=vmro3(:,{2.},:)
  ts1_lat(:,:,1)=vmro3(:,{10.},:)
  ts1_lat(:,:,2)=vmro3(:,{70.},:)
  ts1_lat!0="time"
  ts1_lat!1="lat"
  ts1_lat!2="lev"
  ts1_lat&time=vmro3&time
  ts1_lat&lat=vmro3&lat
  ts1_lat_avg=dim_avg_n_Wrap(ts1_lat,0)
  lev_label=(/"2 hPa","10 hPa","70 hPa"/)

;regions
  region_label=(/"Global","Trop","NH Midlatitudes","SH Midlatitudes","NH Polar","SH Polar"/)
  region_max=(/90.,10.,60.,-35.,90.,-60./)
  region_min=(/-90.,-10.,35.,-60.,60.,-90./)

  ts_regions=new((/dimsizes(ts_avg&time),dimsizes(plev),dimsizes(region_label)/),typeof(ts_avg))
  do ir=0,dimsizes(region_label)-1
       ts_regions(:,:,ir)=area_operations(vmro31,region_max(ir),region_min(ir),0,360,"average",True)
  end do
  ts_regions!2="region"
  ts_region_avg=new((/3,dimsizes(plev),dimsizes(region_label)/),typeof(ts_avg))
  ts_region_avg(0,:,0:3)=dim_avg_n_Wrap(ts_regions(:,:,0:3),0)
  ts_region_avg(1,:,0:3)=dim_avg_n_Wrap(ts_regions(0:179,:,0:3),0)
  ts_region_avg(2,:,0:3)=dim_avg_n_Wrap(ts_regions(180:,:,0:3),0)
  FMA=new((/33,dimsizes(plev)/),typeof(ts_avg))
  SON=new((/33,dimsizes(plev)/),typeof(ts_avg))
  do il=0,dimsizes(plev)-1
     FMA(:,il) = month_to_season(ts_regions(:,il,4),"FMA")
     SON(:,il)= month_to_season(ts_regions(:,il,5),"SON")
  end do
  ts_region_avg(0,:,4)=dim_avg_n_Wrap(FMA,0)
  ts_region_avg(1,:,4)=dim_avg_n_Wrap(FMA(0:14,:),0)
  ts_region_avg(2,:,4)=dim_avg_n_Wrap(FMA(15:,:),0)
  ts_region_avg(0,:,5)=dim_avg_n_Wrap(SON,0)
  ts_region_avg(1,:,5)=dim_avg_n_Wrap(SON(0:14,:),0)
  ts_region_avg(2,:,5)=dim_avg_n_Wrap(SON(15:,:),0)
  seas=(/"ANN","ANN","ANN","ANN","FMA","SON"/)
  ts_region_avg!0="line"
  ts_region_avg!1="plev"
  ts_region_avg!2="region"
  ts_region_avg&plev=plev
;***********************time serie************************************************
  plot=new(12,"graphic")
  wks1 = gsn_open_wks("ps","ts_region_plev_time"+"CDS-OZONE-L3-LP")
  gsn_define_colormap(wks1,"MPL_RdYlGn")
  gsn_reverse_colormap(wks1)
  k=0
  do ir=0,5
     res2 = True                                    ; plot mods desired
     res2@gsnDraw                  = False            ; don't draw yet
     res2@gsnFrame                 = False            ; don't advance frame yet

     res2@cnFillOn         = True               ; color Fill
     res2@cnFillMode       = "RasterFill"       ; Raster Mode
     res2@cnLinesOn        =  False             ; Turn off contour lines
     res2@cnLineLabelsOn   =  False
     res2@cnLinesOn        =  False
     res2@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     res2@cnMinLevelValF       = 0.5               ; set min contour level
     res2@cnMaxLevelValF       = 10.               ; set max contour level
     res2@cnLevelSpacingF      = 0.5               ; set contour spacing
    ; res2@cnLevelSelectionMode = "ExplicitLevels"
    ; res2@cnLevels = (/0.01,0.1,0.3,0.5,0.7,1.,1.5,2.,3.,4.,5.,6.,7.,8.,9.,10./)

     res2@lbLabelBarOn = False
     res2@tmYRMode             = "Automatic"
     res2@tiMainString         = region_label(ir)    ; plot title
     fmt = "%3.2f"
;  res@gsnLeftString      = "mean: "+sprintf(fmt,avg_avail)
     res2@gsnRightString     = "[ppmv]"
     res2@vpHeightF= 0.4                    ; change aspect ratio of plot
     res2@vpWidthF = 0.6
     res2@cnInfoLabelOn=False
     resTick                = True
     ;resTick@ttmFormat      = "%N/%Y"
     resTick@ttmFormat      = "%Y"
     resTick@ttmAxis        = "XB"
     resTick@ttmMajorStride = 36
     resTick@ttmMinorStride = 12
     time_axis_labels( ts_regions&time,  res2, resTick )

     ts_regions&plev=plev
     plot(k) =gsn_csm_pres_hgt(wks1,ts_regions(plev|:,time|:,region|ir), res2)  ; create plot
     k=k+1
     res21=True
     res21@gsnDraw                  = False            ; don't draw yet
     res21@gsnFrame                 = False            ; don't advance frame yet
     res21@vpHeightF= 0.4                    ; change aspect ratio of plot
     res21@vpWidthF = 0.3
     res21@tiXAxisString = "O3 [ppmv]"                 ; x axis title
     res21@trYMinF      = 1
     res21@trYMaxF      = 250
     res21@trXMinF      = 0.
     res21@trXMaxF      = 10.
     res21@tiMainString         = seas(ir)
     res21@tmYLMode        = "Explicit"
     res21@tmYLValues      = (/250,200,150,100,70,50,30,10/)
     res21@tmYLLabels      = (/"250","200","150","100","70","50","30","10"/)
     res21@tmXBMode        = "Explicit"
     res21@tmXBValues      = (/0.,1.,2.,3.,4.,5.,6.,7.,8.,9./)
     res21@tmXBLabels      = (/"0.","1.","2.","3.","4.","5.","6.","7.","8.","9."/)
     res21@tmXBLabelFontHeightF=0.015
     res21@xyYStyle          = "Log"
     res21@trYReverse =True
     res21@xyLineColors     = (/"NavyBlue"/)
     ;res21@xyMonoLineThickness =False
     res21@xyLineThickness = 3
     plot(k) = gsn_csm_xy(wks1,ts_region_avg(line|0,plev|:,region|ir),ts_region_avg&plev, res21)  ; create plot
     k=k+1
end do

   resP=True
   resP@gsnPanelLabelBar    = True
   resP@gsnPanelXWhiteSpacePercent =0.5
   gsn_panel(wks1,plot,(/3,4/),resP)

;***********************time serie 2************************************************
  delete(res2)
  delete(res21)
  plota=new(12,"graphic")
  wksa = gsn_open_wks("ps","ts2_region_plev_time"+"CDS-OZONE-L3-LP") 
  gsn_define_colormap(wksa,"MPL_RdYlGn")
  gsn_reverse_colormap(wksa)
  k=0
  do ir=0,5
     res2 = True                                    ; plot mods desired
     res2@tiMainString         = region_label(ir)    ; plot title
     fmt = "%3.2f"
     res2@gsnDraw                  = False            ; don't draw yet
     res2@gsnFrame                 = False            ; don't advance frame yet
     res2@vpHeightF= 0.4                    ; change aspect ratio of plot
     res2@vpWidthF = 0.6
     res2@tiYAxisString = "O3 [ppmv]"                 ; x axis title
     res2@tmYLValues      = (/0.,1.,2.,3.,4.,5.,6.,7.,8.,9./)
     res2@tmYLLabels      = (/"0.","1.","2.","3.","4.","5.","6.","7.","8.","9."/)
   
     resTick                = True
     ;resTick@ttmFormat      = "%N/%Y"
     resTick@ttmFormat      = "%Y"
     resTick@ttmAxis        = "XB"
     resTick@ttmMajorStride = 36
     resTick@ttmMinorStride = 12
     time_axis_labels( ts_regions&time,  res2, resTick )

     ts_regions&plev=plev
     res2@xyDashPattern  = 0 
     res2@xyLineColors     = (/"Black","Green","Red"/)
     res2@xyMonoLineThickness =False
     res2@xyLineThickness = 2.     
     time1=ts_regions&time
     res2@trXMinF      = time1(0) 
     res2@trXMaxF      = time1(dimsizes(time1)-1)
     if (ir.eq.0)then
        res2@pmLegendDisplayMode    = "Always"            ; turn on legend
     end if

     if (ir.eq.1)then
        delete(res2@pmLegendDisplayMode)
     end if

     res2@pmLegendSide           = "Top"               ; Change location of
     res2@pmLegendParallelPosF   = 0.7                  ; move units right
     res2@pmLegendOrthogonalPosF = -0.6                ; move units down
     res2@pmLegendWidthF         = 0.2                ; Change width and
     res2@pmLegendHeightF        = 0.18                ; height of legend.
     res2@lgPerimOn              = False               ; turn off box around
     res2@lgLabelFontHeightF     = .02                 ; label font height
     res2@xyExplicitLegendLabels = (/"70 hPa","10 hPa","2 hPa"/)         ;
     sel_lev=(/70.,10.,2./)
     ta_regionsa=ts_regions({plev|sel_lev},time|:,region|ir)
     ta_regionsa(0,:)=rmAnnCycle1D(ta_regionsa(0,:))
     ta_regionsa(1,:)=rmAnnCycle1D(ta_regionsa(1,:))
     ta_regionsa(2,:)=rmAnnCycle1D(ta_regionsa(2,:))

     plota(k) =gsn_csm_xy(wksa,ts_regions&time,ta_regionsa, res2)  ; create plot
     k=k+1
     res21=True
     res21@gsnDraw                  = False            ; don't draw yet
     res21@gsnFrame                 = False            ; don't advance frame yet
     res21@vpHeightF= 0.4                    ; change aspect ratio of plot
     res21@vpWidthF = 0.3
     res21@tiXAxisString = "O3 [ppmv]"                 ; x axis title
     res21@tiYAxisString = "Pressure (hPa)" 
     res21@trYMinF      = 1
     res21@trYMaxF      = 250
     res21@trXMinF      = 0.
     res21@trXMaxF      = 10.
     res21@tiMainString         = seas(ir)+" "+ region_label(ir)
     res21@tmYLMode        = "Explicit"
     res21@tmYLValues      = (/250,200,150,100,70,50,30,10/)
     res21@tmYLLabels      = (/"250","200","150","100","70","50","30","10"/)
     res21@tmXBMode        = "Explicit"
     res21@tmXBValues      = (/0.,1.,2.,3.,4.,5.,6.,7.,8.,9./)
     res21@tmXBLabels      = (/"0.","1.","2.","3.","4.","5.","6.","7.","8.","9."/)
     res21@tmXBLabelFontHeightF=0.015
     res21@xyYStyle          = "Log"
     res21@trYReverse =True
     res21@xyLineColors     = (/"NavyBlue"/)
     res21@xyLineThickness = 3
     plota(k) = gsn_csm_xy(wksa,ts_region_avg(line|0,plev|:,region|ir),ts_region_avg&plev, res21)  ; create plot
     k=k+1
end do

   resP=True
   resP@gsnPanelLabelBar    = True
   resP@gsnPanelXWhiteSpacePercent =0.1 
   gsn_panel(wksa,plota((/0,2,4,6,8,10/)),(/3,2/),resP)
   gsn_panel(wksa,plota((/1,3,5,7,9,11/)),(/3,2/),resP)
;*********èèèè3333333time series++++++++++++

;***********************time serie 2************************************************
delete(res2)
delete(res21)
delete(resTick)
plotb=new(12,"graphic")
  wksb = gsn_open_wks("ps","ts3_region_plev_time"+"CDS-OZONE-L3-LP")
  gsn_define_colormap(wksa,"MPL_RdYlGn")
  gsn_reverse_colormap(wksa)
k=0
do ir=0,5
  res2 = True                                    ; plot mods desired
  res2@tiMainString         = region_label(ir)    ; plot title
  fmt = "%3.2f"
  res2@gsnDraw                  = False            ; don't draw yet
  res2@gsnFrame                 = False            ; don't advance frame yet
  res2@vpHeightF= 0.4                    ; change aspect ratio of plot
  res2@vpWidthF = 0.6
  res2@tiYAxisString = "O3 [ppmv]"                 ; x axis title
  res2@tmYLValues      = (/0.,1.,2.,3.,4.,5.,6.,7.,8.,9./)
  res2@tmYLLabels      = (/"0.","1.","2.","3.","4.","5.","6.","7.","8.","9."/)
  sel_lev=(/2.,10.,70./)
  anni=month_to_annual(ts_regions({plev|2.},time|:,region|ir),1)
  ts_regions_ann=new((/3,33/),"float")
  ts_regions_ann(0,:)=month_to_annual(ts_regions({plev|2.},time|:,region|ir),1)
  ts_regions_ann(1,:)=month_to_annual(ts_regions({plev|10.},time|:,region|ir),1)
  ts_regions_ann(2,:)=month_to_annual(ts_regions({plev|70.},time|:,region|ir),1)
  ts_regions_ann!0="lev"
  ts_regions_ann!1="year"
  ts_regions_ann&year=ispan(1985,2017,1)




   resTick                = True
   resTick@ttmFormat      = "%Y"
   resTick@ttmAxis        = "XB"
   resTick@ttmMajorStride = 3
   resTick@ttmMinorStride = 1

   ts_regions&plev=plev
   res2@xyLineColors     = (/"Red","Green","Black"/)
   res2@xyMonoLineThickness =False
   res2@xyLineThickness = 3
   time1=ts_regions&time
   res2@trXMinF      = 1985
   res2@trXMaxF      = 2017
   if (ir.eq.5)then
     res2@pmLegendDisplayMode    = "Always"            ; turn on legend
   end if
   res2@pmLegendSide           = "Top"               ; Change location of
   res2@pmLegendParallelPosF   = -0.3                  ; move units right
   res2@pmLegendOrthogonalPosF = -0.6                ; move units down
   res2@pmLegendWidthF         = 0.2                ; Change width and
   res2@pmLegendHeightF        = 0.3                ; height of legend.
   res2@lgPerimOn              = False               ; turn off box around
   res2@lgLabelFontHeightF     = .02                 ; label font height
   res2@xyExplicitLegendLabels = (/"2 hPa","10 hPa","70 hPa"/)         ;
   plotb(k) =gsn_csm_xy(wksb,ispan(1985,2017,1),ts_regions_ann, res2)  ; create plot
   k=k+1
   res21=True
   res21@gsnDraw                  = False            ; don't draw yet
   res21@gsnFrame                 = False            ; don't advance frame yet
   res21@vpHeightF= 0.4                    ; change aspect ratio of plot
   res21@vpWidthF = 0.3
   res21@tiXAxisString = "O3 [ppmv]"                 ; x axis title
   res21@trYMinF      = 1
   res21@trYMaxF      = 250
   res21@trXMinF      = 0.
   res21@trXMaxF      = 10.
   res21@tiMainString         = seas(ir)
   res21@tmYLMode        = "Explicit"
   res21@tmYLValues      = (/250,200,150,100,70,50,30,10/)
   res21@tmYLLabels      = (/"250","200","150","100","70","50","30","10"/)
   res21@tmXBMode        = "Explicit"
   res21@tmXBValues      = (/0.,1.,2.,3.,4.,5.,6.,7.,8.,9./)
   res21@tmXBLabels      = (/"0.","1.","2.","3.","4.","5.","6.","7.","8.","9."/)
   res21@tmXBLabelFontHeightF=0.015
   res21@xyYStyle          = "Log"
   res21@trYReverse =True
   res21@xyLineColors     = (/"NavyBlue"/)
   res21@xyLineThickness = 3
   plotb(k) = gsn_csm_xy(wksb,ts_region_avg(line|0,plev|:,region|ir),ts_region_avg&plev, res21)  ; create plot
   k=k+1
end do

  resP=True
  resP@gsnPanelLabelBar    = True
  resP@gsnPanelXWhiteSpacePercent =0.5
  gsn_panel(wksb,plotb,(/3,4/),resP)
;************************************************
  plot3=new(6,"graphic")
  wks3 = gsn_open_wks("png","zon_"+"CDS-OZONE-L3-LP")       

  gsn_define_colormap(wks3,"MPL_RdYlGn")
  gsn_reverse_colormap(wks3)
  res = True                                    ; plot mods desired
  res@gsnDraw                  = False            ; don't draw yet
  res@gsnFrame                 = False            ; don't advance frame yet
  res@lbLabelBarOn = False
  res@cnLineLabelsOn   =  False
  res@cnLinesOn        =  False 
  res@cnFillOn         = True               ; color Fill
  res@cnFillMode       = "RasterFill"       ; Raster Mode
  res@cnLinesOn        =  False             ; Turn off contour lines
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = 0.5               ; set min contour level
  res@cnMaxLevelValF       = 10.               ; set max contour level
  res@cnLevelSpacingF      = 0.5               ; set contour spacing
 ; res@cnLevelSelectionMode = "ExplicitLevels"
  ;res@cnLevels = (/0.01,0.1,0.3,0.5,0.7,1.,1.5,2.,3.,4.,5.,6.,7.,8.,9.,10./)


  fmt = "%3.2f"
  res@vpHeightF= 0.4                    ; change aspect ratio of plot
  res@vpWidthF = 0.7
  res@tmYLMode        = "Explicit"
  res@tmYLValues      = (/-60,-30,0,30,60/)
  res@tmYLLabels      = (/"60S","30S","Eq","30N","60N"/)
  res@tmYLLabelFontHeightF=0.015
  res@tmXBLabelFontHeightF=0.015
  res@cnInfoLabelOn=False
  resTick                = True
  ;resTick@ttmFormat      = "%N/%Y"
  resTick@ttmFormat      = "%Y"
  resTick@ttmAxis        = "XB"
  resTick@ttmMajorStride = 36
  resTick@ttmMinorStride = 12
  res@gsnRightString     = "[ppmv]"
; Set resources necessary to customize Y axis labels
  time_axis_labels( ts1_lat&time,  res, resTick )
  res@tiMainString         = lev_label(0) 
  plot3(0) = gsn_csm_lat_time(wks3,ts1_lat(lat|:,time|:,lev|0), res)  ; create plot

  res@tiMainString         = lev_label(1) 
  plot3(2) = gsn_csm_lat_time(wks3,ts1_lat(lat|:,time|:,lev|1), res)  ; create plot
  res@tiMainString         = lev_label(2) 
  plot3(4) = gsn_csm_lat_time(wks3,ts1_lat(lat|:,time|:,lev|2), res)  ; create plot
  delete(res21)
  res21=True
  res21@gsnDraw                  = False            ; don't draw yet
  res21@gsnFrame                 = False            ; don't advance frame yet
 ; res21@tiXAxisString = "Tropospheric Ozone Column [DU]"           ; y axis title
  res21@vpHeightF= 0.4                    ; change aspect ratio of plot
  res21@vpWidthF = 0.3
  ;res@tiXAxisString = "Time"                 ; x axis title
  res21@trYMinF      = -90.
  res21@trYMaxF      = 90.
  res21@trXMinF      = 0.
  res21@trXMaxF      = 10.
  res21@tmYLMode        = "Explicit"
  res21@tmYLValues      = (/-60,-30,0,30,60/)
  res21@tmYLLabels      = (/"60S","30S","Eq","30N","60N"/)
  res21@tmYLLabelFontHeightF=0.015
  res21@tmXBMode        = "Explicit"
  res21@tmXBValues      = (/0.,2.,4.,6.,8.,10./)
  res21@tmXBLabels      = (/"0.","2.","4.","6.","8.","10."/)
  res21@tmXBLabelFontHeightF=0.015

;  res1@gsnAttachPlotsXAxis=True
  res21@xyLineColor      = "NavyBlue"
  res21@xyLineThicknessF  =3.5
  res21@tiXAxisString = "O3 [ppm]"
  res21@tiMainString = "lev_label(0)"
  plot3(1) = gsn_csm_xy(wks3,ts1_lat_avg(lat|:,lev|0),ts1_lat_avg&lat, res21)  ; create plot
  res21@tiMainString = "lev_label(1)"
  plot3(3) = gsn_csm_xy(wks3,ts1_lat_avg(lat|:,lev|1),ts1_lat_avg&lat, res21)
  res21@tiMainString = "lev_label(2)" 
  plot3(5) = gsn_csm_xy(wks3,ts1_lat_avg(lat|:,lev|2),ts1_lat_avg&lat, res21)
  resP=True
  resP@gsnPanelLabelBar    = True
  resP@gsnPanelXWhiteSpacePercent =0.5
  gsn_panel(wks3,plot3,(/3,2/),resP)

;************************************************
  delete(res)
  wks4 = gsn_open_wks("png","lat_lev"+"CDS-OZONE-L3-LP")                ; send graphics to PNG file
; Set resources necessary to customize Y axis labels
  media=dim_avg_n_Wrap(vmro3,0)

  gsn_define_colormap(wks4,"MPL_RdYlGn")
  gsn_reverse_colormap(wks4)
  res = True                                    ; plot mods desired
  ;res@gsnDraw                  = False            ; don't draw yet
  ;res@gsnFrame                 = False            ; don't advance frame yet
  ;res@lbLabelBarOn = False
   res@cnLineLabelsOn   =  False
  res@cnLinesOn        =  False
 res@cnFillOn         = True               ; color Fill
  res@cnFillMode       = "RasterFill"       ; Raster Mode
  res@cnLinesOn        =  False             ; Turn off contour lines
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = 0.5               ; set min contour level
  res@cnMaxLevelValF       = 10.               ; set max contour level
  res@cnLevelSpacingF      = 0.5               ; set contour spacing
  ;res@cnLevelSelectionMode = "ExplicitLevels"
  ;res@cnLevels = (/0.01,0.1,0.3,0.5,0.7,1.,1.5,2.,3.,4.,5.,6.,7.,8.,9.,10./)

  fmt = "%3.2f"
  res@vpHeightF= 0.5                    ; change aspect ratio of plot
  res@vpWidthF = 0.5
   res@gsnRightString     = "[ppmv]"
  res@tmXBMode        = "Explicit"
  res@tmXBValues      = (/-60,-30,0,30,60/)
  res@tmXBLabels      = (/"60S","30S","Eq","30N","60N"/)
  res@tmXBLabelFontHeightF=0.015
  res@tmXBLabelFontHeightF=0.015
  res@cnInfoLabelOn=False
 
  p5=gsn_csm_pres_hgt(wks4,media(plev|:,lat|:), res)  ; create plot

;****************************************************
;************************************************
  delete(res)
  plot6=new(12,"graphic")
  month_label=(/"Jannuary","February","March","April","May","June","July","August","September","October","November","December"/)
  wks5 = gsn_open_wks("png","seasonal_lat_lev"+"CDS-OZONE-L3-LP")                ; send graphics to PNG file
; Set resources necessary to customize Y axis labels
 media_seas=clmMonTLL(vmro3)

  gsn_define_colormap(wks5,"MPL_RdYlGn")
  gsn_reverse_colormap(wks5)
  do im=0,11
    res = True                                    ; plot mods desired
    res@gsnDraw                  = False            ; don't draw yet
    res@gsnFrame                 = False            ; don't advance frame yet
    res@lbLabelBarOn = False
    res@cnLineLabelsOn   =  False
    res@cnLinesOn        =  False
    res@cnFillOn         = True               ; color Fill
    res@cnFillMode       = "RasterFill"       ; Raster Mode
    res@cnLinesOn        =  False             ; Turn off contour lines
    res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
    res@cnMinLevelValF       = 0.5               ; set min contour level
    res@cnMaxLevelValF       = 10.               ; set max contour level
    res@cnLevelSpacingF      = 0.5               ; set contour spacing
  ;  res@cnLevelSelectionMode = "ExplicitLevels"
  ;  res@cnLevels = (/0.01,0.1,0.3,0.5,0.7,1.,1.5,2.,3.,4.,5.,6.,7.,8.,9.,10./) 
    res@gsnRightString     = "[ppmv]"
    fmt = "%3.2f"
    res@vpHeightF= 0.5                    ; change aspect ratio of plot
    res@vpWidthF = 0.5
    res@tmXBMode        = "Explicit"
    res@tmXBValues      = (/-60,-30,0,30,60/)
    res@tmXBLabels      = (/"60S","30S","Eq","30N","60N"/)
    res@tmXBLabelFontHeightF=0.015
    res@tmXBLabelFontHeightF=0.015
    res@cnInfoLabelOn=False
    res@tiMainString         = month_label(im)
     plot6(im)=gsn_csm_pres_hgt(wks5,media_seas(plev|:,lat|:,month|im), res)  ; create plot
  end do
  resP=True
  resP@gsnPanelLabelBar    = True
  resP@gsnPanelXWhiteSpacePercent =0.5
  gsn_panel(wks5,plot6,(/4,3/),resP)
;*******************************trend***********************************
 printVarSummary(vmro3)
  dim=dimsizes(vmro3)
  print(dim)
 
 opt = False
    avail1=new((/dim(1),dim(2)/),typeof(vmro3))
     do ilev=0,dim(1)-1
      do ilat=0,dim(2)-1
          avail1(ilev,ilat)=num(.not.ismissing(vmro3(:,ilev,ilat)))/tofloat(dim(0))
      end do
    end do

  vmro32=rmMonAnnCycTLL(vmro3)   ;remove annual cycle
  trend1  = trend_manken(vmro32, opt, 0)  ; ===> px(2,nlat,mlon)
  trend1(1,:,:)=(/trend1(1,:,:)*12*10/)
  trend1(0,:,:)=where(avail1.gt.0.5,trend1(0,:,:),trend1@_FillValue)
  trend1(1,:,:)=where(avail1.gt.0.5,trend1(1,:,:),trend1@_FillValue)                                   ; if meta data is desired
  copy_VarCoords(vmro32(0,:,:),trend1(0,:,:))
  trend1!0= "prob_trend"               ; ===> size 2
  wkst = gsn_open_wks("png","trend_1_"+"CDS-OZONE-L3-LP")         ; send graphics to PNG file
  gsn_define_colormap(wkst,"NEO_div_vegetation_c")
  plott=new(2,"graphic")
  delete(res)
  res = True                                    ; plot mods desired
  res@gsnDraw                  = False            ; don't draw yet
  res@gsnFrame                 = False            ; don't advance frame yet

  res@cnFillOn         = True               ; color Fill
  res@cnFillMode       = "RasterFill"       ; Raster Mode
  res@cnLinesOn        =  False             ; Turn off contour lines
  res@cnLineLabelsOn       =False
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = -0.5               ; set min contour level
  res@cnMaxLevelValF       = 0.5              ; set max contour level
  res@cnLevelSpacingF      = 0.02               ; set contour spacing
  res@tmYRMode             = "Automatic"
  ;res@cnLabelBarEndStyle="ExcludeOuterBoxes"
   res@tiMainString         = "Trend 1985-2017"    ; plot title
   fmt = "%3.2f"
   res@gsnLeftString      = "[ppmv/decade]"
   ;res@gsnRightString     = "mean: "+sprintf(fmt,avg_clim1)+" std:"+sprintf(fmt,avg_std1)           ; plot center string
   res@lbOrientation = "vertical"
   ;res@lbLabelBarOn = False
   plott(0) = gsn_csm_pres_hgt(wkst,trend1(prob_trend|1,plev|:,lat|:), res)  ; create plot
   res@cnMinLevelValF       = 0               ; set min contour level
   res@cnMaxLevelValF       = 1               ; set max contour level
   res@cnLevelSpacingF      = 0.05               ; set contour spacing
   res@tiMainString         = "Trend Significance"
   res@gsnLeftString      = "[1]"
   plott(1) = gsn_csm_pres_hgt(wkst,trend1(prob_trend|0,plev|:,lat|:), res)
   delete(resP)
    resP=True
  ;resP@gsnPanelLabelBar    = True                ; add common colorbar
   gsn_panel(wkst,plott,(/2,1/),resP)
;*******************************trend 2***********************************
   plott_red=new(2,"graphic")

   red_vmro3=vmro3(180:,:,:)
   dim=dimsizes(red_vmro3)

   opt = False
   avail1=new((/dim(1),dim(2)/),typeof(vmro3))
   do ilev=0,dim(1)-1
      do ilat=0,dim(2)-1
          avail1(ilev,ilat)=num(.not.ismissing(red_vmro3(:,ilev,ilat)))/tofloat(dim(0))
      end do
   end do

  red_vmro32=rmMonAnnCycTLL(red_vmro3)   ;remove annual cycle
  red_trend1  = trend_manken(red_vmro32, opt, 0)  ; ===> px(2,nlat,mlon)
  red_trend1(1,:,:)=(/red_trend1(1,:,:)*12*10/)
  red_trend1(0,:,:)=where(avail1.gt.0.5,red_trend1(0,:,:),red_trend1@_FillValue)
  red_trend1(1,:,:)=where(avail1.gt.0.5,red_trend1(1,:,:),red_trend1@_FillValue)                                   ; if meta data is desired
  copy_VarCoords(red_vmro32(0,:,:),red_trend1(0,:,:))
  
  red_trend1&plev=plev
  red_trend1!0= "prob_trend"               ; ===> size 2
  wkst1 = gsn_open_wks("png","trend_2000_2017_"+"CDS-OZONE-L3-LP")         ; send graphics to PNG file
  gsn_define_colormap(wkst1,"NEO_div_vegetation_c")
  ;plott_red=new(2,"graphic")
  delete(res)
  res = True                                    ; plot mods desired
  res@gsnDraw                  = False            ; don't draw yet
  res@gsnFrame                 = False            ; don't advance frame yet

  res@cnFillOn         = True               ; color Fill
  res@cnFillMode       = "RasterFill"       ; Raster Mode
  res@cnLinesOn        =  False             ; Turn off contour lines
  res@cnLineLabelsOn       =False
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = -0.3               ; set min contour level
  res@cnMaxLevelValF       = 0.3              ; set max contour level
  res@cnLevelSpacingF      = 0.02               ; set contour spacing
  res@tmYRMode             = "Automatic"
  ;res@cnLabelBarEndStyle="ExcludeOuterBoxes"
  res@tiMainString         = "Trend 2000-2017 "    ; plot title
  fmt = "%3.2f"
  res@gsnLeftString      = "[ppmv/decade]"
  ;res@gsnRightString     = "mean: "+sprintf(fmt,avg_clim1)+" std:"+sprintf(fmt,avg_std1)           ; plot center string
  res@lbOrientation = "vertical"
  ;res@lbLabelBarOn = False
  plott_red(0) = gsn_csm_pres_hgt(wkst1,red_trend1(prob_trend|1,plev|:,lat|:), res)  ; create plot
  res@cnMinLevelValF       = 0               ; set min contour level
  res@cnMaxLevelValF       = 1               ; set max contour level
  res@cnLevelSpacingF      = 0.05               ; set contour spacing
  res@tiMainString         = "Trend Significance"
  res@gsnLeftString      = "[1]"
  plott_red(1) = gsn_csm_pres_hgt(wkst1,red_trend1(prob_trend|0,plev|:,lat|:), res)
  delete(resP)
  resP=True
  ;resP@gsnPanelLabelBar    = True                ; add common colorbar
  gsn_panel(wkst1,plott_red,(/2,1/),resP)



end
