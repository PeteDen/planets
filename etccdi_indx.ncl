;---------------------------------------------------------------
; - This script was written by Xu (Pete) Deng in Mar 2020       
; - University of New South Wales, Canberra                     
;---------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

    model_n = ""
    model_ns = ""
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    rip_n = ""
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    fvar_n0 = "sftlf_fx_"
    fvar_n = "tasmax_day_"
    fvar_n_min = "tasmin_day_"

    var1_arr = (/"txx_y","txn_y","tnx_y","tnn_y","dtr_y","fd_y","tr_y","su_y"/)

;======================================================
    y_period = "1950_2014"
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    latS = -45.5
    ;-39.5
    latN = -9.5
    lonW = 109.5
    lonE = 155.5

    n_lon = 47
    n_lat = 37
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    s_year = 1950
    e_year = 2014
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    n_year = e_year-s_year+1

    s_byear = 1961
    e_byear = 1990
    n_byear = e_byear-s_byear+1


    year_arr1 = ispan(s_year,e_year,1)
    year_arr1@calendar = "365"
    ;year_arr1@calendar = "365"

    considerFebL = "considered"
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    byear_arr1 = ispan(s_byear,e_byear,1)
    ;byear_arr1@calendar = "standard"
    byear_arr1@calendar = year_arr1@calendar

    k_temp = 273.15

;===================================================== 
    indir = ""
    indir_min = ""
    indir1 = ""
    
    outdir = ""
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;=====================================================
    f1_t = addfile(indir+fvar_n+model_n+"_historical"+rip_n+".nc", "r")
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    t1_t = f1_t->tasmax(:,{latS:latN},{lonW:lonE})

    dims1_t = dimsizes(t1_t)
    print(dims1_t)
    ntime1_t = dims1_t(0)
    nlat1_t = dims1_t(1)
    nlon1_t = dims1_t(2)

;==================data processing=================
    files1 = systemfunc("ls "+indir+fvar_n+model_n+"_historical"+rip_n+"*"+".nc")
    print(files1)

    files1_min = systemfunc("ls "+indir_min+fvar_n_min+model_n+"_historical"+rip_n+"*"+".nc")
    print(files1_min)

    f0 = addfile(indir1+fvar_n0+model_n+"_historical"+rip_n+".nc","r")
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    sftlf1 = f0->sftlf({latS:latN},{lonW:lonE})
    ;land_area_fraction
    printVarSummary(sftlf1)


    n_lyear = num(isleapyear(year_arr1))
    print("=================number of leap years: "+n_lyear+"=================")

    n_days = 366*n_lyear+365*(n_year-n_lyear)

    t1_maxt = new((/n_days,nlat1_t,nlon1_t/),typeof(t1_t),t1_t@_FillValue)
    printVarSummary(t1_maxt)

    t1_mint = new((/n_days,nlat1_t,nlon1_t/),typeof(t1_t),t1_t@_FillValue)
    printVarSummary(t1_mint)


    file_s = 5
    ;;;;;;;;;;;;;;;;;;;;;;;;
    file_e = dimsizes(files1)-1
    s_indxt = 0

    do i = file_s, file_e
        ;print(files1(i))
        ;print(files1_min(i))

        f1 = addfile(files1(i),"r")
        f1_min = addfile(files1_min(i),"r")

        time1 = f1->time

        dims1 = dimsizes(time1)
        ;print(dims1)
        ntime1 = dims1(0)

        ut1 = cd_calendar(time1, 0)
        ;print(ut1(n_days-1,:))

        if (i.eq.file_s) then 
            s_indx = min(ind(ut1(:,0).eq.s_year))
            e_indx = ntime1-1
        else if (i.eq.file_e) then
            s_indx = 0
            e_indx = max(ind(ut1(:,0).eq.e_year))
        else
            s_indx = 0
            e_indx = ntime1-1
        end if 
        end if 

        t1 = f1->tasmax(s_indx:e_indx,{latS:latN},{lonW:lonE})
        t1_min = f1_min->tasmin(s_indx:e_indx,{latS:latN},{lonW:lonE})
        ;time_1 = f1_min->time

        e_indxt = s_indxt+(e_indx-s_indx)

        t1_maxt(s_indxt:e_indxt,:,:) = t1
        t1_mint(s_indxt:e_indxt,:,:) = t1_min

        s_indxt = e_indxt+1

        delete(time1)
        delete(dims1)
        delete(ntime1)
        delete(ut1)
        delete(t1)
        delete(t1_min)
    end do 

    t1 = t1_maxt
    t1_min = t1_mint

    printVarSummary(t1)
    printVarSummary(t1_min)

    ut1 = cd_calendar(t1&time, 0)

;================regriding=================
    klon = fspan(lonW,lonE,n_lon)
    klat = fspan(latS,latN,n_lat)

    t1_aus= linint2_Wrap(t1&lon,t1&lat,t1,False,klon,klat,0) ;(ntim,mlat,mlon)
    t1_min_aus= linint2_Wrap(t1_min&lon,t1_min&lat,t1_min,False,klon,klat,0) ;(ntim,mlat,mlon)

    lmask1= linint2_Wrap(sftlf1&lon,sftlf1&lat,sftlf1,False,klon,klat,0) ;(ntim,mlat,mlon)
    printVarSummary(t1_aus)

;================masking the ocean================
    t1_ausr = t1_aus
    t1_ausr = mask(t1_aus, (lmask1.ge.75.),True)
    
    t1_min_ausr = t1_min_aus
    t1_min_ausr = mask(t1_min_aus, (lmask1.ge.75.),True)

    t1_ausr = t1_ausr - k_temp
    t1_min_ausr = t1_min_ausr - k_temp
    dtr1_ausr = t1_ausr - t1_min_ausr

;===================================================
    dims1r = dimsizes(t1_ausr)
    print(dims1r)
    nt1r = dims1r(0)
    nlat1r = dims1r(1)
    nlon1r = dims1r(2)

;=========txx, txn, tnx, tnn, dtr, su, tr and fd=========
    txx1_y = new((/n_year,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)
    txn1_y = new((/n_year,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)
    tnx1_y = new((/n_year,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)
    tnn1_y = new((/n_year,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)
    dtr1_y = new((/n_year,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)

    fd1_y = new((/n_year,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)
    tr1_y = new((/n_year,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)
    su1_y = new((/n_year,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)

    do ii = 0, n_year-1
        s_y1 = min(ind(ut1(:,0).eq.s_year+ii))
        e_y1 = max(ind(ut1(:,0).eq.s_year+ii))

        txx1_y(ii,:,:) = dim_max_n_Wrap(t1_ausr(s_y1:e_y1,:,:),0)
        txn1_y(ii,:,:) = dim_min_n_Wrap(t1_ausr(s_y1:e_y1,:,:),0)
        tnx1_y(ii,:,:) = dim_max_n_Wrap(t1_min_ausr(s_y1:e_y1,:,:),0)
        tnn1_y(ii,:,:) = dim_min_n_Wrap(t1_min_ausr(s_y1:e_y1,:,:),0)

        dtr1_y(ii,:,:) = dim_avg_n_Wrap(dtr1_ausr(s_y1:e_y1,:,:),0)

        fd1_y(ii,:,:) = dim_num_n(t1_min_ausr(s_y1:e_y1,:,:).lt.0., 0)
        tr1_y(ii,:,:) = dim_num_n(t1_min_ausr(s_y1:e_y1,:,:).gt.20., 0)
        su1_y(ii,:,:) = dim_num_n(t1_ausr(s_y1:e_y1,:,:).gt.25., 0)
    end do 
    ;printVarSummary(fd1_y)
    ;printMinMax(fd1_y, True)

    fd1_y = mask(fd1_y, (lmask1.ge.75.),True)
    tr1_y = mask(tr1_y, (lmask1.ge.75.),True)
    id1_y = mask(id1_y, (lmask1.ge.75.),True)
    su1_y = mask(su1_y, (lmask1.ge.75.),True)

;==============generating ascii files===============
    asciiwrite(outdir+model_ns+"_"+"aus"+y_period+"_"+var1_arr(0)+".txt",txx1_y)
    asciiwrite(outdir+model_ns+"_"+"aus"+y_period+"_"+var1_arr(1)+".txt",txn1_y)
    asciiwrite(outdir+model_ns+"_"+"aus"+y_period+"_"+var1_arr(2)+".txt",tnx1_y)
    asciiwrite(outdir+model_ns+"_"+"aus"+y_period+"_"+var1_arr(3)+".txt",tnn1_y)
    asciiwrite(outdir+model_ns+"_"+"aus"+y_period+"_"+var1_arr(4)+".txt",dtr1_y)

    asciiwrite(outdir+model_ns+"_"+"aus"+y_period+"_"+var1_arr(5)+".txt",fd1_y)
    asciiwrite(outdir+model_ns+"_"+"aus"+y_period+"_"+var1_arr(6)+".txt",tr1_y)
    asciiwrite(outdir+model_ns+"_"+"aus"+y_period+"_"+var1_arr(7)+".txt",su1_y)    

    print("=================end=================")
end 

