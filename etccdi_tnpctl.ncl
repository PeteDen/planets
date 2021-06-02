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

;================================================
    var_tn10p = "tn10p_y"
    var_tn90p = "tn90p_y"
    var_csdi = "csdi_y"

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
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    n_year = e_year-s_year+1

    s_byear = 1961
    e_byear = 1990
    n_byear = e_byear-s_byear+1
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


    year_arr1 = ispan(s_year,e_year,1)
    ;year_arr1@calendar = "standard"
    year_arr1@calendar = "365"

    considerFebL = "considered"
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    byear_arr1 = ispan(s_byear,e_byear,1)
    ;byear_arr1@calendar = "standard"
    byear_arr1@calendar = year_arr1@calendar

    k_temp = 273.15

    n_win = 5*n_byear
    pth_per1 = 10.
    pth_per2 = 90.
    
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
    print(ut1(n_days-1,:))

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

;===================================================
    dims1r = dimsizes(t1_ausr)
    print(dims1r)
    nt1r = dims1r(0)
    nlat1r = dims1r(1)
    nlon1r = dims1r(2)

;=======================contruct the 5-day window in out-base period=======================
    win_outb1 = new((/365,n_win,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)
    printVarSummary(win_outb1)

    do ii = 0, 365-1
        do jj=0, n_byear-1 

        indx_b = min(ind(ut1(:,0).eq.byear_arr1(jj)))
        s_ind = 5*jj
        e_ind = 5*jj+4

        if ((isleapyear(byear_arr1(jj))).and.(ii.gt.58)) then

            win_outb1(ii,s_ind:e_ind,:,:) = t1_min_ausr(indx_b-1+ii:indx_b+3+ii,:,:)
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        else
            win_outb1(ii,s_ind:e_ind,:,:) = t1_min_ausr(indx_b-2+ii:indx_b+2+ii,:,:)

        end if
        end do 
    end do 

    printVarSummary(win_outb1)

;=======================calculate the 10th and 90th percentile=======================
    indx_outb1  = dim_pqsort_n(win_outb1, 2, 1)

    ppp1 = pth_per1/100.
    ppp2 = pth_per2/100.

    jjj1 = toint(floor(ppp1*(n_win+1)))
    jjj2 = toint(floor(ppp2*(n_win+1)))

    fff1 = ppp1*(n_win+1)-jjj1
    fff2 = ppp2*(n_win+1)-jjj2

    print(ppp1+"=========="+jjj1+"=========="+fff1+"==========")
    print(ppp2+"=========="+jjj2+"=========="+fff2+"==========")


    yj1 = win_outb1(:,jjj1-1,:,:)

    yjp1 = win_outb1(:,jjj1,:,:)

    per_out1 = (1-fff1)*yj1 + fff1*yjp1


    yj2 = win_outb1(:,jjj2-1,:,:)

    yjp2 = win_outb1(:,jjj2,:,:)

    per_out2 = (1-fff2)*yj2 + fff2*yjp2


    delete(yj1)
    delete(yjp1)
    delete(yj2)
    delete(yjp2)
    delete(win_outb1)
    printVarSummary(per_out1)

;=================for the day Feb 29th in out-base period===============
    if ((byear_arr1@calendar.eq."standard").and.(considerFebL.eq."considered")) then

        n_lyear = num(isleapyear(byear_arr1))
    
        win_outb1_leap = new((/n_win,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)
        ;printVarSummary(win_outb1_leap)
    
        n_win_leap = n_lyear*5+(n_byear-n_lyear)*4
    
        do ii=59,59
        do jj=0, n_byear-1
    
            indx_b = min(ind(ut1(:,0).eq.s_byear+jj)) 
            ;indx_b = min(ind(ut1(:,0).eq.byear_arr1(jj)))
            s_ind = 5*jj
            e_ind = 5*jj+4
    
            if ((isleapyear(byear_arr1(jj)))) then
                win_outb1_leap(s_ind:e_ind,:,:) = t1_min_ausr(indx_b-2+ii:indx_b+2+ii,:,:)
                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            else
                win_outb1_leap(s_ind:s_ind+1,:,:) = t1_min_ausr(indx_b-2+ii:indx_b-1+ii,:,:)
                win_outb1_leap(s_ind+2:s_ind+2,:,:) = t1_t@_FillValue
                win_outb1_leap(s_ind+3:e_ind,:,:) = t1_min_ausr(indx_b+ii:indx_b+1+ii,:,:)
            end if
        end do
        end do 
    
        indx_outb1_leap  = dim_pqsort_n(win_outb1_leap, 2, 0)
    
        jjj1_leap = toint(floor(ppp1*(n_win_leap+1)))
    
        fff1_leap = ppp1*(n_win_leap+1)-jjj1_leap

        yj1_leap = win_outb1_leap(jjj1_leap-1,:,:)
        yjp1_leap = win_outb1_leap(jjj1_leap,:,:)

        per_out1_leap = (1-fff1_leap)*yj1_leap + fff1_leap*yjp1_leap


        jjj2_leap = toint(floor(ppp2*(n_win_leap+1)))
    
        fff2_leap = ppp2*(n_win_leap+1)-jjj2_leap
        
        yj2_leap = win_outb1_leap(jjj2_leap-1,:,:)
        yjp2_leap = win_outb1_leap(jjj2_leap,:,:)

        per_out2_leap = (1-fff2_leap)*yj2_leap + fff2_leap*yjp2_leap


        delete(yj1_leap)
        delete(yjp1_leap)
        delete(yj2_leap)
        delete(yjp2_leap)
        delete(win_outb1_leap)
        printVarSummary(per_out1_leap)

    end if

;=======================get initial results for out-base period=======================
    t1_min_ausr1 = t1_min_ausr
    t1_min_ausr2 = t1_min_ausr

    do ii = 0, n_year-1
        s_y1 = min(ind(ut1(:,0).eq.s_year+ii))
        e_y1 = max(ind(ut1(:,0).eq.s_year+ii))

        if ((isleapyear(year_arr1(ii)))) then
            t1_min_ausr1(s_y1:s_y1+58,:,:) = t1_min_ausr1(s_y1:s_y1+58,:,:) - per_out1(0:58,:,:)
            t1_min_ausr1(s_y1+59,:,:) = t1_min_ausr1(s_y1+59,:,:) - per_out1_leap(:,:)
            t1_min_ausr1(s_y1+60:e_y1,:,:) = t1_min_ausr1(s_y1+60:e_y1,:,:) - per_out1(59:,:,:)
            

            t1_min_ausr2(s_y1:s_y1+58,:,:) = t1_min_ausr2(s_y1:s_y1+58,:,:) - per_out2(0:58,:,:)
            t1_min_ausr2(s_y1+59,:,:) = t1_min_ausr2(s_y1+59,:,:) - per_out2_leap(:,:)
            t1_min_ausr2(s_y1+60:e_y1,:,:) = t1_min_ausr2(s_y1+60:e_y1,:,:) - per_out2(59:,:,:)
            
        else
            t1_min_ausr1(s_y1:e_y1,:,:) = t1_min_ausr1(s_y1:e_y1,:,:) - per_out1(:,:,:)
            
            t1_min_ausr2(s_y1:e_y1,:,:) = t1_min_ausr2(s_y1:e_y1,:,:) - per_out2(:,:,:)
        end if
    end do

    delete(per_out1)
    delete(per_out2)
    if (byear_arr1@calendar.eq."standard") then
        delete(per_out1_leap)
        delete(per_out2_leap)
    end if
    printVarSummary(t1_min_ausr1)


    tn10p1_o = where(t1_min_ausr1.lt.0, 1, 0) 

    tn90p1_o = where(t1_min_ausr2.gt.0, 1, 0) 

    printVarSummary(tn10p1_o)
    printMinMax(tn10p1_o, True)
    printVarSummary(tn90p1_o)
    printMinMax(tn90p1_o, True)
    delete(t1_min_ausr1)
    delete(t1_min_ausr2)


    tn10p1_y = new((/n_year,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)
    tn90p1_y = new((/n_year,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)
    csdi1_y = new((/n_year,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)

    tn10p1_y_per = new((/n_year,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)
    tn90p1_y_per = new((/n_year,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)

    do ii = 0, n_year-1
        s_y1 = min(ind(ut1(:,0).eq.s_year+ii))
        e_y1 = max(ind(ut1(:,0).eq.s_year+ii))

        tn10p1_arr = dim_numrun_n(tn10p1_o(s_y1:e_y1,:,:), 0, 0)
        tn10p1_y(ii,:,:) = tn10p1_arr(0,:,:)
        aaa = tofloat(e_y1-s_y1+1)
        tn10p1_y_per(ii,:,:) = (tn10p1_arr(0,:,:)/aaa)*100.

        tn90p1_arr = dim_numrun_n(tn90p1_o(s_y1:e_y1,:,:), 0, 0)
        tn90p1_y(ii,:,:) = tn90p1_arr(0,:,:)
        tn90p1_y_per(ii,:,:) = (tn90p1_arr(0,:,:)/aaa)*100.

        csdi1_arr = dim_numrun_n(tn10p1_o(s_y1:e_y1,:,:), 1, 0)
        wgt = fspan(6.,(e_y1-s_y1+1)*1.,(e_y1-s_y1+1-5))
        csdi1_y(ii,:,:) = dim_sum_wgt_n_Wrap(csdi1_arr(5:,:,:), wgt, 1, 0)


        delete(tn10p1_arr)
        delete(tn90p1_arr)
        delete(csdi1_arr)
        delete(wgt)
    end do

    delete(tn90p1_o)
    ;delete(tn10p1_o)
    printVarSummary(tn10p1_y)
    printVarSummary(csdi1_y)

;===========================================
    csdi1_y1 = csdi1_y
    csdi1_y2 = csdi1_y

    tn10p1_y1 = tn10p1_y
    tn90p1_y1 = tn90p1_y
    tn10p1_y_per1 = tn10p1_y_per
    tn90p1_y_per1 = tn90p1_y_per

    csdi1_y_n = csdi1_y2
    csdi1_y_n = mask(csdi1_y2, (lmask1.ge.75.),True)

    tn10p1_y_n = tn10p1_y1
    tn10p1_y_n = mask(tn10p1_y1, (lmask1.ge.75.),True)

    tn90p1_y_n = tn90p1_y1
    tn90p1_y_n = mask(tn90p1_y1, (lmask1.ge.75.),True)

    tn10p1_y_pern = tn10p1_y_per1
    tn10p1_y_pern = mask(tn10p1_y_per1, (lmask1.ge.75.),True)

    tn90p1_y_pern = tn90p1_y_per1
    tn90p1_y_pern = mask(tn90p1_y_per1, (lmask1.ge.75.),True)

    delete(csdi1_y2)
    delete(csdi1_y_n)
    delete(tn10p1_y1)
    delete(tn10p1_y_n)
    delete(tn90p1_y1)
    delete(tn90p1_y_n)
    delete(tn10p1_y_per1)
    delete(tn10p1_y_pern)
    delete(tn90p1_y_per1)
    delete(tn90p1_y_pern)


    do kk = 0, nlon1r-1
        do jj = 0, nlat1r-1
            do ii = 0, n_year-1-1

                e_y1 = max(ind(ut1(:,0).eq.s_year+ii))

                if (.not.(all(tn10p1_o(e_y1-5:e_y1,jj,kk).eq.1.))) then 

                    if (all(tn10p1_o(e_y1-4:e_y1+1,jj,kk).eq.1.)) then 
                        csdi1_y1(ii,jj,kk) = 5 + csdi1_y1(ii,jj,kk)

                    else if (all(tn10p1_o(e_y1-3:e_y1+2,jj,kk).eq.1.)) then 
                        csdi1_y1(ii,jj,kk) = 4 + csdi1_y1(ii,jj,kk)

                    else if (all(tn10p1_o(e_y1-2:e_y1+3,jj,kk).eq.1.)) then
                        csdi1_y1(ii,jj,kk) = 3 + csdi1_y1(ii,jj,kk)
                        
                    else if (all(tn10p1_o(e_y1-1:e_y1+4,jj,kk).eq.1.)) then
                        csdi1_y1(ii,jj,kk) = 2 + csdi1_y1(ii,jj,kk)

                    else if (all(tn10p1_o(e_y1:e_y1+5,jj,kk).eq.1.)) then
                        csdi1_y1(ii,jj,kk) = 1 + csdi1_y1(ii,jj,kk)
                    end if 
                    end if
                    end if 
                    end if 
                    end if 

                end if 

            end do 
        end do 
    end do


    do kk = 0, nlon1r-1
        do jj = 0, nlat1r-1
            do ii = 1, n_year-1

                s_y1 = min(ind(ut1(:,0).eq.s_year+ii))

                if (.not.(all(tn10p1_o(s_y1:s_y1+5,jj,kk).eq.1.))) then 

                    if (all(tn10p1_o(s_y1-1:s_y1+4,jj,kk).eq.1.)) then 
                        csdi1_y1(ii,jj,kk) = 5 + csdi1_y1(ii,jj,kk)

                    else if (all(tn10p1_o(s_y1-2:s_y1+3,jj,kk).eq.1.)) then 
                        csdi1_y1(ii,jj,kk) = 4 + csdi1_y1(ii,jj,kk)

                    else if (all(tn10p1_o(s_y1-3:s_y1+2,jj,kk).eq.1.)) then
                        csdi1_y1(ii,jj,kk) = 3 + csdi1_y1(ii,jj,kk)
                        
                    else if (all(tn10p1_o(s_y1-4:s_y1+1,jj,kk).eq.1.)) then
                        csdi1_y1(ii,jj,kk) = 2 + csdi1_y1(ii,jj,kk)

                    else if (all(tn10p1_o(s_y1-5:s_y1,jj,kk).eq.1.)) then
                        csdi1_y1(ii,jj,kk) = 1 + csdi1_y1(ii,jj,kk)
                    end if 
                    end if
                    end if 
                    end if 
                    end if 

                end if 

            end do 
        end do 
    end do

    csdi1_y_n = csdi1_y1
    csdi1_y_n = mask(csdi1_y1, (lmask1.ge.75.),True)

;    asciiwrite(outdir+model_ns+"_"+"aus"+y_period+"_"+var_csdi+"org_nbs.txt",csdi1_y_n)
;    delete(csdi1_y1)
;    delete(csdi1_y_n)

;=======================in-base period=======================
    win_inb1 = new((/n_byear,n_byear-1,365,n_win,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)

    do ii = 0, n_byear-1
        indx_inb1 = ind(byear_arr1.ne.s_byear+ii)
        ;print(indx_inb1)
        
        do jj= 0, n_byear-1-1
            indx_inb2 = array_append_record(indx_inb1,indx_inb1(jj),0)
    
            do kk = 0, 365-1
                do mm = 0, n_byear-1

                    indx_b = min(ind(ut1(:,0).eq.byear_arr1(indx_inb2(mm))))
    
                    s_ind = 5*mm
                    e_ind = 5*mm+4
    
                    if ((isleapyear(byear_arr1(indx_inb2(mm)))).and.(kk.gt.58)) then
                        win_inb1(ii,jj,kk,s_ind:e_ind,:,:) = t1_min_ausr(indx_b-1+kk:indx_b+3+kk,:,:)
                        ;some metadata can be passed into the win_inb1
                        ;;;;;;;;;;;;;;;;;;;;;;;;;;;
                    else 
                        win_inb1(ii,jj,kk,s_ind:e_ind,:,:) = t1_min_ausr(indx_b-2+kk:indx_b+2+kk,:,:)
                    end if

                end do
            end do 
        end do
    end do

    printVarSummary(win_inb1)

;===================================================
    indx_inb1s  = dim_pqsort_n(win_inb1, 2, 3)
                
    yj_inb1 = win_inb1(:,:,:,jjj1-1,:,:)
    ;printVarSummary(yj)
    yjp_inb1 = win_inb1(:,:,:,jjj1,:,:)
    per_ind1 = (1-fff1)*yj_inb1 + fff1*yjp_inb1


    yj_inb2 = win_inb1(:,:,:,jjj2-1,:,:)
    ;printVarSummary(yj)
    yjp_inb2 = win_inb1(:,:,:,jjj2,:,:)
    per_ind2 = (1-fff2)*yj_inb2 + fff2*yjp_inb2


    delete(yj_inb1)
    delete(yjp_inb1)
    delete(yj_inb2)
    delete(yjp_inb2)
    delete(win_inb1)
    printVarSummary(per_ind1)

;=================for the day Feb 29th in in-base period===============
    if ((byear_arr1@calendar.eq."standard").and.(considerFebL.eq."considered")) then

        win_inb1_leap = new((/n_byear,n_byear-1,n_win,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)

        per_ind1_leap = new((/n_byear,n_byear-1,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)

        per_ind2_leap = new((/n_byear,n_byear-1,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)
        

        do mm = 0,n_byear-1
            indx_inb1_leap = ind(byear_arr1.ne.s_byear+mm)
        
            do nn= 0, n_byear-1-1
                indx_inb2_leap = array_append_record(indx_inb1_leap,indx_inb1_leap(nn),0)

                n_lyear = num(isleapyear(byear_arr1(indx_inb2_leap)))

                n_win_leap = n_lyear*5+(n_byear-n_lyear)*4
    
                do ii=59,59
                    do jj=0, n_byear-1
    
                        indx_b = min(ind(ut1(:,0).eq.byear_arr1(indx_inb2_leap(jj)))) 
                        s_ind = 5*jj
                        e_ind = 5*jj+4
                        
                        if (isleapyear(byear_arr1(indx_inb2_leap(jj)))) then
                            win_inb1_leap(mm,nn,s_ind:e_ind,:,:) = t1_min_ausr(indx_b-2+ii:indx_b+2+ii,:,:)
                            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                        else
                            win_inb1_leap(mm,nn,s_ind:s_ind+1,:,:) = t1_min_ausr(indx_b-2+ii:indx_b-1+ii,:,:)
                            win_inb1_leap(mm,nn,s_ind+2:s_ind+2,:,:) = t1_t@_FillValue
                            win_inb1_leap(mm,nn,s_ind+3:e_ind,:,:) = t1_min_ausr(indx_b+ii:indx_b+1+ii,:,:)
                        end if

                    end do
                end do 

                indx_inb1_leaps  = dim_pqsort_n(win_inb1_leap(mm,nn,:,:,:), 2, 0)

                jjj_inb1_leap = toint(floor(ppp1*(n_win_leap+1)))

                fff_inb1_leap = ppp1*(n_win_leap+1)-jjj_inb1_leap

                yj_inb1_leap = win_inb1_leap(mm,nn,jjj_inb1_leap-1,:,:)
                yjp_inb1_leap = win_inb1_leap(mm,nn,jjj_inb1_leap,:,:)

                per_ind1_leap(mm,nn,:,:) = (1-fff_inb1_leap)*yj_inb1_leap + fff_inb1_leap*yjp_inb1_leap


                jjj_inb2_leap = toint(floor(ppp2*(n_win_leap+1)))

                fff_inb2_leap = ppp2*(n_win_leap+1)-jjj_inb2_leap

                yj_inb2_leap = win_inb1_leap(mm,nn,jjj_inb2_leap-1,:,:)
                yjp_inb2_leap = win_inb1_leap(mm,nn,jjj_inb2_leap,:,:)

                per_ind2_leap(mm,nn,:,:) = (1-fff_inb2_leap)*yj_inb2_leap + fff_inb2_leap*yjp_inb2_leap

                delete(indx_inb1_leaps)                
                delete(yj_inb1_leap)
                delete(yjp_inb1_leap)
                delete(yj_inb2_leap)
                delete(yjp_inb2_leap)

            end do 
        end do 
        delete(win_inb1_leap)
        printVarSummary(per_ind1_leap)
    end if

;=======================get initial results for in-base period=======================
    tn10p1_indy = new((/n_byear,n_byear-1,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)
    tn90p1_indy = new((/n_byear,n_byear-1,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)
    tn10p1_ob = new((/n_byear-1,nt1r,nlat1r,nlon1r/),integer)
    tn90p1_ob = new((/n_byear-1,nt1r,nlat1r,nlon1r/),integer)
    csdi1_indy = new((/n_byear,n_byear-1,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)
    
    tn10p1_indy_per = new((/n_byear,n_byear-1,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)
    tn90p1_indy_per = new((/n_byear,n_byear-1,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)

    do jj = 0, n_byear-1-1

        t1_min_ausr1 = t1_min_ausr
        t1_min_ausr2 = t1_min_ausr

    do ii = 0, n_byear-1
    
        s_y1 = min(ind(ut1(:,0).eq.s_byear+ii))
        e_y1 = max(ind(ut1(:,0).eq.s_byear+ii))

        if ((isleapyear(byear_arr1(ii)))) then
            t1_min_ausr1(s_y1:s_y1+58,:,:) = t1_min_ausr1(s_y1:s_y1+58,:,:) - per_ind1(ii,jj,0:58,:,:)
            t1_min_ausr1(s_y1+59,:,:) = t1_min_ausr1(s_y1+59,:,:) - per_ind1_leap(ii,jj,:,:)
            t1_min_ausr1(s_y1+60:e_y1,:,:) = t1_min_ausr1(s_y1+60:e_y1,:,:) - per_ind1(ii,jj,59:,:,:)
            
            t1_min_ausr2(s_y1:s_y1+58,:,:) = t1_min_ausr2(s_y1:s_y1+58,:,:) - per_ind2(ii,jj,0:58,:,:)
            t1_min_ausr2(s_y1+59,:,:) = t1_min_ausr2(s_y1+59,:,:) - per_ind2_leap(ii,jj,:,:)
            t1_min_ausr2(s_y1+60:e_y1,:,:) = t1_min_ausr2(s_y1+60:e_y1,:,:) - per_ind2(ii,jj,59:,:,:)
        else

            t1_min_ausr1(s_y1:e_y1,:,:) = t1_min_ausr1(s_y1:e_y1,:,:) - per_ind1(ii,jj,:,:,:)
            
            t1_min_ausr2(s_y1:e_y1,:,:) = t1_min_ausr2(s_y1:e_y1,:,:) - per_ind2(ii,jj,:,:,:)
        end if

        tn10p1_ob(jj,s_y1:e_y1,:,:) = where(t1_min_ausr1(s_y1:e_y1,:,:).lt.0, 1, 0)

        tn90p1_ob(jj,s_y1:e_y1,:,:) = where(t1_min_ausr2(s_y1:e_y1,:,:).gt.0, 1, 0)
    
            tn10p1_arr = dim_numrun_n(tn10p1_ob(jj,s_y1:e_y1,:,:), 0, 0)
            tn10p1_indy(ii,jj,:,:) = tn10p1_arr(0,:,:)
            aaa = tofloat(e_y1-s_y1+1)
            tn10p1_indy_per(ii,jj,:,:) = (tn10p1_arr(0,:,:)/aaa)*100.

            tn90p1_arr = dim_numrun_n(tn90p1_ob(jj,s_y1:e_y1,:,:), 0, 0)
            tn90p1_indy(ii,jj,:,:) = tn90p1_arr(0,:,:)
            tn90p1_indy_per(ii,jj,:,:) = (tn90p1_arr(0,:,:)/aaa)*100.

    
            csdi1_arr = dim_numrun_n(tn10p1_ob(jj,s_y1:e_y1,:,:), 1, 0)
            wgt = fspan(6.,(e_y1-s_y1+1),(e_y1-s_y1+1-5))
            csdi1_indy(ii,jj,:,:) = dim_sum_wgt_n_Wrap(csdi1_arr(5:,:,:), wgt, 1, 0)
            
            delete(tn10p1_arr)
            delete(tn90p1_arr)
            delete(csdi1_arr)
            delete(wgt)

    end do
    end do 

    delete(per_ind1)
    delete(per_ind2)
    if (byear_arr1@calendar.eq."standard") then
        delete(per_ind1_leap)
        delete(per_ind2_leap)
    end if
    delete(tn90p1_ob)
    printVarSummary(t1_min_ausr)

;============CSDI for inter-year============
    do kk = 0, nlon1r-1
        do jj = 0, nlat1r-1
            do ii = 0, 9

                e_y1 = max(ind(ut1(:,0).eq.s_year+ii))

                if (.not.(all(tn10p1_o(e_y1-5:e_y1,jj,kk).eq.1.))) then 

                    if (all(tn10p1_o(e_y1-4:e_y1+1,jj,kk).eq.1.)) then 
                        csdi1_y(ii,jj,kk) = 5 + csdi1_y(ii,jj,kk)

                    else if (all(tn10p1_o(e_y1-3:e_y1+2,jj,kk).eq.1.)) then 
                        csdi1_y(ii,jj,kk) = 4 + csdi1_y(ii,jj,kk)

                    else if (all(tn10p1_o(e_y1-2:e_y1+3,jj,kk).eq.1.)) then
                        csdi1_y(ii,jj,kk) = 3 + csdi1_y(ii,jj,kk)
                        
                    else if (all(tn10p1_o(e_y1-1:e_y1+4,jj,kk).eq.1.)) then
                        csdi1_y(ii,jj,kk) = 2 + csdi1_y(ii,jj,kk)

                    else if (all(tn10p1_o(e_y1:e_y1+5,jj,kk).eq.1.)) then
                        csdi1_y(ii,jj,kk) = 1 + csdi1_y(ii,jj,kk)
                    end if 
                    end if
                    end if 
                    end if 
                    end if 

                end if 

            end do 
        end do 
    end do

    do kk = 0, nlon1r-1
        do jj = 0, nlat1r-1
            do ii = 41, n_year-1-1

                e_y1 = max(ind(ut1(:,0).eq.s_year+ii))

                if (.not.(all(tn10p1_o(e_y1-5:e_y1,jj,kk).eq.1.))) then 

                    if (all(tn10p1_o(e_y1-4:e_y1+1,jj,kk).eq.1.)) then 
                        csdi1_y(ii,jj,kk) = 5 + csdi1_y(ii,jj,kk)

                    else if (all(tn10p1_o(e_y1-3:e_y1+2,jj,kk).eq.1.)) then 
                        csdi1_y(ii,jj,kk) = 4 + csdi1_y(ii,jj,kk)

                    else if (all(tn10p1_o(e_y1-2:e_y1+3,jj,kk).eq.1.)) then
                        csdi1_y(ii,jj,kk) = 3 + csdi1_y(ii,jj,kk)
                        
                    else if (all(tn10p1_o(e_y1-1:e_y1+4,jj,kk).eq.1.)) then
                        csdi1_y(ii,jj,kk) = 2 + csdi1_y(ii,jj,kk)

                    else if (all(tn10p1_o(e_y1:e_y1+5,jj,kk).eq.1.)) then
                        csdi1_y(ii,jj,kk) = 1 + csdi1_y(ii,jj,kk)
                    end if 
                    end if
                    end if 
                    end if 
                    end if 

                end if 

            end do 
        end do 
    end do

    do kk = 0, nlon1r-1
        do jj = 0, nlat1r-1
            do ii = 1, 10

                s_y1 = min(ind(ut1(:,0).eq.s_year+ii))

                if (.not.(all(tn10p1_o(s_y1:s_y1+5,jj,kk).eq.1.))) then 

                    if (all(tn10p1_o(s_y1-1:s_y1+4,jj,kk).eq.1.)) then 
                        csdi1_y(ii,jj,kk) = 5 + csdi1_y(ii,jj,kk)

                    else if (all(tn10p1_o(s_y1-2:s_y1+3,jj,kk).eq.1.)) then 
                        csdi1_y(ii,jj,kk) = 4 + csdi1_y(ii,jj,kk)

                    else if (all(tn10p1_o(s_y1-3:s_y1+2,jj,kk).eq.1.)) then
                        csdi1_y(ii,jj,kk) = 3 + csdi1_y(ii,jj,kk)
                        
                    else if (all(tn10p1_o(s_y1-4:s_y1+1,jj,kk).eq.1.)) then
                        csdi1_y(ii,jj,kk) = 2 + csdi1_y(ii,jj,kk)

                    else if (all(tn10p1_o(s_y1-5:s_y1,jj,kk).eq.1.)) then
                        csdi1_y(ii,jj,kk) = 1 + csdi1_y(ii,jj,kk)
                    end if 
                    end if
                    end if 
                    end if 
                    end if 

                end if 

            end do 
        end do 
    end do

    do kk = 0, nlon1r-1
        do jj = 0, nlat1r-1
            do ii = 42, n_year-1

                s_y1 = min(ind(ut1(:,0).eq.s_year+ii))

                if (.not.(all(tn10p1_o(s_y1:s_y1+5,jj,kk).eq.1.))) then 

                    if (all(tn10p1_o(s_y1-1:s_y1+4,jj,kk).eq.1.)) then 
                        csdi1_y(ii,jj,kk) = 5 + csdi1_y(ii,jj,kk)

                    else if (all(tn10p1_o(s_y1-2:s_y1+3,jj,kk).eq.1.)) then 
                        csdi1_y(ii,jj,kk) = 4 + csdi1_y(ii,jj,kk)

                    else if (all(tn10p1_o(s_y1-3:s_y1+2,jj,kk).eq.1.)) then
                        csdi1_y(ii,jj,kk) = 3 + csdi1_y(ii,jj,kk)
                        
                    else if (all(tn10p1_o(s_y1-4:s_y1+1,jj,kk).eq.1.)) then
                        csdi1_y(ii,jj,kk) = 2 + csdi1_y(ii,jj,kk)

                    else if (all(tn10p1_o(s_y1-5:s_y1,jj,kk).eq.1.)) then
                        csdi1_y(ii,jj,kk) = 1 + csdi1_y(ii,jj,kk)
                    end if 
                    end if
                    end if 
                    end if 
                    end if 

                end if 

            end do 
        end do 
    end do

;===================base period===================
    do mm = 0, n_byear-1-1
    do kk = 0, nlon1r-1
        do jj = 0, nlat1r-1
            do ii = 0, n_byear-1-1

                e_y1 = max(ind(ut1(:,0).eq.s_byear+ii))

                if (.not.(all(tn10p1_ob(mm,e_y1-5:e_y1,jj,kk).eq.1.))) then 

                    if ((all(tn10p1_ob(mm,e_y1-4:e_y1+1,jj,kk).eq.1.))) then 
                        csdi1_indy(ii,mm,jj,kk) = 5 + csdi1_indy(ii,mm,jj,kk)

                    else if ((all(tn10p1_ob(mm,e_y1-3:e_y1+2,jj,kk).eq.1.))) then 
                        csdi1_indy(ii,mm,jj,kk) = 4 + csdi1_indy(ii,mm,jj,kk)

                    else if ((all(tn10p1_ob(mm,e_y1-2:e_y1+3,jj,kk).eq.1.))) then 
                        csdi1_indy(ii,mm,jj,kk) = 3 + csdi1_indy(ii,mm,jj,kk)

                    else if ((all(tn10p1_ob(mm,e_y1-1:e_y1+4,jj,kk).eq.1.))) then 
                        csdi1_indy(ii,mm,jj,kk) = 2 + csdi1_indy(ii,mm,jj,kk)

                    else if ((all(tn10p1_ob(mm,e_y1:e_y1+5,jj,kk).eq.1.))) then 
                        csdi1_indy(ii,mm,jj,kk) = 1 + csdi1_indy(ii,mm,jj,kk)
                    end if 
                    end if
                    end if 
                    end if 
                    end if 

                end if 
            end do 
        end do 
    end do
    end do

    do mm = 0, n_byear-1-1
    do kk = 0, nlon1r-1
        do jj = 0, nlat1r-1
            do ii = 1, n_byear-1

                s_y1 = min(ind(ut1(:,0).eq.s_byear+ii))

                if (.not.(all(tn10p1_ob(mm,s_y1:s_y1+5,jj,kk).eq.1.))) then 

                    if ((all(tn10p1_ob(mm,s_y1-1:s_y1+4,jj,kk).eq.1.))) then 
                        csdi1_indy(ii,mm,jj,kk) = 5 + csdi1_indy(ii,mm,jj,kk)

                    else if ((all(tn10p1_ob(mm,s_y1-2:s_y1+3,jj,kk).eq.1.))) then 
                        csdi1_indy(ii,mm,jj,kk) = 4 + csdi1_indy(ii,mm,jj,kk)

                    else if ((all(tn10p1_ob(mm,s_y1-3:s_y1+2,jj,kk).eq.1.))) then 
                        csdi1_indy(ii,mm,jj,kk) = 3 + csdi1_indy(ii,mm,jj,kk)

                    else if ((all(tn10p1_ob(mm,s_y1-4:s_y1+1,jj,kk).eq.1.))) then 
                        csdi1_indy(ii,mm,jj,kk) = 2 + csdi1_indy(ii,mm,jj,kk)

                    else if ((all(tn10p1_ob(mm,s_y1-5:s_y1,jj,kk).eq.1.))) then 
                        csdi1_indy(ii,mm,jj,kk) = 1 + csdi1_indy(ii,mm,jj,kk)
                    end if 
                    end if
                    end if 
                    end if 
                    end if 

                end if 
            end do 
        end do 
    end do
    end do


;===============================================
    csdi1_y1960 = new((/n_byear-1,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)
    csdi1_y1991 = new((/n_byear-1,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)

    csdi1_y1961 = new((/n_byear-1,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)
    csdi1_y1990 = new((/n_byear-1,nlat1r,nlon1r/),typeof(t1_t),t1_t@_FillValue)

    print(num(ismissing(csdi1_y1960)))
    csdi1_y1960 = 0.0
    csdi1_y1991 = 0.0

    csdi1_y1961 = 0.0
    csdi1_y1990 = 0.0
    print(num(ismissing(csdi1_y1960)))

    do mm = 0, n_byear-1-1
    do kk = 0, nlon1r-1
        do jj = 0, nlat1r-1
            do ii = 10, 10

                e_y1 = max(ind(ut1(:,0).eq.s_year+ii))

                if (.not.(all(tn10p1_o(e_y1-5:e_y1,jj,kk).eq.1.))) then 

                    if ((all(tn10p1_o(e_y1-4:e_y1,jj,kk).eq.1.)).and.(all(tn10p1_ob(mm,e_y1+1:e_y1+1,jj,kk).eq.1.))) then 
                        csdi1_y1960(mm,jj,kk) = 5 

                    else if ((all(tn10p1_o(e_y1-3:e_y1,jj,kk).eq.1.)).and.(all(tn10p1_ob(mm,e_y1+1:e_y1+2,jj,kk).eq.1.))) then 
                        csdi1_y1960(mm,jj,kk) = 4 

                    else if ((all(tn10p1_o(e_y1-2:e_y1,jj,kk).eq.1.)).and.(all(tn10p1_ob(mm,e_y1+1:e_y1+3,jj,kk).eq.1.))) then
                        csdi1_y1960(mm,jj,kk) = 3 
                        
                    else if ((all(tn10p1_o(e_y1-1:e_y1,jj,kk).eq.1.)).and.(all(tn10p1_ob(mm,e_y1+1:e_y1+4,jj,kk).eq.1.))) then
                        csdi1_y1960(mm,jj,kk) = 2 

                    else if ((all(tn10p1_o(e_y1:e_y1,jj,kk).eq.1.)).and.(all(tn10p1_ob(mm,e_y1+1:e_y1+5,jj,kk).eq.1.))) then
                        csdi1_y1960(mm,jj,kk) = 1 
                    end if 
                    end if
                    end if 
                    end if 
                    end if 

                end if 

            end do 
        end do 
    end do
    end do 


    do mm = 0, n_byear-1-1
    do kk = 0, nlon1r-1
        do jj = 0, nlat1r-1
            do ii = 41, 41

                s_y1 = min(ind(ut1(:,0).eq.s_year+ii))

                if (.not.(all(tn10p1_o(s_y1:s_y1+5,jj,kk).eq.1.))) then 

                    if ((all(tn10p1_o(s_y1:s_y1+4,jj,kk).eq.1.)).and.(all(tn10p1_ob(mm,s_y1-1:s_y1-1,jj,kk).eq.1.))) then 
                        csdi1_y1991(mm,jj,kk) = 5 

                    else if ((all(tn10p1_o(s_y1:s_y1+3,jj,kk).eq.1.)).and.(all(tn10p1_ob(mm,s_y1-2:s_y1-1,jj,kk).eq.1.))) then 
                        csdi1_y1991(mm,jj,kk) = 4 

                    else if ((all(tn10p1_o(s_y1:s_y1+2,jj,kk).eq.1.)).and.(all(tn10p1_ob(mm,s_y1-3:s_y1-1,jj,kk).eq.1.))) then
                        csdi1_y1991(mm,jj,kk) = 3 
                        
                    else if ((all(tn10p1_o(s_y1:s_y1+1,jj,kk).eq.1.)).and.(all(tn10p1_ob(mm,s_y1-4:s_y1-1,jj,kk).eq.1.))) then
                        csdi1_y1991(mm,jj,kk) = 2 

                    else if ((all(tn10p1_o(s_y1:s_y1,jj,kk).eq.1.)).and.(all(tn10p1_ob(mm,s_y1-5:s_y1-1,jj,kk).eq.1.))) then
                        csdi1_y1991(mm,jj,kk) = 1 
                    end if 
                    end if
                    end if 
                    end if 
                    end if 

                end if 

            end do 
        end do 
    end do
    end do 


    do mm = 0, n_byear-1-1
    do kk = 0, nlon1r-1
        do jj = 0, nlat1r-1
            do ii = 29, 29

                e_y1 = max(ind(ut1(:,0).eq.s_byear+ii))

                if (.not.(all(tn10p1_ob(mm,e_y1-5:e_y1,jj,kk).eq.1.))) then 

                    if ((all(tn10p1_ob(mm,e_y1-4:e_y1,jj,kk).eq.1.)).and.(all(tn10p1_o(e_y1+1:e_y1+1,jj,kk).eq.1.))) then 
                        csdi1_y1990(mm,jj,kk) = 5 

                    else if ((all(tn10p1_ob(mm,e_y1-3:e_y1,jj,kk).eq.1.)).and.(all(tn10p1_o(e_y1+1:e_y1+2,jj,kk).eq.1.))) then 
                        csdi1_y1990(mm,jj,kk) = 4 

                    else if ((all(tn10p1_ob(mm,e_y1-2:e_y1,jj,kk).eq.1.)).and.(all(tn10p1_o(e_y1+1:e_y1+3,jj,kk).eq.1.))) then 
                        csdi1_y1990(mm,jj,kk) = 3 

                    else if ((all(tn10p1_ob(mm,e_y1-1:e_y1,jj,kk).eq.1.)).and.(all(tn10p1_o(e_y1+1:e_y1+4,jj,kk).eq.1.))) then 
                        csdi1_y1990(mm,jj,kk) = 2 

                    else if ((all(tn10p1_ob(mm,e_y1:e_y1,jj,kk).eq.1.)).and.(all(tn10p1_o(e_y1+1:e_y1+5,jj,kk).eq.1.))) then 
                        csdi1_y1990(mm,jj,kk) = 1 
                    end if 
                    end if
                    end if 
                    end if 
                    end if 

                end if 
            end do 
        end do 
    end do
    end do


    do mm = 0, n_byear-1-1
    do kk = 0, nlon1r-1
        do jj = 0, nlat1r-1
            do ii = 0, 0

                s_y1 = min(ind(ut1(:,0).eq.s_byear+ii))

                if (.not.(all(tn10p1_ob(mm,s_y1:s_y1+5,jj,kk).eq.1.))) then 

                    if ((all(tn10p1_ob(mm,s_y1:s_y1+4,jj,kk).eq.1.)).and.(all(tn10p1_o(s_y1-1:s_y1-1,jj,kk).eq.1.))) then 
                        csdi1_y1961(mm,jj,kk) = 5 

                    else if ((all(tn10p1_ob(mm,s_y1:s_y1+3,jj,kk).eq.1.)).and.(all(tn10p1_o(s_y1-2:s_y1-1,jj,kk).eq.1.))) then 
                        csdi1_y1961(mm,jj,kk) = 4 

                    else if ((all(tn10p1_ob(mm,s_y1:s_y1+2,jj,kk).eq.1.)).and.(all(tn10p1_o(s_y1-3:s_y1-1,jj,kk).eq.1.))) then 
                        csdi1_y1961(mm,jj,kk) = 3 

                    else if ((all(tn10p1_ob(mm,s_y1:s_y1+1,jj,kk).eq.1.)).and.(all(tn10p1_o(s_y1-4:s_y1-1,jj,kk).eq.1.))) then 
                        csdi1_y1961(mm,jj,kk) = 2 

                    else if ((all(tn10p1_ob(mm,s_y1:s_y1,jj,kk).eq.1.)).and.(all(tn10p1_o(s_y1-5:s_y1-1,jj,kk).eq.1.))) then 
                        csdi1_y1961(mm,jj,kk) = 1 
                    end if 
                    end if
                    end if 
                    end if 
                    end if 

                end if 
            end do 
        end do 
    end do
    end do


    delete(tn10p1_o)
    delete(tn10p1_ob)
    printVarSummary(csdi1_indy)

    print(num(ismissing(csdi1_y1960)))

;======================================================
    tn10p1_indym = dim_avg_n(tn10p1_indy,1)
    delete(tn10p1_indy)
    tn90p1_indym = dim_avg_n(tn90p1_indy,1)
    delete(tn90p1_indy)

    tn10p1_indym_per = dim_avg_n(tn10p1_indy_per,1)
    delete(tn10p1_indy_per)
    tn90p1_indym_per = dim_avg_n(tn90p1_indy_per,1)
    delete(tn90p1_indy_per)

    csdi1_indym = dim_avg_n(csdi1_indy,1)
    delete(csdi1_indy)
    csdi1_y1960m = dim_avg_n(csdi1_y1960,0)
    csdi1_y1961m = dim_avg_n(csdi1_y1961,0)
    csdi1_y1990m = dim_avg_n(csdi1_y1990,0)
    csdi1_y1991m = dim_avg_n(csdi1_y1991,0)

    print(num(ismissing(csdi1_y1960m)))
    print(num(ismissing(csdi1_y(10,:,:))))
    

    tn10p1_y(11:40,:,:) = tn10p1_indym
    tn90p1_y(11:40,:,:) = tn90p1_indym

    tn10p1_y_per(11:40,:,:) = tn10p1_indym_per
    tn90p1_y_per(11:40,:,:) = tn90p1_indym_per

    csdi1_y(11:40,:,:) = csdi1_indym

    csdi1_y(10,:,:) = csdi1_y(10,:,:) + csdi1_y1960m(:,:)
    csdi1_y(11,:,:) = csdi1_y(11,:,:) + csdi1_y1961m(:,:)
    csdi1_y(40,:,:) = csdi1_y(40,:,:) + csdi1_y1990m(:,:)
    csdi1_y(41,:,:) = csdi1_y(41,:,:) + csdi1_y1991m(:,:)

    print(num(ismissing(csdi1_y(10,:,:))))

;==============Masking again=============
    tn10p1_y_n = tn10p1_y
    tn10p1_y_n = mask(tn10p1_y, (lmask1.ge.75.),True)
    
    tn90p1_y_n = tn90p1_y
    tn90p1_y_n = mask(tn90p1_y, (lmask1.ge.75.),True)

    tn10p1_y_pern = tn10p1_y_per
    tn10p1_y_pern = mask(tn10p1_y_per, (lmask1.ge.75.),True)
    
    tn90p1_y_pern = tn90p1_y_per
    tn90p1_y_pern = mask(tn90p1_y_per, (lmask1.ge.75.),True)

    csdi1_y_n = csdi1_y
    csdi1_y_n = mask(csdi1_y, (lmask1.ge.75.),True)

    delete(tn10p1_indym)
    delete(tn90p1_indym)
    delete(csdi1_indym)
    printVarSummary(tn10p1_y)

;==============generating ascii files===============
;    asciiwrite(outdir+model_ns+"_"+"aus"+y_period+"_"+var_tn10p+"_num.txt",tn10p1_y_n)
;    asciiwrite(outdir+model_ns+"_"+"aus"+y_period+"_"+var_tn90p+"_num.txt",tn90p1_y_n)
    asciiwrite(outdir+model_ns+"_"+"aus"+y_period+"_"+var_csdi+".txt",csdi1_y_n)

    asciiwrite(outdir+model_ns+"_"+"aus"+y_period+"_"+var_tn10p+".txt",tn10p1_y_pern)
    asciiwrite(outdir+model_ns+"_"+"aus"+y_period+"_"+var_tn90p+".txt",tn90p1_y_pern)

    print("===================end of the code===================")
end

