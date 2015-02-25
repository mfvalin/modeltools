#!/bin/ksh93
Delta=${1:-1month}
CurrentDate=${2:-$(grep Step_runstrt_S sps.cfg | sed -e "s/[.].*//" -e "s/.*'//")}
DaTe="$(date -d${CurrentDate}+${Delta} +%Y%m%d)"
Where=$(true_path ${0})
What=${Where%.sh}_new
echo INFO: setting data to ${DaTe} from ${CurrentDate} in ${What}
Date1=$(date -d${DaTe} +%s)
Date2=$(date -d${DaTe}+${Delta#-} +%s)
Nsteps=$(((Date2-Date1)/900))
typeset -Z6 Nhours
Nhours=$(((Date2-Date1)/3600))
echo INFO: there are ${Nsteps} timesteps in this integration ending at $(date -d${DaTe}+${1:-1month} +%Y%m%d) file=pm${DaTe}000000.....${Nhours}h
cat <<EOT >${What}
version=100

@grid_cfgs
Grd_typ_S     = 'LU'
Grd_dx        =    0.25
Grd_dy          =    0.25
Grd_ni        =  144   
Grd_nj          =  115
Grd_iref      =   82   
Grd_jref        =   60
Grd_latr      =    0.0 
Grd_lonr        =  180.0
Grd_xlat1     =   57.5 
Grd_xlon1       = -130.
Grd_xlat2     =    0.  
Grd_xlon2       =  -40.

@time_cfgs
Step_runstrt_S = '${DaTe}.000000'
Step_dt    = 900.
Step_gstat = 12
Step_bkup  = -1
Step_total = ${Nsteps}

@levels_cfgs
Lvl_typ_S = 'HU'
Lvl_ptop_8 = 10.0
Lvl_NoTopThL_L = .false.
Lvl_Tlift_L = .false.

@sps_cfgs
ip1a = 93423264
int_accu_S = 'CONST'
adapt_L = .FALSE.
lapserate = 0.0065
read_hu_L = .TRUE.
@

&physics_cfgs
 schmsol = 'ISBA' ,

 glaciers   = 'DYN'      ,
 gl_dyn     = .true.     ,
 gl_bins    = 8          ,
 gl_snlays  = 1          ,
 gl_ilays   = 2          ,
 gl_prdfac  = 2.0        ,
 gl_prscalf = .true.     ,

 zta=2.0,
 zua=10.0,
 LIMSNODP  = .true.,
 ICEMELT  = .true.,
/
EOT
