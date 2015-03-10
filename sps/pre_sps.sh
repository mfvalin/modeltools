#!/bin/ksh93
#
[[ -r exper.cfg ]] || { echo "ERROR: cannot find $(pwd -P)/exper.cfg" ; exit 1 ; }
source ./exper.cfg
#
Delta=${exper_delta:-1month}
exper_current_date=${exper_current_date:-${exper_start_date}}
CurrentDate=${exper_current_date}
EndDate="$(date -d${CurrentDate}+${Delta} +%Y%m%d)"
echo "INFO: running from ${CurrentDate} to ${EndDate}, (${Delta})"
#
# get initial_condition file
#
rm -f Data/Input/anal     # get rid of old file
for Target in ${exper_anal1} ${exper_anal2}
do
  [[ -f ${Target}_${CurrentDate} ]] && cp ${Target}_${CurrentDate} Data/Input/anal && echo "INFO: using ${Target}_${CurrentDate}" && break
done
[[ -f Data/Input/anal ]] || { echo "ERROR: cound not find initial conditions file for ${CurrentDate}" ; exit 1 ; }
#
# files for this month (better be available or else !!)
#
#ls -l ${exper_depot1}/*_${CurrentDate%??} ${exper_depot2}/*_${CurrentDate%??}
Nfiles=0
for Target in $(ls -1 ${exper_depot1}/*_${CurrentDate%??} ${exper_depot2}/*_${CurrentDate%??} )
do
  Target2=${Target##*/}
  if [[ -r Data/Input/inrep/${Target2} ]] ; then
    echo "INFO: ${Target2} found in Data/Input/inrep/"
  else
    cp ${Target} Data/Input/inrep/.
  fi
  ((Nfiles=Nfiles+1))
done
if ((Nfiles == 0)) then
  echo "ERROR: no more driving data"
  exit 1
fi
#ls -l Data/Input/inrep/*_${CurrentDate%??}
#
# files for next month if available
#
#ls -l ${exper_depot2}/*_${EndDate%??} ${exper_depot2}/*_${EndDate%??}
for Target in $(ls -1 ${exper_depot1}/*_${EndDate%??} ${exper_depot2}/*_${EndDate%??} )
do
  Target2=${Target##*/}
  if [[ -r Data/Input/inrep/${Target2} ]] ; then
    echo "INFO: ${Target2} found in Data/Input/inrep/"
  else
    cp ${Target} Data/Input/inrep/.
  fi
done
#ls -l Data/Input/inrep/*_${CurrentDate%??}
#
grep -v exper_current_date exper.cfg >exper.cfg_new
mv exper.cfg_new exper.cfg
#echo "exper_current_date=${EndDate}" >>exper.cfg
echo "exper_current_date=${exper_current_date}" >>exper.cfg
#
Date1=$(date -d${CurrentDate}GMT0 +%s)
Date2=$(date -d${EndDate}GMT0 +%s)
Nsteps=$(((Date2-Date1)/900))
typeset -Z6 Nhours
Nhours=$(((Date2-Date1)/3600))
echo INFO: performing ${Nsteps} timesteps in ${Nhours} hours integration   file=pm${DaTe}000000-??-??_000000h
find . -mindepth 3 -maxdepth 3 -name sps.cfg -exec cp {} sps.cfg_old \;
sed -e "s/Step_runstrt_S =[^.]*/Step_runstrt_S = '${CurrentDate}/" -e "s/Step_total.*/Step_total = ${Nsteps}/" <sps.cfg_old >sps.cfg_new
find . -mindepth 3 -maxdepth 3 -name sps.cfg -exec mv sps.cfg_new {} \;
exit 0
#
