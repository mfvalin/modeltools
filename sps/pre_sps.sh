#!/bin/bash
#
[[ -r exper.cfg ]] || { echo "ERROR: cannot find $(pwd -P)/exper.cfg" ; exit 1 ; }
source ./exper.cfg
#
[[ -n "${exper_current_date}" ]] || { echo "ERROR: exper_current_date not set" ; exit 1 ; }
#
Delta=${exper_delta:-1month}
#
StepStartDate=${exper_current_date}
StepEndDate="$(date -d${StepStartDate}+${Delta} +%Y%m%d)"
echo "INFO: running from ${StepStartDate} to ${StepEndDate}, (${Delta})"
#
# get initial conditions files for this month (better be available or else !!)
#
rm -f Data/Input/anal     # get rid of old file
for Target in ${exper_anal1} ${exper_anal2}
do
  [[ -f ${Target}_${StepStartDate} ]] && cp ${Target}_${StepStartDate} Data/Input/anal && echo "INFO: using ${Target}_${StepStartDate}" && break
done
[[ -f Data/Input/anal ]] || { echo "ERROR: cound not find initial conditions file for ${StepStartDate}" ; exit 1 ; }
#
# get driving data files for this month (better be available or else !!)
#
#ls -l ${exper_depot1}/*_${StepStartDate%??} ${exper_depot2}/*_${StepStartDate%??}
Nfiles=0
for Target in $(ls -1 ${exper_depot1}/*_${StepStartDate%??} ${exper_depot2}/*_${StepStartDate%??} )
do
  Target2=${Target##*/}
  if [[ -r Data/Input/inrep/${Target2} ]] ; then
    echo "INFO: ${Target2} found in Data/Input/inrep/"
  else
    cp ${Target} Data/Input/inrep/.     # file not already there,  copy it
    echo "INFO: using ${Target}"
  fi
  ((Nfiles=Nfiles+1))
done
if ((Nfiles == 0));then
  echo "ERROR: no driving data found for ${StepStartDate%??}"
  exit 1
fi
#ls -l Data/Input/inrep/*_${StepStartDate%??}
#
# get driving data files for next month (if available and not there)
#
#ls -l ${exper_depot2}/*_${StepEndDate%??} ${exper_depot2}/*_${StepEndDate%??}
for Target in $(ls -1 ${exper_depot1}/*_${StepEndDate%??} ${exper_depot2}/*_${StepEndDate%??} )
do
  Target2=${Target##*/}
  if [[ -r Data/Input/inrep/${Target2} ]] ; then
    echo "INFO: ${Target2} found in Data/Input/inrep/"
  else
    cp ${Target} Data/Input/inrep/.    # file not already there,  copy it
    echo "INFO: using ${Target}"
  fi
done
#ls -l Data/Input/inrep/*_${StepStartDate%??}
#
Date1=$(date -d${StepStartDate}GMT0 +%s)
Date2=$(date -d${StepEndDate}GMT0 +%s)
Nsteps=$(((Date2-Date1)/900))
Nhours=$(((Date2-Date1)/3600))
echo INFO: performing ${Nsteps} timesteps in ${Nhours} hours integration   file=pm${DaTe}000000-??-??_000000h
#
# fix Step_runstrt_S and Step_total in sps.cfg
#
find . -mindepth 3 -maxdepth 3 -name sps.cfg -exec cp {} sps.cfg_old \;
sed -e "s/Step_runstrt_S =[^.]*/Step_runstrt_S = '${StepStartDate}/" -e "s/Step_total.*/Step_total = ${Nsteps}/" <sps.cfg_old >sps.cfg_new
find . -mindepth 3 -maxdepth 3 -name sps.cfg -exec mv sps.cfg_new {} \;
rm -f sps.cfg_old
exit 0
#
