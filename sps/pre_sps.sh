#!/bin/ksh93
#
[[ -r exper.cfg ]] || { echo "ERROR: cannot find $(pwd -P)/exper.cfg" ; exit 1 ; }
source ./exper.cfg
#
Delta=${exper_delta:-1month}
exper_current_date=${exper_current_date:-${exper_start_date}}
#
grep -v exper_current_date exper.cfg >exper.cfg_new
mv exper.cfg_new exper.cfg
# make sure that there is a value for exper_current_date in configuration file
echo "exper_current_date=${exper_current_date}" >>exper.cfg
#
StepStartDate=${exper_current_date}
StepEndDate="$(date -d${StepStartDate}+${Delta} +%Y%m%d)"
echo "INFO: running from ${StepStartDate} to ${StepEndDate}, (${Delta})"
#
# initial conditions and driving data files for this month (better be available or else !!)
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
  fi
  ((Nfiles=Nfiles+1))
done
if ((Nfiles == 0)) then
  echo "ERROR: no driving data found for ${StepStartDate%??}"
  exit 1
fi
#ls -l Data/Input/inrep/*_${StepStartDate%??}
#
# get driving data files for next month if available
#
#ls -l ${exper_depot2}/*_${StepEndDate%??} ${exper_depot2}/*_${StepEndDate%??}
for Target in $(ls -1 ${exper_depot1}/*_${StepEndDate%??} ${exper_depot2}/*_${StepEndDate%??} )
do
  Target2=${Target##*/}
  if [[ -r Data/Input/inrep/${Target2} ]] ; then
    echo "INFO: ${Target2} found in Data/Input/inrep/"
  else
    cp ${Target} Data/Input/inrep/.    # file not already there,  copy it
  fi
done
#ls -l Data/Input/inrep/*_${StepStartDate%??}
#
Date1=$(date -d${StepStartDate}GMT0 +%s)
Date2=$(date -d${StepEndDate}GMT0 +%s)
Nsteps=$(((Date2-Date1)/900))
typeset -Z6 Nhours
Nhours=$(((Date2-Date1)/3600))
echo INFO: performing ${Nsteps} timesteps in ${Nhours} hours integration   file=pm${DaTe}000000-??-??_000000h
find . -mindepth 3 -maxdepth 3 -name sps.cfg -exec cp {} sps.cfg_old \;
sed -e "s/Step_runstrt_S =[^.]*/Step_runstrt_S = '${StepStartDate}/" -e "s/Step_total.*/Step_total = ${Nsteps}/" <sps.cfg_old >sps.cfg_new
find . -mindepth 3 -maxdepth 3 -name sps.cfg -exec mv sps.cfg_new {} \;
rm -f sps.cfg_old
exit 0
#