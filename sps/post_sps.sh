#!/bin/bash
#
[[ -r exper.cfg ]] || { echo "ERROR: cannot find $(pwd -P)/exper.cfg" ; exit 1 ; }
source ./exper.cfg
[[ -n ${exper_current_date} ]] || { echo "ERROR: exper_current_date not found in exper.cfg" ; exit 1 ; }
#
source functions_sps.dot
#
# extension is normally used when running the same driving data over and over (perpetual year or the like)
#
Extension=""
((exper_current_year>0)) && Extension="$(printf '_%3.3d' ${exper_current_year})"
((exper_current_year>0)) && Extension2="$(printf '_%3.3d' $((exper_current_year+1)))"
#
Delta=${exper_delta:-1month}
FileDate=${exper_current_date}
#
CurrentDate=$(date -d${FileDate}+${Delta} +%Y%m%d)                 # bump exper_current_date by Delta
CurrentCmcStamp=$(r.date ${CurrentDate})
echo "INFO: looking for output valid at ${CurrentDate}, CMC stamp = ${CurrentCmcStamp}"
#
Date2=$(date -d${CurrentDate}GMT0 +%s)
Date1=$(date -d${FileDate}GMT0 +%s)
Nhours=$(((Date2-Date1)/3600))     # number of hours for this integration
#
# patch together the output file (bemol)
#
rm -f OUT/${exper}_${FileDate}
echo bemol -src OUT/*/pm${FileDate}000000-??-??_000000h -dst OUT/${exper}_${FileDate}
ls -l OUT/*/pm${FileDate}000000-??-??_000000h
bemol -src OUT/*/pm${FileDate}000000-??-??_000000h -dst OUT/${exper}_${FileDate} >/dev/null
#
# copy output file to archive after making it read only (background copy)
#
mkdir -p ${exper_archive}/${exper}
chmod 444 OUT/${exper}_${FileDate}
cp OUT/${exper}_${FileDate} ${exper_archive}/${exper}/${exper}_${FileDate}${Extension} &
#
# extract last time frame to be used as part of initial conditions for next integration
#
rm -f OUT/${exper}_anal
cat <<EOT >editfst1.dir
 DESIRE(-1,-1,-1,-1,-1,${Nhours},-1))
 DESIRE(-1,['>>','^^'],-1,-1,-1,-1,-1)
 ZAP('A',-1,-1,-1,-1,-1,-1)
 END
EOT
echo editfst -s OUT/${exper}_${FileDate} -d OUT/${exper}_anal -i editfst1.dir
ls -l OUT/${exper}_${FileDate}
editfst -s OUT/${exper}_${FileDate} -d OUT/${exper}_anal -i editfst1.dir
rm -f editfst1.dir
#
# get dynamic initial conditions for next integration from driving data file
#
cat <<EOT >editfst2.dir
 DESIRE(-1,-1,-1,${CurrentCmcStamp},-1,-1,-1)
 DESIRE(-1,['>>','^^'],-1,-1,-1,-1,-1)
 ZAP('A',-1,-1,-1,-1,-1,-1)
 END
EOT
echo editfst -s Data/Input/inrep/*${FileDate%??}  -d OUT/${exper}_anal -i editfst2.dir
ls -l Data/Input/inrep/*${FileDate%??}
editfst -s Data/Input/inrep/*${FileDate%??}  -d OUT/${exper}_anal -i editfst2.dir
rm -f editfst2.dir
#
# we can now get rid of previous driving data files
#
echo rm Data/Input/inrep/*${FileDate%??}
rm Data/Input/inrep/*${FileDate%??}
#
wait   # for copy to archive completion
rm -f OUT/${exper}_${FileDate}
rm -f Data/Input/anal
set -x
cp OUT/${exper}_anal ${exper_archive}/${exper}/anal_depart_${CurrentDate}${Extension}
chmod 444 ${exper_archive}/${exper}/anal_depart_${CurrentDate}${Extension}
rm -f Data/Input/anal_${CurrentDate}
cp OUT/${exper}_anal Data/Input/anal_${CurrentDate}
rm -f OUT/${exper}_anal
#
[[ -r ${exper_archive}/${exper}/anal_depart_${CurrentDate}${Extension} ]] || \
   { echo "ERROR: failed to create initial conditions file anal_depart_${CurrentDate}${Extension}" ; exit 1; }
#
# update exper_current_date
#
update_cfg exper.cfg exper_current_date ${CurrentDate}
#
# if CurrentDate is january 1st, we have done a year of integration, bump down exper_cycle_year
# bump extension for initial conditions file
#
if [[ ${CurrentDate} == *0101 ]] ; then
  mv ${exper_archive}/${exper}/anal_depart_${CurrentDate}${Extension} ${exper_archive}/${exper}/anal_depart_${CurrentDate}${Extension2}
  update_cfg exper.cfg exper_cycle_year $((${exper_cycle_year:-999999}-1))
  update_cfg exper.cfg exper_current_year $((${exper_current_year:--1}+1))
fi
