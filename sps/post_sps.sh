#!/bin/ksh93
#
[[ -r exper.cfg ]] || { echo "ERROR: cannot find $(pwd -P)/exper.cfg" ; exit 1 ; }
source ./exper.cfg
#
[[ -n "${exper_current_date}" ]] || { echo "ERROR: exper_current_date not set" ; exit 1 ; }
#
Extension=""
((exper_cycle_year>0)) && Extension="$(printf '_%3.3d' ${exper_cycle_year})"
#
Delta=${exper_delta:-1month}
# exper_current_date MUST be defined at this point
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
# extract last time frame to be used as part of initial conditions for next run
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
#
# get dynamic initial conditions for next run from pilot files
#
cat <<EOT >editfst2.dir
 DESIRE(-1,-1,-1,${CurrentCmcStamp},-1,-1,-1)
 ZAP('A',-1,-1,-1,-1,-1,-1)
 END
EOT
echo editfst -s Data/Input/inrep/*${FileDate%??}  -d OUT/${exper}_anal -i editfst2.dir
ls -l Data/Input/inrep/*${FileDate%??}
editfst -s Data/Input/inrep/*${FileDate%??}  -d OUT/${exper}_anal -i editfst2.dir
#
# we can now get rid of previous pilot files
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
cp OUT/${exper}_anal Data/Input/anal
rm -f OUT/${exper}_anal
#
[[ -r ${exper_archive}/${exper}/anal_depart_${CurrentDate}${Extension} ]] || { echo "ERROR: failed to create initial conditions file anal_depart_${CurrentDate}${Extension}" ; exit 1; }
#
# update exper_current_date
#
grep -v exper_current_date exper.cfg >exper.cfg_new
mv exper.cfg_new exper.cfg
echo "exper_current_date=${CurrentDate}" >>exper.cfg
###################################################################################
#  if we processed month no 12, bump year count (exper_year_count)
###################################################################################
#
