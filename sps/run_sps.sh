#!/bin/bash
[[ "$1" == -v* ]] && set -x
#
if [[ -d storage_model ]] ; then
  storage_model=$(readlink -e storage_model)
fi
grep -q storage_model exper.cfg || echo "storage_model=${storage_model}" >>./exper.cfg
export storage_model
#
[[ -r exper.cfg ]] || { echo "ERROR: cannot find $(pwd -P)/exper.cfg" ; exit 1 ; }
source ./exper.cfg
#
# make sure that there is a value for exper_current_date, exper_fold_date and storage_model in configuration file
#
[[ -z ${exper_current_date} ]] && exper_current_date=${exper_start_date} && echo "exper_current_date=${exper_current_date}" >>exper.cfg
#
[[ -z ${exper_fold_date} ]] && exper_fold_date="$(date -d${exper_end_date}+1year  +%Y%m%d)" && echo "exper_fold_date=${exper_fold_date}" >>./exper.cfg
#
[[ -d "${storage_model}" ]] || { echo "ERROR: ${storage_model} does not exist" ; exit 1 ; }
#
while true
do
  source ./exper.cfg
  #
  if [[ "${exper_current_date}" == "${exper_end_date}" ]] ; then
    echo "INFO: last date reached: ${exper_end_date}"
    exit 0
  fi
  #
  ./pre_sps.sh  || { echo "ERROR: pre_sps failed" ; exit 1 ; }
  #
  echo "INFO: sps.ksh ${exper_cpu_config}"
  sps.ksh ${exper_cpu_config} >sps_${exper_current_date:-${exper_start_date}}.lst 2>&1 || { echo "ERROR: sps.ksh failed" ; exit 1 ; }
  #
  ./post_sps.sh  || { echo "ERROR: post_sps failed" ; exit 1 ; }
  #
done
