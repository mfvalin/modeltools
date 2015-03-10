#!/bin/bash
set -x
#
export storage_model=$(readlink -e storage_model)
echo storage_model=${storage_model}
[[ -d "${storage_model}" ]] || { echo "ERROR: $storage_model} does not exist" ; exit 1 ; }
#
while [[ yes == yes ]]
do
  [[ -r exper.cfg ]] || { echo "ERROR: cannot find $(pwd -P)/exper.cfg" ; exit 1 ; }
  source ./exper.cfg
  #
  if [[ "${exper_current_date}" == "${exper_end_date}" ]] ; then
    echo "INFO: last date reached: ${exper_end_date}"
    exit 0
  fi
  #
  ./pre_sps.sh  || { echo "ERROR: pre_sps failed" ; exit 1 ; }
  #
  echo TMPDIR=${TMPDIR}
  sps.ksh --verbosity=error --ptopo=${exper_cpu_config:-1x1x1} >sps_${exper_current_date:-${exper_start_date}}.lst 2>&1 || { echo "ERROR: sps.ksh failed" ; exit 1 ; }
  #
  ./post_sps.sh  || { echo "ERROR: post_sps failed" ; exit 1 ; }
  #
done
