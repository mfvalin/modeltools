#!/bin/bash
PBS_NODEFILE=${PBS_NODEFILE:-${GECOSHEP_HOSTS_FILE:-${GECOSHEP_HOSTFILE}}}  # takes care of PBS/TORQUE and Grid Engine
# =============================================================================================
function usage {
echo " usage: ${0##*/} [-debug level] [-instances number_of|@map_file] [-cpumap cpu_specification] -name stream_name -maxidle n_seconds -queues spec_1 spec_2 ... spec_N"
echo "        options are POSITIONAL, queue spec of the form: [instance#:]queue_name "
}
# =============================================================================================
function map_instance_nodes {
#  usage: map_instance_nodes JobNodeFile N1 C1 [N2 C2] ..... [Ni Ci]
JobNodeFile="$1"
shift
MAP="$@"
#
typeset -a HostList
typeset -a HostPop
typeset -a HostTiles
HostNumber=0
set -- $(uniq -c $JobNodeFile | xargs)
while [[ -n "${2}" ]] ; do
  ((HostNumber=HostNumber+1))
  HostList[$HostNumber]=$2
  HostPop[$HostNumber]=$1
  HostTiles[$HostNumber]=0
  shift; shift;
done
set -- ${MAP}
while [[ -n $2 ]] ; do
  if [[ "$1" == all || "$1" == ALL ]] ; then
    for i in $(seq ${HostNumber}) ; do [[ ${HostTiles[$i]} == 0 ]] && HostTiles[$i]=$2 ; done
  else
    HostTiles[$1]=$2
  fi
  shift; shift;
done
#
for j in  $(seq ${HostNumber}) ; do
  for i in $(seq ${HostTiles[$j]}) ; do echo ${HostList[$j]} ; done
  shift; shift;
done
}  # end of function map_instance_nodes
# =============================================================================================
function InstanceEcho {
  echo "$@" | sed "s/^/ I-0${CurrentInstance}: /"
}
# =============================================================================================
function InstanceCat {
  cat "$@" | sed "s/^/ I-0${CurrentInstance}: /"
}
# =============================================================================================
#
#============================= PASS 1 and 2  =============================
#
MyProgramName="$0"
BATCH_MPI_CPUS=${BATCH_MPI_CPUS:-1}
OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
((BATCH_MPI_CPUS=BATCH_MPI_CPUS*OMP_NUM_THREADS))
OMP_NUM_THREADS=1
#
JobNodeFiles="${HOME}/.job_queues/JobNodeFiles"
mkdir -p ${JobNodeFiles}
JobNodeFilesBase=${JobNodeFiles}/$(hostname)_$$
#
DebugFlag=0
[[ "$1" == -debug ]] && DebugFlag=$2 && shift && shift
#
[[ "$1" == -pass2 ]] && Pass2Flag=YES  && shift
if [[ "$1" != -instances && -z ${Pass2Flag} ]] ; then  # -instances allowed only on PASS 1
  set -- -instances 1 "$@"
  echo SETTING INSTANCES to 1
fi
#
#============================= PASS 1 only  =============================
if [[ "$1" == -instances ]]
then
  [[ -z "$2" ]] && usage && exit 1

  if [[ -f "${2#@}" ]] ; then  # we have a map file
    MapFile="${2#@}"
    Instances="$(cat $MapFile | wc -l)"
    if ((DebugFlag>0)) ; then
      echo ==== using $(sort -u ${PBS_NODEFILE} | wc -l) nodes , ${Instances} Instances ====
      cat $MapFile
      echo =======================
    fi
  else                         # we have a number of instances
    MapFile=""
    ( ((Instances=$2)) ) 2>/dev/null || { echo "ERROR: bad number of instances specified : $2" ; exit 1; }
  fi
  shift ; shift
  MyArguments="$@"
  ((BATCH_MPI_CPUS>=Instances)) || ((Instances=BATCH_MPI_CPUS))
#
  if [[ -f "$MapFile" ]] ; then  # we have a map file
    for Instance in $(seq $Instances) ; do
      MAP="$(cat $MapFile | sed -n ${Instance},${Instance}p | sed 's/#.*//' )"
      map_instance_nodes ${PBS_NODEFILE} ${MAP:-ALL 1} >${JobNodeFilesBase}.${Instance}
    done
  else                         # we have a number of instances
    ((BATCH_MPI_CPUS=BATCH_MPI_CPUS/Instances))
    ((Start=1))
    ((Limit=Start+BATCH_MPI_CPUS-1))
    for Instance in $(seq $Instances) ; do
      ((DebugFlag>1)) && echo Start=$Start, Limit=$Limit
      cat ${PBS_NODEFILE} | sed -n "${Start},${Limit}p" >${JobNodeFilesBase}.${Instance}
      ((Start=Start+BATCH_MPI_CPUS))
      ((Limit=Start+BATCH_MPI_CPUS-1))
    done
  fi
fi
#
if [[ -z "${Pass2Flag}" ]] ; then # pass1
  echo "INFO: starting ${Instances} instance(s) of ${MyProgramName}"
  export CurrentInstance=0
#
((DebugFlag>5)) && /sb/home/valin/test/SUPERJOBS/job_instance_job # petit test
#
  if ((DebugFlag>0)) ; then
    InstanceEcho ================================================================
  fi
#
  for Instance in $(seq $Instances) ; do
    CurrentInstance=${Instance}
#   environment variables PARALLEL_NODEFILE and CurrentInstance set by PASS 1, used by PASS 2
    if ((DebugFlag>0)) ; then
      InstanceEcho "PARALLEL_NODEFILE=${JobNodeFilesBase}.${Instance} ${MyProgramName} -debug ${DebugFlag} -pass2 ${MyArguments}"
    fi
    PARALLEL_NODEFILE=${JobNodeFilesBase}.${Instance} \
        ${MyProgramName} -debug ${DebugFlag} -pass2 ${MyArguments} 1>${ListingFile}.${CurrentInstance} 2>&1 &
  done
  CurrentInstance=0
  ((DebugFlag>0)) && InstanceEcho ================================================================
#
  wait
  echo INFO: all instances done
#
#
  exit 0
#============================== PASS 1 ends here ==============================
fi
#============================= PASS 2 begins here =============================
sleep $((CurrentInstance))
#PBS_NODEFILE=${PARALLEL_NODEFILE:-${PBS_NODEFILE}}
export BATCH_MPI_CPUS="$(wc -l <${PARALLEL_NODEFILE})"
export OMP_NUM_THREADS=1
#
#((DebugFlag>5)) && /sb/home/valin/test/SUPERJOBS/job_instance_job # petit test
#
if ((DebugFlag>0)) ; then
  InstanceEcho ================================================================
  InstanceEcho "$0 -debug ${DebugFlag} $@"
  ((DebugFlag>1)) && InstanceCat ${PARALLEL_NODEFILE}
fi
#
# collect optional -cpumap option:   -cpumap [instance:][mincpus],[maxcpus]
#               no instance: all instances
#               no mincpus : 1
#               no maxcpus : number of cpus available to the instance
# not specifying the option is equivalent to "-cpumap ,"
#
export MinCpusInJob=1                    # do not execute a job that uses less than this number of cpus
export MaxCpusInJob=${BATCH_MPI_CPUS}    # do not execute a job that uses more than this number of cpus
if [[ "$1" == -cpumap ]] ; then          # process optional -cpumap argument
  shift
  while [[ -n "$1" && "$1" != -* ]] ; do
    [[ "$1" != *:* || "${1%%:*}" == ${CurrentInstance} ]] && CpusInJob=${1##*:}
    shift
  done
  MinCpusInJob=${CpusInJob%%,*}
  MinCpusInJob=${MinCpusInJob:-1}
  MaxCpusInJob=${CpusInJob##*,}
  MaxCpusInJob=${MaxCpusInJob:-${BATCH_MPI_CPUS}}
fi
#
# collect mandatory name option
#
[[ "$1" != -name ]] && echo "ERROR: -name option is mandatory, found '$1' instead" && usage && exit 1
[[ -z "$2" ]] && echo ERROR: no argument found after -name option && usage && exit 1
StreamProcessor="${2}_${JOB_ID}.${CurrentInstance:-1}"
shift ; shift
#
# collect mandatory maxidle option
#
[[ "$1" != -maxidle ]] && echo "ERROR: -maxidle option is mandatory, found '$1' instead" && usage && exit 1
[[ -z "$2" ]] && echo "ERROR: no argument found after -maxidle option" && usage && exit 1
MaxIdle=$2
shift ; shift
#
# collect streams from mandatory queues option
#
[[ "$1" != -queues ]] && echo "ERROR: -queues option is mandatory, found '$1' instead" && usage && exit 1
[[ -z "$2" ]] && echo "ERROR: no argument found after -queues option" && usage && exit 1
shift
for stream in $* ; do
  [[ "${stream}" != *:* || "${stream%%:*}" == ${CurrentInstance} ]] && Streams="${Streams} ${stream##*:}"
done
InstanceEcho "INSTANCE ${CurrentInstance}: BATCH_MPI_CPUS=$BATCH_MPI_CPUS OMP_NUM_THREADS=$OMP_NUM_THREADS StreamProcessor=$StreamProcessor MaxIdle=$MaxIdle Streams=$Streams #Cpus=${MinCpusInJob}->${MaxCpusInJob}"
#
cd ${HOME}/.job_queues                       || exit 1
printf "${Streams}" >${HOME}/.job_queues/.active_${StreamProcessor}.queues || exit 1
#
# launch monitor
#
[[ -x ${0%/*}/u.job-monitor ]] && ${0%/*}/u.job-monitor &
sleep 10  # first sleep cycle
#
IdleSince=$(date +%s)
TimeNow=$(date +%s)
((IdleTime=TimeNow-IdleSince))
while ((IdleTime<MaxIdle))
do
    printf "#CurrentInstance=${CurrentInstance}\nMinCpusInJob=${MinCpusInJob}\nMaxCpusInJob=${MaxCpusInJob}\nMaxIdle=${MaxIdle}\nTotalIdleTime=${IdleTime}\n" \
           > ${HOME}/.job_queues/.active_${StreamProcessor}
    ((ProcessedJobs=0))
    for Stream in $(cat ${HOME}/.job_queues/.active_${StreamProcessor}.queues)   # loop over streams
    do
	[[ -d ${HOME}/.job_queues/${Stream} ]] || continue
	cd ${HOME}/.job_queues/${Stream}       || continue
        touch .active_${StreamProcessor}       || continue
        [[ -d .active_${StreamProcessor}.d ]]  || mkdir .active_${StreamProcessor}.d || continue

	export StreamDir=`pwd -P`
	export StreamFlagFile=${StreamDir}/.active_${StreamProcessor}
	export StreamFlagDir=${StreamFlagFile}.d

	if [[ -r $StreamFlagFile ]]  # try to get a job from stream if stream is active
	then
            [[ -r ${HOME}/.job_queues/.active_${StreamProcessor} ]] || break 10  # stop signal from user
	    TimeNow=$(date +%s)
	    ((TimeLeft=JobStartTime+JobTimeLimit-TimeNow))
	    ((TimeLeft<75)) && \
              InstanceEcho "INFO: less than 75 seconds left, exiting" && \
              break 10   # too little time left, break and cleanup
	    ${0%/*}/pick_and_run_work+  $TimeLeft  && \
              ((ProcessedJobs=ProcessedJobs+1)) && \
              IdleSince=$(date +%s) && \
              break  # got a job, try again to pick from highest queue
	fi

    done

    if((ProcessedJobs==0))  # no job found in any stream, run a filler job
    then
      TimeNow=$(date +%s)
      ((IdleTime=TimeNow-IdleSince))
      ((nsleeps==0)) && InstanceEcho "IDLE=$IdleTime, MAXIDLE=$MaxIdle" 1>&2
      ((TimeLeft=JobStartTime+JobTimeLimit-TimeNow))
#      ((TimeLeft>=90))                && time run_simu_opt30_nompi
      ((TimeLeft>=90))                && sleep 10
      ((TimeLeft>=15 && TimeLeft<90)) && sleep 10
      ((nsleeps=(nsleeps+1)%30))
    fi

    TimeNow=$(date +%s)
    ((IdleTime=TimeNow-IdleSince))

    [[ -r ${HOME}/.job_queues/.active_${StreamProcessor} ]] || break   # stop signal from user if file has been removed
    . ${HOME}/.job_queues/.active_${StreamProcessor}                   # update MaxIdle, MinCpusInJob in case user changed the value
    ((TimeLeft<75)) && InstanceEcho "INFO: less than 75 seconds left, exiting" && break      # too little time left, break and cleanup
done # end of while ((IdleTime<MaxIdle))
#
rm -f ${HOME}/.job_queues/.active_${StreamProcessor}   # cleanup at end of streams
rm -f ${HOME}/.job_queues/.active_${StreamProcessor}.queues
for Stream in $Streams
do
    [[ -d ${HOME}/.job_queues/$Stream ]] || continue
    InstanceEcho ==== JOBS run from stream $Stream ====
    InstanceCat ${HOME}/.job_queues/${Stream}/.active_${StreamProcessor}
    rm -f ${HOME}/.job_queues/${Stream}/.active_${StreamProcessor}
    rmdir ${HOME}/.job_queues/${Stream}/.active_${StreamProcessor}.d
done
rm ${PARALLEL_NODEFILE}
InstanceEcho "INFO: instance no ${CurrentInstance} done, removing ${PARALLEL_NODEFILE}"