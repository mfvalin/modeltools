#!/bin/bash 
# intercept ld call, add rpath argument to make sure all needed shared libraries will be found at runtime
allParams=("$@")     # collect arguments to ld
# for i in "${allParams[@]}" ; do
#   echo $i
# done
# 
TRUE_LD=/usr/bin/ld
n=0
nParams=${#allParams[*]}
static=0
LibPath=()
LibList=()
RpathList='-rpath $ORIGIN/../lib -rpath $ORIGIN/../lib64'
unset PATHUSED
declare -A PATHUSED
PATHUSED["/usr/lib"]=1         # ignore /usr/lib and /usr/lib64 as they are part of default rpath
PATHUSED["/usr/lib64"]=1
while ((n < nParams)); do
    p=${allParams[n]}
    p2=${allParams[$((n+1))]}
    if [[ "${p:0:3}" == -L/ ]]; then
	LibPath+=("${p:2}")
    elif [[ "$p" == -L ]]; then
	LibPath+=("${p2}")
	((n = n + 1))
    elif [[ "$p" == -l && ${static} == 0 ]]; then
	LibList+=("${p2}")
	((n = n + 1))
    elif [[ "${p:0:2}" == -l && ${static} == 0 ]]; then
	LibList+=("${p:2}")
    elif [[ "$p" == -dynamic-linker ]]; then     # ignore library name that follows -dynamic-linker
	((n = n + 1))
    elif [[ "$p" == -Bstatic ]]; then            # static mode, ignore library names that follow
	static=1
    elif [[ "$p" == -static ]]; then             # static mode, ignore library names that follow
	static=1
    elif [[ "$p" == -Bdynamic ]]; then           # dynamic mode, collect library names that follow
	static=0
    elif [[ "$p" == -rpath ]]; then              # explicit rpath, use canonical path if absolute path used
        path="${p2}"
        [[ ${p2} == /* ]] && path="$(readlink -f ${p2})" && allParams[$((n+1))]="${path}"
        PATHUSED["${path}"]=1 
        [[ -n ${LD_WRAPPER_DEBUG} ]] && echo "EXPLICIT rpath = ${path} (${p2})"
    elif [[ "$p" =~ ^[^-].*\.so($|\.) ]]; then
	# direct reference to a shared library, so add its canonical path to rpath.
	path="$(readlink -e $(dirname "$p"))";
	PATHUSED["${path}"]=1
	RpathList="${RpathList} -rpath ${path}"
    fi
    ((n = n + 1))
done
[[ -n ${LD_WRAPPER_DEBUG} ]] &&  echo LibPath=${LibPath[@]}
[[ -n ${LD_WRAPPER_DEBUG} ]] &&  echo LibList=${LibList[@]}
#
unset LIBS_SEEN
declare -A LIBS_SEEN
for i in ${LibPath[@]}; do
  path="$(readlink -f ${i})"
  [[ -n ${LD_WRAPPER_DEBUG} ]] &&    [[ -n ${PATHUSED["$path"]} ]] && echo "<<< IGNORING $path >>>" && continue
  [[ -n ${PATHUSED["$path"]} ]] && continue
  for j in ${LibList[@]}; do
      foundlib=${LIBS_SEEN["$j"]}
      if [ -z "$foundlib" -a -f "$i/lib$j.so" ]; then
	  RpathList="${RpathList} -rpath ${path}"
	  LIBS_SEEN["$j"]=1
	  break
      fi
  done
  PATHUSED["${path}"]=1
done
[[ -n ${LD_WRAPPER_DEBUG} ]] &&  echo RpathList=${RpathList[@]}
[[ -n ${LD_WRAPPER_DEBUG} ]] &&  set -x
exec ${TRUE_LD} "${allParams[@]}" ${RpathList[@]}   # and now call the effective ld