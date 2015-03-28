#!/bin/bash
if [[ "$1" == --mod* ]] ; then
  IsModule=YES
  [[ "$1" == *=*  ]] && ModuleName=${1#*=}
  shift
fi
StrSize=""
if [[ "$1" == -L* ]] ; then
  StrSize=${1#-L}
  shift
fi
EnumName=${1}
shift
ModuleName=${ModuleName:-mod_${EnumName}}
[[ -z $1 ]] && echo "EROR: list of names is empty" && exit 1
#
StrSize=${StrSize:-$(for Name in $* ; do echo $Name ; done | awk '{ print length($0); }' | sort -n | tail -1)}
#
[[ -f ${EnumName}.inc ]] && { echo "ERROR: ${EnumName}.inc exists" ; exit 1 ; }
Printf=$(which printf)
#
cat <<EOT >> ${EnumName}.inc
ENUM, BIND(C)
  ENUMERATOR :: ${EnumName}_first=0
EOT
for Name in $*
do
  echo "  ENUMERATOR :: ${Name} " >>${EnumName}.inc
done
cat <<EOT >> ${EnumName}.inc
  ENUMERATOR :: ${EnumName}_last
END ENUM
EOT
if [[ -n ${IsModule} ]] ; then
  rm -f ${EnumName}.F90
  [[ "${ModuleName}" != NoNe ]] && echo "module ${ModuleName}" >>${EnumName}.F90
  echo " include '${EnumName}.inc'" >>${EnumName}.F90
  echo "character(len=${StrSize}), DIMENSION(${EnumName}_first:${EnumName}_last) :: ${EnumName}_table = [  &" >>${EnumName}.F90
  ((names=0))
  ((MaxPerLine=72/(StrSize+4)))
  for Name in "" $*
  do
    ${Printf} " '%-${StrSize}s'," "${Name}" >>${EnumName}.F90
    ((names=names+1))
    if ((names==MaxPerLine)) ; then
      ((names=0))
      echo " &" >>${EnumName}.F90
    fi
  done
  ${Printf} " '%${StrSize}s'" "" >>${EnumName}.F90
  echo "]" >>${EnumName}.F90
  [[ "${ModuleName}" != NoNe ]] && echo "end module ${ModuleName}">>${EnumName}.F90
fi