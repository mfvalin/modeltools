#!/bin/bash
if [[ "$1" == --mod* ]] ; then
  Extension=F90
  IsModule=YES
  [[ "$1" == *=*  ]] && ModuleName=${1#*=}
  shift
else
  Extension=inc
fi
StrSize=""
if [[ "$1" == -L* ]] ; then
  StrSize=${1#-L}
  shift
fi
EnumName=${1:-liste}
shift
[[ -z $1 ]] && echo "EROR: list of names is empty" && exit 1
#
StrSize=${StrSize:-$(for Name in $* ; do echo $Name ; done | awk '{ print length($0); }' | sort -n | tail -1)}
#
[[ -f ${EnumName}.${Extension} ]] && { echo "ERROR: ${EnumName}.${Extension} exists" ; exit 1 ; }
Printf=$(which printf)
#
ModuleName=${ModuleName:-mod_${EnumName}}
[[ -n ${IsModule} ]] && echo "module ${ModuleName}" >>${EnumName}.${Extension}
cat <<EOT >> ${EnumName}.${Extension}
ENUM, BIND(C)
  ENUMERATOR :: ${EnumName}_first=0
EOT
for Name in $*
do
  echo "  ENUMERATOR :: ${Name} " >>${EnumName}.${Extension}
done
cat <<EOT >> ${EnumName}.${Extension}
  ENUMERATOR :: ${EnumName}_last
END ENUM
character(len=${StrSize}), DIMENSION(${EnumName}_first:${EnumName}_last) :: ${EnumName}_table = [  &
EOT
((names=0))
((MaxPerLine=72/(StrSize+4)))
for Name in "" $*
do
  ${Printf} " '%-${StrSize}s'," "${Name}" >>${EnumName}.${Extension}
  ((names=names+1))
  if ((names==MaxPerLine)) ; then
    ((names=0))
    echo " &" >>${EnumName}.${Extension}
  fi
done
${Printf} " '%${StrSize}s'" "" >>${EnumName}.${Extension}
echo "  ]" >>${EnumName}.${Extension}
[[ -n ${IsModule} ]] && echo "end module ${ModuleName}" >>${EnumName}.${Extension}
