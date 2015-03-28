#!/bin/bash
if [[ "$1" == --mod* ]] ; then
  Extension=F90
  IsModule=YES
  shift
else
  Extension=inc
fi
EnumName=${1:-liste}
shift
StrSize=${1:-8}
shift
[[ -f ${EnumName}.${Extension} ]] && { echo "ERROR: ${EnumName}.${Extension} exists" ; exit 1 ; }
Printf=$(which printf)
#
[[ -n ${IsModule} ]] && echo "module mod_${EnumName}" >>${EnumName}.${Extension}
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
for Name in "" $*
do
  ${Printf} " '%-${StrSize}s'," "${Name}" >>${EnumName}.${Extension}
  ((names=names+1))
  if ((names==4)) ; then
    ((names=0))
    echo " &" >>${EnumName}.${Extension}
  fi
done
${Printf} " '%${StrSize}s'" "" >>${EnumName}.${Extension}
echo "  ]" >>${EnumName}.${Extension}
[[ -n ${IsModule} ]] && echo "end module mod_${EnumName}" >>${EnumName}.${Extension}
