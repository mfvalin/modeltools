alias shortcut='. s.ssmuse.dot'
_shortcut_compfn()
{
 local cur prev opts cmd
 cmd="${1##*/}"
 COMPREPLY=()
 cur="${COMP_WORDS[COMP_CWORD]}"
 COMPREPLY=( $(compgen -W "$(_complete_shortcut ${cur} )" ) )
 [[ $BASH_VERSION == 4* ]] && compopt -o nospace
 return 0
}
complete -F _shortcut_compfn shortcut

_echo_dir()
{
  _item=${1}
  [[ -d ${_item} ]] && echo ${_item}/
  [[ -f ${_item} ]] && echo ${_item}
}

_shortcut_possibilities()
{
for j in $( echo ${SSM_SHORTCUT_PATH} ${HOME}/my_ssm_domains  $(r.unified_setup-cfg -local || echo $ARMNLIB)/data/ssm_domains $MODULEPATH | tr ':' ' ' )
do
  if [[ -d $j ]] ; then
    cd $j
    for i in ${1}* ; do
     _echo_dir "$i"
     done
  fi
done
}

_complete_shortcut()
{
Liste="$(_shortcut_possibilities "${1:-[a-zA-Z0-9]}" | sed -e 's/[.]sh$//g' -e 's/[.]bndl$//' -e 's://:/:' | grep -v '[*]$' | sort -u)"
echo ${Liste:-${1:-[a-zA-Z0-9]}}
}
