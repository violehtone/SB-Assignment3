#! /bin/bash

itpfile=$1

nstart=`grep -n "atoms" $itpfile | grep "\[" | grep "\]" | cut -d\: -f1`
taglinelist=`grep -n "\[" $itpfile | grep "\]" | cut -d\: -f1`
for i in $taglinelist
do
  if [ $i -gt $nstart ]; then
    nfin=$i
    break
  fi
done
nstart=$[$nstart+1]
for i in 0 1 2 3 4 5
do
  line=`sed -n "$[$nfin-$i]"p $itpfile`
  if [ "`echo $line|awk '{print $1}'| grep "[0-9]"`" != "" ];then
    ncutlines=$i
    break
  fi
done
nfin=$[$nfin-$ncutlines]

qtot=`head -$nfin $itpfile | tail -$[$nfin-$nstart+1] | awk '{sum+=$7;print sum}' | tail -1`
echo "Total charge in \"$itpfile\": $qtot"

exit