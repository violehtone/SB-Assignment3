#! /bin/bash
# calculates the center of mass of the backbone atoms of a given CG protein structfile (.gro format)

structfile=$1
lstructfile=`wc -l $structfile|awk '{print $1}'`
ncutlines_bot=1
ncutlines_top=2
tmpstruct="tmpstruct.gro"
head -$[$lstructfile-$ncutlines_bot] $structfile | tail -$[$lstructfile-$ncutlines_bot-$ncutlines_top]> $tmpstruct
indexfile=$2
xyzfile="com.xyz"
distfile="com.dist"

lines=(`grep -n "\[" $indexfile | grep "backbone" | awk '{print $1}' | cut -d\: -f1`)
llines=${#lines[@]}
names=(`grep -n "\[" $indexfile | grep "backbone" | awk '{print $2}'`)
lnames=${#names[@]}
lastline=`wc -l $indexfile | awk '{print $1}'`
if [ $llines -eq $lnames ]; then
  echo $llines > $xyzfile
  echo "Generated by com.sh on" `date` "NOTE: COORDINATES IN ANGSTROMS!!">> $xyzfile
  for ((i=0;i<$llines;i++))
  do
    start[${i}]=$[${lines[$i]}+1]
    flaglines=`grep -n "\[" $indexfile | cut -d\: -f1`
    flaglines="$flaglines $lastline"
    for j in $flaglines
    do
      if [ $j -gt ${lines[$i]} ];then
        if [ $j -lt $lastline ];then
          fin[$i]=$[$j-1]
        else
          fin[$i]=$j
        fi
        break
      fi
    done
    BNlist=`head -${fin[$i]} $indexfile | tail -$[${fin[$i]}-${start[$i]}+1]`
    for j in $BNlist
    do
      sed -n ${j}p $tmpstruct
    done |\
    awk '{n++;comx+=$4;comy+=$5;comz+=$6;print comx,comy,comz,n}' |\
    tail -1 |\
    awk '{print "N",10*$1/$4,10*$2/$4,10*$3/$4}' >> $xyzfile
  done
fi

echo "Distances between COMs:" > $distfile
echo "=======================" >> $distfile
comx=(`tail -$llines $xyzfile | awk '{print $2}'`)
comy=(`tail -$llines $xyzfile | awk '{print $3}'`)
comz=(`tail -$llines $xyzfile | awk '{print $4}'`)
length=${#comz[@]}
for((i=0;i<$[$length-1];i++))
do
  for((j=$[$i+1];j<$length;j++))
  do
    dx2=`echo ${comx[$i]} ${comx[$j]} | awk '{print ($2-$1)^2}'`
    dy2=`echo ${comy[$i]} ${comy[$j]} | awk '{print ($2-$1)^2}'`
    dz2=`echo ${comz[$i]} ${comz[$j]} | awk '{print ($2-$1)^2}'`
    dr=`echo $dx2 $dy2 $dz2 | awk '{print sqrt($1+$2+$3)*0.1}'`
    echo "COM distance between groups" "\""${names[$i]}"\"" "and" "\""${names[$j]}"\"" "=" $dr "nm" >> $distfile
  done
done

rm $tmpstruct

exit