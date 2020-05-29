#! /bin/bash

l=$1
tclfile="box.tcl"

echo "draw line {0 0 0} {$l 0 0}"  > $tclfile
echo "draw line {$l 0 0} {$l $l 0}" >> $tclfile
echo "draw line {$l $l 0} {0 $l 0}" >> $tclfile
echo "draw line {0 $l 0} {0 0 0}" >> $tclfile

echo "draw line {0 0 0} {0 0 $l}"  >> $tclfile
echo "draw line {$l 0 0} {$l 0 $l}" >> $tclfile
echo "draw line {$l $l 0} {$l $l $l}" >> $tclfile
echo "draw line {0 $l 0} {0 $l $l}" >> $tclfile

echo "draw line {0 0 $l} {$l 0 $l}"  >> $tclfile
echo "draw line {$l 0 $l} {$l $l $l}" >> $tclfile
echo "draw line {$l $l $l} {0 $l $l}" >> $tclfile
echo "draw line {0 $l $l} {0 0 $l}" >> $tclfile

exit