#!/usr/bin/perl -w

if (@ARGV != 4) {
     print STDERR "Wrong number of arguments.
     Usage: $0 input lower_cutoff higher_cutoff force_constant\n\n";
     exit(0);
}

open(FILE, $ARGV[0]) || die "Can't open $ARGV[0]: $!";
$cutoff1 = $ARGV[1];
$cutoff2 = $ARGV[2];
$forcek = $ARGV[3];

while ($line = <FILE>) {
    if ($line =~ /(\d+)\s+(\d+)\s+(\d+)\s+([0-9.e+-]+)/ ) {
       $atom1 = $1;
       $atom2 = $2;
       $type  = $3;
       $dist = $4;
       if ($dist < $cutoff2 && $dist > $cutoff1) {
          printf "%d \t %d \t %d \t %.4f \t %.0f\n", $atom1,$atom2,$type,$dist,$forcek;
       }
    }
}

close(FILE);


