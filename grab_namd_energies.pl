#!/usr/bin/perl
#



$inp = $ARGV[0];
$out = $ARGV[1];

$deltaT = 2E-6;

open(INP,"grep 'ENERGY: ' $inp|");
open(OUT,">$out");
	$count = 0;
	while ($line = <INP>) {
		
		if ($count > 1) {
			chomp $line;
			@array = split ' ', $line;
			$step = $array[1];
			$elecE = $array[6];
			$vdwE = $array[7];
			$potE = $array[13];
			printf OUT "%10.5f %10.5f %10.5f %10.5f\n", ($step+1)*$deltaT, $elecE, $vdwE, $potE; 
		}
		$count++

	}

close(OUT);
close(INP);


