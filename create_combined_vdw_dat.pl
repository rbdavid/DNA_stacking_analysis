#!/usr/bin/perl
#
#

$nBases = 15;

$timeSteps = `wc -l base1_base2/base1_base2.energies.dat`;
print $timeSteps;

for ($i=0;$i<$timeSteps;$i++) {
	for ($base1=1;$base1<=$nBases;$base1++) {
		$vdwE[$i][$base1] = 0;
	}
}


for ($base1=1;$base1<$nBases;$base1++) {

	for ($base2=$base1+1;$base2<=$nBases;$base2++) {

		open(INP,"<base$base1\_base$base2/base$base1\_base$base2.energies.dat");
		$step=0;
		while ($line=<INP>) {
			chomp $line;
			@array = split ' ', $line;
			$time[$step] = $array[0];
			$vdwE[$step][$base1] += $array[2];
			$vdwE[$step][$base2] += $array[2];
			$step++;
		}


		close(INP);


	}
}

open(OUT,">combined_vdw.dat");
for ($i=0;$i<$timeSteps;$i++) {
	printf OUT "%10.5f ", $time[$i];
        for ($base1=1;$base1<=$nBases;$base1++) {
                printf OUT "%10.5f", $vdwE[$i][$base1];
        }
	print OUT "\n";
}
close(OUT);



