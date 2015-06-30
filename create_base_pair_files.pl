#!/usr/bin/perl
#
#

$nBases = 15;
#$nBases = 3;

@baseAtomNames = qw(N9 C8 H8 N7 C5 C6 N6 H61 H62 N1 C2 H2 N3 C4 N2 H21 H22);

for ($base1=1;$base1<$nBases;$base1++) {

	for ($base2=$base1+1;$base2<=$nBases;$base2++) {
		if ($base1==1 && $base2==2) {
		} else {	
			system("mkdir base$base1\_base$base2");
			system("sed -e s/base1_base2/base$base1\_base$base2/g < base1_base2/pairenergy.cfg > base$base1\_base$base2/pairenergy.cfg");
		}	
			# now we change pdb file selection
			open(INP,"grep ATOM ../nucleic_ions.pdb|");
			open(OUT,">base$base1\_base$base2/base$base1\_base$base2.pdb");
		
			while($line=<INP>) {
				chomp $line;
				@array = split ' ', $line;
				$atom = $array[0];
				$atomNumber = $array[1];
				$atomName = $array[2];
				$resName = $array[3];
				$chain = $array[4];
				$resid = $array[5];
				$x = $array[6];
				$y = $array[7];
				$z = $array[8];
				if ($resid == $base1) {
					$match = "false";
					$i=0;
					while ($match eq "false" && $i<@baseAtomNames) {
						if ($baseAtomNames[$i] eq $atomName) {
							$match = "true";
						}
						$i++;
					}
					if ($match eq "true") {
						$beta = 1;
					} else {
						$beta = 0;
					}
				} elsif ($resid == $base2) {	
					$match = "false";
					$i=0;
					while ($match eq "false" && $i<@baseAtomNames) {
						if ($baseAtomNames[$i] eq $atomName) {
							$match = "true";
						}
						$i++;
					}
					if ($match eq "true") {
						$beta = 2;
					} else {
						$beta = 0;
					}
				} else {
					$beta = 0;
				}
				printf OUT "%4s %6d %4s %3s %1s %3d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", $atom, $atomNumber, $atomName, $resName, $chain, $resid, $x, $y, $z, 0.00, $beta;
			}

			close(INP);
			close(OUT);


			# run namd
			system("cd base$base1\_base$base2;namd2 +p2 pairenergy.cfg > base$base1\_base$base2.log; ../grab_namd_energies.pl base$base1\_base$base2.log base$base1\_base$base2.energies.dat");
	}
}

system("./create_combined_vdw_dat.pl")



