#!/usr/bin/perl
############################### This program is free software; you can redistribute it and/or modify
############################### it under the terms of the GNU General Public License as published by
############################### the Free Software Foundation; either version 2 of the License, or
############################### (at your option) any later version.
###############################
############################### This program is distributed in the hope that it will be useful,
############################### but WITHOUT ANY WARRANTY; without even the implied warranty of
############################### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
############################### GNU General Public License for more details.
###############################
############################### You should have received a copy of the GNU General Public License
############################### along with this program; if not, write to the Free Software
############################### Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
###############################
###############################
###############################
############################### Copyright 2019 Damon P. Little
############################### MODULES
use strict;
############################### PARSE OPTIONS
my $ax = '20';
my $d = '';
my $f = 0;
my $i = 1;
my $g = 1;
my $m = 100;
my $n = '';
my %o = {};
my $r = 0;
for(my $k = $#ARGV; $k >= 0; $k--){
	if($ARGV[$k] eq '-a'){
		if($ARGV[$k+1] =~ m/^(2|3|4|5|6|6dso|6kgb|6sr|8|10|11|12|15|18|20)$/){
			$ax = $ARGV[$k+1];
		}
	} elsif($ARGV[$k] eq '-d'){
		if(-e $ARGV[$k+1]){
			$d = $ARGV[$k+1];
			$d =~ s*/$**;
		} else{
			print(STDERR "\nInput directory (-d) could not be read!\n");
		}
	} elsif($ARGV[$k] eq '-f'){
		$f = 1;
	} elsif($ARGV[$k] eq '-i'){
		$i = 0;
	} elsif($ARGV[$k] eq '-g'){
		$g = 0;
	} elsif($ARGV[$k] eq '-m'){
		$ARGV[$k+1] =~ tr/0123456789//cd;
		if(($ARGV[$k+1] > 0) && ($ARGV[$k+1] <= 100)){
			$m = $ARGV[$k+1]
		}
	} elsif($ARGV[$k] eq '-n'){
		$n = $ARGV[$k+1];
	} elsif($ARGV[$k] eq '-o'){
		my @x = split(/,/, $ARGV[$k+1]);
		for(my $j = $#x; $j >= 0; $j--){
			$o{cleanName($x[$j])} = 1;
		}
	} elsif($ARGV[$k] eq '-r'){
		$r = 1;
	}
}
$m /= 100;
my @infile = ();
if(length($d)){
	opendir(INDIR, $d) or die("Cannot open directory: $d");
	@infile = sort({$b cmp $a} readdir(INDIR));
	closedir(INDIR);
	if($#infile == -1){
		print(STDERR "Input directory (-d) does not contain any files!\n");
	}
}
############################### START
if(($#infile != -1) && length($n)){
	my @files = ('tnt', 'cna', 'phy', 'par');
	my $s = "find $d/ -type f | xargs grep -h '>' | tr ' ' '_' | tr '/' '_' | tr '\\\\' '_' | tr '-' '_' ";
	if($f == 0){
		$s .= "| perl -F'#' -lane '{print(\$F[0])}'"
	}
	$s .=	"| perl -lane '{\$F[0]=~tr/[A-z][0-9]\\_\\.//cd;print(\$F[0])}' | sort -u -r";
	my @species = split('\n', qx`$s`);
	my @speciesFNV = ();
	my $max = 0;
	for(my $k = $#species; $k >= 0; $k--){
		$species[$k] = cleanName($species[$k]);
		$speciesFNV[$k] = fnv32a($species[$k]);
		push(@files, $speciesFNV[$k]);
		if(length($species[$k]) > $max){
			$max = length($species[$k]);
		}
	}
	my @matrixSpeciesTNT = ();
	my @matrixSpeciesPHY = ();
	for(my $k = $#species; $k >= 0; $k--){
		if(exists($o{$species[$k]})){
			$matrixSpeciesTNT[$k] = '@OUTGROUP_' . $species[$k] . (' ' x ($max-length($species[$k])+1));
		} else {
			$matrixSpeciesTNT[$k] = '@INGROUP_' . $species[$k] . (' ' x ($max-length($species[$k])+2));
		}
		$matrixSpeciesPHY[$k] = $species[$k] . (' ' x ($max-length($species[$k])+2));
	}
	my @terminals = ();
	for(my $x = $#infile; $x >= 0; $x--){ ### not done as a subshell because file number restrictions: my $p = qx`find $i/ -type f | xargs grep -c -h '>' | awk '{if(\$1>3){print $#species+1-\$1}}' | sort -n | head -n \$(find $d/ -type f | xargs grep -c -h '>' | awk 'BEGIN{x=0}{if(\$1>3){x++}}END{print int(x*$m)}') | tail -1`; $p = $#species+1-$p;
		if($infile[$x] !~ m/^\.+$/){
			my $t = qx`grep -c '>' $d/$infile[$x]`;
			$t =~ tr/0123456789//cd;
			if($t > 3){
				push(@terminals, $t);
			}
		}
	}
	@terminals = sort({$a <=> $b} @terminals);
	my $p = $terminals[int((1-$m)*($#terminals+1))];
	my %files = {};
	my %buffers = {};
	for(my $k = $#files; $k >= 0; $k--){
		open($files{$files[$k]}, '>' . $files[$k] . '-temporary-do-not-touch') or die("Could not open $files[$k]-temporary-do-not-touch!");
		$buffers{$files[$k]} = '';
	}
	$buffers{'cna'} = ";\ncc - .;\n";
	if(scalar(keys(%o)) > 1){
		$buffers{'cna'} .= "outgroup[OUTGROUP;\n";
	}
	$buffers{'cna'} .= "proc/;\n\n\n#\n\$\n;\ncn ";
	my $AA = 0;
	my $first = 1;
	my $ncharTNT = 0;
	my $ncharFAS = 0;
	my $ncharPHY = 0;
	for(my $x = $#infile; $x >= 0; $x--){ ############################### READ INFILE(S) PARSE AND WRITE TO TEMPORARY FILES
		if(($infile[$x] =~ m/^\.{1,2}$/) || ($infile[$x] =~ m/^\.DS_Store$/)){
			next;
		}
		my $aa = 0;
		my $file = $infile[$x];
		$file =~ tr/[A-z][0-9]\.//cd;
		$file = substr($file, 0, index($file, '.'));
		my @indelNames = ();
		my $matrix = {}; ### ->{terminal}->[sequence]
		my $sequence = '';
		my $species = '';
		open(INFILE, '<', $d . '/' .$infile[$x]) or die("Could not open $infile[$x]!\n");
		while(my $line = <INFILE>){
			$line =~ s/[\r\n]+$//;
			if(length($line)){
				if($line =~ m/^>(?:(.+)#(.+)|([^#]+))$/){
					if(length($species) && length($sequence)){
						if(insert($matrix, $sequence, $species, $ax)){
							$aa = 1;
						}
					}
					if($f == 0){
						$species = $1;
					} else {
						$species = $3
					}
					$sequence = '';
				} else {
					$sequence .= $line;
				}
			}
		}
		if(length($species) && length($sequence)){
			if(insert($matrix, $sequence, $species, $ax)){
				$aa = 1;
			}
		} else {
			print(STDERR "input file '$infile[$x]' skipped...\n");
			next;
		}
		if($aa == 0){
			$ax = '20';
		}
		close(INFILE);
		my $y = scalar(keys(%{$matrix}));
		if(scalar(keys(%{$matrix})) < $p){
			print(STDERR "input file '$infile[$x]' skipped because it has too much missing data (m = $m; p = $p; t = $y)...\n");
			next;
		}
		my $sequenceCharacters = $#{$matrix->{$species}};
		for(my $k = $#species; $k >= 0; $k--){
			if(exists($matrix->{$species[$k]})){
				if($#{$matrix->{$species[$k]}} != $sequenceCharacters){
					die("Sequences must aligned. Error in $infile[$x]: $species[$k] = $#{$matrix->{$species[$k]}} and matrix = $sequenceCharacters!\n");
				}
			}
		}
		if($i == 1){
			my $ends->[0][0] = ();
			my %indels = ();
			for(my $k = $#species; $k >= 0; $k--){
				if(exists($matrix->{$species[$k]})){
					my $lastFive = 0;
					my $lastThree = 0;
					my $inGap = 0;
					for(my $j = $sequenceCharacters; $j >= 0; $j--){ ### FIND INDELS
						if(($matrix->{$species[$k]}->[$j] eq '-') && ($inGap == 0)){ ### gap 3'
							$lastThree = $j;
							$inGap = 1;
						}
						if(($matrix->{$species[$k]}->[$j] ne '-') && ($inGap == 1)){ ### gap 5'
							$inGap = 0;
							$lastFive = $j+1;
							if(($lastFive != 0) && ($lastThree != $sequenceCharacters)){ ### gaps not at ends
								$indels{$lastFive . '-' . $lastThree} = $lastFive;
							}
						}
					}
					for(my $j = 0; $j <= $sequenceCharacters; $j++){ ### FIND ENDS
						if($matrix->{$species[$k]}->[$j] =~ m/A|B|C|D|E|F|G|H|I|K|L|M|N|O|P|Q|R|S|T|U|V|W|X|Y|Z/){
							$ends->[$k][0] = $j;
							last;
						}
					}
					for(my $j = $sequenceCharacters; $j >= 0; $j--){ ### FIND ENDS
						if($matrix->{$species[$k]}->[$j] =~ m/A|B|C|D|E|F|G|H|I|K|L|M|N|O|P|Q|R|S|T|U|V|W|X|Y|Z/){
							$ends->[$k][1] = $j;
							last;
						}
					}
				}
			}
			if(scalar(keys(%indels))){
				@indelNames = sort({$indels{$a} <=> $indels{$b}} keys(%indels));
			}
			undef(%indels);
			for(my $g = $#indelNames; $g >= 0; $g--){ ### SCORE INDELS
				(my $fiveEnd, my $threeEnd) = split(/-/, $indelNames[$g]);
				for(my $k = $#species; $k >= 0; $k--){
					if(exists($matrix->{$species[$k]})){
						my $score = '?';
						if(($fiveEnd > $ends->[$k][0]) && ($threeEnd < $ends->[$k][1])){
							for(my $j = $fiveEnd; $j <= $threeEnd; $j++){
								if($matrix->{$species[$k]}->[$j] =~ m/A|B|C|D|E|F|G|H|I|K|L|M|N|O|P|Q|R|S|T|U|V|W|X|Y|Z/){
									$score = 0;
									last;
								}
							}
							if($score eq '?'){
								if(($matrix->{$species[$k]}->[$fiveEnd-1] ne '-') && ($matrix->{$species[$k]}->[$threeEnd+1] ne '-')){
									$score = 1;
								} else {
									$score = '-';
								}
							}
						}
						$matrix->{$species[$k]}->[$sequenceCharacters+$g+1] = $score;
					}
				}
			}
		}
		my @informative = ();
		my $informativeSequence = 0;
		my $informativeIndel = 0;
		for(my $j = $#{$matrix->{$species}}; $j >= 0; $j--){ ### FIND INFORMATIVE
			my %score = ();
			my $states = 0;
			for(my $k = $#species; $k >= 0; $k--){
				if(exists($matrix->{$species[$k]})){
					if(($matrix->{$species[$k]}->[$j] ne '?') && ($matrix->{$species[$k]}->[$j] ne '-')){
						if(!exists($score{$matrix->{$species[$k]}->[$j]})){
							$score{$matrix->{$species[$k]}->[$j]} = 0;
							$states++;
						}
						$score{$matrix->{$species[$k]}->[$j]} += 1;
					}
				}
			}
			if($states > 1){
				my $minSteps = $states-1;
				my $maxSteps = 0;
				my @steps = sort({$score{$b} <=> $score{$a}} keys(%score));
				for(my $i = $#steps; $i > 0; $i--){
					$maxSteps += $score{$steps[$i]};
				}
				if(($maxSteps-$minSteps) > 0){
					$informative[$j] = 1;
					if($j <= $sequenceCharacters){
						$informativeSequence++;
					} else {
						$informativeIndel++;
					}
				} else {
					$informative[$j] = 0;
				}
			} else {
				$informative[$j] = 0;
			}
		}
		if($informativeSequence){
			if($aa == 1){
				if($ax =~ m/^(2|3|4|5|6|6dso|6kgb|6sr|8|10|11|12|15)$/){
					$buffers{'tnt'} .= "&[num]\n";
				} else {
					$buffers{'tnt'} .= "&[prot nogaps]\n";
				}
			} else {
				$buffers{'tnt'} .= "&[dna nogaps]\n";
			}
		}
		if($first == 0){
			$buffers{'phy'} .= "\n";
		}
		my @gap = (0, 0);
		for(my $k = $#species; $k >= 0; $k--){ ### PRINT SEQUENCE DATA FOR TNT, FAS, AND PHY; CALCULATE GENE ABSENCE/PRESENCE (GAP) INFORMATION
			if($informativeSequence){
				$buffers{'tnt'} .= $matrixSpeciesTNT[$k];
				print({$files{'tnt'}} $buffers{'tnt'});
				$buffers{'tnt'} = '';
			}
			if($first == 1){
				$buffers{$speciesFNV[$k]} .= '>' . $species[$k];
				$buffers{'phy'} .= $matrixSpeciesPHY[$k];
			}
			print({$files{'phy'}} $buffers{'phy'});
			$buffers{'phy'} = '';
			if(exists($matrix->{$species[$k]})){
				for(my $j = 0; $j <= $sequenceCharacters; $j++){
					if($informative[$j] == 1){
						$buffers{'tnt'} .= $matrix->{$species[$k]}->[$j];
					}
					if(($j % 80) == 0){
						$buffers{$speciesFNV[$k]} .= "\n";
					}
					$buffers{$speciesFNV[$k]} .= $matrix->{$species[$k]}->[$j];
					$buffers{'phy'} .= $matrix->{$species[$k]}->[$j];
				}
				$gap[1] += 1;
			} else {
				if($informativeSequence){
					$buffers{'tnt'} .= ('?' x $informativeSequence);
				}
				for(my $j = 0; $j <= $sequenceCharacters; $j++){
					if(($j % 80) == 0){
						$buffers{$speciesFNV[$k]} .= "\n";
					}
					$buffers{$speciesFNV[$k]} .= '-';
				}
				$buffers{'phy'} .= ('-' x ($sequenceCharacters+1));
				$gap[0] += 1;
			}
			if($informativeSequence){
				$buffers{'tnt'} = char2num($buffers{'tnt'}, $ax, 0) . "\n";
			}
			$buffers{$speciesFNV[$k]} .= "\n";
			$buffers{'phy'} = char2num($buffers{'phy'}, $ax, 1) . "\n";
		}
		my $informativeGAP = 0;
		if(($gap[0] > 1) && ($gap[1] > 1)){
			if($gap[0] < $gap[1]){
				$informativeGAP = $gap[0]-1;
			} else {
				$informativeGAP = $gap[1]-1;
			}
		}
		if($first == 1){
			$first = 0;
		}
		my $stateNames = " A C G T /;\n";
		if($aa == 1){
			if($ax eq '20'){
				$stateNames = " A C D E F G H I K L M N O P Q R S T U V W Y /;\n";
			} elsif($ax eq '2'){
				$stateNames = " LVIMCAGSTPFYW EDNQKRH /;\n";
			} elsif($ax eq '3'){
				$stateNames = " LASGVTIPMC EKRDNQH FYW /;\n";
			} elsif($ax eq '4'){
				$stateNames = " LVIMC AGSTP FYW EDNQKRH /;\n";
			} elsif($ax eq '5'){
				$stateNames = " LVIMC ASGTP FYW EDNQ KRH /;\n";
			} elsif($ax eq '6'){
				$stateNames = " LVIM ASGT PHC FYW EDNQ KR /;\n";
			} elsif($ax eq '6dso'){
				$stateNames = " AGPST DENQ HKR MIVL WFY C /;\n";
			} elsif($ax eq '6kgb'){
				$stateNames = " AGPS DENQHKRT MIL W FY CV /;\n";
			} elsif($ax eq '6sr'){
				$stateNames = " APST DENG QKR MIVL WC FYH /;\n";
			} elsif($ax eq '8'){
				$stateNames = " LVIMC AG ST P FYW EDNQ KR H /;\n";
			} elsif($ax eq '10'){
				$stateNames = " LVIM C A G ST P FYW EDNQ KR H /;\n";
			} elsif($ax eq '11'){
				$stateNames = " KREDQN C G H ILV M F Y W P STA /;\n";
			} elsif($ax eq '12'){
				$stateNames = " LVIM C A G ST P FY W EQ DN KR H /;\n";
			} elsif($ax eq '15'){
				$stateNames = " LVIM C A G S T P FY W E D N Q KR H /;\n";
			} elsif($ax eq '18'){
				$stateNames = " LM VI C A G S T P F Y W E D N Q K R H /;\n";
			}
		}
		for(my $j = 0; $j <= $sequenceCharacters; $j++){ ### PRINT SEQUENCE PARTITION DATA FOR CNA
			if($informative[$j] == 1){
				$buffers{'cna'} .= '{' . $ncharTNT . ' ' . $file . '_sequence_' . ($j+1) . $stateNames;
				$ncharTNT++;
			}
		}
		if($aa == 0){ ### PRINT SEQUENCE PARTITION DATA FOR PAR
			if($r == 1){
				$buffers{'par'} .= 'DNAX';
			} else {
				$buffers{'par'} .= 'GTR+G';
			}
		} else {
			if($r == 1){
				if($ax eq '20'){
					$buffers{'par'} .= 'GTR';
				} elsif($ax eq '2'){
					$buffers{'par'} .= 'BIN';
				} else {
					$buffers{'par'} .= 'MULTI';
				}
			} else {
				if($ax eq '20'){
					$buffers{'par'} .= 'PROTGTR+G';
				} elsif($ax eq '2'){
					$buffers{'par'} .= 'BIN+G';
				} elsif($ax eq '3'){
					$buffers{'par'} .= 'MULTI3_GTR+G';
				} elsif($ax eq '4'){
					$buffers{'par'} .= 'MULTI4_GTR+G';
				} elsif($ax eq '5'){
					$buffers{'par'} .= 'MULTI5_GTR+G';
				} elsif(($ax eq '6') || ($ax eq '6dso') || ($ax eq '6kgb') || ($ax eq '6sr')){
					$buffers{'par'} .= 'MULTI6_GTR+G';
				} elsif($ax eq '8'){
					$buffers{'par'} .= 'MULTI8_GTR+G';
				} elsif($ax eq '10'){
					$buffers{'par'} .= 'MULTI10_GTR+G';
				} elsif($ax eq '11'){
					$buffers{'par'} .= 'MULTI11_GTR+G';
				} elsif($ax eq '12'){
					$buffers{'par'} .= 'MULTI12_GTR+G';
				} elsif($ax eq '15'){
					$buffers{'par'} .= 'MULTI15_GTR+G';
				} elsif($ax eq '18'){
					$buffers{'par'} .= 'MULTI18_GTR+G';
				}
			}
		}
		$buffers{'par'} .= ', ' . $file . '_sequence = ' . ($ncharPHY+1) . '-' . ($ncharPHY+$sequenceCharacters+1) . "\n";
		$ncharFAS += $sequenceCharacters+1;
		$ncharPHY += $sequenceCharacters+1;
		if($i == 1){
			if($informativeIndel){
				$buffers{'tnt'} .= "&[num]\n";
			}
			if($#indelNames >= 0){
				$buffers{'phy'} .= "\n";
			}
			for(my $k = $#species; $k >= 0; $k--){ ### PRINT INDEL DATA FOR TNT AND PHY
				if($informativeIndel){
					$buffers{'tnt'} .= $matrixSpeciesTNT[$k];
				}
				if(exists($matrix->{$species[$k]})){
					for(my $j = 0; $j <= $#indelNames; $j++){
						if($informative[$j+$sequenceCharacters+1] == 1){
							$buffers{'tnt'} .= $matrix->{$species[$k]}->[$j+$sequenceCharacters+1];
						}
						$buffers{'phy'} .= $matrix->{$species[$k]}->[$j+$sequenceCharacters+1];
					}
				} else {
					if($informativeIndel){
						$buffers{'tnt'} .= ('?' x $informativeIndel);
					}
					if($#indelNames >= 0){
						$buffers{'phy'} .= ('?' x ($#indelNames+1));
					}
				}
				if($informativeIndel){
					$buffers{'tnt'} .= "\n";
				}
				if($#indelNames >= 0){
					$buffers{'phy'} .= "\n";
				}
			}
			for(my $g = 0; $g <= $#indelNames; $g++){ ### PRINT INDEL PARTITION DATA FOR CNA
				if($informative[$g+$sequenceCharacters+1] == 1){
					(my $fiveEnd, my $threeEnd) = split(/-/, $indelNames[$g]);
					$buffers{'cna'} .= '{' . $ncharTNT . ' ' . $file . '_indel_' . ($fiveEnd+1) . '—' . ($threeEnd+1) . " absent present /;\n";
					$ncharTNT++;
				}
			}
			if($#indelNames >= 0){ ### PRINT INDEL PARTITION DATA FOR PAR
				$ncharPHY++;
				my $start = $ncharPHY;
				$ncharPHY += $#indelNames;
				if($r == 1){
					$buffers{'par'} .= 'BIN, ';
				} else {
					$buffers{'par'} .= 'BIN+G, ';
				}
				$buffers{'par'} .= $file . '_indel = ' . $start . '-' . $ncharPHY . "\n";
			}
		}
		if($g == 1){ ### PRINT GENE ABSENCE/PRESENCE (GAP) CHARACTER
			if($informativeGAP){
				$buffers{'tnt'} .= "&[num]\n";
			}
			$buffers{'phy'} .= "\n";
			for(my $k = $#species; $k >= 0; $k--){ ### PRINT FOR TNT AND PHY
				if(exists($matrix->{$species[$k]})){
					if($informativeGAP){
						$buffers{'tnt'} .= $matrixSpeciesTNT[$k] . "1\n";
					}
					$buffers{'phy'} .= "1\n";
				} else {
					if($informativeGAP){
						$buffers{'tnt'} .= $matrixSpeciesTNT[$k] . "0\n";
					}
					$buffers{'phy'} .= "0\n";
				}
			}
			if($informativeGAP){ ### PRINT GAP PARTITION DATA FOR CNA
				$buffers{'cna'} .= '{' . $ncharTNT . ' ' . $file . " absent present /;\n";
				$ncharTNT++;
			}
			$ncharPHY++;
			if($r == 1){  ### PRINT GAP PARTITION DATA FOR PAR
				$buffers{'par'} .= 'BIN, ' . $file . ' = ' . $ncharPHY . '-' . $ncharPHY . "\n";
			} else {
				$buffers{'par'} .= 'BIN+G, ' . $file . ' = ' . $ncharPHY . '-' . $ncharPHY . "\n";
			}
		}
		for(my $k = $#files; $k >= 0; $k--){
			if(length($buffers{$files[$k]}) > 10000){
				print({$files{$files[$k]}} $buffers{$files[$k]});
				$buffers{$files[$k]} = '';
			}
		}
		if($aa == 1){
			$AA = 1;
		}
	}
	for(my $k = $#files; $k >= 0; $k--){
		print({$files{$files[$k]}} $buffers{$files[$k]});
		$buffers{$files[$k]} = '';
		close($files{$files[$k]});
	}
	open(TNTFILE, '>', $n . '.tnt'); ### ADD TNT HEADERS
	my $ram = 0;
	if($AA == 1){
		if($ax =~ m/^(2|3|4|5|6|6dso|6kgb|6sr|8)$/){
			$ram = $ncharTNT;
		} elsif($ax =~ m/^(10|11|12|15)$/){
			$ram = $ncharTNT*2;
		} else {
			$ram = $ncharTNT*4;
		}
	} else {
		$ram = $ncharTNT;
	}
	$ram = int($ram/1024)+2048; ### B -> mB
	$buffers{'tnt'} .= "mxram $ram;\n";
	$buffers{'tnt'} .= "taxname+64;\n";
	$buffers{'tnt'} .= "taxonomy=;\n";
	if($AA == 1){
		if($ax =~ m/^(2|3|4|5|6|6dso|6kgb|6sr|8)$/){
			$buffers{'tnt'} .= "nstates 8;\n";
		} elsif($ax =~ m/^(10|11|12|15)$/){
			$buffers{'tnt'} .= "nstates 16;\n";
		} else {
			$buffers{'tnt'} .= "nstates 32;\n";
		}
	} else {
		$buffers{'tnt'} .= "nstates 8;\n";
	}
	$buffers{'tnt'} .= "xread\n";
	my ($day, $month, $year) = (localtime)[3,4,5];
	my @months = ('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December');
	$buffers{'tnt'} .= "'Data processed on " . ($year+1900) . ' ' . $months[$month] . ' ' . $day . ' using BAD2matrix.pl ' . join(' ', @ARGV);
	if($i == 1){
		$buffers{'tnt'} .= "\nindel characters coded using the \"simple gap coding\" method of SIMMONS AND OCHOTERENA. 2000. Gaps as characters in sequence-based phylogenetic analysis. Systematic Biology 49: 369-381. DOI 10.1080/10635159950173889";
	}
	if($ax =~ m/^(2|3|4|5|6|8|10|12|15|18)$/){
		$buffers{'tnt'} .= "\namino acid states reduced to $ax following Murphy, Wallqvist, and Levy. 2000. Simplified amino acid alphabets for protein fold recognition and implications for folding. Protein Engineering 13: 149–52. DOI 10.1093/protein/13.3.149";
	} elsif($ax eq '6dso'){
		$buffers{'tnt'} .= "\namino acid states reduced to 6 following Dayhoff, Schwartz, and Orcutt. 1978. A model of evolutionary change in proteins. In Atlas of Protein Sequence and Structure, Dayhoff ed. pp. 345–352";
	} elsif($ax eq '6kgb'){
		$buffers{'tnt'} .= "\namino acid states reduced to 6 following Kosiol, Goldman, and Buttimore. 2004. A new criterion and method for amino acid classification. Journal of Theoretical Biology 228:97–106. DOI 10.1016/j.jtbi.2003.12.010";
	} elsif($ax eq '6sr'){
		$buffers{'tnt'} .= "\namino acid states reduced to 6 following Susko and Roger. 2007. On reduced amino acid alphabets for phylogenetic inference. Molecular Biology and Evolution 24: 2139–2150. DOI 10.1093/molbev/msm144";
	} elsif($ax eq '11'){
		$buffers{'tnt'} .= "\namino acid states reduced to 11 following Buchfink, Xie, and Huson. 2015. Fast and sensitive protein alignment using DIAMOND. Nature Methods 12: 59–60. DOI 10.1038/nmeth.3176";
	}
	$buffers{'tnt'} .= "'\n";
	$buffers{'tnt'} .= $ncharTNT . ' ' . ($#species+1) . "\n";
	print(TNTFILE $buffers{'tnt'});
	close(TNTFILE);
	qx`cat tnt-temporary-do-not-touch cna-temporary-do-not-touch >> $n.tnt`;
	qx`find . -type f -name '0x*-temporary-do-not-touch' | xargs cat > $n.fas`;
	open(PHYFILE, '>', $n . '.phy'); ### ADD PHY HEADERS
	print(PHYFILE ($#species+1) . ' ' . $ncharPHY . "\n");
	close(PHYFILE);
	qx`cat phy-temporary-do-not-touch >> $n.phy`;
	rename('par-temporary-do-not-touch', $n . '.part');
	for(my $k = $#files; $k >= 0; $k--){
		unlink($files[$k] . '-temporary-do-not-touch');
	}
	############################### OPTIONS
} else {
	print(STDERR "\nA Perl script for merging and translating FASTA alignments into TNT extended\n");
	print(STDERR "XREAD, FastTree FASTA, and RAxML/IQ-Tree extended PHYLIP with indel characters\n");
	print(STDERR "(optionally, .tnt & .phy) coded using the 'simple' gap coding method of Simmons\n");
	print(STDERR "and Ochoterena (2000; Gaps as characters in sequence–based phylogenetic\n");
	print(STDERR "analysis. Systematic Biology 49: 369–381. DOI 10.1080/10635159950173889), and\n");
	print(STDERR "with coded gene content (absence/presence) characters (optionally, .tnt & .phy).\n");
	print(STDERR "This script is slower than 2matrix.pl due to more disk use, but will not run out\n");
	print(STDERR "of RAM (hopefully). The RAxML .part file is either in original RAxML (-r) or\n");
	print(STDERR "RAxML-NG format.\n\n");
	print(STDERR "USAGE:\tBAD2matrix.pl [ -a 2|3|4|5|6|6dso|6kgb|6sr|8|10|11|12|15|18|20 ]\n");
	print(STDERR "\t-d <directory> [ -f ] [ -g ] [ -i ] [ -m int ] -n <root-name>\n");
	print(STDERR "\t[ -o speciesA,speciesB... ] [ -r ]\n\n");
	print(STDERR "OPTIONS:\n-a\tNumber of amino acid states (default = $ax). Reduction with option '6dso'\n");
	print(STDERR "\tfollows Dayhoff, Schwartz, and Orcutt (1978. A model of evolutionary\n");
	print(STDERR "\tchange in proteins. In Atlas of Protein Sequence and Structure, Dayhoff\n");
	print(STDERR "\ted. pp. 345–352); option '6kgb' follows Kosiol, Goldman, and Buttimore\n");
	print(STDERR "\t(2004. A new criterion and method for amino acid classification. Journal\n");
	print(STDERR "\tof Theoretical Biology 228:97–106. DOI 10.1016/j.jtbi.2003.12.010);\n");
	print(STDERR "\toption '6sr' follows Susko and Roger (2007. On reduced amino acid\n");
	print(STDERR "\talphabets for phylogenetic inference. Molecular Biology and Evolution 24:\n");
	print(STDERR "\t2139–2150. DOI 10.1093/molbev/msm144); option '11' follows Buchfink, Xie,\n");
	print(STDERR "\tand Huson (2015. Fast and sensitive protein alignment using DIAMOND.\n");
	print(STDERR "\tNature Methods 12: 59–60. DOI 10.1038/nmeth.3176); and all other options\n");
	print(STDERR "\tfollow Murphy, Wallqvist, and Levy (2000. Simplified amino acid alphabets\n");
	print(STDERR "\tfor protein fold recognition and implications for folding. Protein\n");
	print(STDERR "\tEngineering 13: 149–52. DOI 10.1093/protein/13.3.149).\n");
	print(STDERR "-d\tSpecifies a directory of aligned FASTA files. By default, names should\n");
	print(STDERR "\tuse the OrthologID naming convention ('>species#sequenceID').\n");
	print(STDERR "\tCharacters other than letters, numbers, periods, and underscores will be\n");
	print(STDERR "\tdeleted. Use -f for an alternate naming convention.\n");
	print(STDERR "-f\tUse full FASTA names. Characters other than letters, numbers, periods,\n");
	print(STDERR "\tand underscores will be deleted.\n");
	print(STDERR "-g\tDo NOT code gene content.\n");
	print(STDERR "-i\tDo NOT code indels.\n");
	print(STDERR "-m\tRetain the upper x percentile of genes in the distribution of missing\n");
	print(STDERR "\tsequences (default = $m; i.e. include all genes).\n");
	print(STDERR "-n\t<root-name> for output files.\n");
	print(STDERR "-o\tOutgroup(s) for rooting trees.\n");
	print(STDERR "-r\tOutput original RAxML formatted .part file.\n\n");
}
sub char2num {
	my $sequence = $_[0];
	my $ax = $_[1];
	my $phy = $_[2];
	if($ax eq '20'){
		return($sequence);
	} elsif($ax eq '2'){
		$sequence =~ tr/FEX\-/01?\-/;
	} elsif($ax eq '3'){
		$sequence =~ tr/IEFX\-/012?\-/;
	} elsif($ax eq '4'){
		$sequence =~ tr/IPFEX\-/0123?\-/;
	} elsif($ax eq '5'){
		$sequence =~ tr/IPFEKX\-/01234?\-/;
	} elsif($ax eq '6'){
		$sequence =~ tr/IAPFEKX\-/012345?\-/;
	} elsif($ax eq '6dso'){
		$sequence =~ tr/PEHIFCX\-/012345?\-/;
	} elsif($ax eq '6kgb'){
		$sequence =~ tr/PEIWFCX\-/012345?\-/;
	} elsif($ax eq '6sr'){
		$sequence =~ tr/PEQICFX\-/012345?\-/;
	} elsif($ax eq '8'){
		$sequence =~ tr/IASPFEKHX\-/01234567?\-/;
	} elsif($ax eq '10'){
		$sequence =~ tr/ICAGSPFEKHX\-/0123456789?\-/;
	} elsif($ax eq '11'){
		$sequence =~ tr/ECGHIMFYWPAX\-/0123456789A?\-/;
	} elsif($ax eq '12'){
		$sequence =~ tr/ICAGSPFWEDKHX\-/0123456789AB?\-/;
	} elsif($ax eq '15'){
		$sequence =~ tr/ICAGSTPFWEDNQKHX\-/0123456789ABCDE?\-/;
	} elsif(($ax eq '18') && ($phy == 1)){
		$sequence =~ tr/LICAGSTPFYWEDNQKRHX\-/0123456789ABCDEFGH?\-/;
	}
	return($sequence);
}
sub cleanName {
	my $x = $_[0];
	$x =~ tr/ \/\\\-/\_\_\_\_/;
	$x =~ tr/[A-z][0-9]\_\.//cd;
	return($x);
}
sub insert {
	my $matrix = $_[0];
	my $sequence = uc($_[1]);
	my $species = cleanName($_[2]);
	my $ax = $_[3];
	$sequence =~ tr/\?/\-/;
	$sequence =~ tr/ABCDEFGHIKLMNOPQRSTUVWXYZ\-//cd;
	if($ax ne '20'){
		if($ax eq '2'){
			$sequence =~ tr/LVIMCAGSTPFYW EDNQKRH                 \-/FFFFFFFFFFFFF EEEEEEE                 \-/;
		} elsif($ax eq '3'){
			$sequence =~ tr/LASGVTIPMC EKRDNQH FYW                \-/IIIIIIIIII EEEEEEE FFF                \-/;
		} elsif($ax eq '4'){
			$sequence =~ tr/LVIMC AGSTP FYW EDNQKRH               \-/IIIII PPPPP FFF EEEEEEE               \-/;
		} elsif($ax eq '5'){
			$sequence =~ tr/LVIMC ASGTP FYW EDNQ KRH              \-/IIIII PPPPP FFF EEEE KKK              \-/;
		} elsif($ax eq '6'){
			$sequence =~ tr/LVIM ASGT PHC FYW EDNQ KR             \-/IIII AAAA PPP FFF EEEE KK             \-/;
		} elsif($ax eq '6dso'){
			$sequence =~ tr/AGPST DENQ HKR MIVL WFY C             \-/PPPPP EEEE HHH IIII FFF C             \-/;
		} elsif($ax eq '6kgb'){
			$sequence =~ tr/AGPS DENQHKRT MIL W FY CV             \-/PPPP EEEEEEEE III W FF CC             \-/;
		} elsif($ax eq '6sr'){
			$sequence =~ tr/APST DENG QKR MIVL WC FYH             \-/PPPP EEEE QQQ IIII CC FFF             \-/;
		} elsif($ax eq '8'){
			$sequence =~ tr/LVIMC AG ST P FYW EDNQ KR H           \-/IIIII AA SS P FFF EEEE KK H           \-/;
		} elsif($ax eq '10'){
			$sequence =~ tr/LVIM C A G ST P FYW EDNQ KR H         \-/IIII C A G SS P FFF EEEE KK H         \-/;
		} elsif($ax eq '11'){
			$sequence =~ tr/KREDQN C G H ILV M F Y W P STA        \-/EEEEEE C G H III M F Y W P AAA        \-/;
		} elsif($ax eq '12'){
			$sequence =~ tr/LVIM C A G ST P FY W EQ DN KR H       \-/IIII C A G SS P FF W EE DD KK H       \-/;
		} elsif($ax eq '15'){
			$sequence =~ tr/LVIM C A G S T P FY W E D N Q KR H    \-/IIII C A G S T P FF W E D N Q KK H    \-/;
		} elsif($ax eq '18'){
			$sequence =~ tr/LM VI C A G S T P F Y W E D N Q K R H \-/LL II C A G S T P F Y W E D N Q K R H \-/;
		}
	}
	if(length($sequence) && length($species)){
		my @residues = split(//, $sequence);
		for(my $k = $#residues; $k >= 0; $k--){
			$matrix->{$species}->[$k] = $residues[$k];
		}
	}
	return($sequence =~ m/E|F|I|L|O|P|Q|X|Z/);
}
############################### FOWLER–NOLL–VO 32-BIT HASH
sub fnv32a {
	my @data = split(//, $_[0]);
	my $hash = 0x811c9dc5;
	for(my $k = 0; $k <= $#data; $k++){
		$hash = $hash^ord($data[$k]);
		$hash = ($hash*0x01000193)%4294967296; ### 2**32
	}
	return(sprintf("0x%X", $hash));
}
exit(0);

# raxmlHPC -f c -q A_min.part -m MULTIGAMMA -p 1 -n x -s A_min.phy -K GTR
# raxml-ng --check --model A_min.part --msa A_min.phy
