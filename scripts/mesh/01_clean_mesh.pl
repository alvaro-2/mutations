open(IN_FILE, "../../raw_data/d2021.bin");
open(OUT_FILE, ">../../raw_data/d2021_processed.tsv");

#perl 01_clean_mesh.pl
#https://nlmpubs.nlm.nih.gov/projects/mesh/2021/
#head target screen
#*NEWRECORD
#RECTYPE = C
#NM = N-acetylglucosaminylasparagine
#UI = C000591739
#NM_TH = OMIM (2013)
#MN = L01.559.598.400.556.131
#variables texto/numero
$ui = "";
$mh = "";
#variables arreglos/lista
@t = ();
@cr = ();

print OUT_FILE "ui\tmesh_name\ttree\tcross_reference\n";
while( $linea = <IN_FILE> ) {
	chomp $linea;
	if ($linea =~ /^\*NEWRECORD$/){	
		if ($ui ne ""){
			$p1 = join(",", @t);
			$p2 = join(",", @cr);
			print OUT_FILE "$ui\t$mh\t$p1\t$p2\n";
		}
		$ui = "";
		$mh = "";
		@t = ();
		@cr = ();
	}else{
		if ($linea =~ /^MH\s+=\s+(.*)$/){	
			$mh = $1;
		}elsif ($linea =~ /^UI\s+=\s+(.*)$/){	
			$ui = $1;
		}elsif ($linea =~ /^MH_TH\s+=\s+(.*)$/){	
			push @cr, $1;
		}elsif ($linea =~ /^MN\s+=\s+(.*)$/){	
			push @t, $1;
		}
	}		
}

close IN_FILE;
close OUT_FILE;
exit;
