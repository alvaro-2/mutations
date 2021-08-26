$file_n = $ARGV[0]; #CosmicMutantExport.tsv
$file_o = $ARGV[1]; #CosmicMutantExport_sel.tsv
open(ALL_MUT, "" . $file_n);
open(SEL_MUT, ">" . $file_o);
print SEL_MUT "accesion\tgene_name\tenst\thgnc_id\tid_sample\tPrimary_site\tSite_subtype_1\tSite_subtype_2\tSite_subtype_3\tPrimary_histology\tHistology_subtype_1\tHistology_subtype_2\tHistology_subtype_3\tLEGACY_MUTATION_ID\tMUTATION_ID\tCDS\tAA\tchromosome\tstart\tend\tsomatic_status\tpubmed\n";
while($line = <ALL_MUT>) {
	chomp $line;
	@cols = split(/\t/, $line);
	#remove rows without mutations in coding part "p.?";
    	if (($cols[24] eq '38') and ($cols[19] ne "c.?") and ($cols[20] ne "p.?")) {
    		#kept 0, 1, 3, 5, 7-14, 17-18, 19, 20, 25, 31, 32
    		#0-1: accesion, gene_name, 2: enst, 3: hgnc, 4: id_sample, 
    		#5-12: primary site, histology and subhistologys,
    		#13: legacy mutation, 14: mutation id, 
		#15: cds, 16: aa, 17-19:chromosome, start, end, 20:somatic status, 21: pubmed
    		#gene name notation    		
        	@ttt = split(/(_ENST)|\./, $cols[0]);
    		$cols[0] = $cols[0] . "\t". $ttt[0];	
    		#remove version transcript;
		@ttt = split(/\./, $cols[1]); #ensembl_enst.version
		$cols[1] = $ttt[0];	
		#replace : and - by \t to get cromosome, start and end position
    		$cols[25] =~ s/:|-/\t/g;
    		#remove the columns 2, 4, 6, 15-16, 21-24, 26-30, 33-39;
		splice @cols, 33, 7;
		splice @cols, 26, 5;
		splice @cols, 21, 4;
		splice @cols, 15, 2;
		splice @cols, 6, 1;		
		splice @cols, 4, 1;
		splice @cols, 2, 1;			
    		my $new_line = join("\t", @cols);
		$new_line = $new_line . "\n";
		print SEL_MUT $new_line;
    	}
}
close ALL_MUT;
close SEL_MUT;
exit;
