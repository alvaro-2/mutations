$file_n = $ARGV[0]; #CosmicMutantExport.tsv
$file_o = $ARGV[1]; #CosmicMutantExport_sel.tsv
open(ALL_MUT, "" . $file_n);
open(SEL_MUT, ">" . $file_o);
print SEL_MUT "enst\tPrimary_site\tSite_subtype_1\tSite_subtype_2\tSite_subtype_3\tPrimary_histology\tHistology_subtype_1\tHistology_subtype_2\tHistology_subtype_3\tGENOMIC_MUTATION_ID\tCDS\tAA\tchromosome\tstart\tend\tpubmed\n";
while($line = <ALL_MUT>) {
	chomp $line;
	@cols = split(/\t/, $line);
	#remove rows without mutations in coding part "p.?";
    	if (($cols[24] eq '38') and ($cols[19] ne "c.?") and ($cols[20] ne "p.?")) {
    		#kept 1, 7-14, 16, 19, 20, 25, 32
    		#0: enst, 1-8: primary site, histologys and subhistologys,
    		#9: genomic mutation id, 10: cds, 11: aa, 
		#12-14:chromosome, start, end, 15: pubmed
    		
    		#remove version transcript;
		@ttt = split(/\./, $cols[1]); #ensembl_enst.version
		$cols[1] = $ttt[0];	
		#replace : and - by \t to get cromosome, start and end position
    		$cols[25] =~ s/:|-/\t/g;
    		#remove the columns 0, 2-6, 15, 17-18, 21-24, 26-31, 33-39;
		splice @cols, 33, 7;
		splice @cols, 26, 6;
		splice @cols, 21, 4;
		splice @cols, 17, 2;
		splice @cols, 15, 1;	
		splice @cols, 2, 5;	
		splice @cols, 0, 1;		
    		my $new_line = join("\t", @cols);
		$new_line = $new_line . "\n";
		print SEL_MUT $new_line;
    	}
}
close ALL_MUT;
close SEL_MUT;
exit;
