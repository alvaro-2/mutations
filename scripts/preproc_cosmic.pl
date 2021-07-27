$file_n = $ARGV[0]; #CosmicMutantExport.tsv
open(ALL_MUT, "" . $file_n . ".tsv");
open(SEL_MUT, ">" . $file_n . "_sel.tsv");
print SEL_MUT "id_sample\tGENOMIC_MUTATION_ID\tLEGACY_MUTATION_ID\tMUTATION_ID\tCDS\tAA\tconsequence\tgenomic\tsomatic_status\tpubmed\tensembl_ensp\tensembl_enst\n";
while($line = <ALL_MUT>) {
	chomp $line;
	@cols = split(/\t/, $line);
	#remove rows without muttaions in coding part ;
	#remove the non canonic transcript
	@ttt = split(/(_ENST)|\./, $cols[0]);	
    	if ((scalar(@ttt) == 1) and ($cols[24] eq '38') and ($cols[19] ne "c.?") and ($cols[20] ne "p.?")) {
    		#remove the columns 0-4, 6-15, 22-24, 26-30, 33-36, 39-39;
    		splice @cols, 39, 1;
		splice @cols, 33, 4;
		splice @cols, 26, 5;
		splice @cols, 22, 3;
		splice @cols, 6, 10;
		splice @cols, 0, 5;
		#0: id_sample, 1-3: notation of mutation, 4: cds, 5: aa, 6:consequence
		#7:genomic, 8:somatic status, 9: PUBmed, 10: protein, 11:transcript
		#remove version transcript;
		@ttt = split(/\./, $cols[10]); #ensembl_ensp.version:p.XXX
		$cols[10] = $ttt[0];	
		#remove notation in protein;	
		@ttt = split(/\./, $cols[11]);
		$cols[11] = $ttt[0];			
    		my $new_line = join("\t", @cols);
		$new_line = $new_line . "\n";
		print SEL_MUT $new_line;
    	}
}
close ALL_MUT;
close SEL_MUT;
exit;
