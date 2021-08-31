$file_n = $ARGV[0]; #../raw_data/uniprot_all_data_associated_to_mlo.txt
$file_o = $ARGV[1]; #../raw_data/uniprot_all_data_associated_to_mlo_preproc.tsv
open(ALL_MUT, $file_n);
open(SEL_MUT, ">". $file_o);
print SEL_MUT "id_protein\tuniprot_name\tuniprot_status\tlength\tuniprot_acc\tuniprot_other_accesions\tprotein_names\tflags\tgene_name\tgene_name_synonyms\tisoform_canonic\tisoforms\tccds_canonic\thgnc_id\thgnc_name\tgene_id\trefseq_canonic\tensembl_canonic\tchromosome\tsequence\n";

$isof_name = "";
@d = ("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "");
#$d[0] = uniprot_name
#$d[1] = uniprot_status
#$d[2] = length
#$d[3] = uniprot_primary
#$d[4] = uniprot_other_acccessions, can be more than one...
#$d[5] = protein_names
#$d[6] = flags
#$d[7] = gene_name
#$d[8] = gene_name_synomyms
#$d[9] = isoform_canonic
#$d[10] = isoforms
#$d[11] = ccds_canonic
#$d[12] = hgnc_id
#$d[13] = hgnc_name
#$d[14] = gene_id
#$d[15] = refseq_canonic
#$d[16] = ensembl_canonic
#$d[17] = chromosome
#$d[18] = sequence

$protein_id = 1;
$find_iso = "";
while($line = <ALL_MUT>) {
	chomp $line;
	if ($d[9] ne ""){ #buscar por isoforma canonica
        	$find_iso = "\\[" . $d[9] . "\\]";
	}
	if ($line =~ /^ID\s{3}([^\s]+)\s+([^\s]+);\s+(\d+)\s+AA\.$/){
        	$d[0] = $1;
        	$d[1] = $2;
        	$d[2] = $3;  	
    }elsif ($line =~ /^AC\s{3}(.*)$/){    
        #puede haber varias lineas
        $d[3] = $d[3] . $1; 
    }elsif ($line =~ /^DE\s{3}RecName:\s*Full=(.+);$/){    
        #si es un reviewed
        $d[5] = $1;  
    }elsif ($line =~ /^DE\s{3}SubName:\s*Full=(.+);$/){ 
        #si es un unreviewed
        $d[5] = $1;                    
    }elsif ($line =~ /^DE\s{12}Short=(.+);$/){    
        $d[5] = $d[5] . " (" . $1 . ")";  
    }elsif ($line =~ /^DE\s{3}AltName:\s*Full=(.+);$/){    
        $d[5] = $d[5] . " (" . $1 . ")";      
    }elsif ($line =~ /^DE\s{3}Flags:\s*(.+)$/){
        $d[6] = $1; 
    }elsif ($line =~ /^GN\s{3}(.+)$/){
        #hay lineas con "and" cuando tienen mas de un gene name
        $d[7] = $d[7] . $1;  
    }elsif ($line =~ /^CC\s{9}IsoId=(.+);\s*Sequence=(.+);$/){
        $isof_name = $1;
        @y = split(",", $isof_name, 2);
        $isof_name = $y[0];
        	$isof_name =~ s/ //g;
        $i_s = $2;
        @y = split(",", $i_s, 2);
        $i_s = $y[0];        
        	$i_s =~ s/ //g;
        	if ($i_s eq "Displayed"){ 
            	$d[9] = $isof_name;
        	}
        	$d[10] = ($d[10] eq "" ? "" : $d[10]. ";"). $isof_name;
    }elsif ($line =~ /^CC\s{9}IsoId=(.+);$/){
        $isof_name = $1;
        @y = split(",", $isof_name, 2);
        $isof_name = $y[0];
        	$isof_name =~ s/ //g;
        	$d[10] = ($d[10] eq "" ? "" : $d[10]. ";"). $isof_name;    	
    }elsif ($line =~ /^CC\s{9}Sequence=(.+);$/){
        $i_s = $1;
        @y = split(",", $i_s, 2);
        $i_s = $y[0];
        	$i_s =~ s/ //g;
        if ($i_s eq "Displayed"){ 
            	$d[9] = $isof_name;
        	}
    }elsif ($line =~ /^DR\s{3}CCDS;\s*(.+)$find_iso$/){  
        $d[11] =  ($d[11] eq "" ? "" : $d[11]. ";"). $1;  
    }elsif ($line =~ /^DR\s{3}HGNC;\s*(.+)$/){  
        $d[12] =  ($d[12] eq "" ? "" : $d[12]. ";"). $1;  
    }elsif ($line =~ /^DR\s{3}GeneID;\s*(.+)$/){  
        $d[14] =  ($d[14] eq "" ? "" : $d[14]. ";"). $1;      
    }elsif ($line =~ /^DR\s{3}RefSeq;\s*(.+)$find_iso$/){  
        $d[15] =  ($d[15] eq "" ? "" : $d[15]. ";"). $1;  
    }elsif ($line =~ /^DR\s{3}Ensembl;\s*(.+)$find_iso$/){  
        $d[16] =  ($d[16] eq "" ? "" : $d[16]. ";"). $1;     
    }elsif ($line =~ /^DR\s{3}Proteomes;.+Chromosome\s*([A-Z0-9]{1,2}).*$/){
        #DR   Proteomes; UP000005640; Chromosome 19.  
        $d[17] =  ($d[17] eq "" ? "" : $d[17]. ";"). $1;  
    }elsif ($line =~ /^DR\s{3}Proteomes;.+Mitochondrion.*$/){
        $d[17] =  ($d[17] eq "" ? "" : $d[17]. ";"). "MT";          
    }elsif ($line =~ /^\s+(.+)$/){
    #sequence 
        $d[18] =  $d[18] . $1;    	
	}elsif ($line =~ /^\/\//){
        	$d[3] =~ s/ //g; #uniprot_primary
        	$d[3] =~ s/;$//g; 
        	@y = split(";", $d[3], 2);
        	$d[3] = $y[0]; #uniprot_others     	
        	$d[4] = (scalar(@y)== 2 ? $y[1] : "");
        	
        	$d[5] =~ s/\s*{[^{]+}//g; #protein_names
        	
        	$d[6] =~ s/;$//g; #flags
        	
        	$d[7] =~ s/ //g;#gene_name
        	$d[7] =~ s/and//g;
        	$d[7] =~ s/\s*{[^{]+}//g;
        	$d[7] =~ s/ORFNames=([^;]+);//g;
        	$d[8] = $d[7];
        	$d[7] =~ s/Synonyms=([^;]+);//g;
        	$d[7] =~ s/Name=//g;
        	$d[7] =~ s/;$//g; 
        	$d[8] =~ s/Name=([^;]+);//g; #gene_name_synomyms
        	$d[8] =~ s/Synonyms=//g;
        	$d[8] =~ s/;$//g;  
        	
        	$d[11] =~ s/ //g; #ccds_canonic
        	$d[11] =~ s/-\.;?//g;
        	$d[11] =~ s/;$//g; 
        	
        $d[12] =~ s/ //g; #hgnc_id
        $d[13] = $d[12]; #hgnc_name
        $d[12] =~ s/;([-_A-Za-z0-9]+)\.//g; #con esto se quedan solo los HGNC:\d+
        $d[12] =~ s/HGNC://g; #con esto se quedan solo los numeros sin HGNC
        $d[13] =~ s/HGNC:(\d+);//g; #con esto se quedan los hgnc_symbol        
        $d[13] =~ s/\.//g;     
           
        $d[14] =~ s/ //g; #gene_id
        	$d[14] =~ s/-\.;?//g;
        	$d[14] =~ s/;$//g; 
        	
        $d[15] =~ s/ //g; #refseq_canonic
        	$d[15] =~ s/\.$//g; 
        	
        $d[16] =~ s/ //g; #ensembl_canonic
        	$d[16] =~ s/\.$//g; 
        	
        	$d[18] =~ s/ //g; #sequence
        	
        	my $new_line = join("\t", @d);
		$new_line = $protein_id . "\t" . $new_line . "\n";
		print SEL_MUT $new_line;
		@d = ("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "");
		$find_iso = "";
		$isof_name = "";
		$protein_id = $protein_id + 1;
		
    	}
}

close ALL_MUT;
close SEL_MUT;
exit;
