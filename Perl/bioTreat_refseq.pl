use strict;
use warnings;
use Bio::Perl;
use Bio::SeqIO;
use Bio::DB::GenBank;
use Bio::DB::Query::GenBank;
use Bio::DB::SwissProt;
use Bio::DB::EntrezGene;
use Bio::DB::Taxonomy;
#use IO::String;
use Bio::SeqFeatureI;
use threads;
use List::MoreUtils;
use List::Util;


############
#import personal modules
use lib "/home/yuan/eRNA/bin";
use func_basic;
use func_common;
use func_data;
use func_bioperl;
use func_bioseq;
use func_phip;

#################
#main
#################
sub main{
	my($lib, $task)=@_;
	$lib='' unless $lib;
	#
	my %attr;
	#print "###read peptides\n";
	$attr{'lib'}=$lib;
	$attr{'task'}=$task;
	$attr{'out_dir'}='/home/yuan/phip/ref_seq/';
	$attr{'annot_dir'}='/home/yuan/data_preparation/phip_annot/';
	$attr{'task_dir'}=$attr{'annot_dir'}.$lib.'_'.$task.'/';
	mkdir($attr{'task_dir'}, 0755) unless (-d $attr{'task_dir'} or $task eq 'export');
	
	#peptides file
	$attr{'phip_fa'}=$attr{'annot_dir'}.'T7Pep2_VirScan.fa';
	$attr{'phip_human1'}=$attr{'annot_dir'}.'Misc_human_annot.txt';
	$attr{'phip_virus1'}=$attr{'annot_dir'}.'VirScanSupplement_annot.txt';
	$attr{'phip_toxome1'}=$attr{'annot_dir'}.'input_Toxome1.txt';
	$attr{'phip_toxome2'}=$attr{'annot_dir'}.'input_Toxome2.txt';
	$attr{'phip_provirome'}=$attr{'annot_dir'}.'input_ProkaryoticVirome.csv';
	$attr{'phip_allergome1'}=$attr{'annot_dir'}.'input_Allergome1.txt';
	$attr{'phip_allergome2'}=$attr{'annot_dir'}.'input_Allergome2.txt';
	$attr{'phip_mouse'}=$attr{'annot_dir'}.'mouse-annot.tsv';
	#
	$attr{'dependent_file'}=$attr{'out_dir'}.'virus_dependent_peptides.csv';
	
	#extract control peptide
	my $control_pep_pointer=sub_phip::extract_control_pep($attr{phip_fa});
	$attr{'control_pep_pointer'}=$control_pep_pointer;

	#exported columns in annot.txt
	my @pro_cols=('pep_rank', 'UniProt_acc', 'pep_dna', 'pep_aa', 'pro_len', 'pro_motifs', 'GO',);
	$attr{pro_cols_pointer}=\@pro_cols;

	######################################
	if ($attr{'lib'} eq 'human'){
		human(\%attr);
	}elsif ($attr{'lib'} eq 'virus'){
		virus(\%attr);
	}elsif ($attr{'lib'} eq 'toxome'){
		toxome(\%attr);
	}elsif ($attr{'lib'} eq 'provirome'){
		provirome(\%attr);
	}elsif ($attr{'lib'} eq 'allergome'){
		allergome(\%attr);
	}elsif ($attr{'lib'} eq 'mouse'){
        mouse(\%attr);
    }
	#
}


######################################
###human
sub human{
	my($basic_attr_pointer)=@_;
	my %basic_attr=%$basic_attr_pointer;
	
	#read peptides
	my $attr_pointer =sub_phip::A_extract_human(\%basic_attr);
	my %attr=%$attr_pointer;
	#query accession
	my $query_accs_pointer=sub_phip::A_judge_query_acc($attr_pointer);
	my @query_accs=@$query_accs_pointer;

	#
	if($basic_attr{'task'} eq 'GenBank'){
		print "###query1: GenBank\n";
		sub_common::multi_threading(\&sub_bioperl::query_GenBank, 8, \@query_accs, $basic_attr{task_dir});
	}elsif($basic_attr{'task'} eq 'UniProt'){
		print "###query 2: SwissProt query\n";
		sub_common::multi_threading(\&sub_bioperl::query_SwissProt, 8, \@query_accs, $basic_attr{task_dir});
	}elsif($basic_attr{'task'} eq 'ProMotifs'){
		print "###query 5: motifs searching\n";
		sub_phip::C_protein_motifs($attr_pointer);
	}elsif($basic_attr{'task'} eq 'EntrezGene'){
		print "###query 3: EntrezGene query\n";
		sub_common::multi_threading(\&sub_bioperl::query_EntrezGene, 8, \@query_accs, $basic_attr{task_dir});
	}elsif($basic_attr{'task'} eq 'RefGenome'){
		print "###query 4: newest release of NCBI genome\n";
		my $gff_file=$basic_attr{'annot_dir'}.'ref_GRCh38.p7_top_level.gff3';
		sub_phip::human_RefGenome_query(\@query_accs, $basic_attr{task_dir}, $gff_file);
	}elsif($basic_attr{'task'} eq 'PPI'){
		print "###query 5: human PPI searching\n";
		sub_phip::human_PPI_query($attr_pointer);
	}elsif($basic_attr{'task'} eq 'autoantigen'){
		print "###query 6: human: integrate autoantigen\n";
		sub_phip::human_autoantigen($attr_pointer);
	}elsif($basic_attr{'task'} eq 'export'){
		#combination
		my @annot_dirs=('UniProt', 'EntrezGene', 'PPI', 'RefGenome', 'GenBank', 'ProMotifs'); 
		@annot_dirs=map($basic_attr{annot_dir}.'human_'.$_.'/', @annot_dirs);
		my $human_pointer=sub_phip::human_combine_annot($attr{annot_pointer}, \@annot_dirs);
		
		print "export human annotations\n";
		my @human_cols=@{$basic_attr{pro_cols_pointer}};
		my @human_add=('gene_symbol', 'product', 'transcript_id', 'chromosome', 'pos_start', 'pos_end', 'seq_len', 'map', 
				'mol_type', 'organism', 'Ensembl', 'gene_synonym', 'GeneID', 'GI', 'protein_id', 'note', 'Exon count', 
				'HGNC', 'HPRD', 'MIM', 'GO_Entrez','GO_BP','GO_CC','GO_MF', 'KEGG', 'InterPro', 'PPI', 'autoantigen',);
		push(@human_cols, @human_add);
		sub_phip::D_export_fa($human_pointer, \@human_cols, $attr{out_dir}.'T7Pep2_human', $attr{pros_pointer});
		sub_phip::D_update_fa($attr{control_pep_pointer}, \@human_cols, $attr{out_dir}.'T7Pep2_human');
	}
	
	#match the order peptides between fa, annot.txt, and old matrix.
	#my $RC_file='/home/yuan/data_1/results_phip/phipseq5_human/statistics/lowRC.txt';
	#my $human_peps_pointer=sub_phip::read_pep_order($RC_file);
	#my $fa_file='/home/yuan/phip/ref_seq/T7Pep2_human.fa';
	#sub_phip::order_pep_fa($human_peps_pointer, $fa_file);
	#my $annot_file='/home/yuan/phip/ref_seq/T7Pep2_human_annot.txt';
	#sub_phip::order_pep_annot($human_peps_pointer, $annot_file);
}

######################################
#provirome annotations
sub provirome{
	my($basic_attr_pointer)=@_;
	my %basic_attr=%$basic_attr_pointer;
	my @provirome_cols=@{$basic_attr{pro_cols_pointer}};
	
	#1: extract peptides
	my $attr_pointer=sub_phip::A_extract_provirome($basic_attr_pointer);
	#query accession
	my $query_accs_pointer=sub_phip::A_judge_query_acc($attr_pointer);
	my @query_accs=@$query_accs_pointer;
	
	#2:
	if($basic_attr{'task'} eq 'overlap'){
		print "###query1: peptides overlapping including all provirome peptides\n";
		sub_phip::B_pep_overlapping($attr_pointer);
	}elsif($basic_attr{'task'} eq 'UniProt'){
		print "###query 2: SwissProt query of provirome\n";
		sub_common::multi_threading(\&sub_bioperl::query_SwissProt, 8, \@query_accs, $basic_attr{task_dir});
	}elsif($basic_attr{'task'} eq 'ProMotifs'){
		print "###query 3: motifs searching\n";
		sub_phip::C_protein_motifs($attr_pointer);
	}elsif($basic_attr{'task'} eq 'export'){
		#combine annotations
		my $comb_pointer=sub_phip::D_combine_annot($attr_pointer, 'UniProt', 'ProMotifs');
		#3: export into text
		my @pro_cols=@{$basic_attr{pro_cols_pointer}};
		my @pro_add = ('RefSeq', 'taxon_id', 'taxon_specie', 'taxon_genus', 'taxon_superkingdom', 
					'gene_symbol', 'gene_synonym', 'product', 'description', 'Pfam', 'EMBL', 'InterPro',);
		push(@pro_cols, @pro_add);
		sub_phip::D_export_fa($comb_pointer, \@pro_cols, $basic_attr{out_dir}.'ProVirome');
	}
}
######################################
#toxome annotations
sub allergome{
	my($basic_attr_pointer)=@_;
	my %basic_attr=%$basic_attr_pointer;

	#1: extract peptides
	my $attr_pointer=sub_phip::A_extract_allergome($basic_attr_pointer);
	my %attr=%$attr_pointer;
	#query accession
	my $query_accs_pointer=sub_phip::A_judge_query_acc($attr_pointer);
	my @query_accs=@$query_accs_pointer;
	
	#2: extract annotations from database
	if($basic_attr{'task'} eq 'overlap'){
		print "###query1: search overlapped peptides\n";
		sub_phip::B_pep_overlapping($attr_pointer);
	}elsif($basic_attr{'task'} eq 'UniProt'){
		print "###query 2: SwissProt query\n";
		sub_common::multi_threading(\&sub_bioperl::query_SwissProt, 8, \@query_accs, $basic_attr{task_dir});
	}elsif($basic_attr{'task'} eq 'ProMotifs'){
		print "###query 3: motifs searching\n";
		sub_phip::C_protein_motifs($attr_pointer);
	}elsif($basic_attr{'task'} eq 'export'){
		#combine annotations
		my $comb_pointer=sub_phip::D_combine_annot($attr_pointer, 'UniProt', 'ProMotifs');
		#3: export into text
		my @pro_cols=@{$basic_attr{pro_cols_pointer}};
		my @pro_add = ('RefSeq', 'taxon_id', 'taxon_specie', 'taxon_genus', 'taxon_superkingdom', 
					'gene_symbol', 'gene_synonym', 'product', 'description', 'Pfam', 'EMBL', 'InterPro',);
		push(@pro_cols, @pro_add);
		sub_phip::D_export_fa($comb_pointer, \@pro_cols,$basic_attr{out_dir}.'Allergome');
	}
}

######################################
#toxome annotations
sub toxome{
	my($basic_attr_pointer)=@_;
	my %basic_attr=%$basic_attr_pointer;

	#1: extract peptides
	my $attr_pointer=sub_phip::A_extract_toxome($basic_attr_pointer);
	#query accession
	my $query_accs_pointer=sub_phip::A_judge_query_acc($attr_pointer);
	my @query_accs=@$query_accs_pointer;
	
	#2: extract annotations from database
	if($basic_attr{'task'} eq 'overlap'){
		print "###query1: search overlapped peptides\n";
		sub_phip::B_pep_overlapping($attr_pointer);
	}elsif($basic_attr{'task'} eq 'UniProt'){
		print "###query 2: SwissProt query\n";
		sub_common::multi_threading(\&sub_bioperl::query_SwissProt, 8, \@query_accs, $basic_attr{task_dir});
	}elsif($basic_attr{'task'} eq 'ProMotifs'){
		print "###query 3: motifs searching\n";
		sub_phip::C_protein_motifs($attr_pointer);
	}elsif($basic_attr{'task'} eq 'export'){
		#combine annotations
		my $comb_pointer=sub_phip::D_combine_annot($attr_pointer, 'UniProt', 'ProMotifs');
		#3: export into text
		my @pro_cols=@{$basic_attr{pro_cols_pointer}};
		my @pro_add = ('RefSeq', 'taxon_id', 'taxon_specie', 'taxon_genus', 'taxon_superkingdom', 
				'gene_symbol', 'gene_synonym', 'product', 'description', 'Pfam', 'EMBL', 'InterPro',);
		push(@pro_cols, @pro_add);
		sub_phip::D_export_fa($comb_pointer, \@pro_cols, $basic_attr{out_dir}.'Toxome');
	}
}

######################################
##virus annotations###
sub virus{
	my($basic_attr_pointer)=@_;
	my %basic_attr=%$basic_attr_pointer;

	
	#1: get basic virus annotations form phip_fa
	my $attr_pointer=sub_phip::A_extract_virus(\%basic_attr);
	my %attr=%$attr_pointer;
	#query accession
	my $query_accs_pointer=sub_phip::A_judge_query_acc($attr_pointer);
	my @query_accs=@$query_accs_pointer;
	
	#################
	if($basic_attr{'task'} eq 'overlap'){
		print "###query1: peptides overlapping including all virus peptides\n";
		sub_phip::B_pep_overlapping($attr_pointer);
	}elsif($basic_attr{'task'} eq 'UniProt'){
		print "###query 2: SwissProt query of virus\n";
		#sub_common::multi_threading(\&sub_bioperl::query_SwissProt, 8, \@query_accs, $basic_attr{task_dir});

		#my @virus_cols = ('taxon_superkingdom',  'taxon_family', 'taxon_genus', 'taxon_specie', );
		#sub_phip::export_annot_table(\@virus_cols, $attr{annot_dir}.'virus_UniProt/', $attr{annot_dir});

		#before run this step: edit personal taxonomy
		sub_phip::personal_taxonomy_v2($attr{annot_dir});
	}elsif($basic_attr{'task'} eq 'InterTaxon'){
		####################
		print "###query 3: add overlapped pep stretches of virus by taxonomy rank\n";
		#family+subfamily for old data
		sub_phip::virusS3_interpep_query($attr_pointer, 'phip_taxon', 'InterTaxon', \@query_accs);
		print "###query 4: intertaxon\n";
		#total counting of overlapped peptides between taxonomic rank
		sub_phip::virusS4_interpep_stat($attr_pointer, 'phip_taxon', 'InterTaxon', \@query_accs);
	}elsif($basic_attr{'task'} eq 'ProMotifs'){
		print "###query 5: motifs searching\n";
		sub_phip::C_protein_motifs($attr_pointer);
	}elsif($basic_attr{'task'} eq 'export'){
		#combine annotations
		my $comb_pointer=sub_phip::D_combine_annot($attr_pointer, 'UniProt', 'InterTaxon', 'ProMotifs');
		
		print "export virus annotations\n";
		my @virus_cols=@{$basic_attr{pro_cols_pointer}};
		my @virus_add = ('RefSeq', 'GeneID', 'taxon_id', 'taxon_specie', 'taxon_genus',  'taxon_family', 
				'taxon_phip', 'InterTaxon', 'description', 'Pfam', 'EMBL', 'InterPro',);
		push(@virus_cols, @virus_add);
		sub_phip::D_export_fa($comb_pointer, \@virus_cols, $basic_attr{out_dir}.'VirScan', $attr{pros_pointer});
		sub_phip::D_update_fa($attr{control_pep_pointer}, \@virus_cols, $basic_attr{out_dir}.'VirScan');
	}

	#match the order peptides between fa, annot.txt, and old matrix.
	#$RC_file='/home/yuan/data_1/results_phip/phipseq5_virus/statistics/lowRC.txt';
	#my $virus_peps_pointer=sub_phip::read_pep_order($RC_file);
	#$fa_file='/home/yuan/phip/ref_seq/T7Pep2_virus.fa';
	#sub_phip::order_pep_fa($virus_peps_pointer, $fa_file);
	#$annot_file='/home/yuan/phip/ref_seq/T7Pep2_virus_annot.txt';
	#sub_phip::order_pep_annot($virus_peps_pointer, $annot_file);

	#
	#my $in_file='/home/yuan/mysql_pre/JHU_projects/input/20160208_ranked_HIV_peptides.csv';
	#my $out_file='/home/yuan/phip/ref_seq/20160208_ranked_HIV_peptides.txt';
	#sub_phip::treat_intra_peptides($in_file, $out_file);

}

######################################
#mouse annotations
sub mouse{
    my($basic_attr_pointer)=@_;
    my %basic_attr=%$basic_attr_pointer;
    print $basic_attr{'task'};
    
    #1: extract peptides
    my $attr_pointer=sub_phip::A_extract_mouse($basic_attr_pointer);
    #query accession
    my $query_accs_pointer=sub_phip::A_judge_query_acc($attr_pointer);
    my @query_accs=@$query_accs_pointer;
    
    #2: extract annotations from database
    if($basic_attr{'task'} eq 'UniProt'){
        print "###query 2: SwissProt query\n";
        sub_common::multi_threading(\&sub_bioperl::query_SwissProt, 8, \@query_accs, $basic_attr{task_dir});
    }elsif($basic_attr{'task'} eq 'ProMotifs'){
        print "###query 3: motifs searching\n";
        sub_phip::C_protein_motifs($attr_pointer);
    }elsif($basic_attr{'task'} eq 'export'){
        #combine annotations
        my $comb_pointer=sub_phip::D_combine_annot($attr_pointer, 'UniProt', 'ProMotifs');
        #3: export into text
        my @pro_cols=@{$basic_attr{pro_cols_pointer}};
        my @pro_add = ('RefSeq', 'taxon_id', 'taxon_specie', 'taxon_genus', 'taxon_superkingdom', 
                'gene_symbol', 'gene_synonym', 'product', 'description', 'Pfam', 'EMBL', 'InterPro',);
        push(@pro_cols, @pro_add);
        sub_phip::D_export_fa($comb_pointer, \@pro_cols, $basic_attr{out_dir}.$basic_attr{'lib'});
    }
}

############################################################################
############################################################################
#main running
############################################################################
############################################################################
#$ARGV[0]=human, virus,toxome,allergome,toxome,mouse
#$ARGV[1]=GenBank, UniProt,EntrezGene,RefGenome, ProMotifs,overlap,PPI, autoantigen, other, export
my $lib= ($ARGV[0]) ? $ARGV[0] : "";
my $task= ($ARGV[1]) ? $ARGV[1] : "";
main($lib, $task);


#others
	#########local_blast
	#print "Local blast\n";
	#my %par=('program'=>'blastp', 'database'=>'/home/yuan/phip/blast/human-representative-cds', 'expect'=>0.01, 
	#			'identity_percent'=>0, 'infile'=>'/home/yuan/phip/ref_seq/VirScan_v1_AA.fa', 
	#			'outfile'=>$out_dir.'VirScan_blast_virus_protein.txt',); 
	#sub_bioperl::local_blast(\%par);
	#sub_phip::collapse_matrix($out_dir.'VirScan_blast_virus_protein.txt', $out_dir.'VirScan_blast_virus_specie.txt', $out_dir.'VirScan_blast_virus_unalignment.txt'); 
	
print "ok\n";