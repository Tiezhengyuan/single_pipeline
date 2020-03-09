#!/usr/bin/perl                                                                                   
use warnings;
use strict;
use Bio::Tools::Run::StandAloneBlast;

my $blast_obj=Bio::Tools::Run::StandAloneBlast->new(-program=>'blastn', 
		-database=>'/home/yuan/phip/blast/VirScan_v1',  -expect => 0.01, 
		-outfile=>'/home/yuan/phip/blast/test.bls');
#input is fasta file
my $seqio=Bio::SeqIO->new(-file=>"/home/yuan/phip/ref_seq/test.fa", -format=>"Fasta");

while(my $seq_obj=$seqio->next_seq()){
	my $report_obj=$blast_obj->blastall($seq_obj);
	while( my $result_obj = $report_obj->next_result ) {     
		my $query_name=$result_obj->query_name;
		my $query_acc=$result_obj->query_accession;
		my $query_desc=$result_obj->query_description;
		while( my $hit= $result_obj->next_hit ) {      
			my $hit_name=$hit->name;
			my $hits_num=$hit->num_hsps;
			printf("Query name: %s\t Query accession:%s\t Description: %s\t Hits_num:%s\n",  $query_name, $query_acc, $query_desc, $hits_num);
			while( my $hsp = $hit->next_hsp ) {     
				my $identity_percent= $hsp->percent_identity;
				if ( $identity_percent > 90) {             
					printf("Hit name:%s\t Length:%s\t Percent_id:%s\t Score:%s\n", $hit_name,$hsp->length('total'),  $identity_percent, $hsp->score);         
				}       
			}
		}
	}
}
print "ok\n";



