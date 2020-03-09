#!/usr/bin/perl
use strict;
use warnings;
use Bio::Perl;
use Bio::SeqIO;
#use Bio::DB::Query::GenBank;
#use IO::String;
use Bio::SeqFeatureI;


#notice: protein.gbk should be revised before running
#%s/gene_synonym="/synonym="/
#%s/db_xref="GeneID:/GeneID="/

#notice: protein.gbk should be revised before running
#%s/gene_synonym="/synonym="/
#%s/db_xref="GeneID:/GeneID="/
#%s/db_xref="miRBase:/miRBase="/
#%s/miRBase="MIMAT/miRBase_MIMAT_No="MIMAT/
#%s/miRBase="MI/miRBase_MI_No="MI/
#%s/db_xref="GI:/GI="/

##############################################################
#extract information from rna.gbk
#The table include 15columns
sub rna_gbk{#1
  my ($variables_pointer)=@_;
  my %variables=%$variables_pointer;
  $variables{result_dir}=$variables{result_dir}.'RNA/';
  mkdir($variables{result_dir}, 0755) unless -d $variables{result_dir};
  
  my %rna_feature;
  open my ($LOG), ">", $variables{result_dir}.$variables{specie}.'_RNA.log' or die;
  print $LOG "\nExtract RNA information from $variables{genome_dir}.$variables{rna_gbk}:\n\n";
  open my ($RNA_txt), ">", $variables{result_dir}.$variables{specie}.'_RNA_info.txt' or die;

  my $in = Bio::SeqIO->new(-file => $variables{rna_gbk}, -format => 'genbank');
  while ( my $seq = $in->next_seq() ) { #2 sequence circyling
    my %rna_info=(display_id=>"NA", accession=>"NA", description=>"NA", rna_seq =>"NA", 
                            rna_len =>0, chromosome =>"NA", rna_map =>"NA", mol_type =>"NA", );
    $rna_info{display_id}=$seq->display_id();
    $rna_info{GI}=$seq->primary_id();
    $rna_info{accession}=$seq->accession();
    $rna_info{description}=$seq->desc();
    $rna_info{rna_seq}=$seq->seq();  
    $rna_info{rna_len}=$seq->length();
    #print "$rna_info{accession}\n";

    foreach my $feat ( $seq->get_SeqFeatures() ) { #3 features circyling
        my %feature;
        my $tag_name=$feat->primary_tag();
        $feature{name}=$tag_name;
        $feature{start}=$feat->start;
        $feature{end}=$feat->end;
        $feature{strand}=$feat->strand;
        $feature{seq}=substr($rna_info{rna_seq}, $feature{start}-1, $feature{end}-$feature{start}+1);
        $feature{gene}='NA';
        #get feature's info
         foreach my $tag ( $feat->get_all_tags() ) {  #4
           my $tag_value=join("", $feat->get_tag_values($tag));
            if($tag=~/db_xref/){
              my($a, $b)=split(":", $tag_value);
              $feature{$a}=$b;
            }
            else{
                $feature{$tag}=$tag_value;
            }
            
         }#4
         
         #get pseudo-genes
         if ($tag_name=~/misc_RNA/ and exists $feature{pseudo}){
            $feature{name}='pseudo';
        }
        #extract source_features
        if ( $tag_name=~/source/ ) { #4
          $rna_info{chromosome}=$feature{chromosome} if exists $feature{chromosome};
          $rna_info{rna_map}=$feature{map} if exists $feature{map};
          $rna_info{mol_type}=$feature{mol_type} if exists $feature{mol_type};
          #rna sequences information and source features 
          print $RNA_txt join("\t", $rna_info{accession}, $rna_info{display_id}, $rna_info{description}, 
                                          $rna_info{rna_len}, $rna_info{chromosome}, $rna_info{rna_map}, 
                                          $rna_info{mol_type}, $rna_info{rna_seq}), "\n";
          #change feature name
          $feature{mol_type}=~s/\s/_/g;
          $feature{name}=$feature{mol_type};
        }#4
        elsif ($tag_name=~/ncRNA/) {#4
          $feature{name}=$feature{ncRNA_class};
        }#4
        my $feature_index=$rna_info{accession}.':'.$rna_info{rna_len}.':'.$feature{gene}.':'.$feature{start}.'_'.$feature{end};
        $rna_feature{$feature{name}}->{$feature_index}=$feature{seq};
        
        if ($tag_name=~/CDS/){#4
          #5'UTR
          if($feature{start}>10){
            my $UTR5_start=1; 
            my $UTR5_end=$feature{start}-1;
            my $UTR5_seq=substr($rna_info{rna_seq}, $UTR5_start-1, $UTR5_end-$UTR5_start+1);
            my $feature_index=$rna_info{accession}.':'.$rna_info{rna_len}.':'.$feature{gene}.':'.$UTR5_start.'_'.$UTR5_end; 
            $rna_feature{'5UTR'}->{$feature_index}=$UTR5_seq;
          }
          #3'UTR
          if ($rna_info{rna_len}-$feature{end}>10){
             my $UTR3_start=$feature{end}+1;
             my $UTR3_end=$rna_info{rna_len};
             my $UTR3_seq=substr($rna_info{rna_seq}, $UTR3_start-1, $UTR3_end-$UTR3_start+1);
             my $feature_index=$rna_info{accession}.':'.$rna_info{rna_len}.':'.$feature{gene}.':'.$UTR3_start.'_'.$UTR3_end; 
             $rna_feature{'3UTR'}->{$feature_index}=$UTR3_seq;
          }
       }#4
     }#3 features circyling
  } #2 sequence circyling
  
  print "export sequences of features\n";
  foreach my $feature_name(keys %rna_feature){
    print "$feature_name\n";
    my $pointer=$rna_feature{$feature_name};
    my %hash=%$pointer;
    open my $OUT, ">", $variables{result_dir}.$variables{specie}.'_RNA_'.$feature_name.'.fa' or die; 
    foreach my $feature_index(keys %hash){
      print $OUT ">$feature_name:$feature_index\n", "$hash{$feature_index}\n";
    }
    close($OUT);
    my $num=keys %hash;
    print $LOG "$feature_name=$num\n";
  }

  close($LOG);
}#1



############################################################
#extract information from protein.gbk
#The table include 15columns
sub protein_gbk{#1
  my ($variables_pointer)=@_;
  my %variables=%$variables_pointer;

  my $num=0;
  my $CDS_num=0;
  my $protein_num=0;
  print "\nExtract protein information from $variables{genome_dir}.$variables{protein_gbk}:\n\n";
  my $in = Bio::SeqIO->new(-file => $variables{genome_dir}.$variables{protein_gbk}, -format => 'genbank');
  while ( my $seq = $in->next_seq() ) { #2
    my %pro_info=(
      display_id=>"NA", GI=>"NA", accession=>"NA", secondary_acc=>"NA", 
      description=>"NA", pro_seq=>"NA", pro_len=>0, chromosome=>"NA", 
      pro_map        =>"NA",      pro_note       =>"NA",       pro_product    =>"NA",       EC_number      =>"NA", 
      pro_weight     =>0, pro_start      =>0,      pro_end        =>0,
      gene           =>"NA",       gene_synonym   =>"NA",       coded_by       =>"NA",       geneid         =>"NA",
      CDS_note       =>"NA",       CDS_start      =>0,      CDS_end        =>0,
        );

    $num++;    
    $pro_info{display_id}=$seq->display_id();
    $pro_info{GI}=$seq->primary_id();
    $pro_info{accession}=$seq->accession();
    $pro_info{description}=$seq->desc();
    $pro_info{pro_seq}=$seq->seq();  #protein sequence
    $pro_info{pro_len}=$seq->length();
    my @secondary_acc=$seq->get_secondary_accessions();   
    $pro_info{secondary_acc}=join(",", @secondary_acc) if @secondary_acc>0;
    
    my (@chromosome, @map, @product, @pro_note, @CDS_note, @MW, @gene, @gene_synonym, @coded_by, @geneid, @EC_number);
    foreach my $feat ( $seq->get_SeqFeatures() ) { #3
      my $tag_name=$feat->primary_tag();
      #extract source_features
      if ( $tag_name=~/Source/ ) { #4
        foreach my $tag ( $feat->get_all_tags() ) {  #5
          @chromosome = $feat->get_tag_values($tag) if $tag=~/chromosome/;
          @map = $feat->get_tag_values($tag) if $tag=~/map/;
        }#5
        $pro_info{chromosome}=join("", @chromosome) if @chromosome>0;
        $pro_info{pro_map}=join("", @map) if @map>0;
      }#4
           
      #extract protein_features
      if ( $tag_name=~/Protein/ ) { #4
        foreach my $tag ( $feat->get_all_tags() ) {  #5
          @product = $feat->get_tag_values($tag) if $tag=~/product/;
          @pro_note = $feat->get_tag_values($tag) if $tag=~/note/;
          @MW = $feat->get_tag_values($tag) if $tag=~/calculated_mol_wt/;
          @EC_number = $feat->get_tag_values($tag) if $tag=~/EC_number/;
          $pro_info{pro_start}=$feat->start;
          $pro_info{pro_end}=$feat->end;
                }#5
       $protein_num++;
       $pro_info{pro_product}=join("", @product) if @product>0; 
       $pro_info{pro_note}=join("", @pro_note) if @pro_note>0;
       $pro_info{pro_weight}=$MW[0] if @MW>0;
       $pro_info{EC_number}=join("", @EC_number) if @EC_number>0;      
           }#4
 
     #extract CDS_features
     if ( $tag_name=~/CDS/ ) { #4
        foreach my $tag ( $feat->get_all_tags() ) {  #5
          @gene = $feat->get_tag_values($tag) if $tag=~/gene/;
          @gene_synonym = $feat->get_tag_values($tag) if $tag=~/synonym/;
          @coded_by = $feat->get_tag_values($tag) if $tag=~/coded_by/;
          @geneid = $feat->get_tag_values($tag) if $tag=~/GeneID/;
          @CDS_note = $feat->get_tag_values($tag) if $tag=~/note/;
          $pro_info{CDS_start}=$feat->start;
          $pro_info{CDS_end}=$feat->end;
                }#5
       $CDS_num++; 
       $pro_info{gene}=join("", @gene) if @gene>0;
       $pro_info{gene_synonym}=join("", @gene_synonym) if @gene_synonym>0;  
       $pro_info{gene_synonym}=~s/;\s/,/g;  
       $pro_info{CDS_note}=join("", @CDS_note) if @CDS_note>0;
       $pro_info{coded_by}=join("", @coded_by) if @coded_by>0;
       $pro_info{geneid}=join("", @geneid) if @geneid>0;      
           }#4        
         
        }#3

   print "$num:\t", "protein=$protein_num\t", "CDS=$CDS_num\t", "$pro_info{gene_synonym}\n";
   #summary features
   print "$pro_info{accession}\t", "$pro_info{description}\t", "$pro_info{display_id}\t", "$pro_info{GI}\t", "$pro_info{secondary_acc}\t", "$pro_info{pro_len}\t", "$pro_info{pro_seq}\t"; 
   #source features 
   print "$pro_info{chromosome}\t", "$pro_info{pro_map}\t"; 
   #protein features
   print "$pro_info{EC_number}\t", "$pro_info{pro_start}\t", "$pro_info{pro_end}\t", "$pro_info{pro_note}\t", "$pro_info{pro_product}\t", "$pro_info{pro_weight}\t";
   #CDS features
   print  "$pro_info{CDS_note}\t", "$pro_info{coded_by}\t", "$pro_info{gene}\t", "$pro_info{gene_synonym}\t", "$pro_info{geneid}\t";
   } #2 
   
   print "The items number: $num\n";   
   print "The number of proteins: $protein_num\n";
   print "The number of CDS regions: $CDS_num\n";
}#1


############################################################
#extract information from genome.gbk
sub genome_gbk{#1
  my $variables_pointer=$_[0];
  my %variables=%$variables_pointer;
  my $chromosomes_pointer=$_[1];
  my @chromosomes=@$chromosomes_pointer;

  my (%num_out, @feature_names);
  open my ($OUT), ">", $variables{result_dir}.$variables{genome_log} or die;
  print $OUT "\nExtract genome information:\n\n";
  open my ($OUT1), ">", $variables{result_dir}.$variables{genome} or die;
  open my ($OUT2), ">", $variables{result_dir}.$variables{genome_feature} or die;
  open my ($OUT3), ">", $variables{result_dir}.$variables{genome_fasta} or die;
  open my ($OUT4), ">", $variables{result_dir}.$variables{genome_feature_fasta} or die;
  open my ($OUT5), ">", $variables{result_dir}.$variables{genome_gene_fasta} or die;
  open my ($OUT6), ">", $variables{result_dir}.$variables{genome_exon_fasta} or die;
  foreach my $chr_file(@chromosomes){#2 chromosome circyling
    my $chr_name=substr($chr_file, 0, 6); 
    my $in_obj = Bio::SeqIO->new(-file => $variables{genome_dir}.$chr_file, -format => 'genbank');
    while ( my $seq = $in_obj->next_seq() ) { #3 gbk circyling
      my %info=( display_id=>"NA",  accession=>"NA",  secondary_acc=>"NA", description=>"NA", sequence=>"NA",
                 seq_len=>0, GI=>"NA", chromosome=>"NA", mol_type=>"NA");
      $info{display_id}=$seq->display_id();
      $info{GI}=$seq->primary_id();
      $info{accession}=$seq->accession();
      $info{description}=$seq->desc();
      $info{sequence}=$seq->seq();  
      $info{seq_len}=$seq->length();
      $info{secondary_acc}=join(",", $seq->get_secondary_accessions());

      #extract feature information   
      foreach my $feat ( $seq->get_SeqFeatures() ) { #4 feature circyling
        my %tag_info=(tag_name =>"NA", tag_start=>0, tag_end=>0, tag_strand=>"NA",
                      tag_seq=>"NA", tag_seq_len=>0, tag_exon_seq=>"NA", tag_exon_len=>0,    
                      tag_gene=>"NA", tag_synonym=>"NA", tag_product=>"NA", tag_note=>"NA", tag_ncRNA_class=>"NA",
                      tag_transcript_id=>"NA", tag_db_xref=>"NA", tag_protein_id=>"NA", tag_gene_id=>"NA", tag_GI=>"NA",
                      );
      
        $tag_info{tag_name}=$feat->primary_tag();
        push(@feature_names, $tag_info{tag_name}) unless List::Util::first {$tag_info{tag_name} eq $_} @feature_names;
        $num_out{$chr_name}->{$tag_info{tag_name}}++;
        $num_out{'total'}->{$tag_info{tag_name}}++;

        $tag_info{tag_start}=$feat->start;
        $tag_info{tag_end}=$feat->end;
        $tag_info{tag_strand}=$feat->strand;
        $tag_info{tag_seq}=$feat->seq()->seq; # intron and exon
        $tag_info{tag_seq_len}=$feat->seq()->length;
        $tag_info{tag_exon_seq}=$feat->spliced_seq()->seq;  #exon only
        $tag_info{tag_exon_len}=$feat->spliced_seq()->length;
        #print "$len\t$tag_info{tag_seq_len}\t$tag_info{tag_exon_len}\n";
      
        if ($tag_info{tag_name}=~/source/) {#5
          foreach my $tag ( $feat->get_all_tags() ) {  #6
            $info{chromosome}=join("", $feat->get_tag_values($tag)) if $tag=~/chromosome/;
            $info{mol_type}=join("", $feat->get_tag_values($tag)) if $tag=~/mol_type/;
          }#6
          print $OUT1 "$info{accession}\t", "$info{display_id}\t", "$info{GI}\t", "$info{secondary_acc}\t", "$info{chromosome}\t", "$info{mol_type}\t", "$info{seq_len}\t", "$info{description}\t", "$info{sequence}\n";
          print $OUT3 ">genomeDNA:", "chr$info{chromosome}#", "$info{accession}\n", "$info{sequence}\n"; 
        }#5
        else{#5
          foreach my $tag ( $feat->get_all_tags() ) {  #6
            $tag_info{tag_note}=join("", $feat->get_tag_values($tag)) if $tag=~/note/;
            $tag_info{tag_gene}=join("", $feat->get_tag_values($tag)) if $tag eq "gene";
            $tag_info{tag_synonym}=join("", $feat->get_tag_values($tag)) if $tag eq "gene_synonym";
            $tag_info{tag_product}=join("", $feat->get_tag_values($tag)) if $tag=~/product/;
            $tag_info{tag_protein_id}=join("", $feat->get_tag_values($tag)) if $tag=~/protein_id/;
            $tag_info{tag_transcript_id}=join("", $feat->get_tag_values($tag)) if $tag=~/transcript_id/;
            $tag_info{tag_ncRNA_class}=join("", $feat->get_tag_values($tag)) if $tag=~/ncRNA_class/;
            if ($tag=~/db_xref/){#7
              my @db_xref=$feat->get_tag_values($tag);
              $tag_info{tag_db_xref}=join(",", @db_xref);
              foreach (@db_xref){
                $tag_info{tag_gene_id}=$_ if $_=~/GeneID/;
                $tag_info{tag_GI}=$_ if $_=~/GI/;
              }
            }#7
          }#6 

          #tag features 
          print $OUT2 "$info{accession}\t", "$tag_info{tag_name}\t", "$tag_info{tag_start}\t", "$tag_info{tag_end}\t", "$tag_info{tag_seq_len}\t", "$tag_info{tag_exon_len}\t", "$tag_info{tag_strand}\t", "$tag_info{tag_gene}\t", "$tag_info{tag_synonym}\t", "$tag_info{tag_transcript_id}\t", "$tag_info{tag_protein_id}\t", "$tag_info{tag_gene_id}\t", "$tag_info{tag_GI}\t", "$tag_info{tag_note}\t", "$tag_info{tag_product}\t", "$tag_info{tag_db_xref}\t", "$tag_info{tag_ncRNA_class}\t", "$tag_info{tag_seq}\t", "$tag_info{tag_exon_seq}\n";
          print "$info{accession}\t", "$tag_info{tag_name}\t", "$tag_info{tag_start}\t", "$tag_info{tag_end}\t", "$tag_info{tag_seq_len}\t", "$tag_info{tag_exon_len}\t",  "$tag_info{tag_gene}\t", "$tag_info{tag_ncRNA_class}\n";
print "$tag_info{tag_seq_len}\n" if $tag_info{tag_seq_len}>8388607;
          print $OUT4 ">$tag_info{tag_name}", ":chr$info{chromosome}#", "$info{accession}#", "$tag_info{tag_gene}\n", "$tag_info{tag_seq}\n";
          print $OUT5 ">$tag_info{tag_name}", ":chr$info{chromosome}#", "$info{accession}#", "$tag_info{tag_gene}\n", "$tag_info{tag_seq}\n" if $tag_info{tag_name} eq "gene";
          print $OUT6 ">$tag_info{tag_name}", ":chr$info{chromosome}#", "$info{accession}#", "$tag_info{tag_gene}\n", "$tag_info{tag_seq}\n" if $tag_info{tag_name} eq "exon";

        }#5
      }#4 feature circyling

      print "$chr_name:\t";
      foreach (@feature_names){
        $num_out{$chr_name}->{$_}=0 unless exists $num_out{$chr_name}->{$_};
        print "$_=$num_out{$chr_name}->{$_}\t";
      }
      print "\n";

    } #3 gbk circyling
  }#2 chromosome circyling

  @feature_names=sort @feature_names;
  push(@feature_names, "total");
  my $out_title=join("\t", "chromosome", @feature_names);
  print $OUT "$out_title\n";
  foreach my $chr(sort(keys %num_out)) {#2
    print $OUT "$chr\t";
    foreach (@feature_names){
      $num_out{$chr}->{$_}=0 unless exists $num_out{$chr}->{$_};
      print $OUT "$_=$num_out{$chr}->{$_}\t";
    }
    print $OUT "\n";
  } #2
  close($OUT1);   
  close($OUT2); 
  close($OUT3); 
  close($OUT4); 
  close($OUT5); 
  close($OUT6);
  close($OUT);

  return(\%num_out);
}#1

############################################################
#extract assembled genome information from *.gbs and *.fa
sub assembled_genome{#1
  my $variables_pointer=$_[0];
  my %variables=%$variables_pointer;
  my $chromosomes_pointer=$_[1];
  my @chromosomes=@$chromosomes_pointer;

  my (%num_out, @feature_names);
  open my ($OUT), ">", $variables{result_dir}.$variables{genome_log} or die;
  print $OUT "\nExtract genome information:\n\n";
  open my ($OUT1), ">", $variables{result_dir}.$variables{genome} or die;
  open my ($OUT2), ">", $variables{result_dir}.$variables{genome_feature} or die;
  open my ($OUT3), ">", $variables{result_dir}.$variables{genome_fasta} or die;
  open my ($OUT4), ">", $variables{result_dir}.$variables{genome_feature_fasta} or die;
  open my ($OUT5), ">", $variables{result_dir}.$variables{genome_gene_fasta} or die;
  open my ($OUT6), ">", $variables{result_dir}.$variables{genome_exon_fasta} or die;
  foreach my $chr_name(@chromosomes){#2 chromosome circyling

    #extract information from fasta file
    my %info=( accession=>"NA",  description=>"NA", sequence=>"NA",  seq_len=>0, GI=>"NA", chromosome=>"NA",);
    my $seqio_obj = Bio::SeqIO->new(-file => $variables{genome_dir}."Assembled_chromosomes/seq/hs_ref_GRCh37.p10_".$chr_name.".fa", -format => 'fasta');
    my $seq_obj=$seqio_obj->next_seq();
    my @displayid=split(/\|/, $seq_obj->display_id());
    $info{GI}=$displayid[1];
    $info{accession}=$displayid[3];
    $info{description}=$seq_obj->desc();
    $info{sequence}=$seq_obj->seq();  
    $info{seq_len}=$seq_obj->length();

    #extract information from fasta file
    my $in_obj = Bio::SeqIO->new(-file => $variables{genome_dir}."Assembled_chromosomes/gbs/hs_ref_GRCh37.p10_".$chr_name.".gbs", -format => 'genbank');
    while ( my $seq = $in_obj->next_seq() ) { #3 gbs circyling


      #extract feature information   
      foreach my $feat ( $seq->get_SeqFeatures() ) { #4 feature circyling
        my %tag_info=(tag_name =>"NA", tag_start=>0, tag_end=>0, tag_strand=>"NA",
                      tag_seq=>"NA", tag_seq_len=>0, tag_exon_seq=>"NA", tag_exon_len=>0,    
                      tag_gene=>"NA", tag_synonym=>"NA", tag_product=>"NA", tag_note=>"NA", tag_ncRNA_class=>"NA",
                      tag_transcript_id=>"NA", tag_db_xref=>"NA", tag_protein_id=>"NA", tag_gene_id=>"NA", tag_GI=>"NA",
                      );
      
        $tag_info{tag_name}=$feat->primary_tag();
        push(@feature_names, $tag_info{tag_name}) unless List::Util::first {$tag_info{tag_name} eq $_} @feature_names;
        $num_out{$chr_name}->{$tag_info{tag_name}}++;
        $num_out{'total'}->{$tag_info{tag_name}}++;

        $tag_info{tag_start}=$feat->start;
        $tag_info{tag_end}=$feat->end;
        $tag_info{tag_strand}=$feat->strand;
        $tag_info{tag_seq}=substr($info{sequence}, $tag_info{tag_start}, $tag_info{tag_end}-$tag_info{tag_start}); 
        $tag_info{tag_seq_len}=length($tag_info{tag_seq});

        print "$tag_info{tag_strand}\t$tag_info{tag_start}:$tag_info{tag_end}\n\n\n" if $tag_info{tag_end}<$tag_info{tag_start};
      
        if ($tag_info{tag_name}=~/source/) {#5
          foreach my $tag ( $feat->get_all_tags() ) {  #6
            $info{chromosome}=join("", $feat->get_tag_values($tag)) if $tag=~/chromosome/;
            $info{mol_type}=join("", $feat->get_tag_values($tag)) if $tag=~/mol_type/;
          }#6
          print $OUT1 "$info{accession}\t", "$info{GI}\t", "$info{chromosome}\t", "$info{mol_type}\t", "$info{seq_len}\t", "$info{description}\t", "$info{sequence}\n";
          print $OUT3 ">assembledDNA:", "chr$info{chromosome}#", "$info{accession}\n", "$info{sequence}\n"; 
        }#5
        else{#5
          foreach my $tag ( $feat->get_all_tags() ) {  #6
            $tag_info{tag_note}=join("", $feat->get_tag_values($tag)) if $tag=~/note/;
            $tag_info{tag_gene}=join("", $feat->get_tag_values($tag)) if $tag eq "gene";
            $tag_info{tag_synonym}=join("", $feat->get_tag_values($tag)) if $tag eq "gene_synonym";
            $tag_info{tag_product}=join("", $feat->get_tag_values($tag)) if $tag=~/product/;
            $tag_info{tag_protein_id}=join("", $feat->get_tag_values($tag)) if $tag=~/protein_id/;
            $tag_info{tag_transcript_id}=join("", $feat->get_tag_values($tag)) if $tag=~/transcript_id/;
            $tag_info{tag_ncRNA_class}=join("", $feat->get_tag_values($tag)) if $tag=~/ncRNA_class/;
            if ($tag=~/db_xref/){#7
              my @db_xref=$feat->get_tag_values($tag);
              $tag_info{tag_db_xref}=join(",", @db_xref);
              foreach (@db_xref){
                $tag_info{tag_gene_id}=$_ if $_=~/GeneID/;
                $tag_info{tag_GI}=$_ if $_=~/GI/;
              }
            }#7
          }#6 

          #tag features 
          print $OUT2 "$info{accession}\t", "$tag_info{tag_name}\t", "$tag_info{tag_start}\t", "$tag_info{tag_end}\t", "$tag_info{tag_seq_len}\t", "$tag_info{tag_exon_len}\t", "$tag_info{tag_strand}\t", "$tag_info{tag_gene}\t", "$tag_info{tag_synonym}\t", "$tag_info{tag_transcript_id}\t", "$tag_info{tag_protein_id}\t", "$tag_info{tag_gene_id}\t", "$tag_info{tag_GI}\t", "$tag_info{tag_note}\t", "$tag_info{tag_product}\t", "$tag_info{tag_db_xref}\t", "$tag_info{tag_ncRNA_class}\t", "$tag_info{tag_seq}\t", "$tag_info{tag_exon_seq}\n";
          print $OUT4 ">$tag_info{tag_name}", ":chr$info{chromosome}#", "$info{accession}#", "$tag_info{tag_gene}\n", "$tag_info{tag_seq}\n";
          print $OUT5 ">$tag_info{tag_name}", ":chr$info{chromosome}#", "$info{accession}#", "$tag_info{tag_gene}\n", "$tag_info{tag_seq}\n" if $tag_info{tag_name} eq "gene";
          print $OUT6 ">$tag_info{tag_name}", ":chr$info{chromosome}#", "$info{accession}#", "$tag_info{tag_gene}\n", "$tag_info{tag_seq}\n" if $tag_info{tag_name} eq "exon";

        }#5
      }#4 feature circyling

      print "$chr_name:\t";
      foreach (@feature_names){
        $num_out{$chr_name}->{$_}=0 unless exists $num_out{$chr_name}->{$_};
        print "$_=$num_out{$chr_name}->{$_}\t";
      }
      print "\n";

    } #3 gbs circyling

  }#2 chromosome circyling

  @feature_names=sort @feature_names;
  push(@feature_names, "total");
  my $out_title=join("\t", "chromosome", @feature_names);
  print $OUT "$out_title\n";
  foreach my $chr(sort(keys %num_out)) {#2
    print $OUT "$chr\t";
    foreach (@feature_names){
      $num_out{$chr}->{$_}=0 unless exists $num_out{$chr}->{$_};
      print $OUT "$_=$num_out{$chr}->{$_}\t";
    }
    print $OUT "\n";
  } #2
  close($OUT1);   
  close($OUT2); 
  close($OUT3); 
  close($OUT4); 
  close($OUT5); 
  close($OUT6);
  close($OUT);

  return(\%num_out);
}#1

#######################
sub gnomon_mRNA{
  my $variables_pointer=$_[0];
  my %variables=%$variables_pointer;

  open my($OUT), ">", "/home/yuan/mysql_pre/human/RNA/human_predicted_mRNA.fa" or die;
  my $in_obj = Bio::SeqIO->new(-file => $variables{gnomon_mRNA_fa}, -format => 'fasta');
  while ( my $seq_obj = $in_obj->next_seq() ) { #3 gbs circyling
    my $seq=$seq_obj->seq();
    my $displayid=$seq_obj->display_id();
    $displayid=~s/gnl\|//;
    $displayid=~s/\.m//;
    $displayid=~s/\|/_/;
    print $OUT ">$displayid\n", "$seq\n";
  }
  close($OUT)
}
########
sub FLJ_cDNA{
  my $variables_pointer=$_[0];
  my %variables=%$variables_pointer;

  open my($OUT), ">", "/home/yuan/mysql_pre/human/RNA/human_FLJ_cDNA.fa" or die;
  my $in_obj = Bio::SeqIO->new(-file => $variables{FLJ_fa}, -format => 'fasta');
  while ( my $seq_obj = $in_obj->next_seq() ) { #3 gbs circyling
    my $seq=$seq_obj->seq();
    my $displayid=$seq_obj->display_id();
    print $OUT ">$displayid\n", "$seq\n";
    #print "$displayid\n";
  }
  close($OUT)
}
#######################################################
sub LncRNA{#1
  my $variables_pointer=$_[0];
  my %variables=%$variables_pointer;

  my $num=0;
  my $in_obj = Bio::SeqIO->new(-file => $variables{LncRNA}, -format => 'fasta');
  open my($OUT1), ">", "/home/yuan/mysql_pre/EmiRNA/bowtie/human_lncRNA.fa" or die;
  while ( my $seq = $in_obj->next_seq() ) { #2
    my ($nc_id, $acc, $lnc_class, $species, $lnc_name, $lnc_ref, $transcript_id, $url, $cpcScore, $cnci)=split(",", $seq->display_id());
    my $lnc_seq=$seq->seq();
    print $OUT1 ">lncRNA:$acc\n", "$lnc_seq\n";
    $num++;
    print "Lnc_RNA=$num\n";
  }#2
  close($OUT1);
}#1


#######################################################
sub GenBank_RNA{#1
  my $variables_pointer=$_[0];
  my %variables=%$variables_pointer;

  my (@RNAs, %RNA_hash);
  open my($OUT1), ">", $variables{result_dir}.$variables{query_RNA_specie}.'.txt' or die;
  my $query_obj = Bio::DB::Query::GenBank->new(-db => "nucleotide", -query => $variables{db_query});
  my $gb_obj=Bio::DB::GenBank->new;
  my $stream_obj=$gb_obj->get_Stream_by_query($query_obj);
  while ( my $seq_obj = $stream_obj->next_seq() ) { #2
    my %info=( display_id=>"NA",  display_name=>"NA", accession=>"NA",  secondary_acc=>"NA", description=>"NA", 
                        sequence=>"NA", seq_len=>0, GI=>"NA", ncRNA_class=>"NA", product=>"NA");

    $info{display_id}=$seq_obj->display_id;
    $info{display_name}=$seq_obj->display_name;
    $info{GI}=$seq_obj->primary_id;
    $info{accession}=$seq_obj->accession;
    $info{secondary_acc}=$seq_obj->get_secondary_accessions if $seq_obj->get_secondary_accessions != 0;
    $info{description}=$seq_obj->desc;
    $info{sequence}=$seq_obj->seq;
    my $seq=$info{sequence};
    $info{seq_len}=$seq_obj->length;
    foreach my $feat ( $seq_obj->get_SeqFeatures() ) { #3
      my $tag_name=$feat->primary_tag();
      if ($tag_name=~/$variables{query_RNA_specie}/) {#4
        foreach my $tag ( $feat->get_all_tags() ) {  #5
                print "$tag_name\t$tag\t$info{description}\n";
          if ($tag=~/product/){#6
            $info{product}=join("", $feat->get_tag_values($tag)) ;
            push (@RNAs, $info{product}) unless List::Util::first {$info{product} eq $_} @RNAs;
            if ($RNA_hash{$seq}) {  $RNA_hash{$seq} .= ",".$info{accession};    }
            else {  $RNA_hash{$seq}=$info{product}.','.$info{accession}; }
          }#6
        }#5

        print $OUT1 join("\t", $info{display_id}, $info{display_name}, $info{GI}, $info{accession}, $info{secondary_acc}, 
                                     $info{product}, $info{description}, $info{seq_len}, $info{sequence}), "\n"; 
      }#4

    }#3
  }#2
  close($OUT1);
  my $num=@RNAs;
  print "the number of RNAs from GeneBank is $num.\n";

  my $uniq_num=0;
  open my($OUT2), ">", $variables{result_dir}.$variables{query_RNA_specie}.'.fa' or die;
  foreach my $uniq_seq(keys %RNA_hash){
    print $OUT2 ">".$variables{query_RNA_specie}.":$RNA_hash{$uniq_seq}\n", "$uniq_seq\n";
    $uniq_num++;
  }
  close($OUT2);
  print "the unique number of RNAs is $uniq_num.\n";
}#1

################
sub siRNA{
  my $variables_pointer=$_[0];
  my %variables=%$variables_pointer;

  my (%seq_hash, @len);
  open my($IN), "<", $variables{result_dir}.'sirnadb_050915.txt' or die;
  while(<$IN>){
    chomp($_);
    my @items=split("\t", $_);
    my $seq=$items[1];
    my $accession=$items[3];
    if (exists $seq_hash{$seq}){
      $seq_hash{$seq} .= ','.$accession;
    }
    else{
        $seq_hash{$seq} = $accession;
    }
    push(@len, length($seq));
    print "$items[1]\t$items[3]\n";
  }
  close($IN);
  my $ave_len=(List::Util::sum @len) /@len;
print "Length on average:$ave_len\n";
  #
  open my($OUT), ">", $variables{result_dir}.'human_siRNA.fa' or die;
  foreach(keys %seq_hash){
      print $OUT ">siRNA:$seq_hash{$_}\n", "$_\n";
  }
  close($OUT);
}

######################################################
#main program

#H. sapiens genome data extraction
my %variables=(
   'genome_dir'=>"/home/yuan/data_preparation/ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/",
   'result_dir'=>"/home/yuan/mysql_pre/human/", 'specie'=>"human",
   'rna_gbk'=>"/home/yuan/data_preparation/ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/RNA/rna.gbk", 
   'protein_gbk'=>"/home/yuan/data_preparation/ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/protein/protein.gbk",
   'gnomon_mRNA_fa'=>"/home/yuan/data_preparation/ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/RNA/Gnomon_mRNA.fsa", 
   'FLJ_fa'=>"/home/yuan/mysql_pre/human/V3_FLJ_FLJPJ_nq.fasta",
);


# extract protein information from protein.gbk
#protein_gbk(\%variables);


# extract rna information from rna.gbk
#rna_gbk(\%variables);

#predicted transcripts
#gnomon_mRNA(\%variables);

#FLJ cDNA
#FLJ_cDNA(\%variables);


#extract genome information
my @chromosomes=qw(
CHR_01/hs_ref_GRCh37.p10_chr1.gbk  CHR_10/hs_ref_GRCh37.p10_chr10.gbk  CHR_19/hs_ref_GRCh37.p10_chr19.gbk
CHR_02/hs_ref_GRCh37.p10_chr2.gbk  CHR_11/hs_ref_GRCh37.p10_chr11.gbk  CHR_20/hs_ref_GRCh37.p10_chr20.gbk
CHR_03/hs_ref_GRCh37.p10_chr3.gbk  CHR_12/hs_ref_GRCh37.p10_chr12.gbk  CHR_21/hs_ref_GRCh37.p10_chr21.gbk
CHR_04/hs_ref_GRCh37.p10_chr4.gbk  CHR_13/hs_ref_GRCh37.p10_chr13.gbk  CHR_22/hs_ref_GRCh37.p10_chr22.gbk
CHR_05/hs_ref_GRCh37.p10_chr5.gbk  CHR_14/hs_ref_GRCh37.p10_chr14.gbk  CHR_MT/hs_ref_GRCh37.p10_chrMT.gbk
CHR_06/hs_ref_GRCh37.p10_chr6.gbk  CHR_15/hs_ref_GRCh37.p10_chr15.gbk  CHR_Un/hs_ref_GRCh37.p10_chrUn.gbk
CHR_07/hs_ref_GRCh37.p10_chr7.gbk  CHR_16/hs_ref_GRCh37.p10_chr16.gbk  CHR_X/hs_ref_GRCh37.p10_chrX.gbk
CHR_08/hs_ref_GRCh37.p10_chr8.gbk  CHR_17/hs_ref_GRCh37.p10_chr17.gbk  CHR_Y/hs_ref_GRCh37.p10_chrY.gbk
CHR_09/hs_ref_GRCh37.p10_chr9.gbk  CHR_18/hs_ref_GRCh37.p10_chr18.gbk   );

#genome_gbk(\%variables, \@chromosomes);


#extract assembled_genome information
my @assembled_chromosomes=qw(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY);
$variables{genome_log}="human_assembled_genome.log";
$variables{genome}="human_assembled_genome_DNA.txt";
$variables{genome_feature}="human_assembled_genome_feature.txt",  
$variables{genome_fasta}="human_assembled_genome_DNA.fasta",
$variables{genome_feature_fasta}="human_assembled_genome_feature.fasta",
$variables{genome_gene_fasta}="human_assembled_genome_gene.fasta",
$variables{genome_exon_fasta}="human_assembled_genome_exon.fasta",
#assembled_genome(\%variables, \@assembled_chromosomes);

#extract LncRNAs information for NONCODE database
$variables{LncRNA}="/data/NONCODE/human/human.fa";

#Homo sapiens lncRNAs:From NONCODE 3.0 length >= 200
#Scripture Reconstruction LincRNAs By Gutterman [ Reference PMID:20436462,19182780 ]
#Scripture Reconstruction LincRNAs By Luo [ unpublished ]
#>ncid | accn | class | organism | name | ref | transcriptID | url | cpcScore | cnci
#sequence
#LncRNA(\%variables);


#$variables{db_query}="Homo sapiens[organism] AND piRNA[TITL] AND 0:5000[SLEN]";
#$variables{query_RNA_specie}='piRNA';
#GenBank_RNA(\%variables);

#tRNA GenBank
#$variables{db_query}="Homo sapiens[organism] AND tRNA[keywords] AND 0:500[SLEN]";
#$variables{query_RNA_specie}='tRNA';
#GenBank_RNA(\%variables);

#tRNA GenBank
#siRNA(\%variables);


###############33
#maize 

# extract rna information from rna.gbk
%variables=(#'genome_dir'=>'/home/yuan/backup_1/data_preparation/ftp.ncbi.nlm.nih.gov/genomes/Zea_mays/', 
    'result_dir'=>'/home/yuan/data_preparation/NCBI/maize/', 'specie'=>'maize', 
    'rna_gbk'=>'/home/yuan/data_preparation/NCBI/maize/rna.gbk');
#rna_gbk(\%variables);

#piwiRNA
$variables{db_query}="Zea mays[organism] AND piRNA[TITL] AND 0:5000[SLEN]";
$variables{query_RNA_specie}='piRNA';
GenBank_RNA(\%variables);
print "ok\n";
