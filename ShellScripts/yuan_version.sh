#!/bin/bash

#Tiezheng Yuan 20180914
#check version of installed software at Linux VM (CentOS 7)

#load $PATH
source ~/.bashrc


#software list
arr_sw=(Albacore Artemis augustus bamtools barrnap bcftools bedtools bowtie1 bowtie2 BUSCO canu centrifuge Circos hisat2 hmmer infernal japsa minced miniasm minimap2 miRDeep2 miR-PREFeR.py NanoFilt NanoPlot parallelGNU pilon Prodigal prokka quast racon randfold samtools SPAdes Streamformatics stringtie tabix tbl2asn Unicycler Vienna-RNAfold)
arr_ver=('read_fast5_basecaller.py --version' 'art -h' 'augustus --version' 'bamtools --version' 'barrnap -v' 'bcftools --version' 'bedtools --version' 'bowtie --version' 'bowtie2 --version' 'run_BUSCO.py -v' 'canu --version' 'centrifuge --version' 'circos -v' 'hisat2 --version' 'hmmalign -h' 'cmbuild -h' 'jsa' 'minced --version' 'miniasm -V' 'minimap2 -V' 'miRDeep2.pl' 'mirprefer' 'NanoFilt --version' 'NanoPlot -v' 'parallel -h' 'pilon -v' 'prodigal -v' 'prokka -v' 'quast.py -v' 'racon --version' 'randfold' 'samtools --version' 'spades.py -v' 'npm -v' 'stringtie --version' 'tabix --version' 'tbl2asn' 'unicycler --version' 'RNAfold --version')

#return version
for (( i=0; i<${#arr_sw[@]}; i++));do
 echo '#################################################################'
 echo -e "\t\t###\t$(($i+1)): ${arr_sw[i]}\t###"
 echo '#################################################################'
 #return version of software
 ${arr_ver[i]}
 echo -e '_____________\n\n'
done


#end

