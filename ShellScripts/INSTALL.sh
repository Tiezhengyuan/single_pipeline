#!/usr/bin/env bash

#
install_file=$(readlink -f "$0")
#echo $install_file
home_dir=$(dirname "$install_file")
echo "$home_dir"



#################
echo '1. Set soft links of python classes for importing python modules'
python_dir='/anaconda/envs/py35/lib/python3.5/site-packages/altriaseq'
printf '\tDefault python library directory: %s\n' $python_dir
ln -sfn $home_dir $python_dir
printf '...\n...\n...\n'


###########
echo '2. Set shortcuts for launching Python scripts:'
bin_dir='/usr/local/bin/'
printf '\tDefault executable directory: %s\n' $bin_dir
#
bin_names=(tools mRNASeq GenomeAssembly GenomeComparison NanoporeSeq miRNASeq Quant)
file_names=(main.py RNAseq.py GENseq.py COMseq.py NANseq.py MIRseq.py Quant.py)
for i in "${!bin_names[@]}"
do
  main_file=$home_dir'/src/'${file_names[$i]}
  main_ln=$bin_dir'Altria_'${bin_names[$i]}
  printf '\t...%s: %s\n' ${file_names[$i]} $main_ln
  #echo $main_ln
  ln -sfn $main_file $main_ln
  chmod +x $main_ln
done
printf '\n'

#
bin_names=(ballgown)
file_names=(mrna_ballgown.R)
for i in "${!bin_names[@]}"
do
  main_file=$home_dir'/R/'${file_names[$i]}
  main_ln=$bin_dir'Altria_'${bin_names[$i]}
  printf '\t...%s: %s\n' ${file_names[$i]} $main_ln
  #echo $main_ln
  ln -sfn $main_file $main_ln
  chmod +x $main_ln
done
printf '\n'

################
echo 'Installation is done. Great!!!'
#end
