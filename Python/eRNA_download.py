# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 15:09:03 2017

@author: yuan
"""
#standard lib
import sys

#in-house lib
import myDownload
import myGenome
import myIO
import mySystem

################################################################################
class download_annot:
    def __init__(self, par):
        self.par=par

#download genome annotations
    def genome_annot(self):
        if par['web_site']=='ENSEML':
            #1: select data type
            data_types=['dna_fa','genome_CDS','gtf','gff', 'protein']
            self.par['data_type']=mySystem.system().select_key(data_types)
            #2: download
            if par['data_type']=='dna_fa':
                myDownload.ensembl(par['specie'], par['dir_out']).download_dna()
            else:
                myDownload.ensembl(par['specie'], par['dir_out']).download_annot(par['data_type'])
        elif par['web_site']=='NCBI':
            #1: select data type
            data_types=['dna_fa','RNA', 'gff', 'protein']
            self.par['data_type']=mySystem.system().select_key(data_types)
            #2: donwload
            if par['data_type']=='dna_fa':
                myDownload.NCBI(par['specie'], par['dir_out']).download_dna()
            else:
                myDownload.NCBI(par['specie'], par['dir_out']).download_annot(par['data_type'])
            
            
#match fasta with gtf        
    def match_fasta(self):
        files=myIO.dir_os(self.par['dir_out']).incrusive_files()
        #select a fasta file
        fa_files=filter(lambda x: x.endswith(('.fa','.fasta')), files)
        self.par['match_fa']=mySystem.system().select_key(fa_files)
        #select a gtf or gff file
        gtf_files=filter(lambda x: x.endswith(('.gtf','.gff3')), files)
        self.par['match_gtf']=mySystem.system().select_key(gtf_files)
        
        #match
        if par['web_site']=='ENSEML':
            myGenome.genome(par['match_fa']).match_ensembl_fa(par['match_gtf']);
        elif par['web_site']=='NCBI':
            myGenome.genome(par['match_fa']).match_ncbi_fa(par['match_gtf']);

#
    def uniprot_idmapping(self):
        #download the idmapping file
        local_file=myDownload.uniprot(par['dir_out']).download_idmapping()
                    
#################################################################################
if __name__=="__main__":
    #initiate dictionary saving parameters
    par={'in_out':'Continue'};
    annot=download_annot(par)
    ########################################
    #
   
    #2: download dir
    par['dir_home']=myIO.dir_os('/home/yuan/data_preparation/').stdin_dir('Enter the directory path storing downloads files')
    print par['dir_home']
    while(par['in_out']=='Continue'):
        #2:select ftp or web site
        web_sites=['NCBI','ENSEML', 'UniProt']
        par['web_site']=mySystem.system().select_key(web_sites, 'Select public database')
        par['dir_out'] = par['dir_home']+par['web_site']+'/'
        #1: select file types
        if par['web_site'] in ['NCBI','ENSEML']:
            operations=['Genome annotation', 'match fasta and gtf']
            par['operations']=mySystem.system().select_key(operations, 'What is your operations')
        elif par['web_site']=='UniProt':
            par['operations']='UniProt idmapping'        
            
        #download genome annotations
        if par['operations']=='Genome annotation':
            #2: select specie
            species=['human','mouse','rat','maize']
            par['specie']=mySystem.system().select_key(species, 'Select specie of genome')
            #4 launch download
            print par['dir_out']
            download_annot(par).genome_annot()
        #match
        elif par['operations']=='match fasta and gtf':
            download_annot(par).match_fasta()
        #
        elif par['operations']=='UniProt idmapping':
            download_annot(par).uniprot_idmapping()
            
        #decide going on or stop
        #6: other operations
        par['in_out']=mySystem.system().select_key(['Continue','Exit'], 'Next')

    
    print par
    print "\n\nGreat! It is done!\n\n"
#end