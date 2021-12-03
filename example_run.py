import geneplexus
import numpy as np
import argparse


'''
For the jobname need to figure out how to save it with unique idntifier
'''

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',
                    default = '/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/data/disease_gene_files/Alzheimers_Disease.txt',
                    type = str,
                    help = 'file path to the gene functions')
parser.add_argument('-fl','--file_loc',
                    default = 'HPCC',
                    type = str,
                    help = 'options local, HPCC, cloud')
parser.add_argument('-n','--net_type',
                    default = 'BioGRID',
                    type = str,
                    help = 'options are BioGRID, STRING-EXP, STRING, GIANT-TN')
parser.add_argument('-f','--features',
                    default = 'Embedding',
                    type = str,
                    help = 'options are Embedding, Adjacency, Influence')
parser.add_argument('-g','--GSC',
                    default = 'GO',
                    type = str,
                    help = 'options are GO, DisGeNet')
parser.add_argument('-j','--jobname',
                    default = 'mynameargparse_AZ',
                    type = str,
                    help = 'sting to use for jobname')
parser.add_argument('-s','--fp_save',
                    default = 'results/',
                    type = str,
                    help = 'The path to save the file to (needs the / in it for now)')
args = parser.parse_args()
#### in the above file_loc not being used right now ####


'''
Example of just reading in the gene by user
'''
# input_genes = np.loadtxt('input_genes.txt',dtype=str,delimiter=', ')
# input_genes = [item.strip("'") for item in input_genes]
# # geneplexus.validate_genes(input_genes)
# geneplexus.run_model(input_genes,'BioGRID','GO','Embedding','mytestjob')



'''
Example allowing user to use the command line
'''
input_genes = geneplexus.read_input_file(args.input) # maybe add delimter as argparse argument
### The read_input_file needs some work and not same as Doug's ###
geneplexus.run_model(input_genes,args.net_type,args.GSC,args.features,args.jobname,args.fp_save)
