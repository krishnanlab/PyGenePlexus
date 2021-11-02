import argparse
import utls
import numpy

parser = argparse.ArgumentParser()
parser.add_argument('-t','--task',
                    default = 'validate',
                    type = str,
                    help = 'options are validate or results')
parser.add_argument('-i','--input',
                    default = 'input_genes.txt',
                    type = str,
                    help = 'file path to the gene functions')
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

args = parser.parse_args()


if args.task == 'validate':
    input_genes = np.loadtxt(args.input,dtype=str,delimiter=', ')
    input_genes = [item.strip("'") for item in input_genes_txtfile]
    convert_IDs, df_convert_out = utls.intial_ID_convert(input_genes)
    pos_genes_in_net, genes_not_in_net, net_genes = utls.get_genes_in_network(convert_IDs,net_type)
    df_convert_out = utls.make_validation_df(df_convert_out,pos_genes_in_net)