import argparse
import geneplexus2
import numpy as np
import pandas as pd

# Main thing to discuss is how to handle the file_loc problem
# Also, when to save objects? Probbaly not have that be something we do for the user?

###################################################################################################################

# The first step is for a user to load up a set of genes as a list
# Need to figure out best way for doing multiple sets at once

input_genes = np.loadtxt('input_genes.txt',dtype=str,delimiter=', ')
input_genes = [item.strip("'") for item in input_genes]

###################################################################################################################

myclass = geneplexus2.GenePlexus()
myclass.load_genes(input_genes)
convert_IDs, df_convert_out = myclass.intial_ID_convert()
df_convert_out, table_summary, input_count = myclass.make_validation_df()
myclass.set_params('BioGRID','Embedding','GO')
pos_genes_in_net, genes_not_in_net, net_genes = myclass.get_genes_in_network()
print(pos_genes_in_net)


    

    
