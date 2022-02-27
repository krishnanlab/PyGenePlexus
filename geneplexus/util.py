def read_gene_list(path):
    return [gene.strip("'") for gene in open(path, "r").read().split(", ")]
