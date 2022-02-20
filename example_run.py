import time

import numpy as np

from geneplexus import geneplexus


t = time.perf_counter()
input_genes = np.loadtxt("input_genes.txt", dtype=str, delimiter=", ")
input_genes = [item.strip("'") for item in input_genes]
geneplexus.download_select_data("/Users/christophermancuso/Documents/DataSets/from_Azure2/", tasks="NetworkGraph")
print(f"Took {time.perf_counter() - t:.3f} sec to run.")
