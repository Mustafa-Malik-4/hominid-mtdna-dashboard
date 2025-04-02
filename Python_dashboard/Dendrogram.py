
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from Bio import Entrez, SeqIO
import edlib

### Provide email address to NCBI
Entrez.email = "mrmalik@ucdavis.edu"

def fetch_sequences(accessions):
    ### Fetch sequence from NCBI with accession number
    handle = Entrez.efetch(db="nucleotide", id=",".join(accessions), rettype="fasta")
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    return {accessions[i]: str(records[i].seq) for i in range(len(records))}

names = ['European Neanderthal', 'European Human', 'Altai Denisovan', 'Sumatran Orangutan',
         'Bornean Orangutan', 'Western Lowland Gorilla', 'Eastern Lowland Gorilla', 'Bonobo',
         'Chimp Schweinfurthii', 'Chimp Ellioti', 'Chimp Verus', 'Chimp Troglodytes']

numbers = ['NC_011137', 'NC_012920', 'FN673705', 'NC_002083',
           'NC_001646', 'NC_011120', 'KM242275', 'NC_001644',
           'JF727201', 'KM679417', 'JF727217', 'JF727180']

sequences = fetch_sequences(numbers)
hominids = {names[i]: sequences[numbers[i]] for i in range(len(names))}

### Using the dictionary, we iteratively run edit distance on every combination of species.
### First set up matrix for edit distances
ds_mat = np.zeros((len(names), len(names)))

for i in range(len(names)):
    for j in range(i + 1, len(names)):
        result = edlib.align(hominids[names[i]], hominids[names[j]], task='distance')
        edit_ds = result['editDistance']
        ds_mat[i, j] = ds_mat[j, i] = edit_ds

print(ds_mat)

np.savetxt('ds_mat.csv', ds_mat, delimiter=',')

### Doing the hierarchical clustering again using functions from the 'sch' library:

condensed_dist_mat = sch.distance.squareform(ds_mat)
linkage_matrix = sch.linkage(condensed_dist_mat, method='average')

### Dendrogram
plt.figure(figsize=(10, 7))
sch.dendrogram(linkage_matrix, labels=names)
plt.title('Hominid Species Mitochondrial DNA Edit Distance')
plt.xlabel('Species')
plt.xticks(rotation=90)
plt.ylabel('Distance')
plt.show()
