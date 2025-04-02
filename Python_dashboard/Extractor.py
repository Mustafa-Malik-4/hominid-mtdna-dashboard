
import pandas as pd
from Bio import Entrez, SeqIO

### Provide email address to NCBI
Entrez.email = "mrmalik@ucdavis.edu"

def fetch_sequence(accession):
    ### Fetch sequence from NCBI with accession number
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return record

### Create two lists for two different ways to identify hominids in the code that follows.

names = ['German Neanderthal','Russian Neanderthal', 'European Human',
         'Mountain Gorilla Rwanda', 'Chimp Troglodytes', 'Puti Orangutan',
         'Jari Orangutan', 'Western Lowland Gorilla', 'Eastern Lowland Gorilla',
         'Chimp Schweinfurthii', 'Chimp Vellerosus', 'Chimp Verus', 'Altai Denisovan']

numbers = ['AF011222', 'AF254446', 'X90314', 'AF089820', 'AF176766',
           'AF451972', 'AF451964', 'AY079510', 'AF050738', 'AF176722',
           'AF315498', 'AF176731', 'FN673705']

names2 = ['European Neanderthal', 'European Human', 'Altai Denisovan', 'Sumatran Orangutan',
         'Bornean Orangutan', 'Western Lowland Gorilla', 'Eastern Lowland Gorilla', 'Bonobo',
         'Chimp Schweinfurthii', 'Chimp Ellioti', 'Chimp Verus', 'Chimp Troglodytes']


numbers2 = ['NC_011137','NC_012920','FN673705','NC_002083',
            'NC_001646','NC_011120','KM242275','NC_001644',
            'JF727201','KM679417','JF727217','JF727180']

### Now create a dataframe using the accession numbers to retrieve the sequences.

hominids = dict()
index = 0
for accession_number in numbers2:
    sequence_record = fetch_sequence(accession_number)
    sequence = str(sequence_record.seq)
    hominids[names2[index]] = sequence.strip()
    index += 1

df = pd.DataFrame(list(hominids.items()), columns=['Species', 'mtDNA'])

lengths = []
for seq in df['mtDNA']:
     length = len(seq)
     lengths.append(length)
df['Length'] = lengths

df.to_csv('/Users/mustafam/Desktop/dashboard_test/python_dash/mtDNA.csv', index=False)