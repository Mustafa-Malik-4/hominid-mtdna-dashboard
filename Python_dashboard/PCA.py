
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from collections import Counter
import matplotlib.pyplot as plt

df = pd.read_csv('mtDNA.csv')

# Function to compute k-mer frequencies
def kmer_freq(seq, k=3):
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
    counts = Counter(kmers)
    total_kmers = sum(counts.values())
    return {kmer: counts[kmer] / total_kmers for kmer in counts}

# Set k-mer length 
k = 3

# Compute k-mer frequencies for each genome
df['Kmer_Freq'] = df['mtDNA'].apply(lambda seq: kmer_freq(seq, k))

# Collect all k-mers for consistent feature space
all_kmers = set()
for kmer_counts in df['Kmer_Freq']:
    all_kmers.update(kmer_counts.keys())

# Convert into a feature matrix
kmer_features = pd.DataFrame([{kmer: kmer_dict.get(kmer, 0) for kmer in all_kmers} for kmer_dict in df['Kmer_Freq']])

# Standardize
scaler = StandardScaler()
kmer_features_scaled = scaler.fit_transform(kmer_features)

# PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(kmer_features_scaled)

# PCA dataframe
pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
pca_df['Species'] = df['Species']

# Plotting
plt.scatter(pca_df['PC1'], pca_df['PC2'], marker='o')
for i, species in enumerate(pca_df['Species']):
    plt.text(pca_df.loc[i, 'PC1'], pca_df.loc[i, 'PC2'], species, fontsize=8)
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("PCA for NIH Hominid mtDNA Sequences")
plt.show()