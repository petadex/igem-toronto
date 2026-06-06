import pandas as pd
import numpy as np

SEED = 67
np.random.seed(SEED)

def build_fragments(centroids, length, replicas):
    seq_lengths = centroids['sequence'].str.len().to_numpy()
    frag_lengths = (seq_lengths * length).astype(int)
    ids = centroids['30pid_superfamily_id'].to_numpy()
    sequences = centroids['sequence'].to_numpy()

    rows = []
    for seq, frag_length, sid, seq_len in zip(sequences, frag_lengths, ids, seq_lengths):
        starts = np.random.randint(0, seq_len - frag_length, size=replicas)
        for j, start in enumerate(starts):
            rows.append({
                'id': sid,
                'slug': f'{sid}_replicate_{j}_start_{start}',
                'sequence': seq[start:start + frag_length],
            })

    fragments_df = pd.DataFrame(rows)
    output_path = f'../data/fragments_length_{length}_replicas_{replicas}.csv'
    fragments_df.to_csv(output_path, index=False)
    print(f'Fragments saved to {output_path}')
    return True

if __name__ == "__main__":
    centroids = pd.read_csv('../superfamily_clusters_with_sequences.csv')
    for i in range(9):
        length = 0.1 + i * 0.1
        build_fragments(centroids, round(length, 1), 5)