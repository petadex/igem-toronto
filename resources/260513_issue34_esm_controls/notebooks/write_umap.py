import numpy as np
import pandas as pd
import psycopg2
import umap
import matplotlib.pyplot as plt

control_embeds = np.load('../data/family_embedding_controls.esm2-150M.d640.n517840.npz')
real_embeds = np.load('../data/family_embeddings.esm2-150M.d640.n64730.npz')

merged_embeds = np.concatenate([control_embeds["embeddings"], real_embeds["embeddings"]])

reducer = umap.UMAP()
embedding_2d = reducer.fit_transform(merged_embeds)

df = pd.DataFrame(embedding_2d, columns=['x', 'y'])
df.to_csv('../data/umap_embeddings.csv', index=False)