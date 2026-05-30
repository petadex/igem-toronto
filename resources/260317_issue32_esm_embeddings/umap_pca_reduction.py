import numpy as np
import pandas as pd
from umap import UMAP
import plotly.express as px
from sklearn.decomposition import PCA

# 1. Load your ESM-2 embeddings and sequence names
embeddings = np.load('C:/Users/Pixel/OneDrive/Documents/PETase_Stuff/science_seqs/t36_embeddings.npy', allow_pickle=True)
seq_names = np.load('C:/Users/Pixel/OneDrive/Documents/PETase_Stuff/science_seqs/t36_names.npy', allow_pickle=True)

# 2. Load your activity data from TSV
activity_df = pd.read_csv("C:/Users/Pixel/OneDrive/Documents/PETase_Stuff/science_seqs/label_activity_scaled.tsv", sep="\t",
                         header=None)
activity_df.columns = ['sequence_name', 'activity']

# Apply padding to solve complex array shape problem
max_length = max(len(e) for e in embeddings)
print(f"Shape of loaded embeddings: {len(embeddings)} sequences found")
print(f"Maximum sequence length: {max_length}")

embeddings_padded = []
for e in embeddings:
    padded_embedding = np.pad(e, ((0, max_length - len(e)), (0, 0)), mode='constant', constant_values=0)
    embeddings_padded.append(padded_embedding)

embeddings_padded = np.array(embeddings_padded)
print(f"Shape of padded embeddings: {embeddings_padded.shape}")

# Reshape for PCA
embeddings_flat = embeddings_padded.reshape(embeddings_padded.shape[0], -1)
print(f"Shape of flattened embeddings: {embeddings_flat.shape}")

# 3. Apply PCA to reduce to 50 dimensions
print("Applying PCA to reduce dimensions to 50...")
pca = PCA(n_components=50)
pca_result = pca.fit_transform(embeddings_flat)
print(f"Shape after PCA: {pca_result.shape}")

# Print explained variance information
explained_variance = pca.explained_variance_ratio_
print(f"Total variance explained by 50 PCs: {sum(explained_variance)*100:.2f}%")
print(f"Variance explained by first PC: {explained_variance[0]*100:.2f}%")
print(f"Variance explained by last PC: {explained_variance[-1]*100:.2f}%")

# 4. Apply UMAP to the PCA results
print("Applying UMAP on PCA results...")
umap_reducer = UMAP(n_components=3, random_state=42)
umap_embeddings = umap_reducer.fit_transform(pca_result)

# 5. Create result dataframe
result_df = pd.DataFrame({
    'sequence_name': seq_names,
    'UMAP_1': umap_embeddings[:, 0],
    'UMAP_2': umap_embeddings[:, 1],
    'UMAP_3': umap_embeddings[:, 2]
})

# 6. Merge with activity data
result_df = result_df.merge(activity_df, on='sequence_name')

# 7. Save the final output
output_df = result_df[['sequence_name', 'UMAP_1', 'UMAP_2', 'UMAP_3', 'activity']]
output_df.to_csv("protein_pca_umap_3d_with_activity.csv", index=False)

# 8. Create an interactive 3D visualization with Plotly
fig = px.scatter_3d(
    result_df,
    x='UMAP_1',
    y='UMAP_2',
    z='UMAP_3',
    color='activity',
    color_continuous_scale='viridis',
    opacity=0.8,
    hover_name='sequence_name',
    hover_data={'UMAP_1': True, 'UMAP_2': True, 'UMAP_3': True, 'activity': ':.4f'},
    title='Interactive 3D UMAP of PCA-Reduced Embeddings Colored by Activity'
)

# Customize the layout
fig.update_layout(
    scene=dict(
        xaxis_title='UMAP Dimension 1',
        yaxis_title='UMAP Dimension 2',
        zaxis_title='UMAP Dimension 3',
        bgcolor='white'
    ),
    margin=dict(r=20, l=10, b=10, t=30),
    coloraxis_colorbar=dict(title="Activity"),
    width=900,
    height=700
)

# Save as interactive HTML
fig.write_html("interactive_pca_umap_3d.html")

# Display the figure
fig.show()