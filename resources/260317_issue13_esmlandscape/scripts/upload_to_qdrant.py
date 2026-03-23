from qdrant_client import QdrantClient
from qdrant_client.models import VectorParams, Distance
import numpy as np
import pandas as pd
import uuid
import math

def upload_embeddings_to_qdrant():
    try:
        client = QdrantClient(host="localhost", port=6333)
    except:
        print("Failed to connect to Qdrant. Make sure it's running on localhost:6333.")
        return
    # Load the embeddings and family IDs
    data = np.load('./data/embeddings.npz')
    family_ids = data['family_ids']
    embeddings = data['embeddings']
    representatives_df = pd.read_csv('./data/representatives.csv')

    # Prepare the data for Qdrant
    points = []
    for i in range(len(family_ids)):
        family_id = family_ids[i]
        embedding = embeddings[i].tolist()  # Convert to list for JSON serialization
        representative_info = representatives_df[representatives_df['family_id'] == family_id]
        if not representative_info.empty:
            representative_seq = str(representative_info['sequence'].iloc[0])
        else:
            representative_seq = "Unknown"
        
        point = {
            "id": str(uuid.uuid4()),  # Generate a unique ID for each point
            "vector": embedding,
            "payload": {
                "family_id": int(family_id),
                "representative_seq": representative_seq
            }
        }
        points.append(point)

    # Upload the points to Qdrant
    try:
        if not client.collection_exists(collection_name="esm_embeddings"):
            client.create_collection(
                collection_name="esm_embeddings",
                vectors_config=VectorParams(
                    size=embeddings.shape[1],
                    distance=Distance.COSINE
                )
            )
        total_batches = 1000
        batch_size = math.ceil(len(points) / total_batches) if points else 0

        uploaded_points = 0
        for batch_num in range(total_batches):
            start = batch_num * batch_size
            end = min(start + batch_size, len(points))
            batch_points = points[start:end]

            if not batch_points:
                continue

            client.upsert(
                collection_name="esm_embeddings",
                points=batch_points
            )
            uploaded_points += len(batch_points)
            print(f"Uploaded batch {batch_num + 1}/{total_batches} ({len(batch_points)} points).")

        print(f"Successfully uploaded {uploaded_points} points to Qdrant in {total_batches} batches.")
    except Exception as e:
        print(f"Failed to upload points to Qdrant: {e}")

if __name__ == "__main__":
    upload_embeddings_to_qdrant()