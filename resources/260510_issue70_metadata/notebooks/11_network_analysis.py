# ---
# jupyter:
#   jupytext:
#     formats: py:light
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Notebook 11 — Network Analysis (NetworkX)
#
# Models the plastic biodegradation landscape as a **bipartite network**:
# - Nodes: organisms (blue) + plastic types (orange)
# - Edges: reported degradation activity
#
# Computes graph metrics (degree centrality, betweenness, clustering) to identify
# *hub organisms* (broad degraders) and *bridge organisms* (connecting otherwise
# separate plastic-degradation communities). Also builds an organism–organism
# similarity network based on shared plastic substrates.

# +
import sys
from pathlib import Path as _P
sys.path.insert(0, str(_P(__file__).parent.parent))
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from collections import defaultdict

from src.data_loader import load_all
from src.analysis import plastic_co_occurrence

data = load_all()
df = data['plasticdb']
organisms = data['organisms']
# -

# ## 1. Build the bipartite organism–plastic network

B = nx.Graph()
top_plastics = df['plastic'].value_counts().head(20).index.tolist()
df_top = df[df['plastic'].isin(top_plastics)].dropna(subset=['organism', 'plastic'])

organism_nodes = set()
plastic_nodes = set()
for _, row in df_top.iterrows():
    org = row['organism']
    pla = row['plastic']
    B.add_node(org, bipartite=0, node_type='organism')
    B.add_node(pla, bipartite=1, node_type='plastic')
    if not B.has_edge(org, pla):
        B.add_edge(org, pla, weight=1)
    else:
        B[org][pla]['weight'] += 1
    organism_nodes.add(org)
    plastic_nodes.add(pla)

print(f"Bipartite network:")
print(f"  Organism nodes: {len(organism_nodes):,}")
print(f"  Plastic nodes:  {len(plastic_nodes)}")
print(f"  Edges:          {B.number_of_edges():,}")
print(f"  Is bipartite:   {nx.is_bipartite(B)}")

# ## 2. Degree centrality — hub organisms

degree_dict = dict(B.degree())
organism_degrees = {n: d for n, d in degree_dict.items() if n in organism_nodes}
top_organisms = sorted(organism_degrees, key=organism_degrees.get, reverse=True)[:20]

print("\nTop 20 hub organisms (highest degree = most plastic types):")
for org in top_organisms:
    plastics_for_org = [n for n in B.neighbors(org)]
    print(f"  {org:<45s}: degree={organism_degrees[org]:2d}  "
          f"plastics=[{', '.join(sorted(plastics_for_org)[:4])}{'...' if len(plastics_for_org)>4 else ''}]")

# ## 3. Plastic node degrees

plastic_degrees = {n: d for n, d in degree_dict.items() if n in plastic_nodes}
print("\nPlastic types by number of unique degrading organisms (top 15):")
for pla in sorted(plastic_degrees, key=plastic_degrees.get, reverse=True)[:15]:
    print(f"  {pla:<10s}: {plastic_degrees[pla]:4d} organisms")

# ## 4. Visualise bipartite network (top 80 organisms by degree)

top80_orgs = sorted(organism_degrees, key=organism_degrees.get, reverse=True)[:80]
subgraph_nodes = set(top80_orgs) | plastic_nodes
B_sub = B.subgraph(subgraph_nodes)

pos = {}
plastic_list = sorted(plastic_nodes)
org_subset = sorted(top80_orgs)
for i, pla in enumerate(plastic_list):
    pos[pla] = (0, i * 2)
for i, org in enumerate(org_subset):
    pos[org] = (5, i * (len(plastic_list) * 2 / len(org_subset)))

fig, ax = plt.subplots(figsize=(14, 18))
org_color   = '#2E86AB'
pla_color   = '#F18F01'
node_colors = [org_color if B_sub.nodes[n]['node_type'] == 'organism' else pla_color
               for n in B_sub.nodes()]
node_sizes  = [20 + 5 * B_sub.degree(n) if B_sub.nodes[n]['node_type'] == 'organism' else 400
               for n in B_sub.nodes()]
edge_weights = [B_sub[u][v].get('weight', 1) for u, v in B_sub.edges()]
nx.draw_networkx(
    B_sub, pos=pos, ax=ax,
    node_color=node_colors, node_size=node_sizes,
    edge_color='gray', alpha=0.6, width=0.5,
    font_size=6, with_labels=True,
    labels={n: n if n in plastic_nodes else ' '.join(n.split()[:2])
            for n in B_sub.nodes()},
)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=org_color, label='Organism'),
                   Patch(color=pla_color, label='Plastic Type')])
ax.set_title('Bipartite Organism–Plastic Network\n(top 80 organisms by degree)', fontsize=13)
ax.axis('off')
plt.tight_layout()
plt.savefig('outputs/figures/11_bipartite_network.png', dpi=120, bbox_inches='tight')
plt.show()

# ## 5. Projected organism–organism similarity network

from networkx.algorithms import bipartite
organism_graph = bipartite.weighted_projected_graph(B, organism_nodes)
print(f"\nOrganism–organism projection:")
print(f"  Nodes: {organism_graph.number_of_nodes():,}")
print(f"  Edges: {organism_graph.number_of_edges():,}")

# Keep only top-degree organisms for clarity
top100 = sorted(organism_degrees, key=organism_degrees.get, reverse=True)[:100]
G_sub = organism_graph.subgraph(top100)

# Betweenness centrality — bridge species
betweenness = nx.betweenness_centrality(G_sub, weight='weight', normalized=True)
top_bridges = sorted(betweenness, key=betweenness.get, reverse=True)[:15]
print("\nTop 15 bridge organisms (betweenness centrality):")
for org in top_bridges:
    print(f"  {org:<45s}: {betweenness[org]:.5f}")

# ## 6. Network community detection (greedy modularity)

try:
    from networkx.algorithms.community import greedy_modularity_communities
    communities = list(greedy_modularity_communities(G_sub, weight='weight'))
    print(f"\nCommunity detection: {len(communities)} communities found")
    for i, comm in enumerate(sorted(communities, key=len, reverse=True)[:5]):
        print(f"  Community {i+1}: {len(comm)} organisms  "
              f"| {', '.join(list(comm)[:4])}{'...' if len(comm)>4 else ''}")
    modularity = nx.algorithms.community.quality.modularity(
        G_sub, communities, weight='weight'
    )
    print(f"  Modularity Q = {modularity:.4f}")
except Exception as e:
    communities = []
    print(f"Community detection note: {e}")

# ## 7. Plastic co-occurrence network

co = plastic_co_occurrence(organisms)
top12_pla = co.sum().nlargest(12).index.tolist()
co_sub_arr = co.loc[top12_pla, top12_pla].values.astype(float)
np.fill_diagonal(co_sub_arr, 0)
co_sub = pd.DataFrame(co_sub_arr, index=top12_pla, columns=top12_pla)

G_pla = nx.from_pandas_adjacency(co_sub)
pos_pla = nx.spring_layout(G_pla, seed=42, weight='weight')
edge_widths = [G_pla[u][v]['weight'] / 30 for u, v in G_pla.edges()]
node_sizes  = [500 + 20 * G_pla.degree(n, weight='weight') for n in G_pla.nodes()]

fig, ax = plt.subplots(figsize=(10, 8))
nx.draw_networkx(
    G_pla, pos=pos_pla, ax=ax,
    node_color='#44BBA4', node_size=node_sizes,
    edge_color='gray', width=edge_widths, alpha=0.7,
    font_size=9, font_weight='bold',
)
ax.set_title('Plastic Co-occurrence Network\n(edge weight = shared degrading organisms)', fontsize=13)
ax.axis('off')
plt.tight_layout()
plt.savefig('outputs/figures/11_plastic_cooccurrence_network.png', dpi=150, bbox_inches='tight')
plt.show()

# ## 8. Graph summary statistics

print("\n=== Network Summary ===")
print(f"Bipartite (org+plastic):  nodes={B.number_of_nodes():,}  edges={B.number_of_edges():,}  "
      f"density={nx.density(B):.5f}")
print(f"Organism projection:      nodes={organism_graph.number_of_nodes():,}  "
      f"edges={organism_graph.number_of_edges():,}  density={nx.density(organism_graph):.5f}")
print(f"Plastic co-occurrence:    nodes={G_pla.number_of_nodes()}  "
      f"edges={G_pla.number_of_edges()}  density={nx.density(G_pla):.4f}")
avg_clustering = nx.average_clustering(G_sub, weight='weight')
print(f"Avg clustering (organism projection top-100): {avg_clustering:.4f}")
