"""
Visualize a ceRNA miRNA-mRNA interaction network, with miRNA and mRNA nodes
colored separately and (optionally) a subset of nodes labeled.

Design note on node labeling
----------------------------
``plot_network_graph`` always labels a highlighted node with its own
identifier, taken directly from the supplied graph; there is no mechanism
for substituting a different node's label or value, and no node or edge is
added by the plotting code. The network passed in should come only from
``analysis.network_centrality.build_interaction_network`` (or an
equivalent function that builds edges purely from the interaction table
and expression correlation), with no manually added nodes/edges.

Expected input
--------------
A ``networkx.Graph`` such as one produced by
``analysis.network_centrality.build_interaction_network``, whose miRNA
nodes are identifiable by a shared prefix (default: ``"hsa-mir"``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import networkx as nx

from config import FIGURE_DPI
from utils import apply_default_style

MIRNA_COLOR = "#E78AC3"
MRNA_COLOR = "#87CEEB"


def plot_network_graph(
    graph: nx.Graph,
    output_path: Path,
    mirna_prefix: str = "hsa-mir",
    highlight_nodes: Sequence[str] = (),
    layout_seed: int = 42,
) -> None:
    """Draw the network with miRNA/mRNA nodes colored separately.

    Every node in ``highlight_nodes`` is labeled with its own identifier --
    there is no way to display a different label at a node's position, by
    design (see the module docstring).
    """
    apply_default_style()

    mirna_nodes = [n for n in graph.nodes() if n.startswith(mirna_prefix)]
    mrna_nodes = [n for n in graph.nodes() if not n.startswith(mirna_prefix)]

    fig = plt.figure(figsize=(16, 12))
    pos = nx.kamada_kawai_layout(graph)

    nx.draw_networkx_nodes(
        graph, pos, nodelist=mirna_nodes, node_color=MIRNA_COLOR, alpha=0.95,
        node_size=600, edgecolors="black", linewidths=0.8,
    )
    nx.draw_networkx_nodes(
        graph, pos, nodelist=mrna_nodes, node_color=MRNA_COLOR, alpha=0.95,
        node_size=200, edgecolors="black", linewidths=0.8,
    )
    nx.draw_networkx_edges(graph, pos, alpha=0.8)

    for node in highlight_nodes:
        if node not in pos:
            continue
        x, y = pos[node]
        plt.text(
            x, y, node, fontsize=22, ha="center", va="center",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, boxstyle="round,pad=0.2"),
        )

    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor=MIRNA_COLOR, edgecolor=MIRNA_COLOR, label="miRNA"),
        Patch(facecolor=MRNA_COLOR, edgecolor=MRNA_COLOR, label="mRNA"),
    ]
    plt.legend(handles=legend_elements, loc="upper center", bbox_to_anchor=(0.7, 0.95), fontsize=20, frameon=True)

    plt.axis("off")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="png", dpi=FIGURE_DPI, transparent=True)
    plt.close(fig)
