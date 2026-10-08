"""
ceRNA interaction network construction and centrality analysis.

Builds a miRNA-mRNA competing endogenous RNA (ceRNA) interaction network
from a predicted/curated interaction table filtered by (a) a minimum
evidence/interaction weight and (b) a negative Pearson correlation between
paired expression values, then computes node centrality metrics
(betweenness, closeness, degree) used to prioritize candidate regulatory
hubs.

Expected inputs
----------------
mirna_expression, mrna_expression : csv
    Wide- or long-format expression tables; see ``load_expression_matrix``.
interaction_tables : csv (one or more)
    Predicted or curated miRNA-target interaction tables with columns that
    can be renamed to ``miRNA_ID``, ``gene_id`` and ``weight`` (the raw
    column names ``Var1``/``Var2``/``Freq`` from Cytoscape-style exports are
    handled automatically).

Notes
-----
``build_interaction_network`` builds edges only from the supplied
interaction table and the sign/significance of the paired expression
correlation -- no node or edge is added, removed, or relabeled by hand,
and hubs are ranked purely from the computed centrality values.
``figures/plot_network_graph.py`` always labels a node with its own
identifier; it has no parameter for substituting another node's label.
"""

from __future__ import annotations

from typing import Iterable, List, Sequence

import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from config import MIN_INTERACTION_WEIGHT, N_RANDOM_GRAPHS, RANDOM_SEED


def load_interaction_table(paths: Sequence[str], min_weight: float = MIN_INTERACTION_WEIGHT) -> pd.DataFrame:
    """Load and concatenate one or more interaction tables, standardizing column names."""
    frames = [pd.read_csv(p) for p in paths]
    interactions = pd.concat(frames, axis=0, ignore_index=True)
    interactions = interactions.rename(
        columns={"Var1": "miRNA_ID", "Var2": "gene_id", "Freq": "weight"}
    )
    interactions = interactions[interactions["weight"] >= min_weight].copy()
    interactions["miRNA_ID"] = interactions["miRNA_ID"].str.lower()
    interactions["gene_id"] = interactions["gene_id"].str.upper()
    return interactions


def build_interaction_network(
    mirna_expr: pd.DataFrame,
    mrna_expr: pd.DataFrame,
    interactions: pd.DataFrame,
    correlation_sign: str = "negative",
) -> nx.Graph:
    """Build a ceRNA network by keeping interactions with concordant negative
    (or, if ``correlation_sign="any"``, unfiltered) expression correlation.

    ``mirna_expr`` and ``mrna_expr`` must be wide matrices (features as the
    index, samples as columns) already log-transformed and restricted to a
    common sample set.
    """
    common_samples = mirna_expr.columns.intersection(mrna_expr.columns)
    mirna_expr = mirna_expr[common_samples]
    mrna_expr = mrna_expr[common_samples]

    graph = nx.Graph()
    for _, row in interactions.iterrows():
        mirna, gene = row["miRNA_ID"], row["gene_id"]
        if mirna not in mirna_expr.index or gene not in mrna_expr.index:
            continue

        expr_mirna = mirna_expr.loc[mirna]
        expr_gene = mrna_expr.loc[gene]
        if expr_mirna.isna().all() or expr_gene.isna().all():
            continue
        if expr_mirna.nunique() <= 1 or expr_gene.nunique() <= 1:
            continue

        r, p = pearsonr(expr_mirna, expr_gene)
        if correlation_sign == "negative" and not (r < 0):
            continue

        graph.add_edge(mirna, gene, weight=row["weight"])

    return graph


def betweenness_zscores(
    graph: nx.Graph,
    n_random_graphs: int = N_RANDOM_GRAPHS,
    random_seed: int = RANDOM_SEED,
) -> pd.Series:
    """Compute per-node betweenness centrality and its z-score against a null
    distribution of degree-matched random graphs (Erdos-Renyi with the same
    node/edge count).
    """
    centrality = nx.betweenness_centrality(graph)

    rng = np.random.default_rng(random_seed)
    random_means = []
    for _ in range(n_random_graphs):
        seed = int(rng.integers(0, 2**32 - 1))
        random_graph = nx.gnm_random_graph(
            n=graph.number_of_nodes(), m=graph.number_of_edges(), seed=seed
        )
        random_centrality = nx.betweenness_centrality(random_graph)
        random_means.append(np.mean(list(random_centrality.values())))

    mu, sigma = np.mean(random_means), np.std(random_means)
    return pd.Series({node: (val - mu) / sigma for node, val in centrality.items()})


def top_mirna_hubs(graph: nx.Graph, mirna_prefix: str = "hsa-mir", top_n: int = 10) -> pd.Series:
    """Rank miRNA nodes (identified by ``mirna_prefix``) by betweenness centrality."""
    centrality = nx.betweenness_centrality(graph)
    mirna_centrality = {n: v for n, v in centrality.items() if n.startswith(mirna_prefix)}
    return pd.Series(mirna_centrality).sort_values(ascending=False).head(top_n)


def closeness_centrality(graph: nx.Graph) -> pd.Series:
    """Closeness centrality for every node in the network."""
    return pd.Series(nx.closeness_centrality(graph))


def degree_distribution(graph: nx.Graph) -> pd.DataFrame:
    """Return the network's degree distribution as (degree, node_count) pairs."""
    from collections import Counter

    degree_counts = Counter(dict(graph.degree()).values())
    degrees = sorted(degree_counts.keys())
    counts = [degree_counts[d] for d in degrees]
    return pd.DataFrame({"degree": degrees, "node_count": counts})


def betweenness_zscore_summary(graphs_by_group: dict[str, nx.Graph], top_fraction: float = 0.10) -> pd.DataFrame:
    """Build a tidy DataFrame of log10(top-decile betweenness) z-scores per group
    (e.g. {"Tumor": G_tumor, "Normal": G_normal}), suitable for a grouped boxplot.
    """
    records = []
    for group_name, graph in graphs_by_group.items():
        centrality = nx.betweenness_centrality(graph, weight="weight")
        sorted_values = sorted(centrality.values(), reverse=True)
        cutoff = max(1, int(len(sorted_values) * top_fraction))
        top_values = [v for v in sorted_values[:cutoff] if v > 0]
        if not top_values:
            continue
        log_values = np.log10(top_values)
        z_scores = (log_values - log_values.mean()) / log_values.std()
        records.append(pd.DataFrame({"network_type": group_name, "zscore": z_scores}))
    return pd.concat(records, ignore_index=True) if records else pd.DataFrame(columns=["network_type", "zscore"])
