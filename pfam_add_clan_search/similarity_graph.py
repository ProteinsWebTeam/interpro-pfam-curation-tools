"""
Build a Pfam-Pfam similarity network from Foldseek results and:
1) Suggest additions of Pfams to existing clans (based on neighbour agreement).
2) Discover *de novo* groups among Pfams with no clan (potential new clans).

INPUT format (foldseekpfam_map_all.tsv):
    query, target, target_clan, fident, alnlen, evalue, bitscore

Outputs (written to outdir):
    nodes.csv                 (Cosmograph nodes, includes 'clan' and 'proposed_group')
    edges.csv                 (Cosmograph edges, includes source/target clan and same_clan)
    clan_candidates.tsv       (Pfams w/ NoClan, suggested existing clan)
    proposed_groups.tsv       (Groups among NoClan Pfams, by connected components + modularity)
    summary.txt               (Quick stats)
"""

import argparse
from collections import Counter, defaultdict
from pathlib import Path
import os
import numpy as np
import pandas as pd
import networkx as nx

def infer_weight(df: pd.DataFrame, mode: str) -> pd.Series:
    if mode == 'neglog10_e':
        return -np.log10(df['evalue'].clip(lower=np.finfo(float).tiny))
    elif mode == 'bitscore':
        return df['bitscore']
    elif mode == 'bitscore_per_res':
        return df['bitscore'] / df['alnlen'].replace(0, np.nan)
    else:
        raise ValueError(f"Unknown weight mode: {mode}")


def build_graph(df: pd.DataFrame, weight_col: str) -> nx.Graph:
    G = nx.Graph()
    for _, r in df.iterrows():
        u = r['query']
        v = r['target']
        w = float(r[weight_col])
        u_clan = 'NoClan'
        v_clan = str(r.get('target_clan', '') or 'NoClan')

        if u not in G:
            G.add_node(u, clan=u_clan)
        if v not in G:
            G.add_node(v, clan=v_clan)

        if G.has_edge(u, v):
            if w > G[u][v]['weight']:
                G[u][v].update(
                    weight=w,
                    evalue=float(r['evalue']),
                    bitscore=float(r['bitscore']),
                    alnlen=int(r['alnlen'])
                )
        else:
            G.add_edge(u, v,
                       weight=w,
                       evalue=float(r['evalue']),
                       bitscore=float(r['bitscore']),
                       alnlen=int(r['alnlen']))
    return G


def backfill_clans_from_targets(df: pd.DataFrame, G: nx.Graph) -> None:
    tmap = (df[['target', 'target_clan']]
            .dropna()
            .query("target_clan != ''")
            .drop_duplicates(subset=['target'], keep='first')
            .set_index('target')['target_clan']
            .to_dict())
    for n in G.nodes:
        if G.nodes[n].get('clan', 'NoClan') == 'NoClan':
            if n in tmap:
                G.nodes[n]['clan'] = tmap[n]


def propose_existing_clan_candidates(G: nx.Graph, min_neigh: int, agreement: float) -> pd.DataFrame:
    records = []
    for node in G.nodes:
        if G.nodes[node].get('clan', 'NoClan') != 'NoClan':
            continue
        neigh = list(G.neighbors(node))
        ann = [G.nodes[n]['clan'] for n in neigh if G.nodes[n]['clan'] != 'NoClan']
        if len(ann) < min_neigh:
            continue
        counts = Counter(ann)
        top_clan, top_count = counts.most_common(1)[0]
        frac = top_count / len(ann)
        if frac >= agreement:
            clan_neighs = [n for n in neigh if G.nodes[n]['clan'] == top_clan]
            best_neigh = max(clan_neighs, key=lambda n: G[node][n]['weight'])
            best_weight = G[node][best_neigh]['weight']
            support_weight = sum(G[node][n]['weight'] for n in clan_neighs)
            records.append({
                'pfam': node,
                'suggested_clan': top_clan,
                'supporting_neighbours': top_count,
                'total_annotated_neighbours': len(ann),
                'agreement_fraction': round(frac, 3),
                'support_weight_sum': round(support_weight, 3),
                'example_target': best_neigh,
                'example_weight': round(best_weight, 3),
            })
    return pd.DataFrame(records)


def group_noclan_pfams(G: nx.Graph, min_group_size: int = 3) -> pd.DataFrame:
    nodes_nc = [n for n in G.nodes if G.nodes[n].get('clan', 'NoClan') == 'NoClan']
    H = G.subgraph(nodes_nc).copy()
    groups = []
    gid = 0
    for comp in nx.connected_components(H):
        comp = list(comp)
        if len(comp) < min_group_size:
            continue
        gid += 1
        W = sum(H[u][v]['weight'] for u, v in H.subgraph(comp).edges)
        groups.append({
            'proposed_group_id': f'NCOMP_{gid:03d}',
            'method': 'connected_component',
            'size': len(comp),
            'total_weight': round(float(W), 3),
            'nodes': ';'.join(sorted(comp)),
        })
    return pd.DataFrame(groups)


def annotate_nodes_with_groups(nodes_df: pd.DataFrame, groups_df: pd.DataFrame) -> pd.DataFrame:
    node2groups = defaultdict(list)
    for _, r in groups_df.iterrows():
        gid = r['proposed_group_id']
        for n in str(r['nodes']).split(';'):
            node2groups[n].append(gid)
    nodes_df['proposed_group'] = [','.join(node2groups.get(n, [])) for n in nodes_df['id']]
    return nodes_df


def analyse_foldseek_results(input, sep='\t', weight='bitscore_per_res', min_neigh=1, agreement=0.7, min_group_size=2):

    outdir = os.path.dirname(os.path.abspath(input))
    nodes_file = os.path.join(outdir, "pfam_nodes_all.csv")
    edges_file = os.path.join(outdir, "pfam_edges_all.csv")
    candidates_file = os.path.join(outdir, "clan_candidates_all.csv")
    proposed_groups_file = os.path.join(outdir, "proposed_groups_all.csv")
    summary_file = os.path.join(outdir, "summary.txt")

    cols = ['query', 'target', 'target_clan', 'fident', 'alnlen', 'evalue', 'bitscore']
    df = pd.read_csv(input, sep=sep, names=cols, header=0)

    for c in ['fident', 'alnlen', 'evalue', 'bitscore']:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    df = df.dropna(subset=['alnlen', 'evalue', 'bitscore'])
    df = df[df['alnlen'] > 0]

    df['weight_calc'] = infer_weight(df, weight)
    G = build_graph(df, 'weight_calc')
    backfill_clans_from_targets(df, G)

    cand_df = propose_existing_clan_candidates(G, min_neigh, agreement)
    groups_df = group_noclan_pfams(G, min_group_size)

    nodes_df = pd.DataFrame({
        'id': list(G.nodes),
        'clan': [G.nodes[n]['clan'] for n in G.nodes],
        'degree': [G.degree(n) for n in G.nodes],
    })
    nodes_df = annotate_nodes_with_groups(nodes_df, groups_df)
    nodes_df.to_csv(nodes_file, index=False)

    edge_rows = []
    for u, v in G.edges:
        su = G.nodes[u]['clan']
        sv = G.nodes[v]['clan']
        same = (su != 'NoClan') and (sv != 'NoClan') and (su == sv)
        edge_rows.append({
            'source': u,
            'target': v,
            'weight': G[u][v]['weight'],
            'evalue': G[u][v]['evalue'],
            'bitscore': G[u][v]['bitscore'],
            'alnlen': G[u][v]['alnlen'],
            'source_clan': su,
            'target_clan': sv,
            'same_clan': same
        })
    edges_df = pd.DataFrame(edge_rows)
    edges_df.to_csv(edges_file, index=False)

    cand_df.to_csv(candidates_file, sep=',', index=False)
    groups_df.to_csv(proposed_groups_file, sep=',', index=False)

    with open(summary_file, 'w') as fh:
        fh.write(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges\n")
        fh.write(f"Existing-clan candidates: {len(cand_df)}\n")
        fh.write(f"NoClan proposed groups: {len(groups_df)}\n")

    print(f"Done. Wrote results to {outdir}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True, help='Foldseek filtered file')
    ap.add_argument('--sep', default='\t', help='Field separator: "tab", "\t", ","')
    ap.add_argument('--weight', default='bitscore_per_res',
                    choices=['bitscore_per_res', 'neglog10_e', 'bitscore'])
    ap.add_argument('--min-neigh', type=int, default=1,
                    help='Min annotated neighbours for suggesting existing clan')
    ap.add_argument('--agreement', type=float, default=0.7,
                    help='Neighbour clan agreement fraction')
    ap.add_argument('--min-group-size', type=int, default=2,
                    help='Min size for a NoClan proposed group')
    args = ap.parse_args()

    analyse_foldseek_results(args.input, args.sep, args.weight, args.min_neigh, args.agreement, args.min_group_size)
