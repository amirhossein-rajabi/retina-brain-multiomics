# Network, regulatory and enrichment analyses

These steps were run with web-based tools. This folder contains their input and the
exported results, so that the downstream tables and figures can be traced without
re-running the tools.

| File | Tool | Content |
|---|---|---|
| `shared_genes_148.csv` | — | The 148 shared protein-coding genes (97 genomic, 52 transcriptomic; SEMA3C in both) with disease membership, source and direction of change. Input to STRING and g:Profiler. |
| `string_interactions.tsv` | STRING v12.5 | Interactions among the shared genes (full network, *Homo sapiens*, combined score ≥ 0.400, text mining excluded, no added interactors): 65 connected proteins, 92 edges. |
| `gephi_node_metrics_modules.csv` | Gephi 0.11.3 | Per-gene degree, weighted degree, betweenness, closeness, eigenvector centrality, clustering coefficient, Louvain module (resolution 1.0, weighted) and hub status for the 65 connected genes. |
| `omicsnet_tf_gene.csv` | OmicsNet 2.0 (TRRUST) | Transcription factor–target interactions of the 65 network genes. |
| `omicsnet_mirna_gene.csv` | OmicsNet 2.0 (miRTarBase v9.0) | miRNA–target interactions; `Id1` is the miRNA name and `Name1` its miRTarBase accession. |
| `omicsnet_metabolite_protein.csv` | OmicsNet 2.0 (KEGG) | Metabolite–protein interactions as exported, before filtering. |
| `metabolite_filtering.csv` | — | Every metabolite–protein interaction with its network module, filtering result and metabolite class. |
| `gprofiler_shared_genes.csv` | g:Profiler e114_eg62_p19_9dc92ac8 | Over-representation analysis of the 148 shared genes. |
| `gprofiler_modules.csv` | g:Profiler e114_eg62_p19_9dc92ac8 | Over-representation analysis of each network module with ≥ 4 genes (M1–M8); only modules with significant terms are listed. |

FKBP1C (genomic layer; AD and glaucoma) has no STRING interaction at the chosen
settings and no annotation in the g:Profiler sources, so it does not change the
network or the enrichment results.

## Settings

**STRING.** Active interaction sources: experiments, databases, co-expression,
neighbourhood, gene fusion, co-occurrence; text mining off. Minimum combined
score 0.400.

**Gephi.** Undirected graph weighted by the STRING combined score. Betweenness and
closeness were computed on the unweighted graph; closeness is interpreted within
the largest connected component. Hub genes rank in the top seven (≈ 10 % of 65
connected genes) by degree, betweenness and closeness simultaneously.

**OmicsNet.** The 65 connected genes were submitted by symbol. Three TF–gene
interactions whose target lies outside the network (RUNX2 → TRPS1, TWIST1, TWIST2)
were excluded from the analysis.

**Metabolite filtering.** KEGG links a protein to every substrate and product of its
reactions. Three classes of compounds were removed:

1. currency metabolites and cofactors (H₂O, H⁺, O₂, ATP, ADP, NAD⁺, NADH, NADP⁺,
   NADPH);
2. non-specific KEGG entries that denote a compound class (e.g. Protein, RH, Fatty
   acid);
3. xenobiotics and drugs (e.g. tamoxifen, aflatoxins, benzo[a]pyrene, caffeine,
   lidocaine).

Endogenous substrates and products of the cytochrome P450 enzymes (steroid
hormones, eicosanoids, retinoids, melatonin) were retained. Of 259 interactions,
93 interactions between 54 metabolites and 8 proteins were retained.

**g:Profiler.** Unordered queries; sources GO (MF, BP, CC), KEGG, Reactome and
WikiPathways; custom background of 24,238 genes (HGNC protein-coding genes plus
all genes tested in the five transcriptomic datasets), annotated genes only; g:SCS
correction, adjusted *P* < 0.05; terms with 10–500 genes.
