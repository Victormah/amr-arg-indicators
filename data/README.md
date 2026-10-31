# Input data

The pipeline needs three files in this folder: `metadata.tsv` and `ARG.tsv` (public, not
redistributed here; download as described below) and the CARD mapping table (included).

| File | Content | Source |
|---|---|---|
| `metadata.tsv` | Sample metadata (host label, project, country, read counts) | 214K-metagenome resource (Martiny et al., 2022) |
| `ARG.tsv` | ARG alignment counts (ResFinder references) per sequencing run | 214K-metagenome resource (Martiny et al., 2022) |
| `aro.tsv` (optional) | CARD antibiotic resistance ontology; needed only to rebuild the mapping table | CARD, https://card.mcmaster.ca/download (release used: downloaded 23 October 2026) |
| `AMRgenelist_mapped_CARD_ALL.csv` | ResFinder gene names mapped to CARD terms (ARO accession, name, drug class, mechanism, gene family) | **included in this repository** (CARD free-text descriptions removed; they are not used in any calculation). Can be rebuilt with `scripts/build_card_mapping.R` |

**214K-metagenome resource**

- Martiny H-M, Munk P, Brinch C, Aarestrup FM, Petersen TN. A curated data resource of 214K
  metagenomes for characterization of the global antimicrobial resistome. *PLoS Biology*
  2022;20:e3001792. https://doi.org/10.1371/journal.pbio.3001792
- Data (version 3): Zenodo, https://doi.org/10.5281/zenodo.7553776. Download the metadata table
  and the ResFinder ARG table and save them here as `metadata.tsv` and `ARG.tsv`.

**Smaller alternative to `ARG.tsv`.** `ARG.tsv` is large (about 425 MB). The pipeline also accepts
`ARG_eligible_5sources.tsv`, the same table restricted to the 2,758 samples labeled bird, compost,
feces, food or wastewater metagenome. It gives identical results and loads much faster; it can be
made from `ARG.tsv` by keeping rows whose sample appears in `metadata.tsv` with one of these five
host labels.

Only samples whose metadata host field is exactly one of "bird metagenome", "compost metagenome",
"feces metagenome", "food metagenome" or "wastewater metagenome" are used.
