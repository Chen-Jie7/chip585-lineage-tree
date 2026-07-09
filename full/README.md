# Full lineage tree — every expressed family, valid hits vs parent artifacts

*Companion to `lineage_tracing_explained.md`. This one covers the **full** tree (all
expressed designs, not just the top-15 per target) and the specific question: of the extra
sequences we did **not** order, how many are real hits vs. byproducts of an abundant parent?*

Interactive view: `lineage_full_tree.html` (self-contained; published as an Artifact).

---

## The node model you asked for

A **node = (design, mutation, barcode)** — a specific protein genotype observed on one of the
27 shared panel barcodes. A **family** = one ordered WT design plus all its accidental point
mutants. Within a family the **WT (mutation = "") is the trunk**, each distinct mutation is a
**branch**, and each panel barcode is a leaf observation of that genotype. Every node carries its
read count and barcode count in each of the three pools: **expression → Sort 1 → Sort 2**.

The 27 barcodes are a **fixed shared panel** (pool/sample indices reused across all 12k designs),
*not* per-design tags — so barcodes do not define families; the protein sequence does. Barcodes
that sit one nucleotide off a panel member (~8% of reads) are sequencing errors and are
**collapsed onto their panel member** before counting.

**Two barcode counts per genotype.** Because raw barcode strings include those sequencing errors,
each genotype carries two counts in every pool:
- **raw barcodes** — distinct barcode strings exactly as sequenced.
- **true-panel barcodes (of 27)** — how many members of the real panel those map to. Each observed
  barcode is assigned to its **nearest panel member, ties broken by panel index** (this is a
  *forced* assignment: every barcode maps to something, regardless of distance). The count is
  bounded by 27 and dedupes error-barcodes, so a genotype seen on `20 raw / 15 panel` barcodes was
  on 15 distinct panel positions.

  > **Caveat on the forced mapping.** Node *identity* still collapses only clean, unambiguous 1-nt
  > error barcodes onto the panel — that part is conservative. But the reported panel *count* uses
  > nearest-member-wins even for far/ambiguous barcodes. A distinct off-panel population exists here
  > (a `GATTCCGA`-type cluster ~3 nt from the panel, ~6% of reads — the H3 spike in the distance
  > histogram, not sequencing error of any single panel member). Under the forced rule those are
  > folded onto their lowest-index nearest panel member (e.g. `GATTCCGA → GACACTGA`), which can
  > overstate a genotype's true panel presence. Example: the `r5_470…` V45I mutant's expression/S1
  > reads live entirely on that off-panel cluster; it reports `4 raw / 2 panel` at expression, where
  > the 2 panel members are the tie-broken nearest, not barcodes V45I is genuinely an error of. The
  > **classification heuristic still runs on the honest raw barcode footprint**, so this reporting
  > choice does not change valid-vs-artifact calls.

**Barcode sharing with the WT parent.** For each mutant we also report how many of its true-panel
barcodes are also carried by its **same-design WT** (the parent is unambiguous — a mutant of design
D descends from D's wild-type). A mutant whose panel barcodes are *all* shared with an abundant WT,
with no independent sort signal, is the bleed-through artifact. (Note: a real affinity-matured
binder like DM2·S13Y also co-occurs on the parent's barcodes — there the *enrichment* is what
proves it real, not barcode independence.)

## Coverage — what actually expressed

| | count | of 12,000 |
|---|---:|---:|
| Ordered designs | 12,000 | 100% |
| **Expressed** (≥1 read in the expression pool) | **8,800** | **73%** |
| Never detected | 3,200 | 27% |

Of the 8,800 expressed families, **6,663 (76%) are WT-only** — the sequencing never picked up any
mutant of them. Only **2,141 families carry at least one mutant branch**, and those are where the
valid-vs-artifact question lives.

## The headline answer: extra sequences are mostly parent artifacts

Across all 8,800 families there are **34,612 mutant genotypes** (the "extra sequences on top of
what we ordered"). Classified by whether each one behaves independently of its WT parent:

| Class | count | share | meaning |
|---|---:|---:|---|
| **Valid extra hit** | **11** | 0.03% | ≥2 true-panel barcodes incl. ≥1 the WT lacks, >100 Sort-2 reads & >100 indep. reads — or cleared binder cutoffs |
| Ambiguous | 5,644 | 16% | some panel independence, below the read/barcode floors |
| **Parent artifact** | **28,957** | **84%** | occupies only panel barcodes the WT also has — no independent panel position |

So the direct answer to *"how much of these extra sequences are valid hits rather than byproduct
artifacts of the lineage parent"* is: **~84% are artifacts, ~16% are ambiguous, and only 11 mutant
genotypes (0.03%) are confident valid hits** — 2 of which are the already-known affinity-matured
binders `DM2 · S13Y` and `DM2 · L52F`, the rest dominated by the recurrent `V45I` variant (seen
independently across all six targets, reaching 1.19M Sort-2 reads in estradiol). Judging independence
on the 27-panel footprint rather than raw barcodes is what makes this honest — see the A39S example below.

## How a mutant is judged valid vs. artifact

The bleed-through artifact (caveat #2 in the original explainer) has a clear signature at
(design, mutation, barcode) granularity: a mutant whose reads appear **only on barcodes where its
abundant WT parent is also present**, at levels consistent with a low per-read error rate, is
almost certainly sequencing error of the parent — not a real variant. The heuristic per mutant:

- **`frac_on_parent_bc`** — fraction of the mutant's barcodes that also carry the WT (high ⇒ suspicious).
- **`mut_only_bc`** — barcodes carrying the mutant but not the WT (its own independent footprint).
- **`independent_sort_reads`** — reads the mutant accumulates through S1/S2 on those own barcodes.

Classification — **independence is judged on the true-panel footprint, not raw barcodes** (raw
strings fragment the off-panel error cluster and can make a mutant look falsely independent):
- **binder** — clears the enrichment cutoffs (r1 ≥ 5×, r2 ≥ 2×) → strongest evidence, valid.
- **likely_valid** — holds **≥1 true-panel barcode its WT parent does not occupy**, sits on **≥2 of
  the 27 true-panel barcodes**, and has **>100 Sort-2 reads** with **>100 reads on the independent
  panel positions** → grew independently of the trunk with real read support.
- **likely_artifact** — **every** panel barcode it maps to is also the WT's (`panel_only = 0`); it
  holds no panel position of its own → bleed-through of the abundant parent.
- **ambiguous** — has some panel independence but below the 50-read floor.

> **Why panel-level, not raw.** A mutant like `vi_rank258_path3_tbg · A39S` sits entirely on the
> off-panel `GATTCCGA`-cluster, which *raw* comparison made look "unique" (its raw strings differed
> from the WT's), so it was mis-called valid. But those barcodes force-map to the *same* 2 panel
> positions the WT already occupies — `panel_only = 0` — so it's correctly an artifact. Judging
> independence on the 27-panel fixes this whole class of false positives.

The standout valid hit is **`V45I` on `r5_470_ttr_ef2_…`**: it appears across all six targets with
only ~25% of its barcodes shared with the parent, and in estradiol it reaches **1.19M Sort-2 reads**
— a genuine independent variant, not bleed-through.

## Files produced

| File | What it is |
|---|---|
| `build_lineage_full.py` | the pipeline: collapse error-barcodes, build nodes, classify, emit outputs |
| `build_lineage_viz.py` | renders the self-contained interactive HTML from the tree JSON |
| `lineage_full_tree.html` | the interactive tree (this is the visualization) |
| `ngs_counts_out/lineage_full_nodes.csv.gz` | **every** (ligand, design, mut, barcode) node with expr/S1/S2 reads — 592,223 rows |
| `ngs_counts_out/lineage_mutant_calls.csv` | one row per mutant genotype with the artifact heuristic + classification |
| `ngs_counts_out/lineage_family_summary.csv` | one row per family: WT survival, mutant counts, per-pool totals |
| `ngs_counts_out/lineage_full_tree.json` | pruned tree (top-40 families + all valid-hit families per target) driving the HTML |
| `ngs_counts_out/lineage_coverage.json` | the headline coverage numbers |

## Caveats

1. The interactive tree shows the **top-40 families per target plus every family with a valid
   hit**; within a family it shows WT + all valid mutants + a sample of artifacts, with a hidden
   count. The **full data is in the CSVs** — nothing is discarded, only hidden from the picture.
2. "Valid" here means **independent of the parent**, not necessarily **enriched**. A few
   likely_valid mutants have their own barcodes but actually *dropped* through Sort 2 (e.g. an
   `E109D` that fell from 68k to 66 reads). The drawer shows per-round enrichment so you can tell
   independent-and-real from independent-but-dying.
3. Same shared-panel caveat as before: barcode counts are an independent-observation signal, not
   per-molecule tags.
4. **Read floors are pre-applied.** The upstream counter (`count_ngs_mutations.py`) already
   dropped barcode rows below ~5 reads and required ≥2-barcode support at the sequence level, so
   the per-node barcode counts here are **floored, not raw** — very low-abundance observations of
   a genotype were removed before this build ever saw them. (The counter also confirmed zero
   `ambiguous` rows in these files, so no design assignments were tie-broken.)
