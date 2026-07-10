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

**Padding mutations are merged out first.** Each displayed protein is
`prefix_pad + core + suffix_pad`, where the pads are flexible GS-linker sequence that does not touch
the ligand. A mutation in the padding cannot change binding, so before anything else we **collapse
any read whose only differences from WT lie in the padding onto the core genotype** (summing counts).
Using `design_unpadded_mapping.tsv`, a 1-based mutation at position `pos` is padding iff
`pos ≤ len(prefix_pad)` or `pos > len(sequence) − len(suffix_pad)` (verified against real WT
sequences for all 12,000 designs; the per-row sequence length is authoritative — `padded_len` and
`unpadded_seq` are off-by-one for 482 "terminal-overwritten" designs and are not used). Surviving
**core mutations are renumbered to unpadded coordinates** (`pos − len(prefix_pad)`), so what was
`S13Y` on the padded sequence appears as its core-relative position. This removed **~8,500 mutant
genotypes (25%)** that were padding-only, and notably reclassified `EST · S15F` — previously a
"valid hit" — as a **padding mutation** (position 15 sits in an 18-residue prefix pad): its 8,344
Sort-2 reads correctly fold into that design's WT.

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
  > overstate a genotype's true panel presence. Example: the `r5_470…` V25I mutant's expression/S1
  > reads live entirely on that off-panel cluster; it reports `4 raw / 2 panel` at expression, where
  > the 2 panel members are the tie-broken nearest, not barcodes V25I is genuinely an error of. The
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

Of the 8,800 expressed families, **7,023 (80%) are WT-only** after the padding merge — the
sequencing never picked up any *core* mutant of them. Only **1,781 families carry at least one core
mutant branch**, and those are where the valid-vs-artifact question lives.

## The headline answer: extra sequences are mostly parent artifacts

Across all 8,800 families, after merging out padding mutations there are **26,127 mutant genotypes**
(the "extra sequences on top of what we ordered"; down from 34,612 before the padding merge).
Classified by whether each one behaves independently of its WT parent *and* enriched through selection:

| Class | count | share | meaning |
|---|---:|---:|---|
| **Valid extra hit** | **3** | 0.01% | independent of the WT parent, enriched through selection, AND not just a low-% read shadow of the parent |
| Ambiguous | 4,275 | 16% | some panel independence but failed the selection, read, or parent-% floors |
| **Parent artifact** | **21,849** | **84%** | occupies only the WT's panel barcodes, or is <5% of the parent's reads on shared barcodes |

So the direct answer to *"how much of these extra sequences are valid hits rather than byproduct
artifacts of the lineage parent"* is: **~84% are artifacts, ~16% are ambiguous, and only 3 mutant
genotypes are confident valid hits** (positions in core/unpadded coordinates): `DM2 · L37F` (a real
affinity-matured binder), `EST · V25I` (74% of the estradiol Sort-2 pool), and `EF2 · A127V`.
(`EST · S15F`, a former "hit", turned out to be a **padding mutation** and folded into its WT.)

`DM2 · S13Y` — previously called a binder — is now correctly demoted: although it enriched, on the
barcodes it shares with its (binder) WT parent it is only **1.57%** of the parent's reads, a level
consistent with a fixed sequencing-error rate. Its apparent enrichment is most likely the *parent's*
enrichment bleeding onto shared barcodes. `L37F` survives the same test at **27%** of its parent —
too abundant to be error bleed, so a genuine variant.

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
- **binder** — clears the enrichment cutoffs (r1 ≥ 5×, r2 ≥ 2×) **and** is >5% of its WT parent's
  reads on the barcodes they share → strongest evidence, valid. (The parent-% test guards against a
  mutant that enriched only because the *parent* enriched and bled onto shared barcodes — a mutant
  at ~1% of an abundant binder parent is error-rate-consistent and is demoted to artifact.)
- **likely_valid** — must clear **all** of: (a) **itself enriched through selection** — net
  enrichment from expression to Sort 2 with a strong final-round gain (the *same* selection standard
  the WT/binders pass, not a rigid round-1 ≥5× — binders often dip in round 1 then explode); (b)
  holds **≥1 true-panel barcode its WT parent does not occupy**; (c) sits on **≥2 of the 27 true-panel
  barcodes**; (d) **>100 Sort-2 reads** and **>100 reads on the independent panel positions**.
  A mutant that merely rides an abundant parent but *shrinks* each round fails (a) and is not valid.
- **likely_artifact** — **every** panel barcode it maps to is also the WT's (`panel_only = 0`); it
  holds no panel position of its own → bleed-through of the abundant parent.
- **ambiguous** — has some panel independence but below the 50-read floor.

> **Why panel-level, not raw.** A mutant like `vi_rank258_path3_tbg · A22S` sits entirely on the
> off-panel `GATTCCGA`-cluster, which *raw* comparison made look "unique" (its raw strings differed
> from the WT's), so it was mis-called valid. But those barcodes force-map to the *same* 2 panel
> positions the WT already occupies — `panel_only = 0` — so it's correctly an artifact. Judging
> independence on the 27-panel fixes this whole class of false positives.

The standout valid hit is **`V25I` on `r5_470_ttr_ef2_…`**: it appears across all six targets with
only ~25% of its barcodes shared with the parent, and in estradiol it reaches **1.27M Sort-2 reads**
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
2. **Pool depths are uneven — always compare % of pool, not raw reads.** The pools were sequenced
   to different matched depths (expression ≈ 2.14M matched reads; each Sort 1 ≈ 1.0–1.6M; each
   Sort 2 ≈ 0.67–1.6M), and **Sort 1 is shallower than expression for every ligand**. So a genotype
   can show *fewer raw reads* in a pool while being a *larger share* of it. Every metric uses the
   genotype's **percentage of that pool's matched reads** (`reads / pool_total`), where `pool_total`
   is the sum of all per-genotype reads — so **within a pool the genotype percentages sum to exactly
   100%**. (The raw sequencer total is *not* used as the denominator: ~31–59% of raw reads never
   matched a design, so shares would sum to <100% and a "% of pool" would be incoherent.) Enrichment
   ratios (`enrich_s1`, `enrich_s2`, net) are ratios of these fractions; the flow diagram sizes dots
   by **% of pool** and the drawer shows a `% of pool` column next to raw reads. Example: `EST · V25I`
   climbs 0.048% (expr) → 0.040% (S1) → **74.30% (Sort 2)** — it becomes three-quarters of the whole
   estradiol Sort-2 pool.
3. Same shared-panel caveat as before: barcode counts are an independent-observation signal, not
   per-molecule tags.
4. **Read floors are pre-applied.** The upstream counter (`count_ngs_mutations.py`) already
   dropped barcode rows below ~5 reads and required ≥2-barcode support at the sequence level, so
   the per-node barcode counts here are **floored, not raw** — very low-abundance observations of
   a genotype were removed before this build ever saw them. (The counter also confirmed zero
   `ambiguous` rows in these files, so no design assignments were tie-broken.)
