# Mouse cancer cell-line companion app — research plan

_Written 2026-05-11 (Stockholm) after a web sweep of available mouse cancer cell-line resources. The goal: figure out whether a Correlate-style app for mouse cell lines is feasible, and if so what data sources to build on._

## ⚑ Architecture decision (2026-05-12)

**Scope tightened to a standalone Mouse Cell Line Browser app — NOT a fork of Correlate V2.**

The new app is its own product: same visual language and CLB layout as V2 (so it feels like a sibling alongside Correlate / CoExpress), but no dependency on the V2 codebase, no shared web data, no analysis modes (network / correlation / cluster / AI export), no menu link from V2 to mouse or vice versa. Stand-alone repo, stand-alone deploy, stand-alone landing page.

Working name: **MouseCLB** (placeholder; final name TBD).

**Why standalone and not a fork:**
- Roughly 60 % of V2's code drives analysis flows (network rendering, correlation engine, AI export, drug-response section, CRISPR-dependency view, gate-comparison plots). None of that applies for mouse — those data layers don't exist at scale (CRISPR / PRISM) or aren't in the MCCA download (drug response, fusions). Forking would just mean carrying dead code.
- The CLB itself is a self-contained part of V2 — the detail card, Wiki modal, collections panel, sort, list, oncoprint-style filters. Lifting just that as a new app is cleaner than maintaining one bloated codebase with feature flags.
- Separate repo means the mouse app can iterate on schema / data format without touching V2.
- The mouse data has a different reference genome (mm10), different gene-ID convention (ENSMUSG / MGI), different lineage ontology, different driver gene set. Mixing into V2 would entangle these.

**What carries over (visual / interaction language):**
- The 3-block detail card (Identity / specific features / Wiki button).
- The Wiki modal with its sectioned layout, executive summary, distribution histograms, large green section headers.
- The oncoprint-style collections panel with ✓/✗ include/exclude, the active-filter chip strip.
- The sort dropdown with inline values + caption.
- The plain-English low/medium/high tier descriptors on numeric metrics.
- The gene-hover MyGene-style tooltips (using MGI for mouse).
- The export buttons on the genome-metric histograms (PNG / SVG / CSV).

**What gets dropped vs V2:**
- Network / correlation / cluster tabs (the whole `tab-network`, `tab-correlations`, `tab-clusters` machinery).
- Mutation-analysis tab, including hotspot picker, Welch's t-test, compare-by-tissue / hotspot / fusion tables.
- AI JSON export.
- Synonyms / orthologs mode.
- Gene effect distribution / scatter inspect modals.
- PRISM Drug response section in the Wiki (no data).
- CRISPR dependencies section in the Wiki (no data).
- Druggable dependencies / Pathway dependencies interpretive panels (depend on CRISPR).

**What needs new design** (different in mouse):
- Lineage / tissue ontology. MCCA's `Tissue` and `MouseModel` columns instead of Oncotree.
- Driver-gene panels. Need mouse-orthologue versions of the cancer pathway panels (TP53 → Trp53, KRAS → Kras, etc.) via the MGI homology table.
- Hotspot definitions. Mouse cancer hotspots aren't catalogued the same way — Trp53 R172H is the canonical mouse-tumour hotspot vs human TP53 R175H. May need a small curated mouse-hotspot list.
- Authentication. STR profiles are sparse for mouse (141 lines) — the STR section becomes optional with a Cellosaurus link as the primary handle.

**File layout (suggested):**
```
~/Documents/MouseCLB/
  index.html              # single-page app
  app.js                  # all logic; mirrors V2 patterns but smaller
  web_data/
    metadata.json         # cell-line metadata (lineage, tissue, sex, strain, model type, …)
    mutations.json        # hotspot + damaging matrices (derived from MCCA-Mutations xlsx)
    cn.json               # gene-level CN matrix (segments → gene overlay)
    expression.bin.gz     # log2-like expression matrix (from VST batch-corrected txt)
    expression_genes.json # gene-symbol map (ENSMUSG → symbol via MGI)
    pathway_panels.json   # mouse-orthologue pathway gene panels
    facs_markers.json     # mouse surface-antigen panel (orthologue of V2 FACS list)
    cellosaurus_rrid.json # MCCA-ID → Cellosaurus CVCL_ mapping
  scripts/
    process_mcca_metadata.py
    process_mcca_mutations.py
    process_mcca_cn.py
    process_mcca_expression.py
    process_pathway_panels.py
    fetch_mgi_homology.py
  README.md
  DESIGN_NOTES.md
```

The sections under this point describe the original research sweep — kept for reference but the architecture above supersedes the "fork V2" plan in the "Suggested architecture" section below.

---


## TL;DR

- **Feasibility: yes, but smaller.** A new mouse companion app is realistic for **mutations + transcriptomics + curated metadata + (separately) immuno-oncology overlay**.
- **The single most important resource is brand-new (2025/26): the Mouse Cancer Cell Line Atlas (MCCA)** — 590 cell lines across 22 lineages / 46 cancer types, with multilayered genomic, transcriptomic, morphological and histopathological data. Published in _Nature Genetics_ ("Mapping cancer evolution with the mouse cell line atlas"). Portal at **www.mcca.tum.de** (was unreachable from here when I checked — likely transient TUM infrastructure issue).
- **No mouse equivalent of CRISPR dependency screens** (DepMap-style) and **no mouse equivalent of PRISM-scale drug screens**. Those are the two layers we would have to drop, or scrape from individual papers.
- **Best secondary resource: TISMO** (Tumor Immune Syngeneic MOuse) — ~49 syngeneic cell lines + 68 in-vivo tumour models with RNA-seq, immune infiltration, cytokine and checkpoint-blockade treatment context. http://tismo.cistrome.org. Tightly focused on immuno-oncology; ideal as an immunology overlay for the subset of cell lines that overlap with MCCA.
- **Cellosaurus** indexes 31,079 mouse cell lines but only 141 have STR profiles. Good as a registry / authentication link, weak as a data source.
- **Taconic Syngeneic Cell Line Reference Database** — useful as a curated catalogue of strain compatibility, immunogenicity rating, metastatic potential. Free, web-only table; would need to scrape.

## Coverage map: what we have for humans vs what we can get for mouse

| Data layer (human DepMap → Correlate V2) | Mouse equivalent | Status |
| --- | --- | --- |
| Hotspot + damaging mutation matrix | **MCCA WES/WGS** (assumed; paper says "genomics") | Likely available |
| Copy number | **MCCA** (likely WGS-derived) | Likely available |
| Gene expression (TPM) | **MCCA** RNA-seq (and **TISMO** for syngeneic subset) | Available |
| CRISPR gene effect (Chronos) | **None at scale** | Major gap — individual papers only |
| PRISM drug response | **None at scale** | Major gap — individual papers only |
| Fusions (clinical + raw) | Likely from MCCA WGS | Need to verify |
| Cell-line metadata (lineage, subtype, sex, age) | MCCA clinical metadata | Available |
| STR authentication | **Cellosaurus** — only 141/31,079 mouse lines | Sparse |
| Global signatures (WGD, CIN, aneuploidy, MSI) | Likely derivable from MCCA WGS | Need to compute |
| Inferred functional loss (CN + mut + expr) | Re-derivable using the same logic from MCCA inputs | Build ourselves |
| Lehmann TNBC / disease subtyping | Different in mouse — would need a separate effort | Out of scope initially |
| MyGene.info hover info | MGI for mouse | Different API, easy to add |

## Primary candidate resources

### 1. Mouse Cancer Cell Line Atlas (MCCA) — _Nature Genetics_ 2025/26
- **What:** 590 mouse cancer cell lines, 22 lineages, 46 cancer types. Genomics + transcriptomics + morphology + histopathology + clinical metadata. Computational methods provided for inferring mouse immunophenotype from genomic data (to guide immunocompetent transplantation).
- **Why this is the right backbone:** explicitly designed as a multi-omic mouse cell-line atlas, with cross-species comparison to human cancers as a core use case. Closest thing in spirit to DepMap for mouse.
- **Portal:** www.mcca.tum.de (interactive web portal — needs verification of bulk-download paths).
- **Paper:** [Mapping cancer evolution with the mouse cell line atlas — _Nature Genetics_](https://www.nature.com/articles/s41588-026-02589-9)
- **Open questions:** raw file formats, license terms, whether CRISPR-style dependency data is included (paper summary doesn't mention it), whether per-cell-line drug screens are included.

### 2. TISMO — Tumor Immune Syngeneic MOuse
- **What:** 605 in-vitro RNA-seq samples from 49 syngeneic cell lines (23 cancer types) + 1,518 in-vivo samples from 68 syngeneic tumour models (19 cancer types). Treatment metadata covers IFNγ / IFNβ / TNFα and anti-PD1 / anti-PDL1 / anti-CTLA4. Six immune-cell-infiltration deconvolution algorithms applied.
- **Why useful:** the immuno-oncology overlay. Perfect for "predicted ICB response" features that mirror our human Immunology category.
- **Portal:** http://tismo.cistrome.org (R-Shiny, free)
- **Paper:** [TISMO: syngeneic mouse tumor database to model tumor immunity and immunotherapy response — _NAR_](https://academic.oup.com/nar/article/50/D1/D1391/6371975)
- **Data:** downloadable from the portal — expression matrices, immune infiltration estimations, metadata.
- **Caveat:** no mutations as primary layer. Combine with MCCA for those.

### 3. Cellosaurus
- **What:** 31,079 mouse cell lines indexed. Each entry has provenance, parent / derivative relationships, contamination flags, RRID.
- **Use:** authentication link from each cell-line page (same as the V2 STR-authentication block), and the source-of-truth for cell-line identity / RRID lookups.
- **Caveat:** only 141 mouse lines have STR profiles (vs almost-universal coverage for human DepMap lines). Need to make the absence of an STR profile a non-error state.
- **Reference:** [Cellosaurus](https://www.cellosaurus.org/)

### 4. Taconic Syngeneic Cell Line Reference Database
- **What:** curated table of syngeneic mouse cancer cell lines with strain compatibility, immunogenicity rating (---/+++ scale), metastatic potential, and PubMed citations.
- **Use:** quick "what mouse strain do I transplant this into" context. Could be a sidebar in the per-cell-line view.
- **Format:** web table, would need to scrape. Free.
- **Reference:** [Taconic Syngeneic Cell Line Reference Database](https://www.taconic.com/resources/databases/syngeneic-cell-line-reference-database)

### 5. MGI Gene Expression Database (GXD)
- **What:** mouse gene expression atlas at the gene-level, not cell-line-level.
- **Use:** **MyGene.info equivalent for mouse genes** — gene-info popups on hover in the app. Different API endpoint (MGI), same UX.
- **Reference:** [MGI GXD](https://www.informatics.jax.org/expression.shtml)

### 6. MMHCdb (Mouse Models of Human Cancer Database) at JAX
- **What:** curated knowledgebase of mouse models (inbred strains, GEMMs, PDXs, diversity panels). Links to GEO for genomic data.
- **Use:** mostly for in-vivo tumour models, not cell lines per se — but useful for cross-referencing a cell line to the GEMM it derived from.
- **Reference:** [MMHCdb](https://tumor.informatics.jax.org/)

### 7. MTB (Mouse Tumor Biology) at JAX
- **Same group as MMHCdb.** Confirmed: **does not currently include cell-line data** (in-vivo tumours only). Not useful for our app.
- **Reference:** [MTB](https://www.jax.org/jax-mice-and-services/preclinical-research-services/mouse-tumor-biology-database)

## Major gaps (and what to do about them)

1. **No mouse CRISPR dependency screen at DepMap scale.** Genome-wide pooled CRISPR has been done in syngeneic models (e.g. PTPN2 paper, mGeCKOa library) but as individual studies, not a unified ~1000-cell-line resource. **Mitigation:** drop the "Top uniquely essential genes" / "Pathway dependencies unique to this line" / "Oncogene addiction (mutation × CRISPR dependency)" features in the mouse app. Or scrape individual papers and compose a small CRISPR view; expect &lt; 30 lines covered.

2. **No mouse PRISM-scale drug screen.** Many individual drug-screen papers exist (NCI panels, single-lab efforts) but no unified resource. **Mitigation:** drop the PRISM Drug response section in the mouse app, or rebuild it from a single-paper supplement (limited coverage).

3. **Sparse STR profiles.** Only 141 of 31,079 lines have STR. **Mitigation:** STR section becomes "if available, show — else just the Cellosaurus link" for authentication.

4. **No mouse-specific TNBC / Lehmann-style molecular subtypes.** Mouse mammary tumour models have their own subtypings (e.g. MMTV-PyMT model classes) but it's not a Lehmann analogue. **Mitigation:** initially skip subtyping; revisit if MCCA provides their own subtyping.

## Suggested architecture

Build a **third sibling app** alongside Correlate V1 (gene effect on humans), CoExpress (gene expression on humans), and Correlate V2 (V1 + expression correlates). Call it something like **Correlate-M** or **MouseCorrelate**.

Reuse the V2 codebase as a starting point and:

- **Keep:** detail card, collections panel (oncoprint-style include/exclude), Wiki modal, executive summary, mutation analysis mode, expression analysis mode, lineage / subtype filters, genome distribution histograms.
- **Adapt:**
  - Replace DepMap human data with MCCA mouse data.
  - MyGene.info → MGI for gene-info popups.
  - Cellosaurus link → mouse-RRID-aware.
  - Hotspot / damaging definitions — mouse uses a different reference (GRCm39 mm39); driver-gene panel should be the mouse orthologues of the human cancer pathway panels.
  - Lineage / subtype labels — mouse uses MTB ontology not Oncotree.
- **Drop:**
  - CRISPR-dependency section + all derived features (oncogene-addiction collections, "Top uniquely essential" Wiki block, dependency-based pathway status).
  - PRISM Drug response section.
  - Lehmann TNBC subtypes.
- **Add:**
  - TISMO immunology overlay: immune-cell infiltration estimates and IFN / ICB treatment-response data for the syngeneic subset (~49 lines).
  - Cross-species comparison view: pick a mouse line, see the closest human CCLE / DepMap line by orthologous-mutation overlap. MCCA explicitly enables this; use their computational methods.

## Step-by-step implementation plan (rough)

### Phase 0 — Verify data access (1 day)
- Get www.mcca.tum.de working — likely transient. If portal is unstable long-term, fall back to the paper's supplementary tables or the underlying repository (TUM datacenter or GEO/SRA accessions referenced in the paper).
- Document exact file formats, schema, license terms.
- Identify the equivalent of DepMap's `Model.csv` (cell-line metadata), `OmicsSomaticMutationsMatrixHotspot.csv`, `OmicsCNGeneWGS.csv`, `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`, `OmicsFusionFiltered.csv` etc. — note that mouse may bundle them differently.

### Phase 1 — Data processing pipeline (3–5 days)
- New `scripts/` directory mirroring V2's processing scripts.
- `process_mcca_metadata.py` — build `cellLineMetadata.json`.
- `process_mcca_mutations.py` — build hotspot + damaging matrices.
- `process_mcca_cn.py` — build CN matrix; reuse V2's relative-CN format if possible.
- `process_mcca_expression.py` — build expression matrix (RNA-seq TPM).
- `process_mcca_fusions.py` — build fusion data; clinical-fusion validation needs a mouse-driver-fusion panel (curate from literature — most mouse models use the human ortholog of the human fusion, e.g. PML-RARA → mouse Pml-Rara).
- `process_tismo.py` — load TISMO RNA-seq + immune infiltration as a per-cell-line layer for the subset that overlaps MCCA by cell-line name (match Cellosaurus RRIDs).
- `process_mouse_pathway_panels.py` — convert human `_WIKI_PATHWAYS()` gene panels to mouse orthologues via MGI homology table.

### Phase 2 — App fork + minimal viable view (3–5 days)
- Fork V2 repo, point at new `web_data/` directory.
- Strip CRISPR-dependency code paths (replace with empty / "data not available for mouse" notes).
- Strip PRISM code paths.
- Strip Lehmann TNBC code paths.
- Update Wiki executive summary's variant-aware driver mapping for mouse — many oncogenes have mouse-specific hotspots (e.g. mouse Hras G12 vs human KRAS G12; mouse Trp53 R172H vs human TP53 R175H).
- Adapt MyGene.info to MGI gene-info endpoint.

### Phase 3 — TISMO immunology overlay (1–2 days)
- For each syngeneic cell line in MCCA that overlaps TISMO: add an "Immuno-oncology context" Wiki block showing TISMO immune-cell infiltration estimates and ICB-response group.
- Collections: "ICB-responder profile (TISMO)", "ICB-resistant profile (TISMO)".

### Phase 4 — Cross-species comparison view (1–2 days)
- For each mouse cell line, compute and display its closest human DepMap line by orthologous-mutation overlap.
- Link out to the human Correlate V2 detail page for that line.

### Phase 5 — Polish (1–2 days)
- Authentication block: Cellosaurus link + caveat about STR sparsity.
- Mouse-specific catalog descriptions in the collections panel (replace human-only therapy mappings).
- Standard QA: try ~10 well-known mouse lines (B16-F10, 4T1, CT26, MC38, EMT6, LL/2, MOC1, Pan02, B16, EL4) and verify the views read sensibly.

### Phase 6 — Open questions to defer (revisit later)
- Whether to attempt to build a "Top depleted genes" CRISPR view from individual mouse-tumour papers. Most likely not worth the effort unless several large datasets emerge.
- Mouse mammary tumour subtyping (PyMT class, etc.) — only if it turns out to be useful to a critical user.
- Cross-species drug-response prediction: hard problem; outside the scope of this initial build.

## Caveats and risks

- **MCCA portal stability:** I couldn't reach www.mcca.tum.de via the redirect from here. If the portal is intermittently down, depend on the published supplementary tables and the GEO/SRA deposit instead.
- **License / re-use:** the paper is _Nature Genetics_, which usually means CC-BY for the article and "freely available for non-commercial research use" for data — confirm before redistributing.
- **Schema drift:** MCCA may not use the same per-file structure as DepMap. The processing scripts will need real engineering — don't assume drop-in.
- **TISMO syngeneic subset is small** (49 lines). The immunology overlay will only fire for a fraction of the panel; design the UI accordingly (no big empty section when TISMO data isn't available).
- **Mouse-vs-human gene orthology:** not always 1:1 (e.g. mouse Trp53 ↔ human TP53 with the standard convention; mouse Kras ↔ human KRAS; but several gene families differ in copy number or have mouse-specific paralogs). Use MGI homology table as the ground truth.

## Bottom line

A mouse companion app is worth doing. **MCCA is the unlock** — without it the mouse landscape was too fragmented for a unified app, but with 590 well-characterised cell lines and harmonised multi-omics it's now achievable. Expect ~2 weeks for a minimum-viable version, mostly spent on the data-processing pipeline rather than UI (the UI inherits from V2). The CRISPR dependency and drug-response features should be cut — they don't exist for mouse at scale and faking them with partial data would be worse than not having them.

Next concrete step: get MCCA portal back up (or its supplementary data) and verify the actual file formats / license terms.

---

### References

- [Mapping cancer evolution with the mouse cell line atlas — _Nature Genetics_ 2025/26](https://www.nature.com/articles/s41588-026-02589-9)
- [TISMO: syngeneic mouse tumor database — _NAR_](https://academic.oup.com/nar/article/50/D1/D1391/6371975)
- [TISMO portal](http://tismo.cistrome.org)
- [Cellosaurus](https://www.cellosaurus.org/)
- [Cellosaurus FAQ on mouse STR coverage](https://www.cellosaurus.org/faq)
- [ICLAC Guide to Mouse Cell Line Authentication](https://iclac.org/wp-content/uploads/ICLAC_Guide-to-Mouse-Cell-Line-Authentication_02-Mar-2023.pdf)
- [MMHCdb (Mouse Models of Human Cancer Database)](https://tumor.informatics.jax.org/)
- [MTB (Mouse Tumor Biology) — note: no cell-line data](https://www.jax.org/jax-mice-and-services/preclinical-research-services/mouse-tumor-biology-database)
- [Taconic Syngeneic Cell Line Reference Database](https://www.taconic.com/resources/databases/syngeneic-cell-line-reference-database)
- [MGI Gene Expression Database (GXD)](https://www.informatics.jax.org/expression.shtml)
- [Charles River Cancer Model Database](https://www.criver.com/cancer-model-database)
- [GEMiCCL: genotype + expression of cancer cell lines](https://academic.oup.com/database/article/doi/10.1093/database/bay041/4991663)
- [Single-cell atlas of the TIME across syngeneic murine models — PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12660277/)
- [A disease model resource reveals core principles of tissue-specific cancer evolution — _Nature_](https://www.nature.com/articles/s41586-026-10187-2)
