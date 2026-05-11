# Cell Line Browser Redesign — Design Notes

## Goal
Reduce the cell line detail card from ~15 stacked sections to 3 compact blocks. 
Move all deep-dive content into a new "Wiki" modal that opens from the card.

## Card structure (3 blocks)
1. Identity — name, ACH ID, sex symbol, tissue · subtype
2. Specific features — 5–8 lines, each a characterizing alteration with its 
   functional consequence (e.g. "BRAF p.V600E → MAPK pathway active")
3. Wiki button — full-width, opens the wiki modal

## Feature-line rules (deterministic, no LLM)
A "specific feature" line is generated for each of:
- Hotspot driver mutations (BRAF V600E, KRAS G12*, etc.)
- Clinically relevant fusions
- LOF of tumor suppressors in functional-loss list (PTEN, TP53, RB1, APC, etc.)
- Focal amplifications and deep deletions
- MSI-high
- TMB-high (>10 mut/Mb)
- MHC-I loss (HLA-A/B/C or B2M alterations)
- WGD and/or high CIN (single combined line)

Each line: glyph + alteration + arrow + consequence.
Glyphs: ● filled = GOF/driver, ○ open = LOF, ◆/◇ = CN event.

## Wiki modal
Opens as a full-screen modal overlay. Sticky left-side TOC. Sections:
1. Overview (auto-generated paragraph + key facts)
2. Identity & origin
3. Genome (ploidy, WGD, aneuploidy, CIN)
4. Pathway status (detailed view of what the card summarized)
5. Mutations (Hotspots / Damaging / All — tabbed)
6. Copy number (Focal / Arm-level)
7. Fusions (Clinically relevant / All raw)
8. Gene dependencies (Most depleted / Enriched / Uniquely depleted / Uniquely 
   enriched vs filtered — each with Enrichr link)

## Implementation order
Prompt 1: Skeleton + feature-line data model
Prompt 2: Identity block
Prompt 3: Specific features block
Prompt 4: Wiki modal scaffold + TOC
Prompt 5: Wiki content sections (migrate existing card sections)
Prompt 6: Polish, accessibility, mobile
