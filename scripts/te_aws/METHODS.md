# Methods: retroelement transcription across the CCLE panel

Draft methods text and the parameter record behind it. Numbers marked
**[pending]** are filled in when the full run completes.

## Data source

Cell line RNA-seq alignments were obtained from the Cancer Dependency Map
(DepMap) Cancer Cell Line Encyclopedia (CCLE) collection hosted on the AWS
Registry of Open Data (`s3://depmap-omics-ccle`, public, no controlled-access
application required). The collection holds 1,032 RNA-seq alignment files
(975 aligned to hg19, 57 to hg38). Analysis was restricted to hg19 alignments
whose cell line was also present in the CRISPR gene-effect cohort used by
Correlate (DepMap 26Q1; 1,208 lines), giving **669 cell lines**. All
coordinates below are hg19.

## Retroelement panel

Element coordinates were taken from the UCSC RepeatMasker annotation for hg19
(`rmsk`). Three families were included:

| Family | Selection rule | n |
|---|---|---|
| LINE-1 | `L1HS` and the longest `L1PA2` copies, length >= 5,800 bp | 237 |
| HERV-K | `HERVK-int` >= 5,000 bp and `LTR5_Hs` >= 900 bp | 350 |
| SVA | `SVA_E`, `SVA_F`, length >= 1,000 bp | 163 |

The length thresholds select full-length elements: truncated copies, which make
up the large majority of retroelement sequence in the genome, cannot initiate
autonomous transcription and contribute predominantly passenger signal.

Elements overlapping any RefSeq gene, extended by 2 kb on each side, were
excluded. A retroelement inside an expressed gene accumulates reads from the
host transcript regardless of its own promoter activity, and this is the
dominant confounder in repeat quantification.

## Controls

Two control sets were carried through the identical pipeline (30 additional
loci, 780 in total):

1. **Passenger control** (n=25): full-length L1HS copies located *inside*
   RefSeq genes. These are expected to track host-gene transcription. If
   intergenic element signal correlated strongly with this set, the assay would
   be measuring bulk transcription rather than element activity.
2. **Library control** (n=5): exons of *ACTB*, *GAPDH*, *RPL13A*, *TBP* and
   *B2M*, used to confirm comparable library composition across lines.

## Quantification

Each alignment file was streamed once, without local storage, and reads were
assigned to panel intervals by a single coordinate-ordered pass. Reads flagged
unmapped, secondary, supplementary or duplicate were discarded. Two counts were
retained per locus: all remaining reads, and reads with MAPQ >= 20 ("unique").
Unique counts are used throughout; multi-mapping reads are intrinsic to repeat
families and their inclusion inflates counts in proportion to family size rather
than to activity.

Library size was taken as the total number of retained mapped reads in the same
pass, and counts were expressed as counts per million (CPM). An element was
called active at CPM > 0.5.

Analysis used pysam 0.24.0 (htslib) on Python 3.12.

## Validation

Two independent checks were performed before the full run.

**Method agreement.** DepMap publishes BAM indices for 58 of the 1,032 RNA
files. For those, counts were obtained a second way, by indexed random access to
each locus rather than by streaming the whole file, on separate hardware. Across
5 cell lines and 1,643 loci with non-zero counts, 1,620 (98.6%) agreed exactly
and the remaining 23 fell within tolerance (<= 2 reads or <= 5%); none
disagreed beyond it.

Total library size is not comparable between the two paths and was excluded from
the comparison: the indexed route reads its total from the BAM index, which
counts every mapped record including secondary and supplementary alignments,
whereas the streaming pass counts primary non-duplicate reads. RNA-seq
multi-mapping makes this a 31-50% difference by construction. The streaming
definition is used as the CPM denominator throughout.

**Specificity.** In a 49-cell-line pilot, summed intergenic element signal
correlated only weakly with the genic-L1 passenger control (Pearson r = +0.27)
and with housekeeping expression (r = +0.17), indicating that element signal is
not a proxy for transcriptional output in general. Active element counts ranged
from 0 to 47 per line (median 11). In the most active line, the single strongest
locus accounted for 19% of element signal and the top five for 42%, so the
signal is distributed across elements rather than driven by one outlier locus.

## Limitations

- Coverage is 669 of 1,208 CRISPR-cohort lines (55%); the remainder have no
  hg19 CCLE RNA-seq alignment in the public collection.
- The assay measures element **transcription**. Retrotransposition (new somatic
  insertions) cannot be inferred from RNA and would require whole-genome
  sequencing.
- Alu was not included. At ~1.1 million copies, locus-level quantification is
  not tractable by this approach, and family-level counting is dominated by
  intronic passenger signal. The relevant Alu ligand for MDA5 is inverted-repeat
  pairs in 3'UTRs, which would require a separately curated interval set.
- Most libraries are unstranded (984 of 1,032), so antisense element
  transcription cannot be separated from sense transcription.

## Reproducibility

Panel definition, quantification and validation code:
`scripts/te_aws/` in the Correlate repository. The panel is distributed as
`te_panel.json`; the quantifier (`te_quant_stream.py`) takes that file and a
list of alignment URLs and is deterministic given both.
