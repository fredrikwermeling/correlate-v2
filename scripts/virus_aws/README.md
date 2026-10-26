# Viral reads across CCLE: cloud (or laptop) run

Looks for virus RNA in every CRISPR-cohort cell line, using the reads in the
public DepMap CCLE RNA-seq BAMs that did not match the human genome. It uses the
same 669 hg19 BAMs and the same `s3://depmap-omics-ccle` bucket as `../te_aws`.

## What was checked first (2026-10-26)

- **The unmapped reads are in the files.** The BAMs are Broad's 2013 TopHat
  1.4.1 alignments plus Picard MarkDuplicates. They are coordinate-sorted and end
  in a block of unmapped reads (`check_one.py`). In RT4 the block is 1.47 GB,
  9% of the file. In Raji, MJ and A549 it is 10 to 15%.
- **It finds known viruses.** The last 25 MB of five files were tested against
  the Cellosaurus records:

  | Cell line | Recorded | Reads in ~330k unmapped | Genome breadth |
  |---|---|---|---|
  | Raji | EBV | EBV 309 | 0.52 |
  | NAMALWA | EBV | EBV 138 | 0.28 |
  | PLC/PRF/5 | HBV | HBV 848 | 0.44 |
  | MJ | HTLV-1 | HTLV-1 37,079 | 1.00 |
  | RT4 | none | none of these; **XMRV 45, breadth 0.5** | |

  The mouse retrovirus XMRV checks out against the literature too:
  - 22Rv1 (the known XMRV producer): 112,926 reads, 35% of its unmapped reads,
    breadth 0.94.
  - VCaP: 97,316 reads, breadth 0.92.
  - PC-3, DU 145 and LNCaP: 0 to 6 reads (negative).

  RT4's 45 XMRV reads are therefore far below a real infection. Low-level hits
  like these are more likely carry-over from another library sequenced in the
  same run than an infection of the cell line.

## How it works

`virus_scan.py` finds where the unmapped block starts with a binary search of
HTTP range reads, then downloads only that block. Each read is cut into 25-mers
(every 4th position) and looked up against every 25-mer of the reference genomes
(`references.json`, both strands). The panel has 18 references: EBV, KSHV,
HSV-1, CMV, HBV, HCV, HTLV-1, HIV-1, HPV16/18/45, MCPyV, BKPyV, SV40,
adenovirus 5, XMRV, MoMLV and Mycoplasma hyorhinis.

- A read counts for a virus when 2 or more of its 25-mers belong to that
  reference alone. 25-mers shared between two references (such as XMRV and
  MoMLV) are dropped, and so are low-complexity ones.
- **Breadth** is the share of 50 equal bins of the genome that the reads touched.
  An infection covers the genome. A vector fragment (SV40 promoter, MLV LTRs) or a
  shared conserved stretch (the one-bin mycoplasma hits seen in every line) does
  not. `analyze.py` calls a virus only when both of these hold:
  - at least 10 reads per million unmapped reads;
  - breadth of 0.2 or more.
- **Not scanned:** unmapped reads whose mate mapped. They sit inside the mapped
  part of the file. Some integration-junction reads are among them.
- **Mycoplasma:** the libraries are polyA-selected and bacterial mRNA has no
  polyA tail. So a mycoplasma hit is informative, but no hit is weak evidence.

## Two ways to run

**Quick screen from a laptop: about 1 hour, 17 GB.** It reads the last 25 MB of
each file (about 330,000 unmapped reads per cell line). That is enough for
strongly expressed viruses and contaminants. A weak signal can be missed.

```
python3 fetch_refs.py
python3 virus_scan.py lines_verify.json counts_verify.json --sample-mb 25 --workers 5
python3 verify_controls.py counts_verify.json          # expect VERIFIED
python3 virus_scan.py lines_all.json counts_sample.json --sample-mb 25 --workers 5
python3 analyze.py counts_sample.json
```

**Full run in AWS: about 1.2 TB of unmapped blocks.** It reads every unmapped
read. Set up as in `../te_aws/README.md`:

- region **us-east-1**, where the bucket is, so the transfer is free and fast;
- `m6i.4xlarge`, Ubuntu, default disk.

```
scp -i KEY.pem -r virus_aws ubuntu@INSTANCE:~/
ssh -i KEY.pem ubuntu@INSTANCE
cd virus_aws && tmux new -s v
bash setup_and_run.sh auto     # verify the five controls, then full only if they pass
# Ctrl-b d to detach; bash status.sh to check; tmux attach -t v to return
scp -i KEY.pem ubuntu@INSTANCE:~/virus_aws/counts_all.json .
```

**Terminate the instance afterwards.** Expect a few hours and roughly $5 to 15.
The run can be resumed: rerunning skips finished cell lines and retries failed ones.

## Files

| File | What it is |
|---|---|
| `check_one.py` | Step one: the aligner, and whether and how big the unmapped block is |
| `fetch_refs.py`, `references.json` | Reference genomes from NCBI into `refs/` |
| `virus_scan.py` | The scan (numpy only, no pysam) |
| `verify_controls.py` | Raji, NAMALWA, PLC/PRF/5, MJ positive; RT4 negative |
| `analyze.py` | Calls per virus, compared with Cellosaurus (`web_data/virus_status.json`) |
| `lines_all.json`, `lines_verify.json` | The 669 hg19 BAMs (from `../te_aws`) and the five controls |
| `lines_extra.json` | 50 more CRISPR-cohort cell lines from the bucket's `metadata/inventory.csv`: 41 hg38 STAR alignments and 9 hg19 files `../te_aws` missed. The unmapped block works the same way whatever the genome build. Together 719 of the 1,208 cohort cell lines; the rest have no public RNA file. |
