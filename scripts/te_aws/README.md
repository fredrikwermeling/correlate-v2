# Retroelement quantification across CCLE — cloud run

Counts transcription of LINE-1, HERV-K and SVA elements in every CRISPR-cohort
cell line, straight off the public DepMap CCLE BAMs, and produces one small
table for the Cell Line Browser.

## Why this runs in AWS and not on a laptop

The BAMs are public (`s3://depmap-omics-ccle`, no credentials), but DepMap
published an index for only 58 of 1,032 RNA files. Without an index there is no
random access, so each BAM must be read once end to end: 669 files, 9.7 TB.
Over a home connection that is weeks. Inside AWS **us-east-1**, where the bucket
lives, it is a few hours and the transfer is free.

Nothing is written to disk — each BAM is streamed and discarded.

## Steps

1. **Account.** Ask KI IT whether the university already has an AWS agreement;
   institutional billing avoids a personal card. Otherwise any AWS account works.

2. **Launch an instance** — this is the only setting that matters:
   - Region: **us-east-1 (N. Virginia)**. Same region as the data, or the run is
     slow and no longer free.
   - Type: `m6i.4xlarge` (16 vCPU, 64 GB) is ample. `m6i.2xlarge` also fine.
   - Image: Ubuntu 22.04 or 24.04.
   - Disk: default 30 GB. Nothing large is stored.

3. **Copy this folder to the instance** and connect to it:
   ```
   scp -i YOUR_KEY.pem -r te_aws ubuntu@INSTANCE:~/
   ssh -i YOUR_KEY.pem ubuntu@INSTANCE
   ```

4. **Verify first — do not skip this.** Runs 5 cell lines already measured on
   the laptop and compares per-locus counts:
   ```
   cd te_aws && bash setup_and_run.sh verify
   ```
   Expect `VERIFIED - streaming reproduces the pilot`. If it reports MISMATCH,
   stop: the two access paths disagree and the full run would be meaningless.

5. **Full run**, inside tmux so it survives a dropped connection:
   ```
   tmux new -s te
   bash setup_and_run.sh full
   #  detach: Ctrl-b then d      reattach: tmux attach -t te
   ```
   669 cell lines. Resumable — rerunning skips whatever is already in
   `counts_all.json`, so an interruption costs one file, not the job.

6. **Bring back one file** and shut the instance down:
   ```
   scp -i YOUR_KEY.pem ubuntu@INSTANCE:~/te_aws/counts_all.json .
   ```
   **Terminate the instance in the AWS console** — it bills per hour until you do.

## Cost and time

| | |
|---|---|
| Data streamed | 9.7 TB (free, same region) |
| Instance | ~$0.77/hour for m6i.4xlarge |
| Expected wall time | ~6–12 h at 12 workers |
| **Expected total** | **roughly $10–40** |

Set a billing alert. The only way this gets expensive is leaving the instance
running after the job finishes.

## What comes back

`counts_all.json` — per cell line: total mapped reads (library size) and, for
each of 780 panel loci, all-read and unique-read (MAPQ>=20) counts. A few MB.
That is the input for the Cell Line Browser layer.

## The panel

`te_panel.json`, 780 loci, hg19, built from UCSC RepeatMasker:

- 237 LINE-1 — full-length (>=5.8 kb) intergenic L1HS plus the longest L1PA2
- 350 HERV-K — full-length HERVK-int proviruses and LTR5_Hs
- 163 SVA — SVA_E and SVA_F
- 30 controls — 25 L1HS copies *inside* genes (passenger transcription) and
  5 housekeeping genes (library health)

Elements overlapping any RefSeq gene (±2 kb) are excluded, because a repeat
inside an expressed gene picks up reads from the host transcript rather than
from its own promoter. The controls exist to prove that separation held: in the
49-line pilot, intergenic element signal correlated only r=+0.27 with genic-L1
passenger signal and r=+0.17 with housekeeping expression.

## Files

| File | Purpose |
|---|---|
| `setup_and_run.sh` | installs deps, runs verify or full |
| `te_quant_stream.py` | the counter; one streaming pass per BAM, no index |
| `verify_against_pilot.py` | streaming vs laptop pilot, per locus |
| `te_analyze.py` | library health, passenger control, dynamic range |
| `te_panel.json` | the 780 loci |
| `lines_verify.json` / `lines_all.json` | 5-line check / 669-line target |
| `pilot_reference_counts.json` | laptop pilot results, the reference |
