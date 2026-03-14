# InterProScan Local Installation on Amazon Linux 2

This document provides a full walkthrough for running InterProScan locally on Amazon Linux 2 for protein FASTA input at moderate scale.

---

## Prerequisites

InterProScan requires Java 11+ and Python 3. Check what you have:

```bash
java -version
python3 --version
```

If Java is missing:

```bash
sudo amazon-linux-extras enable corretto11
sudo yum install -y java-11-amazon-corretto
```

---

## Download and Install InterProScan

InterProScan provides a self-contained tarball with all databases bundled. This is the recommended approach — do not use conda, as the packaged version is often outdated.

```bash
mkdir -p ~/interproscan && cd ~/interproscan

# Download the latest release (check https://github.com/ebi-pf-team/interproscan/releases for current version)
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.72-103.0/interproscan-5.72-103.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.72-103.0/interproscan-5.72-103.0-64-bit.tar.gz.md5

# Verify integrity
md5sum -c interproscan-5.72-103.0-64-bit.tar.gz.md5

# Extract
tar -xzf interproscan-5.72-103.0-64-bit.tar.gz
cd interproscan-5.72-103.0
```

The tarball is large (~20 GB with all member databases), so make sure your instance has sufficient disk. A 50–100 GB EBS volume is reasonable for moderate-scale work.

---

## Initial Setup

InterProScan needs to index its bundled lookup tables before first use:

```bash
python3 setup.py -f interproscan.properties
```

This step decompresses and prepares the H2 database files used by some analyses. It can take 10–20 minutes.

---

## Test the Installation

Run against the bundled test sequences first:

```bash
./interproscan.sh -i test_all_appl.fasta -f tsv -dp
```

The `-dp` flag disables the pre-calculated match lookup (forcing local computation), which is useful to confirm your local install works correctly. Expect this to take a few minutes and produce output files in the current directory.

---

## Running on Your Own Protein FASTA

A typical command for moderate-scale protein analysis:

```bash
./interproscan.sh \
  -i /path/to/your/proteins.fasta \
  -o /path/to/output/results \
  -f TSV,GFF3 \
  -goterms \
  -pa \
  -cpu 8 \
  -dp
```

Key flags:

| Flag | Meaning |
|---|---|
| `-i` | Input FASTA (protein sequences) |
| `-o` | Output file prefix |
| `-f` | Output formats; TSV is easiest to parse downstream |
| `-goterms` | Map matched domains to GO terms |
| `-pa` | Include pathway annotations (MetaCyc, KEGG, Reactome) |
| `-cpu` | Number of threads; set to your instance's vCPU count |
| `-dp` | Disable pre-computed lookup; forces local computation |
| `-appl` | (optional) restrict to specific databases — see below |

---

## Restricting to Specific Analyses

Running all member databases is slow. For PETase-type work you likely care most about:

```bash
-appl Pfam,TIGRFAM,CDD,SUPERFAMILY,Gene3D,PANTHER
```

This cuts runtime significantly. You can list all available applications with:

```bash
./interproscan.sh --list-analyses
```

---

## Practical Considerations at Thousands-of-Sequences Scale

**Chunking input.** InterProScan can struggle with very large single FASTA files. Split your input into chunks of ~5,000 sequences each and run in parallel or sequentially:

```bash
# Split with seqkit (install via conda or download binary)
seqkit split2 -s 5000 proteins.fasta -o chunks/
```

Then loop over chunks:

```bash
for f in chunks/*.fasta; do
  base=$(basename $f .fasta)
  ./interproscan.sh -i $f -o results/${base} -f TSV -goterms -cpu 8 -dp
 done
```

**Memory.** The default JVM heap is 2 GB. If you have a memory-rich instance, increase it in `interproscan.properties`:

```
jvm.maximum.memory.size=8096m
```

**Nohup / tmux.** These runs will take hours. Run inside a tmux session or with nohup so SSH disconnection doesn't kill the job:

```bash
tmux new -s iprscan
# then run your command inside
```

**Output.** TSV output has one domain hit per line with columns: protein ID, MD5, length, analysis, accession, description, start, stop, e-value, status, date, IPR accession, IPR description, GO terms. This is straightforward to parse with pandas or awk for downstream filtering.

---

## Common Issues on Amazon Linux 2

- **libz / libstdc++ errors:** Install `sudo yum install -y libstdc++ zlib` if you see linker errors from embedded binaries (hmmer, signalp, etc.)
- **Perl missing:** Some member databases use Perl scripts; install with `sudo yum install -y perl`
- **Disk space:** The full installation with data is ~25 GB; monitor with `df -h` before starting
