#!/usr/bin/env bash
#
# setup_vm.sh -- install everything run_pipeline.sh needs on a fresh Ubuntu VM.
# Run once.  Must be in us-east-1 (RDS + S3 are in region).
#
set -euo pipefail

sudo apt-get update
# hmmer -> hmmsearch ; awscli -> aws ; pv -> extract progress bar (optional)
sudo apt-get install -y hmmer awscli python3-pip pv

# seqkit: static binary, not in apt
if ! command -v seqkit >/dev/null; then
    wget -q https://github.com/shenwei356/seqkit/releases/latest/download/seqkit_linux_amd64.tar.gz
    tar xf seqkit_linux_amd64.tar.gz
    sudo mv seqkit /usr/local/bin/
    rm -f seqkit_linux_amd64.tar.gz
fi

# python deps: psycopg2 (fetch_clusters SQL), openpyxl (oPools xlsx), matplotlib (figures)
pip3 install --quiet psycopg2-binary openpyxl matplotlib

echo "--- verifying ---"
for t in aws seqkit hmmsearch python3; do
    if command -v "$t" >/dev/null; then echo "ok:      $t"; else echo "MISSING: $t"; fi
done
python3 -c "import psycopg2, openpyxl, matplotlib; print('ok:      python deps (psycopg2, openpyxl, matplotlib)')"
echo "setup complete."
