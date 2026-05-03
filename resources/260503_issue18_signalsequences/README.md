# Issue 18: Signal Sequences

This folder contains the resources and results for issue 18: predicting signal sequences. Is a work in progress(currently contains results of small-scale analysis, not yet database wide)

- fix_fasta.py: when downloaded, the fasta sequences had a formatting issue (specifically line breaks appeared literally as \n), so this script was used to fix the formatting
- sequences.fasta: the raw fastaa sequences (602 sequences total)
- sequences_fixed.fasta: the fastaa sequences properly formatted using fix_fasta.py
- other_summary.csv: the signal sequences prediction results of all the sequences acquired from using the "other" option in SignalP6.0
- eukarya_summary.csv: the signal sequences prediction results of all the sequences acquired from using the "eukarya" option in SignalP6.0
- results_69EFE3B10031...: the localization prediction of the first 301 sequences from sequences_fixed.fasta
- results_69EFF9660038...: the localization prediction of the latter 301 sequences from sequences_fixed.fasta
