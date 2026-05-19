import pandas as pd

# run every script in /scripts
scripts = [
    "generate_rand_fragments.py",
    "generate_rand_matched.py",
    "generate_rand_empirical.py",
    "generate_rand_95th.py",
    "generate_shuffle.py",
    "generate_rand_uniprot.py",
]

for script in scripts:
    print(f"Running {script}...")
    exec(open(f"./scripts/{script}").read())