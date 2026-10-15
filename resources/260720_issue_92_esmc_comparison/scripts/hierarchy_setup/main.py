import subprocess
import sys

# pull the huge orfs file from S3
subprocess.run([sys.executable, "./get_orfs.sh"])

# start by running format_data.py, which fetches data from S3 and catalytic orfs file, and saves it locally
subprocess.run([sys.executable, "./format_data.py"])

# run build_hierarchy.py, which builds the hierarchy of 30% superfamilies and 60% families, and saves it locally
subprocess.run([sys.executable, "./build_hierarchy.py"])