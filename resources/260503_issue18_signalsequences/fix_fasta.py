from pathlib import Path

input_path = Path(r"C:\Users\jiang\OneDrive\Desktop\PETAdex\sequences.fasta")
output_path = Path(r"C:\Users\jiang\OneDrive\Desktop\PETAdex\sequences_fixed.fasta")

text = input_path.read_text(encoding="utf-8")

# Convert literal "\n" text into actual newlines
fixed = text.replace("\\n", "\n")

# Optional cleanup: normalize line endings and remove extra blank lines
lines = [line.rstrip() for line in fixed.splitlines()]
cleaned = "\n".join(line for line in lines if line.strip()) + "\n"

output_path.write_text(cleaned, encoding="utf-8")

print(f"Fixed FASTA written to: {output_path}")