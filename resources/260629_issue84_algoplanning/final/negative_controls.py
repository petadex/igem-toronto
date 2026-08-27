"""Inject deliberate faults and confirm the harness catches each one."""
import json, os, shutil, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable                  # whatever interpreter runs this
VALIDATOR = os.path.join(HERE, "validate_design.py")
if len(sys.argv) < 2:
    sys.exit(f"usage: python {os.path.basename(sys.argv[0])} <run_dir>")
SRC = sys.argv[1]

def mutate(tag, fn):
    tmp = tempfile.mkdtemp(prefix="neg_")
    dst = os.path.join(tmp, "run")
    shutil.copytree(SRC, dst)
    S = json.load(open(os.path.join(dst, "summary.json")))
    fn(S)
    json.dump(S, open(os.path.join(dst, "summary.json"), "w"), indent=2)
    p = subprocess.run([PY, VALIDATOR, dst],
                       capture_output=True, text=True)
    fails = [l.strip() for l in p.stdout.splitlines() if "[FAIL]" in l]
    print(f"--- {tag}: exit {p.returncode}")
    for f in fails[:3]:
        print(f"      {f}")
    if not fails:
        print("      *** NOT CAUGHT ***")
    shutil.rmtree(tmp, ignore_errors=True)

def site(S):                       # plant a BsmBI site
    o = S["fragments"][0][0]["oligo"]
    S["fragments"][0][0]["oligo"] = o[:9] + "CGTCTC" + o[15:]

def uncover(S):                    # break a codon so cores stop being encoded
    o = S["fragments"][0][0]["oligo"]
    S["fragments"][0][0]["oligo"] = "TGG" + o[3:]

def junction(S):                   # break overhang uniformity
    u = S["fragments"][0][1]
    S["fragments"][0][1]["oligo"] = u["oligo"][:-3] + "AAA"

def stop(S):                       # plant a stop codon
    o = S["fragments"][0][0]["oligo"]
    S["fragments"][0][0]["oligo"] = o[:6] + "TAA" + o[9:]

def frame(S):                      # break the reading frame
    S["fragments"][0][0]["oligo"] = S["fragments"][0][0]["oligo"][:-1]

def cap(S):                        # violate the library cap
    S["args"]["max_library"] = 10

mutate("forbidden site planted", site)
mutate("codon changed (coverage)", uncover)
mutate("junction codon changed", junction)
mutate("stop codon planted", stop)
mutate("reading frame broken", frame)
mutate("library cap violated", cap)
