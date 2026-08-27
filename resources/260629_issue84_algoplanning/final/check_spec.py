"""Check a design against the wet-lab spec, clause by clause, in its wording."""
import json, sys, itertools
IUPAC = {"A":"A","C":"C","G":"G","T":"T","R":"AG","Y":"CT","S":"CG","W":"AT",
         "K":"GT","M":"AC","B":"CGT","D":"AGT","H":"ACT","V":"ACG","N":"ACGT"}
BASES="TCAG"; AAS="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
TAB={a+b+c:aa for (a,b,c),aa in zip([(x,y,z) for x in BASES for y in BASES
                                     for z in BASES],AAS)}
STOPS={c for c,a in TAB.items() if a=="*"}
def rc(s): return s.translate(str.maketrans("ACGT","TGCA"))[::-1]
def exp(t): return ["".join(x) for x in itertools.product(*(IUPAC[b] for b in t))]
def cod(s): return [s[i:i+3] for i in range(0,len(s),3)]
def may(seq, sites):
    for st in sites:
        for i in range(len(seq)-len(st)+1):
            if all(st[k] in IUPAC.get(seq[i+k],"") for k in range(len(st))):
                return True
    return False

SITE, RCSITE = "CGTCTC", "GAGACG"
FIRST_J, LAST_J = "CGGA", "GGTG"

for run in sys.argv[1:]:
    S=json.load(open(run+"/summary.json")); F=S["fragments"]; cuts=S.get("cuts",[])
    print(f"\n{'='*72}\n{run.split('/')[-1]}   K={S['recommended_K']}  cuts={cuts}")
    olig=[(f+1,i+1,u["oligo"]) for f,us in enumerate(F) for i,u in enumerate(us)]

    # -- 1. no BsmBI site anywhere inside a fragment ----------------------- #
    bad=[f"frag{f}#{i}" for f,i,o in olig if may(o,(SITE,RCSITE))]
    print(f"\n[1] BsmBI {SITE}/{RCSITE} inside any fragment")
    print(f"    {len(olig)} oligos scanned (all IUPAC expansions): "
          f"{'VIOLATIONS '+str(bad) if bad else 'none -- PASS'}")
    xj=[]
    for f in range(len(F)-1):
        for i,a in enumerate(F[f]):
            for j,b in enumerate(F[f+1]):
                if may(a["oligo"][-6:]+b["oligo"][:6],(SITE,RCSITE)):
                    xj.append(f"frag{f+1}#{i+1}|frag{f+2}#{j+1}")
    print(f"    across every fragment junction: "
          f"{'VIOLATIONS '+str(xj[:3]) if xj else 'none -- PASS'}")

    # -- 2. no in-frame stop codons ---------------------------------------- #
    print(f"\n[2] in-frame stop codons")
    frame=[f"frag{f}#{i} len {len(o)}" for f,i,o in olig if len(o)%3]
    print(f"    reading frame (every oligo a whole number of codons): "
          f"{'BROKEN '+str(frame) if frame else 'PASS'}")
    st=[f"frag{f}#{i} codon {ci} ({t})" for f,i,o in olig
        for ci,t in enumerate(cod(o)) if set(exp(t))&STOPS]
    print(f"    stops within fragments (all expansions): "
          f"{'VIOLATIONS '+str(st[:3]) if st else 'none -- PASS'}")
    # do junctions create a NEW codon that could be a stop?
    offs={len(o)%3 for _f,_i,o in olig}
    print(f"    junctions fall on codon boundaries (so no new codon is formed): "
          f"{'YES -- PASS' if offs=={0} else 'NO -- CHECK'}")

    # -- 3. reserved overhangs ---------------------------------------------- #
    print(f"\n[3] reserved junctions {FIRST_J} (backbone->frag1) / "
          f"{LAST_J} (fragK->backbone)")
    ohs=[]
    for f in range(len(cuts)):
        cl=cod(F[f][0]["oligo"])[-1]; cr=cod(F[f+1][0]["oligo"])[0]
        ohs.append((cuts[f], cl[1:]+cr[:2]))
    if not ohs: print("    no internal junctions (K=1)")
    for p,oh in ohs:
        clash=[r for r in (FIRST_J,LAST_J) if oh==r or rc(oh)==r]
        print(f"    internal overhang at col {p}: {oh} (rc {rc(oh)})  "
              f"{'REUSES '+str(clash) if clash else 'distinct from both -- PASS'}")
    # stronger than the spec: 1-mismatch cross-reactivity with the reserved pair
    def conf(a,b):
        return a==b or sum(1 for x,y in zip(a,rc(b)) if x!=y)<=1
    near=[f"{oh} vs {r}" for _p,oh in ohs for r in (FIRST_J,LAST_J) if conf(oh,r)]
    print(f"    (beyond spec) within 1 mismatch of a reserved overhang: "
          f"{near if near else 'none -- PASS'}")

    # -- 4. universal primer tails ------------------------------------------ #
    print(f"\n[4] BsmBI sites + universal primer tails on each end")
    f1=F[0][0]["oligo"]; fk=F[-1][0]["oligo"]
    print(f"    emitted fragment 1 begins {f1[:12]}... ; fragment {len(F)} "
          f"ends ...{fk[-12:]}")
    print(f"    contains a BsmBI site: {may(f1,(SITE,))} / {may(fk,(SITE,))}")
    print(f"    -> NOT PRESENT: the designer emits CODING SEQUENCE ONLY.")
    print(f"       No BsmBI sites, no spacers, no primer tails, and no {FIRST_J}/"
          f"{LAST_J} backbone overhangs are in these sequences.")
