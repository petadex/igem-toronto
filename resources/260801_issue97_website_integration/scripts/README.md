# Scripts

Small helpers for website integration checks. They do not replace petadex.io.

## probe_structure_s3.sh

Range-GET probes on Alex’s structure bucket for a known ORF id (default `4981589`).

```bash
chmod +x probe_structure_s3.sh
./probe_structure_s3.sh
# optional:
./probe_structure_s3.sh 4981589
```

Expected (Aug 2026):

| Path | Typical status |
|------|----------------|
| `esmfold2-centroids/test2/...` | 206 (readable) |
| `esmfold2-centroids/60pid/...` | 403 until Dennis opens ACL |
| `esmfold2-centroids/60pid-msa/...` | 403 until Dennis opens ACL |

206 / 200 = object is publicly readable. 403 = ACL still closed. 404 = key missing (or wrong ORF id).
