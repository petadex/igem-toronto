# Angela website wiring (issue #97)

Sara wired Angela’s SignalP6 table onto petadex.io corpus ORF pages, plus a PID cluster nav prototype.

## Code

- Branch: `feat/angela-sequence-annotations` on https://github.com/sarapr06/petadex.io
- Docs: `frontend/docs/angela-sequence-annotations.md`, `frontend/docs/angela-followups.md`

## Status

| Piece | Status |
|-------|--------|
| SignalP6 (`signalp6_orf_predictions`) | Live on `/sequence/orf/:id` |
| DeepLoc / biochem | Soft stubs; waiting on DB update + Angela load |
| PID parent path + pin centroids | Live on `/cluster/:level/:id` |
| Child enumerate 30→60→90 | Deferred (needs clustering index or map) |
| Ancestral AA | Waiting on Thomas running Angela’s scripts |

## Verify

- `/sequence/orf/294247546` (SignalP SP hit)
- `/cluster/90/1124517` (parent path + pin tray)
