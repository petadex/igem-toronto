# igem-toronto

## Documentation

Analysis is documented in Jupyter notebooks under `notebooks/`. Each notebook corresponds to a GitHub issue and follows the naming convention:

```
notebooks/YYMMDD_issue<N>_<topic>.ipynb
files/YYMMDD_issue<N>_<topic>/   ← working files, scripts, figures
```

### Notebook structure

Each notebook contains:
- **Header** — lead, issue link, start/complete dates, file paths, S3 archive location
- **Introduction** — biological or technical motivation
- **Objectives & Methods** — what was done and how
- **Results & Discussion** — findings, coverage stats, follow-up questions
- **Analysis cells** — code that reads data, generates figures, and saves them to `files/<notebook>/`

### Environment setup

A setup cell at the bottom of each notebook creates a local `.venv` and registers an `igem-toronto` Jupyter kernel. Run it once, switch to that kernel, then re-run the notebook. Dependencies are listed in `requirements.txt`.

### Data

Large input files are not committed. They live in `files/<notebook>/` locally and are archived to S3 (`s3://petadex/`). Each notebook documents the exact `aws s3 cp` commands used to fetch and upload data, with access dates.