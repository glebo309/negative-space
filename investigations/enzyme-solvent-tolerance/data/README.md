# Data

`raw/measurements.csv.gz` contains the frozen OrganiZymeDB input used for this analysis, downloaded on 2026-09-11.

- [OrganiZymeDB dataset and column reference](https://organizymedb.org/download)
- [Current complete CSV](https://organizymedb.org/download/measurements.csv)
- [OrganiZymeDB v1.0 record and DOI](https://doi.org/10.5281/zenodo.21874133)

The live database may change as measurements are added or corrected. The compressed file preserves the exact input behind the reported results. Extract it before running the analysis:

```sh
gzip -dk data/raw/measurements.csv.gz
shasum -a 256 data/raw/measurements.csv
```

SHA-256:

`400823131c996a3899d95cdb0c82e37a6aa5833d32f0daa23e15a696045b3a06`

OrganiZymeDB identifies Romain Debruyne and Fabrizio Pucci at Université Libre de Bruxelles as its creators and licenses the data under [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/). Commercial use requires contacting the database authors.

The snapshot contains public scientific records and protein sequences. It remains third-party data and is not covered by the repository's MIT software license.

Database-provided structures are not included. `python3 -m enzsolv.run_surface` retrieves the referenced structures and verifies full-chain sequence identity before calculating accessibility features.
