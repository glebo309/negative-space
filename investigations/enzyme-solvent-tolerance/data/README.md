# Data

`raw/measurements.csv.gz` contains a compressed frozen snapshot from [OrganiZymeDB](https://organizymedb.org/), downloaded on 2026-09-11. Extract it before running the analysis:

```sh
gzip -dk data/raw/measurements.csv.gz
shasum -a 256 data/raw/measurements.csv
```

SHA-256:

`400823131c996a3899d95cdb0c82e37a6aa5833d32f0daa23e15a696045b3a06`

OrganiZymeDB identifies Romain Debruyne and Fabrizio Pucci at Université Libre de Bruxelles as its creators and licenses the data under [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/). Commercial use requires contacting the database authors.

The analysis preserves this exact snapshot because the database may continue to change. The snapshot contains public scientific records and protein sequences. It remains third-party data and is not covered by the repository's MIT software license.

Database-provided structures are not included. `python3 -m enzsolv.run_surface` retrieves the referenced structures and verifies full-chain sequence identity before calculating accessibility features.
