# Protein surface-residue calculations

Completed 2026-09-11 for 289 of 295 cohort sequences (97.97% coverage).
The accepted structures cover 114,170 residues and 5,675 of the 5,798 measurements.
The subsequent predictive comparison is complete; see
`../surface_benchmark/REPORT.md`. It found no reliable surface advantage on the
primary matched-solvent ranking outcome.

## Outputs

- `residues/seqXXXX.csv`: one row per sequence residue with amino acid, 1-based
  sequence position, absolute accessible area (square angstroms), relative
  accessible area, polar/apolar areas, CA confidence and exposure flags.
- `features.csv`: 27 whole-structure surface descriptors and 27 corresponding
  confidence-restricted descriptors, keyed by sequence ID.
- `manifest.json`: one entry for every requested protein, including sequence and
  structure hashes, coverage status, structure source, confidence and errors.
- `summary.json`: cohort, method and software metadata.

## Definitions and checks

FreeSASA 2.2.1 calculates Lee-Richards accessibility with a 1.4 angstrom probe and
its default classifier/reference areas. RSA is relative accessible area; values
above one are retained. Exposure flags use RSA at least 0.10, 0.20 and 0.30. The
primary descriptor threshold is 0.20.

Every accepted structure chain matches its complete input sequence exactly.
Raw input sequences and downloaded structures were not edited. All accepted
structures are single-chain models. Calculations use standard protein heavy
atoms in the first provided model, excluding ligands and waters. These are not
verified biological assemblies. Surface descriptors therefore describe the
provided predicted conformation, not necessarily the enzyme's solution assembly.

Main descriptors include exposed-residue fraction, mean RSA, accessible area per
residue and frequencies of the 20 amino acids among exposed residues. Area-weighted
acidic (DE), basic (KR), hydrophobic (AVILMFWY) and uncharged polar (STNQ) fractions
sum the accessible areas of those residue groups across the whole chain, divided
by total chain accessible area. These groups do not exhaust all amino acids.

CA B-factors are interpreted as pLDDT for the database-provided predicted models.
Confidence-restricted summaries use residues with pLDDT at least 70, while keeping
all residues in the physical accessibility calculation. The median protein mean
pLDDT is 93.19. Three proteins have no residues reaching 70:

| Sequence | Mean pLDDT | Maximum pLDDT |
| --- | ---: | ---: |
| seq0041 | 33.38 | 64.88 |
| seq0099 | 24.39 | 37.00 |
| seq0170 | 30.01 | 38.88 |

Their confidence-restricted descriptors are intentionally missing, not zero.
Thus 286 proteins have confidence-restricted surface summaries. Their low-confidence
counterparts remain explicitly flagged in the manifest. Confidence filtering and
exposure-threshold sensitivity must be considered before interpreting prediction.

Validation checked all accepted residue sequences, consecutive positions,
nonnegative finite areas, confidence bounds, exposure flags and normalized amino
acid frequencies. An independent BioPython Shrake-Rupley calculation on cached
protein 17 gave total area 8,043.62 square angstroms versus 7,973.54 from FreeSASA,
a difference of 0.87%. This is a one-protein sanity check using different default
atomic radii, not a proof that every predicted structure is accurate.

## Six unresolved proteins

| Sequence | Database protein | Reason |
| --- | ---: | --- |
| seq0007 | 171 | Structure has 253 extra N-terminal residues; all 390 target residues otherwise match |
| seq0014 | 219 | Structure has 28 extra N-terminal residues; all 388 target residues otherwise match |
| seq0032 | 41 | Structure has one extra N-terminal residue; all 316 target residues otherwise match |
| seq0147 | 154 | Published page has an empty structure-file attribute |
| seq0209 | 155 | Published page has an empty structure-file attribute |
| seq0250 | 86 | Equal length, but one amino-acid mismatch among 568 residues |

These were not silently trimmed or substituted. Extra modeled segments can
change calculated exposure; a differing residue is not the exact target protein.

## Reproduction and next comparison

```bash
PYTHONPATH=src python3 -m enzsolv.run_surface --workers 2
```

Completed structures and catalogue pages are cached. Network requests retry
transient failures up to three attempts. The first fetch was interrupted by a
timeout; the resumed run completed all candidate requests and calculations.
Twenty-four tests pass.

Compare conditions-only, sequence and surface models on the same structure-covered
rows and saved grouping assignments. The full-cohort baseline alone is not a fair
comparator for the smaller structure-covered subset. Keep quality flags out of
unexamined biological interpretations and report the three low-confidence cases.

Method reference: https://freesasa.github.io/python/classes.html
Structures: https://organizymedb.org/
