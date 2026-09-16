# Can protein sequence and surface composition predict enzyme solvent tolerance?

| Field | Value |
| --- | --- |
| Stage | completed |
| Outcome | practical predictor claim not supported; underlying biological hypothesis unresolved |
| Started | 2026-09-11 |
| Completed | 2026-09-12 |
| Source dataset | [OrganiZymeDB download page](https://organizymedb.org/download) |
| Live dataset | [complete measurements CSV](https://organizymedb.org/download/measurements.csv) |
| Dataset record | [OrganiZymeDB v1.0](https://doi.org/10.5281/zenodo.21874133) |
| Frozen input | [`data/raw/measurements.csv.gz`](data/raw/measurements.csv.gz), downloaded 2026-09-11 |
| Raw SHA-256 | `400823131c996a3899d95cdb0c82e37a6aa5833d32f0daa23e15a696045b3a06` |

## Short result

Across the unfamiliar-family holdouts, protein sequence representations and global surface descriptors did not reliably improve organic-solvent ranking over models based on the recorded experimental conditions. The planned computational gate for prospective bench validation was not met, so I stopped development of the general predictor.

## Question

Given an enzyme sequence, predicted structure and proposed organic cosolvent, can a model rank which solvents and concentration windows the enzyme is most likely to tolerate?

The intended use was solvent selection for a previously unseen enzyme. Absolute retained activities measured in different laboratories were not treated as directly comparable targets.

## Brief context

Organic-solvent tolerance is influenced by protein sequence, surface chemistry, hydration, charge, packing, flexibility and solvent access. Published measurements, however, combine different enzymes, solvents, substrates, assays, incubation conditions and laboratories. A model can appear successful by learning concentration or publication patterns without learning transferable protein-specific solvent tolerance.

The investigation therefore tested whether explicit protein information improves a conditions-only comparator under grouped holdouts and a ranking metric aligned with solvent selection.

## Original biological hypothesis

Enzymes that remain stable in organic solvents contain transferable sequence and surface-composition features that distinguish them from solvent-sensitive enzymes. In particular, the identity and distribution of solvent-exposed residues may help an enzyme retain its structure and activity during solvent exposure. By comparing tolerant and sensitive enzymes, it should therefore be possible to learn a sequence or surface signature that predicts the solvent stability of an unseen enzyme.

The relevant residue pattern may depend on the solvent, concentration, temperature, exposure time, and enzyme family.

## Operational predictive hypothesis

If the biological signal is sufficiently consistent and represented in the available data, frozen protein-language-model representations and explicit solvent-accessible surface descriptors should improve matched-condition solvent ranking over experimental conditions alone when evaluated on unfamiliar protein groups and separated publications.

The computational experiment tested this predictive version. It did not isolate the underlying biological mechanism.

## Experiment

### Data

- Complete downloaded resource: 28,173 measurements
- Incubation-stability subset: 9,091 measurements
- Conservative analysis cohort: 5,798 measurements from 295 exact sequences
- Exact structure-covered cohort: 5,675 measurements from 289 sequences
- Primary surface-ranking evaluation: 763 matched-solvent series from 270 proteins
- Surface calculation: 114,170 residues and 27 global surface descriptors

Only 15 sequences appeared in multiple publications before structural filtering. Curator comments flagged assay or control uncertainty in 3,470 of the 5,798 conservative-cohort rows. No cross-paper repeats matched all required complete numerical conditions.

### Models

- Experimental conditions and solvent identity
- Whole-sequence amino-acid composition
- Frozen ESM2 35M and ESM2 650M representations
- Twenty-seven explicit surface descriptors
- ESM650 plus surface descriptors
- Extra Trees as the primary learner
- A supplementary fixed-alpha interaction ridge model

### Primary evaluation

The primary metric was mean within-protein Spearman correlation for ranking solvents measured under matching recorded conditions. At least three solvents and varying outcomes were required. Models were compared on identical rows using five held-out folds. Sequence identity groups were separated, and the strictest design also separated publications.

Paired 95 percent intervals used 5,000 component-bootstrap draws conditional on the fitted models and stated grouping. They were not treated as equivalence tests or independent biological replications.

### Controls

- Training-only preprocessing and imputation
- Group and publication leakage checks
- Identical evaluation rows for paired models
- Multiple surface-exposure thresholds
- Confidence-restricted surface features
- Additional learner seeds
- Supplementary predictor family
- Planted-signal, null-feature and shuffled-feature controls
- Independent reproduction of decisive calculations

## Results

### Surface benchmark

| Holdout | Conditions | Composition | Surface | Surface minus conditions, paired 95 percent interval |
| --- | ---: | ---: | ---: | --- |
| 30 percent identity plus publication | 0.181 | 0.214 | 0.195 | +0.014, -0.030 to +0.063 |
| 30 percent identity | 0.212 | 0.229 | 0.205 | -0.007, -0.054 to +0.040 |
| 40 percent identity | 0.215 | 0.226 | 0.206 | -0.009, -0.054 to +0.036 |
| Exact sequence | 0.284 | 0.267 | 0.244 | -0.041, -0.083 to +0.001 |

All 18 surface-versus-conditions intervals crossed zero across the main analysis, exposure and confidence sensitivities, and additional seeds.

### Sequence and surface combination

In the strictest publication-linked holdout, ESM650 scored 0.199 and ESM650 plus surface scored 0.204 against 0.181 for conditions. Adding surface to ESM650 changed the score by +0.006 with a paired interval from -0.012 to +0.022. Neither ESM650 nor the combined model reliably exceeded conditions in all four primary holdouts.

### Final decision diagnostics

After explicit metadata controls, the full unfamiliar-family plus publication comparison gave 0.2199 for ESM650 with metadata and 0.1910 for metadata alone. The difference was +0.0289 with an interval from -0.0150 to +0.0756.

A narrower familiar-family composition analysis retained a modest hint after removing a harmful family-label feature: 0.2263 versus 0.2021, difference +0.0242 with an interval from -0.0089 to +0.0824. Best-solvent selection was 26.60 percent versus 25.02 percent, a difference of 1.57 percentage points with an interval from -3.50 to +6.23. This did not establish practical solvent-selection value.

Strong planted protein-by-solvent signals were recovered at correlations around 0.79 under both strict fold designs. Null and shuffled-feature controls did not show a clear advantage. The pipeline could recover a simple strong signal; its sensitivity to weaker or more local biological relationships remains unknown.

The final computational check passed 60 tests in the public package. These tests describe software and artifact coverage, not scientific replication.

## Conclusion

The heterogeneous literature measurements, available metadata, frozen sequence representations, and global surface descriptors did not establish a useful general solvent-selection model. I therefore stopped predictor development and did not launch the proposed robotic validation panel.

The original biological hypothesis remains unresolved. The tested representations did not extract a transferable signal strong enough for solvent selection from this dataset.

## Limits of the conclusion

The experiment did not test individual surface residues or local structural mechanisms directly. It also did not establish equivalence between the protein-aware and conditions-only models. A modest signal remained in the familiar-family composition analysis, and a standardized prospective dataset may produce a different result.

The analysis supports a narrow conclusion about the tested representations and literature-derived labels. It does not support a general claim that protein sequence or surface chemistry is unrelated to solvent tolerance.

## Why the investigation stopped

The predefined gate for prospective bench validation was not met. Searching more model architectures against the same labels would not resolve the main limitation, which is the heterogeneity and incomplete context of the measurements. I would reopen the project for a more comparable independent dataset, a defensible reconstruction of missing experimental context, or a mechanism-specific target and representation.

## Reproduction

Requirements:

- Python 3.11 or newer
- MMseqs2 on `PATH`, or its executable specified through `MMSEQS_BIN`
- Optional ESM and surface dependencies defined in `pyproject.toml`

```sh
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -e '.[dev,esm,surface]'
gzip -dk data/raw/measurements.csv.gz
shasum -a 256 data/raw/measurements.csv
python3 -m pytest -q
PYTHONPATH=src python3 -m enzsolv.run_benchmark
PYTHONPATH=src python3 -m enzsolv.run_esm
PYTHONPATH=src python3 -m enzsolv.run_surface --workers 2
PYTHONPATH=src python3 -m enzsolv.run_surface_benchmark --phase all
PYTHONPATH=src python3 -m enzsolv.run_interactions
PYTHONPATH=src python3 -m enzsolv.run_diagnostic_controls
PYTHONPATH=src python3 -m enzsolv.run_metadata_diagnostics
PYTHONPATH=src python3 -m enzsolv.run_composition_ablation
```

The checksum printed above must match the digest at the top of this page. The frozen input is documented in [data/README.md](data/README.md). Compact final reports are preserved under `reports/`. Generated structures, model caches and large out-of-fold prediction tables are intentionally omitted and can be reconstructed by the commands above.

## Next decisive experiment

A stronger test would use several unrelated enzymes measured prospectively against the same solvent panel, concentrations, incubation conditions, and activity assay.

The current dataset does not justify another round of broad model search.
