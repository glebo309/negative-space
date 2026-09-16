# Publication and prior-art scope check

Checked 2026-09-11 using targeted searches and primary resource/publisher pages.
This is a bounded novelty check, not a systematic review or proof of absence.

## Resource provenance

The [OrganiZymeDB Zenodo record](https://zenodo.org/records/21874133) lists Romain
Debruyne and Fabrizio Pucci as creators. It was created on 10 August 2026. That is
a repository-record date, not evidence that no one previously studied prediction
of solvent tolerance. The record remains embargoed pending the associated
manuscript, with an explicitly placeholder release deadline. The public website
provides the measurements used here. Its [citation section](https://organizymedb.org/about)
still has placeholder publication fields at this check.

The website itself explicitly acknowledges uncertain assay/control compositions
and intentionally distinguishes stability assays from combined activity/stability.
Our audit quantifies those existing annotations in our selected cohort. It should
not be portrayed as discovering an error the database creators concealed.
The live About page now describes ColabFold/AlphaFold2 models; a search-index
snippet returned an older AlphaFold3 description. Our calculations use the cached
published coordinate files, with exact sequence checks, and do not depend on an
unverified generator-version assertion.

## Related work already exists

- Dachuri, Lee and Jang, 2018, [Organic Solvent-Tolerant Esterase from Sphingomonas
  glacialis Based on Amino Acid Composition Analysis](https://pubmed.ncbi.nlm.nih.gov/30176710/),
  DOI 10.4014/jmb.1806.06032, already uses composition to identify a tolerant
  esterase. A narrower composition result here cannot be framed as the discovery
  that residue frequencies can contain solvent-tolerance information.
- Ma, Li, Zhang and Xu, 2023, [Methanol tolerance upgrading of Proteus mirabilis
  lipase by machine learning-assisted directed evolution](https://journal.hep.com.cn/smab/EN/1160585065489032134),
  DOI 10.1007/s43393-023-00179-y. Models trained on 266 combinatorial mutants
  guided a restricted library in one lipase. This is direct ML-assisted solvent
  tolerance engineering, but not transfer across unfamiliar enzyme families.
- [Prediction of the solvent affecting site and the computational design of
  stable Candida antarctica lipase B in a hydrophilic organic solvent](https://pubmed.ncbi.nlm.nih.gov/23178554/)
  used molecular dynamics and surface-site redesign, followed by mutation tests.
  The broad idea of surface-guided solvent tolerance design is not new.
- [ESM-Ezy: a deep learning strategy for the mining of novel multicopper oxidases
  with superior properties](https://pmc.ncbi.nlm.nih.gov/articles/PMC11972304/)
  uses an ESM-based mining strategy with experimental property testing, including
  solvent tolerance. It is adjacent protein-language-model prior art, not the
  same matched-solvent, family-and-publication-held-out supervised benchmark.

## Defensible publication boundary

The potentially distinctive question is whether protein representations improve
solvent choice across independent enzymes and publications in a heterogeneous
literature resource. Neither a recent database nor a negative run guarantees a
paper. The completed work can support an auditable transfer benchmark and a
measurement-context limitations analysis. It does not yet establish a useful
predictor, a universal absence of sequence effects, or a new surface mechanism.

A negative computational note would need a substantive benchmark contribution,
careful positioning against the resource authors' forthcoming work, clear
reproducibility, and a venue interested in such results. Those are publication
development choices, not reasons to initiate a model-validation bench panel when
the model has not passed its computational gate. No authors were contacted.
