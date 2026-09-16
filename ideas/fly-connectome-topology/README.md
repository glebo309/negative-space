# Can fly-connectome topology provide a useful neural-network architecture?

| Field | Value |
| --- | --- |
| Stage | idea |
| Outcome | broad version superseded by prior work and scope |
| Investigated | 2026-09-15 to 2026-09-16 |
| Experiment performed | no |
| Empirical result | none |

## Short result

The proposed broad benchmark was not performed. Literature and repository searches showed that reduced fly connectomes had already been used as reservoirs, connectome topology had already been compared with randomized controls, and current projects were pursuing multi-task null-graph benchmarks. Stronger controls also showed that some reported advantages disappear after matching initialization and directed degree sequence.

The remaining defensible questions require a substantially larger causal or embodied study and are not a small hobby project.

## Question

Can a reduced topology extracted from the complete male fruit-fly central nervous system outperform conventional recurrent architectures or matched artificial graphs at the same computational budget?

## Hypothesis

A carefully reduced fly-derived topology will exhibit task-dependent advantages in memory, temporal prediction, robustness or sample efficiency compared with artificial graphs and conventional recurrent networks when size, edge count, weights, dynamics, training and compute are controlled.

## Proposed experiment

Use the same recurrent update rule for every graph and initially train only a linear readout. Compare a MaleCNS-derived topology against directed degree-preserving rewires, weight and sign shuffles, density-matched random graphs, spatially constrained nulls, artificial modular graphs and standard echo-state reservoirs.

Test several reduction sizes and repeat the analysis across functionally distinct fly regions, aligned and deliberately mismatched tasks, multiple input-output placements and at least two principled coarsening methods.

An embodied extension would independently vary neural topology, body morphology and anatomical versus shuffled sensor-motor mappings.

## Prior-art result

- [Morra and Daley (2022)](https://arxiv.org/abs/2201.09359) imposed a female fly connectome on an echo-state network.
- [Morra et al. (2023)](https://arxiv.org/abs/2306.01885) studied multifunctionality in a lateral-horn reservoir.
- [Costi et al. (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12109256/) tested reduced fly-connectome reservoirs and topology-weight variants.
- A [connectome-constrained visual model](https://www.nature.com/articles/s41586-024-07939-3) combined measured fly wiring with task optimization.
- [FLYNN](https://arxiv.org/abs/2607.00025) used fly-derived topology for robot navigation.
- [Dhiman (2026)](https://arxiv.org/abs/2604.04033) found that several apparent advantages largely disappeared under stronger controls.
- The public [Connectome Architecture Benchmark](https://github.com/tejasnaladala/connectome-bpu) implements multiple connectomes, tasks and null graph families. Its corrected complete result was pending at the time of review.

## Results

No computational experiment was performed.

The result is a scope decision. The simple reduced-connectome-versus-random-network experiment is not novel. Fair attribution requires controls for degree, geometry, initialization, input-output placement, weights and coarse-graining. A defensible project would need causal feature attribution or a region-by-task and body-topology interaction.

## Conclusion

Do not implement the broad 1,024-node fly-topology benchmark as a publication-oriented side project. The useful record is the original hypothesis, the fair design it would require and the narrower causal questions that remain.

This does not show that biological topology lacks computational value. It records that no experiment was run and that the broad question was already occupied.
