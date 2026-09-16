# Can fly-connectome topology provide a useful neural-network architecture?

| Field | Value |
| --- | --- |
| Stage | idea |
| Outcome | broad version already covered by prior work and too large for the intended side project |
| Investigated | 2026-09-15 to 2026-09-16 |
| Experiment performed | no |
| Empirical result | none |

## Summary

I considered reducing the complete male fruit-fly connectome to roughly 1,024 nodes and using its topology as a recurrent neural network. The plan was to compare it with standard recurrent architectures and random graphs at the same size and computational budget.

I did not run the benchmark. Several groups had already tested fly-derived reservoirs, connectome-constrained networks, and randomized topology controls. More recent work also showed that some apparent advantages disappear when the control graph preserves the directed input and output degree of every node.

The remaining question is narrower: which measured features of a connectome contribute to a specific computation after degree, geometry, dynamics, input placement, and coarse-graining have been controlled? Answering that would require a larger causal study than the project I had in mind.

## Original hypothesis

A reduced fly-derived topology may improve memory, temporal prediction, robustness, or sample efficiency on some tasks when compared with artificial graphs and conventional recurrent networks under matched conditions.

## Proposed test

Every graph would use the same recurrent update rule, with only a linear readout trained in the first experiment. The MaleCNS-derived graph would be compared with:

- directed degree-preserving rewires;
- weight and sign shuffles;
- density-matched random graphs;
- spatially constrained null graphs;
- artificial modular graphs;
- standard echo-state reservoirs.

The comparison would cover several network sizes, multiple input and output placements, and at least two coarse-graining methods. A region-by-task analysis would test whether topology taken from a particular biological subsystem performs differently on aligned and mismatched tasks. An embodied extension would vary the neural topology, body morphology, and sensor-motor mapping independently.

## Prior work found during the check

- [Morra and Daley (2022)](https://arxiv.org/abs/2201.09359) imposed a female fly connectome on an echo-state network.
- [Morra et al. (2023)](https://arxiv.org/abs/2306.01885) studied multifunctionality in a lateral-horn reservoir.
- [Costi et al. (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12109256/) tested reduced fly-connectome reservoirs and topology-weight variants.
- A [connectome-constrained visual model](https://www.nature.com/articles/s41586-024-07939-3) combined measured fly wiring with task optimization.
- [FLYNN](https://arxiv.org/abs/2607.00025) used fly-derived topology for robot navigation.
- [Dhiman (2026)](https://arxiv.org/abs/2604.04033) reported that several apparent advantages largely disappeared under stronger controls.
- The public [Connectome Architecture Benchmark](https://github.com/tejasnaladala/connectome-bpu) implements multiple connectomes, tasks, and null graph families. Its corrected complete result was pending when I checked it.

## Outcome

No computational result exists. The literature check showed that the original connectome-versus-random-network comparison was no longer a useful standalone project. A new experiment would need causal feature attribution or a region-by-task analysis with matched null graphs.

## Why I stopped

The broad benchmark had already been attempted in several forms, while the controls required for a useful new result expanded the project far beyond a small computational side project. I stopped at the design and prior-art stage.

Biological topology may still contain useful computational biases. This entry records why the broad test was not run and what a more informative study would need to isolate.
