# Short-Time QMC for the Spinless t-V1-V2 Model

This repository contains an open-source implementation of the **Short-Time QMC** algorithm, a short-time version of determinantal quantum Monte Carlo (DQMC) introduced in [arXiv:2410.18854](https://arxiv.org/abs/2410.18854).

## Model

We study the spinless fermion model with nearest-neighbor (V1) and next-nearest-neighbor (V2) interactions on a square lattice.

The Hamiltonian is:

    H = -t ∑⟨i,j⟩ c†_i c_j + V1 ∑⟨i,j⟩ n_i n_j + V2 ∑⟨⟨i,j⟩⟩ n_i n_j

where:
- ⟨i,j⟩ denotes nearest-neighbor pairs,
- ⟨⟨i,j⟩⟩ denotes next-nearest-neighbor pairs,
- c†_i and c_j are fermionic creation and annihilation operators,
- n_i = c†_i c_i is the number operator.

This model exhibits a fermion sign problem when V1 < 0 or V2 > 0.

## Method

Our implementation uses the **Short-Time QMC** framework, which allows controlled simulation of regimes with sign problems, as introduced in:

- [arXiv:2410.18854](https://arxiv.org/abs/2410.18854)

For a review of the physical properties of the spinless t-V1-V2 model, see:

- [arXiv:1609.01161](https://arxiv.org/pdf/1609.01161)

## License

MIT
