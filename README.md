# Short-Time QMC for the Spinless \( t\text{-}V_1\text{-}V_2 \) Model

This repository contains an open-source implementation of the **Short-Time QMC** algorithm, a short-time version of the determinantal quantum Monte Carlo (DQMC) framework proposed in [arXiv:2410.18854](https://arxiv.org/abs/2410.18854).

## Model

We study the spinless fermion model with nearest-neighbor (V₁) and next-nearest-neighbor (V₂) interactions, defined by the Hamiltonian:

\[
H = -t \sum_{\langle ij \rangle} c_i^\dagger c_j + V_1 \sum_{\langle ij \rangle} n_i n_j + V_2 \sum_{\langle\langle ij \rangle\rangle} n_i n_j
\]

where:
- \( \langle ij \rangle \): nearest neighbors  
- \( \langle\langle ij \rangle\rangle \): next-nearest neighbors  
- \( c_i^\dagger \), \( c_j \): fermionic creation/annihilation operators  
- \( n_i = c_i^\dagger c_i \): number operator

This model suffers from the fermion sign problem when \( V_1 < 0 \) or \( V_2 > 0 \).

For a comprehensive review of the physical properties of the spinless \( t\text{-}V_1\text{-}V_2 \) model, refer to:

- [arXiv:1609.01161](https://arxiv.org/pdf/1609.01161)

## Method

Our implementation is based on the **Short-Time QMC** framework introduced in:

- [arXiv:2410.18854](https://arxiv.org/abs/2410.18854): A new state-of-the-art QMC framework designed to mitigate the sign problem in challenging interaction regimes.