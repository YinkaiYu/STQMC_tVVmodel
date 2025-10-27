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

For a review of the physical properties of the spinless t-V1-V2 model, see:

- [arXiv:1609.01161](https://arxiv.org/pdf/1609.01161)


## Method

This code implements the **Short-Time QMC** algorithm, a state-of-the-art quantum Monte Carlo method introduced in [arXiv:2410.18854](https://arxiv.org/abs/2410.18854). Unlike traditional DQMC approaches, which often suffer from severe fermion sign problems in strongly interacting regimes, Short-Time QMC leverages nonequilibrium critical dynamics to **preempt the sign problem** and enable simulations that are otherwise inaccessible.

The key idea is to perform simulations in the short-time regime of imaginary-time evolution, where the sign problem is much less severe. By combining this with universal scaling theory of relaxation dynamics, Short-Time QMC allows accurate determination of quantum critical points and critical exponents, all while retaining the numerical exactness of QMC.

This framework represents a significant advancement in tackling the fermion sign problem, and opens new possibilities for studying interacting fermion systems near quantum criticality.


## Subroutine Introduction

### main subroutine of Short-Time QMC

- `SuNF`: main subroutine of Short-Time QMC

### Variable Definitions
- `blockc`: define and allocate variables.
- `block_obs`: define and allocate observables.
- `salph`: define and assign values to intermediate variables.

### Set Up Hopping
- `sli`: define functions related to lattice structure.
- `setH`: set up the hopping matrix.
- `sthop`: set up the exponential of the hopping matrix.

### Field Updates
- `upgradeV1`: update auxiliary field for V1 interaction.
- `upgradeV2`: update auxiliary field for V2 interaction.

### Matrix Operations
- `mmuur`: multiply UR matrix by interaction trotter exponential.
- `mmthr`: multiply UR matrix by hopping trotter exponential.
- `mmuul`: multiply UL matrix by interaction trotter exponential.
- `mmthl`: multiply UL matrix by hopping trotter exponential.
- `mmuurm1`: divide UR matrix by interaction trotter exponential.
- `mmthrm1`: divide UR matrix by hopping trotter exponential.
- `mmuulm1`: divide UL matrix by interaction trotter exponential.
- `mmthlm1`: divide UL matrix by hopping trotter exponential.

### Stabilization
- `ortho`: perform SVD orthogonalization of UR matrix.

### Time-Dependent Measurements
- `dyn`: compute time-dependent single-particle Green's function.
- `propr`, `proprm1`: propagate Green's function from τ to τ+1.

### Measurement and Output
- `calcgr`: calculate equal-time single-particle Green's function.
- `obser`: evaluate equal-time observables.
- `obsert`: evaluate unequal-time observables.
- `preq`: output equal-time observables.
- `prtau`: output unequal-time observables.

## License

MIT
