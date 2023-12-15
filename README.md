# Computing-algebraic-degrees-of-phylogenetic-varieties

Within this repository, you'll find code featuring functions dedicated to computing algebraic degrees, including EDdeg, gEDdeg, and MLdeg. These calculations apply to both affine and projective phylogenetic varieties associated with diverse trees and models within the Phylogenetic framework. The computations are executed using HomotopyContinuation.jl package and Oscar.jl.

👉 `PhylogeneticData.jl` the primary goal of this file is to generate the matrix representing the monomial parametrization of the (affine and projective) variety, given the tree and the model.

👉 `Algebraic_Invariants_functions.jl` includes functions for computing algebraic degrees, utilizing both the `solve()` function and `monodromy_solve()` in `HomotopyContinuation.jl`.

👉 `Symbolic_Computations.jl` includes functions for computing the EDdeg and gEDdeg using the `msolve` functionality in `Oscar.jl`.

## Examples

You can find various examples at 👉 `Benchmarks.jl`

## 👷 Development team

- Luis David Garcia Puente <lgarciapuente@coloradocollege.edu>
- Marina Garrote-López <marinagarrotelopez@gmail.com>
- Elima Shehu <elima.shehu@mis.mpg.de>
