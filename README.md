# sparse-tensor-cp
This package is based on the work in the following article:
> Pengfei Huang, Minru Bai.
> An ideal-sparse generalized moment problem reformulation for completely positive tensor decomposition exploiting maximal cliques of multi-hypergraphs.
> (https://arxiv.org/abs/2505.15056)

This code is designed for the completely positive tensor decomposition problem via a semidefinite relaxation method. It checks whether a symmetric tensor is completely positive and, if so, provides a completely positive decomposition.

### Package requirements:
- CSV
- DataFrames 
- DynamicPolynomials
- HomotopyContinuation 
- JuMP 
- Mosek 
- MosekTools 
- MultivariateMoments 
- MultivariatePolynomials 
- RowEchelon 
- TSSOS = "1.3.4"
- TypedPolynomials
## Getting started
To generate and run all examples in Julia, one can 
``` julia
include("src\\demo.jl")
```

## Acknowledgement
- This code builds upon "ju-cp-rank" (https://github.com/JAndriesJ/ju-cp-rank).
> M. Korda, M. Laurent, V. Magron, and A. Steenkamp. 
> Exploiting ideal-sparsity in the generalized moment problem with application to matrix factorization ranks.
> Mathematical Programming, 2023. https://doi.org/10.1007/s10107-023-01993-x.

Our code is mainly focused on tensor cases, with an emphasis on how to generate the maximal cliques of the associated multi-hypergraph ("src\\tensorGraphs.jl") and handle tensor indices when the order is greater than two.

- We are grateful for the work of Ni Guyan and Li Ying to deal with subscripts of tensors.
> Ni Guyan and Li Ying.
> A Semidefinite Relaxation Method for Partially Symmetric Tensor Decomposition.
> Mathematics of Operations Research, 2022, 10.1287/moor.2021.1231.
