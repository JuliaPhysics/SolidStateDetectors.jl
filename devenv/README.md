# Julia development environments

This directory contains the Julia project environments `cpu` and `cuda`
that can be used when developing SolidStateDetectors with Julia >= v1.11.
They contain all direct, test and doc-gen dependencies of SolidStateDetectors,
plus BenchmarkTools and Cthulhu.

Note: These environments can't be used with Julia versions <= v1.10, as they
use a `[sources]` section in the `Project.toml` to ensure SolidStateDetectors
is loaded from the local source directory, this Pkg feature was introduced
in Julia v1.11.
