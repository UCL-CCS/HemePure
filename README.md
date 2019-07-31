# HemePure
Build dependencies before attempting to build HemePure.

## DEPENDENCIES #
1) Create `dep/build/`.
2) In `dep/build/` run `ccmake -B. -H../' or 'ccmake ..`.
3) Configure using CMake.
4) Run `make` in `dep/build/`.

## SOURCE #
1) Create `src/build/`.
2) In `src/build/` run `ccmake -B. -H../` or `ccmake ..`.
3) Configure using CMake.
4) Run `make` in `src/build/`.

## CASES #
- Two simple examples are included:
  - pressure-driven flow in a pipe.
  - pressure-driven flow in a wye joint (provided at two resolutions).
- Run with: `mpirun -np xx <HemePure binary> -in input.xml`.
- Case should run with the default build options.
- Ensure that you are running with enough MPI ranks to satisfy HEMELB_READING_GROUP_SIZE.

## DEVELOPMENT #
A few parts of the code are under heavy development.

Although the current version of ALL has been fully incorporated, ALL itself is undergoing development. As a result, it is likely that the implementation in `src/geometry/decomposition/BasicDecomposition.cc` (in the function `BasicDecomposition::DecomposeBlock()`) will need to be modified. Once development is complete, ALL should be enabled through cmake. After rotation of the geometry (necessary for ALL decomposition), HemeLB and ALL interact through

`points.push_back(p);`

where

```
#ifdef HEMELB_USE_GMYPLUS
	ALL_Point<double> p(3,block_coords,blockWeights.at(blockNumber));
#else
	ALL_Point<double> p(3,block_coords,1);
#endif
```

This is dependant on block weighting (included in `.gmy+` files). To use `BasicDecomposition::DecomposeBlock()`, uncomment the call in `GeometryReader.cc`. Note, geometry rotation is necessary so that layers of sites and partition boundary surfaces are not parallel. Rotation is performed in `BasicDecomposition::RotateAndAllocate`.

The new shared memory implementation (found in `GeometryReader.cc`) only shares `principalProcForEachBlock`. Further savings can be achieved by sharing other common variables, e.g. `blockInformation`. Also note that only the `MPI_WIN_UNIFIED` memory model is supported (see https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node278.htm#Node278). **If the machine in use does not support this, it will not be possible to use this implementation**. EPCC ARCHER seems to have issues with the current implementation. Further investigation is necessary.

Verification checks have yet to be reintroduced following conversion of sequential containers `std::vector` to associative containers `std::unordered_map`. Checks are performed during loading and decomposition. They are currently commented out and disabled at compile time. Cross-reference with master branch of HemeLB. Checks will need to be rewritten.

## NOTES #
- Compilation of some dependencies is disabled:
  - CppUnit: compilation of unit tests is disabled.
- This version of HemePure is in continuous development - some features are in testing.
- Disable optimised decomposition for memory-intensive runs. ParMETIS is resource hungry.
