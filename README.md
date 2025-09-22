# HemePure
HemePure is a modified version of the HemeLB code that improved memory usage and scaling behaviour. It makes use of the lattice Boltzmann method to solve macroscopic blood flow in 3D vascular geometries. The CPU version (this repository) has demonstrated strong scaling behaviour to hundreds of thousands of cores with a sufficiently sized simulation domain. The related [HemePure-GPU](https://github.com/UCL-CCS/HemePure-GPU) allows for execution of the code on a variety of GPU architectures with strong scaling shown on tens of thousands of GPU cards.

This version of the code was developed from the full HemeLB code in late 2018 and has since seen further developments. Publications using some form HemePure:
* Xue, X., Athawale, T. M., McCullough, J. W. S., Lo, S. C., Zacharoudiou, I., Joo, B., ... & Coveney, P. V. (2025). An Uncertainty Visualization Framework for Large-Scale Cardiovascular Flow Simulations: A Case Study on Aortic Stenosis. arXiv preprint arXiv:2508.15420.
* Benemerito, I.,  McCullough, J., Narracott, A., Coveney, P. V.  & Marzo, A. (2025). Comparison of Navier-Stokes and lattice Boltzmann solvers for subject-specific modelling of intracranial aneurysms. Computers in Biology and Medicine, Vol. 197 Part B, 111050.
* Lo, S. C., Zingaro, A., McCullough, J. W. S., Xue, X., Gonzalez-Martin, P., Joo, B., ... & Coveney, P. V. (2025). A multi-component, multi-physics computational model for solving coupled cardiac electromechanics and vascular haemodynamics. Computer Methods in Applied Mechanics and Engineering, 446, 118185.
* Xue, X., McCullough, J. W. S., Lo, S. C., Zacharoudiou, I., Joó, B., & Coveney, P. V. (2024, June). The lattice Boltzmann based large eddy simulations for the stenosis of the aorta. In International Conference on Computational Science (pp. 408-420). Cham: Springer Nature Switzerland.
* McCullough, J. W. S. , & Coveney, P. V. (2024). Uncertainty quantification of the lattice Boltzmann method focussing on studies of human-scale vascular blood flow. Nature Scientific Reports 14, 11317.
* Lo, S. C., McCullough, J. W. S., Xue, X., & Coveney, P. V. (2024). Uncertainty quantification of the impact of peripheral arterial disease on abdominal aortic aneurysms in blood flow simulations. Journal of the Royal Society Interface, 21(213), 20230656.
* McCullough, J. W. S. & Coveney, P. V. (2023). High resolution simulation of basilar artery infarct and flow within the circle of Willis, Scientific Reports 13, 21665.
* Zacharoudiou, I., McCullough, J. W. S. & Coveney, P. V. (2023). Development and performance of a HemeLB GPU code for human-scale blood flow simulation. Computer Physics Communications, 282, 108548.
* Lo, S. C., McCullough, J. W. S. & Coveney, P. V. (2022). Parametric Analysis of an Efficient Boundary Condition to Control Outlet Flow Rates in Large Arterial Networks. Scientific Reports 12, 19092.
* McCullough, J. W. S. & Coveney, P. V. (2021). An efficient, localised approach for the simulation of elastic blood vessels using the lattice Boltzmann method. Scientific Reports 11, 24260.
* McCullough, J. W. S. & Coveney, P. V. (2021). High fidelity blood flow in patient-specific arteriovenous fistula. Scientific Reports 11, 22301.
* McCullough, J. W. S., Richardson, R. A., Patronis, A., Halver, R., Marshall, R., Ruefenacht, M., Wylie, B. J. N., Odaker, T., Wiedemann, M., Lloyd, B., Neufeld, E., Sutmann, G., Skjellum, A., Kranzlmüller, D., & Coveney, P. V. (2020). Towards blood flow in the virtual human: efficient self-coupling of HemeLB. Journal of the Royal Society Interface Focus 11, 20190119.

## Features #
The CPU version of HemePure can be executed with the following functionality. Some must be specified at the compilation of the `hemepure` executable. Simulations are conducted using a D3Q19 lattice stencil. 

Collision kernels (compile time):
* LBGK - Single relaxation time
* TRT - Two relaxation time
* LBGK + LES - Inclusion of large eddy simulation approximation to assist in higher Re flow modelling.

Inlet/Outlet boundary conditions (compile time):
* Pressure
  - Sinusoidal profile (constant pressure enabled using)
  - Transient profile
  - Windkessel
  - Sponge layer (outlets) - Acts on a pressure outlet but modifies the viscosity near the outlet to increase stability of the simulation (defined with collision kernel).
* Velocity
   - Constant magnitude with parabolic profile for circular inlets
   - Transient profile with parabolic profile for circular inlets
   - Transient profile with Poiseuille-like profile for non-circular inlets

Wall boundary conditions (compile time):
* Bounceback - simple rigid walls
* BFL - Modified bounceback approach with
* GZS - Extrapolation based wall condition
* GZSElastic - Method for approximating the effect of elastic walls

Intrinsics for vectorisation, may improve performance on CPU (compile time):
* SSE3
* AVX
* AVX512

Data output (run time):
- Extraction of data at point, line, plane, inlets, outlets, wall surface, surfaces within a sphere, whole domain.
- Checkpoint restart from written data file

## Compilation #
Build dependencies before attempting to build HemePure. Once dependencies are built, they do not need to be recompiled for (re)compilation of the source code. The following steps can be followed to build the dependency and source code mannually, the FullBuild.sh file collects these into a single location. This file may need to be modified to reflect the settings, defaults and software available on a given machine. This file can be modified to provide alternative functionality and boundary conditions.

For benchmarking, it is recommended to use the SRCbuild_Benchmark function in FullBuild.sh to compile the source code. This enables pressure boundary conditions at inlets and outlets of the domain, and enables more reliable simulation performance and stability for evaluating performance on a given computer.

### DEPENDENCIES #
1) Create `dep/build/`.
2) In `dep/build/` run `ccmake -B. -H../' or 'ccmake ..`.
3) Configure using CMake.
4) Run `make` in `dep/build/`.

### SOURCE #
1) Create `src/build/`.
2) In `src/build/` run `ccmake -B. -H../` or `ccmake ..`.
3) Configure using CMake.
4) Run `make` in `src/build/`.

## EXAMPLE CASES #
The [cases](cases) folder provides some example input files for running jobs with HemePure. Some are simple simulations, some help illustrate the usage of particular features of the code.
Execution should occur with: `mpirun -np xx <HemePure binary> -in input.xml -out results`, ensure that you are running with enough MPI ranks to satisfy HEMELB_READING_GROUP_SIZE + 1 (by default, at least 3 ranks). (`mpirun` may need to be replaced with an alternative call based on machine settings). Note that the simulation will not proceed if a `results` folder with the same name already exists in the run directory - delete the existing folder, or specify a new folder in which to store results, prior to execution.

Two simple examples are included:
* [Pipe](case/pipe): pressure-driven flow in a pipe.
* [Bifurcation](cases/bifurcation): pressure-driven flow in a Y-bifurcation (provided at two resolutions).

Further examples provide extra input domains or illustrate extra functionality:
* Data extraction
* Checkpointing
* Elastic walls
* LES and sponge layer outlets
* Windkessel outlets
  
## DEVELOPMENT #
A few parts of the code are under heavy development.

<!--
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
-->

The new shared memory implementation (found in `GeometryReader.cc`) only shares `principalProcForEachBlock`. Further savings can be achieved by sharing other common variables, e.g. `blockInformation`. Also note that only the `MPI_WIN_UNIFIED` memory model is supported (see https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node278.htm#Node278). **If the machine in use does not support this, it will not be possible to use this implementation**. EPCC ARCHER seems to have issues with the current implementation. Further investigation is necessary.

Verification checks have yet to be reintroduced following conversion of sequential containers `std::vector` to associative containers `std::unordered_map`. Checks are performed during loading and decomposition. They are currently commented out and disabled at compile time. Cross-reference with master branch of HemeLB. Checks will need to be rewritten.

## NOTES #
- Compilation of some dependencies is disabled:
  - CppUnit: compilation of unit tests is disabled.
- This version of HemePure is in continuous development - some features are in testing.
- Disable optimised decomposition for memory-intensive runs. ParMETIS is resource hungry.
