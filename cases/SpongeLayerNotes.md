**Sponge Layer Implementations**

How to use:
- In the compiler option, set HEMELB_KERNEL to “LBGKSL”, “TRTSL”, or “LBGKLESSL”.
- In the input file, add a “sponge_layer” block under “initialconditions”.

Example input block:
```
<initialconditions>
	<sponge_layer>
		<viscosity_ratio value="1000" units="dimensionless"/>
		<width value="0.04" units="m"/>
		<lifetime value="99999" units="lattice"/>
	</sponge_layer>
</initialconditions>
```

Pull requests: #25, #27, #33 in UCL_CCS/HemePure (possibly more)

Used in: 
* https://doi.org/10.1016/j.cma.2025.118185
* https://doi.org/10.1007/978-3-031-63775-9_30
* https://doi.org/10.48550/arXiv.2508.15420


Remarks:
- “width”: This parameter actually corresponds to the radius of the sponge layer. The current naming is potential misleading and should be revised. (The terminology originated from an earlier plan to implement sponge layers of different shapes.)
- “lifetime”: This parameter specifies the duration of the sponge layer from the beginning of the simulation. The sponge layer begins to decay halfway through its lifetime and disappears completely at the end (see line 127 in https://github.com/UCL-CCS/HemePure/blob/master/src/lb/kernels/LBGKSpongeLayer.h). This feature is found to be impractical in real simulations as instabilities will still occur once the sponge layer disappears in cases where such a layer is required. It was implemented solely to verify that the flow can return to its original state.
