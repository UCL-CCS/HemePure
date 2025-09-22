**Two-Element Windkessel Model**

How to use:
- In the compiler option, set HEMELB_INLET_BOUNDARY or HEMELB_OUTLET_BOUNDARY to NASHZEROTHORDERPRESSUREIOLET or YANGPRESSUREIOLET.
- In the input file, select “pressure” or “yangpressure” as the boundary condition type and “WK2” as the subtype.

Example input block:
```
<outlet>
	<condition subtype="WK2" type="pressure">
		<R value="4.18e+09" units="kg/m^4*s"/>
		<C value="2.82e-10" units="m^4*s^2/kg"/>
		<radius value="0.00427" units="m"/>
		<area value="5.6e-05" units="m^2"/>
	</condition>
	<normal units="dimensionless" value="(0.555,0.308,-0.773)"/>
	<position units="lattice" value="(131.4,290.8,2817.1)"/>
</outlet>
```
Pull requests: #2, #7, #8, #13, #15, #19 in UCL_CCS/HemePure (possibly more)

Used in: 
* https://doi.org/10.1038/s41598-022-21923-9
* https://doi.org/10.1098/rsif.2023.0656


===========================================================================================


**Three-Element Windkessel Model** 

How to use:
- In the compiler option, set HEMELB_INLET_BOUNDARY or HEMELB_OUTLET_BOUNDARY to NASHZEROTHORDERPRESSUREIOLET or YANGPRESSUREIOLET.
- In the input file, select “pressure” or “yangpressure” as the boundary condition type and “WK3” as the subtype.

Example input block:
```
<outlet>
	<condition subtype="WK3" type="pressure">
		<Rc value="2e6" units="kg/m^4*s"/>
		<Rp value="1e9" units="kg/m^4*s"/>
		<Cp value="1.5e-8" units="m^4*s^2/kg"/>
		<radius value="0.02" units="m"/>
		<area value="0.00126" units="m^2"/>
	</condition>
	<normal units="dimensionless" value="(0,0,-1)"/>
	<position units="lattice" value="(216.5,216.5,4273)"/>
</outlet>
```

Pull requests: #33 in UCL_CCS/HemePure

Used in: 
* https://doi.org/10.1016/j.cma.2025.118185

