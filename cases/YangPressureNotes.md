**Yang Pressure Boundary Condition** 

How to use: 
* In the compiler option, set HEMELB_INLET_BOUNDARY or HEMELB_OUTLET_BOUNDARY to YANGPRESSUREIOLET.
* In the input file, specify “yangpressure” as the boundary condition type for an inlet or outlet. The choice of subtype remains the same as for “pressure”.

Example input block:
```
<inlet>
    <condition subtype="cosine" type="yangpressure">
      <amplitude units="mmHg" value="0.0"/>
      <mean units="mmHg" value="0.01"/>
      <phase units="rad" value="0.0"/>
      <period units="s" value="1"/>
      <radius value="0.001" units="m"/>
      <area value="3.14e-06" units="m^2"/>
    </condition>
  <normal units="dimensionless" value="(0,-0.342,0.940)"/>
  <position units="lattice" value="(13,46.6,6.42)"/>
</inlet>
```

Pull requests: #4, #28 in UCL_CCS/HemePure

Used in: 
* Sharp Lo’s PhD thesis

Remarks:

- This implementation has been verified against laminar flow in a cylinder tilted at an arbitrary angle, demonstrating second-order accuracy for both velocity and pressure. Data available at https://doi.org/10.5522/04/c.7811009

- However, challenges may arise when applying it to realistic vascular geometries with irregular boundary planes. The issue stems from the possible lack of sufficient lattice sites along a given direction to support interpolation/extrapolation.
