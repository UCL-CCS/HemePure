
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_IOLETS_INOUTLETWOMERSLEYELASTICVELOCITY_H
#define HEMELB_LB_IOLETS_INOUTLETWOMERSLEYELASTICVELOCITY_H
#include "lb/iolets/InOutLetVelocity.h"
#include <complex>

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      /**
       * This class implements the InOutLetVelocity interface. It imposes a Womersley velocity profile
       * at a given inlet/outlet of a simulation of laminar flow in a cylinder with elastic walls of average radius R and zero
       * average sinusoidal pressure gradient given by
       *
       *   pressureGradientAmplitude * sin(2*pi*t/period)
       *
       * where t is the current time. See physics validation paper for the analytical expression of
       * velocity implemented by GetVelocity(x, t) which is also a function of pressureGradientAmplitude,
       * radius, period, womersleyNumber and elastic wall parameters (see PhD thesis by Figueroa 2006 - A coupled-momentum method to model blood flow and vessel deformation in human arteries: Applications in disease research and simulation-based medical planning (Stanford)).
       *
       * Use at both ends (?) with length paramter set according to length of cylinder modelled
       */
      class InOutLetWomersleyElasticVelocity : public InOutLetVelocity
      {
        public:

          /**
           * Returns a copy of the current iolet. The caller is responsible for freeing that memory.
           *
           * @return copy of the current iolet
           */
          InOutLet* Clone() const;

          /**
           * Get Womersley velocity for a given time and position.
           *
           * @param x lattice site position
           * @param t time
           * @return velocity
           */
          LatticeVelocity GetVelocity(const LatticePosition& x, const LatticeTimeStep t) const;

          /**
           * Get the amplitude of the zero average pressure gradient sine wave imposed.
           *
           * @return pressure gradient amplitude
           */
          const LatticePressureGradient& GetPressureGradientAmplitude() const;

          /**
           * Set the amplitude of the zero average pressure gradient sine wave imposed.
           *
           * @param pressGradAmp pressure gradient amplitude
           */
          void SetPressureGradientAmplitude(const LatticePressureGradient& pressGradAmp);

          /**
           * Get the period of the zero average pressure gradient sine wave imposed
           *
           * @return sinusoidal pressure profile period
           */
          const LatticeTime& GetPeriod() const;

          /**
           * Set the period of the zero average pressure gradient sine wave imposed
           *
           * @param per sinusoidal pressure profile period
           */
          void SetPeriod(const LatticeTime& per);

          /**
           * Get the Womersley number characterising the pulsatile flow simulation
           *
           * @return Womersley number
           */
          const Dimensionless& GetWomersleyNumber() const;

          /**
           * Set the Womersley number characterising the pulsatile flow simulation
           *
           * @param womNumber Womersley number
           */
          void SetWomersleyNumber(const Dimensionless& womNumber);
            
	  /**
           * Get the Poisson Ratio characterising the elastic wall
           *
           * @return Poisson Ratio
           */
          const Dimensionless& GetPoissonRatio() const;

          /**
           * Set the Poisson Ratio characterising the elastic wall
           *
           * @param poisson PoissonRatio
           */
          void SetPoissonRatio(const Dimensionless& poisson);

  
	  /**
           * Get the Wall Youngs Modulus characterising the elastic wall
           *
           * @return Wall Youngs Modulus
           */
          const LatticePressure& GetWallYoungsModulus() const;

          /**
           * Set the Wall Youngs Modulsu characterising the elastic wall
           *
           * @param wallYoungs Wall Youngs Modulus
           */
          void SetWallYoungsModulus(const LatticePressure& wallYoungs);

  
	  /**
           * Get the Axial Position of boundary
           *
           * @return Axial Position
           */
          const LatticeDistance& GetAxialPosition() const;

          /**
           * Set the Axial Position of boundary
           *
           * @param axialPos Axial Position
           */
          void SetAxialPosition(const LatticeDistance& axialPos);

        private:
          typedef std::complex<double> Complex;
          static const Complex i;
          static const Complex iPowThreeHalves;
          LatticePressureGradient pressureGradientAmplitude; ///< See class documentation
          LatticeTime period; ///< See class documentation
          double womersleyNumber; ///< See class documentation
          double poisson; ///< See class documentation
	  LatticeDistance axialPos;
	  LatticePressure wallYoungs;
      };
    }
  }
}
#endif // HEMELB_LB_IOLETS_INOUTLETWOMERSLEYELASTICVELOCITY_H
