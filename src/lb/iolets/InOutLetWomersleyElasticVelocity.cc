
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "lb/iolets/InOutLetWomersleyElasticVelocity.h"
#include "configuration/SimConfig.h"
#include "util/Bessel.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      const InOutLetWomersleyElasticVelocity::Complex InOutLetWomersleyElasticVelocity::i = Complex(0, 1);
      const InOutLetWomersleyElasticVelocity::Complex InOutLetWomersleyElasticVelocity::iPowThreeHalves =
          pow(i, 1.5);

      InOutLet* InOutLetWomersleyElasticVelocity::Clone() const
      {
        InOutLet* copy = new InOutLetWomersleyElasticVelocity(*this);
        return copy;
      }

      LatticeVelocity InOutLetWomersleyElasticVelocity::GetVelocity(const LatticePosition& x,
                                                             const LatticeTimeStep t) const
      {
        LatticePosition displ = x - position;
        LatticeDistance z = displ.Dot(normal);
        Dimensionless r = sqrt(displ.GetMagnitudeSquared() - z * z);

        double omega = 2.0 * PI / period;
        LatticeDensity density = 1.0;

	Complex lambda = iPowThreeHalves * womersleyNumber;

	Complex g = 2.0 * util::BesselJ1ComplexArgument(lambda) / ( lambda * util::BesselJ0ComplexArgument(lambda));

	//Build terms of frequency equation, we will ASSUME that the wall thickness is 0.1*Radius and rho_wall == rho_fluid
	//Quadratic equation of form Ax^2 + Bx + C = 0
	double thick = 0.1;
	Complex A = (g - 1.0) * (poisson * poisson - 1);
       	Complex B = thick * (g - 1.0) + (2.0 * poisson - 0.5) * g - 2.0;
	Complex C = 2.0*thick + g;

	Complex rt1 = (-B + sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
	Complex rt2 = (-B - sqrt(B * B - 4.0 * A * C)) / (2.0 * A);

	Complex nu;
        if ((rt1.real() > rt2.real()) || ((rt1.real() == rt2.real()) && (rt1.imag() > rt2.imag())))
	{
		nu = rt1;
	}
	else
	{
		nu = rt2;
	}

	Complex M = (2.0 + nu * (2.0 * poisson - 1.0)) / (nu * (2.0 * poisson - g));

	double cR = 1.0 / std::real(1.0 / sqrt(thick * wallYoungs / ((1.0 - poisson * poisson) * nu)));

        Complex besselNumer = util::BesselJ0ComplexArgument(lambda * r / radius);
        Complex besselDenom = util::BesselJ0ComplexArgument(lambda);
       
	LatticeSpeed velocityMagnitude = std::real(i * pressureGradientAmplitude / (density * omega)
            * (1.0 - M * besselNumer / besselDenom) * exp(i * omega * (double(t) - axialPos / cR)));

	//std::cout << "omega " << omega << ", t " << t << ", axial " << axialPos << ", cR " << cR << ", exp " << exp(i * omega * (double(t) - axialPos / cR)) << std::endl;

//	std::cout << "r/radius = " << r/radius << ", lambda " << lambda << ", M " << M << ", bessel division " << besselNumer/besselDenom << std::endl;
        return normal * velocityMagnitude;
      }

      const LatticePressureGradient& InOutLetWomersleyElasticVelocity::GetPressureGradientAmplitude() const
      {
        return pressureGradientAmplitude;
      }

      void InOutLetWomersleyElasticVelocity::SetPressureGradientAmplitude(const LatticePressureGradient& pressGradAmp)
      {
        pressureGradientAmplitude = pressGradAmp;
      }

      const LatticeTime& InOutLetWomersleyElasticVelocity::GetPeriod() const
      {
        return period;
      }

      void InOutLetWomersleyElasticVelocity::SetPeriod(const LatticeTime& per)
      {
        period = per;
      }

      const Dimensionless& InOutLetWomersleyElasticVelocity::GetWomersleyNumber() const
      {
        return womersleyNumber;
      }

      void InOutLetWomersleyElasticVelocity::SetWomersleyNumber(const Dimensionless& womNumber)
      {
        womersleyNumber = womNumber;
      }

      const Dimensionless& InOutLetWomersleyElasticVelocity::GetPoissonRatio() const
      {
        return poisson;
      }

      void InOutLetWomersleyElasticVelocity::SetPoissonRatio(const Dimensionless& pr)
      {
        poisson = pr;
      }

      const LatticePressure& InOutLetWomersleyElasticVelocity::GetWallYoungsModulus() const
      {
        return wallYoungs;
      }

      void InOutLetWomersleyElasticVelocity::SetWallYoungsModulus(const LatticePressure& wallYoungMod)
      {
        wallYoungs = wallYoungMod;
      }

      const LatticeDistance& InOutLetWomersleyElasticVelocity::GetAxialPosition() const
      {
        return axialPos;
      }

      void InOutLetWomersleyElasticVelocity::SetAxialPosition(const LatticeDistance& axPosn)
      {
        axialPos = axPosn;
      }
    }
  }
}
