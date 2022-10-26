
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_LBMPARAMETERS_H
#define HEMELB_LB_LBMPARAMETERS_H

#include <cmath>
#include "constants.h"
#include <vector>
#include <cassert>

namespace hemelb
{
  namespace lb
  {
    enum StressTypes
    {
      VonMises = 0,
      ShearStress = 1,
      IgnoreStress = 2
    };

    struct LbmParameters
    {
      public:
        LbmParameters(PhysicalTime timeStepLength, PhysicalDistance voxelSize)
        {
          Update(timeStepLength, voxelSize);
        }

        void Update(PhysicalTime timeStepLength, PhysicalDistance voxelSizeMetres)
        {
          timestep = timeStepLength;
          voxelSize = voxelSizeMetres;
          tau = 0.5
              + (timeStepLength * BLOOD_VISCOSITY_Pa_s / BLOOD_DENSITY_Kg_per_m3)
                  / (Cs2 * voxelSize * voxelSize);

          omega = -1.0 / tau;
          beta = -1.0 / (2.0 * tau);
        }

        PhysicalTime GetTimeStep() const
        {
          return timestep;
        }

        PhysicalDistance GetVoxelSize() const
        {
          return voxelSize;
        }

        distribn_t GetOmega() const
        {
          return omega;
        }

        distribn_t GetTau() const
        {
          return tau;
        }

        distribn_t GetBeta() const
        {
          return beta;
        }

        // A parameter used in the relaxation of some collision operators.
        void SetRelaxationParameter(const distribn_t& param)
        {
          relaxationParameter = param;
        }
        distribn_t GetRelaxationParameter() const
        {
          return relaxationParameter;
        }

        StressTypes StressType;
        distribn_t ElasticWallStiffness;
        distribn_t BoundaryVelocityRatio;

      private:
        PhysicalTime timestep;
        PhysicalDistance voxelSize;
        distribn_t omega;
        distribn_t tau;
        distribn_t beta; ///< Viscous dissipation in ELBM
        distribn_t relaxationParameter;
    };
  }
}

#endif //HEMELB_LB_LBMPARAMETERS_H
