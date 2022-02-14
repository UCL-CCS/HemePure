
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "lb/iolets/InOutLetWK.h"
#include "lb/iolets/BoundaryComms.h"
#include "configuration/SimConfig.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      InOutLetWK::InOutLetWK() :
          InOutLet(), density(1.0), densityNew(1.0), radius(1.0),
          resistance(1.0), capacitance(1.0)
      {
      }

      InOutLet* InOutLetWK::Clone() const
      {
        InOutLetWK* copy = new InOutLetWK(*this);
        return copy;
      }

      InOutLetWK::~InOutLetWK()
      {
      }

      void InOutLetWK::DoComms(const BoundaryCommunicator& boundaryComm, const LatticeTimeStep timeStep)
      {
        if (comms->GetNumProcs() == 1) return;

        comms->Receive(&density);
        comms->Send(&densityNew);
        comms->WaitAllComms();
      }

      LatticeDistance InOutLetWK::GetDistance(const LatticePosition& x) const
      {
        LatticePosition displ = x - position;
        LatticeDistance z = displ.Dot(normal);
        return std::sqrt(displ.GetMagnitudeSquared() - z * z);
      }

      distribn_t InOutLetWK::GetScaleFactor(const LatticePosition& x) const
      {
        // Q = vLocal (0.5pi a**2)(a**2/(a**2 - r**2)
        // where r is the distance from the centreline
        // a is the radius of the circular iolet
        LatticePosition displ = x - position;
        LatticeDistance z = displ.Dot(normal);
        Dimensionless rFactor = (radius * radius)/(radius * radius - (displ.GetMagnitudeSquared() - z * z) );
        return 0.5 * PI * radius * radius * rFactor;
      }
    }
  }
}
