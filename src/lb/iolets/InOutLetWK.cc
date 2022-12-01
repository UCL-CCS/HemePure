
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
          InOutLet(), radius(1.0), area(1.0), resistance(1.0), capacitance(1.0),
          density(1.0), densityNew(1.0), flowRate(0.0), flowRateNew(0.0),
          siteCount(0)
      {
      }

      InOutLet* InOutLetWK::Clone() const
      {
        InOutLet* copy = new InOutLetWK(*this);
        return copy;
      }

      InOutLetWK::~InOutLetWK()
      {
      }

      void InOutLetWK::DoComms(const BoundaryCommunicator& boundaryComm, const LatticeTimeStep timeStep)
      {
        /**
         * Here the send and receive requests are placed. The message is received at or before the wait
         * barrier set by BoundaryValues::FinishReceive().
         */
        comms->Receive(&density);
        comms->Send(&densityNew);

        // Here the reductions are blocking communications; they have to be made non-blocking
        const BoundaryCommunicator& bcComm = comms->GetCommunicator();
        flowRate = bcComm.Reduce(flowRateNew, MPI_SUM, bcComm.GetBCProcRank());
        siteCount = bcComm.Reduce(siteCount, MPI_SUM, bcComm.GetBCProcRank());
        if (siteCount != 0)
        {
          flowRate = flowRate * area / siteCount;
        }
        flowRateNew = 0.0;
        siteCount = 0;
      }

      LatticeDistance InOutLetWK::GetDistanceSquared(const LatticePosition& x) const
      {
        LatticePosition displ = x - position;
        LatticeDistance z = displ.Dot(normal);
        return displ.GetMagnitudeSquared() - z * z;
      }

      distribn_t InOutLetWK::GetScaleFactor(const LatticePosition& x) const
      {
        // Q = vLocal (0.5pi a**2)(a**2/(a**2 - r**2)
        // where r is the distance from the centreline
        // a is the radius of the circular iolet
        Dimensionless rFactor = (radius * radius) / (radius * radius - GetDistanceSquared(x));
        return 0.5 * PI * radius * radius * rFactor;
      }

      void InOutLetWK::DoPreStreamCoupling(const site_t& siteID,
                                           const LatticeTimeStep& timeStep,
                                           const LatticeVector& sitePos,
                                           const LatticeDensity& density,
                                           const LatticeVelocity& velocity)
      {
        if (siteID == centreSiteID)
				{
					LatticePressure pressure = GetPressure(0); // the argument is dummy
					distribn_t R0 = resistance, C0 = capacitance;

					// Explicit integration scheme
					//LatticePressure pressureNew = (1.0/C0)*flowRate + (1.0 - 1.0/(R0*C0))*pressure;

					// Semi-implicit integration scheme
					LatticePressure pressureNew = R0/(1.0 + R0*C0) * (flowRate + C0*pressure);

					densityNew = pressureNew / Cs2;
				}
        flowRateNew += velocity.Dot(-normal);
        siteCount ++;
      }

      void InOutLetWK::DoPostStreamCoupling(const site_t& siteID,
                                            const LatticeTimeStep& timeStep,
                                            const LatticeVector& sitePos)
      {
        if (siteID == centreSiteID)
        {
          density = densityNew;
        }
      }
    }
  }
}
