
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "lb/iolets/InOutLetWK3.h"
#include "lb/iolets/BoundaryComms.h"
#include "configuration/SimConfig.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      InOutLetWK3::InOutLetWK3() :
          InOutLet(), radius(1.0), area(1.0), resistance_c(1.0), resistance_p(1.0), capacitance(1.0),
          density(1.0), densityNew(1.0), flowRateOld(0.0), flowRate(0.0), flowRateNew(0.0),
          siteCount(0)
      {
      }

      InOutLet* InOutLetWK3::Clone() const
      {
        InOutLet* copy = new InOutLetWK3(*this);
        return copy;
      }

      InOutLetWK3::~InOutLetWK3()
      {
      }

      void InOutLetWK3::DoComms(const BoundaryCommunicator& boundaryComm, const LatticeTimeStep timeStep)
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

      LatticeDistance InOutLetWK3::GetDistanceSquared(const LatticePosition& x) const
      {
        LatticePosition displ = x - position;
        LatticeDistance z = displ.Dot(normal);
        return displ.GetMagnitudeSquared() - z * z;
      }

      distribn_t InOutLetWK3::GetScaleFactor(const LatticePosition& x) const
      {
        // Q = vLocal (0.5pi a**2)(a**2/(a**2 - r**2)
        // where r is the distance from the centreline
        // a is the radius of the circular iolet
        Dimensionless rFactor = (radius * radius) / (radius * radius - GetDistanceSquared(x));
        return 0.5 * PI * radius * radius * rFactor;
      }

      void InOutLetWK3::DoPreStreamCoupling(const site_t& siteID,
                                           const LatticeTimeStep& timeStep,
                                           const LatticeVector& sitePos,
                                           const LatticeDensity& density,
                                           const LatticeVelocity& velocity)
      {
        if (siteID == centreSiteID)
				{
					LatticePressure pressure = GetPressure(0); // the argument is dummy
					const distribn_t Rc = resistance_c, Rp = resistance_p, Cp = capacitance;

          /* To solve for P in the ODE: P + Rp*Cp*dP/dt = (Rc+Rp)*Q + Rc*Rp*Cp*dQ/dt
           * - P is approximated as P^{n+1}
           * - Q^{n} is not available, so Q is approximated as Q^{n} = 2*Q^{n-1} - Q^{n-2}
           * - Similarly, dQ/dt|^{n} = (Q^{n-1}-Q^{n-2})/dt
           * - (1+Rp*Cp/dt)*P^{n+1} = (Rp*Cp/dt)*P^{n} + (Rc+Rp)*Q^{n} + (Rc*Rp*Cp/dt)*(Q^{n}-Q^{n-1})
           * - dt = 1 in lattice unit
           */
          LatticePressure pressureNew = (Rp*Cp*pressure + (Rc+Rp)*(2.0*flowRate - flowRateOld) 
                                          + Rc*Rp*Cp*(flowRate - flowRateOld)) / (1.0 + Rp*Cp);

					densityNew = pressureNew / Cs2;
				}
        flowRateNew += velocity.Dot(-normal);
        siteCount ++;
      }

      void InOutLetWK3::DoPostStreamCoupling(const site_t& siteID,
                                            const LatticeTimeStep& timeStep,
                                            const LatticeVector& sitePos)
      {
        if (siteID == centreSiteID)
        {
          density = densityNew;
          flowRateOld = flowRate;
        }
      }
    }
  }
}
