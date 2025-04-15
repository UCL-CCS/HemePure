
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_IOLETS_INOUTLETWK2_H
#define HEMELB_LB_IOLETS_INOUTLETWK2_H
#include "lb/iolets/InOutLet.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      class InOutLetWK2 : public InOutLet
      {
        public:
          InOutLetWK2();
          virtual ~InOutLetWK2();
          virtual InOutLet* Clone() const;

          bool IsCommsRequired() const
          {
            return true;
          }

          void DoComms(const BoundaryCommunicator& boundaryComm, const LatticeTimeStep timeStep);

          virtual LatticeDensity GetDensityMin() const
          {
            return 1.0;
          }

          virtual LatticeDensity GetDensityMax() const
          {
            return 1.0;
          }

          LatticeDensity GetDensity(LatticeTimeStep time_step) const
          {
            return density;
          }

          void SetDensity(const LatticeDensity& d)
          {
            density = d;
          }
          
          LatticeDensity GetDensityNew(LatticeTimeStep time_step) const
          {
            return densityNew;
          }

          void SetDensityNew(const LatticeDensity& d)
          {
            densityNew = d;
          }

          virtual void Reset(SimulationState &state)
          {
            //pass;
          }
          
	        /**
           * Note that the radius and area for these are specified in LATTICE UNITS in the XML file.
           * This is indeed a horrible hack.
           * @return
           */
          const LatticeDistance& GetRadius() const
          {
            return radius;
          }
          
          void SetRadius(const LatticeDistance& r)
          {
            radius = r;
          }

          const distribn_t& GetArea() const
          {
            return area;
          }

          void SetArea(const distribn_t& a)
          {
            area = a;
          }
        
          virtual distribn_t GetScaleFactor(const LatticePosition& x) const;
          
          LatticeDistance GetDistanceSquared(const LatticePosition& x) const;
          
          const distribn_t& GetResistance() const
          {
            return resistance;
          }

          void SetResistance(const distribn_t& R)
          {
            resistance = R;
          }

          const distribn_t& GetCapacitance() const
          {
            return capacitance;
          }

          void SetCapacitance(const distribn_t& C)
          {
            capacitance = C;
          }

          void DoPreStreamCoupling(const site_t& siteID,
                                   const LatticeTimeStep& timeStep,
                                   const LatticeVector& sitePos,
                                   const LatticeDensity& density,
                                   const LatticeVelocity& velocity);

          void DoPostStreamCoupling(const site_t& siteID,
                                    const LatticeTimeStep& timeStep,
                                    const LatticeVector& sitePos);

        protected:
          LatticeDistance radius;
          distribn_t area, resistance, capacitance;
          LatticeDensity density, densityNew;
          distribn_t flowRate, flowRateNew;
          site_t siteCount;
      };
    }
  }
}
#endif // HEMELB_LB_IOLETS_INOUTLETWK2_H
