
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_IOLETS_INOUTLETWK_H
#define HEMELB_LB_IOLETS_INOUTLETWK_H
#include "lb/iolets/InOutLet.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      class InOutLetWK : public InOutLet
      {
        public:
          InOutLetWK();
          virtual ~InOutLetWK();
          virtual InOutLet* Clone() const;

          bool IsCommsRequired() const
          {
            return true;
          }

          void DoComms(const BoundaryCommunicator& boundaryComm, const LatticeTimeStep timeStep);

          virtual LatticeDensity GetDensityMin() const
          {
            //pass;
          }

          virtual LatticeDensity GetDensityMax() const
          {
            //pass;
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
           * Note that the radius and max speed for these are specified in LATTICE UNITS in the XML file.
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
        
          virtual distribn_t GetScaleFactor(const LatticePosition& x) const;
          
          LatticeDistance GetDistance(const LatticePosition& x) const;
          
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

        protected:
          LatticeDensity density, densityNew;
          LatticeDistance radius;
          distribn_t resistance;
          distribn_t capacitance;
      };
    }
  }
}
#endif // HEMELB_LB_IOLETS_INOUTLETWK_H
