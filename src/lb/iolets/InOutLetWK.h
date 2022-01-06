
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

	  LatticeDensity GetDensity(unsigned long time_step) const
	  {
		  return density;
	  }

	  void SetDensity(const LatticeDensity& d)
	  {
		  density = d;
	  }
	  
	  LatticeDensity GetDensityNew(unsigned long time_step) const
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

          const LatticeDensity& GetDensityMean() const
          {
            return densityMean;
          }
          void SetDensityMean(const LatticeDensity& rho)
          {
            densityMean = rho;
          }

          const LatticeDensity& GetDensityAmp() const
          {
            return densityAmp;
          }
          void SetDensityAmp(const LatticeDensity& rho)
          {
            densityAmp = rho;
          }
          
	  LatticeDensity GetDensityMin() const
          {
            return (densityMean - densityAmp);
          }
          LatticeDensity GetDensityMax() const
          {
            return (densityMean + densityAmp);
          }
        
	  distribn_t GetQtScaleFactor(const LatticePosition& x) const;
	  
	  distribn_t GetDistance(const LatticePosition& x) const;
	  
	  const distribn_t& GetRwk() const
	  {
		  return rwk;
          }
	  void SetRwk(const distribn_t r)
	  {
		  rwk = r;
	  }

	  const distribn_t& GetCwk() const
	  {
		  return cwk;
          }
	  void SetCwk(const distribn_t c)
	  {
		  cwk = c;
	  }

        private:
          LatticeDistance radius;
          distribn_t rwk;
          distribn_t cwk;
          LatticeDensity density;
          LatticeDensity densityNew;
          LatticeDensity densityMean;
          LatticeDensity densityAmp;

      };
    }
  }
}
#endif // HEMELB_LB_IOLETS_INOUTLETWK_H
