
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_IOLETS_INOUTLETFILEWK_H
#define HEMELB_LB_IOLETS_INOUTLETFILEWK_H
#include "lb/iolets/InOutLet.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      class InOutLetFileWK : public InOutLet
      {
        public:
          InOutLetFileWK();
          virtual ~InOutLetFileWK();
          virtual InOutLet* Clone() const;
 
	  LatticeDensity GetDensity(unsigned long time_step) const
	  {
		  return density;
	  }

	  void SetDensity(const LatticeDensity& d)
	  {
		  density = d;
	  }
 
          const std::string& GetFilePath()
          {
            return wkWeightsFilePath;
          }
          void SetFilePath(const std::string& path)
          {
            wkWeightsFilePath = path;
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
          
	  const distribn_t& GetArea() const
          {
            return area;
          }
          void SetArea(const distribn_t& r)
          {
            area = r;
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

          bool useWeightsFromFile;

          void Initialise(const util::UnitConverter* unitConverter);
        
	private:
          LatticeDistance radius;
          distribn_t rwk;
          distribn_t cwk;
          distribn_t area;
          LatticeDensity density;
          LatticeDensity densityMean;
          LatticeDensity densityAmp;

          std::map<std::vector<int>, double> weights_table;
          std::string wkWeightsFilePath;
          const util::UnitConverter* units;
      };
    }
  }
}
#endif // HEMELB_LB_IOLETS_INOUTLETWK_H
