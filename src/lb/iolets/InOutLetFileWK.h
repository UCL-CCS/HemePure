
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_IOLETS_INOUTLETFILEWK_H
#define HEMELB_LB_IOLETS_INOUTLETFILEWK_H
#include "lb/iolets/InOutLetWK2.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      class InOutLetFileWK : public InOutLetWK2
      {
        public:
          InOutLetFileWK();
          virtual ~InOutLetFileWK();
          virtual InOutLet* Clone() const;
          
          virtual distribn_t GetScaleFactor(const LatticePosition& x) const;

          const std::string& GetFilePath()
          {
            return wkWeightsFilePath;
          }

          void SetFilePath(const std::string& path)
          {
            wkWeightsFilePath = path;
          }
          
          bool useWeightsFromFile;

          void Initialise(const util::UnitConverter* unitConverter);
        
        private:
          std::map<std::vector<int>, double> weights_table;
          std::string wkWeightsFilePath;
          const util::UnitConverter* units;
      };
    }
  }
}
#endif // HEMELB_LB_IOLETS_INOUTLETWK_H
