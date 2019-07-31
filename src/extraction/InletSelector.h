
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_INLETSELECTOR_H
#define HEMELB_EXTRACTION_INLETSELECTOR_H

#include "extraction/GeometrySelector.h"

namespace hemelb
{
  namespace extraction
  {
    /**
     * Class for selecting inlets.
     */
    class InletSelector : public GeometrySelector
    {
      protected:
        /**
         * Returns true for all locations.
         *
         * @param data
         * @param location
         * @return
         */
        bool IsWithinGeometry(const extraction::IterableDataSource& data, const util::Vector3D<site_t>& location);
    };
  }
}

#endif /* HEMELB_EXTRACTION_INLETSELECTOR_H */
