
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_EXTRACTION_SURFACEWITHINSPHERESELECTOR_H
#define HEMELB_EXTRACTION_SURFACEWITHINSPHERESELECTOR_H

#include "extraction/SphereGeometrySelector.h"

namespace hemelb
{
  namespace extraction
  {

    /**
     * This class defines a geometry selector for a point in the geometry surface. By geometry surface we mean the
     * original surface passed into the setup tool (a .stl file possibly). The class will select any lattice site
     * known to be a wall (in LatticeData terms) within a sphere.
     */
    class SurfaceWithinSphereSelector : public SphereGeometrySelector
    {
      public:
        /**
         * Constructor makes a sphere geometry object about a given point with given radius.
         * @param point
         * @param radius
         */
        SurfaceWithinSphereSelector(const util::Vector3D<float>& point, float radius);

      protected:
        /**
         * Returns true if the given location is within the selection.
         * @param data data source abstraction
         * @param location lattice site coordinates to evaluate the selector on
         * @return whether location is within the selection
         */
        bool IsWithinGeometry(const extraction::IterableDataSource& data, const util::Vector3D<site_t>& location);
    };

  } /* namespace extraction */
} /* namespace hemelb */
#endif /* HEMELB_EXTRACTION_SURFACEWITHINSPHERESELECTOR_H */
