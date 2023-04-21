
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_SPHEREGEOMETRYSELECTOR_H
#define HEMELB_EXTRACTION_SPHEREGEOMETRYSELECTOR_H

#include "extraction/GeometrySelector.h"

namespace hemelb
{
  namespace extraction
  {
    /**
     * Selects a geometry that forms a squat cylinder with height 0.5 lattice units
     * around a planar circle with specified normal and centre, and optionally specified
     * radius (assumed to be infinite when absent).
     */
    class SphereGeometrySelector : public GeometrySelector
    {
      public:
        /**
         * Constructor makes a sphere geometry object about a given point with given radius.
         * @param point
         * @param radius
         */
        SphereGeometrySelector(const util::Vector3D<float>& point, float radius);

        /**
         * Returns a point that lies on the sphere.
         * @return
         */
        const util::Vector3D<float>& GetPoint() const;

        /**
         * Returns the radius of the sphere.
         */
        float GetRadius() const;

      protected:
        /**
         * Returns true for any location within 0.5 lattice units of the sphere.
         *
         * @param data
         * @param location
         * @return
         */
        bool IsWithinGeometry(const extraction::IterableDataSource& data, const util::Vector3D<site_t>& location);

      private:
        /**
         * A point on the sphere.
         */
        const util::Vector3D<float> spherePoint;

        /**
         * The radius around the spherePoint to select. Radius <= 0 is taken as infinite.
         */
        const float radius;
    };
  }
}

#endif /* HEMELB_EXTRACTION_SPHEREGEOMETRYSELECTOR_H */
