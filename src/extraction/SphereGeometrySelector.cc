
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "extraction/SphereGeometrySelector.h"

namespace hemelb
{
  namespace extraction
  {
    SphereGeometrySelector::SphereGeometrySelector(const util::Vector3D<float>& point,
                                                 float radius) :
      spherePoint(point), radius(radius)
    {

    }

    const util::Vector3D<float>& SphereGeometrySelector::GetPoint() const
    {
      return spherePoint;
    }

    float SphereGeometrySelector::GetRadius() const
    {
      return radius;
    }

    bool SphereGeometrySelector::IsWithinGeometry(const extraction::IterableDataSource& data,
                                                 const util::Vector3D<site_t>& location)
    {
      util::Vector3D<float> coords = util::Vector3D<float>(location) * data.GetVoxelSize() + data.GetOrigin();
      const float distance = (coords - spherePoint).GetMagnitude();
      return distance <= radius;
    }
  }
}
