
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "extraction/SurfaceWithinSphereSelector.h"

namespace hemelb
{
  namespace extraction
  {
    SurfaceWithinSphereSelector::SurfaceWithinSphereSelector(const util::Vector3D<float>& point,
                                                 float radius) :
        SphereGeometrySelector(point, radius)
    {
    }

    bool SurfaceWithinSphereSelector::IsWithinGeometry(const extraction::IterableDataSource& data,
                                                const util::Vector3D<site_t>& location)
    {
      return data.IsWallSite(location) & SphereGeometrySelector::IsWithinGeometry(data, location);
    }

  } /* namespace extraction */
} /* namespace hemelb */
