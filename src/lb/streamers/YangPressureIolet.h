
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_LB_STREAMERS_YANGPRESSUREIOLET_H
#define HEMELB_LB_STREAMERS_YANGPRESSUREIOLET_H

#include "lb/streamers/StreamerTypeFactory.h"
#include "lb/streamers/YangPressureDelegate.h"
#include "lb/streamers/SimpleBounceBackDelegate.h"
#include "lb/streamers/BouzidiFirdaousLallemandDelegate.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<class CollisionType>
      struct YangPressureIolet
      {
          typedef IoletStreamerTypeFactory<CollisionType, YangPressureDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct YangPressureIoletSBB
      {
          typedef WallIoletStreamerTypeFactory<CollisionType, SimpleBounceBackDelegate<CollisionType> ,
              YangPressureDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct YangPressureIoletBFL
      {
          typedef WallIoletStreamerTypeFactory<CollisionType, BouzidiFirdaousLallemandDelegate<CollisionType> ,
              YangPressureDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct YangPressureIoletGZS
      {
          typedef WallIoletStreamerTypeFactory<CollisionType, GuoZhengShiDelegate<CollisionType> ,
              YangPressureDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct YangPressureIoletGZSE
      {
          typedef WallIoletStreamerTypeFactory<CollisionType, GuoZhengShiElasticWallDelegate<CollisionType> ,
              YangPressureDelegate<CollisionType> > Type;
      };

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_YANGPRESSUREIOLET_H */
