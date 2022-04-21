
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_LB_STREAMERS_GRINBERGKARNIADAKISWKIOLET_H
#define HEMELB_LB_STREAMERS_GRINBERGKARNIADAKISWKIOLET_H

#include "lb/streamers/StreamerTypeFactory.h"
#include "lb/streamers/GrinbergKarniadakisWKDelegate.h"
#include "lb/streamers/SimpleBounceBackDelegate.h"
#include "lb/streamers/BouzidiFirdaousLallemandDelegate.h"
#include "lb/streamers/JunkYangFactory.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<class CollisionType>
      struct GrinbergKarniadakisWKIolet
      {
          typedef IoletStreamerTypeFactory<CollisionType, GrinbergKarniadakisWKDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct GrinbergKarniadakisWKIoletSBB
      {
          typedef WallIoletStreamerTypeFactory<CollisionType, SimpleBounceBackDelegate<CollisionType> ,
              GrinbergKarniadakisWKDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct GrinbergKarniadakisWKIoletBFL
      {
          typedef WallIoletStreamerTypeFactory<CollisionType, BouzidiFirdaousLallemandDelegate<CollisionType> ,
              GrinbergKarniadakisWKDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct GrinbergKarniadakisWKIoletGZS
      {
          typedef WallIoletStreamerTypeFactory<CollisionType, GuoZhengShiDelegate<CollisionType> ,
              GrinbergKarniadakisWKDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct GrinbergKarniadakisWKIoletGZSE
      {
          typedef WallIoletStreamerTypeFactory<CollisionType, GuoZhengShiElasticWallDelegate<CollisionType> ,
              GrinbergKarniadakisWKDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct GrinbergKarniadakisWKIoletJY
      {
          typedef JunkYangFactory<CollisionType, GrinbergKarniadakisWKDelegate<CollisionType> > Type;
      };

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_GRINBERGKARNIADAKISWKIOLET_H */
