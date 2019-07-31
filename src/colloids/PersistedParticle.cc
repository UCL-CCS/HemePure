
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "colloids/PersistedParticle.h"
#include "constants.h"

namespace hemelb
{
	namespace colloids
	{
		PersistedParticle::PersistedParticle(io::xml::Element& xml)
		{
			// assume we are currently at a <subgridParticle> node
			xml.GetAttributeOrThrow("ParticleId", particleId);
			xml.GetAttributeOrThrow("Radius", smallRadius_a0);

			io::xml::Element initPosElem = xml.GetChildOrThrow("initialPosition");
			initPosElem.GetAttributeOrThrow("x", globalPosition.x);
			initPosElem.GetAttributeOrThrow("y", globalPosition.y);
			initPosElem.GetAttributeOrThrow("z", globalPosition.z);

			lastCheckpointTimestep = creationTimestep = 0;
			markedForDeletionTimestep = SITE_OR_BLOCK_SOLID;
		};
	}
}
