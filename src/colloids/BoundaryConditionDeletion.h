
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_COLLOIDS_BOUNDARYCONDITIONDELETION_H
#define HEMELB_COLLOIDS_BOUNDARYCONDITIONDELETION_H

#include "colloids/BoundaryConditions.h"

namespace hemelb
{
	namespace colloids
	{
		class BoundaryConditionDeletion : public BoundaryCondition
		{
			public:
				static BoundaryCondition* ReadFromXml(io::xml::Element& xml)
				{
					LatticeDistance effectiveRange = -0.25;
					xml.GetAttributeOrThrow("effectiveRange", effectiveRange);

					return new BoundaryConditionDeletion(effectiveRange);
				}

				virtual const bool DoSomethingToParticle(
						Particle& particle,
						const std::vector<LatticePosition> particleToWallVectors)
				{
					// TODO: does not do *beyond* just *within* activation distance of boundary
					//LatticeDistance distance = wallNormal.GetMagnitudeSquared();
					log::Logger::Log<log::Trace, log::OnePerCore>(
							"*** In BoundaryConditionDeletion::DoSomethingToParticle for particleId: %lu ***\n",
							particle.GetParticleId());
					return false;//distance < (effectiveRange * effectiveRange);
				}

			protected:
				BoundaryConditionDeletion(LatticeDistance effectiveRange) : effectiveRange(effectiveRange) { };
				LatticeDistance effectiveRange;
		};

		class BoundaryConditionFactoryDeletion : public BoundaryConditionFactory<BoundaryConditionDeletion> { };
	}
}
#endif /* HEMELB_COLLOIDS_BOUNDARYCONDITIONDELETION_H */
