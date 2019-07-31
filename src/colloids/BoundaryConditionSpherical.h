
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_COLLOIDS_BOUNDARYCONDITIONSPHERICAL_H
#define HEMELB_COLLOIDS_BOUNDARYCONDITIONSPHERICAL_H

#include "colloids/BoundaryConditions.h"

namespace hemelb
{
	namespace colloids
	{
		class BoundaryConditionSpherical : public BoundaryCondition
		{
			public:
				static BoundaryCondition* ReadFromXml(io::xml::Element& xml)
				{
					LatticePosition sphereCentre;
					LatticeDistance sphereRadius;

					io::xml::Element radiusElem = xml.GetChildOrThrow("sphereRadius");
					configuration::GetDimensionalValue(radiusElem, "lattice", sphereRadius);

					io::xml::Element centreElem = xml.GetChildOrThrow("sphereCentre");
					centreElem.GetAttributeOrThrow("x", sphereCentre.x);
					centreElem.GetAttributeOrThrow("y", sphereCentre.y);
					centreElem.GetAttributeOrThrow("z", sphereCentre.z);

					return new BoundaryConditionSpherical(sphereCentre, sphereRadius);
				}

				virtual const bool DoSomethingToParticle(
						Particle& particle,
						const std::vector<LatticePosition> particleToWallVectors)
				{
					const bool IsInSphere = (
							pow(particle.GetGlobalPosition().x - sphereCentre.x,2) +
							pow(particle.GetGlobalPosition().y - sphereCentre.y,2) +
							pow(particle.GetGlobalPosition().z - sphereCentre.z,2) < pow(sphereRadius,2));

					return IsInSphere;
				}

			protected:
				BoundaryConditionSpherical(const LatticePosition sphereCentre, const LatticeDistance sphereRadius) :
					sphereCentre(sphereCentre), sphereRadius(sphereRadius) { };

				const LatticePosition sphereCentre;
				const LatticeDistance sphereRadius;
		};

		class BoundaryConditionFactorySpherical : public BoundaryConditionFactory<BoundaryConditionSpherical> { };
	}
}
#endif /* HEMELB_COLLOIDS_BOUNDARYCONDITIONSPHERICAL_H */
