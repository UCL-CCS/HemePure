
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_COLLOIDS_BOUNDARYCONDITIONLUBRICATION_H
#define HEMELB_COLLOIDS_BOUNDARYCONDITIONLUBRICATION_H

#include "colloids/BoundaryConditions.h"
#include "io/xml/XmlAbstractionLayer.h"

namespace hemelb
{
	namespace colloids
	{
		class BoundaryConditionLubrication : public BoundaryCondition
		{
			public:
				static BoundaryCondition* ReadFromXml(io::xml::Element& xml)
				{
					LatticeDistance effectiveRange = 1.0;
					xml.GetAttributeOrThrow("effectiveRange", effectiveRange);
					return new BoundaryConditionLubrication(effectiveRange);
				}

				virtual const bool DoSomethingToParticle(
						Particle& particle,
						const std::vector<LatticePosition> particleToWallVectors)
				{
					const bool keep = true;

					log::Logger::Log<log::Trace, log::OnePerCore>(
							"*** In BoundaryConditionLubrication::DoSomethingToParticle for particleId: %lu ***\n",
							particle.GetParticleId());

					LatticeVelocity velocity;
					velocity *= 0.0;

					for (std::vector<LatticePosition>::const_iterator iter = particleToWallVectors.begin();
							iter != particleToWallVectors.end();
							iter++)
					{
						const LatticePosition particleToWallVector = *iter;
						const LatticeDistance separation_h = particleToWallVector.GetMagnitude()
							- particle.GetRadius(false);

						log::Logger::Log<log::Trace, log::OnePerCore>(
								"*** In BoundaryConditionLubrication::DoSomethingToParticle - wall vector: {%g,%g,%g}, mag: %g, particle radius: %g, separation_h: %g\n",
								particleToWallVector.x,
								particleToWallVector.y,
								particleToWallVector.z,
								particleToWallVector.GetMagnitude(),
								particle.GetRadius(false),
								separation_h);

						// from tenCate2002, Eq. (10)
						if (separation_h <= effectiveRange) {
							velocity += -(
									 particleToWallVector.GetNormalised()
									*particleToWallVector.GetNormalised().Dot(particle.GetVelocity()))
									*((separation_h-effectiveRange)/(separation_h*effectiveRange))
									*particle.GetRadius(false);

							// used to plot lubrication force in 'cube' case
							//LatticeVelocity velocity_p = particle.GetVelocity();
							//LatticeForceVector dimless;
							//dimless.x = velocity.x/velocity_p.x+1.0;
							//dimless.y = velocity.y/velocity_p.y+1.0;
							//dimless.z = velocity.z/velocity_p.z+1.0;
							//log::Logger::Log<log::Info, log::OnePerCore>("{%g %g}", separation_h/particle.GetRadius(false), dimless.z);
						}
					} particle.UpdateVelocity(velocity);
					return keep;
				}

			protected:
				BoundaryConditionLubrication(const LatticeDistance effectiveRange) : effectiveRange(effectiveRange) {};
				LatticeDistance effectiveRange;
		};

		class BoundaryConditionFactoryLubrication : public BoundaryConditionFactory<BoundaryConditionLubrication> { };
	}
}
#endif /* HEMELB_COLLOIDS_BOUNDARYCONDITIONLUBRICATION_H */
