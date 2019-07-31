
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <unistd.h>

#include "extraction/InletSelector.h"

namespace hemelb
{
	namespace extraction
	{
		bool InletSelector::IsWithinGeometry(
				const extraction::IterableDataSource& data,
				const util::Vector3D<site_t>& location)
		{
			//if (site.GetSiteType() == hemelb::geometry::INLET_TYPE)
			//	log::Logger::Log<log::Info, log::OnePerCore>("site %i has pressure", site);

			if (data.IsInletSite(location))
			{
				//{
				//	int i = 0;
				//	char hostname[256];
				//	gethostname(hostname, sizeof(hostname));
				//	printf("PID %d (%d) on %s ready for attach\n", getpid(), 1, hostname);
				//	fflush(stdout);
				//	while (0 == i)
				//		sleep(5);
				//}
				return true;
			}
		}
	}
}
