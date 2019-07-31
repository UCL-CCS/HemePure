
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_FORMATS_COLLOIDS_H
#define HEMELB_IO_FORMATS_COLLOIDS_H

#include "io/formats/formats.h"

namespace hemelb
{
	namespace io
	{
		namespace formats
		{
			namespace colloids
			{
				enum
				{
					/* Identify colloid files
					 * ASCII for 'col', then EOF
					 * Combined magic number is
					 * hex:    68 6c 62 21 63 6f 6c 04
					 * ascii:  h  l  b  !  c  o  l EOF
					 */
					MagicNumber = 0x636f6c04
				};
				enum
				{
					VersionNumber = 1
				};
				enum
				{
					MagicLength = 12
				};
				enum
				{
					HeaderLength = 24
				};
				enum
				{
					RecordLength = 88
				  //RecordLength = 64
				  //RecordLength = 40
				};
			}
		}
	}
}
#endif // HEMELB_IO_FORMATS_COLLOIDS_H
