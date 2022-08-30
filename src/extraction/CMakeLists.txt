# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

add_library(hemelb_extraction
	GeometrySelector.cc
	StraightLineGeometrySelector.cc
	LocalPropertyOutput.cc
	IterableDataSource.cc
	PlaneGeometrySelector.cc
	SphereGeometrySelector.cc
	PropertyActor.cc
	PropertyWriter.cc
	InletSelector.cc
	OutletSelector.cc
	WholeGeometrySelector.cc
	LbDataSourceIterator.cc
	GeometrySurfaceSelector.cc
	SurfaceWithinSphereSelector.cc
	SurfacePointSelector.cc
        LocalDistributionInput.cc)
