// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "SimulationMaster.h"
#include "configuration/SimConfig.h"
#include "extraction/PropertyActor.h"
#include "extraction/LbDataSourceIterator.h"
#include "io/writers/xdr/XdrFileWriter.h"
#include "util/utilityFunctions.h"
#include "geometry/GeometryReader.h"
#include "geometry/LatticeData.h"
#include "util/fileutils.h"
#include "log/Logger.h"
#include "lb/HFunction.h"
#include "io/xml/XmlAbstractionLayer.h"
#include "colloids/ColloidController.h"
#include "net/BuildInfo.h"
#include "net/IOCommunicator.h"
#include "colloids/BodyForces.h"
#include "colloids/BoundaryConditions.h"

#include <map>
#include <limits>
#include <cstdlib>

/**
 * Constructor for the SimulationMaster class
 *
 * Initialises member variables including the network topology
 * object.
 */

double hemelb::tau_;

SimulationMaster::SimulationMaster(hemelb::configuration::CommandLine & options, const hemelb::net::IOCommunicator& ioComm) :
	ioComms(ioComm), timings(ioComm), build_info(), communicationNet(ioComm) {
	timings[hemelb::reporting::Timers::total].Start();

	latticeData = NULL;

	colloidController = NULL;
	latticeBoltzmannModel = NULL;
	propertyDataSource = NULL;
	propertyExtractor = NULL;
	simulationState = NULL;
	stepManager = NULL;
	netConcern = NULL;
	neighbouringDataManager = NULL;
	imagesPerSimulation = options.NumberOfImages();

	fileManager = new hemelb::io::PathManager(options, IsCurrentProcTheIOProc(), GetProcessorCount());
	simConfig = hemelb::configuration::SimConfig::New(fileManager->GetInputFile());
	unitConverter = &simConfig->GetUnitConverter();
	monitoringConfig = simConfig->GetMonitoringConfiguration();

	fileManager->SaveConfiguration(simConfig);
	Initialise();

	if (IsCurrentProcTheIOProc()) {
		reporter = new hemelb::reporting::Reporter(
			fileManager->GetReportPath(),
			fileManager->GetInputFile());
		reporter->AddReportable(&build_info);

		if (monitoringConfig->doIncompressibilityCheck) {
			reporter->AddReportable(incompressibilityChecker);
		}

		reporter->AddReportable(&timings);
		reporter->AddReportable(latticeData);
		reporter->AddReportable(simulationState);
	}
}

/**
 * Destructor for the SimulationMaster class.
 *
 * Deallocates dynamically allocated memory to contained classes.
 */
SimulationMaster::~SimulationMaster() {

	delete latticeData;
	delete colloidController;
	delete latticeBoltzmannModel;
	delete inletValues;
	delete outletValues;
	delete propertyExtractor;
	delete propertyDataSource;
	delete stabilityTester;
	delete entropyTester;
	delete simulationState;
	delete incompressibilityChecker;
	delete neighbouringDataManager;

	delete simConfig;
	delete fileManager;

	if (IsCurrentProcTheIOProc()) {
		delete reporter;
	}

	delete stepManager;
	delete netConcern;
}

/**
 * Returns true if the current processor is the dedicated I/O
 * processor.
 */
bool SimulationMaster::IsCurrentProcTheIOProc() {
	return ioComms.OnIORank();
}

/**
 * Returns the number of processors involved in the simulation.
 */
int SimulationMaster::GetProcessorCount() {
	return ioComms.Size();
}

/**
 * Initialises various elements of the simulation
 */
void SimulationMaster::Initialise() {

	simulationState = new hemelb::lb::SimulationState(
		simConfig->GetTimeStepLength(),
		simConfig->GetTotalTimeSteps());

	timings[hemelb::reporting::Timers::latDatInitialise].Start();

	// Use a reader to read in the file.
	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("INITIALISE");
	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("----------");
	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("--> loading input and decomposing geometry");
	hemelb::geometry::GeometryReader reader(
		latticeType::GetLatticeInfo(),
		timings, ioComms);
	hemelb::geometry::Geometry readGeometryData =
		reader.LoadAndDecompose(simConfig->GetDataFilePath());

	// Create a new lattice based on that info and return it.
	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("--> lattice data");
	latticeData = new hemelb::geometry::LatticeData(latticeType::GetLatticeInfo(),
			readGeometryData,
			ioComms);

	timings[hemelb::reporting::Timers::latDatInitialise].Stop();

	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("--> neighbouring data manager");
	neighbouringDataManager = new hemelb::geometry::neighbouring::NeighbouringDataManager(*latticeData,
				latticeData->GetNeighbouringData(),
				communicationNet);

	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("--> lattice-Boltzmann model");
	latticeBoltzmannModel = new hemelb::lb::LBM<latticeType>(simConfig,
			&communicationNet,
			latticeData,
			simulationState,
			timings,
			neighbouringDataManager);

	hemelb::lb::MacroscopicPropertyCache& propertyCache = latticeBoltzmannModel->GetPropertyCache();

	hemelb::tau_ = latticeBoltzmannModel->GetLbmParams()->GetTau();

	if (simConfig->HasColloidSection()) {
		hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("--> colloid section present");

		timings[hemelb::reporting::Timers::colloidInitialisation].Start();
		hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("----> loading colloid configuration");
		std::string colloidConfigPath = simConfig->GetColloidConfigPath();
		hemelb::io::xml::Document xml(colloidConfigPath);

		hemelb::colloids::BodyForces::InitBodyForces(xml);
		hemelb::colloids::BoundaryConditions::InitBoundaryConditions(latticeData, xml);

		hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("----> initialising colloid controller");
		colloidController =
			new hemelb::colloids::ColloidController(*latticeData,
					*simulationState,
					simConfig,
					readGeometryData,
					xml,
					propertyCache,
					//latticeBoltzmannModel->GetLbmParams(),
					fileManager->GetColloidPath(),
					ioComms,
					timings);
		timings[hemelb::reporting::Timers::colloidInitialisation].Stop();
	}

	stabilityTester = new hemelb::lb::StabilityTester<latticeType>(latticeData,
			&communicationNet,
			simulationState,
			timings,
			monitoringConfig);
	entropyTester = NULL;

	if (monitoringConfig->doIncompressibilityCheck) {
		incompressibilityChecker = new hemelb::lb::IncompressibilityChecker <
		hemelb::net::PhasedBroadcastRegular<> > (latticeData,
				&communicationNet,
				simulationState,
				latticeBoltzmannModel->GetPropertyCache(),
				timings);
	} else {
		incompressibilityChecker = NULL;
	}

	inletValues = new hemelb::lb::iolets::BoundaryValues(hemelb::geometry::INLET_TYPE,
			latticeData,
			simConfig->GetInlets(),
			simulationState,
			ioComms,
			*unitConverter);

	outletValues = new hemelb::lb::iolets::BoundaryValues(hemelb::geometry::OUTLET_TYPE,
			latticeData,
			simConfig->GetOutlets(),
			simulationState,
			ioComms,
			*unitConverter);

	latticeBoltzmannModel->Initialise(inletValues, outletValues, unitConverter);
        latticeBoltzmannModel->SetInitialConditions(ioComms); //JM Checkpoint addition
	neighbouringDataManager->ShareNeeds();
	neighbouringDataManager->TransferNonFieldDependentInformation();

	propertyDataSource =
		new hemelb::extraction::LbDataSourceIterator(latticeBoltzmannModel->GetPropertyCache(),
				*latticeData,
				ioComms.Rank(),
				*unitConverter);

	if (simConfig->PropertyOutputCount() > 0) {

		for (unsigned outputNumber = 0; outputNumber < simConfig->PropertyOutputCount(); ++outputNumber) {
			simConfig->GetPropertyOutput(outputNumber)->filename = fileManager->GetDataExtractionPath()
					+ simConfig->GetPropertyOutput(outputNumber)->filename;
		}

		propertyExtractor = new hemelb::extraction::PropertyActor(*simulationState,
				simConfig->GetPropertyOutputs(),
				*propertyDataSource,
				timings, ioComms);
	}

	imagesPeriod = OutputPeriod(imagesPerSimulation);

	stepManager = new hemelb::net::phased::StepManager(2,
			&timings,
			hemelb::net::separate_communications);
	netConcern = new hemelb::net::phased::NetConcern(communicationNet);
	stepManager->RegisterIteratedActorSteps(*neighbouringDataManager, 0);

	if (colloidController != NULL) {
		stepManager->RegisterIteratedActorSteps(*colloidController, 1);
	}

	stepManager->RegisterIteratedActorSteps(*latticeBoltzmannModel, 1);

	stepManager->RegisterIteratedActorSteps(*inletValues, 1);
	stepManager->RegisterIteratedActorSteps(*outletValues, 1);
	stepManager->RegisterIteratedActorSteps(*stabilityTester, 1);

	if (entropyTester != NULL) {
		stepManager->RegisterIteratedActorSteps(*entropyTester, 1);
	}

	if (monitoringConfig->doIncompressibilityCheck) {
		stepManager->RegisterIteratedActorSteps(*incompressibilityChecker, 1);
	}

	if (propertyExtractor != NULL) {
		stepManager->RegisterIteratedActorSteps(*propertyExtractor, 1);
	}

	stepManager->RegisterCommsForAllPhases(*netConcern);

	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("-------------------");
	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("INITIALISE FINISHED");
}

unsigned int SimulationMaster::OutputPeriod(unsigned int frequency) {
	if (frequency == 0) {
		return 1000000000;
	}

	unsigned long roundedPeriod = simulationState->GetTotalTimeSteps() / frequency;
	return hemelb::util::NumericalFunctions::max(1U, (unsigned int) roundedPeriod);
}

void SimulationMaster::HandleActors() {
	stepManager->CallActions();
}

void SimulationMaster::OnUnstableSimulation() {
	LogStabilityReport();
	hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("ABORTING :: time step length: %f",
			simulationState->GetTimeStepLength());
	Finalise();
	Abort();
}

/**
 * Begin the simulation.
 */
void SimulationMaster::RunSimulation() {
	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("SIMULATION STARTING");
	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("-------------------");
	timings[hemelb::reporting::Timers::simulation].Start();

	while (simulationState->GetTimeStep() <= simulationState->GetTotalTimeSteps()) {
		DoTimeStep();

		if (simulationState->IsTerminating()) {
			break;
		}
	}

	timings[hemelb::reporting::Timers::simulation].Stop();
	Finalise();
}

void SimulationMaster::Finalise() {
	timings[hemelb::reporting::Timers::total].Stop();
	timings.Reduce();

	if (IsCurrentProcTheIOProc()) {
		reporter->FillDictionary();
		reporter->Write();
	}

	// DTMP: Logging output on communication as debug output for now.
	hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("sync points: %lld, bytes sent: %lld",
			communicationNet.SyncPointsCounted,
			communicationNet.BytesSent);

	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("-------------------");
	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("SIMULATION FINISHED");
}

void SimulationMaster::DoTimeStep() {
	bool writeImage = ((simulationState->GetTimeStep() % imagesPeriod) == 0) ?
					true :
					false;

	if (simulationState->GetTimeStep() % 200 == 0) {
		hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("time step %07i :: write_image_to_disk %i",
				simulationState->GetTimeStep(),
				writeImage);
		LogStabilityReport();
	}

	RecalculatePropertyRequirements();

	HandleActors();

	if (simulationState->GetStability() == hemelb::lb::Unstable) {
		OnUnstableSimulation();
	}

	// If the user requested to terminate converged steady flow simulations, mark
	// simulation to be finished at the end of the current timestep.
	if ((simulationState->GetStability() == hemelb::lb::StableAndConverged)
			&& monitoringConfig->convergenceTerminate) {
		LogStabilityReport();
		simulationState->SetIsTerminating(true);
	}

	// Colloid output
	//if ((simulationState->GetTimeStep() % 100 == 0) && colloidController != NULL)
	//	colloidController->OutputInformation(simulationState->GetTimeStep());

	if (simulationState->GetTimeStep() % FORCE_FLUSH_PERIOD == 0 && IsCurrentProcTheIOProc()) {
		fflush(NULL);
	}

	simulationState->Increment();
}

void SimulationMaster::RecalculatePropertyRequirements() {
	// Get the property cache & reset its list of properties to get.
	hemelb::lb::MacroscopicPropertyCache& propertyCache = latticeBoltzmannModel->GetPropertyCache();

	propertyCache.ResetRequirements();

	if (monitoringConfig->doIncompressibilityCheck) {
		propertyCache.densityCache.SetRefreshFlag();
		propertyCache.velocityCache.SetRefreshFlag();
	}

	// If extracting property results, check what's required by them.
	if (propertyExtractor != NULL) {
		propertyExtractor->SetRequiredProperties(propertyCache);
	}
}

/**
 * Called on error to abort the simulation and pull-down the MPI environment.
 */
void SimulationMaster::Abort() {
	// This gives us something to work from when we have an error - we get the rank
	// that calls abort, and we get a stack-trace from the exception having been thrown.
	hemelb::log::Logger::Log<hemelb::log::Critical, hemelb::log::Singleton>("ABORTING");
	hemelb::net::MpiEnvironment::Abort(1);

	exit(1);
}

void SimulationMaster::LogStabilityReport() {
	if (monitoringConfig->doIncompressibilityCheck
			&& incompressibilityChecker->AreDensitiesAvailable()) {
		hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("time step %07i :: tau: %.6f, max_relative_press_diff: %.3f, Ma: %.3f, max_vel_phys: %e",
				simulationState->GetTimeStep(),
				latticeBoltzmannModel->GetLbmParams()->GetTau(),
				incompressibilityChecker->GetMaxRelativeDensityDifference(),
				incompressibilityChecker->GetGlobalLargestVelocityMagnitude()
				/ hemelb::Cs,
				unitConverter->ConvertVelocityToPhysicalUnits(incompressibilityChecker->GetGlobalLargestVelocityMagnitude()));
	}

	if (simulationState->GetStability() == hemelb::lb::StableAndConverged) {
		hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("time step %07i :: steady flow simulation converged",
				simulationState->GetTimeStep());
	}
}
