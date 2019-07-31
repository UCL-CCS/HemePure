
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cstdarg>
#include <sys/time.h>
#include <sys/resource.h>
#include <string.h>

#include "util/utilityFunctions.h"
#include "net/mpi.h"
#include "log/Logger.h"

namespace hemelb
{
	namespace log
	{
		const LogLevel Logger::currentLogLevel = HEMELB_LOG_LEVEL;
		// Use negative value to indicate uninitialised.
		int Logger::thisRank = -1;
		double Logger::startTime = -1.0;

		void Logger::Init()
		{	// If Logger uninitialised
			if (thisRank < 0)
			{
				// Check that MPI is ready
				if (net::MpiEnvironment::Initialized())
				{
					thisRank = net::MpiCommunicator::World().Rank();
				}
				startTime = util::myClock();
			}
		}

		template<>
			void Logger::LogInternal<OnePerCore>(std::string format, std::va_list args)
			{
				std::stringstream output;
				// Set the fill digit to be 0, so the integer 1 renders as 0000001
				output.fill('0');
				output << "[Rank " << std::setw(7) << thisRank << ", " << std::setiosflags(std::ios::scientific)
					<< std::setw(8) << (util::myClock() - startTime) << " s";

#ifdef HAVE_RUSAGE
				rusage usage;
				getrusage(RUSAGE_SELF, &usage);
				output << ", " << std::setw(14) << usage.ru_maxrss;
#endif
				output << " kB] :: " << format << '\n';

				std::string overFormat(output.str());
				std::vprintf(overFormat.c_str(), args);
			}

		template<>
			void Logger::LogInternal<Singleton>(std::string format, std::va_list args)
			{
				// Get the the current process' status file from the proc filesystem
				FILE* procfile = fopen("/proc/self/status", "r");

				long to_read = 8192; char buffer[to_read];
				int read = fread(buffer, sizeof(char), to_read, procfile);
				fclose(procfile);

				bool found_vmrss = false;
				char* search_result;
				long np = net::MpiCommunicator::World().Size();
				long vmrss_kb;
				long vmrss_kb_pp[np];

				// Look through proc status contents line by line
				char delims[] = "\n";
				char* line = strtok(buffer, delims);
				while (line != NULL && !found_vmrss) {
					search_result = strstr(line, "VmRSS:");
					if (search_result != NULL) {
						sscanf(line, "%*s %ld", &vmrss_kb);
						found_vmrss = true;
					} line = strtok(NULL, delims);
				}

				// Gather all measurements and send to rank 0
				MPI_Gather(
						&vmrss_kb,    1, MPI_UNSIGNED_LONG,
						&vmrss_kb_pp, 1, MPI_UNSIGNED_LONG,
						0, net::MpiCommunicator::World());

				if (thisRank == 0) {
					// Only print if rank 0
					long vmrss_kb_global = 0;
					for (int i = 0; i < np; i++)
						vmrss_kb_global += vmrss_kb_pp[i];

					std::stringstream output;
					// Set the fill digit to be 0, so the integer 1 renders as 0000001
					output.fill('0');
					output << "[Rank " << std::setw(7) << thisRank << ", " << std::setiosflags(std::ios::scientific)
						<< std::setw(8) << (util::myClock() - startTime) << " s";

#ifdef HAVE_RUSAGE
					rusage usage;
					getrusage(RUSAGE_SELF, &usage);
					if (found_vmrss) output << ", " << std::setw(14) << vmrss_kb_global;
					else			 output << ", " << std::setw(14) << usage.ru_maxrss;
					//output << ", " << std::setw(7) << usage.ru_maxrss;
#endif
					output << " kB] :: " << format << '\n';

					std::string overFormat(output.str());
					std::vprintf(overFormat.c_str(), args);
				}
			}

//		template<>
//			void Logger::LogInternal<Singleton>(std::string format, std::va_list args)
//			{
//				if (thisRank == 0)
//				{
//#ifdef HAVE_RUSAGE
//					rusage usage;
//					getrusage(RUSAGE_SELF, &usage);
//#endif
//					char lead[30];
//					std::sprintf(lead, "![%06.1fs, %07ldkB]", util::myClock() - startTime, usage.ru_maxrss);
//					std::string newFormat = std::string(lead);
//					std::vprintf(newFormat.append(format).append("\n").c_str(), args);
//				}
//			}
	}
}
