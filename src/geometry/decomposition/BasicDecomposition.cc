
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/decomposition/BasicDecomposition.h"
#include "net/mpi.h"

#include <algorithm>
#include <fstream>

#define ALL_HISTOGRAM_DEFAULT_WIDTH 1.0

namespace hemelb
{
	namespace geometry
	{
		namespace decomposition
		{
			BasicDecomposition::BasicDecomposition(const Geometry& geometry,
					const lb::lattices::LatticeInfo& latticeInfo,
					const net::MpiCommunicator& communicator,
					const std::unordered_map<site_t, std::pair<uint16_t, uint16_t> >& blockInformation,
					const std::unordered_map<site_t, uint16_t>& blockWeights) :
				geometry(geometry), latticeInfo(latticeInfo), communicator(communicator),
				blockInformation(blockInformation), blockWeights(blockWeights)
			{
			}

			// Very dumb decomposition. For testing only!
			void BasicDecomposition::DecomposeDumb(
					std::unordered_map<site_t, proc_t>& procAssignedToEachBlock,
					sitedata_t nonEmptyBlocks)
			{
				site_t solidBlock = 0;
				for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
				{
					if (blockInformation.find(block) != blockInformation.end()) {
						procAssignedToEachBlock[block] = solidBlock/
							(double) nonEmptyBlocks*communicator.Size();
						solidBlock++;
					}
				}
			}

			// Rotate and allocate (owned locally?) the geometry.
			void BasicDecomposition::RotateAndAllocate(
					const double (&m_rot)[3][3][3], const double (&phi)[3],
					const double (&r_min)[3], const double (&r_max)[3],
					std::vector<ALL_Point<double>>& points,
					std::unordered_map<site_t, proc_t>& procAssignedToEachBlock,
					const std::vector<double>& l, const std::vector<double>& u,
					const int local_rank)
			{

				double x[3], z[3];
				double block_coords[3];

				for (site_t blockI = 0; blockI < geometry.GetBlockDimensions().x; blockI++)
					for (site_t blockJ = 0; blockJ < geometry.GetBlockDimensions().y; blockJ++)
						for (site_t blockK = 0; blockK < geometry.GetBlockDimensions().z; blockK++)
						{
							const site_t blockNumber = geometry.GetBlockIdFromBlockCoordinates(
									blockI,
									blockJ,
									blockK);

							// Does block contain fluid sites?
							if (blockInformation.find(blockNumber) == blockInformation.end())
								continue;

							block_coords[0] = blockI;
							block_coords[1] = blockJ;
							block_coords[2] = blockK;

							z[0] = block_coords[0] + 0.5;
							z[1] = block_coords[1] + 0.5;
							z[2] = block_coords[2] + 0.5;

							for (int j = 0; j < 3; j++) {
								x[j] = 0.0;
								for (int i = 0; i < 3; i++) {
									x[j] = x[j] + m_rot[0][i][j]*z[i];
								}
							}

							z[0] = x[0];
							z[1] = x[1];
							z[2] = x[2];

							for (int j = 0; j < 3; j++) {
								x[j] = 0.0;
								for (int i = 0; i < 3; i++) {
									x[j] = x[j] + m_rot[1][i][j]*z[i];
								}
							}

							z[0] = x[0];
							z[1] = x[1];
							z[2] = x[2];

							for (int j = 0; j < 3; j++) {
								x[j] = 0.0;
								for (int i = 0; i < 3; i++) {
									x[j] = x[j] + m_rot[2][i][j]*z[i];
								}
							}

							for (int k = 0; k < 3; k++)
								block_coords[k] = x[k] - r_min[k];

							if (
									block_coords[0] >= l.at(0) && block_coords[0] <= u.at(0) &&
									block_coords[1] >= l.at(1) && block_coords[1] <= u.at(1) &&
									block_coords[2] >= l.at(2) && block_coords[2] <= u.at(2))
							{
#ifdef HEMELB_USE_GMYPLUS
								ALL_Point<double> p(3,block_coords,blockWeights.at(blockNumber));
#else
								ALL_Point<double> p(3,block_coords,1);
#endif
								procAssignedToEachBlock[blockNumber] = local_rank;
								points.push_back(p);
							}
						}
			}

			// Some sort of space partitioning (inspired by sparta/balance_grid.cpp).
			void BasicDecomposition::DecomposeBlock(
					std::unordered_map<site_t, proc_t>& procAssignedToEachBlock, int noderank)
			{
				const int sys_dim = 3;

				double phi[3] = {};
				double r_min[3] = {}, r_max[3] = {};
				double x[3] = {}, z[3] = {};
				double m_rot[3][3][3] = {};

				std::fill_n(r_min, 3, std::numeric_limits<double>::max());

				phi[0] = 2.0*acos(-1.0)*0.123;
				phi[1] = 2.0*acos(-1.0)*0.823;
				phi[2] = 2.0*acos(-1.0)*0.323;

				m_rot[0][0][0] =  1.0;
				m_rot[0][1][1] =  cos(phi[0]);
				m_rot[0][1][2] =  sin(phi[0]);
				m_rot[0][2][1] = -sin(phi[0]);
				m_rot[0][2][2] =  cos(phi[0]);

				m_rot[1][0][0] =  cos(phi[1]);
				m_rot[1][0][2] = -sin(phi[1]);
				m_rot[1][1][1] =  1.0;
				m_rot[1][2][0] =  sin(phi[1]);
				m_rot[1][2][2] =  cos(phi[1]);

				m_rot[2][0][0] =  cos(phi[2]);
				m_rot[2][0][1] =  sin(phi[2]);
				m_rot[2][1][0] = -sin(phi[2]);
				m_rot[2][1][1] =  cos(phi[2]);
				m_rot[2][2][2] =  1.0;

				double block_coords[sys_dim];
				for (site_t blockI = 0; blockI < geometry.GetBlockDimensions().x; blockI++)
					for (site_t blockJ = 0; blockJ < geometry.GetBlockDimensions().y; blockJ++)
						for (site_t blockK = 0; blockK < geometry.GetBlockDimensions().z; blockK++)
						{
							const site_t blockNumber = geometry.GetBlockIdFromBlockCoordinates(
									blockI,
									blockJ,
									blockK);

							// Does block contain fluid sites?
							if (blockInformation.find(blockNumber) == blockInformation.end())
								continue;

							block_coords[0] = blockI;
							block_coords[1] = blockJ;
							block_coords[2] = blockK;

							z[0] = block_coords[0] + 0.5;
							z[1] = block_coords[1] + 0.5;
							z[2] = block_coords[2] + 0.5;

							for (int j = 0; j < 3; j++) {
								x[j] = 0.0;
								for (int i = 0; i < 3; i++) {
									x[j] = x[j] + m_rot[0][i][j]*z[i];
								}
							}

							z[0] = x[0];
							z[1] = x[1];
							z[2] = x[2];

							for (int j = 0; j < 3; j++) {
								x[j] = 0.0;
								for (int i = 0; i < 3; i++) {
									x[j] = x[j] + m_rot[1][i][j]*z[i];
								}
							}

							z[0] = x[0];
							z[1] = x[1];
							z[2] = x[2];

							for (int j = 0; j < 3; j++) {
								x[j] = 0.0;
								for (int i = 0; i < 3; i++) {
									x[j] = x[j] + m_rot[2][i][j]*z[i];
								}
							}

							for (int k = 0; k < 3; k++) {
								if (x[k] < r_min[k]) r_min[k] = x[k];
								if (x[k] > r_max[k]) r_max[k] = x[k];
							}
						}

#ifdef HEMELB_USE_GMYPLUS
				bool weighted_points = true;
#else
				bool weighted_points = false;
#endif
				double gamma = 2048.0;

				double box_size[sys_dim] = {};
				box_size[0] = std::fabs(r_max[0]-r_min[0])*1.2;
				box_size[1] = std::fabs(r_max[1]-r_min[1])*1.2;
				box_size[2] = std::fabs(r_max[2]-r_min[2])*1.2;

				std::vector<double> sys_size(6);
				sys_size.at(0) = 0.0;
				sys_size.at(1) = box_size[0];
				sys_size.at(2) = 0.0;
				sys_size.at(3) = box_size[1];
				sys_size.at(4) = 0.0;
				sys_size.at(5) = box_size[2];

				// setup of vector of points on each MPI rank
				std::vector<double> dummy(sys_dim);
				std::vector<ALL_Point<double>> points;
				int max_neighbors, loc_neighbors;

				// setup of cartesian communicator
				MPI_Comm cart_comm;
				int n_ranks;
				int local_rank;
				int local_coords[sys_dim];
				int periods[sys_dim];
				int glb_dim[sys_dim];

				// domain sizes
				std::vector<double> ds(sys_dim);
				// domain boundaries
				std::vector<double> l(sys_dim);
				std::vector<double> u(sys_dim);

				double min_ratio = 1.01;
				int min_step = -1;

				for (int i = 0; i < sys_dim; ++i)
				{
					periods[i] = 1;
					glb_dim[i] = 0;
				}

				// get number of total ranks
				MPI_Comm_size(communicator, &n_ranks);

				// get distribution into number of dimensions
				MPI_Dims_create(n_ranks, sys_dim, glb_dim);

				// create cartesian MPI communicator
				MPI_Cart_create(communicator, sys_dim, glb_dim, periods, 1, &cart_comm);

				// get local coordinates
				MPI_Cart_get(cart_comm, sys_dim, glb_dim, periods, local_coords);

				// get local rank
				MPI_Cart_rank(cart_comm, local_coords, &local_rank);

				// calculate domain extension
				for (int i = 0; i < sys_dim; ++i)
				{
					ds.at(i) = box_size[i]/(double)glb_dim[i];
					l.at(i) =        local_coords[i]  * ds.at(i);
					u.at(i) = (1.0 + local_coords[i]) * ds.at(i);
				}

				// generate vertices (equal to outline (for tensor only))
				ALL_Point<double>* tp = new ALL_Point<double>(2);
				std::vector<ALL_Point<double>>* tv = new std::vector<ALL_Point<double>>(2);
				tv->at(0) = *tp;
				delete tv;
				delete tp;

				int nvertices = 2;

				ALL_Point<double> lp(l);
				ALL_Point<double> up(u);
				std::vector<ALL_Point<double>>     vertices(nvertices,lp);
				std::vector<ALL_Point<double>> new_vertices(nvertices,lp);
				std::vector<ALL_Point<double>> old_vertices(nvertices,lp);

				vertices.at(0) = lp;
				vertices.at(1) = up;

				RotateAndAllocate(
						m_rot, phi, r_min, r_max,
						points, procAssignedToEachBlock,
						l, u, local_rank);

				procAssignedToEachBlock.clear();

				int n_points = points.size();
				int max_particles = 1;

				double* send;
				double* recv;

				MPI_Allreduce(&n_points,&max_particles,1,MPI_INT,MPI_MAX,cart_comm);

				max_particles = (int)std::ceil((double)max_particles * 1.5);

				if (local_rank == 0)
					log::Logger::Log<log::Info, log::OnePerCore>(
							"----> maximum number of sites on any process: %06i", max_particles);
				int max_neig = 27;

				recv = new double[max_neig * (sys_dim+1) * max_particles];
				send = new double[max_neig * (sys_dim+1) * max_particles];

				// find neighbors on cartesian communicator
				int l_neig[sys_dim];
				int r_neig[sys_dim];
				int self;
				for (int i = 0; i < sys_dim; ++i)
				{
					MPI_Cart_shift(cart_comm,i,+1,&self,&r_neig[i]);
					MPI_Cart_shift(cart_comm,i,-1,&self,&l_neig[i]);
					if (local_coords[i] == 0)
						l_neig[i] = MPI_PROC_NULL;
					if (local_coords[i] == glb_dim[i] - 1)
						r_neig[i] = MPI_PROC_NULL;
				}

				double d_min, d_max, d_ratio;
				double n_local;
				if (!weighted_points)
				{
					n_local = (double)n_points;
				}
				else
				{
					n_local = 0.0;
					for (auto p = points.begin(); p != points.end(); ++p)
						n_local += p->get_weight();
				}

				double n_total;
				MPI_Allreduce(&n_local,&n_total,1,MPI_DOUBLE,MPI_SUM,cart_comm);

				double avg_work = (double)n_total/(double)n_ranks;
				double n_min, n_max;
				MPI_Allreduce(&n_local,&n_min,1,MPI_DOUBLE,MPI_MIN,cart_comm);
				MPI_Allreduce(&n_local,&n_max,1,MPI_DOUBLE,MPI_MAX,cart_comm);

				d_min = n_min/avg_work;
				d_max = n_max/avg_work;
				d_ratio = (d_max - d_min)/(d_max + d_min);

				min_ratio = d_ratio;
				min_step = 0;

				// get initial work (simply number of points if unweighted)
				double total_points;
				MPI_Allreduce(&n_local,&total_points,1,MPI_DOUBLE,MPI_SUM,cart_comm);

				double n_total_points = total_points;
				double limit_efficiency = 0.5;

				// create ALL object
				ALL<double,double> lb_obj(sys_dim,vertices,gamma);
				for (int i_loop = 0; i_loop < 3; ++i_loop)
				{
					if (local_rank == 0)
						log::Logger::Log<log::Info, log::OnePerCore>(
								"----> [%03i] total work: %f", i_loop, n_total_points);

					MPI_Barrier(cart_comm);

					if (d_ratio < limit_efficiency)
					{
						gamma *= 2.0;
						limit_efficiency /= 2.0;
					}

					std::vector<double> work;
					std::vector<int> n_bins(3,-1);

					double histogram_width = ALL_HISTOGRAM_DEFAULT_WIDTH/gamma;

					{
						// work -> number of points on domain
						double lb(-1.0);
						double ub(-1.0);
						double overlap(0.0);
						int d = 2 - i_loop % 3;

						// compute number of bins in each direction
						lb = std::ceil(lp.x(d)/histogram_width)*histogram_width;
						ub = std::ceil(up.x(d)/histogram_width)*histogram_width;
						n_bins.at(d) = (ub - lb)/histogram_width;

						work = std::vector<double>(n_bins.at(d), 0.0);

						// compute histogram of work load
						for (auto p : points)
						{
							int idx = (int)(((p.x(d) - lb)/histogram_width));
							if (idx >= 0)
							{
								if (!weighted_points)
									work.at(idx) += 1.0;
								else
									work.at(idx) += p.get_weight();
							}
							else
							{
								if (!weighted_points)
									overlap += 1.0;
								else
									overlap += p.get_weight();
							}
						}

						// exchange overlapping workload (histograms might overlap
						// over the domain boundaries
						int rank_left, rank_right;
						MPI_Cart_shift(cart_comm,0,1,&rank_left,&rank_right);

						MPI_Request sreq, rreq;
						MPI_Status ssta, rsta;

						double recv_work;

						MPI_Isend(&overlap,
								1,
								MPI_DOUBLE,
								rank_left,
								0,
								cart_comm,
								&sreq);
						MPI_Irecv(&recv_work,
								1,
								MPI_DOUBLE,
								rank_right,
								0,
								cart_comm,
								&rreq);
						MPI_Wait(&sreq,&ssta);
						MPI_Wait(&rreq,&rsta);

						if (local_coords[d] != glb_dim[d] - 1)
							work.at(n_bins.at(d) - 1) += recv_work;
					}

					MPI_Barrier(cart_comm);

					int n_points_global = 0;
					MPI_Reduce(
							&n_points,
							&n_points_global,
							1,MPI_INT,MPI_SUM,0,communicator);

					lb_obj.set_work(work);
					lb_obj.set_communicator(cart_comm);
					lb_obj.setup(
							ALL_LB_t::HISTOGRAM);
					lb_obj.set_method_data(
							ALL_LB_t::HISTOGRAM,n_bins.data());
					lb_obj.set_sys_size(
							ALL_LB_t::HISTOGRAM, sys_size);
					lb_obj.set_vertices(vertices);
					lb_obj.balance(ALL_LB_t::HISTOGRAM);

					new_vertices = lb_obj.get_result_vertices();
					old_vertices = vertices;
					vertices = new_vertices;

					lp = vertices.at(0);
					up = vertices.at(vertices.size()-1);

					// determine current dimension
					int curr_dim = 2 - (i_loop % 3);

					MPI_Comm comm_col;
					// create temporary communicator to exchange borders
					switch (curr_dim)
					{
						case 0:
							{
								// x-plane
								MPI_Comm_split(cart_comm,
										local_coords[1]+local_coords[2]*glb_dim[1],
										local_coords[0],
										&comm_col);
								break;
							}
						case 1:
							{
								// y-plane
								MPI_Comm_split(cart_comm,
										local_coords[0]+local_coords[2]*glb_dim[0],
										local_coords[1],
										&comm_col);
								break;
							}
						case 2:
							{
								// z-plane
								MPI_Comm_split(cart_comm,
										local_coords[0]+local_coords[1]*glb_dim[0],
										local_coords[2],
										&comm_col);
								break;
							}
					}

					// vector to collect the borders into
					std::vector<double> borders(2*glb_dim[curr_dim]);
					std::vector<double> old_borders(2*glb_dim[curr_dim]);
					std::vector<double> local_borders(2);
					std::vector<double> old_border(2);

					local_borders.at(0) = vertices.at(0).x(curr_dim);
					local_borders.at(1) = vertices.at(1).x(curr_dim);
					old_border.at(0) = old_vertices.at(0).x(curr_dim);
					old_border.at(1) = old_vertices.at(1).x(curr_dim);

					int size, rank;
					MPI_Comm_rank(comm_col, &rank);
					MPI_Comm_size(comm_col, &size);

					// collect borders
					MPI_Allgather(local_borders.data(),
							2,
							MPI_DOUBLE,
							borders.data(),
							2,
							MPI_DOUBLE,
							comm_col);

					// collect borders
					MPI_Allgather(old_border.data(),
							2,
							MPI_DOUBLE,
							old_borders.data(),
							2,
							MPI_DOUBLE,
							comm_col);

					// compare old domains with new domains
					std::vector<int> send_neig;
					std::vector<int> recv_neig;

					for (int n = 0; n < glb_dim[curr_dim]; ++n)
					{
						if (
								(borders.at(2*n) <= old_border.at(0) &&
								 borders.at(2*n+1) > old_border.at(0)) ||
								(borders.at(2*n) <= old_border.at(1) &&
								 borders.at(2*n+1) >= old_border.at(1)) ||
								(borders.at(2*n) > old_border.at(0) &&
								 borders.at(2*n+1) < old_border.at(1)))
							send_neig.push_back(n);
						if (
								(old_borders.at(2*n) <= local_borders.at(0) &&
								 old_borders.at(2*n+1) > local_borders.at(0)) ||
								(old_borders.at(2*n) <= local_borders.at(1) &&
								 old_borders.at(2*n+1) >= local_borders.at(1)) ||
								(old_borders.at(2*n) > local_borders.at(0) &&
								 old_borders.at(2*n+1) < local_borders.at(1)))
							recv_neig.push_back(n);
					}

					// vectors to send and received points
					std::vector<std::vector<double>> send_vec(send_neig.size());

					// all points are sorted into the correct vectors
					for (auto p : points)
					{
						for (int n = 0; n < send_neig.size(); ++n)
						{
							int neig_id = send_neig.at(n);
							if (
									(borders.at(2*neig_id) <= p.x(curr_dim) &&
									 borders.at(2*neig_id+1) > p.x(curr_dim)))
							{
								send_vec.at(n).push_back(p.x(0));
								send_vec.at(n).push_back(p.x(1));
								send_vec.at(n).push_back(p.x(2));
								send_vec.at(n).push_back(p.get_weight());
							}
						}
					}

					// clear old point vector
					points.clear();
					n_points = 0;

					std::vector<int> n_send(send_neig.size());
					for (int n = 0; n < send_neig.size(); ++n)
						n_send.at(n) = send_vec.at(n).size();

					// communicate number of particles to be send to neighbors
					std::vector<int> n_recv(recv_neig.size());

					std::vector<MPI_Request> sreq(send_neig.size());
					std::vector<MPI_Request> rreq(recv_neig.size());
					std::vector<MPI_Status> ssta(send_neig.size());
					std::vector<MPI_Status> rsta(recv_neig.size());

					for (int n = 0; n < send_neig.size(); ++n)
					{
						MPI_Isend(n_send.data()+n,
								1,
								MPI_INT,
								send_neig.at(n),
								2000,
								comm_col,
								sreq.data()+n);
					}
					for (int n = 0; n < recv_neig.size(); ++n)
					{
						MPI_Irecv(n_recv.data()+n,
								1,
								MPI_INT,
								recv_neig.at(n),
								2000,
								comm_col,
								rreq.data()+n);
					}
					MPI_Waitall(send_neig.size(),sreq.data(),ssta.data());

					int recvs = 0;

					std::vector<std::vector<double>> recv_vec(recv_neig.size());
					while (recvs < recv_neig.size())
					{
						int idx;
						MPI_Waitany(recv_neig.size(),rreq.data(),&idx,rsta.data());
						recv_vec.at(idx).resize(n_recv.at(idx));
						recvs++;
					}

					// send points from old domains to new domains
					for (int n = 0; n < send_neig.size(); ++n)
					{
						MPI_Isend(send_vec.at(n).data(),
								send_vec.at(n).size(),
								MPI_DOUBLE,
								send_neig.at(n),
								3000,
								comm_col,
								sreq.data()+n);
					}
					for (int n = 0; n < recv_neig.size(); ++n)
					{
						MPI_Irecv(recv_vec.at(n).data(),
								recv_vec.at(n).size(),
								MPI_DOUBLE,
								recv_neig.at(n),
								3000,
								comm_col,
								rreq.data()+n);
					}
					MPI_Waitall(send_neig.size(),sreq.data(),ssta.data());

					recvs = 0;
					while (recvs < recv_neig.size())
					{
						int idx;
						MPI_Waitany(recv_neig.size(),rreq.data(),&idx,rsta.data());
						for (int p = 0; p < recv_vec.at(idx).size(); p+=4)
						{
							ALL_Point<double> tmp(3,
									recv_vec.at(idx).data()+p,
									recv_vec.at(idx).at(p+3)
									);
							points.push_back(tmp);
						}
						recvs++;
					}
					n_points = points.size();

					MPI_Comm_free(&comm_col);

					// calculate quality parameters
					if (!weighted_points)
					{
						n_local = (double)n_points;
					}
					else
					{
						n_local = 0.0;
						for (auto p = points.begin(); p != points.end(); ++p)
							n_local += p->get_weight();
					}

					MPI_Allreduce(&n_local,&n_total_points,1,MPI_DOUBLE,MPI_SUM,cart_comm);

					MPI_Allreduce(&n_local,&n_total,1,MPI_DOUBLE,MPI_SUM,cart_comm);
					avg_work = n_total/(double)n_ranks;
					MPI_Allreduce(&n_local,&n_min,1,MPI_DOUBLE,MPI_MIN,cart_comm);
					MPI_Allreduce(&n_local,&n_max,1,MPI_DOUBLE,MPI_MAX,cart_comm);
					d_min = n_min/avg_work;
					d_max = n_max/avg_work;
					d_ratio = (d_max - d_min)/(d_max + d_min);

					double loc_diff, tot_diff;
					loc_diff = (n_local - avg_work) * (n_local - avg_work);
					MPI_Reduce(&loc_diff,&tot_diff,1,MPI_DOUBLE,MPI_SUM,0,cart_comm);

					d_min = n_min/avg_work;
					d_max = n_max/avg_work;
					d_ratio = (d_max - d_min)/(d_max + d_min);

					if (local_rank == 0)
					{
						log::Logger::Log<log::Info, log::OnePerCore>(
								"----> [%03i] d_min: %f, d_max: %f, d_ratio: %f ", i_loop,
								d_min, d_max, d_ratio);

						if (d_ratio < min_ratio)
						{
							min_ratio = d_ratio;
							min_step  = i_loop;
						}
					}

					if (std::abs(n_total_points - total_points) > 1e-6)
					{
						log::Logger::Log<log::Error, log::OnePerCore>("Work has been lost!");
					}
				}

				for (int i = 0; i < sys_dim; ++i)
				{
					l[i] = lp.x(i);
					u[i] = up.x(i);
				}

				RotateAndAllocate(
						m_rot, phi, r_min, r_max,
						points, procAssignedToEachBlock,
						l, u, local_rank);

				sitedata_t* blocks_local;
				sitedata_t* blocks_neigh;

				int displs[n_ranks] = {};
				int n_recv[n_ranks] = {};
				int neighs[n_ranks] = {};

				int n_send = procAssignedToEachBlock.size();

				MPI_Allgather(&n_send    ,1,MPI_INT,&n_recv[0],1,MPI_INT,cart_comm);
				MPI_Allgather(&local_rank,1,MPI_INT,&neighs[0],1,MPI_INT,cart_comm);

				displs[0] = 0;
				for (int i = 1; i < n_ranks; ++i)
					displs[i] = displs[i-1] + n_recv[i-1];

				blocks_local = new sitedata_t[procAssignedToEachBlock.size()];
				blocks_neigh = new sitedata_t[       blockInformation.size()]; sitedata_t i = 0;
				for (std::unordered_map<site_t, proc_t>::const_iterator iter = procAssignedToEachBlock.begin();
						iter != procAssignedToEachBlock.end(); ++iter)
					blocks_local[i++] = iter->first;

				MPI_Allgatherv(
						&blocks_local[0],
						n_send,
						MPI_UNSIGNED_LONG,
						&blocks_neigh[0],
						n_recv,
						displs,
						MPI_UNSIGNED_LONG,
						cart_comm);

				if (noderank == 1)
				{
					int j = 1; proc_t rank = neighs[0];
					for (int i = 0; i < blockInformation.size(); ++i)
					{
						if (i == displs[j])
							rank = neighs[j++];
						procAssignedToEachBlock[blocks_neigh[i]] = rank;
					}

					std::unordered_map<proc_t, proc_t> test1;
					std::unordered_map<proc_t, proc_t> test2;
					for (std::unordered_map<site_t, proc_t>::const_iterator iter = procAssignedToEachBlock.begin();
							iter != procAssignedToEachBlock.end(); ++iter)
					{
						// Check that we have the correct number of entries.
						if (test1.find(iter->second) == test1.end())
							test1[iter->second] = 0;
						test1.at(iter->second)++;

						// Check for duplicate entries.
						if (test2.find(iter->first) == test2.end())
							test2[iter->first] = 0;
						else std::cout << "duplicate!" << std::endl;
					}
				}

				//{
				//	int i = 0;
				//	char hostname[256];
				//	gethostname(hostname, sizeof(hostname));
				//	printf("PID %d (%d) on %s ready for attach\n", getpid(), communicator.Rank(), hostname);
				//	fflush(stdout);
				//	while (0 == i)
				//		sleep(5);
				//}

				delete blocks_local;
				delete blocks_neigh;
				delete recv;
				delete send;
			}

			void BasicDecomposition::Decompose(std::unordered_map<site_t, proc_t>& procAssignedToEachBlock)
			{
				// Keep a count of the number of non-empty blocks that haven't yet been assigned
				// a processor.
				site_t unvisitedFluidBlockCount = 0;
				for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
				{
					if (blockInformation.find(block) != blockInformation.end())
					{
#ifdef HEMELB_USE_GMYPLUS
						unvisitedFluidBlockCount += blockWeights.at(block);
#else
						unvisitedFluidBlockCount++;
#endif
					}
				}

				DivideBlocks(procAssignedToEachBlock,
						unvisitedFluidBlockCount,
						geometry,
						communicator.Size(),
						blockInformation);
			}

			//void BasicDecomposition::Validate(std::vector<proc_t>& procAssignedToEachBlock)
			//{
			//	log::Logger::Log<log::Debug, log::OnePerCore>("Validating procForEachBlock");

			//	std::vector<proc_t> procForEachBlockRecv = communicator.AllReduce(procAssignedToEachBlock, MPI_MAX);

			//	for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
			//	{
			//		if (procAssignedToEachBlock[block] != procForEachBlockRecv[block])
			//		{
			//			log::Logger::Log<log::Critical, log::OnePerCore>("At least one other proc thought block %li should be on proc %li but we locally had it as %li",
			//					block,
			//					procAssignedToEachBlock[block],
			//					procForEachBlockRecv[block]);
			//		}
			//	}
			//}

			void BasicDecomposition::DivideBlocks(std::unordered_map<site_t, proc_t>& unitForEachBlock,
					site_t unassignedBlocks,
					const Geometry& geometry,
					const proc_t unitCount,
					const std::unordered_map<site_t, std::pair<uint16_t, uint16_t> >& blockInformation)
			{
				// Initialise the unit being assigned to, and the approximate number of blocks
				// required on each unit.
				proc_t currentUnit = 0;

				site_t targetBlocksPerUnit = (site_t) ceil((double) unassignedBlocks / (double) (communicator.Size()));

				// Create a set to monitor whether each block has been assigned yet.
				std::unordered_set<site_t> blockAssigned;

				// Create lists of the current edge of blocks on the current proc and the edge being expanded into
				std::vector<BlockLocation> currentEdge;
				std::vector<BlockLocation> expandedEdge;

				// Domain Decomposition. Pick a site. Set it to the rank we are
				// looking at. Find its neighbours and put those on the same
				// rank, then find the next-nearest neighbours, etc. until we
				// have a completely joined region, or there are enough fluid
				// sites on the rank.  In the former case, start again at
				// another site. In the latter case, move on to the next rank.
				// Do this until all sites are assigned to a rank. There is a
				// high chance of of all sites on a rank being joined.

				site_t blockNumber = -1;
				site_t blocksOnCurrentProc = 0;

				// Iterate over all blocks.
				for (site_t blockCoordI = 0; blockCoordI < geometry.GetBlockDimensions().x; blockCoordI++)
				{
					for (site_t blockCoordJ = 0; blockCoordJ < geometry.GetBlockDimensions().y; blockCoordJ++)
					{
						for (site_t blockCoordK = 0; blockCoordK < geometry.GetBlockDimensions().z; blockCoordK++)
						{
							// Block number is the number of the block we're currently on.
							blockNumber++;

							if (blockInformation.find(blockNumber) == blockInformation.end())
							{
								continue;
							}

							// Alternatively, if this block has already been assigned, move on.
							if (blockAssigned.find(blockNumber) != blockAssigned.end())
							{
								continue;
							}

							// Assign this block to the current unit.
							blockAssigned.insert(blockNumber);
							unitForEachBlock[blockNumber] = currentUnit;

#ifdef HEMELB_USE_GMYPLUS
							blocksOnCurrentProc += blockWeights.at(blockNumber);
#else
							blocksOnCurrentProc++;
#endif

							// Record the location of this initial block.
							currentEdge.clear();
							BlockLocation lNew(blockCoordI, blockCoordJ, blockCoordK);
							currentEdge.push_back(lNew);

							// The subdomain can grow.
							bool isRegionGrowing = true;

							int blockAssignedCount = blockAssigned.size();

							// While the region can grow (i.e. it is not bounded by solids or visited
							// sites), and we need more sites on this particular rank.
							while (blocksOnCurrentProc < targetBlocksPerUnit && isRegionGrowing)
							{
								expandedEdge.clear();

								// Sites added to the edge of the mClusters during the iteration.
								isRegionGrowing = Expand(expandedEdge,
										blockAssigned,
										unitForEachBlock,
										blocksOnCurrentProc,
										currentEdge,
										currentUnit,
										targetBlocksPerUnit);

								// When the new layer of edge sites has been found, swap the buffers for
								// the current and new layers of edge sites.
								currentEdge.swap(expandedEdge);
							}

							//// If we tried to expand and did not get far.
							//// Prevents creation of small isolated regions.
							//if (((blockAssigned.size() - blockAssignedCount) < 15) && (currentUnit > 0))
							//{
							//	for (std::unordered_map<site_t, proc_t>::const_iterator iter = unitForEachBlock.begin();
							//			iter != unitForEachBlock.end(); ++iter)
							//	{
							//		if (iter->second == currentUnit)
							//			// Allocate these blocks to previous unit.
							//			unitForEachBlock.at(iter->first) = currentUnit - 1;
							//	}
							//}

							// If we have enough sites, we have finished.
							if (blocksOnCurrentProc >= targetBlocksPerUnit)
							{
								++currentUnit;

								unassignedBlocks -= blocksOnCurrentProc;
								targetBlocksPerUnit = (site_t) ceil((double) unassignedBlocks / (double) (unitCount - currentUnit));

								blocksOnCurrentProc = 0;
							}
							// If not, we have to start growing a different region for the same rank:
							// region expansions could get trapped.
						} // Block co-ord k
					} // Block co-ord j
				} // Block co-ord i

				//blockNumber = -1;
				//// Check which ranks own the neighbouring blocks of each block.
				//for (site_t blockCoordI = 0; blockCoordI < geometry.GetBlockDimensions().x; blockCoordI++)
				//{
				//	for (site_t blockCoordJ = 0; blockCoordJ < geometry.GetBlockDimensions().y; blockCoordJ++)
				//	{
				//		for (site_t blockCoordK = 0; blockCoordK < geometry.GetBlockDimensions().z; blockCoordK++)
				//		{
				//			// Block number is the number of the block we're currently on.
				//			if (unitForEachBlock.find(++blockNumber) == unitForEachBlock.end())
				//			{
				//				continue;
				//			}

				//			// Create map to record rank of neighbouring blocks.
				//			std::unordered_map<proc_t, int> KFC;
				//			for (site_t neighI = util::NumericalFunctions::max<site_t>(0, blockCoordI - 1); (neighI
				//						<= (blockCoordI + 1)) && (neighI < geometry.GetBlockDimensions().x); ++neighI)
				//			{
				//				for (site_t neighJ = util::NumericalFunctions::max<site_t>(0, blockCoordJ - 1); (neighJ
				//							<= (blockCoordJ + 1)) && (neighJ < geometry.GetBlockDimensions().y); ++neighJ)
				//				{
				//					for (site_t neighK = util::NumericalFunctions::max<site_t>(0, blockCoordK - 1); (neighK
				//								<= (blockCoordK + 1)) && (neighK < geometry.GetBlockDimensions().z); ++neighK)
				//					{
				//						site_t neighBlockId = geometry.GetBlockIdFromBlockCoordinates(neighI, neighJ, neighK);
				//						if (unitForEachBlock.find(neighBlockId) != unitForEachBlock.end())
				//						{
				//							proc_t neighRank = unitForEachBlock.at(neighBlockId);
				//							if (KFC.find(neighRank) == KFC.end())
				//								KFC[neighRank] = 0;
				//							KFC.at(neighRank)++;
				//						}
				//					}
				//				}
				//			}

				//			// If no neighbouring blocks are owned by the rank owning the source block...
				//			// ... or if the number of blocks owned by rank owning source block is < 8.
				//			if (
				//					KFC.find(unitForEachBlock.at(blockNumber)) == KFC.end() ||
				//					KFC.at(  unitForEachBlock.at(blockNumber)) < 8)
				//			{
				//				int maxValue = 0; proc_t maxValueProc = unitForEachBlock.at(blockNumber);
				//				for (std::unordered_map<proc_t, int>::const_iterator iter = KFC.begin();
				//						iter != KFC.end(); ++iter)
				//				{
				//					if (iter->second > maxValue)
				//					{
				//						maxValue     = iter->second;
				//						maxValueProc = iter->first;
				//					}
				//				}
				//				unitForEachBlock.at(blockNumber) = maxValueProc;
				//			}
				//		}
				//	}
				//}

				// To measure quality of distribution.
				std::vector<sitedata_t> totalBlockWeights(communicator.Size(), 0);

				// Iterate over all blocks (again).
				for (site_t blockNumber = 0; blockNumber < geometry.GetBlockCount(); ++blockNumber)
				{
					// Weight of all blocks on partition.
					if (unitForEachBlock.find(blockNumber) != unitForEachBlock.end())
#ifdef HEMELB_USE_GMYPLUS
						totalBlockWeights[unitForEachBlock.at(blockNumber)] += blockWeights.at(blockNumber);
#else
					totalBlockWeights[unitForEachBlock.at(blockNumber)]++;
#endif
				}

				sitedata_t max = *std::max_element(
						totalBlockWeights.begin(), totalBlockWeights.end());
				sitedata_t min = *std::min_element(
						totalBlockWeights.begin(), totalBlockWeights.end());
				double average = std::accumulate(
						totalBlockWeights.begin(), totalBlockWeights.end(), 0)/(double)totalBlockWeights.size();

				double d_min = (double)min / average;
				double d_max = (double)max / average;

				double d_ratio = (d_max - d_min) / (d_max + d_min);

				//communicator.Barrier();
				if (communicator.Rank() == 0)
					log::Logger::Log<log::Info, log::OnePerCore>("----> load distribution: %f", d_ratio);//totalBlockWeights[communicator.Rank()]/average);
			}

			bool BasicDecomposition::Expand(std::vector<BlockLocation>& expansionBlocks,
					std::unordered_set<site_t>& blockAssigned,
					std::unordered_map<site_t, proc_t>& unitForEachBlock,
					site_t &blocksOnCurrentUnit,
					const std::vector<BlockLocation>& edgeBlocks,
					const proc_t currentUnit,
					const site_t blocksPerUnit)
			{
				bool regionExpanded = false;

				// For sites on the edge of the domain (sites_a), deal with the neighbours.
				for (unsigned int edgeBlockId = 0; (edgeBlockId < edgeBlocks.size()) && (blocksOnCurrentUnit < blocksPerUnit); edgeBlockId++)
				{
					const BlockLocation& edgeBlockCoords = edgeBlocks[edgeBlockId];

					for (Direction direction = 1; direction < latticeInfo.GetNumVectors() && blocksOnCurrentUnit < blocksPerUnit; direction++)
					{
						// Record neighbour location.
						BlockLocation neighbourCoords = edgeBlockCoords + latticeInfo.GetVector(direction);

						//{
						//	int i = 0;
						//	char hostname[256];
						//	gethostname(hostname, sizeof(hostname));
						//	printf("PID %d (%d) on %s ready for attach\n", getpid(), communicator.Rank(), hostname);
						//	fflush(stdout);
						//	while (0 == i)
						//		sleep(5);
						//}

						// Move on if neighbour is outside the bounding box.
						if (!geometry.AreBlockCoordinatesValid(neighbourCoords))
						{
							continue;
						}

						site_t neighBlockId = geometry.GetBlockIdFromBlockCoordinates(neighbourCoords.x,
								neighbourCoords.y,
								neighbourCoords.z);

						// Don't use this block if it has no fluid sites, or if it has already been assigned to a processor.
						if (blockInformation.find(neighBlockId) == blockInformation.end() ||
								blockAssigned.find(neighBlockId) != blockAssigned.end())
						{
							continue;
						}

						// Set the rank for a neighbour and update the fluid site counters.
						blockAssigned.insert(neighBlockId);
						unitForEachBlock[neighBlockId] = currentUnit;
#ifdef HEMELB_USE_GMYPLUS
						blocksOnCurrentUnit += blockWeights.at(neighBlockId);
#else
						blocksOnCurrentUnit++;
#endif

						// Neighbour was found, so the region can grow.
						regionExpanded = true;

						// Record the location of the neighbour.
						expansionBlocks.push_back(neighbourCoords);
					}
				}
				return regionExpanded;
			}
		} /* namespace decomposition */
	} /* namespace geometry */
} /* namespace hemelb */

//int px = 0, py = 0, pz = 0;
//int ix, ipx, iy, ipy, iz, ipz;
//
//int nx = geometry.GetBlockDimensions().x;
//int ny = geometry.GetBlockDimensions().y;
//int nz = geometry.GetBlockDimensions().z;
//
//// Returns number of partitions in each direction (minimize surface area).
//Procs2Grid(nx, ny, nz, px, py, pz);
//
//if (px*py*pz != communicator.Size())
//	log::Logger::Log<log::Error, log::OnePerCore>("Bad grid of processors!");
//
//for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
//{
//	if (blockInformation.find(block) == blockInformation.end()) continue;
//
//	site_t blockIJData = block / geometry.GetBlockDimensions().z;
//
//	// Find block x, y, z.
//	ix = blockIJData / geometry.GetBlockDimensions().y;
//	iy = blockIJData % geometry.GetBlockDimensions().y;
//	iz = block       % geometry.GetBlockDimensions().z;
//
//	ipx = ix*px / nx;
//	ipy = iy*py / ny;
//	ipz = iz*pz / nz;
//
//	// Assign rank to block.
//	procAssignedToEachBlock[block] = ipz*px*py + ipy*px + ipx;
//}

//void BasicDecomposition::Procs2Grid(int nx, int ny, int nz,
//		int &px, int &py, int &pz)
//{
//	int upx = px;
//	int upy = py;
//	int upz = pz;
//
//	int nprocs = communicator.Size();
//
//	// Determine cross-sectional areas:
//	// area[0] = xy, area[1] = xz, area[2] = yz
//	double area[3];
//	area[0] = nx*ny;
//	area[1] = nx*nz;
//	area[2] = ny*nz;
//
//	double bestsurf = 2.0 * (area[0]+area[1]+area[2]);
//
//	double surf;
//	int ipx,ipy,ipz,valid;
//	ipx = 1;
//
//	// Loop through all possible factorizations of nprocs.
//	// 'surf' is surface area of a proc sub-domain.
//	while (ipx <= nprocs) {
//		valid = 1;
//		if (upx && ipx != upx) valid = 0;
//		if (nprocs % ipx) valid = 0;
//		if (!valid) {
//			ipx++;
//			continue;
//		}
//
//		ipy = 1;
//		while (ipy <= nprocs/ipx) {
//			valid = 1;
//			if (upy && ipy != upy) valid = 0;
//			if ((nprocs/ipx) % ipy) valid = 0;
//			if (!valid) {
//				ipy++;
//				continue;
//			}
//
//			ipz = nprocs/ipx/ipy;
//			valid = 1;
//			if (upz && ipz != upz) valid = 0;
//			if (!valid) {
//				ipy++;
//				continue;
//			}
//
//			surf = area[0]/ipx/ipy + area[1]/ipx/ipz + area[2]/ipy/ipz;
//			if (surf < bestsurf) {
//				bestsurf = surf;
//				px = ipx; py = ipy; pz = ipz;
//			} ipy++;
//		} ipx++;
//	}
//}
