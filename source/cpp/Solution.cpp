#include<Solution.hpp>
#include<Data.hpp>
#include<EdgeDataHandler.hpp>
#include<FaceDataHandler.hpp>
#include<omp.h>
Solution::Solution()
{
	auto max_thread = omp_get_max_threads();
}

Solution::~Solution()
{
}

GeometryResult Solution::ReadGeometry(const std::string& file_name)
{
	return geometry.Read(file_name);
}

GeometryResult Solution::WriteSolution(const std::string& file_name)
{
	GeometryResult result;
	auto vertex_list = geometry.GetVertexList();
	auto face_list = geometry.GetFaceList();
	std::ofstream myfile;
	auto file = file_name + ".dat";
	myfile.open(file);
	myfile << "TITLE = \"title\"\n";
	myfile << "VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"P\"\n";
	myfile << "ZONE N = " << vertex_list.size() << ", E = " << face_list.size() << ", DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE\n";
	myfile << "VARLOCATION = ([3,4,5] = CELLCENTERED)" << "\n";
	for (auto it : vertex_list)
	{
		auto coord = it->GetPositionVector();
		myfile << coord[0] << " ";
	}
	myfile << "\n";
	for (auto it : vertex_list)
	{
		auto coord = it->GetPositionVector();
		myfile << coord[1] << "\n";
	}
	myfile << "\n";
	for (auto it : face_list)
	{
		myfile << it->u << "\n";
	}
	myfile << "\n";

	for (auto it : face_list)
	{
		myfile << it->v << "\n";
	}
	myfile << "\n";

	for (auto it : face_list)
	{
		myfile << it->p << "\n";
	}
	myfile << "\n";

	for (auto it : face_list)
	{
		myfile << it->GetHalfEdge()[0]->GetStart()->GetId() << " " << it->GetHalfEdge()[1]->GetStart()->GetId() << " " << it->GetHalfEdge()[2]->GetStart()->GetId() << "\n";
	}

	myfile.close();

	result.success = true;
	return result;
}

void Solution::SetBoundaryCondition(std::string name, const int& type, const double& value)
{
	BoundaryCondition bc;
	bc.SetBoundaryCondition(type, value);
	geometry.GetPhysicalgroupByName(name)->SetBoundaryCondition(bc);
}

void Solution::SetTimeStep(const double& dt)
{
	parameter.dt = dt;
}

void Solution::SetDensity(const double& density)
{
	parameter.density = density;
}

void Solution::SetViscosity(const double& viscosity)
{
	parameter.viscosity = viscosity;
}

void Solution::SetMaxTimeStep(const int& max_timestep)
{
	parameter.max_timestep = max_timestep;
}

Geometry& Solution::GetGeometry()
{
	return geometry;
}

void Solution::SolveMomentum()
{
	auto& faces = geometry.GetFaceList();

	LIS_INT n = faces.size();

	lis_matrix_create(0, &A_us);
	lis_matrix_create(0, &A_p );
	lis_matrix_create(0, &A_u );

	lis_vector_create(0, &b_us);
	lis_vector_create(0, &b_vs);
	lis_vector_create(0, &b_u );
	lis_vector_create(0, &b_v );
	lis_vector_create(0, &b_p );
	lis_vector_create(0, &x   );

	lis_matrix_set_size(A_us, 0, n);
	lis_matrix_set_size(A_p , 0, n); 
	lis_matrix_set_size(A_u , 0, n);

	lis_vector_set_size(b_us, 0, n);
	lis_vector_set_size(b_vs, 0, n);
	lis_vector_set_size(b_u , 0, n);
	lis_vector_set_size(b_v , 0, n);
	lis_vector_set_size(b_p , 0, n);
	lis_vector_set_size(x   , 0, n);

	lis_solver_create(&solver_momentum);
	lis_solver_create(&solver_pressure);
	lis_solver_create(&solver_velocity);

	lis_solver_set_option("-i bicrstab -p ilu", solver_momentum);
	lis_solver_set_option("-i bicrstab -p ilu", solver_pressure);
	lis_solver_set_option("-i bicrstab -p ilu", solver_velocity);
	lis_solver_set_option("-tol 1.0e-12", solver_momentum);
	lis_solver_set_option("-tol 1.0e-12", solver_pressure);
	lis_solver_set_option("-tol 1.0e-12", solver_velocity);
	//lis_solver_set_option("-maxiter 10000", solver_momentum);
	//lis_solver_set_option(const_cast<char*> (("-omp_num_threads " + std::to_string(max_thread)).c_str()), solver_momentum);
	//lis_solver_set_option(const_cast<char*> (("-omp_num_threads " + std::to_string(max_thread)).c_str()), solver_pressure);
	//lis_solver_set_option(const_cast<char*> (("-omp_num_threads " + std::to_string(max_thread)).c_str()), solver_velocity);


	auto physical_groups = geometry.GetPhysicalGroup();

	//setting A mat
	//----------------------------------------------------------------------------------------
	for (auto face : faces)
	{
		auto cell_volume = face->GetArea();
		face->central_term = parameter.density * cell_volume / parameter.dt;
		auto diagonal_term = 0.0;
		auto& half_edges = face->GetHalfEdge();

		for (auto half_edge : half_edges)
		{
			auto neighbour_half_edge = half_edge->GetNeighbourHalfEdge();
			auto neighbour_face = neighbour_half_edge->GetFace();

			if (neighbour_face != nullptr)
			{
				auto face_area = half_edge->GetEdgeVector().Abs();
				auto cell_centeroid_distance = (face->GetCentroid() - neighbour_face->GetCentroid()).Abs();
				auto coefficient = parameter.viscosity * face_area / cell_centeroid_distance;
				face->central_term += coefficient;
				diagonal_term = -coefficient;
				lis_matrix_set_value(LIS_INS_VALUE, face->GetId(), neighbour_face->GetId(), diagonal_term, A_us);
			}
		}
		lis_matrix_set_value(LIS_INS_VALUE, face->GetId(), face->GetId(), face->central_term, A_us);
	}

	for (auto& physical_group : physical_groups)
	{
		auto bc = physical_group.second->GetBoundaryCondition();
		if (bc.GetType() == WALL || bc.GetType() == INFLOW)
		{
			auto boundary_edges = physical_group.second->GetEdges();
			for (auto boundary_edge : boundary_edges)
			{
				Face* face = nullptr;
				HalfEdge* he = nullptr;
				if (boundary_edge->GetHalfEdge(0)->GetFace() == nullptr)
				{
					he = boundary_edge->GetHalfEdge(1);
					face = he->GetFace();
				}
				else
				{
					he = boundary_edge->GetHalfEdge(0);
					face = he->GetFace();
				}

				auto cell_centeroid_to_wall_distance = abs((boundary_edge->GetCenter() - face->GetCentroid()) * (he->GetNormal()));
				auto face_area = he->GetEdgeVector().Abs();
				auto coefficient = parameter.viscosity * face_area / cell_centeroid_to_wall_distance;
				face->central_term += coefficient;
				lis_matrix_set_value(LIS_INS_VALUE, face->GetId(), face->GetId(), face->central_term, A_us);


				if (bc.GetType() == INFLOW)
				{
					he->mass_flux = parameter.density * (bc.GetValue() * he->GetEdgeVector().Abs() * he->GetNormal().GetDx());
				}
			}
		}
	}
	lis_matrix_set_type(A_us, LIS_MATRIX_CSR);
	lis_matrix_assemble(A_us);
	
	//end setting A_us mat
	//----------------------------------------------------------------------------------------


	//setting A_p mat
	//----------------------------------------------------------------------------------------
	for (auto face : faces)
	{
		auto cell_volume = face->GetArea();
		face->central_term = 0.0;
		auto diagonal_term = 0.0;
		auto& half_edges = face->GetHalfEdge();

		for (auto half_edge : half_edges)
		{
			auto neighbour_half_edge = half_edge->GetNeighbourHalfEdge();
			auto neighbour_face = neighbour_half_edge->GetFace();

			if (neighbour_face != nullptr)
			{
				auto face_area = half_edge->GetEdgeVector().Abs();
				auto cell_centeroid_distance = (face->GetCentroid() - neighbour_face->GetCentroid()).Abs();
				auto coefficient = face_area / cell_centeroid_distance;
				face->central_term -= coefficient;
				diagonal_term = coefficient;
				lis_matrix_set_value(LIS_INS_VALUE, face->GetId(), neighbour_face->GetId(), diagonal_term, A_p);
			}
		}
		lis_matrix_set_value(LIS_INS_VALUE, face->GetId(), face->GetId(), face->central_term, A_p);
	}

	lis_matrix_set_type(A_p, LIS_MATRIX_CSR);
	lis_matrix_assemble(A_p);
	//end setting A_p mat
	//----------------------------------------------------------------------------------------


	//setting A_u mat
//----------------------------------------------------------------------------------------
	for (auto face : faces)
	{
		auto cell_volume = face->GetArea();
		face->central_term = parameter.density * cell_volume / parameter.dt;
		auto diagonal_term = 0.0;
		auto& half_edges = face->GetHalfEdge();

		for (auto half_edge : half_edges)
		{
			auto neighbour_half_edge = half_edge->GetNeighbourHalfEdge();
			auto neighbour_face = neighbour_half_edge->GetFace();

			if (neighbour_face != nullptr)
			{
				auto face_area = half_edge->GetEdgeVector().Abs();
				auto cell_centeroid_distance = (face->GetCentroid() - neighbour_face->GetCentroid()).Abs();
				auto coefficient = parameter.viscosity * face_area / cell_centeroid_distance;
				face->central_term += coefficient;
				face->central_term += (1.0 - parameter.gamma) * std::max(half_edge->mass_flux, 0.0);
				face->central_term += parameter.gamma * ((half_edge->mass_flux*neighbour_face->GetArea())/
					                                    (face->GetArea()+ neighbour_face->GetArea()));

				diagonal_term = -coefficient;
				diagonal_term -= (1.0 - parameter.gamma) * std::max(-half_edge->mass_flux, 0.0);
				diagonal_term += parameter.gamma * ((half_edge->mass_flux * face->GetArea()) /
					                               (face->GetArea() + neighbour_face->GetArea()));

				lis_matrix_set_value(LIS_INS_VALUE, face->GetId(), neighbour_face->GetId(), diagonal_term, A_u);
			}
		}
		lis_matrix_set_value(LIS_INS_VALUE, face->GetId(), face->GetId(), face->central_term, A_u);
	}

	for (auto physical_group : physical_groups)
	{
		auto bc = physical_group.second->GetBoundaryCondition();
		if (bc.GetType() == WALL || bc.GetType() == INFLOW)
		{
			auto boundary_edges = physical_group.second->GetEdges();
			for (auto boundary_edge : boundary_edges)
			{
				Face* face = nullptr;
				HalfEdge* he = nullptr;
				if (boundary_edge->GetHalfEdge(0)->GetFace() == nullptr)
				{
					he = boundary_edge->GetHalfEdge(1);
					face = he->GetFace();
				}
				else
				{
					he = boundary_edge->GetHalfEdge(0);
					face = he->GetFace();
				}

				auto cell_centeroid_to_wall_distance = abs((boundary_edge->GetCenter() - face->GetCentroid()) * (he->GetNormal()));
				auto face_area = he->GetEdgeVector().Abs();
				auto coefficient = parameter.viscosity * face_area / cell_centeroid_to_wall_distance;
				face->central_term += coefficient;
				lis_matrix_set_value(LIS_INS_VALUE, face->GetId(), face->GetId(), face->central_term, A_u);
			}
		}
		if (bc.GetType() == OUTFLOW)
		{
			auto boundary_edges = physical_group.second->GetEdges();
			for (auto boundary_edge : boundary_edges)
			{
				Face* face = nullptr;
				HalfEdge* he = nullptr;
				if (boundary_edge->GetHalfEdge(0)->GetFace() == nullptr)
				{
					he = boundary_edge->GetHalfEdge(1);
					face = he->GetFace();
				}
				else
				{
					he = boundary_edge->GetHalfEdge(0);
					face = he->GetFace();
				}
				face->central_term += he->mass_flux;
			}
		}
	}

	lis_matrix_set_type(A_u, LIS_MATRIX_CSR);
	lis_matrix_assemble(A_u);

	//end setting A_u mat
	//----------------------------------------------------------------------------------------
;

	for (int int_iter = 0; int_iter < 25; ++int_iter)
	{
		std::cout << "iter: " << int_iter << std::endl;
		//iter part setting b vector
		//----------------------------------------
		for (auto face : faces)
		{
			auto cell_volume = face->GetArea();
			auto right_hand_term_us = face->u * parameter.density * cell_volume / parameter.dt;
			auto right_hand_term_vs = face->v * parameter.density * cell_volume / parameter.dt;
			auto& half_edges = face->GetHalfEdge();

			for (auto half_edge : half_edges)
			{
				auto neighbour_half_edge = half_edge->GetNeighbourHalfEdge();
				auto neighbour_face = neighbour_half_edge->GetFace();

				if (neighbour_face != nullptr)
				{
					auto us_f = (face->us * face->GetArea() + neighbour_face->us * neighbour_face->GetArea()) / (face->GetArea() + neighbour_face->GetArea());
					auto vs_f = (face->vs * face->GetArea() + neighbour_face->vs * neighbour_face->GetArea()) / (face->GetArea() + neighbour_face->GetArea());

					right_hand_term_us -= (parameter.gamma * (us_f * half_edge->mass_flux) + (1.0 - parameter.gamma) *
						((std::max( half_edge->mass_flux, 0.0) * face->us          ) -
						 (std::max(-half_edge->mass_flux, 0.0) * neighbour_face->us)));
					right_hand_term_vs -= (parameter.gamma * (vs_f * half_edge->mass_flux) + (1.0 - parameter.gamma) *
						((std::max( half_edge->mass_flux, 0.0) * face->vs          ) -
						 (std::max(-half_edge->mass_flux, 0.0) * neighbour_face->vs)));
				}
			}

			lis_vector_set_value(LIS_INS_VALUE, face->GetId(), right_hand_term_us, b_us);
			lis_vector_set_value(LIS_INS_VALUE, face->GetId(), right_hand_term_vs, b_vs);
		}

		for (auto physical_group : physical_groups)
		{
			auto bc = physical_group.second->GetBoundaryCondition();

			if (bc.GetType() == INFLOW)
			{
				double value_u = bc.GetValue();
				double value_v = 0.0;//workaround
				auto boundary_edges = physical_group.second->GetEdges();
				for (auto boundary_edge : boundary_edges)
				{
					Face* face = nullptr;
					HalfEdge* he = nullptr;
					if (boundary_edge->GetHalfEdge(0)->GetFace() == nullptr)
					{
						he = boundary_edge->GetHalfEdge(1);
						face = he->GetFace();
					}
					else
					{
						he = boundary_edge->GetHalfEdge(0);
						face = he->GetFace();
					}

					auto cell_centeroid_to_wall_distance = abs((boundary_edge->GetCenter() - face->GetCentroid()) * (he->GetNormal()));

					auto face_area = he->GetEdgeVector().Abs();

					double right_hand_term_us;
					double right_hand_term_vs;
					lis_vector_get_value(b_us, face->GetId(), &right_hand_term_us);
					lis_vector_get_value(b_vs, face->GetId(), &right_hand_term_vs);

					right_hand_term_us += value_u * parameter.viscosity * face_area / cell_centeroid_to_wall_distance;
					right_hand_term_vs += value_v * parameter.viscosity * face_area / cell_centeroid_to_wall_distance;

					right_hand_term_us -= he->mass_flux * value_u;
					right_hand_term_vs -= he->mass_flux * value_v;

					lis_vector_set_value(LIS_INS_VALUE, face->GetId(), right_hand_term_us, b_us);
					lis_vector_set_value(LIS_INS_VALUE, face->GetId(), right_hand_term_vs, b_vs);
				}
			}
			if (bc.GetType() == OUTFLOW)
			{
				auto boundary_edges = physical_group.second->GetEdges();
				for (auto boundary_edge : boundary_edges)
				{
					Face* face = nullptr;
					HalfEdge* half_edge = nullptr;
					if (boundary_edge->GetHalfEdge(0)->GetFace() == nullptr)
					{
						half_edge = boundary_edge->GetHalfEdge(1);
						face = half_edge->GetFace();
					}
					else
					{
						half_edge = boundary_edge->GetHalfEdge(0);
						face = half_edge->GetFace();
					}
					double right_hand_term_us;
					double right_hand_term_vs;
					lis_vector_get_value(b_us, face->GetId(), &right_hand_term_us);
					lis_vector_get_value(b_vs, face->GetId(), &right_hand_term_vs);

					auto us_f = face->us;
					auto vs_f = face->vs;

					right_hand_term_us -= (us_f * half_edge->mass_flux);
					right_hand_term_vs -= (vs_f * half_edge->mass_flux);
					lis_vector_set_value(LIS_INS_VALUE, face->GetId(), right_hand_term_us, b_us);
					lis_vector_set_value(LIS_INS_VALUE, face->GetId(), right_hand_term_vs, b_vs);
				}
			}
		}


		//solving (*) velocity
		//----------------------------------------
		lis_solve(A_us, b_us, x, solver_momentum);
		for (auto face : faces)
		{
			double us_value;
			lis_vector_get_value(x, face->GetId(), &us_value);
			face->us = us_value;
		}

		lis_solve(A_us, b_vs, x, solver_momentum);
		for (auto face : faces)
		{
			double vs_value;
			lis_vector_get_value(x, face->GetId(), &vs_value);
			face->vs = vs_value;
		}

		//update mass flux
		//-------------------------------------------------------
		for (auto edge : geometry.GetEdgeList())
		{
			HalfEdge* half_edge[2]{ nullptr, nullptr };
			Face* face[2]{ nullptr, nullptr };
			half_edge[0] = edge->GetHalfEdge(0);
			half_edge[1] = edge->GetHalfEdge(1);
			face[0] = half_edge[0]->GetFace();
			face[1] = half_edge[1]->GetFace();

			auto us_f = (face[0]->us * face[1]->GetArea() + face[1]->us * face[0]->GetArea()) / (face[0]->GetArea() + face[1]->GetArea());
			auto vs_f = (face[0]->vs * face[1]->GetArea() + face[1]->vs * face[0]->GetArea()) / (face[0]->GetArea() + face[1]->GetArea());
			half_edge[0]->mass_flux = parameter.density * (us_f * half_edge[0]->GetEdgeVector().Abs() * half_edge[0]->GetNormal().GetDx()) +
				parameter.density * (vs_f * half_edge[0]->GetEdgeVector().Abs() * half_edge[0]->GetNormal().GetDy());

			half_edge[1]->mass_flux = -half_edge[0]->mass_flux;
		}

		for (auto physical_group : physical_groups)
		{
			auto bc = physical_group.second->GetBoundaryCondition();

			if (bc.GetType() == INFLOW)
			{
				double value_u = bc.GetValue();
				double value_v = 0.0;//workaround
				auto boundary_edges = physical_group.second->GetEdges();
				for (auto boundary_edge : boundary_edges)
				{
					Face* face = nullptr;
					HalfEdge* he = nullptr;
					if (boundary_edge->GetHalfEdge(0)->GetFace() == nullptr)
					{
						he = boundary_edge->GetHalfEdge(1);
					}
					else
					{
						he = boundary_edge->GetHalfEdge(0);
					}
					he->mass_flux = parameter.density * (value_u * he->GetEdgeVector().Abs() * he->GetNormal().GetDx()) +
						parameter.density * (value_v * he->GetEdgeVector().Abs() * he->GetNormal().GetDy());

				}
			}
			if (bc.GetType() == OUTFLOW)
			{
				auto boundary_edges = physical_group.second->GetEdges();
				for (auto boundary_edge : boundary_edges)
				{
					Face* face = nullptr;
					HalfEdge* he = nullptr;
					if (boundary_edge->GetHalfEdge(0)->GetFace() == nullptr)
					{
						he = boundary_edge->GetHalfEdge(1);
						face = he->GetFace();
					}
					else
					{
						he = boundary_edge->GetHalfEdge(0);
						face = he->GetFace();
					}

					double value_u = face->us;
					double value_v = face->vs;

					he->mass_flux = parameter.density * (value_u * he->GetEdgeVector().Abs() * he->GetNormal().GetDx()) +
						parameter.density * (value_v * he->GetEdgeVector().Abs() * he->GetNormal().GetDy());
				}
			}
			if (bc.GetType() == WALL)
			{
				auto boundary_edges = physical_group.second->GetEdges();
				for (auto boundary_edge : boundary_edges)
				{
					Face* face = nullptr;
					HalfEdge* he = nullptr;
					if (boundary_edge->GetHalfEdge(0)->GetFace() == nullptr)
					{
						he = boundary_edge->GetHalfEdge(1);
					}
					else
					{
						he = boundary_edge->GetHalfEdge(0);
					}

					he->mass_flux = 0.0;
				}
			}
		}

		for (auto face : faces)
		{
			auto& half_edges = face->GetHalfEdge();

			double mass_flux = 0.0;
			for (auto half_edge : half_edges)
			{
				mass_flux += half_edge->mass_flux;
			}
			auto val = mass_flux / parameter.dt;
			lis_vector_set_value(LIS_INS_VALUE, face->GetId(), val, b_p);
		}

		lis_solve(A_p, b_p, x, solver_pressure);


		double value;
		lis_vector_get_value(x, 0, &value);

		for (auto face : faces)
		{
			double ps_value;
			lis_vector_get_value(x, face->GetId(), &ps_value);
			face->ps = ps_value - value;
		}

		//updated mass flux
		for (auto edge : geometry.GetEdgeList())
		{
			HalfEdge* half_edge[2]{ nullptr, nullptr };
			Face* face[2]{ nullptr, nullptr };
			half_edge[0] = edge->GetHalfEdge(0);
			half_edge[1] = edge->GetHalfEdge(1);
			face[0] = half_edge[0]->GetFace();
			face[1] = half_edge[1]->GetFace();

			auto face_area = half_edge[0]->GetEdgeVector().Abs();
			auto cell_centeroid_distance = (face[0]->GetCentroid() - face[1]->GetCentroid()).Abs();

			half_edge[0]->mass_flux -= parameter.dt * face_area * (face[1]->ps - face[0]->ps) / cell_centeroid_distance;
			half_edge[1]->mass_flux = -half_edge[0]->mass_flux;
		}
	}

	//getting pressure
	for (auto face : faces)
	{
		face->p = face->ps;
	}

	// setting b_u vector
	//----------------------------------------
	for (auto face : faces)
	{
		auto cell_volume = face->GetArea();
		auto right_hand_term_u = face->u * parameter.density * cell_volume / parameter.dt;
		auto right_hand_term_v = face->v * parameter.density * cell_volume / parameter.dt;
		auto& half_edges = face->GetHalfEdge();

		for (auto half_edge : half_edges)
		{
			auto neighbour_half_edge = half_edge->GetNeighbourHalfEdge();
			auto neighbour_face = neighbour_half_edge->GetFace();

			if (neighbour_face != nullptr)
			{
				auto p_f = (face->p * neighbour_face->GetArea() + neighbour_face->us * face->GetArea()) / (face->GetArea() + neighbour_face->GetArea());
				
				right_hand_term_u -= (p_f * half_edge->GetEdgeVector().Abs() * half_edge->GetNormal().GetDx());
				right_hand_term_v -= (p_f * half_edge->GetEdgeVector().Abs() * half_edge->GetNormal().GetDy());
			}
		}

		lis_vector_set_value(LIS_INS_VALUE, face->GetId(), right_hand_term_u, b_u);
		lis_vector_set_value(LIS_INS_VALUE, face->GetId(), right_hand_term_v, b_v);
	}

	for (auto physical_group : physical_groups)
	{
		auto bc = physical_group.second->GetBoundaryCondition();

		if (bc.GetType() == INFLOW)
		{
			double value_u = bc.GetValue();
			double value_v = 0.0;//workaround
			auto boundary_edges = physical_group.second->GetEdges();
			for (auto boundary_edge : boundary_edges)
			{
				Face* face = nullptr;
				HalfEdge* he = nullptr;
				if (boundary_edge->GetHalfEdge(0)->GetFace() == nullptr)
				{
					he = boundary_edge->GetHalfEdge(1);
					face = he->GetFace();
				}
				else
				{
					he = boundary_edge->GetHalfEdge(0);
					face = he->GetFace();
				}

				auto cell_centeroid_to_wall_distance = abs((boundary_edge->GetCenter() - face->GetCentroid()) * (he->GetNormal()));

				auto face_area = he->GetEdgeVector().Abs();

				double right_hand_term_u;
				double right_hand_term_v;
				lis_vector_get_value(b_us, face->GetId(), &right_hand_term_u);
				lis_vector_get_value(b_vs, face->GetId(), &right_hand_term_v);

				auto p_f = face->p;

				right_hand_term_u -= (p_f * he->GetEdgeVector().Abs() * he->GetNormal().GetDx());
				right_hand_term_v -= (p_f * he->GetEdgeVector().Abs() * he->GetNormal().GetDy());

				right_hand_term_u += value_u * parameter.viscosity * face_area / cell_centeroid_to_wall_distance;
				right_hand_term_v += value_v * parameter.viscosity * face_area / cell_centeroid_to_wall_distance;

				right_hand_term_u -= he->mass_flux * value_u;
				right_hand_term_v -= he->mass_flux * value_v;

				lis_vector_set_value(LIS_INS_VALUE, face->GetId(), right_hand_term_u, b_u);
				lis_vector_set_value(LIS_INS_VALUE, face->GetId(), right_hand_term_v, b_v);
			}
		}
		if (bc.GetType() == OUTFLOW)
		{
			auto boundary_edges = physical_group.second->GetEdges();
			for (auto boundary_edge : boundary_edges)
			{
				Face* face = nullptr;
				HalfEdge* half_edge = nullptr;
				if (boundary_edge->GetHalfEdge(0)->GetFace() == nullptr)
				{
					half_edge = boundary_edge->GetHalfEdge(1);
					face = half_edge->GetFace();
				}
				else
				{
					half_edge = boundary_edge->GetHalfEdge(0);
					face = half_edge->GetFace();
				}
				double right_hand_term_u;
				double right_hand_term_v;
				lis_vector_get_value(b_us, face->GetId(), &right_hand_term_u);
				lis_vector_get_value(b_vs, face->GetId(), &right_hand_term_v);

				auto p_f = face->p;

				right_hand_term_u -= (p_f * half_edge->GetEdgeVector().Abs() * half_edge->GetNormal().GetDx());
				right_hand_term_v -= (p_f * half_edge->GetEdgeVector().Abs() * half_edge->GetNormal().GetDy());

				lis_vector_set_value(LIS_INS_VALUE, face->GetId(), right_hand_term_u, b_u);
				lis_vector_set_value(LIS_INS_VALUE, face->GetId(), right_hand_term_v, b_v);
			}
		}
		if (bc.GetType() == WALL)
		{
			auto boundary_edges = physical_group.second->GetEdges();
			for (auto boundary_edge : boundary_edges)
			{
				Face* face = nullptr;
				HalfEdge* half_edge = nullptr;
				if (boundary_edge->GetHalfEdge(0)->GetFace() == nullptr)
				{
					half_edge = boundary_edge->GetHalfEdge(1);
					face = half_edge->GetFace();
				}
				else
				{
					half_edge = boundary_edge->GetHalfEdge(0);
					face = half_edge->GetFace();
				}
				double right_hand_term_u;
				double right_hand_term_v;
				lis_vector_get_value(b_us, face->GetId(), &right_hand_term_u);
				lis_vector_get_value(b_vs, face->GetId(), &right_hand_term_v);

				auto p_f = face->p;

				right_hand_term_u -= (p_f * half_edge->GetEdgeVector().Abs() * half_edge->GetNormal().GetDx());
				right_hand_term_v -= (p_f * half_edge->GetEdgeVector().Abs() * half_edge->GetNormal().GetDy());

				lis_vector_set_value(LIS_INS_VALUE, face->GetId(), right_hand_term_u, b_u);
				lis_vector_set_value(LIS_INS_VALUE, face->GetId(), right_hand_term_v, b_v);
			}
		}
	}


	//solving velocity
	//----------------------------------------
	lis_solve(A_u, b_u, x, solver_velocity);
	for (auto face : faces)
	{
		double u_value;
		lis_vector_get_value(x, face->GetId(), &u_value);
		face->u = u_value;
	}

	lis_solve(A_u, b_v, x, solver_velocity);
	for (auto face : faces)
	{
		double v_value;
		lis_vector_get_value(x, face->GetId(), &v_value);
		face->v = v_value;
	}

}

void Solution::TestFeature(Data& data)
{
	//auto& temperature = data.GetFaceData(TEMPERATURE);
	EdgeDataHandler edge_handler(data);
	FaceDataHandler face_handler(data);

	auto& edge_temperature_data = data.GetEdgeData(TEMPERATURE);
	auto& face_temperature_data = data.GetFaceData(TEMPERATURE);
	auto& face_temperature_grad_data = data.GetFaceGradData(TEMPERATURE);

	auto physical_groups = geometry.GetPhysicalGroup();

	auto& faces = geometry.GetFaceList();
	for (auto& face : faces)
	{
		auto face_cg = face->GetCentroid();
		face_temperature_data[face->GetId()] = face_cg[0] * 7.0;
	}

	edge_handler.Update(TEMPERATURE);

	//BoundaryCondition
	for (auto physical_group : physical_groups)
	{
		auto bc = physical_group.second->GetBoundaryCondition();

		if (bc.GetType() == WALL)
		{
			auto boundary_edges = physical_group.second->GetEdges();
			for (auto boundary_edge : boundary_edges)
			{
				auto boundary_edge_id = boundary_edge->GetId();
				edge_temperature_data[boundary_edge_id] = 0;
			}
		}

		if (bc.GetType() == INFLOW)
		{
			auto boundary_edges = physical_group.second->GetEdges();
			for (auto boundary_edge : boundary_edges)
			{
				auto boundary_edge_id = boundary_edge->GetId();
				edge_temperature_data[boundary_edge_id] = 100;
			}
		}

		if (bc.GetType() == OUTFLOW)
		{
			auto boundary_edges = physical_group.second->GetEdges();
			for (auto boundary_edge : boundary_edges)
			{
				auto boundary_edge_id = boundary_edge->GetId();
				edge_temperature_data[boundary_edge_id] = 200;
			}
		}
	}
	face_handler.UpdateGradient(TEMPERATURE);











	//------------------------------------------------------------------------
	auto vertex_list = geometry.GetVertexList();
	auto face_list = geometry.GetFaceList();
	std::ofstream myfile;
	auto file = "temperature.dat";
	myfile.open(file);
	myfile << "TITLE = \"title\"\n";
	myfile << "VARIABLES = \"X\", \"Y\", \"Temperature\", \"TemperatureBack\", \"TemperatureGrad\"\n";
	myfile << "ZONE N = " << vertex_list.size() << ", E = " << face_list.size() << ", DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE\n";
	myfile << "VARLOCATION = ([3,4,5] = CELLCENTERED)" << "\n";
	for (auto it : vertex_list)
	{
		auto coord = it->GetPositionVector();
		myfile << coord[0] << " ";
	}
	myfile << "\n";
	for (auto it : vertex_list)
	{
		auto coord = it->GetPositionVector();
		myfile << coord[1] << "\n";
	}
	myfile << "\n";
	for (auto it : face_list)
	{
		auto face_id = it->GetId();
		myfile << face_temperature_data[face_id] << "\n";
	}
	myfile << "\n";
	for (auto it : face_list)
	{
		auto face_id = it->GetId();
		auto half_edges = it->GetHalfEdge();
		double temperature = 0;
		for (auto& half_edge : half_edges)
		{
			auto edge_id = half_edge->GetParentEdge()->GetId();
			temperature += edge_temperature_data[edge_id] / 3.0;
		}
		myfile << temperature << "\n";
	}
	myfile << "\n";
	for (auto it : face_list)
	{
		auto face_id = it->GetId();
		myfile << face_temperature_grad_data[face_id][0] << "\n";
	}
	myfile << "\n";

	for (auto it : face_list)
	{
		myfile << it->GetHalfEdge()[0]->GetStart()->GetId() << " " << it->GetHalfEdge()[1]->GetStart()->GetId() << " " << it->GetHalfEdge()[2]->GetStart()->GetId() << "\n";
	}

	myfile.close();
	//------------------------------------------------------------------------
}

void Solution::PoissonSolver()
{

}
