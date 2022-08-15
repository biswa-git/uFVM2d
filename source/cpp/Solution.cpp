#include<Solution.hpp>

Solution::Solution()
{
	m_edge_face_map[E_VEL_XS] = F_VEL_XS;
	m_edge_face_map[E_VEL_YS] = F_VEL_YS;
	m_edge_face_map[E_TEMPERATURE] = F_TEMPERATURE;

	m_face_grad_edge_map[FG_TEMPERATURE] = E_TEMPERATURE;
}

Solution::~Solution()
{
}

GeometryResult Solution::ReadGeometry(const std::string& file_name)
{
	auto result = m_geometry.Read(file_name);
	return result;
}

GeometryResult Solution::WriteSolution(const std::string& file_name)
{
	GeometryResult result;
	auto vertex_list = m_geometry.GetVertexList();
	auto face_list = m_geometry.GetFaceList();
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
		//myfile << it->GetFaceData(F_VEL_X) << "\n";
		double val = 0.0;
		for (auto edge : it->GetHalfEdge())
		{
			val += edge->GetParentEdge()->GetEdgeData(E_MASS_FLUX);
		}
		myfile << val << "\n";
	}
	myfile << "\n";

	for (auto it : face_list)
	{
		myfile << it->GetFaceData(F_VEL_Y) << "\n";
	}
	myfile << "\n";

	for (auto it : face_list)
	{
		myfile << it->GetFaceData(F_PRESSURE) << "\n";
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
	BoundaryCondition boundary_condition;
	boundary_condition.SetBoundaryCondition(type, value);
	m_geometry.GetPhysicalgroupByName(name)->SetBoundaryCondition(boundary_condition);
}

void Solution::SetTimeStep(const double& dt)
{
	m_parameter.dt = dt;
}

void Solution::SetDensity(const double& density)
{
	m_parameter.density = density;
}

void Solution::SetViscosity(const double& viscosity)
{
	m_parameter.viscosity = viscosity;
}

void Solution::SetMaxTimeStep(const int& max_timestep)
{
	m_parameter.max_timestep = max_timestep;
}

Geometry& Solution::GetGeometry()
{
	return m_geometry;
}
/*
void Solution::SolveMomentum()
{
	auto& faces = m_geometry.GetFaceList();

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


	auto physical_groups = m_geometry.GetPhysicalGroup();

	//setting A mat
	//----------------------------------------------------------------------------------------
	for (auto face : faces)
	{
		auto cell_volume = face->GetArea();
		face->central_term = m_parameter.density * cell_volume / m_parameter.dt;
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
				auto coefficient = m_parameter.viscosity * face_area / cell_centeroid_distance;
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
				auto coefficient =m_parameter.viscosity * face_area / cell_centeroid_to_wall_distance;
				face->central_term += coefficient;
				lis_matrix_set_value(LIS_INS_VALUE, face->GetId(), face->GetId(), face->central_term, A_us);


				if (bc.GetType() == INFLOW)
				{
					he->mass_flux =m_parameter.density * (bc.GetValue() * he->GetEdgeVector().Abs() * he->GetNormal().GetDx());
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
		face->central_term =m_parameter.density * cell_volume /m_parameter.dt;
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
				auto coefficient =m_parameter.viscosity * face_area / cell_centeroid_distance;
				face->central_term += coefficient;
				face->central_term += (1.0 -m_parameter.gamma) * std::max(half_edge->mass_flux, 0.0);
				face->central_term +=m_parameter.gamma * ((half_edge->mass_flux*neighbour_face->GetArea())/
					                                    (face->GetArea()+ neighbour_face->GetArea()));

				diagonal_term = -coefficient;
				diagonal_term -= (1.0 -m_parameter.gamma) * std::max(-half_edge->mass_flux, 0.0);
				diagonal_term +=m_parameter.gamma * ((half_edge->mass_flux * face->GetArea()) /
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
				auto coefficient =m_parameter.viscosity * face_area / cell_centeroid_to_wall_distance;
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
			auto right_hand_term_us = face->u *m_parameter.density * cell_volume /m_parameter.dt;
			auto right_hand_term_vs = face->v *m_parameter.density * cell_volume /m_parameter.dt;
			auto& half_edges = face->GetHalfEdge();

			for (auto half_edge : half_edges)
			{
				auto neighbour_half_edge = half_edge->GetNeighbourHalfEdge();
				auto neighbour_face = neighbour_half_edge->GetFace();

				if (neighbour_face != nullptr)
				{
					auto us_f = (face->us * face->GetArea() + neighbour_face->us * neighbour_face->GetArea()) / (face->GetArea() + neighbour_face->GetArea());
					auto vs_f = (face->vs * face->GetArea() + neighbour_face->vs * neighbour_face->GetArea()) / (face->GetArea() + neighbour_face->GetArea());

					right_hand_term_us -= (m_parameter.gamma * (us_f * half_edge->mass_flux) + (1.0 -m_parameter.gamma) *
						((std::max( half_edge->mass_flux, 0.0) * face->us          ) -
						 (std::max(-half_edge->mass_flux, 0.0) * neighbour_face->us)));
					right_hand_term_vs -= (m_parameter.gamma * (vs_f * half_edge->mass_flux) + (1.0 -m_parameter.gamma) *
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

					right_hand_term_us += value_u *m_parameter.viscosity * face_area / cell_centeroid_to_wall_distance;
					right_hand_term_vs += value_v *m_parameter.viscosity * face_area / cell_centeroid_to_wall_distance;

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
		for (auto edge : m_geometry.GetEdgeList())
		{
			HalfEdge* half_edge[2]{ nullptr, nullptr };
			Face* face[2]{ nullptr, nullptr };
			half_edge[0] = edge->GetHalfEdge(0);
			half_edge[1] = edge->GetHalfEdge(1);
			face[0] = half_edge[0]->GetFace();
			face[1] = half_edge[1]->GetFace();

			auto us_f = (face[0]->us * face[1]->GetArea() + face[1]->us * face[0]->GetArea()) / (face[0]->GetArea() + face[1]->GetArea());
			auto vs_f = (face[0]->vs * face[1]->GetArea() + face[1]->vs * face[0]->GetArea()) / (face[0]->GetArea() + face[1]->GetArea());
			half_edge[0]->mass_flux =m_parameter.density * (us_f * half_edge[0]->GetEdgeVector().Abs() * half_edge[0]->GetNormal().GetDx()) +
				m_parameter.density * (vs_f * half_edge[0]->GetEdgeVector().Abs() * half_edge[0]->GetNormal().GetDy());

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
					he->mass_flux =m_parameter.density * (value_u * he->GetEdgeVector().Abs() * he->GetNormal().GetDx()) +
						m_parameter.density * (value_v * he->GetEdgeVector().Abs() * he->GetNormal().GetDy());

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

					he->mass_flux =m_parameter.density * (value_u * he->GetEdgeVector().Abs() * he->GetNormal().GetDx()) +
						m_parameter.density * (value_v * he->GetEdgeVector().Abs() * he->GetNormal().GetDy());
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
			auto val = mass_flux /m_parameter.dt;
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
		for (auto edge : m_geometry.GetEdgeList())
		{
			HalfEdge* half_edge[2]{ nullptr, nullptr };
			Face* face[2]{ nullptr, nullptr };
			half_edge[0] = edge->GetHalfEdge(0);
			half_edge[1] = edge->GetHalfEdge(1);
			face[0] = half_edge[0]->GetFace();
			face[1] = half_edge[1]->GetFace();

			auto face_area = half_edge[0]->GetEdgeVector().Abs();
			auto cell_centeroid_distance = (face[0]->GetCentroid() - face[1]->GetCentroid()).Abs();

			half_edge[0]->mass_flux -=m_parameter.dt * face_area * (face[1]->ps - face[0]->ps) / cell_centeroid_distance;
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
		auto right_hand_term_u = face->u *m_parameter.density * cell_volume /m_parameter.dt;
		auto right_hand_term_v = face->v *m_parameter.density * cell_volume /m_parameter.dt;
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

				right_hand_term_u += value_u *m_parameter.viscosity * face_area / cell_centeroid_to_wall_distance;
				right_hand_term_v += value_v *m_parameter.viscosity * face_area / cell_centeroid_to_wall_distance;

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
*/

void Solution::TestFeature()
{
	auto& physical_groups = m_geometry.GetPhysicalGroup();

	auto& faces = m_geometry.GetFaceList();
	for (auto& face : faces)
	{
		auto face_cg = face->GetCentroid();
		face->GetFaceData(F_TEMPERATURE) = face_cg[0] * 7.0;
		face->GetFaceData(F_VEL_XS) = face_cg[0];
		face->GetFaceData(F_VEL_YS) = face_cg[1];
	}

	UpdateEdgeData(E_TEMPERATURE);

	//BoundaryCondition
	for (auto& physical_group : physical_groups)
	{
		auto boundary_condition = physical_group.second->GetBoundaryCondition();

		if (boundary_condition.GetType() == WALL)
		{
			auto& boundary_edges = physical_group.second->GetEdges();
			for (auto boundary_edge : boundary_edges)
			{
				auto boundary_edge_id = boundary_edge->GetId();
				boundary_edge->GetEdgeData(E_TEMPERATURE) = 0;
			}
		}

		if (boundary_condition.GetType() == INFLOW)
		{
			auto boundary_edges = physical_group.second->GetEdges();
			for (auto boundary_edge : boundary_edges)
			{
				auto boundary_edge_id = boundary_edge->GetId();
				boundary_edge->GetEdgeData(E_TEMPERATURE) = 100;
			}
		}

		if (boundary_condition.GetType() == OUTFLOW)
		{
			auto boundary_edges = physical_group.second->GetEdges();
			for (auto boundary_edge : boundary_edges)
			{
				auto boundary_edge_id = boundary_edge->GetId();
				boundary_edge->GetEdgeData(E_TEMPERATURE) = 200;
			}
		}
	}

	UpdateFaceGradient(FG_TEMPERATURE);

	UpdateFlux(E_VEL_XS, E_VEL_YS, E_MASS_FLUX_S);








	//------------------------------------------------------------------------
	auto vertex_list = m_geometry.GetVertexList();
	auto face_list = m_geometry.GetFaceList();
	std::ofstream myfile;
	auto file = "temperature.dat";
	myfile.open(file);
	myfile << "TITLE = \"title\"\n";
	myfile << "VARIABLES = \"X\", \"Y\", \"Temperature\", \"TemperatureBack\", \"TemperatureGrad\", \"US\", \"VS\", \"MassFlux\"\n";
	myfile << "ZONE N = " << vertex_list.size() << ", E = " << face_list.size() << ", DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE\n";
	myfile << "VARLOCATION = ([3,4,5,6,7,8] = CELLCENTERED)" << "\n";
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
		myfile << it->GetFaceData(F_TEMPERATURE) << "\n";
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
			temperature += half_edge->GetParentEdge()->GetEdgeData(E_TEMPERATURE) / 3.0;
		}
		myfile << temperature << "\n";
	}
	myfile << "\n";
	for (auto it : face_list)
	{
		auto face_id = it->GetId();
		myfile << it->GetFaceGradData(FG_TEMPERATURE)[0] << "\n";
	}
	myfile << "\n";
	for (auto it : face_list)
	{
		auto face_id = it->GetId();
		auto half_edges = it->GetHalfEdge();
		double u = 0;
		for (auto& half_edge : half_edges)
		{
			auto edge_id = half_edge->GetParentEdge()->GetId();
			u += half_edge->GetParentEdge()->GetEdgeData(E_VEL_XS)/3;
		}
		myfile << u << "\n";
	}
	myfile << "\n";
	for (auto it : face_list)
	{
		auto face_id = it->GetId();
		auto& half_edges = it->GetHalfEdge();
		double v = 0;
		for (auto half_edge : half_edges)
		{
			auto edge_id = half_edge->GetParentEdge()->GetId();
			v += half_edge->GetParentEdge()->GetEdgeData(E_VEL_YS) / 3;
		}
		myfile << v << "\n";
	}
	myfile << "\n";
	for (auto it : face_list)
	{
		
		auto face_id = it->GetId();
		auto half_edges = it->GetHalfEdge();
		double mass_flux = 0;
		for (auto& half_edge : half_edges)
		{
			auto edge_id = half_edge->GetParentEdge()->GetId();
			mass_flux += half_edge->GetParentEdge()->GetEdgeData(E_MASS_FLUX_S)* double(half_edge->GDirectionCoefficient()) / 2;
		}
		std::cout << it->GetArea() << " ------- " << mass_flux << std::endl;
		myfile << mass_flux << "\n";
	}
	myfile << "\n";

	for (auto it : face_list)
	{
		myfile << it->GetHalfEdge()[0]->GetStart()->GetId() << " " << it->GetHalfEdge()[1]->GetStart()->GetId() << " " << it->GetHalfEdge()[2]->GetStart()->GetId() << "\n";
	}

	myfile.close();
	//------------------------------------------------------------------------
}

void Solution::UpdateEdgeData(const int& edge_data_id)
{
	int face_data_id = m_edge_face_map[edge_data_id];
	std::vector<Face*> face(2, nullptr);

	for (auto edge : m_geometry.GetEdgeList())
	{
		face[0] = edge->GetHalfEdge(0)->GetFace();
		face[1] = edge->GetHalfEdge(1)->GetFace();

		if (face[0] != nullptr && face[1] != nullptr)
		{
			edge->GetEdgeData(edge_data_id) = (face[0]->GetFaceData(face_data_id) * face[1]->GetArea() +
				face[1]->GetFaceData(face_data_id) * face[0]->GetArea()) /
				(face[0]->GetArea() + face[1]->GetArea());
		}
	}
}

void Solution::UpdateFaceGradient(const int& face_grad_data_id)
{
	int edge_data_id = m_face_grad_edge_map[face_grad_data_id];
	for (auto& face : m_geometry.GetFaceList())
	{
		auto& half_edges = face->GetHalfEdge();
		auto volume = face->GetArea();

		face->GetFaceGradData(face_grad_data_id)[0] = 0.0;
		face->GetFaceGradData(face_grad_data_id)[1] = 0.0;
		face->GetFaceGradData(face_grad_data_id)[2] = 0.0;

		for (auto half_edge : half_edges)
		{
			auto edge = half_edge->GetParentEdge();
			auto edge_id = edge->GetId();

			auto area_vector = half_edge->GetNormal() * half_edge->GetEdgeVector().Abs();

			face->GetFaceGradData(face_grad_data_id)[0] += half_edge->GetParentEdge()->GetEdgeData(edge_data_id) * area_vector.GetDx();
			face->GetFaceGradData(face_grad_data_id)[1] += half_edge->GetParentEdge()->GetEdgeData(edge_data_id) * area_vector.GetDy();
			face->GetFaceGradData(face_grad_data_id)[2] += half_edge->GetParentEdge()->GetEdgeData(edge_data_id) * area_vector.GetDz();
		}
		face->GetFaceGradData(face_grad_data_id)[0] /= volume;
		face->GetFaceGradData(face_grad_data_id)[1] /= volume;
		face->GetFaceGradData(face_grad_data_id)[2] /= volume;
	}
}

void Solution::UpdateFlux(const int& e_vel_x_data_id, const int& e_vel_y_data_id, const int& e_flux_data_id)
{
	UpdateEdgeData(e_vel_x_data_id);
	UpdateEdgeData(e_vel_y_data_id);

	for (auto edge : m_geometry.GetEdgeList())
	{
		edge->GetEdgeData(e_flux_data_id) = m_parameter.density * edge->GetEdgeData(e_vel_x_data_id) * edge->GetHalfEdge(0)->GeAreaVector().GetDx() +
			                                m_parameter.density * edge->GetEdgeData(e_vel_y_data_id) * edge->GetHalfEdge(0)->GeAreaVector().GetDy();

	}

	//BOUNDARY_FLUX
	auto& physical_groups = m_geometry.GetPhysicalGroup();
	for (auto& physical_group : physical_groups)
	{
		auto boundary_condition = physical_group.second->GetBoundaryCondition();

		if (boundary_condition.GetType() == INFLOW)
		{
			auto& boundary_edges = physical_group.second->GetEdges();
			for (auto boundary_edge : boundary_edges)
			{
				HalfEdge* half_edge = nullptr;
				if (boundary_edge->GetHalfEdge(0)->GetFace() == nullptr)
				{
					half_edge = boundary_edge->GetHalfEdge(1);
				}
				else
				{
					half_edge = boundary_edge->GetHalfEdge(0);
				}

				double value_u = -boundary_condition.GetValue() * half_edge->GeAreaVector()[0];
				double value_v = -boundary_condition.GetValue() * half_edge->GeAreaVector()[1];

				half_edge->GetParentEdge()->GetEdgeData(e_flux_data_id) = m_parameter.density * value_u * half_edge->GeAreaVector().GetDx() +
					                                                      m_parameter.density * value_v * half_edge->GeAreaVector().GetDy();
			}
		}
		if (boundary_condition.GetType() == OUTFLOW)
		{
			auto& boundary_edges = physical_group.second->GetEdges();
			for (auto boundary_edge : boundary_edges)
			{
				HalfEdge* half_edge = nullptr;
				if (boundary_edge->GetHalfEdge(0)->GetFace() == nullptr)
				{
					half_edge = boundary_edge->GetHalfEdge(1);
				}
				else
				{
					half_edge = boundary_edge->GetHalfEdge(0);
				}

				Face* face = half_edge->GetFace();

				double value_u = face->GetFaceData(e_vel_x_data_id);
				double value_v = face->GetFaceData(e_vel_y_data_id);

				half_edge->GetParentEdge()->GetEdgeData(e_flux_data_id) = m_parameter.density * value_u * half_edge->GeAreaVector().GetDx() +
					                                                    m_parameter.density * value_v * half_edge->GeAreaVector().GetDy();
			}
		}

		if (boundary_condition.GetType() == WALL)
		{
			auto& boundary_edges = physical_group.second->GetEdges();
			for (auto boundary_edge : boundary_edges)
			{
				HalfEdge* half_edge = nullptr;
				if (boundary_edge->GetHalfEdge(0)->GetFace() == nullptr)
				{
					half_edge = boundary_edge->GetHalfEdge(1);
				}
				else
				{
					half_edge = boundary_edge->GetHalfEdge(0);
				}

				half_edge->GetParentEdge()->GetEdgeData(e_flux_data_id) = 0.0;
			}
		}
	}
}

void Solution::CopyEdgeData(const int& source, const int& destination)
{
	auto& edges = m_geometry.GetEdgeList();

	for (auto edge : edges)
	{
		edge->GetEdgeData(destination) = edge->GetEdgeData(source);
	}
}

void Solution::CopyFaceData(const int& source, const int& destination)
{
	auto& faces = m_geometry.GetFaceList();

	for (auto face : faces)
	{
		face->GetFaceData(destination) = face->GetFaceData(source);
	}
}

void Solution::CopyFaceGradData(const int& source, const int& destination)
{
	auto& faces = m_geometry.GetFaceList();

	for (auto face : faces)
	{
		face->GetFaceGradData(destination) = face->GetFaceGradData(source);
	}
}

void Solution::SetEdgeData(const int& data_id, const int& value)
{
	auto& edges = m_geometry.GetEdgeList();

	for (auto edge : edges)
	{
		edge->GetEdgeData(data_id) = value;
	}

	auto& physical_groups = m_geometry.GetPhysicalGroup();
	for (auto& physical_group : physical_groups)
	{
		auto boundary_condition = physical_group.second->GetBoundaryCondition();
		if (boundary_condition.GetType() == WALL || boundary_condition.GetType() == INFLOW)
		{
			auto& boundary_edges = physical_group.second->GetEdges();
			for (auto boundary_edge : boundary_edges)
			{
				boundary_edge->GetEdgeData(data_id) = value;
			}
		}
	}
}

void Solution::SetFaceData(const int& data_id, const int& value)
{
	auto& faces = m_geometry.GetFaceList();

	for (auto face : faces)
	{
		face->GetFaceData(data_id) = value;
	}
}

void Solution::ConstructMomentumMatrix()
{
	auto& faces = m_geometry.GetFaceList();
	m_momentum.SetSize(static_cast<LIS_INT>(faces.size()));
	//setting momentum mat
	//----------------------------------------------------------------------------------------
	for (auto face : faces)
	{
		auto cell_volume = face->GetArea();
		face->GetFaceData(F_CENTRAL_TERM) = m_parameter.density * cell_volume / m_parameter.dt;
		auto diagonal_term = 0.0;
		auto& half_edges = face->GetHalfEdge();

		for (auto half_edge : half_edges)
		{
			auto neighbour_half_edge = half_edge->GetNeighbourHalfEdge();
			auto neighbour_face = neighbour_half_edge->GetFace();

			if (neighbour_face != nullptr)
			{
				auto coefficient = m_parameter.viscosity * half_edge->GetParentEdge()->GetAreaByDistance();
				face->GetFaceData(F_CENTRAL_TERM) += coefficient;
				diagonal_term = -coefficient;
				m_momentum.MatrixSetValue(face->GetId(), neighbour_face->GetId(), diagonal_term);
			}
		}
		m_momentum.MatrixSetValue(face->GetId(), face->GetId(), face->GetFaceData(F_CENTRAL_TERM));
	}

	auto& physical_groups = m_geometry.GetPhysicalGroup();
	for (auto& physical_group : physical_groups)
	{
		auto boundary_condition = physical_group.second->GetBoundaryCondition();
		if (boundary_condition.GetType() == WALL || boundary_condition.GetType() == INFLOW)
		{
			auto& boundary_edges = physical_group.second->GetEdges();
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

				auto coefficient = m_parameter.viscosity * half_edge->GetParentEdge()->GetAreaByDistance();
				face->GetFaceData(F_CENTRAL_TERM) += coefficient;
				m_momentum.MatrixSetValue(face->GetId(), face->GetId(), face->GetFaceData(F_CENTRAL_TERM));
			}
		}

		if (boundary_condition.GetType() == OUTFLOW)
		{
			auto& boundary_edges = physical_group.second->GetEdges();
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
				face->GetFaceData(F_CENTRAL_TERM) += half_edge->GetParentEdge()->GetEdgeData(E_MASS_FLUX_S);
			}
		}
	}
	m_momentum.Prepare();
	//end setting momentum mat
	//----------------------------------------------------------------------------------------

}


void Solution::Solve()
{

	SetFaceData(F_VEL_X, 0);
	SetFaceData(F_VEL_Y, 0);
	SetFaceData(F_PRESSURE, 1);


	SetEdgeData(E_VEL_X, 0);
	SetEdgeData(E_VEL_Y, 0);

	auto& physical_groups = m_geometry.GetPhysicalGroup();

	//Initial codition
	//----------------------------------------------------------------------------------

	for (auto& physical_group : physical_groups)
	{
		auto boundary_conditon = physical_group.second->GetBoundaryCondition();

		if (boundary_conditon.GetType() == INFLOW)
		{
			auto& bundary_half_edges = physical_group.second->GetEdges();
			for (auto edge : bundary_half_edges)
			{
				edge->GetEdgeData(E_VEL_X) = boundary_conditon.GetValue();
				HalfEdge* half_edge = nullptr;
				if (edge->GetHalfEdge(0)->GetFace() == nullptr)
				{
					half_edge = edge->GetHalfEdge(1);
				}
				else
				{
					half_edge = edge->GetHalfEdge(0);
				}
				edge->GetEdgeData(E_MASS_FLUX) = -boundary_conditon.GetValue()* half_edge->GeAreaVector()[0];
			}
		}
	}
	//----------------------------------------------------------------------------------
	WriteSolution("newOutput1");




	//step 1

	CopyEdgeData(E_MASS_FLUX, E_MASS_FLUX_S);

	//step 2





}