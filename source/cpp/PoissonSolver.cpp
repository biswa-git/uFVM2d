#pragma once
#include<PoissonSolver.hpp>

/*
PoissonSolver::PoissonSolver(Geometry& geometry): m_geometry(&geometry)
{
	auto& faces = geometry.GetFaceList();

	lis_matrix_create(0, &A);
	lis_matrix_set_size(A, 0, faces.size());

	lis_vector_create(0, &b);
	lis_vector_set_size(b, 0, faces.size());
	lis_vector_create(0, &x);
	lis_vector_set_size(x, 0, faces.size());

	lis_solver_create(&solver);
	lis_solver_set_option("-i bicgstab -p ilut", solver);
	lis_solver_set_option("-tol 1.0e-15", solver);
	//setting A mat
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
				lis_matrix_set_value(LIS_INS_VALUE, face->GetId(), neighbour_face->GetId(), diagonal_term, A);
			}
		}
		lis_matrix_set_value(LIS_INS_VALUE, face->GetId(), face->GetId(), face->central_term, A);
	}


	auto physical_groups = geometry.GetPhysicalGroup();

	for (auto& physical_group : physical_groups)
	{
		auto boundary_condition = physical_group.second->GetBoundaryCondition();
		if (boundary_condition.GetType() == INFLOW)
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
				auto coefficient = face_area / cell_centeroid_to_wall_distance;
				face->central_term -= coefficient;
				lis_matrix_set_value(LIS_INS_VALUE, face->GetId(), face->GetId(), face->central_term, A);
				lis_vector_set_value(LIS_INS_VALUE, face->GetId(), - coefficient * boundary_condition.GetValue(), b);	
			}
		}
	}

	lis_matrix_set_type(A, LIS_MATRIX_CSR);
	lis_matrix_assemble(A);

	lis_solve(A, b, x, solver);

}
*/

PoissonSolver::PoissonSolver()
{
	lis_solver_create(&m_solver);
	lis_solver_set_option("-i bicgstab -p ilut", m_solver);
	lis_solver_set_option("-tol 1.0e-15", m_solver);
}

PoissonSolver::~PoissonSolver()
{
	lis_solver_destroy(m_solver);
	lis_matrix_destroy(m_A);
	lis_vector_destroy(m_b);
	lis_vector_destroy(m_x);
}


void PoissonSolver::SetSize(const LIS_INT& size)
{
	m_size = size;
	m_temp_b_vector.resize(m_size);

	lis_matrix_create(0, &m_A);
	lis_matrix_set_size(m_A, 0, m_size);

	lis_vector_create(0, &m_b);
	lis_vector_set_size(m_b, 0, m_size);
	lis_vector_create(0, &m_x);
	lis_vector_set_size(m_x, 0, m_size);
}

void PoissonSolver::MatrixSetValue(const int& i, const int& j, const double& value)
{
	lis_matrix_set_value(LIS_INS_VALUE, i, j, value, m_A);
}

void PoissonSolver::BVectorSetValue(const int& i, const double& value)
{
	m_temp_b_vector[i] = value;
}

double PoissonSolver::BVectorGetValue(const int& i)
{
	return m_temp_b_vector[i];
}

void PoissonSolver::XVectorSetValue(const int& i, const double& value)
{
	lis_vector_set_value(LIS_INS_VALUE, i, value, m_x);
}

void PoissonSolver::PrepareMatrix()
{
	lis_matrix_set_type(m_A, LIS_MATRIX_CSR);
	lis_matrix_assemble(m_A);
}

void PoissonSolver::PrepareBVector()
{
	for (int i = 0; i < m_size; ++i)
	{
		lis_vector_set_value(LIS_INS_VALUE, i, m_temp_b_vector[i], m_b);
	}
}
LIS_VECTOR& PoissonSolver::GetResult()
{
	lis_solve(m_A, m_b, m_x, m_solver);
	return m_x;
}