#include<iostream>
#include<Solution.hpp>
#include<set>
#include<iostream>
#include<chrono>
#include<lis.h>
int main()
{
    std::string name ="test";
    std::string file_location = "D:/Project/uFVM/input/";
    std::string file_name = name + ".msh";

    Solution solution;
    solution.ReadGeometry(file_location + file_name);
    
    LIS_MATRIX A;
    LIS_VECTOR b, x;
    LIS_SOLVER solver;
    
    LIS_INT i, n;
    n = 199;

    lis_vector_create(0, &b);
    lis_vector_create(0, &x);
    lis_vector_set_size(b, 0, n); /* or lis_vector_set_size(v,n,0); */
    lis_vector_set_size(x, 0, n); /* or lis_vector_set_size(v,n,0); */

    for (i = 0; i < n; i++)
    {
        lis_vector_set_value(LIS_INS_VALUE, i, 0, b);
        lis_vector_set_value(LIS_INS_VALUE, i, 0, x);
    }
    lis_vector_set_value(LIS_INS_VALUE, n-1, -100, b);
    
    lis_matrix_create(0, &A);
    lis_matrix_set_size(A, 0, n); /* or lis_matrix_set_size(A,n,0); */

    for (i = 0; i < n; i++)
    {
        if (i > 0) lis_matrix_set_value(LIS_INS_VALUE, i, i - 1, 1.0, A);
        if (i < n - 1) lis_matrix_set_value(LIS_INS_VALUE, i, i + 1, 1.0, A);
        lis_matrix_set_value(LIS_INS_VALUE, i, i, -2.0, A);
    }
    lis_matrix_set_type(A, LIS_MATRIX_CSR);
    lis_matrix_assemble(A);
    
    lis_solver_create(&solver);
    lis_solver_set_option("-i bicg -p none", solver);
    lis_solver_set_option("-tol 1.0e-12", solver);
    lis_solve(A, b, x, solver);


    LIS_SCALAR val;
    for (i = 0; i < n; i++)
    {
        lis_vector_get_value(x, i, &val);
        std::cout << "val: " << val << std::endl;
    }

    return 0;
}