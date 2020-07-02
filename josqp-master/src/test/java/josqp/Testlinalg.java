package josqp;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;

import com.quantego.josqp.CSCMatrix;
import com.quantego.josqp.DataGenerator_Lin_Alg;
import com.quantego.josqp.DataGenerator_Lin_Alg.lin_alg_sols_data;
import com.quantego.josqp.LinAlg;
import com.quantego.josqp.OSQP;
import com.quantego.josqp.OSQP.Data;
import com.quantego.josqp.OSQP.Settings;
import java.lang.Math;

public class Testlinalg{
    final static double TESTS_TOL = 0.0001;
    final static int OSQP_TIME_LIMIT_REACHED = -6;
    public double[] csc_to_dns(CSCMatrix matrix)
    {
        int i = 0;
        int j = 0;
        int index = 0;
        double[]A = new double[matrix.m * matrix.n];
        for ( index = 0; index < matrix.Ap[matrix.n]; ++index)
        {
            i = matrix.Ai[index];
            while(matrix.Ap[j+1] <= index)
            {
                j = j+1;
            }
            A[j* (matrix.m) + i] = matrix.Ax[index];
        }
        return A;


    };

    @Test
    public void test_constr_sparse_mat() {

    double []Adns; // Conversion to dense matrix
    lin_alg_sols_data data =  DataGenerator_Lin_Alg.generate_problem_lin_alg_sols_data();

    // Convert sparse to dense
    Adns = csc_to_dns(data.test_sp_matrix_A);

    // Compute norm of the elementwise difference with
    assertTrue((LinAlg.vec_norm_inf_diff(Adns, data.test_sp_matrix_Adns,data.test_sp_matrix_A.m*data.test_sp_matrix_A.n) < TESTS_TOL),"Linear algebra tests: error in constructing sparse/dense matrix!");
 
        
    };
    /*
    public void test_vec_operations() {
        double  norm_inf, vecprod; // normInf;
        double []ew_reciprocal;
        double []add_scaled;
        double []vec_ew_max_vec_test, vec_ew_min_vec_test;
      
        lin_alg_sols_data data = generate_problem_lin_alg_sols_data();
      
      
        // Add scaled
        add_scaled = vec_copy(data.test_vec_ops_v1, data.test_vec_ops_n);
        vec_add_scaled(add_scaled,
                       add_scaled,
                       data.test_vec_ops_v2,
                       data.test_vec_ops_n,
                       data.test_vec_ops_sc);
        assertTrue(
          "Linear algebra tests: error in vector operation, adding scaled vector",
          vec_norm_inf_diff(add_scaled, data.test_vec_ops_add_scaled,
                            data.test_vec_ops_n) < TESTS_TOL);
      
        // Norm_inf of the difference
        assertTrue(
          "Linear algebra tests: error in vector operation, norm_inf of difference",
          Math.abs(vec_norm_inf_diff(data.test_vec_ops_v1,
                                     data.test_vec_ops_v2,
                                     data.test_vec_ops_n) -
                   data.test_vec_ops_norm_inf_diff) <
          TESTS_TOL);
      
        // norm_inf
        norm_inf = vec_norm_inf(data.test_vec_ops_v1, data.test_vec_ops_n);
        assertTrue("Linear algebra tests: error in vector operation, norm_inf",
                  Math.abs(norm_inf - data.test_vec_ops_norm_inf) < TESTS_TOL);
      
        // Elementwise reciprocal
        ew_reciprocal = double[data.test_vec_ops_n];
        vec_ew_recipr(data.test_vec_ops_v1, ew_reciprocal, data.test_vec_ops_n);
        assertTrue(
          "Linear algebra tests: error in vector operation, elementwise reciprocal",
          vec_norm_inf_diff(ew_reciprocal, data.test_vec_ops_ew_reciprocal,data.test_vec_ops_n) < TESTS_TOL);
      
      
        // Vector product
        vecprod = vec_prod(data.test_vec_ops_v1,
                           data.test_vec_ops_v2,
                           data.test_vec_ops_n);
        assertTrue("Linear algebra tests: error in vector operation, vector product",
        Math.(vecprod - data.test_vec_ops_vec_prod) < TESTS_TOL);
      
        // Elementwise maximum between two vectors
        vec_ew_max_vec_test = double[data.test_vec_ops_n];
        vec_ew_max_vec(data.test_vec_ops_v1,
                       data.test_vec_ops_v2,
                       vec_ew_max_vec_test,
                       data.test_vec_ops_n);
        assertTrue(
          "Linear algebra tests: error in vector operation, elementwise maximum between vectors",
          vec_norm_inf_diff(vec_ew_max_vec_test, data.test_vec_ops_ew_max_vec,
                            data.test_vec_ops_n) < TESTS_TOL);  
        // Elementwise minimum between two vectors
        vec_ew_min_vec_test = new double[data.test_vec_ops_n ];
          //c_float *)c_malloc(data.test_vec_ops_n * sizeof(c_float));
        vec_ew_min_vec(data.test_vec_ops_v1,
                       data.test_vec_ops_v2,
                       vec_ew_min_vec_test,
                       data.test_vec_ops_n);
        assertTrue(
          "Linear algebra tests: error in vector operation, elementwise minimum between vectors",
          vec_norm_inf_diff(vec_ew_min_vec_test, data.test_vec_ops_ew_min_vec,
                            data.test_vec_ops_n) < TESTS_TOL);     
        return 0;
    };
    public void test_mat_operations() {
        CSCMatrix Ad, dA; // Matrices used for tests
        // csc *A_ewsq, *A_ewabs;     // Matrices used for tests
        int exitflag = 0;
      
        // c_float trace, fro_sq;
        double []inf_norm_cols_rows_test;
      
      
        lin_alg_sols_data data = generate_problem_lin_alg_sols_data();
      
      
        // Copy matrices
        Ad = copy_csc_mat(data.test_mat_ops_A);
        dA = copy_csc_mat(data.test_mat_ops_A);

      
      
        // Premultiply matrix A
        mat_premult_diag(dA, data.test_mat_ops_d);
        assertTrue(
          "Linear algebra tests: error in matrix operation, premultiply diagonal",
          is_eq_csc(dA, data.test_mat_ops_prem_diag, TESTS_TOL));
      
      
        // Postmultiply matrix A
        mat_postmult_diag(Ad, data.test_mat_ops_d);
        assertTrue(
          "Linear algebra tests: error in matrix operation, postmultiply diagonal",
          is_eq_csc(Ad, data.test_mat_ops_postm_diag, TESTS_TOL));
      
        // Maximum norm over columns
        inf_norm_cols_rows_test = double[data.test_mat_ops_n];
        mat_inf_norm_cols(data.test_mat_ops_A, inf_norm_cols_rows_test);
        assertTrue(
          "Linear algebra tests: error in matrix operation, max norm over columns",
          vec_norm_inf_diff(inf_norm_cols_rows_test, data.test_mat_ops_inf_norm_cols,
                            data.test_mat_ops_n) < TESTS_TOL);
      
        // Maximum norm over rows
        mat_inf_norm_rows(data.test_mat_ops_A, inf_norm_cols_rows_test);
        assertTrue("Linear algebra tests: error in matrix operation, max norm over rows",
                  vec_norm_inf_diff(inf_norm_cols_rows_test,
                                    data.test_mat_ops_inf_norm_rows,
                                    data.test_mat_ops_n) < TESTS_TOL);

        return 0;
      };
      public void test_mat_vec_multiplication() {
        double[]Ax, ATy, Px, Ax_cum, ATy_cum, Px_cum;
      
        lin_alg_sols_data data = generate_problem_lin_alg_sols_data();
      
      
        // Allocate vectors
        Ax  = new double[data.test_mat_vec_m];
        ATy = new double[data.test_mat_vec_n];
        Px  = new double[data.test_mat_vec_n ];
      
      
        // Matrix-vector multiplication:  y = Ax
        mat_vec(data.test_mat_vec_A, data.test_mat_vec_x, Ax, 0);
         assertTrue(
          "Linear algebra tests: error in matrix-vector operation, matrix-vector multiplication",
          vec_norm_inf_diff(Ax, data.test_mat_vec_Ax,
                            data.test_mat_vec_m) < TESTS_TOL);
      
        // Cumulative matrix-vector multiplication:  y += Ax
        Ax_cum = vec_copy(data.test_mat_vec_y, data.test_mat_vec_m);
        mat_vec(data.test_mat_vec_A, data.test_mat_vec_x, Ax_cum, 1);
         assertTrue(
          "Linear algebra tests: error in matrix-vector operation, cumulative matrix-vector multiplication",
          vec_norm_inf_diff(Ax_cum, data.test_mat_vec_Ax_cum,
                            data.test_mat_vec_m) < TESTS_TOL);
      
        // Matrix-transpose-vector multiplication:  x = A'*y
        mat_tpose_vec(data.test_mat_vec_A, data.test_mat_vec_y, ATy, 0, 0);
         assertTrue(
          "Linear algebra tests: error in matrix-vector operation, matrix-transpose-vector multiplication",
          vec_norm_inf_diff(ATy, data.test_mat_vec_ATy,
                            data.test_mat_vec_n) < TESTS_TOL);
      
        // Cumulative matrix-transpose-vector multiplication:  x += A'*y
        ATy_cum = vec_copy(data.test_mat_vec_x, data.test_mat_vec_n);
        mat_tpose_vec(data.test_mat_vec_A, data.test_mat_vec_y, ATy_cum, 1, 0);
         assertTrue(
          "Linear algebra tests: error in matrix-vector operation, cumulative matrix-transpose-vector multiplication",
          vec_norm_inf_diff(ATy_cum, data.test_mat_vec_ATy_cum,
                            data.test_mat_vec_n) < TESTS_TOL);
      
        // Symmetric-matrix-vector multiplication (only upper part is stored)
        mat_vec(data.test_mat_vec_Pu, data.test_mat_vec_x, Px, 0);          // upper
                                                                              // traingular
                                                                              // part
        mat_tpose_vec(data.test_mat_vec_Pu, data.test_mat_vec_x, Px, 1, 1); // lower
                                                                              // traingular
                                                                              // part
                                                                              // (without
                                                                              // diagonal)
         assertTrue(
          "Linear algebra tests: error in matrix-vector operation, symmetric matrix-vector multiplication",
          vec_norm_inf_diff(Px, data.test_mat_vec_Px,
                            data.test_mat_vec_n) < TESTS_TOL);
      
      
        // Cumulative symmetric-matrix-vector multiplication
        Px_cum = vec_copy(data.test_mat_vec_x, data.test_mat_vec_n);
        mat_vec(data.test_mat_vec_Pu, data.test_mat_vec_x, Px_cum, 1);          // upper
                                                                                  // traingular
                                                                                  // part
        mat_tpose_vec(data.test_mat_vec_Pu, data.test_mat_vec_x, Px_cum, 1, 1); // lower
                                                                                  // traingular
                                                                                  // part
                                                                                  // (without
                                                                                  // diagonal)
         assertTrue(
          "Linear algebra tests: error in matrix-vector operation, cumulative symmetric matrix-vector multiplication",
          vec_norm_inf_diff(Px_cum, data.test_mat_vec_Px_cum,
                            data.test_mat_vec_n) < TESTS_TOL);
        return 0;
      };
*/
};

