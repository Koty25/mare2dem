/*******************************************************************************
* Copyright 2019-2020 Intel Corporation.
*
* This software and the related documents are Intel copyrighted  materials,  and
* your use of  them is  governed by the  express license  under which  they were
* provided to you (License).  Unless the License provides otherwise, you may not
* use, modify, copy, publish, distribute,  disclose or transmit this software or
* the related documents without Intel's prior written permission.
*
* This software and the related documents  are provided as  is,  with no express
* or implied  warranties,  other  than those  that are  expressly stated  in the
* License.
*******************************************************************************/

template <typename fp, typename intType>
void cleanup_arrays(std::vector<fp *> &fp_ptr_vec,
                    std::vector<intType *> &int_ptr_vec,
                    cl::sycl::context cxt)
{
    for (int i = 0; i < fp_ptr_vec.size(); i++)
        free(fp_ptr_vec[i], cxt);

    for (int i = 0; i < int_ptr_vec.size(); i++)
        free(int_ptr_vec[i], cxt);
}

// Creating the 3arrays CSR representation (ia, ja, values)
// of stencil-based matrix with size nx=ny=nz
template <typename fp, typename intType>
void generate_sparse_matrix(const intType nx,
                            std::vector<intType, mkl_allocator<intType, 64>> &ia,
                            std::vector<intType, mkl_allocator<intType, 64>> &ja,
                            std::vector<fp, mkl_allocator<fp, 64>> &a)
{
    intType nz = nx, ny = nx;
    intType nnz = 0;
    intType current_row;

    ia[0] = 0;

    for (intType iz = 0; iz < nz; iz++) {
        for (intType iy = 0; iy < ny; iy++) {
            for (intType ix = 0; ix < nx; ix++) {

                current_row = iz * nx * ny + iy * nx + ix;

                for (intType sz = -1; sz <= 1; sz++) {
                    if (iz + sz > -1 && iz + sz < nz) {
                        for (intType sy = -1; sy <= 1; sy++) {
                            if (iy + sy > -1 && iy + sy < ny) {
                                for (intType sx = -1; sx <= 1; sx++) {
                                    if (ix + sx > -1 && ix + sx < nx) {
                                        intType current_column =
                                                current_row + sz * nx * ny + sy * nx + sx;
                                        ja[nnz] = current_column;
                                        if (current_column == current_row) {
                                            a[nnz++] = set_fp_value(fp(26.0), fp(0.0));
                                        }
                                        else {
                                            a[nnz++] = set_fp_value(fp(-1.0), fp(0.0));
                                        }
                                    } // end
                                      // x
                                      // bounds
                                      // test
                                }     // end sx loop
                            }         // end y bounds test
                        }             // end sy loop
                    }                 // end z bounds test
                }                     // end sz loop
                ia[current_row + 1] = nnz;

            } // end ix loop
        }     // end iy loop
    }         // end iz loop
}

template <typename fp, typename intType>
void generate_sparse_matrix(const intType nx, intType *ia, intType *ja, fp *a)
{
    intType nz = nx, ny = nx;
    intType nnz = 0;
    intType current_row;

    ia[0] = 0;

    for (intType iz = 0; iz < nz; iz++) {
        for (intType iy = 0; iy < ny; iy++) {
            for (intType ix = 0; ix < nx; ix++) {

                current_row = iz * nx * ny + iy * nx + ix;

                for (intType sz = -1; sz <= 1; sz++) {
                    if (iz + sz > -1 && iz + sz < nz) {
                        for (intType sy = -1; sy <= 1; sy++) {
                            if (iy + sy > -1 && iy + sy < ny) {
                                for (intType sx = -1; sx <= 1; sx++) {
                                    if (ix + sx > -1 && ix + sx < nx) {
                                        intType current_column =
                                                current_row + sz * nx * ny + sy * nx + sx;
                                        ja[nnz] = current_column;
                                        if (current_column == current_row) {
                                            a[nnz++] = set_fp_value(fp(26.0), fp(0.0));
                                        }
                                        else {
                                            a[nnz++] = set_fp_value(fp(-1.0), fp(0.0));
                                        }
                                    } // end
                                      // x
                                      // bounds
                                      // test
                                }     // end sx loop
                            }         // end y bounds test
                        }             // end sy loop
                    }                 // end z bounds test
                }                     // end sz loop
                ia[current_row + 1] = nnz;

            } // end ix loop
        }     // end iy loop
    }         // end iz loop
}

// Creating the 3arrays CSR representation (ia, ja, values)
// of general random sparse matrix
// with density (0 < density <= 1.0)
// -0.5 <= value < 0.5
template <typename fp, typename intType>
void generate_random_sparse_matrix(const intType nrows,
                                   const intType ncols,
                                   const double density_val,
                                   std::vector<intType, mkl_allocator<intType, 64>> &ia,
                                   std::vector<intType, mkl_allocator<intType, 64>> &ja,
                                   std::vector<fp, mkl_allocator<fp, 64>> &a)
{
    intType nnz = 0;
    ia.push_back(0); // starting index of row0.

    for (intType i = 0; i < nrows; i++) {
        ia.push_back(nnz); // ending index of row_i.
        for (intType j = 0; j < ncols; j++) {
            if ((double)std::rand() / RAND_MAX < density_val) {
                a.push_back(rand_scalar<fp>());
                ja.push_back(j);
                nnz++;
            }
        }

        ia[i + 1] = nnz; // update ending index of row_i.
    }
}

template <typename fp, typename intType>
bool check_result(fp res, fp ref, intType index)
{
    bool check;
    fp bound = std::numeric_limits<fp>::epsilon();
    fp aerr  = std::abs(res - ref);
    fp rerr  = aerr / std::abs(ref);
    check    = (rerr <= bound) || (aerr <= bound);
    if (!check)
        std::cout << "relative error = " << rerr << " absolute error = " << aerr
                  << " limit = " << bound << " in index: " << index << std::endl;
    return check;
}

template <typename fp, typename intType>
bool check_result(fp res, fp ref, intType nFlops, intType index)
{
    bool check;
    fp bound = std::numeric_limits<fp>::epsilon() * nFlops;
    fp aerr  = std::abs(res - ref);
    fp rerr  = aerr / std::abs(ref);
    check    = (rerr <= bound) || (aerr <= bound);
    if (!check)
        std::cout << "relative error = " << rerr << " absolute error = " << aerr
                  << " limit = " << bound << " in index: " << index << std::endl;
    return check;
}
