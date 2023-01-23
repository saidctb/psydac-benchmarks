from pyccel.decorators import types
@types("float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]",
       "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]",
       "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]",
       "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]",
       "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]",
       "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]",
       "int64[:]", "int64[:]", "int64[:]", "int64[:]", "int64[:]",
       "int64[:]", "int64[:]", "int64[:]", "int64[:]", "float64[:,:]",
       "float64[:,:]", "float64[:,:]",
       "int64", "int64", "int64", "int64", "int64", "int64", "int64",
       "int64", "int64", "int64", "int64", "int64", "int64", "int64",
       "int64", "int64", "int64", "int64", "int64", "int64", "int64",
       "int64", "int64", "int64", "int64", "int64", "int64",
       "float64[:,:,:,:,:,:]", "float64[:,:,:,:,:,:]", "float64[:,:,:,:,:,:]",
       "int64[:,:]", "int64[:,:,:]", "int64[:]", "int64[:]",
       "int64[:]", "int64[:]", "int64[:]", "int64[:]", "int64")
def assemble_matrix(global_test_basis_TE_0_1, global_test_basis_TE_0_2, global_test_basis_TE_0_3,
                             global_test_basis_TE_1_1, global_test_basis_TE_1_2, global_test_basis_TE_1_3,
                             global_test_basis_TE_2_1, global_test_basis_TE_2_2, global_test_basis_TE_2_3,
                             global_trial_basis_E_0_1, global_trial_basis_E_0_2, global_trial_basis_E_0_3,
                             global_trial_basis_E_1_1, global_trial_basis_E_1_2, global_trial_basis_E_1_3,
                             global_trial_basis_E_2_1, global_trial_basis_E_2_2, global_trial_basis_E_2_3,
                             global_span_TE_0_1, global_span_TE_0_2, global_span_TE_0_3, global_span_TE_1_1,
                             global_span_TE_1_2, global_span_TE_1_3, global_span_TE_2_1, global_span_TE_2_2,
                             global_span_TE_2_3, global_x1, global_x2, global_x3,
                             test_TE_0_p1, test_TE_0_p2, test_TE_0_p3,
                             test_TE_1_p1, test_TE_1_p2, test_TE_1_p3, test_TE_2_p1,
                             test_TE_2_p2, test_TE_2_p3, trial_E_0_p1, trial_E_0_p2,
                             trial_E_0_p3, trial_E_1_p1, trial_E_1_p2, trial_E_1_p3,
                             trial_E_2_p1, trial_E_2_p2, trial_E_2_p3,
                             n_element_1, n_element_2, n_element_3, k1, k2, k3, pad1, pad2, pad3,
                             g_mat_00, g_mat_11, g_mat_22,
                             coords_from_rank, rank_from_coords,
                             global_thread_starts_1, global_thread_starts_2, global_thread_starts_3,
                             global_thread_ends_1, global_thread_ends_2, global_thread_ends_3,
                             num_threads):

    from numpy import array, zeros, zeros_like, floor
    from pyccel.stdlib.internal.openmp import omp_get_thread_num
    
    #$ omp parallel default(private) &
    #$ shared(coords_from_rank,rank_from_coords,global_thread_starts_1,global_thread_starts_2,global_thread_starts_3,global_thread_ends_1,global_thread_ends_2,global_thread_ends_3,global_test_basis_TE_0_1,global_test_basis_TE_0_2,global_test_basis_TE_0_3,global_test_basis_TE_1_1,global_test_basis_TE_1_2,global_test_basis_TE_1_3,global_test_basis_TE_2_1,global_test_basis_TE_2_2,global_test_basis_TE_2_3,global_trial_basis_E_0_1,global_trial_basis_E_0_2,global_trial_basis_E_0_3,global_trial_basis_E_1_1,global_trial_basis_E_1_2,global_trial_basis_E_1_3,global_trial_basis_E_2_1,global_trial_basis_E_2_2,global_trial_basis_E_2_3,global_span_TE_0_1,global_span_TE_0_2,global_span_TE_0_3,global_span_TE_1_1,global_span_TE_1_2,global_span_TE_1_3,global_span_TE_2_1,global_span_TE_2_2,global_span_TE_2_3,global_x1,global_x2,global_x3,g_mat_00,g_mat_11,g_mat_22) &
    #$ firstprivate(test_TE_0_p1,test_TE_0_p2,test_TE_0_p3,test_TE_1_p1,test_TE_1_p2,test_TE_1_p3,test_TE_2_p1,test_TE_2_p2,test_TE_2_p3,trial_E_0_p1,trial_E_0_p2,trial_E_0_p3,trial_E_1_p1,trial_E_1_p2,trial_E_1_p3,trial_E_2_p1,trial_E_2_p2,trial_E_2_p3,n_element_1,n_element_2,n_element_3,k1,k2,k3,pad1,pad2,pad3,num_threads) 

    l_mat_00 = zeros((4, 5, 5, 7, 9, 9), dtype='float64')
    l_mat_11 = zeros((5, 4, 5, 9, 7, 9), dtype='float64')
    l_mat_22 = zeros((5, 5, 4, 9, 9, 7), dtype='float64')
    thread_id = omp_get_thread_num()

    thread_coords_1 = coords_from_rank[thread_id,0]
    thread_coords_2 = coords_from_rank[thread_id,1]
    thread_coords_3 = coords_from_rank[thread_id,2]

    global_thread_size_1 = 1 + global_thread_ends_1[thread_coords_1] - global_thread_starts_1[thread_coords_1]
    global_thread_size_2 = 1 + global_thread_ends_2[thread_coords_2] - global_thread_starts_2[thread_coords_2]
    global_thread_size_3 = 1 + global_thread_ends_3[thread_coords_3] - global_thread_starts_3[thread_coords_3]

    local_thread_starts_1 = array((global_thread_starts_1[thread_coords_1], global_thread_starts_1[thread_coords_1]+int((1/2)*global_thread_size_1)))
    local_thread_starts_2 = array((global_thread_starts_2[thread_coords_2], global_thread_starts_2[thread_coords_2]+int((1/2)*global_thread_size_2)))
    local_thread_starts_3 = array((global_thread_starts_3[thread_coords_3], global_thread_starts_3[thread_coords_3]+int((1/2)*global_thread_size_3)))

    local_thread_ends_1 = array((global_thread_starts_1[thread_coords_1]+int((1/2)*global_thread_size_1), global_thread_ends_1[thread_coords_1]))
    local_thread_ends_2 = array((global_thread_starts_2[thread_coords_2]+int((1/2)*global_thread_size_2), global_thread_ends_2[thread_coords_2]))
    local_thread_ends_3 = array((global_thread_starts_3[thread_coords_3]+int((1/2)*global_thread_size_3), global_thread_ends_3[thread_coords_3]))

    for local_i_element_1 in range(0, 2, 1):
        for local_i_element_2 in range(0, 2, 1):
            for local_i_element_3 in range(0, 2, 1):
                #$omp barrier
                for i_element_1 in range(local_thread_starts_1[local_i_element_1], local_thread_ends_1[local_i_element_1], 1):
                    span_TE_0_1 = global_span_TE_0_1[i_element_1]
                    span_TE_1_1 = global_span_TE_1_1[i_element_1]
                    span_TE_2_1 = global_span_TE_2_1[i_element_1]
                    for i_element_2 in range(local_thread_starts_2[local_i_element_2], local_thread_ends_2[local_i_element_2], 1):
                        span_TE_0_2 = global_span_TE_0_2[i_element_2]
                        span_TE_1_2 = global_span_TE_1_2[i_element_2]
                        span_TE_2_2 = global_span_TE_2_2[i_element_2]
                        for i_element_3 in range(local_thread_starts_3[local_i_element_3], local_thread_ends_3[local_i_element_3], 1):
                            span_TE_0_3 = global_span_TE_0_3[i_element_3]
                            span_TE_1_3 = global_span_TE_1_3[i_element_3]
                            span_TE_2_3 = global_span_TE_2_3[i_element_3]

                            for i_quad_1 in range(0, 5, 1):
                                for i_quad_2 in range(0, 5, 1):
                                    for i_quad_3 in range(0, 5, 1):

                                        for i_basis_1 in range(0, 4, 1):
                                            TE_0_1 = global_test_basis_TE_0_1[i_element_1,i_basis_1,0,i_quad_1]
                                            for i_basis_2 in range(0, 5, 1):
                                                TE_0_2 = global_test_basis_TE_0_2[i_element_2,i_basis_2,0,i_quad_2]
                                                for i_basis_3 in range(0, 5, 1):
                                                    TE_0_3 = global_test_basis_TE_0_3[i_element_3,i_basis_3,0,i_quad_3]
                                                    for j_basis_1 in range(0, 4, 1):
                                                        E_0_1 = global_trial_basis_E_0_1[i_element_1,j_basis_1,0,i_quad_1]
                                                        for j_basis_2 in range(0, 5, 1):
                                                            E_0_2  = global_trial_basis_E_0_2[i_element_2,j_basis_2,0,i_quad_2]
                                                            for j_basis_3 in range(0, 5, 1):
                                                                E_0_3 = global_trial_basis_E_0_3[i_element_3,j_basis_3,0,i_quad_3]
                                                                TE_0 = TE_0_1*TE_0_2*TE_0_3
                                                                E_0 = E_0_1*E_0_2*E_0_3
                                                                contribution = E_0*TE_0

                                                                l_mat_00[i_basis_1,i_basis_2,i_basis_3,3 - i_basis_1 + j_basis_1,4 - i_basis_2 + j_basis_2,4 - i_basis_3 + j_basis_3] += contribution

                                        for i_basis_1 in range(0, 5, 1):
                                            TE_1_1 = global_test_basis_TE_1_1[i_element_1,i_basis_1,0,i_quad_1]
                                            for i_basis_2 in range(0, 4, 1):
                                                TE_1_2 = global_test_basis_TE_1_2[i_element_2,i_basis_2,0,i_quad_2]
                                                for i_basis_3 in range(0, 5, 1):
                                                    TE_1_3 = global_test_basis_TE_1_3[i_element_3,i_basis_3,0,i_quad_3]
                                                    for j_basis_1 in range(0, 5, 1):
                                                        E_1_1 = global_trial_basis_E_1_1[i_element_1,j_basis_1,0,i_quad_1]
                                                        for j_basis_2 in range(0, 4, 1):
                                                            E_1_2 = global_trial_basis_E_1_2[i_element_2,j_basis_2,0,i_quad_2]
                                                            for j_basis_3 in range(0, 5, 1):
                                                                E_1_3 = global_trial_basis_E_1_3[i_element_3,j_basis_3,0,i_quad_3]

                                                                TE_1 = TE_1_1*TE_1_2*TE_1_3
                                                                E_1 = E_1_1*E_1_2*E_1_3
                                                                contribution = E_1*TE_1
              
                                                                l_mat_11[i_basis_1,i_basis_2,i_basis_3,4 - i_basis_1 + j_basis_1,3 - i_basis_2 + j_basis_2,4 - i_basis_3 + j_basis_3] += contribution

                                        for i_basis_1 in range(0, 5, 1):
                                            TE_2_1 = global_test_basis_TE_2_1[i_element_1,i_basis_1,0,i_quad_1]
                                            for i_basis_2 in range(0, 5, 1):
                                                TE_2_2 = global_test_basis_TE_2_2[i_element_2,i_basis_2,0,i_quad_2]
                                                for i_basis_3 in range(0, 4, 1):
                                                    TE_2_3 = global_test_basis_TE_2_3[i_element_3,i_basis_3,0,i_quad_3]
                                                    for j_basis_1 in range(0, 5, 1):
                                                        E_2_1 = global_trial_basis_E_2_1[i_element_1,j_basis_1,0,i_quad_1]
                                                        for j_basis_2 in range(0, 5, 1):
                                                            E_2_2 = global_trial_basis_E_2_2[i_element_2,j_basis_2,0,i_quad_2]
                                                            for j_basis_3 in range(0, 4, 1):
                                                                E_2_3 = global_trial_basis_E_2_3[i_element_3,j_basis_3,0,i_quad_3]
                                                                TE_2 = TE_2_1*TE_2_2*TE_2_3
                                                                E_2 = E_2_1*E_2_2*E_2_3
                                                                contribution = E_2*TE_2
                                                                        
                                                                l_mat_22[i_basis_1,i_basis_2,i_basis_3,4 - i_basis_1 + j_basis_1,4 - i_basis_2 + j_basis_2,3 - i_basis_3 + j_basis_3] += contribution
                                                            
                                            
                            g_mat_00[pad1 + span_TE_0_1 - test_TE_0_p1:1 + pad1 + span_TE_0_1,pad2 + span_TE_0_2 - test_TE_0_p2:1 + pad2 + span_TE_0_2,pad3 + span_TE_0_3 - test_TE_0_p3:1 + pad3 + span_TE_0_3,:,:,:] += l_mat_00[:,:,:,:,:,:]
                            g_mat_11[pad1 + span_TE_1_1 - test_TE_1_p1:1 + pad1 + span_TE_1_1,pad2 + span_TE_1_2 - test_TE_1_p2:1 + pad2 + span_TE_1_2,pad3 + span_TE_1_3 - test_TE_1_p3:1 + pad3 + span_TE_1_3,:,:,:] += l_mat_11[:,:,:,:,:,:]
                            g_mat_22[pad1 + span_TE_2_1 - test_TE_2_p1:1 + pad1 + span_TE_2_1,pad2 + span_TE_2_2 - test_TE_2_p2:1 + pad2 + span_TE_2_2,pad3 + span_TE_2_3 - test_TE_2_p3:1 + pad3 + span_TE_2_3,:,:,:] += l_mat_22[:,:,:,:,:,:]
                        
    #$ omp end parallel
    return
