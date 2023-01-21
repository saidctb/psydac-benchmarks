from pyccel.decorators import types,pure

@types("float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "int64[:]", "int64[:]", "int64[:]", "int64[:]", "int64[:]", "int64[:]", "float64[:,:]", "float64[:,:]", "float64[:,:]", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "float64[:,:,:]", "float64[:,:,:]", "float64[:,:,:]", "float64[:,:,:,:,:,:]")
def assemble_matrix_2(global_test_basis_v_1, global_test_basis_v_2, global_test_basis_v_3, global_trial_basis_u_1, global_trial_basis_u_2, global_trial_basis_u_3, global_test_basis_mapping_v_1, global_test_basis_mapping_v_2, global_test_basis_mapping_v_3, global_span_v_1, global_span_v_2, global_span_v_3, global_span_mapping_v_1, global_span_mapping_v_2, global_span_mapping_v_3, global_x1, global_x2, global_x3, test_v_p1, test_v_p2, test_v_p3, trial_u_p1, trial_u_p2, trial_u_p3, test_mapping_v_p1, test_mapping_v_p2, test_mapping_v_p3, n_element_1, n_element_2, n_element_3, k1, k2, k3, pad1, pad2, pad3, global_coeffs_x, global_coeffs_y, global_coeffs_z, g_mat_u_v_z6fv8xca):

    from numpy import array, zeros, zeros_like, floor
    from math import sqrt
    coeffs_x = zeros((1 + test_mapping_v_p1, 1 + test_mapping_v_p2, 1 + test_mapping_v_p3), dtype='float64')
    coeffs_y = zeros((1 + test_mapping_v_p1, 1 + test_mapping_v_p2, 1 + test_mapping_v_p3), dtype='float64')
    coeffs_z = zeros((1 + test_mapping_v_p1, 1 + test_mapping_v_p2, 1 + test_mapping_v_p3), dtype='float64')

    l_mat_u_v_z6fv8xca = zeros((3, 3, 3, 5, 5, 5), dtype='float64')
    for i_element_1 in range(0, n_element_1, 1):
        span_v_1 = global_span_v_1[i_element_1]
        span_mapping_v_1 = global_span_mapping_v_1[i_element_1]
        for i_element_2 in range(0, n_element_2, 1):
            span_v_2 = global_span_v_2[i_element_2]
            span_mapping_v_2 = global_span_mapping_v_2[i_element_2]
            for i_element_3 in range(0, n_element_3, 1):
                span_v_3 = global_span_v_3[i_element_3]
                span_mapping_v_3 = global_span_mapping_v_3[i_element_3]
                coeffs_x[:,:,:] = global_coeffs_x[2 + span_mapping_v_1 - test_mapping_v_p1:3 + span_mapping_v_1,2 + span_mapping_v_2 - test_mapping_v_p2:3 + span_mapping_v_2,2 + span_mapping_v_3 - test_mapping_v_p3:3 + span_mapping_v_3]
                coeffs_y[:,:,:] = global_coeffs_y[2 + span_mapping_v_1 - test_mapping_v_p1:3 + span_mapping_v_1,2 + span_mapping_v_2 - test_mapping_v_p2:3 + span_mapping_v_2,2 + span_mapping_v_3 - test_mapping_v_p3:3 + span_mapping_v_3]
                coeffs_z[:,:,:] = global_coeffs_z[2 + span_mapping_v_1 - test_mapping_v_p1:3 + span_mapping_v_1,2 + span_mapping_v_2 - test_mapping_v_p2:3 + span_mapping_v_2,2 + span_mapping_v_3 - test_mapping_v_p3:3 + span_mapping_v_3]
                for i_quad_1 in range(0, 3, 1):
                    for i_quad_2 in range(0, 3, 1):
                        for i_quad_3 in range(0, 3, 1):
                            x       = 0.
                            x_x3    = 0.
                            x_x3x3  = 0.
                            x_x2    = 0.
                            x_x2x3  = 0.
                            x_x2x2  = 0.
                            x_x1    = 0.
                            x_x1x3  = 0.
                            x_x1x2  = 0.
                            x_x1x1  = 0.
                            y       = 0.
                            y_x3    = 0.
                            y_x3x3  = 0.
                            y_x2    = 0.
                            y_x2x3  = 0.
                            y_x2x2  = 0.
                            y_x1    = 0.
                            y_x1x3  = 0.
                            y_x1x2  = 0.
                            y_x1x1  = 0.
                            z       = 0.
                            z_x3    = 0.
                            z_x3x3  = 0.
                            z_x2    = 0.
                            z_x2x3  = 0.
                            z_x2x2  = 0.
                            z_x1    = 0.
                            z_x1x3  = 0.
                            z_x1x2  = 0.
                            z_x1x1  = 0.
                            for i_basis_1 in range(0, 3, 1):
                                mapping_v_1 = global_test_basis_mapping_v_1[i_element_1,i_basis_1,0,i_quad_1]
                                mapping_v_1_x1 = global_test_basis_mapping_v_1[i_element_1,i_basis_1,1,i_quad_1]
                                mapping_v_1_x1x1 = global_test_basis_mapping_v_1[i_element_1,i_basis_1,2,i_quad_1]
                                for i_basis_2 in range(0, 3, 1):
                                    mapping_v_2 = global_test_basis_mapping_v_2[i_element_2,i_basis_2,0,i_quad_2]
                                    mapping_v_2_x2 = global_test_basis_mapping_v_2[i_element_2,i_basis_2,1,i_quad_2]
                                    mapping_v_2_x2x2 = global_test_basis_mapping_v_2[i_element_2,i_basis_2,2,i_quad_2]
                                    for i_basis_3 in range(0, 3, 1):
                                        mapping_v_3 = global_test_basis_mapping_v_3[i_element_3,i_basis_3,0,i_quad_3]
                                        mapping_v_3_x3 = global_test_basis_mapping_v_3[i_element_3,i_basis_3,1,i_quad_3]
                                        mapping_v_3_x3x3 = global_test_basis_mapping_v_3[i_element_3,i_basis_3,2,i_quad_3]
                                        coeff_x = coeffs_x[i_basis_1,i_basis_2,i_basis_3]
                                        coeff_y = coeffs_y[i_basis_1,i_basis_2,i_basis_3]
                                        coeff_z = coeffs_z[i_basis_1,i_basis_2,i_basis_3]
                                        mapping_v_x1x3 = mapping_v_1_x1*mapping_v_2*mapping_v_3_x3
                                        mapping_v_x2 = mapping_v_1*mapping_v_2_x2*mapping_v_3
                                        mapping_v_x1x1 = mapping_v_1_x1x1*mapping_v_2*mapping_v_3
                                        mapping_v_x1x2 = mapping_v_1_x1*mapping_v_2_x2*mapping_v_3
                                        mapping_v_x3x3 = mapping_v_1*mapping_v_2*mapping_v_3_x3x3
                                        mapping_v = mapping_v_1*mapping_v_2*mapping_v_3
                                        mapping_v_x3 = mapping_v_1*mapping_v_2*mapping_v_3_x3
                                        mapping_v_x2x3 = mapping_v_1*mapping_v_2_x2*mapping_v_3_x3
                                        mapping_v_x2x2 = mapping_v_1*mapping_v_2_x2x2*mapping_v_3
                                        mapping_v_x1 = mapping_v_1_x1*mapping_v_2*mapping_v_3
                                        x += mapping_v*coeff_x
                                        x_x3 += mapping_v_x3*coeff_x
                                        x_x3x3 += mapping_v_x3x3*coeff_x
                                        x_x2 += mapping_v_x2*coeff_x
                                        x_x2x3 += mapping_v_x2x3*coeff_x
                                        x_x2x2 += mapping_v_x2x2*coeff_x
                                        x_x1 += mapping_v_x1*coeff_x
                                        x_x1x3 += mapping_v_x1x3*coeff_x
                                        x_x1x2 += mapping_v_x1x2*coeff_x
                                        x_x1x1 += mapping_v_x1x1*coeff_x
                                        y += mapping_v*coeff_y
                                        y_x3 += mapping_v_x3*coeff_y
                                        y_x3x3 += mapping_v_x3x3*coeff_y
                                        y_x2 += mapping_v_x2*coeff_y
                                        y_x2x3 += mapping_v_x2x3*coeff_y
                                        y_x2x2 += mapping_v_x2x2*coeff_y
                                        y_x1 += mapping_v_x1*coeff_y
                                        y_x1x3 += mapping_v_x1x3*coeff_y
                                        y_x1x2 += mapping_v_x1x2*coeff_y
                                        y_x1x1 += mapping_v_x1x1*coeff_y
                                        z += mapping_v*coeff_z
                                        z_x3 += mapping_v_x3*coeff_z
                                        z_x3x3 += mapping_v_x3x3*coeff_z
                                        z_x2 += mapping_v_x2*coeff_z
                                        z_x2x3 += mapping_v_x2x3*coeff_z
                                        z_x2x2 += mapping_v_x2x2*coeff_z
                                        z_x1 += mapping_v_x1*coeff_z
                                        z_x1x3 += mapping_v_x1x3*coeff_z
                                        z_x1x2 += mapping_v_x1x2*coeff_z
                                        z_x1x1 += mapping_v_x1x1*coeff_z

                            temp0 = x_x1*y_x2
                            temp1 = x_x2*y_x3
                            temp2 = x_x3*y_x1
                            temp3 = x_x1*y_x3
                            temp4 = x_x2*y_x1
                            temp5 = x_x3*y_x2
                            temp6 = temp0*z_x3 + temp1*z_x1 + temp2*z_x2 - temp3*z_x2 - temp4*z_x3 - temp5*z_x1
                            temp7 = temp6**2
                            temp8 = temp0 - temp4
                            temp9 = x_x1x3*y_x2
                            temp10 = x_x2*y_x1x3
                            temp11 = -temp10
                            temp12 = x_x1*y_x2x3 - x_x2x3*y_x1
                            temp13 = temp11 + temp12 + temp9
                            temp14 = temp6**(-1)
                            temp15 = -temp1 + temp5
                            temp16 = x_x2*y_x3x3
                            temp17 = x_x2x3*y_x3
                            temp18 = x_x3*y_x2x3 + x_x3x3*y_x2
                            temp19 = -temp16 - temp17 + temp18
                            temp20 = -temp2 + temp3
                            temp21 = x_x1x3*y_x3
                            temp22 = x_x3*y_x1x3
                            temp23 = temp21 - temp22 + x_x1*y_x3x3 - x_x3x3*y_x1
                            temp24 = temp7**(-1)
                            temp25 = y_x3*z_x2x3
                            temp26 = y_x3x3*z_x2
                            temp27 = x_x2*z_x3x3 + x_x2x3*z_x3
                            temp28 = y_x2*z_x3x3 + y_x2x3*z_x3
                            temp29 = x_x3*z_x2x3
                            temp30 = x_x3x3*z_x2
                            temp31 = temp24*(-temp1*z_x1x3 + temp10*z_x3 + temp18*z_x1 + temp21*z_x2 - temp22*z_x2 + temp27*y_x1 - temp28*x_x1 + temp5*z_x1x3 - temp9*z_x3 + x_x1*(temp25 + temp26) - y_x1*(temp29 + temp30) - z_x1*(temp16 + temp17))

                            temp35 = temp14*temp8
                            temp36 = x_x1*z_x3 - x_x3*z_x1
                            temp37 = x_x1x2*z_x3
                            temp38 = x_x3*z_x1x2
                            temp39 = -temp38
                            temp40 = x_x1*z_x2x3 - x_x2x3*z_x1
                            temp41 = temp37 + temp39 + temp40
                            temp42 = x_x2*z_x3
                            temp43 = x_x3*z_x2
                            temp44 = temp42 - temp43
                            temp45 = x_x2x3*z_x2
                            temp46 = x_x3*z_x2x2
                            temp47 = x_x2*z_x2x3 + x_x2x2*z_x3
                            temp48 = -temp45 - temp46 + temp47
                            temp49 = x_x1*z_x2 - x_x2*z_x1
                            temp50 = x_x1x2*z_x2
                            temp51 = temp50 + x_x1*z_x2x2 - x_x2*z_x1x2 - x_x2x2*z_x1
                            temp52 = y_x2x3*z_x2
                            temp53 = y_x3*z_x2x2
                            temp54 = x_x2x3*y_x2 + x_x3*y_x2x2
                            temp55 = y_x2*z_x2x3 + y_x2x2*z_x3
                            temp56 = x_x2*y_x2x3
                            temp57 = x_x2x2*y_x3
                            temp58 = temp24*(-temp1*z_x1x2 - temp37*y_x2 + temp42*y_x1x2 - temp43*y_x1x2 + temp47*y_x1 + temp5*z_x1x2 + temp50*y_x3 + temp54*z_x1 - temp55*x_x1 + x_x1*(temp52 + temp53) - y_x1*(temp45 + temp46) - z_x1*(temp56 + temp57))
                            temp62 = temp14*temp36
                            temp63 = y_x2*z_x3
                            temp64 = y_x3*z_x2
                            temp65 = temp63 - temp64
                            temp66 = y_x1x3*z_x2
                            temp67 = -temp66
                            temp68 = y_x3*z_x1x2
                            temp69 = -temp68
                            temp70 = y_x1x2*z_x3
                            temp71 = y_x2*z_x1x3
                            temp72 = temp70 + temp71
                            temp73 = temp67 + temp69 + temp72
                            temp74 = y_x1*z_x3 - y_x3*z_x1
                            temp75 = y_x1*z_x1x3 + y_x1x1*z_x3 - y_x1x3*z_x1 - y_x3*z_x1x1
                            temp76 = -y_x1*z_x2 + y_x2*z_x1
                            temp77 = -y_x1*z_x1x2 - y_x1x1*z_x2 + y_x1x2*z_x1 + y_x2*z_x1x1
                            temp78 = x_x2*z_x1x3
                            temp79 = temp37 + temp78
                            temp80 = x_x3*y_x1x2
                            temp81 = temp80 + temp9
                            temp82 = x_x1x3*z_x2
                            temp83 = x_x1x2*y_x3
                            temp84 = temp24*(-temp1*z_x1x1 + temp42*y_x1x1 - temp43*y_x1x1 + temp5*z_x1x1 - temp63*x_x1x1 + temp64*x_x1x1 - temp72*x_x1 + temp79*y_x1 + temp81*z_x1 + x_x1*(temp66 + temp68) - y_x1*(temp38 + temp82) - z_x1*(temp10 + temp83))
                            temp88 = temp14*temp65
                            temp89 = x_x1*y_x2x2 + x_x1x2*y_x2 - x_x2*y_x1x2 - x_x2x2*y_x1
                            temp90 = temp54 - temp56 - temp57
                            temp91 = temp12 - temp80 + temp83
                            temp92 = temp14*temp20
                            temp93 = x_x1*z_x3x3 + x_x1x3*z_x3 - x_x3*z_x1x3 - x_x3x3*z_x1
                            temp94 = temp27 - temp29 - temp30
                            temp95 = temp40 - temp78 + temp82
                            temp96 = temp14*temp49
                            temp97 = x_x1*y_x1x2 + x_x1x1*y_x2 - x_x1x2*y_x1 - x_x2*y_x1x1
                            temp98 = temp11 + temp81 - temp83
                            temp99 = x_x1*y_x1x3 + x_x1x1*y_x3 - x_x1x3*y_x1 - x_x3*y_x1x1
                            temp100 = temp14*temp15
                            temp101 = x_x1*z_x1x3 + x_x1x1*z_x3 - x_x1x3*z_x1 - x_x3*z_x1x1
                            temp102 = temp39 + temp79 - temp82
                            temp103 = x_x1*z_x1x2 + x_x1x1*z_x2 - x_x1x2*z_x1 - x_x2*z_x1x1
                            temp104 = temp14*temp44
                            temp105 = -temp25 - temp26 + temp28
                            temp106 = y_x1*z_x3x3 + y_x1x3*z_x3 - y_x3*z_x1x3 - y_x3x3*z_x1
                            temp107 = y_x2x3*z_x1
                            temp108 = y_x1*z_x2x3
                            temp109 = temp107 - temp108 + temp67 + temp71
                            temp110 = temp14*temp76
                            temp111 = -temp52 - temp53 + temp55
                            temp112 = -temp107 + temp108 + temp69 + temp70
                            temp113 = -y_x1*z_x2x2 - y_x1x2*z_x2 + y_x2*z_x1x2 + y_x2x2*z_x1
                            temp114 = temp14*temp74
                            for i_basis_1 in range(0, 3, 1):
                                v_1 = global_test_basis_v_1[i_element_1,i_basis_1,0,i_quad_1]    
                                v_1_x1 = global_test_basis_v_1[i_element_1,i_basis_1,1,i_quad_1]
                                v_1_x1x1 = global_test_basis_v_1[i_element_1,i_basis_1,2,i_quad_1]
                                for i_basis_2 in range(0, 3, 1):
                                    v_2 = global_test_basis_v_2[i_element_2,i_basis_2,0,i_quad_2]
                                    v_2_x2 = global_test_basis_v_2[i_element_2,i_basis_2,1,i_quad_2]
                                    v_2_x2x2 = global_test_basis_v_2[i_element_2,i_basis_2,2,i_quad_2]
                                    for i_basis_3 in range(0, 3, 1):
                                        v_3 = global_test_basis_v_3[i_element_3,i_basis_3,0,i_quad_3]
                                        v_3_x3 = global_test_basis_v_3[i_element_3,i_basis_3,1,i_quad_3]
                                        v_3_x3x3 = global_test_basis_v_3[i_element_3,i_basis_3,2,i_quad_3]
                                        for j_basis_1 in range(0, 3, 1):
                                            u_1 = global_trial_basis_u_1[i_element_1,j_basis_1,0,i_quad_1]
                                            u_1_x1 = global_trial_basis_u_1[i_element_1,j_basis_1,1,i_quad_1]
                                            u_1_x1x1 = global_trial_basis_u_1[i_element_1,j_basis_1,2,i_quad_1]
                                            for j_basis_2 in range(0, 3, 1):
                                                u_2 = global_trial_basis_u_2[i_element_2,j_basis_2,0,i_quad_2]
                                                u_2_x2 = global_trial_basis_u_2[i_element_2,j_basis_2,1,i_quad_2]
                                                u_2_x2x2 = global_trial_basis_u_2[i_element_2,j_basis_2,2,i_quad_2]
                                                for j_basis_3 in range(0, 3, 1):
                                                    u_3 = global_trial_basis_u_3[i_element_3,j_basis_3,0,i_quad_3]
                                                    u_3_x3 = global_trial_basis_u_3[i_element_3,j_basis_3,1,i_quad_3]
                                                    u_3_x3x3 = global_trial_basis_u_3[i_element_3,j_basis_3,2,i_quad_3]

                                                    v = v_1*v_2*v_3
                                                    v_x3 = v_1*v_2*v_3_x3
                                                    v_x3x3 = v_1*v_2*v_3_x3x3
                                                    v_x2 = v_1*v_2_x2*v_3
                                                    v_x2x3 = v_1*v_2_x2*v_3_x3
                                                    v_x2x2 = v_1*v_2_x2x2*v_3
                                                    v_x1 = v_1_x1*v_2*v_3
                                                    v_x1x3 = v_1_x1*v_2*v_3_x3
                                                    v_x1x2 = v_1_x1*v_2_x2*v_3
                                                    v_x1x1 = v_1_x1x1*v_2*v_3
                                                    u = u_1*u_2*u_3
                                                    u_x3 = u_1*u_2*u_3_x3
                                                    u_x3x3 = u_1*u_2*u_3_x3x3
                                                    u_x2 = u_1*u_2_x2*u_3
                                                    u_x2x3 = u_1*u_2_x2*u_3_x3
                                                    u_x2x2 = u_1*u_2_x2x2*u_3
                                                    u_x1 = u_1_x1*u_2*u_3
                                                    u_x1x3 = u_1_x1*u_2*u_3_x3
                                                    u_x1x2 = u_1_x1*u_2_x2*u_3
                                                    u_x1x1 = u_1_x1x1*u_2*u_3
                                                    temp32  = temp8*u_x3
                                                    temp33  = temp15*u_x1
                                                    temp34  = temp20*u_x2
                                                    temp59  = temp36*u_x2
                                                    temp60  = temp44*u_x1
                                                    temp61  = temp49*u_x3
                                                    temp85  = temp65*u_x1
                                                    temp86  = temp74*u_x2
                                                    temp87  = temp76*u_x3
                                                    temp115 = temp8*v_x3
                                                    temp116 = temp15*v_x1
                                                    temp117 = temp20*v_x2
                                                    temp118 = temp36*v_x2
                                                    temp119 = temp44*v_x1
                                                    temp120 = temp49*v_x3
                                                    temp121 = temp65*v_x1
                                                    temp122 = temp74*v_x2
                                                    temp123 = temp76*v_x3
                                                    contribution_v_u_z6fv8xca = sqrt(temp7)*(-temp100*(temp115*temp84 - temp116*temp84 - temp117*temp84 - temp14*(temp15*v_x1x1 + temp98*v_x1) - temp14*(temp20*v_x1x2 + temp99*v_x2) + temp14*(temp8*v_x1x3 + temp97*v_x3)) - temp104*(temp118*temp84 - temp119*temp84 - temp120*temp84 + temp14*(temp101*v_x2 + temp36*v_x1x2) - temp14*(temp102*v_x1 + temp44*v_x1x1) - temp14*(temp103*v_x3 + temp49*v_x1x3)) - temp110*(temp121*temp31 - temp122*temp31 - temp123*temp31 + temp14*(temp105*v_x1 + temp65*v_x1x3) - temp14*(temp106*v_x2 + temp74*v_x2x3) - temp14*(temp109*v_x3 + temp76*v_x3x3)) - temp114*(temp121*temp58 - temp122*temp58 - temp123*temp58 + temp14*(temp111*v_x1 + temp65*v_x1x2) - temp14*(temp112*v_x2 + temp74*v_x2x2) - temp14*(temp113*v_x3 + temp76*v_x2x3)) + temp35*(temp115*temp31 - temp116*temp31 - temp117*temp31 + temp14*(temp13*v_x3 + temp8*v_x3x3) - temp14*(temp15*v_x1x3 + temp19*v_x1) - temp14*(temp20*v_x2x3 + temp23*v_x2)) + temp62*(temp118*temp58 - temp119*temp58 - temp120*temp58 + temp14*(temp36*v_x2x2 + temp41*v_x2) - temp14*(temp44*v_x1x2 + temp48*v_x1) - temp14*(temp49*v_x2x3 + temp51*v_x3)) + temp88*(temp121*temp84 - temp122*temp84 - temp123*temp84 + temp14*(temp65*v_x1x1 + temp73*v_x1) - temp14*(temp74*v_x1x2 + temp75*v_x2) - temp14*(temp76*v_x1x3 + temp77*v_x3)) - temp92*(temp115*temp58 - temp116*temp58 - temp117*temp58 - temp14*(temp15*v_x1x2 + temp90*v_x1) - temp14*(temp20*v_x2x2 + temp91*v_x2) + temp14*(temp8*v_x2x3 + temp89*v_x3)) - temp96*(temp118*temp31 - temp119*temp31 - temp120*temp31 + temp14*(temp36*v_x2x3 + temp93*v_x2) - temp14*(temp44*v_x1x3 + temp94*v_x1) - temp14*(temp49*v_x3x3 + temp95*v_x3)))*(-temp100*(-temp14*(temp15*u_x1x1 + temp98*u_x1) - temp14*(temp20*u_x1x2 + temp99*u_x2) + temp14*(temp8*u_x1x3 + temp97*u_x3) + temp32*temp84 - temp33*temp84 - temp34*temp84) - temp104*(temp14*(temp101*u_x2 + temp36*u_x1x2) - temp14*(temp102*u_x1 + temp44*u_x1x1) - temp14*(temp103*u_x3 + temp49*u_x1x3) + temp59*temp84 - temp60*temp84 - temp61*temp84) - temp110*(temp14*(temp105*u_x1 + temp65*u_x1x3) - temp14*(temp106*u_x2 + temp74*u_x2x3) - temp14*(temp109*u_x3 + temp76*u_x3x3) + temp31*temp85 - temp31*temp86 - temp31*temp87) - temp114*(temp14*(temp111*u_x1 + temp65*u_x1x2) - temp14*(temp112*u_x2 + temp74*u_x2x2) - temp14*(temp113*u_x3 + temp76*u_x2x3) + temp58*temp85 - temp58*temp86 - temp58*temp87) + temp35*(temp14*(temp13*u_x3 + temp8*u_x3x3) - temp14*(temp15*u_x1x3 + temp19*u_x1) - temp14*(temp20*u_x2x3 + temp23*u_x2) + temp31*temp32 - temp31*temp33 - temp31*temp34) + temp62*(temp14*(temp36*u_x2x2 + temp41*u_x2) - temp14*(temp44*u_x1x2 + temp48*u_x1) - temp14*(temp49*u_x2x3 + temp51*u_x3) + temp58*temp59 - temp58*temp60 - temp58*temp61) + temp88*(temp14*(temp65*u_x1x1 + temp73*u_x1) - temp14*(temp74*u_x1x2 + temp75*u_x2) - temp14*(temp76*u_x1x3 + temp77*u_x3) + temp84*temp85 - temp84*temp86 - temp84*temp87) - temp92*(-temp14*(temp15*u_x1x2 + temp90*u_x1) - temp14*(temp20*u_x2x2 + temp91*u_x2) + temp14*(temp8*u_x2x3 + temp89*u_x3) + temp32*temp58 - temp33*temp58 - temp34*temp58) - temp96*(temp14*(temp36*u_x2x3 + temp93*u_x2) - temp14*(temp44*u_x1x3 + temp94*u_x1) - temp14*(temp49*u_x3x3 + temp95*u_x3) + temp31*temp59 - temp31*temp60 - temp31*temp61))

                                                    l_mat_u_v_z6fv8xca[i_basis_1,i_basis_2,i_basis_3,2 - i_basis_1 + j_basis_1,2 - i_basis_2 + j_basis_2,2 - i_basis_3 + j_basis_3] += contribution_v_u_z6fv8xca
                                    

                g_mat_u_v_z6fv8xca[pad1 + span_v_1 - test_v_p1:1 + pad1 + span_v_1,pad2 + span_v_2 - test_v_p2:1 + pad2 + span_v_2,pad3 + span_v_3 - test_v_p3:1 + pad3 + span_v_3,:,:,:] += l_mat_u_v_z6fv8xca[:,:,:,:,:,:]
            
        
    
    return

@types("float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "int64[:]", "int64[:]", "int64[:]", "int64[:]", "int64[:]", "int64[:]", "float64[:,:]", "float64[:,:]", "float64[:,:]", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "float64[:,:,:]", "float64[:,:,:]", "float64[:,:,:]", "float64[:,:,:,:,:,:]")
def assemble_matrix_3(global_test_basis_v_1, global_test_basis_v_2, global_test_basis_v_3, global_trial_basis_u_1, global_trial_basis_u_2, global_trial_basis_u_3, global_test_basis_mapping_v_1, global_test_basis_mapping_v_2, global_test_basis_mapping_v_3, global_span_v_1, global_span_v_2, global_span_v_3, global_span_mapping_v_1, global_span_mapping_v_2, global_span_mapping_v_3, global_x1, global_x2, global_x3, test_v_p1, test_v_p2, test_v_p3, trial_u_p1, trial_u_p2, trial_u_p3, test_mapping_v_p1, test_mapping_v_p2, test_mapping_v_p3, n_element_1, n_element_2, n_element_3, k1, k2, k3, pad1, pad2, pad3, global_coeffs_x, global_coeffs_y, global_coeffs_z, g_mat_u_v_z6fv8xca):

    from numpy import array, zeros, zeros_like, floor
    from math import sqrt
    coeffs_x = zeros((1 + test_mapping_v_p1, 1 + test_mapping_v_p2, 1 + test_mapping_v_p3), dtype='float64')
    coeffs_y = zeros((1 + test_mapping_v_p1, 1 + test_mapping_v_p2, 1 + test_mapping_v_p3), dtype='float64')
    coeffs_z = zeros((1 + test_mapping_v_p1, 1 + test_mapping_v_p2, 1 + test_mapping_v_p3), dtype='float64')

    l_mat_u_v_z6fv8xca = zeros((4, 4, 4, 7, 7, 7), dtype='float64')
    for i_element_1 in range(0, n_element_1, 1):
        span_v_1 = global_span_v_1[i_element_1]
        span_mapping_v_1 = global_span_mapping_v_1[i_element_1]
        for i_element_2 in range(0, n_element_2, 1):
            span_v_2 = global_span_v_2[i_element_2]
            span_mapping_v_2 = global_span_mapping_v_2[i_element_2]
            for i_element_3 in range(0, n_element_3, 1):
                span_v_3 = global_span_v_3[i_element_3]
                span_mapping_v_3 = global_span_mapping_v_3[i_element_3]
                coeffs_x[:,:,:] = global_coeffs_x[3 + span_mapping_v_1 - test_mapping_v_p1:4 + span_mapping_v_1,3 + span_mapping_v_2 - test_mapping_v_p2:4 + span_mapping_v_2,3 + span_mapping_v_3 - test_mapping_v_p3:4 + span_mapping_v_3]
                coeffs_y[:,:,:] = global_coeffs_y[3 + span_mapping_v_1 - test_mapping_v_p1:4 + span_mapping_v_1,3 + span_mapping_v_2 - test_mapping_v_p2:4 + span_mapping_v_2,3 + span_mapping_v_3 - test_mapping_v_p3:4 + span_mapping_v_3]
                coeffs_z[:,:,:] = global_coeffs_z[3 + span_mapping_v_1 - test_mapping_v_p1:4 + span_mapping_v_1,3 + span_mapping_v_2 - test_mapping_v_p2:4 + span_mapping_v_2,3 + span_mapping_v_3 - test_mapping_v_p3:4 + span_mapping_v_3]
                for i_quad_1 in range(0, 4, 1):
                    for i_quad_2 in range(0, 4, 1):
                        for i_quad_3 in range(0, 4, 1):
                            x       = 0.
                            x_x3    = 0.
                            x_x3x3  = 0.
                            x_x2    = 0.
                            x_x2x3  = 0.
                            x_x2x2  = 0.
                            x_x1    = 0.
                            x_x1x3  = 0.
                            x_x1x2  = 0.
                            x_x1x1  = 0.
                            y       = 0.
                            y_x3    = 0.
                            y_x3x3  = 0.
                            y_x2    = 0.
                            y_x2x3  = 0.
                            y_x2x2  = 0.
                            y_x1    = 0.
                            y_x1x3  = 0.
                            y_x1x2  = 0.
                            y_x1x1  = 0.
                            z       = 0.
                            z_x3    = 0.
                            z_x3x3  = 0.
                            z_x2    = 0.
                            z_x2x3  = 0.
                            z_x2x2  = 0.
                            z_x1    = 0.
                            z_x1x3  = 0.
                            z_x1x2  = 0.
                            z_x1x1  = 0.
                            for i_basis_1 in range(0, 4, 1):
                                mapping_v_1 = global_test_basis_mapping_v_1[i_element_1,i_basis_1,0,i_quad_1]
                                mapping_v_1_x1 = global_test_basis_mapping_v_1[i_element_1,i_basis_1,1,i_quad_1]
                                mapping_v_1_x1x1 = global_test_basis_mapping_v_1[i_element_1,i_basis_1,2,i_quad_1]
                                for i_basis_2 in range(0, 4, 1):
                                    mapping_v_2 = global_test_basis_mapping_v_2[i_element_2,i_basis_2,0,i_quad_2]
                                    mapping_v_2_x2 = global_test_basis_mapping_v_2[i_element_2,i_basis_2,1,i_quad_2]
                                    mapping_v_2_x2x2 = global_test_basis_mapping_v_2[i_element_2,i_basis_2,2,i_quad_2]
                                    for i_basis_3 in range(0, 4, 1):
                                        mapping_v_3 = global_test_basis_mapping_v_3[i_element_3,i_basis_3,0,i_quad_3]
                                        mapping_v_3_x3 = global_test_basis_mapping_v_3[i_element_3,i_basis_3,1,i_quad_3]
                                        mapping_v_3_x3x3 = global_test_basis_mapping_v_3[i_element_3,i_basis_3,2,i_quad_3]
                                        coeff_x = coeffs_x[i_basis_1,i_basis_2,i_basis_3]
                                        coeff_y = coeffs_y[i_basis_1,i_basis_2,i_basis_3]
                                        coeff_z = coeffs_z[i_basis_1,i_basis_2,i_basis_3]
                                        mapping_v_x1x3 = mapping_v_1_x1*mapping_v_2*mapping_v_3_x3
                                        mapping_v_x2 = mapping_v_1*mapping_v_2_x2*mapping_v_3
                                        mapping_v_x1x1 = mapping_v_1_x1x1*mapping_v_2*mapping_v_3
                                        mapping_v_x1x2 = mapping_v_1_x1*mapping_v_2_x2*mapping_v_3
                                        mapping_v_x3x3 = mapping_v_1*mapping_v_2*mapping_v_3_x3x3
                                        mapping_v = mapping_v_1*mapping_v_2*mapping_v_3
                                        mapping_v_x3 = mapping_v_1*mapping_v_2*mapping_v_3_x3
                                        mapping_v_x2x3 = mapping_v_1*mapping_v_2_x2*mapping_v_3_x3
                                        mapping_v_x2x2 = mapping_v_1*mapping_v_2_x2x2*mapping_v_3
                                        mapping_v_x1 = mapping_v_1_x1*mapping_v_2*mapping_v_3
                                        x += mapping_v*coeff_x
                                        x_x3 += mapping_v_x3*coeff_x
                                        x_x3x3 += mapping_v_x3x3*coeff_x
                                        x_x2 += mapping_v_x2*coeff_x
                                        x_x2x3 += mapping_v_x2x3*coeff_x
                                        x_x2x2 += mapping_v_x2x2*coeff_x
                                        x_x1 += mapping_v_x1*coeff_x
                                        x_x1x3 += mapping_v_x1x3*coeff_x
                                        x_x1x2 += mapping_v_x1x2*coeff_x
                                        x_x1x1 += mapping_v_x1x1*coeff_x
                                        y += mapping_v*coeff_y
                                        y_x3 += mapping_v_x3*coeff_y
                                        y_x3x3 += mapping_v_x3x3*coeff_y
                                        y_x2 += mapping_v_x2*coeff_y
                                        y_x2x3 += mapping_v_x2x3*coeff_y
                                        y_x2x2 += mapping_v_x2x2*coeff_y
                                        y_x1 += mapping_v_x1*coeff_y
                                        y_x1x3 += mapping_v_x1x3*coeff_y
                                        y_x1x2 += mapping_v_x1x2*coeff_y
                                        y_x1x1 += mapping_v_x1x1*coeff_y
                                        z += mapping_v*coeff_z
                                        z_x3 += mapping_v_x3*coeff_z
                                        z_x3x3 += mapping_v_x3x3*coeff_z
                                        z_x2 += mapping_v_x2*coeff_z
                                        z_x2x3 += mapping_v_x2x3*coeff_z
                                        z_x2x2 += mapping_v_x2x2*coeff_z
                                        z_x1 += mapping_v_x1*coeff_z
                                        z_x1x3 += mapping_v_x1x3*coeff_z
                                        z_x1x2 += mapping_v_x1x2*coeff_z
                                        z_x1x1 += mapping_v_x1x1*coeff_z
                            temp0 = x_x1*y_x2
                            temp1 = x_x2*y_x3
                            temp2 = x_x3*y_x1
                            temp3 = x_x1*y_x3
                            temp4 = x_x2*y_x1
                            temp5 = x_x3*y_x2
                            temp6 = temp0*z_x3 + temp1*z_x1 + temp2*z_x2 - temp3*z_x2 - temp4*z_x3 - temp5*z_x1
                            temp7 = temp6**2
                            temp8 = temp0 - temp4
                            temp9 = x_x1x3*y_x2
                            temp10 = x_x2*y_x1x3
                            temp11 = -temp10
                            temp12 = x_x1*y_x2x3 - x_x2x3*y_x1
                            temp13 = temp11 + temp12 + temp9
                            temp14 = temp6**(-1)
                            temp15 = -temp1 + temp5
                            temp16 = x_x2*y_x3x3
                            temp17 = x_x2x3*y_x3
                            temp18 = x_x3*y_x2x3 + x_x3x3*y_x2
                            temp19 = -temp16 - temp17 + temp18
                            temp20 = -temp2 + temp3
                            temp21 = x_x1x3*y_x3
                            temp22 = x_x3*y_x1x3
                            temp23 = temp21 - temp22 + x_x1*y_x3x3 - x_x3x3*y_x1
                            temp24 = temp7**(-1)
                            temp25 = y_x3*z_x2x3
                            temp26 = y_x3x3*z_x2
                            temp27 = x_x2*z_x3x3 + x_x2x3*z_x3
                            temp28 = y_x2*z_x3x3 + y_x2x3*z_x3
                            temp29 = x_x3*z_x2x3
                            temp30 = x_x3x3*z_x2
                            temp31 = temp24*(-temp1*z_x1x3 + temp10*z_x3 + temp18*z_x1 + temp21*z_x2 - temp22*z_x2 + temp27*y_x1 - temp28*x_x1 + temp5*z_x1x3 - temp9*z_x3 + x_x1*(temp25 + temp26) - y_x1*(temp29 + temp30) - z_x1*(temp16 + temp17))

                            temp35 = temp14*temp8
                            temp36 = x_x1*z_x3 - x_x3*z_x1
                            temp37 = x_x1x2*z_x3
                            temp38 = x_x3*z_x1x2
                            temp39 = -temp38
                            temp40 = x_x1*z_x2x3 - x_x2x3*z_x1
                            temp41 = temp37 + temp39 + temp40
                            temp42 = x_x2*z_x3
                            temp43 = x_x3*z_x2
                            temp44 = temp42 - temp43
                            temp45 = x_x2x3*z_x2
                            temp46 = x_x3*z_x2x2
                            temp47 = x_x2*z_x2x3 + x_x2x2*z_x3
                            temp48 = -temp45 - temp46 + temp47
                            temp49 = x_x1*z_x2 - x_x2*z_x1
                            temp50 = x_x1x2*z_x2
                            temp51 = temp50 + x_x1*z_x2x2 - x_x2*z_x1x2 - x_x2x2*z_x1
                            temp52 = y_x2x3*z_x2
                            temp53 = y_x3*z_x2x2
                            temp54 = x_x2x3*y_x2 + x_x3*y_x2x2
                            temp55 = y_x2*z_x2x3 + y_x2x2*z_x3
                            temp56 = x_x2*y_x2x3
                            temp57 = x_x2x2*y_x3
                            temp58 = temp24*(-temp1*z_x1x2 - temp37*y_x2 + temp42*y_x1x2 - temp43*y_x1x2 + temp47*y_x1 + temp5*z_x1x2 + temp50*y_x3 + temp54*z_x1 - temp55*x_x1 + x_x1*(temp52 + temp53) - y_x1*(temp45 + temp46) - z_x1*(temp56 + temp57))
                            temp62 = temp14*temp36
                            temp63 = y_x2*z_x3
                            temp64 = y_x3*z_x2
                            temp65 = temp63 - temp64
                            temp66 = y_x1x3*z_x2
                            temp67 = -temp66
                            temp68 = y_x3*z_x1x2
                            temp69 = -temp68
                            temp70 = y_x1x2*z_x3
                            temp71 = y_x2*z_x1x3
                            temp72 = temp70 + temp71
                            temp73 = temp67 + temp69 + temp72
                            temp74 = y_x1*z_x3 - y_x3*z_x1
                            temp75 = y_x1*z_x1x3 + y_x1x1*z_x3 - y_x1x3*z_x1 - y_x3*z_x1x1
                            temp76 = -y_x1*z_x2 + y_x2*z_x1
                            temp77 = -y_x1*z_x1x2 - y_x1x1*z_x2 + y_x1x2*z_x1 + y_x2*z_x1x1
                            temp78 = x_x2*z_x1x3
                            temp79 = temp37 + temp78
                            temp80 = x_x3*y_x1x2
                            temp81 = temp80 + temp9
                            temp82 = x_x1x3*z_x2
                            temp83 = x_x1x2*y_x3
                            temp84 = temp24*(-temp1*z_x1x1 + temp42*y_x1x1 - temp43*y_x1x1 + temp5*z_x1x1 - temp63*x_x1x1 + temp64*x_x1x1 - temp72*x_x1 + temp79*y_x1 + temp81*z_x1 + x_x1*(temp66 + temp68) - y_x1*(temp38 + temp82) - z_x1*(temp10 + temp83))
                            temp88 = temp14*temp65
                            temp89 = x_x1*y_x2x2 + x_x1x2*y_x2 - x_x2*y_x1x2 - x_x2x2*y_x1
                            temp90 = temp54 - temp56 - temp57
                            temp91 = temp12 - temp80 + temp83
                            temp92 = temp14*temp20
                            temp93 = x_x1*z_x3x3 + x_x1x3*z_x3 - x_x3*z_x1x3 - x_x3x3*z_x1
                            temp94 = temp27 - temp29 - temp30
                            temp95 = temp40 - temp78 + temp82
                            temp96 = temp14*temp49
                            temp97 = x_x1*y_x1x2 + x_x1x1*y_x2 - x_x1x2*y_x1 - x_x2*y_x1x1
                            temp98 = temp11 + temp81 - temp83
                            temp99 = x_x1*y_x1x3 + x_x1x1*y_x3 - x_x1x3*y_x1 - x_x3*y_x1x1
                            temp100 = temp14*temp15
                            temp101 = x_x1*z_x1x3 + x_x1x1*z_x3 - x_x1x3*z_x1 - x_x3*z_x1x1
                            temp102 = temp39 + temp79 - temp82
                            temp103 = x_x1*z_x1x2 + x_x1x1*z_x2 - x_x1x2*z_x1 - x_x2*z_x1x1
                            temp104 = temp14*temp44
                            temp105 = -temp25 - temp26 + temp28
                            temp106 = y_x1*z_x3x3 + y_x1x3*z_x3 - y_x3*z_x1x3 - y_x3x3*z_x1
                            temp107 = y_x2x3*z_x1
                            temp108 = y_x1*z_x2x3
                            temp109 = temp107 - temp108 + temp67 + temp71
                            temp110 = temp14*temp76
                            temp111 = -temp52 - temp53 + temp55
                            temp112 = -temp107 + temp108 + temp69 + temp70
                            temp113 = -y_x1*z_x2x2 - y_x1x2*z_x2 + y_x2*z_x1x2 + y_x2x2*z_x1
                            temp114 = temp14*temp74
                            for i_basis_1 in range(0, 4, 1):
                                v_1 = global_test_basis_v_1[i_element_1,i_basis_1,0,i_quad_1]    
                                v_1_x1 = global_test_basis_v_1[i_element_1,i_basis_1,1,i_quad_1]
                                v_1_x1x1 = global_test_basis_v_1[i_element_1,i_basis_1,2,i_quad_1]
                                for i_basis_2 in range(0, 4, 1):
                                    v_2 = global_test_basis_v_2[i_element_2,i_basis_2,0,i_quad_2]
                                    v_2_x2 = global_test_basis_v_2[i_element_2,i_basis_2,1,i_quad_2]
                                    v_2_x2x2 = global_test_basis_v_2[i_element_2,i_basis_2,2,i_quad_2]
                                    for i_basis_3 in range(0, 4, 1):
                                        v_3 = global_test_basis_v_3[i_element_3,i_basis_3,0,i_quad_3]
                                        v_3_x3 = global_test_basis_v_3[i_element_3,i_basis_3,1,i_quad_3]
                                        v_3_x3x3 = global_test_basis_v_3[i_element_3,i_basis_3,2,i_quad_3]
                                        for j_basis_1 in range(0, 4, 1):
                                            u_1 = global_trial_basis_u_1[i_element_1,j_basis_1,0,i_quad_1]
                                            u_1_x1 = global_trial_basis_u_1[i_element_1,j_basis_1,1,i_quad_1]
                                            u_1_x1x1 = global_trial_basis_u_1[i_element_1,j_basis_1,2,i_quad_1]
                                            for j_basis_2 in range(0, 4, 1):
                                                u_2 = global_trial_basis_u_2[i_element_2,j_basis_2,0,i_quad_2]
                                                u_2_x2 = global_trial_basis_u_2[i_element_2,j_basis_2,1,i_quad_2]
                                                u_2_x2x2 = global_trial_basis_u_2[i_element_2,j_basis_2,2,i_quad_2]
                                                for j_basis_3 in range(0, 4, 1):
                                                    u_3 = global_trial_basis_u_3[i_element_3,j_basis_3,0,i_quad_3]
                                                    u_3_x3 = global_trial_basis_u_3[i_element_3,j_basis_3,1,i_quad_3]
                                                    u_3_x3x3 = global_trial_basis_u_3[i_element_3,j_basis_3,2,i_quad_3]

                                                    v = v_1*v_2*v_3
                                                    v_x3 = v_1*v_2*v_3_x3
                                                    v_x3x3 = v_1*v_2*v_3_x3x3
                                                    v_x2 = v_1*v_2_x2*v_3
                                                    v_x2x3 = v_1*v_2_x2*v_3_x3
                                                    v_x2x2 = v_1*v_2_x2x2*v_3
                                                    v_x1 = v_1_x1*v_2*v_3
                                                    v_x1x3 = v_1_x1*v_2*v_3_x3
                                                    v_x1x2 = v_1_x1*v_2_x2*v_3
                                                    v_x1x1 = v_1_x1x1*v_2*v_3
                                                    u = u_1*u_2*u_3
                                                    u_x3 = u_1*u_2*u_3_x3
                                                    u_x3x3 = u_1*u_2*u_3_x3x3
                                                    u_x2 = u_1*u_2_x2*u_3
                                                    u_x2x3 = u_1*u_2_x2*u_3_x3
                                                    u_x2x2 = u_1*u_2_x2x2*u_3
                                                    u_x1 = u_1_x1*u_2*u_3
                                                    u_x1x3 = u_1_x1*u_2*u_3_x3
                                                    u_x1x2 = u_1_x1*u_2_x2*u_3
                                                    u_x1x1 = u_1_x1x1*u_2*u_3
                                                    temp32  = temp8*u_x3
                                                    temp33  = temp15*u_x1
                                                    temp34  = temp20*u_x2
                                                    temp59  = temp36*u_x2
                                                    temp60  = temp44*u_x1
                                                    temp61  = temp49*u_x3
                                                    temp85  = temp65*u_x1
                                                    temp86  = temp74*u_x2
                                                    temp87  = temp76*u_x3
                                                    temp115 = temp8*v_x3
                                                    temp116 = temp15*v_x1
                                                    temp117 = temp20*v_x2
                                                    temp118 = temp36*v_x2
                                                    temp119 = temp44*v_x1
                                                    temp120 = temp49*v_x3
                                                    temp121 = temp65*v_x1
                                                    temp122 = temp74*v_x2
                                                    temp123 = temp76*v_x3
                                                    contribution_v_u_z6fv8xca = sqrt(temp7)*(-temp100*(temp115*temp84 - temp116*temp84 - temp117*temp84 - temp14*(temp15*v_x1x1 + temp98*v_x1) - temp14*(temp20*v_x1x2 + temp99*v_x2) + temp14*(temp8*v_x1x3 + temp97*v_x3)) - temp104*(temp118*temp84 - temp119*temp84 - temp120*temp84 + temp14*(temp101*v_x2 + temp36*v_x1x2) - temp14*(temp102*v_x1 + temp44*v_x1x1) - temp14*(temp103*v_x3 + temp49*v_x1x3)) - temp110*(temp121*temp31 - temp122*temp31 - temp123*temp31 + temp14*(temp105*v_x1 + temp65*v_x1x3) - temp14*(temp106*v_x2 + temp74*v_x2x3) - temp14*(temp109*v_x3 + temp76*v_x3x3)) - temp114*(temp121*temp58 - temp122*temp58 - temp123*temp58 + temp14*(temp111*v_x1 + temp65*v_x1x2) - temp14*(temp112*v_x2 + temp74*v_x2x2) - temp14*(temp113*v_x3 + temp76*v_x2x3)) + temp35*(temp115*temp31 - temp116*temp31 - temp117*temp31 + temp14*(temp13*v_x3 + temp8*v_x3x3) - temp14*(temp15*v_x1x3 + temp19*v_x1) - temp14*(temp20*v_x2x3 + temp23*v_x2)) + temp62*(temp118*temp58 - temp119*temp58 - temp120*temp58 + temp14*(temp36*v_x2x2 + temp41*v_x2) - temp14*(temp44*v_x1x2 + temp48*v_x1) - temp14*(temp49*v_x2x3 + temp51*v_x3)) + temp88*(temp121*temp84 - temp122*temp84 - temp123*temp84 + temp14*(temp65*v_x1x1 + temp73*v_x1) - temp14*(temp74*v_x1x2 + temp75*v_x2) - temp14*(temp76*v_x1x3 + temp77*v_x3)) - temp92*(temp115*temp58 - temp116*temp58 - temp117*temp58 - temp14*(temp15*v_x1x2 + temp90*v_x1) - temp14*(temp20*v_x2x2 + temp91*v_x2) + temp14*(temp8*v_x2x3 + temp89*v_x3)) - temp96*(temp118*temp31 - temp119*temp31 - temp120*temp31 + temp14*(temp36*v_x2x3 + temp93*v_x2) - temp14*(temp44*v_x1x3 + temp94*v_x1) - temp14*(temp49*v_x3x3 + temp95*v_x3)))*(-temp100*(-temp14*(temp15*u_x1x1 + temp98*u_x1) - temp14*(temp20*u_x1x2 + temp99*u_x2) + temp14*(temp8*u_x1x3 + temp97*u_x3) + temp32*temp84 - temp33*temp84 - temp34*temp84) - temp104*(temp14*(temp101*u_x2 + temp36*u_x1x2) - temp14*(temp102*u_x1 + temp44*u_x1x1) - temp14*(temp103*u_x3 + temp49*u_x1x3) + temp59*temp84 - temp60*temp84 - temp61*temp84) - temp110*(temp14*(temp105*u_x1 + temp65*u_x1x3) - temp14*(temp106*u_x2 + temp74*u_x2x3) - temp14*(temp109*u_x3 + temp76*u_x3x3) + temp31*temp85 - temp31*temp86 - temp31*temp87) - temp114*(temp14*(temp111*u_x1 + temp65*u_x1x2) - temp14*(temp112*u_x2 + temp74*u_x2x2) - temp14*(temp113*u_x3 + temp76*u_x2x3) + temp58*temp85 - temp58*temp86 - temp58*temp87) + temp35*(temp14*(temp13*u_x3 + temp8*u_x3x3) - temp14*(temp15*u_x1x3 + temp19*u_x1) - temp14*(temp20*u_x2x3 + temp23*u_x2) + temp31*temp32 - temp31*temp33 - temp31*temp34) + temp62*(temp14*(temp36*u_x2x2 + temp41*u_x2) - temp14*(temp44*u_x1x2 + temp48*u_x1) - temp14*(temp49*u_x2x3 + temp51*u_x3) + temp58*temp59 - temp58*temp60 - temp58*temp61) + temp88*(temp14*(temp65*u_x1x1 + temp73*u_x1) - temp14*(temp74*u_x1x2 + temp75*u_x2) - temp14*(temp76*u_x1x3 + temp77*u_x3) + temp84*temp85 - temp84*temp86 - temp84*temp87) - temp92*(-temp14*(temp15*u_x1x2 + temp90*u_x1) - temp14*(temp20*u_x2x2 + temp91*u_x2) + temp14*(temp8*u_x2x3 + temp89*u_x3) + temp32*temp58 - temp33*temp58 - temp34*temp58) - temp96*(temp14*(temp36*u_x2x3 + temp93*u_x2) - temp14*(temp44*u_x1x3 + temp94*u_x1) - temp14*(temp49*u_x3x3 + temp95*u_x3) + temp31*temp59 - temp31*temp60 - temp31*temp61))

                                                    l_mat_u_v_z6fv8xca[i_basis_1,i_basis_2,i_basis_3,3 - i_basis_1 + j_basis_1,3 - i_basis_2 + j_basis_2,3 - i_basis_3 + j_basis_3] += contribution_v_u_z6fv8xca
                                    

                g_mat_u_v_z6fv8xca[pad1 + span_v_1 - test_v_p1:1 + pad1 + span_v_1,pad2 + span_v_2 - test_v_p2:1 + pad2 + span_v_2,pad3 + span_v_3 - test_v_p3:1 + pad3 + span_v_3,:,:,:] += l_mat_u_v_z6fv8xca[:,:,:,:,:,:]
    return

@types("float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "float64[:,:,:,:]", "int64[:]", "int64[:]", "int64[:]", "int64[:]", "int64[:]", "int64[:]", "float64[:,:]", "float64[:,:]", "float64[:,:]", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "float64[:,:,:]", "float64[:,:,:]", "float64[:,:,:]", "float64[:,:,:,:,:,:]")
def assemble_matrix_4(global_test_basis_v_1, global_test_basis_v_2, global_test_basis_v_3, global_trial_basis_u_1, global_trial_basis_u_2, global_trial_basis_u_3, global_test_basis_mapping_v_1, global_test_basis_mapping_v_2, global_test_basis_mapping_v_3, global_span_v_1, global_span_v_2, global_span_v_3, global_span_mapping_v_1, global_span_mapping_v_2, global_span_mapping_v_3, global_x1, global_x2, global_x3, test_v_p1, test_v_p2, test_v_p3, trial_u_p1, trial_u_p2, trial_u_p3, test_mapping_v_p1, test_mapping_v_p2, test_mapping_v_p3, n_element_1, n_element_2, n_element_3, k1, k2, k3, pad1, pad2, pad3, global_coeffs_x, global_coeffs_y, global_coeffs_z, g_mat_u_v_z6fv8xca):

    from numpy import array, zeros, zeros_like, floor
    from math import sqrt
    coeffs_x = zeros((5,5,5), dtype='float64')
    coeffs_y = zeros((5,5,5), dtype='float64')
    coeffs_z = zeros((5,5,5), dtype='float64')

    l_mat_u_v = zeros((5, 5, 5, 9, 9, 9), dtype='float64')
    for i_element_1 in range(0, n_element_1, 1):
        span_v_1 = global_span_v_1[i_element_1]
        span_mapping_v_1 = global_span_mapping_v_1[i_element_1]
        for i_element_2 in range(0, n_element_2, 1):
            span_v_2 = global_span_v_2[i_element_2]
            span_mapping_v_2 = global_span_mapping_v_2[i_element_2]
            for i_element_3 in range(0, n_element_3, 1):
                span_v_3 = global_span_v_3[i_element_3]
                span_mapping_v_3 = global_span_mapping_v_3[i_element_3]
                coeffs_x[:,:,:] = global_coeffs_x[span_mapping_v_1:5 + span_mapping_v_1,span_mapping_v_2:5 + span_mapping_v_2,span_mapping_v_3:5 + span_mapping_v_3]
                coeffs_y[:,:,:] = global_coeffs_y[span_mapping_v_1:5 + span_mapping_v_1,span_mapping_v_2:5 + span_mapping_v_2,span_mapping_v_3:5 + span_mapping_v_3]
                coeffs_z[:,:,:] = global_coeffs_z[span_mapping_v_1:5 + span_mapping_v_1,span_mapping_v_2:5 + span_mapping_v_2,span_mapping_v_3:5 + span_mapping_v_3]
                for i_quad_1 in range(0, 5, 1):
                    for i_quad_2 in range(0, 5, 1):
                        for i_quad_3 in range(0, 5, 1):
                            x       = 0.
                            x_x3    = 0.
                            x_x3x3  = 0.
                            x_x2    = 0.
                            x_x2x3  = 0.
                            x_x2x2  = 0.
                            x_x1    = 0.
                            x_x1x3  = 0.
                            x_x1x2  = 0.
                            x_x1x1  = 0.
                            y       = 0.
                            y_x3    = 0.
                            y_x3x3  = 0.
                            y_x2    = 0.
                            y_x2x3  = 0.
                            y_x2x2  = 0.
                            y_x1    = 0.
                            y_x1x3  = 0.
                            y_x1x2  = 0.
                            y_x1x1  = 0.
                            z       = 0.
                            z_x3    = 0.
                            z_x3x3  = 0.
                            z_x2    = 0.
                            z_x2x3  = 0.
                            z_x2x2  = 0.
                            z_x1    = 0.
                            z_x1x3  = 0.
                            z_x1x2  = 0.
                            z_x1x1  = 0.
                            for i_basis_1 in range(0, 5, 1):
                                mapping_v_1 = global_test_basis_mapping_v_1[i_element_1,i_basis_1,0,i_quad_1]
                                mapping_v_1_x1 = global_test_basis_mapping_v_1[i_element_1,i_basis_1,1,i_quad_1]
                                mapping_v_1_x1x1 = global_test_basis_mapping_v_1[i_element_1,i_basis_1,2,i_quad_1]
                                for i_basis_2 in range(0, 5, 1):
                                    mapping_v_2 = global_test_basis_mapping_v_2[i_element_2,i_basis_2,0,i_quad_2]
                                    mapping_v_2_x2 = global_test_basis_mapping_v_2[i_element_2,i_basis_2,1,i_quad_2]
                                    mapping_v_2_x2x2 = global_test_basis_mapping_v_2[i_element_2,i_basis_2,2,i_quad_2]
                                    for i_basis_3 in range(0, 5, 1):
                                        mapping_v_3 = global_test_basis_mapping_v_3[i_element_3,i_basis_3,0,i_quad_3]
                                        mapping_v_3_x3 = global_test_basis_mapping_v_3[i_element_3,i_basis_3,1,i_quad_3]
                                        mapping_v_3_x3x3 = global_test_basis_mapping_v_3[i_element_3,i_basis_3,2,i_quad_3]
                                        coeff_x = coeffs_x[i_basis_1,i_basis_2,i_basis_3]
                                        coeff_y = coeffs_y[i_basis_1,i_basis_2,i_basis_3]
                                        coeff_z = coeffs_z[i_basis_1,i_basis_2,i_basis_3]
                                        mapping_v_x1x3 = mapping_v_1_x1*mapping_v_2*mapping_v_3_x3
                                        mapping_v_x2 = mapping_v_1*mapping_v_2_x2*mapping_v_3
                                        mapping_v_x1x1 = mapping_v_1_x1x1*mapping_v_2*mapping_v_3
                                        mapping_v_x1x2 = mapping_v_1_x1*mapping_v_2_x2*mapping_v_3
                                        mapping_v_x3x3 = mapping_v_1*mapping_v_2*mapping_v_3_x3x3
                                        mapping_v = mapping_v_1*mapping_v_2*mapping_v_3
                                        mapping_v_x3 = mapping_v_1*mapping_v_2*mapping_v_3_x3
                                        mapping_v_x2x3 = mapping_v_1*mapping_v_2_x2*mapping_v_3_x3
                                        mapping_v_x2x2 = mapping_v_1*mapping_v_2_x2x2*mapping_v_3
                                        mapping_v_x1 = mapping_v_1_x1*mapping_v_2*mapping_v_3
                                        x += mapping_v*coeff_x
                                        x_x3 += mapping_v_x3*coeff_x
                                        x_x3x3 += mapping_v_x3x3*coeff_x
                                        x_x2 += mapping_v_x2*coeff_x
                                        x_x2x3 += mapping_v_x2x3*coeff_x
                                        x_x2x2 += mapping_v_x2x2*coeff_x
                                        x_x1 += mapping_v_x1*coeff_x
                                        x_x1x3 += mapping_v_x1x3*coeff_x
                                        x_x1x2 += mapping_v_x1x2*coeff_x
                                        x_x1x1 += mapping_v_x1x1*coeff_x
                                        y += mapping_v*coeff_y
                                        y_x3 += mapping_v_x3*coeff_y
                                        y_x3x3 += mapping_v_x3x3*coeff_y
                                        y_x2 += mapping_v_x2*coeff_y
                                        y_x2x3 += mapping_v_x2x3*coeff_y
                                        y_x2x2 += mapping_v_x2x2*coeff_y
                                        y_x1 += mapping_v_x1*coeff_y
                                        y_x1x3 += mapping_v_x1x3*coeff_y
                                        y_x1x2 += mapping_v_x1x2*coeff_y
                                        y_x1x1 += mapping_v_x1x1*coeff_y
                                        z += mapping_v*coeff_z
                                        z_x3 += mapping_v_x3*coeff_z
                                        z_x3x3 += mapping_v_x3x3*coeff_z
                                        z_x2 += mapping_v_x2*coeff_z
                                        z_x2x3 += mapping_v_x2x3*coeff_z
                                        z_x2x2 += mapping_v_x2x2*coeff_z
                                        z_x1 += mapping_v_x1*coeff_z
                                        z_x1x3 += mapping_v_x1x3*coeff_z
                                        z_x1x2 += mapping_v_x1x2*coeff_z
                                        z_x1x1 += mapping_v_x1x1*coeff_z

                            temp0 = x_x1*y_x2
                            temp1 = x_x2*y_x3
                            temp2 = x_x3*y_x1
                            temp3 = x_x1*y_x3
                            temp4 = x_x2*y_x1
                            temp5 = x_x3*y_x2
                            temp6 = temp0*z_x3 + temp1*z_x1 + temp2*z_x2 - temp3*z_x2 - temp4*z_x3 - temp5*z_x1
                            temp8 = temp0 - temp4
                            temp9 = x_x1x3*y_x2
                            temp10 = x_x2*y_x1x3
                            temp11 = -temp10
                            temp12 = x_x1*y_x2x3 - x_x2x3*y_x1
                            temp13 = temp11 + temp12 + temp9
                            temp14 = temp6**(-1)
                            temp15 = -temp1 + temp5
                            temp16 = x_x2*y_x3x3
                            temp17 = x_x2x3*y_x3
                            temp18 = x_x3*y_x2x3 + x_x3x3*y_x2
                            temp19 = -temp16 - temp17 + temp18
                            temp20 = -temp2 + temp3
                            temp21 = x_x1x3*y_x3
                            temp22 = x_x3*y_x1x3
                            temp23 = temp21 - temp22 + x_x1*y_x3x3 - x_x3x3*y_x1
                            temp25 = y_x3*z_x2x3
                            temp26 = y_x3x3*z_x2
                            temp27 = x_x2*z_x3x3 + x_x2x3*z_x3
                            temp28 = y_x2*z_x3x3 + y_x2x3*z_x3
                            temp29 = x_x3*z_x2x3
                            temp30 = x_x3x3*z_x2
                            temp31 = temp14*(-temp1*z_x1x3 + temp10*z_x3 + temp18*z_x1 + temp21*z_x2 - temp22*z_x2 + temp27*y_x1 - temp28*x_x1 + temp5*z_x1x3 - temp9*z_x3 + x_x1*(temp25 + temp26) - y_x1*(temp29 + temp30) - z_x1*(temp16 + temp17))

                            temp35 = temp14*temp8
                            temp36 = x_x1*z_x3 - x_x3*z_x1
                            temp37 = x_x1x2*z_x3
                            temp38 = x_x3*z_x1x2
                            temp39 = -temp38
                            temp40 = x_x1*z_x2x3 - x_x2x3*z_x1
                            temp41 = temp37 + temp39 + temp40
                            temp42 = x_x2*z_x3
                            temp43 = x_x3*z_x2
                            temp44 = temp42 - temp43
                            temp45 = x_x2x3*z_x2
                            temp46 = x_x3*z_x2x2
                            temp47 = x_x2*z_x2x3 + x_x2x2*z_x3
                            temp48 = -temp45 - temp46 + temp47
                            temp49 = x_x1*z_x2 - x_x2*z_x1
                            temp50 = x_x1x2*z_x2
                            temp51 = temp50 + x_x1*z_x2x2 - x_x2*z_x1x2 - x_x2x2*z_x1
                            temp52 = y_x2x3*z_x2
                            temp53 = y_x3*z_x2x2
                            temp54 = x_x2x3*y_x2 + x_x3*y_x2x2
                            temp55 = y_x2*z_x2x3 + y_x2x2*z_x3
                            temp56 = x_x2*y_x2x3
                            temp57 = x_x2x2*y_x3
                            temp58 = temp14*(-temp1*z_x1x2 - temp37*y_x2 + temp42*y_x1x2 - temp43*y_x1x2 + temp47*y_x1 + temp5*z_x1x2 + temp50*y_x3 + temp54*z_x1 - temp55*x_x1 + x_x1*(temp52 + temp53) - y_x1*(temp45 + temp46) - z_x1*(temp56 + temp57))

                            temp62 = temp14*temp36
                            temp63 = y_x2*z_x3
                            temp64 = y_x3*z_x2
                            temp65 = temp63 - temp64
                            temp66 = y_x1x3*z_x2
                            temp67 = -temp66
                            temp68 = y_x3*z_x1x2
                            temp69 = -temp68
                            temp70 = y_x1x2*z_x3
                            temp71 = y_x2*z_x1x3
                            temp72 = temp70 + temp71
                            temp73 = temp67 + temp69 + temp72
                            temp74 = y_x1*z_x3 - y_x3*z_x1
                            temp75 = y_x1*z_x1x3 + y_x1x1*z_x3 - y_x1x3*z_x1 - y_x3*z_x1x1
                            temp76 = -y_x1*z_x2 + y_x2*z_x1
                            temp77 = -y_x1*z_x1x2 - y_x1x1*z_x2 + y_x1x2*z_x1 + y_x2*z_x1x1
                            temp78 = x_x2*z_x1x3
                            temp79 = temp37 + temp78
                            temp80 = x_x3*y_x1x2
                            temp81 = temp80 + temp9
                            temp82 = x_x1x3*z_x2
                            temp83 = x_x1x2*y_x3
                            temp84 = temp14*(-temp1*z_x1x1 + temp42*y_x1x1 - temp43*y_x1x1 + temp5*z_x1x1 - temp63*x_x1x1 + temp64*x_x1x1 - temp72*x_x1 + temp79*y_x1 + temp81*z_x1 + x_x1*(temp66 + temp68) - y_x1*(temp38 + temp82) - z_x1*(temp10 + temp83))
                            temp88 = temp14*temp65
                            temp89 = x_x1*y_x2x2 + x_x1x2*y_x2 - x_x2*y_x1x2 - x_x2x2*y_x1
                            temp90 = temp54 - temp56 - temp57
                            temp91 = temp12 - temp80 + temp83
                            temp92 = temp14*temp20
                            temp93 = x_x1*z_x3x3 + x_x1x3*z_x3 - x_x3*z_x1x3 - x_x3x3*z_x1
                            temp94 = temp27 - temp29 - temp30
                            temp95 = temp40 - temp78 + temp82
                            temp96 = temp14*temp49
                            temp97 = x_x1*y_x1x2 + x_x1x1*y_x2 - x_x1x2*y_x1 - x_x2*y_x1x1
                            temp98 = temp11 + temp81 - temp83
                            temp99 = x_x1*y_x1x3 + x_x1x1*y_x3 - x_x1x3*y_x1 - x_x3*y_x1x1
                            temp100 = temp14*temp15
                            temp101 = x_x1*z_x1x3 + x_x1x1*z_x3 - x_x1x3*z_x1 - x_x3*z_x1x1
                            temp102 = temp39 + temp79 - temp82
                            temp103 = x_x1*z_x1x2 + x_x1x1*z_x2 - x_x1x2*z_x1 - x_x2*z_x1x1
                            temp104 = temp14*temp44
                            temp105 = -temp25 - temp26 + temp28
                            temp106 = y_x1*z_x3x3 + y_x1x3*z_x3 - y_x3*z_x1x3 - y_x3x3*z_x1
                            temp107 = y_x2x3*z_x1
                            temp108 = y_x1*z_x2x3
                            temp109 = temp107 - temp108 + temp67 + temp71
                            temp110 = temp14*temp76
                            temp111 = -temp52 - temp53 + temp55
                            temp112 = -temp107 + temp108 + temp69 + temp70
                            temp113 = -y_x1*z_x2x2 - y_x1x2*z_x2 + y_x2*z_x1x2 + y_x2x2*z_x1
                            temp114 = temp14*temp74

                            for i_basis_1 in range(0, 5, 1):
                                v_1 = global_test_basis_v_1[i_element_1,i_basis_1,0,i_quad_1]    
                                v_1_x1 = global_test_basis_v_1[i_element_1,i_basis_1,1,i_quad_1]
                                v_1_x1x1 = global_test_basis_v_1[i_element_1,i_basis_1,2,i_quad_1]
                                for i_basis_2 in range(0, 5, 1):
                                    v_2 = global_test_basis_v_2[i_element_2,i_basis_2,0,i_quad_2]
                                    v_2_x2 = global_test_basis_v_2[i_element_2,i_basis_2,1,i_quad_2]
                                    v_2_x2x2 = global_test_basis_v_2[i_element_2,i_basis_2,2,i_quad_2]
                                    for i_basis_3 in range(0, 5, 1):
                                        v_3 = global_test_basis_v_3[i_element_3,i_basis_3,0,i_quad_3]
                                        v_3_x3 = global_test_basis_v_3[i_element_3,i_basis_3,1,i_quad_3]
                                        v_3_x3x3 = global_test_basis_v_3[i_element_3,i_basis_3,2,i_quad_3]
                                        for j_basis_1 in range(0, 5, 1):
                                            u_1 = global_trial_basis_u_1[i_element_1,j_basis_1,0,i_quad_1]
                                            u_1_x1 = global_trial_basis_u_1[i_element_1,j_basis_1,1,i_quad_1]
                                            u_1_x1x1 = global_trial_basis_u_1[i_element_1,j_basis_1,2,i_quad_1]
                                            for j_basis_2 in range(0, 5, 1):
                                                u_2 = global_trial_basis_u_2[i_element_2,j_basis_2,0,i_quad_2]
                                                u_2_x2 = global_trial_basis_u_2[i_element_2,j_basis_2,1,i_quad_2]
                                                u_2_x2x2 = global_trial_basis_u_2[i_element_2,j_basis_2,2,i_quad_2]
                                                for j_basis_3 in range(0, 5, 1):
                                                    u_3 = global_trial_basis_u_3[i_element_3,j_basis_3,0,i_quad_3]
                                                    u_3_x3 = global_trial_basis_u_3[i_element_3,j_basis_3,1,i_quad_3]
                                                    u_3_x3x3 = global_trial_basis_u_3[i_element_3,j_basis_3,2,i_quad_3]

                                                    v = v_1*v_2*v_3
                                                    v_x3 = v_1*v_2*v_3_x3
                                                    v_x3x3 = v_1*v_2*v_3_x3x3
                                                    v_x2 = v_1*v_2_x2*v_3
                                                    v_x2x3 = v_1*v_2_x2*v_3_x3
                                                    v_x2x2 = v_1*v_2_x2x2*v_3
                                                    v_x1 = v_1_x1*v_2*v_3
                                                    v_x1x3 = v_1_x1*v_2*v_3_x3
                                                    v_x1x2 = v_1_x1*v_2_x2*v_3
                                                    v_x1x1 = v_1_x1x1*v_2*v_3
                                                    u = u_1*u_2*u_3
                                                    u_x3 = u_1*u_2*u_3_x3
                                                    u_x3x3 = u_1*u_2*u_3_x3x3
                                                    u_x2 = u_1*u_2_x2*u_3
                                                    u_x2x3 = u_1*u_2_x2*u_3_x3
                                                    u_x2x2 = u_1*u_2_x2x2*u_3
                                                    u_x1 = u_1_x1*u_2*u_3
                                                    u_x1x3 = u_1_x1*u_2*u_3_x3
                                                    u_x1x2 = u_1_x1*u_2_x2*u_3
                                                    u_x1x1 = u_1_x1x1*u_2*u_3

                                                    temp32  = temp8*u_x3
                                                    temp33  = temp15*u_x1
                                                    temp34  = temp20*u_x2
                                                    temp59  = temp36*u_x2
                                                    temp60  = temp44*u_x1
                                                    temp61  = temp49*u_x3
                                                    temp85  = temp65*u_x1
                                                    temp86  = temp74*u_x2
                                                    temp87  = temp76*u_x3
                                                    temp115 = temp8*v_x3
                                                    temp116 = temp15*v_x1
                                                    temp117 = temp20*v_x2
                                                    temp118 = temp36*v_x2
                                                    temp119 = temp44*v_x1
                                                    temp120 = temp49*v_x3
                                                    temp121 = temp65*v_x1
                                                    temp122 = temp74*v_x2
                                                    temp123 = temp76*v_x3

                                                    contribution = temp6*(-temp100*(temp115*temp84 - temp116*temp84 - temp117*temp84 - temp14*(temp15*v_x1x1 + temp98*v_x1) - temp14*(temp20*v_x1x2 + temp99*v_x2) + temp14*(temp8*v_x1x3 + temp97*v_x3)) - temp104*(temp118*temp84 - temp119*temp84 - temp120*temp84 + temp14*(temp101*v_x2 + temp36*v_x1x2) - temp14*(temp102*v_x1 + temp44*v_x1x1) - temp14*(temp103*v_x3 + temp49*v_x1x3)) - temp110*(temp121*temp31 - temp122*temp31 - temp123*temp31 + temp14*(temp105*v_x1 + temp65*v_x1x3) - temp14*(temp106*v_x2 + temp74*v_x2x3) - temp14*(temp109*v_x3 + temp76*v_x3x3)) - temp114*(temp121*temp58 - temp122*temp58 - temp123*temp58 + temp14*(temp111*v_x1 + temp65*v_x1x2) - temp14*(temp112*v_x2 + temp74*v_x2x2) - temp14*(temp113*v_x3 + temp76*v_x2x3)) + temp35*(temp115*temp31 - temp116*temp31 - temp117*temp31 + temp14*(temp13*v_x3 + temp8*v_x3x3) - temp14*(temp15*v_x1x3 + temp19*v_x1) - temp14*(temp20*v_x2x3 + temp23*v_x2)) + temp62*(temp118*temp58 - temp119*temp58 - temp120*temp58 + temp14*(temp36*v_x2x2 + temp41*v_x2) - temp14*(temp44*v_x1x2 + temp48*v_x1) - temp14*(temp49*v_x2x3 + temp51*v_x3)) + temp88*(temp121*temp84 - temp122*temp84 - temp123*temp84 + temp14*(temp65*v_x1x1 + temp73*v_x1) - temp14*(temp74*v_x1x2 + temp75*v_x2) - temp14*(temp76*v_x1x3 + temp77*v_x3)) - temp92*(temp115*temp58 - temp116*temp58 - temp117*temp58 - temp14*(temp15*v_x1x2 + temp90*v_x1) - temp14*(temp20*v_x2x2 + temp91*v_x2) + temp14*(temp8*v_x2x3 + temp89*v_x3)) - temp96*(temp118*temp31 - temp119*temp31 - temp120*temp31 + temp14*(temp36*v_x2x3 + temp93*v_x2) - temp14*(temp44*v_x1x3 + temp94*v_x1) - temp14*(temp49*v_x3x3 + temp95*v_x3)))*(-temp100*(-temp14*(temp15*u_x1x1 + temp98*u_x1) - temp14*(temp20*u_x1x2 + temp99*u_x2) + temp14*(temp8*u_x1x3 + temp97*u_x3) + temp32*temp84 - temp33*temp84 - temp34*temp84) - temp104*(temp14*(temp101*u_x2 + temp36*u_x1x2) - temp14*(temp102*u_x1 + temp44*u_x1x1) - temp14*(temp103*u_x3 + temp49*u_x1x3) + temp59*temp84 - temp60*temp84 - temp61*temp84) - temp110*(temp14*(temp105*u_x1 + temp65*u_x1x3) - temp14*(temp106*u_x2 + temp74*u_x2x3) - temp14*(temp109*u_x3 + temp76*u_x3x3) + temp31*temp85 - temp31*temp86 - temp31*temp87) - temp114*(temp14*(temp111*u_x1 + temp65*u_x1x2) - temp14*(temp112*u_x2 + temp74*u_x2x2) - temp14*(temp113*u_x3 + temp76*u_x2x3) + temp58*temp85 - temp58*temp86 - temp58*temp87) + temp35*(temp14*(temp13*u_x3 + temp8*u_x3x3) - temp14*(temp15*u_x1x3 + temp19*u_x1) - temp14*(temp20*u_x2x3 + temp23*u_x2) + temp31*temp32 - temp31*temp33 - temp31*temp34) + temp62*(temp14*(temp36*u_x2x2 + temp41*u_x2) - temp14*(temp44*u_x1x2 + temp48*u_x1) - temp14*(temp49*u_x2x3 + temp51*u_x3) + temp58*temp59 - temp58*temp60 - temp58*temp61) + temp88*(temp14*(temp65*u_x1x1 + temp73*u_x1) - temp14*(temp74*u_x1x2 + temp75*u_x2) - temp14*(temp76*u_x1x3 + temp77*u_x3) + temp84*temp85 - temp84*temp86 - temp84*temp87) - temp92*(-temp14*(temp15*u_x1x2 + temp90*u_x1) - temp14*(temp20*u_x2x2 + temp91*u_x2) + temp14*(temp8*u_x2x3 + temp89*u_x3) + temp32*temp58 - temp33*temp58 - temp34*temp58) - temp96*(temp14*(temp36*u_x2x3 + temp93*u_x2) - temp14*(temp44*u_x1x3 + temp94*u_x1) - temp14*(temp49*u_x3x3 + temp95*u_x3) + temp31*temp59 - temp31*temp60 - temp31*temp61))

                                                    l_mat_u_v[i_basis_1,i_basis_2,i_basis_3,3 - i_basis_1 + j_basis_1,3 - i_basis_2 + j_basis_2,3 - i_basis_3 + j_basis_3] += contribution
                g_mat_u_v_z6fv8xca[span_v_1:5 + span_v_1,span_v_2:5 + span_v_2,span_v_3:5 + span_v_3,:,:,:] += l_mat_u_v[:,:,:,:,:,:]

    return
