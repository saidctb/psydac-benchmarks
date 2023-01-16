module hcurl_mass_matrix

  use, intrinsic :: ISO_C_BINDING

  implicit none

  contains

  !........................................
  subroutine assemble_matrix(global_test_basis_TE_0_1, &
      global_test_basis_TE_0_2, global_test_basis_TE_0_3, &
      global_test_basis_TE_1_1, global_test_basis_TE_1_2, &
      global_test_basis_TE_1_3, global_test_basis_TE_2_1, &
      global_test_basis_TE_2_2, global_test_basis_TE_2_3, &
      global_trial_basis_E_0_1, global_trial_basis_E_0_2, &
      global_trial_basis_E_0_3, global_trial_basis_E_1_1, &
      global_trial_basis_E_1_2, global_trial_basis_E_1_3, &
      global_trial_basis_E_2_1, global_trial_basis_E_2_2, &
      global_trial_basis_E_2_3, global_span_TE_0_1, global_span_TE_0_2, &
      global_span_TE_0_3, global_span_TE_1_1, global_span_TE_1_2, &
      global_span_TE_1_3, global_span_TE_2_1, global_span_TE_2_2, &
      global_span_TE_2_3, global_x1, global_x2, global_x3, test_TE_0_p1 &
      , test_TE_0_p2, test_TE_0_p3, test_TE_1_p1, test_TE_1_p2, &
      test_TE_1_p3, test_TE_2_p1, test_TE_2_p2, test_TE_2_p3, &
      trial_E_0_p1, trial_E_0_p2, trial_E_0_p3, trial_E_1_p1, &
      trial_E_1_p2, trial_E_1_p3, trial_E_2_p1, trial_E_2_p2, &
      trial_E_2_p3, n_element_1, n_element_2, n_element_3, k1, k2, k3, &
      pad1, pad2, pad3, g_mat_00, g_mat_11, g_mat_22, coords_from_rank, &
      rank_from_coords, global_thread_starts_1, global_thread_starts_2, &
      global_thread_starts_3, global_thread_ends_1, &
      global_thread_ends_2, global_thread_ends_3, num_threads)

    use omp_lib

    implicit none

    real(C_DOUBLE), intent(in) :: global_test_basis_TE_0_1(0:,0:,0:,0:)
    real(C_DOUBLE), intent(in) :: global_test_basis_TE_0_2(0:,0:,0:,0:)
    real(C_DOUBLE), intent(in) :: global_test_basis_TE_0_3(0:,0:,0:,0:)
    real(C_DOUBLE), intent(in) :: global_test_basis_TE_1_1(0:,0:,0:,0:)
    real(C_DOUBLE), intent(in) :: global_test_basis_TE_1_2(0:,0:,0:,0:)
    real(C_DOUBLE), intent(in) :: global_test_basis_TE_1_3(0:,0:,0:,0:)
    real(C_DOUBLE), intent(in) :: global_test_basis_TE_2_1(0:,0:,0:,0:)
    real(C_DOUBLE), intent(in) :: global_test_basis_TE_2_2(0:,0:,0:,0:)
    real(C_DOUBLE), intent(in) :: global_test_basis_TE_2_3(0:,0:,0:,0:)
    real(C_DOUBLE), intent(in) :: global_trial_basis_E_0_1(0:,0:,0:,0:)
    real(C_DOUBLE), intent(in) :: global_trial_basis_E_0_2(0:,0:,0:,0:)
    real(C_DOUBLE), intent(in) :: global_trial_basis_E_0_3(0:,0:,0:,0:)
    real(C_DOUBLE), intent(in) :: global_trial_basis_E_1_1(0:,0:,0:,0:)
    real(C_DOUBLE), intent(in) :: global_trial_basis_E_1_2(0:,0:,0:,0:)
    real(C_DOUBLE), intent(in) :: global_trial_basis_E_1_3(0:,0:,0:,0:)
    real(C_DOUBLE), intent(in) :: global_trial_basis_E_2_1(0:,0:,0:,0:)
    real(C_DOUBLE), intent(in) :: global_trial_basis_E_2_2(0:,0:,0:,0:)
    real(C_DOUBLE), intent(in) :: global_trial_basis_E_2_3(0:,0:,0:,0:)
    integer(C_INT64_T), intent(in) :: global_span_TE_0_1(0:)
    integer(C_INT64_T), intent(in) :: global_span_TE_0_2(0:)
    integer(C_INT64_T), intent(in) :: global_span_TE_0_3(0:)
    integer(C_INT64_T), intent(in) :: global_span_TE_1_1(0:)
    integer(C_INT64_T), intent(in) :: global_span_TE_1_2(0:)
    integer(C_INT64_T), intent(in) :: global_span_TE_1_3(0:)
    integer(C_INT64_T), intent(in) :: global_span_TE_2_1(0:)
    integer(C_INT64_T), intent(in) :: global_span_TE_2_2(0:)
    integer(C_INT64_T), intent(in) :: global_span_TE_2_3(0:)
    real(C_DOUBLE), intent(in) :: global_x1(0:,0:)
    real(C_DOUBLE), intent(in) :: global_x2(0:,0:)
    real(C_DOUBLE), intent(in) :: global_x3(0:,0:)
    integer(C_INT64_T), value :: test_TE_0_p1
    integer(C_INT64_T), value :: test_TE_0_p2
    integer(C_INT64_T), value :: test_TE_0_p3
    integer(C_INT64_T), value :: test_TE_1_p1
    integer(C_INT64_T), value :: test_TE_1_p2
    integer(C_INT64_T), value :: test_TE_1_p3
    integer(C_INT64_T), value :: test_TE_2_p1
    integer(C_INT64_T), value :: test_TE_2_p2
    integer(C_INT64_T), value :: test_TE_2_p3
    integer(C_INT64_T), value :: trial_E_0_p1
    integer(C_INT64_T), value :: trial_E_0_p2
    integer(C_INT64_T), value :: trial_E_0_p3
    integer(C_INT64_T), value :: trial_E_1_p1
    integer(C_INT64_T), value :: trial_E_1_p2
    integer(C_INT64_T), value :: trial_E_1_p3
    integer(C_INT64_T), value :: trial_E_2_p1
    integer(C_INT64_T), value :: trial_E_2_p2
    integer(C_INT64_T), value :: trial_E_2_p3
    integer(C_INT64_T), value :: n_element_1
    integer(C_INT64_T), value :: n_element_2
    integer(C_INT64_T), value :: n_element_3
    integer(C_INT64_T), value :: k1
    integer(C_INT64_T), value :: k2
    integer(C_INT64_T), value :: k3
    integer(C_INT64_T), value :: pad1
    integer(C_INT64_T), value :: pad2
    integer(C_INT64_T), value :: pad3
    real(C_DOUBLE), intent(inout) :: g_mat_00(0:,0:,0:,0:,0:,0:)
    real(C_DOUBLE), intent(inout) :: g_mat_11(0:,0:,0:,0:,0:,0:)
    real(C_DOUBLE), intent(inout) :: g_mat_22(0:,0:,0:,0:,0:,0:)
    integer(C_INT64_T), intent(in) :: coords_from_rank(0:,0:)
    integer(C_INT64_T), intent(in) :: rank_from_coords(0:,0:,0:)
    integer(C_INT64_T), intent(in) :: global_thread_starts_1(0:)
    integer(C_INT64_T), intent(in) :: global_thread_starts_2(0:)
    integer(C_INT64_T), intent(in) :: global_thread_starts_3(0:)
    integer(C_INT64_T), intent(in) :: global_thread_ends_1(0:)
    integer(C_INT64_T), intent(in) :: global_thread_ends_2(0:)
    integer(C_INT64_T), intent(in) :: global_thread_ends_3(0:)
    integer(C_INT64_T), value :: num_threads
    real(C_DOUBLE), allocatable :: l_mat_00(:,:,:,:,:,:)
    real(C_DOUBLE), allocatable :: l_mat_11(:,:,:,:,:,:)
    real(C_DOUBLE), allocatable :: l_mat_22(:,:,:,:,:,:)
    integer(C_INT64_T) :: thread_id
    integer(C_INT64_T) :: thread_coords_1
    integer(C_INT64_T) :: thread_coords_2
    integer(C_INT64_T) :: thread_coords_3
    integer(C_INT64_T) :: global_thread_size_1
    integer(C_INT64_T) :: global_thread_size_2
    integer(C_INT64_T) :: global_thread_size_3
    integer(C_INT64_T), allocatable :: local_thread_starts_1(:)
    integer(C_INT64_T), allocatable :: local_thread_starts_2(:)
    integer(C_INT64_T), allocatable :: local_thread_starts_3(:)
    integer(C_INT64_T), allocatable :: local_thread_ends_1(:)
    integer(C_INT64_T), allocatable :: local_thread_ends_2(:)
    integer(C_INT64_T), allocatable :: local_thread_ends_3(:)
    real(C_DOUBLE), allocatable :: local_test_basis_TE_0_1(:,:,:,:)
    real(C_DOUBLE), allocatable :: local_test_basis_TE_1_1(:,:,:,:)
    real(C_DOUBLE), allocatable :: local_test_basis_TE_2_1(:,:,:,:)
    real(C_DOUBLE), allocatable :: local_test_basis_TE_0_2(:,:,:,:)
    real(C_DOUBLE), allocatable :: local_test_basis_TE_1_2(:,:,:,:)
    real(C_DOUBLE), allocatable :: local_test_basis_TE_2_2(:,:,:,:)
    real(C_DOUBLE), allocatable :: local_test_basis_TE_0_3(:,:,:,:)
    real(C_DOUBLE), allocatable :: local_test_basis_TE_1_3(:,:,:,:)
    real(C_DOUBLE), allocatable :: local_test_basis_TE_2_3(:,:,:,:)
    real(C_DOUBLE), allocatable :: local_trial_basis_E_0_1(:,:,:,:)
    real(C_DOUBLE), allocatable :: local_trial_basis_E_1_1(:,:,:,:)
    real(C_DOUBLE), allocatable :: local_trial_basis_E_2_1(:,:,:,:)
    real(C_DOUBLE), allocatable :: local_trial_basis_E_0_2(:,:,:,:)
    real(C_DOUBLE), allocatable :: local_trial_basis_E_1_2(:,:,:,:)
    real(C_DOUBLE), allocatable :: local_trial_basis_E_2_2(:,:,:,:)
    real(C_DOUBLE), allocatable :: local_trial_basis_E_0_3(:,:,:,:)
    real(C_DOUBLE), allocatable :: local_trial_basis_E_1_3(:,:,:,:)
    real(C_DOUBLE), allocatable :: local_trial_basis_E_2_3(:,:,:,:)
    integer(C_INT64_T), allocatable :: local_span_TE_0_1(:)
    integer(C_INT64_T), allocatable :: local_span_TE_1_1(:)
    integer(C_INT64_T), allocatable :: local_span_TE_2_1(:)
    integer(C_INT64_T), allocatable :: local_span_TE_0_2(:)
    integer(C_INT64_T), allocatable :: local_span_TE_1_2(:)
    integer(C_INT64_T), allocatable :: local_span_TE_2_2(:)
    integer(C_INT64_T), allocatable :: local_span_TE_0_3(:)
    integer(C_INT64_T), allocatable :: local_span_TE_1_3(:)
    integer(C_INT64_T), allocatable :: local_span_TE_2_3(:)
    integer(C_INT64_T) :: local_i_element_1
    integer(C_INT64_T) :: local_i_element_2
    integer(C_INT64_T) :: local_i_element_3
    integer(C_INT64_T) :: i_element_1
    integer(C_INT64_T) :: span_TE_0_1
    integer(C_INT64_T) :: span_TE_1_1
    integer(C_INT64_T) :: span_TE_2_1
    integer(C_INT64_T) :: i_element_2
    integer(C_INT64_T) :: span_TE_0_2
    integer(C_INT64_T) :: span_TE_1_2
    integer(C_INT64_T) :: span_TE_2_2
    integer(C_INT64_T) :: i_element_3
    integer(C_INT64_T) :: span_TE_0_3
    integer(C_INT64_T) :: span_TE_1_3
    integer(C_INT64_T) :: span_TE_2_3
    integer(C_INT64_T) :: i_basis_1
    integer(C_INT64_T) :: i_basis_2
    integer(C_INT64_T) :: i_basis_3
    integer(C_INT64_T) :: j_basis_1
    integer(C_INT64_T) :: j_basis_2
    integer(C_INT64_T) :: j_basis_3
    real(C_DOUBLE) :: contribution
    integer(C_INT64_T) :: i_quad_1
    real(C_DOUBLE) :: TE_0_1
    real(C_DOUBLE) :: E_0_1
    integer(C_INT64_T) :: i_quad_2
    real(C_DOUBLE) :: TE_0_2
    real(C_DOUBLE) :: E_0_2
    integer(C_INT64_T) :: i_quad_3
    real(C_DOUBLE) :: TE_0_3
    real(C_DOUBLE) :: E_0_3
    real(C_DOUBLE) :: TE_0
    real(C_DOUBLE) :: E_0
    real(C_DOUBLE) :: TE_1_1
    real(C_DOUBLE) :: E_1_1
    real(C_DOUBLE) :: TE_1_2
    real(C_DOUBLE) :: E_1_2
    real(C_DOUBLE) :: TE_1_3
    real(C_DOUBLE) :: E_1_3
    real(C_DOUBLE) :: TE_1
    real(C_DOUBLE) :: E_1
    real(C_DOUBLE) :: TE_2_1
    real(C_DOUBLE) :: E_2_1
    real(C_DOUBLE) :: TE_2_2
    real(C_DOUBLE) :: E_2_2
    real(C_DOUBLE) :: TE_2_3
    real(C_DOUBLE) :: E_2_3
    real(C_DOUBLE) :: TE_2
    real(C_DOUBLE) :: E_2
    real(C_DOUBLE) :: start, endd
    integer(C_INT64_T) :: i_0005
    integer(C_INT64_T) :: i_0006
    integer(C_INT64_T) :: i_0007

    !$omp parallel default(private) shared(coords_from_rank, rank_from_coords, global_thread_starts_1, global_thread_starts_2, global_thread_starts_3, global_thread_ends_1, global_thread_ends_2, global_thread_ends_3, global_test_basis_TE_0_1, global_test_basis_TE_0_2, global_test_basis_TE_0_3, global_test_basis_TE_1_1, global_test_basis_TE_1_2, global_test_basis_TE_1_3, global_test_basis_TE_2_1, global_test_basis_TE_2_2, global_test_basis_TE_2_3, global_trial_basis_E_0_1, global_trial_basis_E_0_2, global_trial_basis_E_0_3, global_trial_basis_E_1_1, global_trial_basis_E_1_2, global_trial_basis_E_1_3, global_trial_basis_E_2_1, global_trial_basis_E_2_2, global_trial_basis_E_2_3, global_span_TE_0_1, global_span_TE_0_2, global_span_TE_0_3, global_span_TE_1_1, global_span_TE_1_2, global_span_TE_1_3, global_span_TE_2_1, global_span_TE_2_2, global_span_TE_2_3, global_x1, global_x2, global_x3, g_mat_00, g_mat_11, g_mat_22) firstprivate(test_TE_0_p1, test_TE_0_p2, test_TE_0_p3, test_TE_1_p1, test_TE_1_p2, test_TE_1_p3, test_TE_2_p1, test_TE_2_p2, test_TE_2_p3, trial_E_0_p1, trial_E_0_p2, trial_E_0_p3, trial_E_1_p1, trial_E_1_p2, trial_E_1_p3, trial_E_2_p1, trial_E_2_p2, trial_E_2_p3, n_element_1, n_element_2, n_element_3, k1, k2, k3, pad1, pad2, pad3, num_threads)   
    allocate(l_mat_00(0:8, 0:8, 0:6, 0:4, 0:4, 0:3))
    l_mat_00 = 0.0
    allocate(l_mat_11(0:8, 0:6, 0:8, 0:4, 0:3, 0:4))
    l_mat_11 = 0.0
    allocate(l_mat_22(0:6, 0:8, 0:8, 0:3, 0:4, 0:4))
    l_mat_22 = 0.0
    thread_id = omp_get_thread_num()

    thread_coords_1 = coords_from_rank(0, thread_id)
    thread_coords_2 = coords_from_rank(1, thread_id)
    thread_coords_3 = coords_from_rank(2, thread_id)

    global_thread_size_1 = 1 + global_thread_ends_1(thread_coords_1) - global_thread_starts_1(thread_coords_1)
    global_thread_size_2 = 1 + global_thread_ends_2(thread_coords_2) - global_thread_starts_2(thread_coords_2)
    global_thread_size_3 = 1 + global_thread_ends_3(thread_coords_3) - global_thread_starts_3(thread_coords_3)

    allocate(local_thread_starts_1(0:1))
    local_thread_starts_1 = [0_C_INT64_T, Int(1.0 /2.0 * global_thread_size_1, C_INT64_T)]

    allocate(local_thread_starts_2(0:1))
    local_thread_starts_2 = [0_C_INT64_T, Int(1.0 / 2.0 * global_thread_size_2, C_INT64_T)]

    allocate(local_thread_starts_3(0:1))
    local_thread_starts_3 = [0_C_INT64_T, Int(1.0 / 2.0 * global_thread_size_3, C_INT64_T)]

    allocate(local_thread_ends_1(0:1))
    local_thread_ends_1 = [Int(1.0 / 2.0 * global_thread_size_1, C_INT64_T), global_thread_size_1]

    allocate(local_thread_ends_2(0:1))
    local_thread_ends_2 = [Int(1.0 / 2.0 * global_thread_size_2, C_INT64_T), global_thread_size_2]

    allocate(local_thread_ends_3(0:1))
    local_thread_ends_3 = [Int(1.0 / 2.0*global_thread_size_3, C_INT64_T), global_thread_size_3]

    allocate(local_test_basis_TE_0_1(0:4, 0:1, 0:3, 0:global_thread_size_1-1))
    allocate(local_test_basis_TE_1_1(0:4, 0:1, 0:4, 0:global_thread_size_1-1))
    allocate(local_test_basis_TE_2_1(0:4, 0:1, 0:4, 0:global_thread_size_1-1))
    local_test_basis_TE_0_1(:,:,:,:) = 0.0
    local_test_basis_TE_1_1(:,:,:,:) = 0.0
    local_test_basis_TE_2_1(:,:,:,:) = 0.0

!    local_test_basis_TE_0_1(:,:,:,:) = global_test_basis_TE_0_1(:,:,:, global_thread_starts_1(thread_coords_1):global_thread_starts_1(thread_coords_1)+size(local_test_basis_TE_0_1, 4, C_INT64_T) - 1)
!    local_test_basis_TE_1_1(:,:,:,:) = global_test_basis_TE_1_1(:,:,:, global_thread_starts_1(thread_coords_1):global_thread_starts_1(thread_coords_1)+size(local_test_basis_TE_1_1, 4, C_INT64_T) - 1)
!    local_test_basis_TE_2_1(:,:,:,:) = global_test_basis_TE_2_1(:,:,:, global_thread_starts_1(thread_coords_1):global_thread_starts_1(thread_coords_1)+size(local_test_basis_TE_2_1, 4, C_INT64_T) - 1)

    allocate(local_test_basis_TE_0_2(0:4, 0:1, 0:4, 0:global_thread_size_2 - 1))
    allocate(local_test_basis_TE_1_2(0:4, 0:1, 0:3, 0:global_thread_size_2 - 1))
    allocate(local_test_basis_TE_2_2(0:4, 0:1, 0:4, 0:global_thread_size_2 - 1))
    local_test_basis_TE_0_2(:,:,:,:) = 0.0
    local_test_basis_TE_1_2(:,:,:,:) = 0.0
    local_test_basis_TE_2_2(:,:,:,:) = 0.0

!    local_test_basis_TE_0_2(:,:,:,:) = global_test_basis_TE_0_2(:,:,:, global_thread_starts_2(thread_coords_2):global_thread_starts_2(thread_coords_2)+size(local_test_basis_TE_0_2, 4,C_INT64_T) - 1)
!    local_test_basis_TE_1_2(:,:,:,:) = global_test_basis_TE_1_2(:,:,:, global_thread_starts_2(thread_coords_2):global_thread_starts_2(thread_coords_2)+size(local_test_basis_TE_1_2, 4,C_INT64_T) - 1)
!    local_test_basis_TE_2_2(:,:,:,:) = global_test_basis_TE_2_2(:,:,:, global_thread_starts_2(thread_coords_2):global_thread_starts_2(thread_coords_2)+size(local_test_basis_TE_2_2, 4,C_INT64_T) - 1)

    allocate(local_test_basis_TE_0_3(0:4, 0:1, 0:4, 0:global_thread_size_3 - 1))
    allocate(local_test_basis_TE_1_3(0:4, 0:1, 0:4, 0:global_thread_size_3 - 1))
    allocate(local_test_basis_TE_2_3(0:4, 0:1, 0:3, 0:global_thread_size_3 - 1))
    local_test_basis_TE_0_3(:,:,:,:) = 0.0
    local_test_basis_TE_1_3(:,:,:,:) = 0.0
    local_test_basis_TE_2_3(:,:,:,:) = 0.0

!    local_test_basis_TE_0_3(:,:,:,:) = global_test_basis_TE_0_3(:,:,:, global_thread_starts_3(thread_coords_3):global_thread_starts_3(thread_coords_3)+size(local_test_basis_TE_0_3, 4,C_INT64_T) - 1)
!    local_test_basis_TE_1_3(:,:,:,:) = global_test_basis_TE_1_3(:,:,:, global_thread_starts_3(thread_coords_3):global_thread_starts_3(thread_coords_3)+size(local_test_basis_TE_1_3, 4,C_INT64_T) - 1)
!    local_test_basis_TE_2_3(:,:,:,:) = global_test_basis_TE_2_3(:,:,:, global_thread_starts_3(thread_coords_3):global_thread_starts_3(thread_coords_3)+size(local_test_basis_TE_2_3, 4,C_INT64_T) - 1)

    allocate(local_trial_basis_E_0_1(0:4, 0:1, 0:3, 0:global_thread_size_1 - 1))
    allocate(local_trial_basis_E_1_1(0:4, 0:1, 0:4, 0:global_thread_size_1 - 1))
    allocate(local_trial_basis_E_2_1(0:4,0:1,0:4,0:global_thread_size_1 - 1))
    local_trial_basis_E_0_1(:,:,:,:) = 0.0
    local_trial_basis_E_1_1(:,:,:,:) = 0.0
    local_trial_basis_E_2_1(:,:,:,:) = 0.0

!    local_trial_basis_E_0_1(:,:,:,:) = global_trial_basis_E_0_1(:,:,:, global_thread_starts_1(thread_coords_1):global_thread_starts_1(thread_coords_1)+size(local_trial_basis_E_0_1, 4,C_INT64_T)-1)
!    local_trial_basis_E_1_1(:,:,:,:) = global_trial_basis_E_1_1(:,:,:, global_thread_starts_1(thread_coords_1):global_thread_starts_1(thread_coords_1)+size(local_trial_basis_E_1_1, 4,C_INT64_T)-1)
!    local_trial_basis_E_2_1(:,:,:,:) = global_trial_basis_E_2_1(:,:,:, global_thread_starts_1(thread_coords_1):global_thread_starts_1(thread_coords_1)+size(local_trial_basis_E_2_1, 4,C_INT64_T)-1)

    allocate(local_trial_basis_E_0_2(0:4,0:1,0:4, 0:global_thread_size_2-1))
    allocate(local_trial_basis_E_1_2(0:4,0:1,0:3, 0:global_thread_size_2-1))
    allocate(local_trial_basis_E_2_2(0:4,0:1,0:4, 0:global_thread_size_2-1))
    local_trial_basis_E_0_2(:,:,:,:) = 0.0
    local_trial_basis_E_1_2(:,:,:,:) = 0.0
    local_trial_basis_E_2_2(:,:,:,:) = 0.0

!    local_trial_basis_E_0_2(:,:,:,:) = global_trial_basis_E_0_2(:,:,:, global_thread_starts_2(thread_coords_2):global_thread_starts_2(thread_coords_2)+size(local_trial_basis_E_0_2, 4,C_INT64_T) - 1)
!    local_trial_basis_E_1_2(:,:,:,:) = global_trial_basis_E_1_2(:,:,:, global_thread_starts_2(thread_coords_2):global_thread_starts_2(thread_coords_2)+size(local_trial_basis_E_1_2, 4,C_INT64_T) - 1)
!    local_trial_basis_E_2_2(:,:,:,:) = global_trial_basis_E_2_2(:,:,:, global_thread_starts_2(thread_coords_2):global_thread_starts_2(thread_coords_2)+size(local_trial_basis_E_2_2, 4,C_INT64_T) - 1)

    allocate(local_trial_basis_E_0_3(0:4,0:1,0:4,0:global_thread_size_3-1))
    allocate(local_trial_basis_E_1_3(0:4,0:1,0:4,0:global_thread_size_3-1))
    allocate(local_trial_basis_E_2_3(0:4,0:1,0:3,0:global_thread_size_3-1))
    local_trial_basis_E_0_3(:,:,:,:) = 0.0
    local_trial_basis_E_1_3(:,:,:,:) = 0.0
    local_trial_basis_E_2_3(:,:,:,:) = 0.0

!    local_trial_basis_E_0_3(:,:,:,:) = global_trial_basis_E_0_3(:,:,:, global_thread_starts_3(thread_coords_3):global_thread_starts_3(thread_coords_3)+size(local_trial_basis_E_0_3, 4,C_INT64_T)-1)
!    local_trial_basis_E_1_3(:,:,:,:) = global_trial_basis_E_1_3(:,:,:, global_thread_starts_3(thread_coords_3):global_thread_starts_3(thread_coords_3)+size(local_trial_basis_E_1_3, 4,C_INT64_T)-1)
!    local_trial_basis_E_2_3(:,:,:,:) = global_trial_basis_E_2_3(:,:,:, global_thread_starts_3(thread_coords_3):global_thread_starts_3(thread_coords_3)+size(local_trial_basis_E_2_3, 4,C_INT64_T)-1)

    allocate(local_span_TE_0_1(0:global_thread_size_1-1))
    allocate(local_span_TE_1_1(0:global_thread_size_1-1))
    allocate(local_span_TE_2_1(0:global_thread_size_1-1))
    local_span_TE_0_1(:) = global_span_TE_0_1(global_thread_starts_1(thread_coords_1):1 + global_thread_ends_1(thread_coords_1)-1)
    local_span_TE_1_1(:) = global_span_TE_1_1(global_thread_starts_1(thread_coords_1):1 + global_thread_ends_1(thread_coords_1)-1)
    local_span_TE_2_1(:) = global_span_TE_2_1(global_thread_starts_1(thread_coords_1):1 + global_thread_ends_1(thread_coords_1)-1)

    allocate(local_span_TE_0_2(0:global_thread_size_2-1))
    allocate(local_span_TE_1_2(0:global_thread_size_2-1))
    allocate(local_span_TE_2_2(0:global_thread_size_2-1))
    local_span_TE_0_2(:) = global_span_TE_0_2(global_thread_starts_2(thread_coords_2):1 + global_thread_ends_2(thread_coords_2)-1)
    local_span_TE_1_2(:) = global_span_TE_1_2(global_thread_starts_2(thread_coords_2):1 + global_thread_ends_2(thread_coords_2)-1)
    local_span_TE_2_2(:) = global_span_TE_2_2(global_thread_starts_2(thread_coords_2):1 + global_thread_ends_2(thread_coords_2)-1)

    allocate(local_span_TE_0_3(0:global_thread_size_3-1))
    allocate(local_span_TE_1_3(0:global_thread_size_3-1))
    allocate(local_span_TE_2_3(0:global_thread_size_3-1))
    local_span_TE_0_3(:) = global_span_TE_0_3(global_thread_starts_3(thread_coords_3):1 + global_thread_ends_3(thread_coords_3)-1)
    local_span_TE_1_3(:) = global_span_TE_1_3(global_thread_starts_3(thread_coords_3):1 + global_thread_ends_3(thread_coords_3)-1)
    local_span_TE_2_3(:) = global_span_TE_2_3(global_thread_starts_3(thread_coords_3):1 + global_thread_ends_3(thread_coords_3)-1)

    start = omp_get_wtime()

    do local_i_element_1 = 0,1,1
      do local_i_element_2 = 0,1,1
        do local_i_element_3 = 0,1,1
          !$omp barrier

          do i_element_1 = local_thread_starts_1(local_i_element_1), local_thread_ends_1(local_i_element_1) - 1, 1
            span_TE_0_1 = local_span_TE_0_1(i_element_1)
            span_TE_1_1 = local_span_TE_1_1(i_element_1)
            span_TE_2_1 = local_span_TE_2_1(i_element_1)
            do i_element_2 = local_thread_starts_2(local_i_element_2), local_thread_ends_2(local_i_element_2) - 1, 1
              span_TE_0_2 = local_span_TE_0_2(i_element_2)
              span_TE_1_2 = local_span_TE_1_2(i_element_2)
              span_TE_2_2 = local_span_TE_2_2(i_element_2)
              do i_element_3 = local_thread_starts_3(local_i_element_3), local_thread_ends_3(local_i_element_3) - 1, 1
                span_TE_0_3 = local_span_TE_0_3(i_element_3)
                span_TE_1_3 = local_span_TE_1_3(i_element_3)
                span_TE_2_3 = local_span_TE_2_3(i_element_3)

                do i_basis_1 = 0,3,1
                  do i_basis_2 = 0,4,1
                    do i_basis_3 = 0,4,1
                      do j_basis_1 = 0,3,1
                        do j_basis_2 = 0,4,1
                          do j_basis_3 = 0,4,1
                            contribution = 0.0
                            do i_quad_1 = 0,4,1
                              TE_0_1 = local_test_basis_TE_0_1(i_quad_1,0, i_basis_1, i_element_1)
                              E_0_1 = local_trial_basis_E_0_1(i_quad_1,0, j_basis_1, i_element_1)
                              do i_quad_2 = 0,4,1
                                TE_0_2 = local_test_basis_TE_0_2(i_quad_2, 0, i_basis_2, i_element_2)
                                E_0_2 = local_trial_basis_E_0_2(i_quad_2, 0, j_basis_2, i_element_2)
                                do i_quad_3 = 0,4,1
                                  TE_0_3 = local_test_basis_TE_0_3( i_quad_3, 0, i_basis_3, i_element_3)
                                  E_0_3 = local_trial_basis_E_0_3(i_quad_3, 0, j_basis_3, i_element_3)
                                  TE_0 = TE_0_1 * TE_0_2 * TE_0_3
                                  E_0 = E_0_1 * E_0_2 * E_0_3
                                  contribution = contribution + E_0 *TE_0
                                end do
                              end do
                            end do
                            l_mat_00(4 - i_basis_3 + j_basis_3, 4 - i_basis_2 + j_basis_2, 3- i_basis_1 +j_basis_1, i_basis_3, i_basis_2, i_basis_1) = contribution
                          end do
                        end do
                      end do
                    end do
                  end do
                end do

                do i_basis_1 = 0,4,1
                  do i_basis_2 = 0,3,1
                    do i_basis_3 = 0,4,1
                      do j_basis_1 = 0,4,1
                        do j_basis_2 = 0,3,1
                          do j_basis_3 = 0, 4,1
                            contribution = 0.0
                            do i_quad_1 = 0,4,1
                              TE_1_1 = local_test_basis_TE_1_1(i_quad_1,0, i_basis_1, i_element_1)
                              E_1_1 = local_trial_basis_E_1_1(i_quad_1,0, j_basis_1, i_element_1)
                              do i_quad_2 = 0,4,1
                                TE_1_2 = local_test_basis_TE_1_2(i_quad_2, 0, i_basis_2, i_element_2)
                                E_1_2 = local_trial_basis_E_1_2(i_quad_2, 0, j_basis_2, i_element_2)
                                do i_quad_3 = 0,4,1
                                  TE_1_3 = local_test_basis_TE_1_3(i_quad_3, 0, i_basis_3, i_element_3)
                                  E_1_3 = local_trial_basis_E_1_3(i_quad_3, 0, j_basis_3, i_element_3)
                                  TE_1 = TE_1_1 * TE_1_2 * TE_1_3
                                  E_1 = E_1_1 * E_1_2 * E_1_3
                                  contribution = contribution + E_1*TE_1
                                end do
                              end do
                            end do
                            l_mat_11(4 - i_basis_3 + j_basis_3, 3 - i_basis_2 + j_basis_2, 4 - i_basis_1 +j_basis_1, i_basis_3, i_basis_2, i_basis_1) = contribution
                          end do
                        end do
                      end do
                    end do
                  end do
                end do

                do i_basis_1 = 0,4,1
                  do i_basis_2 = 0,4,1
                    do i_basis_3 = 0,3,1
                      do j_basis_1 = 0,4,1
                        do j_basis_2 = 0,4,1
                          do j_basis_3 = 0,3,1
                            contribution = 0.0
                            do i_quad_1 = 0,4,1
                              TE_2_1 = local_test_basis_TE_2_1(i_quad_1,0, i_basis_1, i_element_1)
                              E_2_1 = local_trial_basis_E_2_1(i_quad_1,0, j_basis_1, i_element_1)
                              do i_quad_2 = 0,4,1
                                TE_2_2 = local_test_basis_TE_2_2( i_quad_2, 0, i_basis_2, i_element_2)
                                E_2_2 = local_trial_basis_E_2_2(i_quad_2,0, j_basis_2, i_element_2)
                                do i_quad_3 = 0,4,1
                                  TE_2_3 = local_test_basis_TE_2_3( i_quad_3, 0, i_basis_3, i_element_3)
                                  E_2_3 = local_trial_basis_E_2_3( i_quad_3, 0, j_basis_3, i_element_3)
                                  TE_2 = TE_2_1 * TE_2_2 * TE_2_3
                                  E_2 = E_2_1 * E_2_2 * E_2_3
                                  contribution = contribution + E_2*TE_2
                                end do
                              end do
                            end do
                            l_mat_22(3 - i_basis_3 + j_basis_3, 4 - i_basis_2 + j_basis_2, 4 - i_basis_1 + j_basis_1, i_basis_3, i_basis_2, i_basis_1) = contribution
                          end do
                        end do
                      end do
                    end do
                  end do
                end do

                do i_0005 = 0, 1 + pad1 + span_TE_0_1 - (pad1 + span_TE_0_1 - test_TE_0_p1) - 1,1
                  do i_0006 = 0, 1 + pad2 + span_TE_0_2 - (pad2 + span_TE_0_2 - test_TE_0_p2) - 1,1
                    do i_0007 = 0, 1 + pad3 +span_TE_0_3 - (pad3 + span_TE_0_3 - test_TE_0_p3) - 1, 1
                            g_mat_00(:,:,:, i_0007 + (pad3 + span_TE_0_3 - test_TE_0_p3), i_0006 + (pad2 + span_TE_0_2 - test_TE_0_p2), i_0005 + (pad1 + span_TE_0_1 - test_TE_0_p1)) = g_mat_00(:,:,:, i_0007 + (pad3 + span_TE_0_3 - test_TE_0_p3), i_0006 + (pad2 + span_TE_0_2 - test_TE_0_p2), i_0005 + (pad1 + span_TE_0_1 - test_TE_0_p1)) + l_mat_00(:,:,:, i_0007, i_0006, i_0005)
                          end do
                        end do
                      end do

                do i_0005 = 0, 1 + pad1 +span_TE_1_1 - (pad1 + span_TE_1_1 - test_TE_1_p1) - 1,1
                  do i_0006 = 0, 1 + pad2 + span_TE_1_2 - (pad2 + span_TE_1_2 - test_TE_1_p2) - 1,1
                    do i_0007 = 0, 1 + pad3 +span_TE_1_3 - (pad3 + span_TE_1_3 - test_TE_1_p3) - 1,1
                            g_mat_11(:,:,:, i_0007 + (pad3 + span_TE_1_3 - test_TE_1_p3), i_0006 + (pad2 + span_TE_1_2 - test_TE_1_p2), i_0005 + (pad1 + span_TE_1_1 - test_TE_1_p1)) = g_mat_11(:,:,:, i_0007 + (pad3 + span_TE_1_3 -test_TE_1_p3), i_0006 + (pad2 + span_TE_1_2 - test_TE_1_p2), i_0005 + (pad1 + span_TE_1_1 - test_TE_1_p1)) + l_mat_11(:,:,:, i_0007, i_0006, i_0005)
                          end do
                        end do
                      end do

                do i_0005 = 0, 1 + pad1 +span_TE_2_1 - (pad1 + span_TE_2_1 - test_TE_2_p1) - 1, 1
                  do i_0006 = 0, 1 + pad2 +span_TE_2_2 - (pad2 + span_TE_2_2 - test_TE_2_p2) - 1,1
                    do i_0007 = 0, 1 + pad3 +span_TE_2_3 - (pad3 + span_TE_2_3 - test_TE_2_p3) - 1,1
                            g_mat_22(:, :, :, i_0007 + (pad3 + span_TE_2_3 - test_TE_2_p3), i_0006 + (pad2 + span_TE_2_2 - test_TE_2_p2), i_0005 + (pad1 + span_TE_2_1 - test_TE_2_p1)) = g_mat_22(:,:,:, i_0007 + (pad3 + span_TE_2_3 - test_TE_2_p3), i_0006 + (pad2 + span_TE_2_2 - test_TE_2_p2), i_0005 + (pad1 + span_TE_2_1 - test_TE_2_p1)) + l_mat_22(:,:,:, i_0007, i_0006, i_0005)
                          end do
                        end do
                      end do
              end do
            end do
          end do
        end do
      end do
    end do
    endd = omp_get_wtime() 
    print *, "timing : ", endd-start
    !$omp end parallel
    return

  end subroutine assemble_matrix
  !........................................

end module hcurl_mass_matrix
