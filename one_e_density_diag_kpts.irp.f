
 BEGIN_PROVIDER [ double precision, one_e_dm_mo_alpha_diag_kpts, (mo_num_per_kpt,kpt_num,N_states) ]
&BEGIN_PROVIDER [ double precision, one_e_dm_mo_beta_diag_kpts, (mo_num_per_kpt,kpt_num,N_states) ]
  implicit none
  BEGIN_DOC
  ! $\alpha$ and $\beta$ one-body density matrix for each state
  ! $\gamma_{\mu\nu} = \langle\Psi|a_{\nu}^{\dagger}a_{\mu}|\Psi\rangle$
  ! $\gamma_{\mu\nu} = \langle a_{\nu} \Psi|a_{\mu} \Psi\rangle$
  ! $\gamma_{\mu\nu} = \sum_{IJ} c^*_J c_I \langle a_{\nu} I|a_{\mu} J\rangle$
  END_DOC
  !todo: implement for kpts
  integer                        :: j,l,m,k_a
  integer                        :: occ(N_int*bit_kind_size,2)
  double precision               :: ck
  integer                        :: kk,k_shft,ii
  integer(bit_kind)              :: tmp_det(N_int,2)
  integer(bit_kind)              :: tmp_det_kpts(N_int,2)
  integer                        :: n_occ(2)
  double precision, allocatable  :: tmp_a(:,:,:), tmp_b(:,:,:)
  integer                        :: krow, kcol

  PROVIDE psi_det psi_coef_complex

  one_e_dm_mo_alpha_diag_kpts = 0.d0
  one_e_dm_mo_beta_diag_kpts  = 0.d0
  !$OMP PARALLEL DEFAULT(NONE)                                      &
      !$OMP PRIVATE(j,k_a,l,m,occ,ck,&
      !$OMP  tmp_a, tmp_b, n_occ, krow, kcol, tmp_det, kk,&
      !$OMP  tmp_det_kpts,k_shft,ii)&
      !$OMP SHARED(psi_det,psi_coef_complex,N_int,N_states,  &
      !$OMP  one_e_dm_mo_alpha_diag_kpts,one_e_dm_mo_beta_diag_kpts,N_det,&
      !$OMP  mo_num_per_kpt,psi_bilinear_matrix_rows,psi_bilinear_matrix_columns,&
      !$OMP  psi_bilinear_matrix_transp_rows, psi_bilinear_matrix_transp_columns,&
      !$OMP  psi_bilinear_matrix_order_reverse, psi_det_alpha_unique, psi_det_beta_unique,&
      !$OMP  psi_bilinear_matrix_values_complex, psi_bilinear_matrix_transp_values_complex,&
      !$OMP  N_det_alpha_unique,N_det_beta_unique,irp_here,kpt_num,kpts_bitmask)
  allocate(tmp_a(mo_num_per_kpt,kpt_num,N_states), tmp_b(mo_num_per_kpt,kpt_num,N_states) )
  tmp_a = 0.d0
  tmp_b = 0.d0
  !$OMP DO SCHEDULE(dynamic,64)
  do k_a=1,N_det
    krow = psi_bilinear_matrix_rows(k_a)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_columns(k_a)
    ASSERT (kcol <= N_det_beta_unique)

    tmp_det(1:N_int,1) = psi_det_alpha_unique(1:N_int,krow)
    tmp_det(1:N_int,2) = psi_det_beta_unique (1:N_int,kcol)

    ! Diagonal part
    ! -------------

    do kk=1,kpt_num
      k_shft = (kk-1)*mo_num_per_kpt
      do ii=1,N_int
        tmp_det_kpts(ii,1) = iand(tmp_det(ii,1),kpts_bitmask(ii,kk))
        tmp_det_kpts(ii,2) = iand(tmp_det(ii,2),kpts_bitmask(ii,kk))
      enddo
      call bitstring_to_list_ab(tmp_det_kpts, occ, n_occ, N_int)
      do m=1,N_states
        ck = cdabs(psi_bilinear_matrix_values_complex(k_a,m)*psi_bilinear_matrix_values_complex(k_a,m))
        !do l=1,elec_alpha_num_kpts(kk)
        do l=1,n_occ(1)
          j = occ(l,1) - k_shft
          tmp_a(j,kk,m) += ck
        enddo
        do l=1,n_occ(2)
          j = occ(l,2) - k_shft
          tmp_b(j,kk,m) += ck
        enddo
      enddo
    enddo

  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  one_e_dm_mo_alpha_diag_kpts(:,:,:) = one_e_dm_mo_alpha_diag_kpts(:,:,:) + tmp_a(:,:,:)
  one_e_dm_mo_beta_diag_kpts(:,:,:)  = one_e_dm_mo_beta_diag_kpts(:,:,:)  + tmp_b(:,:,:)
  !$OMP END CRITICAL
  deallocate(tmp_a,tmp_b)

  !$OMP END PARALLEL

END_PROVIDER

