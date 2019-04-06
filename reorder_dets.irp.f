program print_h0j
  use bitmasks
  implicit none
  BEGIN_DOC
! Print 1st line of H matrix
  END_DOC
  integer :: i,j,k

  do k=1,N_det
    i = psi_bilinear_matrix_rows(k)
    j = psi_bilinear_matrix_columns(k)
    psi_det(:,1,k) = psi_det_alpha_unique(:,i)
    psi_det(:,2,k) = psi_det_beta_unique(:,j)
    psi_coef(k,:) = psi_bilinear_matrix_values(k,:)
  end do
  call save_wavefunction_unsorted

end program

