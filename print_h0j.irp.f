program print_h0j
  use bitmasks
  implicit none
  BEGIN_DOC
! Print 1st line of H matrix
  END_DOC
  integer :: i,j,k
  integer(bit_kind) :: det2(N_int,2)
  double precision :: hij, E

  print *,  '< 0 | H | j >'
  print *,  '============='
  print *,  ''

  E = 0.d0
  do k=1,N_det
    i = psi_bilinear_matrix_rows(k)
    j = psi_bilinear_matrix_columns(k)
    det2(:,1) = psi_det_alpha_unique(:,i)
    det2(:,2) = psi_det_beta_unique(:,j)

    call i_h_j(psi_det(1,1,1), det2(1,1), N_int,hij)
    print *,  k, psi_bilinear_matrix_values(k,1), hij
    E += psi_bilinear_matrix_values(k,1)*hij
  end do
  E = E/psi_bilinear_matrix_values(1,1)
  print *,  'nuclear_repulsion = ', nuclear_repulsion
  print *,  'E = ', E + nuclear_repulsion
end program

