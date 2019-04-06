program print_pgm
 implicit none

 integer :: i,j
 character*(64) :: fmt

 print '(A)', '-------------------------------'
 print '(A)', 'P2'
 print '(A)', '#'
 print *, n_det_alpha_unique, n_det_beta_unique
 print *, 255

 write(fmt,*) '(',n_det_beta_unique,'(I3,X))'

 integer, external :: f
 do i=1,n_det_alpha_unique
    write(*,fmt)  (f(psi_bilinear_matrix(i,j,1)), j=1,n_det_beta_unique)
 enddo
end program

integer function f(x)
  implicit none
  double precision :: x, df
  df = 255.d0*erf(abs(10.d0*x)**(0.25d0))
  f = int(df)
end

