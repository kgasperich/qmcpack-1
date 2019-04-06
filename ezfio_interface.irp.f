! DO NOT MODIFY BY HAND
! Created by $QP_ROOT/scripts/ezfio_interface/ei_handler.py
! from file /home/abenali/Work/src/qp2/src/QMCPack/EZFIO.cfg


BEGIN_PROVIDER [ double precision, ci_threshold  ]
  implicit none
  BEGIN_DOC
! Threshold on the CI coefficients as in QMCChem
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_qmcpack_ci_threshold(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: ci_threshold ] <<<<< ..'
      call ezfio_get_qmcpack_ci_threshold(ci_threshold)
    else
      print *, 'qmcpack/ci_threshold not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( ci_threshold, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read ci_threshold with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER
