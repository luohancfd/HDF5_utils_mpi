program test
    implicit none
    include "mpif.h"
    integer :: n, ncpu, i, j, proc, ierr, master, myid, tag, comm, request(2)
    integer status(MPI_STATUS_SIZE,2)

    integer :: sendbuf, recvbuf, dest

    comm = MPI_COMM_WORLD
    call MPI_Init(ierr)                       ! starts MPI
    call MPI_Comm_rank(comm, myid, ierr)      ! get current proc ID
    call MPI_Comm_size(comm, ncpu, ierr)         ! get number of pro

    dest = myid + 1
    if (dest == ncpu) dest = 0

    sendbuf = myid
    tag = 0
    call MPI_Isend(sendbuf, 1, MPI_INTEGER, dest, tag, comm, request(1), ierr)

    ! write(*,"(A, I2, A, I3)") "Process i=", myid, " send to", dest

    ! call MPI_Irecv(recvbuf, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, request(2), ierr)
    call MPI_Recv(recvbuf, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, status(:,1), ierr)

    ! call MPI_Waitall(2, request, status, ierr)

    call MPI_Wait(request(1), ierr)
    write(*,"(A, I2, A, I3)") "Process i=", myid, " receive data", recvbuf

    call MPI_Finalize(comm, ierr)

end program
