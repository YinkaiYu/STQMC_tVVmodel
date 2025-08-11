       SUBROUTINE INCONFC(ISEED)
         Use blockc
         Implicit Real (KIND=8) (A-G,O-Z)
         Implicit Integer (H-N)
!#define DEC
         INCLUDE 'mpif.h'
         ! Local
         INTEGER STATUS(MPI_STATUS_SIZE)
         INTEGER, Dimension(:,:,:),   Allocatable :: ITMPV1
         INTEGER, Dimension(:,:,:,:), Allocatable :: ITMPV2
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
         Allocate (ITMPV1(LQ,Nfam,LTROT), ITMPV2(LQ,Norb,Nnext,LTROT))
         IF (IRANK .EQ. 0 ) THEN
            OPEN (UNIT=30,FILE='confin', STATUS='UNKNOWN')
            OPEN (UNIT=10,FILE='seeds', STATUS='UNKNOWN')
         ENDIF
         IF ( IRANK.EQ.0 ) THEN
            WRITE(6,*) 'Number of nodes', ISIZE
            READ(30,*) ISEED
            IF (ISEED.EQ.0) THEN
               READ(10,*) ISEED0
               do N = 1,ISIZE - 1
                  !Setup node I and send data.
                  READ(10,*) ITMP
                  do NT = 1,LTROT
                     do I = 1,LQ
                        do nf = 1,Nfam
                           X = RANF(ITMP)
                           NAUX_V1(I,nf,NT) = 1
                           IF (X.GT.0.5) NAUX_V1(I,nf,NT) = -1
                        enddo
                        do no = 1, Norb
                           do nf = 1,Nnext
                              X = RANF(ITMP)
                              NAUX_V2(I,no,nf,NT) = 1
                              IF (X.GT.0.5) NAUX_V2(I,no,nf,NT) = -1
                           enddo
                        enddo
                     enddo
                  enddo
                  call MPI_SEND(ITMP, 1,MPI_INTEGER, N, &
                       & N,MPI_COMM_WORLD,IERR)
                  call MPI_SEND(NAUX_V1, Nfam*LQ*LTROT, MPI_INTEGER, N, &
                       & N+1024, MPI_COMM_WORLD, IERR)
                  call MPI_SEND(NAUX_V2, Nnext*Norb*LQ*LTROT, MPI_INTEGER, N, &
                       & N+2048, MPI_COMM_WORLD, IERR)
               enddo
               ! Set node zero.
               ISEED = ISEED0
               do NT = 1,LTROT
                  do I = 1,LQ
                     do nf = 1, Nfam
                        X = RANF(ISEED)
                        NAUX_V1(I,nf,NT) = 1
                        if (X.GT.0.5) NAUX_V1(I,nf,NT) = -1
                     enddo
                     do no = 1, Norb
                        do nf = 1, Nnext
                           X = RANF(ISEED)
                           NAUX_V2(I,no,nf,NT) = 1
                           if (X.GT.0.5) NAUX_V2(I,no,nf,NT) = -1
                        enddo
                     enddo
                  enddo
               enddo
            else
               ! Read all confins from NODE 0.
               ! Read for Node 0
               do NT = 1,LTROT
                  do i = 1,LQ
                     do nf = 1, NFAM
                        Read(30,*) NAUX_V1(I,nf,NT)
                     enddo
                     do no = 1, Norb
                        do nf = 1, Nnext
                           Read(30,*) NAUX_V2(I,no,nf,NT)
                        enddo
                     enddo
                  enddo
               enddo
               ! Read for Node others
               do N = 1,ISIZE - 1
                  read(30,*) ITMP
                  call MPI_SEND(ITMP, 1,MPI_INTEGER, N, &
                       & N,MPI_COMM_WORLD,IERR)
                  do NT = 1,LTROT
                     do I = 1,LQ
                        do nf = 1,Nfam
                           READ(30,*) ITMPV1(I,nf,NT)
                        enddo
                        do no = 1, Norb
                           do nf = 1, Nnext
                              READ(30,*) ITMPV2(I,no,nf,NT)
                           enddo
                        enddo
                     enddo
                  enddo
                  call MPI_SEND(ITMPV1, Nfam*LQ*LTROT,MPI_INTEGER, N, &
                       & N+1024,MPI_COMM_WORLD,IERR)
                  call MPI_SEND(ITMPV2, Nnext*Norb*LQ*LTROT,MPI_INTEGER, N, &
                       & N+2048,MPI_COMM_WORLD,IERR)
               enddo
            endif
         else
            call MPI_RECV(ISEED, 1,MPI_INTEGER,0, &
                 & IRANK , MPI_COMM_WORLD,STATUS,IERR)
            call MPI_RECV(NAUX_V1, Nfam*LQ*LTROT, MPI_INTEGER,0, &
                 & IRANK + 1024, MPI_COMM_WORLD,STATUS,IERR)
            call MPI_RECV(NAUX_V2, Nnext*Norb*LQ*LTROT, MPI_INTEGER,0, &
                 & IRANK + 2048, MPI_COMM_WORLD,STATUS,IERR)
         endif
         if (IRANK .EQ. 0 ) then
            CLOSE(30); CLOSE(10)
         endif
         Deallocate (ITMPV1, ITMPV2)
         RETURN
       END SUBROUTINE INCONFC
