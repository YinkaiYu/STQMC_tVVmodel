     SUBROUTINE OUTCONFC(ISEED)
        Use Blockc
        Implicit Real (KIND=8) (A-G,O-Z)
        Implicit Integer (H-N)
!#define DEC
        include 'mpif.h'
! Local
        INTEGER STATUS(MPI_STATUS_SIZE)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
        if (IRANK.EQ.0) THEN
            OPEN (UNIT=35,FILE='confout', STATUS='UNKNOWN')
        endif
        if ( IRANK.NE.0 ) THEN
            CALL MPI_SEND(ISEED, 1,MPI_INTEGER, 0, IRANK,MPI_COMM_WORLD,IERR)
            CALL MPI_SEND(NAUX_V1, Nfam*LQ*LTROT,MPI_INTEGER, 0, IRANK+1024,MPI_COMM_WORLD,IERR)
            CALL MPI_SEND(NAUX_V2, Nnext*Norb*LQ*LTROT,MPI_INTEGER, 0, IRANK+2048,MPI_COMM_WORLD,IERR)
        endif
        if (IRANK.EQ.0) THEN
            WRITE(35,*) ISEED
            do NT = 1,LTROT
                do I = 1,LQ
                    do nf = 1,Nfam
                        Write(35,*) NAUX_V1(I,nf,NT)
                    enddo
                    do no = 1,Norb
                        do nf = 1,Nnext
                            Write(35,*) NAUX_V2(I,no,nf,NT)
                        enddo
                    enddo
                enddo
            enddo
            do N = 1,ISIZE - 1
                CALL MPI_RECV(ISEED, 1,MPI_INTEGER,N, N, MPI_COMM_WORLD,STATUS,IERR)
                CALL MPI_RECV(NAUX_V1,Nfam*LQ*LTROT, MPI_INTEGER,N, N+1024, MPI_COMM_WORLD,STATUS,IERR)
                CALL MPI_RECV(NAUX_V2,Nnext*Norb*LQ*LTROT, MPI_INTEGER,N, N+2048, MPI_COMM_WORLD,STATUS,IERR)
                WRITE(35,*) ISEED
                do NT = 1,LTROT
                    do I = 1,LQ
                        do nf = 1,Nfam
                            Write(35,*) NAUX_V1(I,nf,NT)
                        enddo
                        do no = 1,Norb
                            do nf = 1,Nnext
                                Write(35,*) NAUX_V2(I,no,nf,NT)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        endif
        if (IRANK.EQ.0) THEN
            CLOSE(35)
        endif  
    END SUBROUTINE OUTCONFC
