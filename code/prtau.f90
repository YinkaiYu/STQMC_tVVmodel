      SUBROUTINE PRTAU(NOBST)
        Use Blockc
        Use Block_obs
        Implicit Real (KIND=8) (A-G,O-Z)
        Implicit Integer (H-N)
        COMPLEX (KIND =8) :: Znorm
        COMPLEX (KIND=8), Dimension(:,:,:,:), Allocatable :: Collect2
        complex (kind=8), dimension(:,:), allocatable :: collect3
        Character(16) :: filek,filek1
        Interface
           Subroutine Fourier_Trans_tau(gr,filek)
             Complex (Kind=8), dimension(:,:,:,:) :: gr
             Character (16) :: filek
           end Subroutine Fourier_Trans_tau
        end Interface
!#define DEC
        INCLUDE 'mpif.h'
        INTEGER STATUS(MPI_STATUS_SIZE)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
        !write(6,*) ' In Prtau : ', NOBST
        ZNORM = CMPLX(1.D0,0.D0) / DCMPLX( DBLE(NOBST), 0.D0 )
      !   DEN_tau = ZNORM* DEN_tau
      !   spin_tau = znorm* spin_tau
      !   spinpm_tau = znorm * spinpm_tau
      !   GREEN_tau = ZNORM* GREEN_tau
      !   Green1_tau = znorm* Green1_tau
      !   onspair_tau = znorm * onspair_tau

        
      !    Allocate(Collect2(LQ,Norb,Norb,NTDM+1))
      !    allocate(collect3(LQ,NTDM+1))
      !    N = LQ*Norb*Norb*(NTDM+1)
      !    Collect2 = CMPLX(0.D0,0.D0)
      !    CALL MPI_REDUCE(DEN_tau,Collect2,N,MPI_COMPLEX16,MPI_SUM,&
      !         & 0,MPI_COMM_WORLD,IERR)
      !    DEN_tau = Collect2/CMPLX(DBLE(ISIZE),0.D0)
      !    Collect2 = CMPLX(0.D0,0.D0)
      !    CALL MPI_REDUCE(spin_tau,Collect2,N,MPI_COMPLEX16,MPI_SUM,&
      !         & 0,MPI_COMM_WORLD,IERR)
      !    spin_tau = Collect2/CMPLX(DBLE(ISIZE),0.D0)
      !    Collect2 = CMPLX(0.D0,0.D0)
      !    CALL MPI_REDUCE(spinpm_tau,Collect2,N,MPI_COMPLEX16,MPI_SUM,&
      !         & 0,MPI_COMM_WORLD,IERR)
      !    spinpm_tau = Collect2/CMPLX(DBLE(ISIZE),0.D0)
      !    Collect2 = CMPLX(0.D0,0.D0)
      !    CALL MPI_REDUCE(GREEN_tau,Collect2,N,MPI_COMPLEX16,MPI_SUM,&
      !         & 0,MPI_COMM_WORLD,IERR)
      !    GREEN_tau = Collect2/CMPLX(DBLE(ISIZE),0.D0)
      !    Collect2 = CMPLX(0.D0,0.D0)
      !    CALL MPI_REDUCE(GREEN1_tau,Collect2,N,MPI_COMPLEX16,MPI_SUM,&
      !         & 0,MPI_COMM_WORLD,IERR)
      !    GREEN1_tau = Collect2/CMPLX(DBLE(ISIZE),0.D0)
      !    Collect2 = CMPLX(0.D0,0.D0)
      !    CALL MPI_REDUCE(onspair_tau,Collect2,N,MPI_COMPLEX16,MPI_SUM,&
      !         & 0,MPI_COMM_WORLD,IERR)
      !    onspair_tau = Collect2/CMPLX(DBLE(ISIZE),0.D0)

        
      !    IF (IRANK.EQ.0) THEN
      !       filek = "dentau_tot"
      !       Call Fourier_Trans_tau(den_tau ,filek)
      !       filek = "spintau_tot"
      !       Call Fourier_Trans_tau(spin_tau ,filek)
      !       filek = "spinpmtau_tot"
      !       Call Fourier_Trans_tau(spinpm_tau ,filek)
      !       filek = "gtau_tot"
      !       Call Fourier_Trans_tau(green_tau,filek)
      !       filek = "g1tau_tot"
      !       Call Fourier_Trans_tau(green1_tau,filek)
      !       filek = "onspairtau_tot"
      !       Call Fourier_Trans_tau(onspair_tau,filek)
      !    endif

    END SUBROUTINE PRTAU
    
       
    Subroutine Fourier_Trans_tau(gr,filek)
         Use Blockc
         Use Block_obs
         Implicit Real (KIND=8) (A-G,O-Z)
         Implicit Integer (H-N)
         Complex (Kind=8), dimension(:,:,:,:) :: gr
         Integer :: lp
         Character (16) :: filek
         Real (Kind=8) :: xk_p(2), aimj_p(2)
         Complex (Kind=8), allocatable , dimension(:,:,:,:) :: gk
         allocate (gk(LQ,norb,norb,ntdm+1))
         gk = cmplx(0.d0,0.d0)
         do imj = 1,LQ
            aimj_p = dble(list(imj,1))*a1_p + dble(list(imj,2))*a2_p
            do no = 1,norb
               do no1 = 1,norb
                  do nt = 1,ntdm + 1
                     do nk = 1,LQ
                        xk_p = dble(list(nk,1)-1)*b1_p/dble(NLX) + dble(list(nk,2)-1)*b2_p/dble(NLY)
                        gk(nk,no,no1,nt) = gk(nk,no,no1,nt) + &
                             & exp( cmplx( 0.d0, xk_p(1)*aimj_p(1)+xk_p(2)*aimj_p(2)) ) * gr(imj,no,no1,nt)
                     enddo
                  enddo
               enddo
            enddo
         enddo
         gk = gk/cmplx(LQ,0.d0)
         OPEN (UNIT=20,FILE=filek,STATUS='UNKNOWN', action="write", position="append")
            do nk = 1,LQ
               xk_p = dble(list(nk,1)-1)*b1_p/dble(NLX) + dble(list(nk,2)-1)*b2_p/dble(NLY)
               write(20,*) xk_p(1), xk_p(2)
               do nt = 1,ntdm+1
                  do no1 = 1,norb
                     do no2 = 1,norb
                        write(20,*) gk(nk,no1,no2,nt)
                     enddo
                  enddo
               enddo
            enddo
         close(20)
         deallocate (gk)
    end Subroutine Fourier_Trans_tau
    

