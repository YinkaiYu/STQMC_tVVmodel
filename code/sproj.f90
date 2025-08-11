SUBROUTINE SPROJ(DEGEN,EN_FREE)
   !Sets projector.
   Use Blockc
   Use MyMats
   Implicit Real (KIND=8) (A-G,O-Z)
   Implicit Integer (H-N)
   Complex (Kind=8) :: Z1
!#define DEC
   Interface
      Subroutine SetH(HLP2)
         Complex (Kind=8), Dimension(:,:) :: HLP2
      end Subroutine SetH
   end Interface
   INCLUDE 'mpif.h' 
   INTEGER STATUS(MPI_STATUS_SIZE)
   COMPLEX (Kind=8), Dimension(:,:), Allocatable :: TMP
   real (Kind=8), Dimension(:), Allocatable :: WC
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
   Allocate(TMP(Ndim,Ndim), WC(Ndim))
   PROJ = CMPLX(0.d0,0.d0)
   TMP = CMPLX(0.d0,0.d0)
   write(6,*) 'Itwist=', Itwist
   if (Itwist == 1) then ! DSM
      write(6,*) 'Initial state: DSM (Itwist=1)'
      Call SetH(TMP)
   endif
   if (Itwist == 2) then ! DSM
      write(6,*) 'Initial state: DSM (Itwist=2)'
      Call SetH(TMP)
   endif
   if (Itwist == 3) then ! DSM (QAH twist)
      write(6,*) 'Initial state: DSM (QAH twist)'
      Call SetH(TMP)
   endif
   if (Itwist == 0) then ! CDW
      write(6,*) 'Initial state: CDW (Itwist=0)'
      IseedHop = 3958195
      do ii = 1, Ndim
         random = ranf(IseedHop)
         no = Nlist(ii,3)
         R1 = (-RT1)**dble(no) * ( 1.0d0 + TwistX*(random-0.5d0) ) 
         TMP(ii,ii) = CMPLX( R1, 0.d0)
      enddo
   endif
   if (Itwist == -1) then ! Random
      write(6,*) 'Initial state: Random (Itwist=-1)'
      IseedHop = 3958195
      do ii = 1, Ndim
         random = ranf(IseedHop)
         R1 = RT1 * (random-0.5d0)
         TMP(ii,ii) = CMPLX( R1, 0.d0)
      enddo
   endif
   if (Itwist == -2) then ! Random and cdw=0
      write(6,*) 'Initial state: Random and CDW=0 (Itwist=-2)'
      IseedHop = 3958195
      do i = 1, LQ
         random = ranf(IseedHop)
         R1 = RT1 * (random-0.5d0)
         ix = list(i,1)
         iy = list(i,2)
         ii_A = invNlist(ix,iy,1,1)
         ii_B = invNlist(ix,iy,2,1)
         TMP(ii_A,ii_A) = CMPLX( R1, 0.d0)
         TMP(ii_B,ii_B) = CMPLX( -R1, 0.d0)
      enddo
   endif
   if (Itwist == -3) then ! CM cos wave
      write(6,*) 'Initial state: CM cos wave (Itwist=-3)'
      do ix = 1, NLX
         do iy = 1, NLY
            do no = 1, Norb
               ii = invnlist(ix,iy,no,1)
               R1 = -0.25d0 * (-1.0d0)**dble(no)
               if (MOD(ix-iy,3)==0) R1 = -2.0d0 * R1
               TMP(ii,ii) = CMPLX( R1, 0.d0)
            enddo
         enddo
      enddo
   endif
   if (Itwist == -4) then ! QAH
      write(6,*) 'Initial state: QAH (Itwist=-4)'
      ! t1
      do I = 1,LQ
         ii_0 = L_Bonds(i,0)
         do nf = 1,Nbond
            ii_n = L_Bonds(i,nf)
            Z1 = CMPLX(-RT1,0.D0)
            TMP(ii_0,ii_n)  =  Z1
            TMP(ii_n,ii_0)  = conjg(Z1)
         enddo
      enddo
      ! t2
      do ix = 1, NLX
         do iy = 1, NLY
            do i_sublattice = 1, Norb
               do i_direction = 1, Nnext
                  i = invlist(ix,iy)
                  ii_0 = L_next(i,i_sublattice,0)
                  ii_n = L_next(i,i_sublattice,i_direction)
                  R2 = - iniQAHt2 * (-1.0d0)**dble(i_sublattice+i_direction)
                  TMP(ii_0,ii_n) = CMPLX( 0.d0,  R2)
                  TMP(ii_n,ii_0) = CMPLX( 0.d0, -R2)
               enddo
            enddo
         enddo
      enddo
   endif
   if (Itwist == -5) then ! CM (1308.6094)
      write(6,*) 'Initial state: CM (1308.6094)'
      do ix = 1, NLX
         do iy = 1, NLY
            do no = 1, Norb
               ii = invnlist(ix,iy,no,1)
               R1 = -0.5d0 * (-1.0d0)**dble(no)
               if (MOD(ix-iy,3)==0) R1 = - R1
               TMP(ii,ii) = CMPLX( R1, 0.d0)
            enddo
         enddo
      enddo
   endif
   CALL Diag(TMP,PROJ,WC)
   ! write(6,*) 'matrix of TMP'
   ! WRITE(6, '(A7)',advance='no') "Col/Row" 
   ! DO j = 1, Ndim
   !    WRITE(6,'(7X,I5,8X)', advance='no') j
   ! ENDDO
   ! WRITE(6,*) ! \n
   ! DO i = 1, Ndim
   !    WRITE(6, '(I5,2X)',advance='no') i 
   !    DO j = 1, Ndim
   !       WRITE(6,'( "(", F8.4, ",", F8.4, ") " )', advance='no') REAL(TMP(i,j)), AIMAG(TMP(i,j))
   !    ENDDO
   !    WRITE(6,*) ! \n
   ! ENDDO
   ! write(6,*) 'matrix of PROJ'
   ! WRITE(6, '(A7)',advance='no') "Col/Row" 
   ! DO j = 1, Ndim
   !    WRITE(6,'(7X,I5,8X)', advance='no') j
   ! ENDDO
   ! WRITE(6,*) ! \n
   ! DO i = 1, Ndim
   !    WRITE(6, '(I5,2X)',advance='no') i 
   !    DO j = 1, Ndim
   !       WRITE(6,'( "(", F8.4, ",", F8.4, ") " )', advance='no') REAL(PROJ(i,j)), AIMAG(PROJ(i,j))
   !    ENDDO
   !    WRITE(6,*) ! \n
   ! ENDDO
   en_free = 0.d0
   do i = 1,Ne
      en_free = en_free + wc(i)
   enddo
   en_free = en_free*dble(N_SUN)
   DEGEN = WC(NE+1) - WC(NE)
   IF (IRANK == 0) WRITE(50,*) 'Degen: ', DEGEN

   Deallocate(TMP)
   Deallocate(WC)
   RETURN
END SUBROUTINE SPROJ
