SUBROUTINE SPROJ(DEGEN,EN_FREE)
   !Sets projector.
   Use Blockc
   Use MyMats
   Implicit Real (KIND=8) (A-G,O-Z)
   Implicit Integer (H-N)
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
   if (Itwist == 1) then ! DSM
      Call SetH(TMP)
   endif
   if (Itwist == 2) then ! DSM
      Call SetH(TMP)
   endif
   if (Itwist == 0) then ! CDW
      IseedHop = 3958195
      do ii = 1, Ndim
         random = ranf(IseedHop)
         no = Nlist(ii,3)
         R1 = (-RT1)**dble(no) * ( 1.0d0 + TwistX*(random-0.5d0) ) 
         TMP(ii,ii) = CMPLX( R1, 0.d0)
      enddo
   endif
   if (Itwist == -1) then ! Random
      IseedHop = 3958195
      do ii = 1, Ndim
         random = ranf(IseedHop)
         R1 = RT1 * (random-0.5d0)
         TMP(ii,ii) = CMPLX( R1, 0.d0)
      enddo
   endif
   if (Itwist == -2) then ! Random and cdw=0
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
   if (Itwist == -3) then ! CM
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
   CALL Diag(TMP,PROJ,WC)
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
