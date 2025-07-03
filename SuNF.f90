   Program SuNF
   Use Blockc
   Use Block_obs
! Use Block_tau
   Use MyMats
   Implicit Real (KIND=8) (A-G,O-Z)
   Implicit Integer (H-N)
! T=0 SU(N) t-U-J model. Flux phases t-U-J model. With large-N, Honeycomb lattice.
! L is the linear size of square lattice.
! LQ is the number of unit cells.
! Norb number of orbitals in unit cell.
! Ndim = Norb*LQ
! Trotter break-up. There are 6 families.
! e^(J_3) e^(J_2) e^(J_1) e^(K_3) e^(K_2) e^(K_1) e^(T)
! Since observables are rather expensive to compute we will calculate
! them on a maximum of LTROT_ME time slices.
   Interface
      SUBROUTINE MMUUR(A, NF, NT ,NFLAG)
        COMPLEX (Kind=8), Dimension(:,:) :: A
        INTEGER :: NF, NT, NFLAG
      END SUBROUTINE MMUUR
      SUBROUTINE MMUURM1(A, NF, NTAU, NFLAG)
        COMPLEX (KIND=8), Dimension(:,:) :: A
        Integer :: NF, NTAU, NFLAG
      END SUBROUTINE MMUURM1
      SUBROUTINE MMUUL(A, NF,NTAU,NFLAG)
        COMPLEX (Kind=8), Dimension(:,:) :: A
        Integer :: NF,NTAU,NFLAG
      END SUBROUTINE MMUUL
      SUBROUTINE MMUULM1(A, NF, NTAU, NFLAG)
        COMPLEX (KIND=8), Dimension(:,:) :: A
        Integer :: NF, NTAU, NFLAG
      END SUBROUTINE MMUULM1
      SUBROUTINE ORTHO(A,I)
        COMPLEX (Kind=8), Dimension(:,:) :: A
        INTEGER :: I
      END SUBROUTINE ORTHO
      SUBROUTINE obser(UL, UR, ULRINV,phase)     
        COMPLEX (Kind=8), Dimension(:,:) :: UL, UR, ULRINV
        complex(kind=8) :: phase
      END SUBROUTINE obser
      SUBROUTINE upgradeV1(NTAU,NF,ISEED,UL,UR,ULRINV,phase,NFLAG)
        COMPLEX (Kind=8), Dimension(:,:) :: UL, UR, ULRINV
        complex (kind=8) :: phase
        Integer :: NTAU,NF,ISEED,NFLAG
      END SUBROUTINE upgradeV1
      SUBROUTINE upgradeV2(NTAU,NF,ISEED,UL,UR,ULRINV,phase,NFLAG)
        COMPLEX (Kind=8), Dimension(:,:) :: UL, UR, ULRINV
        complex (kind=8) :: phase
        Integer :: NTAU,NF,ISEED,NFLAG
      END SUBROUTINE upgradeV2
      SUBROUTINE MMTHR(A)
        COMPLEX (Kind=8), Dimension(:,:) :: A
      END SUBROUTINE MMTHR
      SUBROUTINE MMTHRM1(A)
        COMPLEX (Kind=8), Dimension(:,:) :: A
      END SUBROUTINE MMTHRM1
      SUBROUTINE MMTHL(A)
        COMPLEX (Kind=8), Dimension(:,:) :: A
      END SUBROUTINE MMTHL
      SUBROUTINE MMTHLM1(A)
        COMPLEX (Kind=8), Dimension(:,:) :: A
      END SUBROUTINE MMTHLM1
      SUBROUTINE DYN( UST, UL, UR, ULR,ULRINV, XMEAN_DYN, XMAX_DYN)
         Complex (Kind=8), Dimension(:,:,:) :: UST
         Complex (Kind=8), Dimension(:,:) :: UL, UR, ULR, ULRINV
         Real (Kind=8) :: XMEAN_DYN, XMAX_DYN
      END SUBROUTINE DYN
   END Interface
!#define DEC
   INCLUDE 'mpif.h'
   INTEGER STATUS(MPI_STATUS_SIZE)
! Space for storage.
   COMPLEX (KIND=8), Dimension(:,:), Allocatable :: UL ,UR , ULR, ULRINV
   COMPLEX (KIND=8), Dimension(:,:), Allocatable :: UL_1,UR_1, ULRINV_1,ULR1,ULRINV1
   COMPLEX (KIND=8), Dimension(:,:,:),Allocatable :: UST
   COMPLEX (KIND=8) :: DETZ, ALPHA_U, phase
   COMPLEX (Kind=8) :: DET1(2)
   NCON = 0
   CALL MPI_INIT(IERR)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
   IF (IRANK == 0) THEN
      OPEN (UNIT=50,FILE='info',STATUS='UNKNOWN')
      OPEN ( UNIT=20, FILE='paramC_sets',STATUS='UNKNOWN' )
      READ(20,*) BETA, LTROT, NWRAP, RT1, RV1, RV2
      READ(20,*) NBIN, NSWEEP, LTAU, NTDM
      READ(20,*) NLX, NLY, Itwist, TwistX, N_SUN, NE
      CLOSE(20)
   ENDIF
      L_Trot_hop = .false.
      NB_Field = 0
   CALL MPI_BCAST(BETA ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST(RV1 ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST(RV2 ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST(RT1 ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST(Itwist ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST(TwistX ,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST(NLX ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST(NLY ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST(NE ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST(LTROT ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST(NWRAP ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST(N_SUN ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST(NBIN ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST(NSWEEP ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST(LTAU ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST(NTDM ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL Allocate_Blockc
   CALL Allocate_obs( NLX, NLY, Norb, LTROT )
   IF (LTAU == 1) THEN
      CALL Allocate_obs_tau(LQ,Norb,LTROT)
      ALLOCATE( UL_1(NE,Ndim), UR_1(Ndim,NE), ULRINV_1(NE,NE) )
   ENDIF
   NST = LTROT/NWRAP
   ALLOCATE(UL (NE,NDIM),UR(NDIM,NE), ULR(NE,NE), ULRINV(NE,NE), UST(NDIM,NE,NST))
   ALLOCATE(ULR1(NE,NE),ULRINV1(NE,NE))
! Write parameters.
   IF (IRANK == 0) THEN
      WRITE(50,*) '========================='
      WRITE(50,*) 'SU(N) t-U-J, Projector '
      WRITE(50,*) 'Linear lenght    :',NLX,NLY
      WRITE(50,*) 'N: SU(N)         :',N_SUN
      WRITE(50,*) 'Particle #       :',NE
      WRITE(50,*) 'Hopping t        :',RT1
      WRITE(50,*) 'V1 < 0           :',RV1
      WRITE(50,*) 'V2 > 0           :',RV2
      WRITE(50,*) 'Theta            :',BETA
      WRITE(50,*) 'Trotter number   :',LTROT
      WRITE(50,*) '=>Dtau           :',DTAU
      WRITE(50,*) 'N_Ortho          :',NWRAP
      WRITE(50,*) '# Bins           :',NBIN
      WRITE(50,*) '# Sweeps/Bin     :',NSWEEP
      WRITE(50,*) 'Measuring on the ',LTROT_ME,' middle time slices'
      IF (LTAU == 1) THEN
         WRITE(50,*) 'TAU_max            :', DBLE(NTDM)*DTAU
         WRITE(50,*) 'Effective proj.    :', DBLE(NTAUIN)*DTAU
      ENDIF
      WRITE(50,*) '# of families    :', NFAM
      WRITE(50,*) 'Length of fam.   :', LFAM
      WRITE(50,*) 'Itwist           :', Itwist
      WRITE(50,*) 'Twist in x-direction phi/phi_0 : ', TwistX
      if ( L_Trot_hop ) then
         Write(50,*) ' Trotter decomposition for hopping'
      else
         Write(50,*) ' No Trotter decomposition for hopping'
      endif
      WRITE(50,*) '# Nodes         :',ISIZE
   ENDIF
! Check parameters.
   IF (LTAU == 1) THEN
      IF (MOD(NTAUIN,NWRAP).NE.0) THEN
         WRITE(50,*) 'NTAUIN is not a multiple of NOFRE'
         STOP
      ENDIF
      IF (MOD(NTDM,NWRAP).NE.0) THEN
         WRITE(50,*) 'NTDM is not a multiple of NOFRE'
         STOP
      ENDIF
   ENDIF
! Setup lists and hopping matrix.
   CALL SLI ! Setup lattice
   CALL SALPH(Alpha_U) ! Calculate DTAU and ALPHA - precalc factors for updates & measurements
   CALL SPROJ(DEGEN,EN_FREE) ! Setup Projector PROJ and TWF, compute charge gap/degeneracy
   CALL STHOP ! Setup & exponentiate hopping matrix
   IF (IRANK == 0) THEN
      WRITE(50,*) 'Trial: Degen, En:  ', DEGEN, EN_FREE
   ENDIF
   CALL INCONFC(ISEED) ! Setup random sigma configuration or read it from unit 10
   CALL SYSTEM_CLOCK(COUNT_RATE=N_P_SEC)
   CALL SYSTEM_CLOCK(COUNT=ICPU_1)
!-------Fill up right storage and set the phase.
   DO NL = 1,NE
      DO I = 1,Ndim
         UR(I,NL) = PROJ(I,NL)
      ENDDO
   ENDDO
   if ( LTROT > 0 ) then
      DO NT = 1,LTROT
         CALL MMTHR(UR)
         IF (abs(RV1) > Zero) THEN
            DO NF = 1,NFAM
               NFLAG = 1
               CALL MMUUR(UR, NF, NT, NFLAG)
            ENDDO
         ENDIF
         IF (abs(RV2) > Zero) THEN
            DO NF = 1,Nnext
               NFLAG = 2
               CALL MMUUR(UR, NF, NT, NFLAG)
            ENDDO
         ENDIF
         IF (MOD(NT,NWRAP) == 0) THEN
            CALL ORTHO(UR,NCON)
            ! Write in Store (NT,NF+1) with ortho.
            NT_ST = NT/NWRAP
            DO NL = 1,NE
               DO NR = 1,Ndim
                  UST(NR,NL,NT_ST) = UR(NR,NL)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
   end if
! UR is now on time slice LTROT.
! Set UL
   DO NL = 1,NE
      DO I = 1,Ndim
         UL(NL,I) = DCONJG(PROJ(I,NL))
      ENDDO
   ENDDO
! Set phase.
   CALL MMULT (ULR,UL,UR)
   CALL INV (ULR,ULRINV,DET1)
   phase = det1(1)/abs(det1(1))

! Test is OK.
!************** Start Sweeps. *********************************
   XMAXP = 0.D0
   ACC_JK = 0.D0
   ACC_U = 0.D0
   AXMAX_DYN = 0.D0
   XMEAN_DYN= 0.D0
   DO NBC = 1, NBIN
   ! Here, you have the green functions on time slice 1.
   ! Set bin observables to zero.
      CALL Init_Obs
      If (Ltau == 1) Call Init_obs_tau
      OBS = CMPLX(0.d0,0.d0)
      NOBS = 0
      NOBST = 0
      SIGNT = 0.D0
      DO NSW = 1, NSWEEP
      ! Set UL, UST is full with right propagations.
         DO NL = 1,NE
            DO I = 1,Ndim
               UL(NL,I) = DCONJG(PROJ(I,NL))
            ENDDO
         ENDDO
         CALL ORTHO(UR,NCON)
         CALL MMULT(ULR,UL,UR)
         CALL INV (ULR,ULRINV,DET1)
         phase = det1(1)/abs(det1(1))
         if ( LTROT > 0 ) then
            DO NT = LTROT, 1, -1
            ! UR, UL on time slice Ltrot.
               IF (MOD(NT,NWRAP) == 0 ) THEN
                  ! Read UR.
                  NT_ST = NT/NWRAP
                  DO NL = 1,NE
                     DO I = 1,Ndim
                        UR(I,NL) = UST(I,NL,NT_ST)
                     ENDDO
                  ENDDO
                  ! Ortho of UL and recalc of ulrinv
                  IF (NT.NE.LTROT) CALL ORTHO(UL,NCON)
                  CALL MMULT(ULR,UL,UR)
                  CALL INV (ULR,ULRINV,DET1)
                  phase = det1(1)/abs(det1(1))
                  ! Store UL.
                  NT_ST = NT/NWRAP
                  DO NL = 1,NE
                     DO I = 1,Ndim
                        UST(I,NL,NT_ST) = UL(NL,I)
                     ENDDO
                  ENDDO
               ENDIF
               ! Compute observables
               ! IF (NT >= NME_ST .AND. NT <= NME_EN) THEN
               IF ( NT == LTROT/2 ) THEN
                  Call OBSER(UL,UR,ULRINV,phase)
                  NOBS = NOBS + 1
               ENDIF
               ! Propagate UL and UR by multiplying U-matrix
               if (abs(RV2) > Zero) THEN
                  DO NF = Nnext,1,-1
                     NFLAG = 2
                     CALL upgradeV2(NT,NF,ISEED,UL,UR,ULRINV,phase,NFLAG)
                     CALL MMUUL (UL,NF,NT,NFLAG)
                     CALL MMUURM1(UR,NF,NT,NFLAG)
                  ENDDO
               endif
               if (abs(RV1) > Zero) THEN
                  DO NF = NFAM,1,-1
                     NFLAG = 1
                     CALL upgradeV1(NT,NF,ISEED,UL,UR,ULRINV,phase,NFLAG)
                     CALL MMUUL (UL,NF,NT,NFLAG)
                     CALL MMUURM1(UR,NF,NT,NFLAG)
                  ENDDO
               endif
         
               CALL MMTHL (UL)
               CALL MMTHRM1(UR)
            ENDDO
         else
         ! LTROT == 0:  initial state, directly observe
            Call OBSER(UL,UR,ULRINV,phase)
            NOBS = NOBS + 1
            phasetot = phasetot + real(phase)
            ncount = ncount + 1
         end if

      ! Set UR, UST is full with left propagations. UL is on time slice 0.
         UR(1:Ndim,1:NE) = PROJ(1:Ndim,1:NE)
         CALL ORTHO(UL,NCON)
         CALL MMULT(ULR,UL,UR)
         CALL INV (ULR,ULRINV,DET1)
         phase = det1(1)/abs(det1(1))
         if ( LTROT > 0 ) then
            DO NT = 1,LTROT
            ! You start the loop on time slice NT-1
               CALL MMTHR (UR)
               CALL MMTHLM1(UL)
         
               IF (abs(RV1) > Zero) THEN
                  DO NF = 1,NFAM
                     NFLAG = 1
                     CALL MMUUR (UR, NF, NT, NFLAG)
                     CALL MMUULM1 (UL, NF, NT, NFLAG)
                     CALL upgradeV1(NT,NF,ISEED,UL,UR,ULRINV,phase,NFLAG)
                  ENDDO
               ENDIF
               IF (abs(RV2) > Zero) THEN
                  DO NF = 1,Nnext
                     NFLAG = 2
                     CALL MMUUR (UR, NF, NT, NFLAG)
                     CALL MMUULM1 (UL, NF, NT, NFLAG)
                     CALL upgradeV2(NT,NF,ISEED,UL,UR,ULRINV,phase,NFLAG)
                  ENDDO
               ENDIF
               ! IF (NT >= NME_ST .AND. NT <= NME_EN) THEN
               IF ( NT == LTROT/2 ) THEN
                  Call OBSER(UL,UR,ULRINV,phase)
                  NOBS = NOBS + 1
               ENDIF
               IF (MOD(NT,NWRAP) == 0) THEN
                  NT_ST = NT/NWRAP
                  ! Read.
                  DO NL = 1,NE
                     DO NR = 1,Ndim
                        UL(NL,NR) = UST(NR,NL,NT_ST)
                     ENDDO
                  ENDDO
                  CALL ORTHO(UR,NCON)
                  CALL MMULT(ULR,UL,UR)
                  CALL INV (ULR,ULRINV,DET1)
                  phsae = DET1(1)/abs(DET1(1))
                  ! Store Here
                  DO NL = 1,NE
                     DO NR = 1,Ndim
                        UST(NR,NL,NT_ST) = UR(NR,NL)
                     ENDDO
                  ENDDO
               ENDIF
               IF (LTAU == 1 .AND. NT == NTAUIN ) THEN
                  UL_1 = UL
                  UR_1 = UR
                  ULRINV_1 = ULRINV
                  CALL DYN(UST,  UL_1, UR_1, ULR,ULRINV_1, XMEAN_DYN, XMAX_DYN)
                  NOBST = NOBST + 1
               ENDIF
            ENDDO
         else
         ! LTROT == 0:  initial state, directly observe
            CALL ORTHO(UR,NCON)
            CALL MMULT(ULR,UL,UR)
            CALL INV (ULR,ULRINV,DET1)
            phase = det1(1)/abs(det1(1))      
            Call OBSER(UL,UR,ULRINV,phase)
            NOBS = NOBS + 1
            phasetot = phasetot + real(phase)
            ncount = ncount + 1
         end if
      ENDDO
      CALL PREQ (NOBS,NBC)
      IF (LTAU == 1) THEN
         CALL PRTAU(NOBST)
      ENDIF
      IF (DBLE(OBS(30)) > ZERO) ACC_JK = ACC_JK + DBLE(OBS(29)/OBS(30))
      IF (DBLE(OBS(28)) > ZERO) ACC_U = ACC_U + DBLE(OBS(27)/OBS(28))
   ENDDO
   CALL OUTCONFC(ISEED)

!!$ IF (LTAU == 1) XMEAN_DYN = XMEAN_DYN/(DBLE(NOBST)*DBLE(NTDM/NWRAP))
   XMAX_DYN1 = 0.0
   CALL MPI_REDUCE(XMAX_DYN,XMAX_DYN1,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,IERR)
   XMAX_DYN = XMAX_DYN1
   XMEAN_DYN1 = 0.D0
   CALL MPI_REDUCE(XMEAN_DYN,XMEAN_DYN1,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,IERR)
   XMEAN_DYN = XMEAN_DYN1/DBLE(ISIZE)
   CALL MPI_REDUCE(XMAXP,XMAXP_1,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,IERR)
   XMAXP = XMAXP_1
   CALL SYSTEM_CLOCK(COUNT=ICPU_2)
   CPUT = 0.D0
   CPUT = DBLE(ICPU_2 - ICPU_1)/DBLE(N_P_SEC)
   CPUT_1 = 0.D0
   CALL MPI_REDUCE(CPUT,CPUT_1,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,IERR)
   CPUT = CPUT_1/DBLE(ISIZE)
   IF (IRANK == 0) THEN
      WRITE(50,*) ' Max  diff Ph :', XMAXP
      IF (LTAU == 1) WRITE(50,*) ' Dynamics Mean,Max: ', XMEAN_DYN, XMAX_DYN
      WRITE(50,*) ' Accep_J    :',ACC_JK/DBLE(NBIN)
      WRITE(50,*) ' Accep_U    :',ACC_U /DBLE(NBIN)
      WRITE(50,*) 'Tot CPU time:', CPUT
   ENDIF
   CALL MPI_FINALIZE(IERR)
    END PROGRAM SuNF
