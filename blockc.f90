      Module  Blockc 
        
!!        Use Lattices_v3


        Real (Kind=8),save :: BETA , RJ, RHUB,  RT1, DTAU, PI,  TwistX
        Integer,      save :: LTROT, NWRAP, N_SUN, NBIN, NSWEEP,  NORB, NORB1, NDIM, LTAU,  &
             &                LFAM, NFAM, Nbond, nspin, LQ, NLX, NLY, NE, Itwist
        real(kind=8), DIMENSION(:), save :: a1_p(2), a2_p(2), b1_p(2), b2_p(2)

        Integer,  Dimension(:,:), Allocatable, Save :: list(:,:), invlist(:,:), nlist(:,:), invnlist(:,:,:,:)
        Integer,  dimension(:,:), allocatable, save :: lattimj(:,:)
	    INTEGER,  Dimension(:,:), Allocatable, Save :: L_bonds, NSIGL_U
        INTEGER,  Dimension(:,:,:), Allocatable, Save :: NSIGL_K
        INTEGER,  Save ::  NFLIPL(-1:1)
        
        COMPLEX (Kind=8), Dimension(:,:),   Allocatable, Save ::  PROJ, ZKRON
        COMPLEX (Kind=8), Dimension(:,:)  , Allocatable, Save ::  URT_tot, URTM1_tot

        COMPLEX (Kind=8),  Save ::   XSIGMA_U(-1:1), DELTA_U(-1:1), &
             &                       UR_K(2,2), URT_K(2,2), obs(30)

        REAL (Kind=8), Save ::   XSIGP2(-1:1),XSIGM2(-1:1), &
             &                   DELLP2(-1:1), DELLM2(-1:1), ETAL(-2:2), &
             &                   FD(4), ZERO

        Logical :: L_Trot_hop

        Contains

        Subroutine Allocate_Blockc

          Implicit none

          Zero    =  1.0D-10
          PI      =  acos(-1.d0)

          a1_p(1) = 1.d0
          a1_p(2) = 0.d0
          a2_p(1) = 0.5d0
          a2_p(2) = 0.5d0*sqrt(3.0d0)
          b1_p(1) = 2.d0*Pi
          b1_p(2) = -2.d0*Pi/sqrt(3.d0)
          b2_p(1) = 0.d0
          b2_p(2) = 2.d0*Pi*2.d0/sqrt(3.d0)
          
          NORB = 2 
          NORB1 = 2
          LQ   = NLX * NLY  
          NDIM = Norb*LQ
          NFAM   = 3  ! * 2 for current and for kenetic terms.
          LFAM   = LQ
          Nbond = NFam
          nspin = 1
          DTAU   = BETA/DBLE(LTROT)
          
          Allocate( L_Bonds(LQ,0:Nbond), list(LQ,2), invlist(NLX,NLY), nlist(Ndim, 4), invnlist(NLX,NLY,norb,nspin), &
               &    NSIGL_K(LQ,Nfam,LTROT), NSIGL_U(NDIM,LTROT), Lattimj(LQ,LQ)  )
          Allocate ( PROJ(NDIM,NDIM),   ZKRON(NDIM,NDIM) ) 
          Allocate ( URT_tot(NDIM,NDIM), URTM1_tot(NDIM,NDIM) )

        end Subroutine Allocate_Blockc

      end Module Blockc
