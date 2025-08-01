      Module  Blockc 
        
!!        Use Lattices_v3
          interface
               function ranf(iseed) result(r)
                    integer, intent(in) :: iseed
                    real(8) :: r
               end function ranf
          end interface

        Real (Kind=8),save :: BETA , RV1, RV2, RT1, DTAU, PI,  TwistX, iniQAHt2
        Integer,      save :: LTROT, NWRAP, N_SUN, NBIN, NSWEEP,  NORB, NORB1, NDIM, LTAU,  &
             &                LFAM, NFAM, Nbond, Nnext, nspin, LQ, NLX, NLY, NE, Itwist
        real(kind=8), DIMENSION(:), save :: a1_p(2), a2_p(2), b1_p(2), b2_p(2)

        Integer,  Dimension(:,:), Allocatable, Save :: list(:,:), invlist(:,:), nlist(:,:), invnlist(:,:,:,:)
        Integer,  dimension(:,:), allocatable, save :: lattimj(:,:)
	     INTEGER,  Dimension(:,:), Allocatable, Save :: L_bonds
	     INTEGER,  Dimension(:,:,:), Allocatable, Save :: L_next
        INTEGER,  Dimension(:,:,:), Allocatable, Save :: NAUX_V1
        INTEGER,  Dimension(:,:,:,:), Allocatable, Save :: NAUX_V2
        INTEGER,  Save ::  NFLIPL(-1:1)
        
        COMPLEX (Kind=8), Dimension(:,:),   Allocatable, Save ::  PROJ, ZKRON
        COMPLEX (Kind=8), Dimension(:,:)  , Allocatable, Save ::  URT_tot, URTM1_tot

        COMPLEX (Kind=8),  Save ::   XSIGMA_U(-1:1), DELTA_U(-1:1), &
             &                       UR_K(2,2), URT_K(2,2), obs(30)

        REAL (Kind=8), Save ::   XSIGP1(-1:1),XSIGM1(-1:1),XSIGP2(-1:1),XSIGM2(-1:1), &
             &                   DELLP1(-1:1),DELLM1(-1:1),DELLP2(-1:1),DELLM2(-1:1), ETAL(-2:2), &
             &                   FD(4), ZERO, ratio_const(-1:1)

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
          NFAM   = 3 ! nearest neighbors for honeycomb: 3 A->B
          LFAM   = LQ
          Nbond = NFam
          Nnext = 3  ! next nearest neighbors for honeycomb: 3 A->A or B->B
          nspin = 1
          DTAU   = BETA/DBLE(LTROT)
          
          Allocate( L_Bonds(LQ,0:Nbond), L_next(LQ,Norb,0:Nnext), list(LQ,2), invlist(NLX,NLY), nlist(Ndim, 4), invnlist(NLX,NLY,norb,nspin), &
               &    NAUX_V1(LQ,Nfam,LTROT), NAUX_V2(LQ,Norb,Nnext,LTROT), Lattimj(LQ,LQ)  )
          Allocate ( PROJ(NDIM,NDIM),   ZKRON(NDIM,NDIM) ) 
          Allocate ( URT_tot(NDIM,NDIM), URTM1_tot(NDIM,NDIM) )

        end Subroutine Allocate_Blockc

      end Module Blockc
