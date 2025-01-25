    Module Block_obs

      Use Matrix

      ! For equal time
      Integer  :: LTROT_ME,  NME_ST, NME_EN 
      COMPLEX  (Kind=8), Dimension(:,:,:,:), Allocatable, Save  ::  DEN
      real (kind=8), save :: kinetic, potential, density, M2
      real (kind=8), save :: S0_11, S0_12, S0_21, S0_22, Sk_11, Sk_12, Sk_21, Sk_22
      REAL (Kind=8), save :: fermicor11_deltaq, fermicor12_deltaq, fermicor21_deltaq, fermicor22_deltaq
      real (kind=8), save :: phaseTot
      Integer, save :: Ncount

      ! For time displaced
      Integer ::  NTDM, NTAUIN

    Contains
      Subroutine Allocate_obs(NLX,NLY,Norb,LTROT)
        Integer :: NLX, NLY, Norb, Ltrot
        Integer :: RX_min, RX_max, RY_min, RY_max
        RX_max = NLX-1
        RX_min = 1-NLX
        RY_max = NLY-1
        RY_min = 1-NLY
        LTROT_ME = 20
        NME_ST = LTROT/2 - LTROT_ME/2 
        NME_EN = LTROT/2 + LTROT_ME/2 
        allocate ( DEN(RX_min:RX_max, RY_min:RY_max, norb, norb) )
      end Subroutine Allocate_obs
      

      Subroutine Allocate_obs_tau(LQ,Norb,LTROT)
        NTAUIN     = (LTROT - NTDM)/2
      end Subroutine Allocate_obs_tau
        
      Subroutine Init_obs
        DEN   = CMPLX( 0.d0 , 0.d0 )
        M2 = 0.0d0
        S0_11 = 0.0d0
        S0_12 = 0.0d0
        S0_21 = 0.0d0
        S0_22 = 0.0d0
        Sk_11 = 0.0d0
        Sk_12 = 0.0d0
        Sk_21 = 0.0d0
        Sk_22 = 0.0d0
        fermicor11_deltaq = 0.d0
        fermicor12_deltaq = 0.d0
        fermicor21_deltaq = 0.d0
        fermicor22_deltaq = 0.d0
        density = 0.0d0
        kinetic = 0.0d0
        potential = 0.d0
        PhaseTot = 0.d0
        Ncount = 0
      End Subroutine Init_obs

      Subroutine Init_obs_tot
      END Subroutine Init_obs_tot

      Subroutine Init_obs_tau
      End Subroutine Init_obs_tau
      
    end Module Block_obs
    


    