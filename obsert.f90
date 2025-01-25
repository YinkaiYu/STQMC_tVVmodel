    SUBROUTINE obsert(NT, GRT0_TAU, GR0T_TAU, GRTT_TAU, GR00_TAU)  
        Use Blockc
        Use Block_obs
        Implicit Real (KIND=8) (A-G,O-Z)
        Implicit Integer (H-N)
        complex(kind=8),dimension(:,:) :: GR00_tau, GRTT_tau, GR0T_tau, GRT0_tau
        
   
        
        !G0T_(i,j) =  <c^{dagger}(tau)_j c_i>
        !GT0_(i,j) = <c_i(tau) c^{dagger}_j>

        NT1 = NT + 1
        !	Arguments.
        do LX1 = 1,NLX
            do LX2 = 1,NLX
                do LY1 = 1,NLY
                    do LY2 = 1,NLY
                        lx = npbcx(LX1-LX2)
                        ly = npbcy(LY1-LY2)
                        imj = invlist(lx,ly) 
                        do no1 = 1,norb    
                            do no2 = 1,norb          
                                i = invnlist(lx1,ly1,no1,1)         
                                j = invnlist(lx2,ly2,no2,1)   
                                ! green_tau(imj,no1,no2, NT1) = green_tau(imj,no1,no2, NT1) + &
                                !     &   ( GRT0_TAU(i, j) )* dble(N_SUN) 
                                ! green1_tau(imj,no1,no2, NT1) = green1_tau(imj,no1,no2, NT1) + &       
                                !     &   ( GR0T_TAU(j, i) )* dble(N_SUN) 
                                ! onspair_tau(imj,no1,no2, NT1) = onspair_tau(imj,no1,no2, NT1) + &      
                                !     &   ( GRT0_TAU(i, j) * GRT0_TAU(i, j) )* dble(N_SUN)
                                ! spin_tau(imj,no1,no2, NT1) = spin_tau(imj,no1,no2, NT1) + &      
                                !     &   ( GR0T_TAU(j, i) * GRT0_TAU(i, j) )* dble(N_SUN)             
                                ! spinpm_tau(imj,no1,no2, NT1) = spinpm_tau(imj,no1,no2, NT1) + &            
                                !     &   ( GR0T_TAU(j, i) * GRT0_TAU(i, j) )* dble(N_SUN)          
                                ! den_tau(imj,no1,no2,NT1) = den_tau(imj,no1,no2,nt1)  +  &                      
                                !     &    ((1.d0-GRTT_tau(i,i)-0.5d0)*(1.d0-GR00_tau(j,j)-0.5d0))*dble(N_SUN*N_SUN) &
                                !     &    + (GR0T_TAU(j,i)*GRT0_TAU(i,j))* dble(N_SUN)                              
                            enddo                        
                        enddo                  
                    enddo                   
                enddo
            enddo
        enddo
        
             
      END SUBROUTINE OBSERT
      
