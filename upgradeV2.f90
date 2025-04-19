    SUBROUTINE upgradeV(NTAU,NF,ISEED,UL,UR,ULRINV,phase,NFLAG)     
    
    Use Blockc
    Use Block_obs
        
    Implicit Real (KIND=8) (A-G,O-Z)
        
    Implicit Integer (H-N)
        
        
        ! Space for tests


        !Arguments

	COMPLEX (Kind=8), Dimension(:,:) :: UL, UR, ULRINV
    complex(kind=8) :: phase    
    Integer :: NTAU,NF,ISEED,NFLAG
        
        !Local
	COMPLEX (Kind=8), Dimension(:), Allocatable :: VEC1, VEC2, &
             &    VHLP1, UHLP1, VHLP2, UHLP2,  V1, V2, U1, U2
	
	COMPLEX (Kind=8) ::  G44UP, G55UP, G45UP, G54UP, RATIOUP, RATIOTOT, &
             &               DENOM, PHT, PHM, Z1, Z2, Z3, Z4
            
    Allocate (VEC1(NE), VEC2(NE),VHLP1(NE), UHLP1(NE), VHLP2(NE), UHLP2(NE), &
             &    V1(NE), V2(NE), U1(NE), U2(NE) )
   
    NF1 = NF

	ACCM = 0.D0
    do no = 1, Norb
        DO I = 1,LQ

            if( NFLAG == 2 ) then ! V2 > 0, (ni - nj)
                I4 = L_next(I,no,0   ) 
                I5 = L_next(I,no,NF1 ) 
                DEL44   =  DELLP2( NAUX_V2(I,NF1,NTAU) )   
                DEL55   =  DELLM2( NAUX_V2(I,NF1,NTAU) )  
            ENDIF

            G44UP = CMPLX(0.D0,0.D0)
            G45UP = CMPLX(0.D0,0.D0)
            G54UP = CMPLX(0.D0,0.D0)
            G55UP = CMPLX(0.D0,0.D0)
                ! ZGEMV
            DO NL = 1,NE
                VHLP1(NL) = CMPLX(DEL44,0.D0)*UR(I4,NL)
                VHLP2(NL) = CMPLX(DEL55,0.D0)*UR(I5,NL)
                UHLP1(NL) = UL(NL,I4)
                UHLP2(NL) = UL(NL,I5)
            ENDDO
            
            DO NL = 1,NE
                V1(NL) = CMPLX(0.D0,0.D0)
                V2(NL) = CMPLX(0.D0,0.D0)
                U1(NL) = CMPLX(0.D0,0.D0)
                U2(NL) = CMPLX(0.D0,0.D0)
            ENDDO
            DO NL = 1,NE
                DO NL1 = 1,NE
                V1(NL) = V1(NL) + VHLP1(NL1)*ULRINV(NL1,NL)
                V2(NL) = V2(NL) + VHLP2(NL1)*ULRINV(NL1,NL)
                ENDDO
            ENDDO

            DO NL = 1,NE
                G44UP = G44UP + V1(NL)*UHLP1(NL)
                G54UP = G54UP + V2(NL)*UHLP1(NL)
                G45UP = G45UP + V1(NL)*UHLP2(NL)
                G55UP = G55UP + V2(NL)*UHLP2(NL)
            ENDDO

            RATIOUP = (DCMPLX(1.D0,0.D0) + G44UP) * &
                        &    (DCMPLX(1.D0,0.D0) + G55UP) - &
                        &     G45UP*G54UP
                
            RATIOtot =    RATIOUP  

            if( NFLAG == 2 ) then ! V2 > 0, (ni - nj)
                RATIOtot =    RATIOUP
            ENDIF


            RATIO_RE_ABS = abs(RATIOtot)
            Random = RANF(ISEED)

            !    write(6,*) 'RATIOUP', RATIOUP
            !    write(6,*) 'RATIO_RE_ABS', RATIO_RE_ABS
                
                
            IF (RATIO_RE_ABS.GT.Random) THEN  
                ! WRITE(50,*) 'Accepted'   
                ! Upgrade the inverse
                    ACCM = ACCM + 1.D0
                    WEIGHT = SQRT(DBLE(RATIOTOT*DCONJG(RATIOTOT)))
                    phase = phase* ratiotot/dcmplx(weight,0.d0)
                    phasetot = phasetot + real(phase)
                    ncount = ncount + 1
                    
                    DO NL  = 1,NE
                        DO NL1 = 1,NE
                            U1(NL) = U1(NL) + ULRINV(NL,NL1)*UHLP1(NL1)
                            U2(NL) = U2(NL) + ULRINV(NL,NL1)*UHLP2(NL1)
                        ENDDO
                    ENDDO

                    Z1 =  CMPLX(1.D0,0.D0)/(CMPLX(1.D0,0.D0) + G55UP)
                    Z2 =  G54UP*Z1
                    Z3 =  G45UP*Z1
                    Z4 =  DCMPLX(1.D0,0.D0) + G44UP - G45UP*G54UP*Z1
                    Z4 =  DCMPLX(1.D0,0.D0)/Z4

                    DO NL = 1,NE
                        UHLP1(NL) = U2(NL)
                        VHLP1(NL) = V2(NL)*Z1
                        UHLP2(NL) = Z4*( U1(NL) - U2(NL)*Z2 )
                        VHLP2(NL) =      V1(NL) - V2(NL)*Z3
                    ENDDO
                
                    DO NL1 = 1,NE
                        DO NL2 = 1,NE
                            ULRINV(NL1,NL2) = ULRINV(NL1,NL2) &
                            &        -  UHLP1(NL1)*VHLP1(NL2)&
                            &        -  UHLP2(NL1)*VHLP2(NL2)
                        ENDDO
                    ENDDO
                    ! Upgrade  URUP
                    DO NL = 1,NE
                        UR(I4,NL) = UR(I4,NL) + &    
                        &     DCMPLX(DEL44,0.D0)*UR(I4,NL)
                        UR(I5,NL) = UR(I5,NL) +  &      
                        &      DCMPLX(DEL55,0.D0)*UR(I5,NL)   
                    ENDDO
                    !	      Flip:
                    if( NFLAG == 2 ) then ! V2 > 0, (ni - nj)
                        NAUX_V2(I,Nf1,NTAU) = - NAUX_V2(I,Nf1,NTAU)    
                    endif
            endif
        ENDDO
    enddo
	OBS(29) = OBS(29) + CMPLX(ACCM/DBLE(LFAM),0.D0)
	OBS(30) = OBS(30) + CMPLX(1.D0,0.D0)
    
    Deallocate (VEC1, VEC2,VHLP1, UHLP1, VHLP2, UHLP2,V1, V2, U1, U2 )
	 
END SUBROUTINE upgradeV
