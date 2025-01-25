	SUBROUTINE SALPH(alpha_U)
!       Calculate DTAU and ALPHA

          Use Blockc
          Implicit Real (KIND=8) (A-G,O-Z)
          Implicit Integer (H-N)

          REAL (Kind=8) ::  XTH(4)
          COMPLEX (Kind = 8) :: Z, ALPHA_U
          COMPLEX (Kind = 8) :: UHLP(4,4)


          !For J term.
          ETAL(-2)= - SQRT(2.D0 * ( 3.D0 + SQRT(6.D0) ) )
          ETAL(-1)= - SQRT(2.D0 * ( 3.D0 - SQRT(6.D0) ) )
          ETAL( 1)=   SQRT(2.D0 * ( 3.D0 - SQRT(6.D0) ) )
          ETAL( 2)=   SQRT(2.D0 * ( 3.D0 + SQRT(6.D0) ) )
          
          ! For flipping 
          NFLIPL(-1) =  1
          NFLIPL( 1) =  -1
          
          ALPHA  =  acosh(exp(RJ*dtau/2.d0))
          
          ! Same for current and for kinetic.
          XSIGP2(-1) = EXP( 1.D0* ALPHA * dble(-1))
          XSIGP2( 1) = EXP( 1.D0* ALPHA * dble(1) )
          
          XSIGM2(-1) = EXP(-1.D0* ALPHA * dble(-1))
          XSIGM2( 1) = EXP(-1.D0* ALPHA * dble(1) )
          
          DO NL = -1,1
             IF (NL.NE.0) THEN
                 IF(NL==1) NLN = -1
                 IF(NL==-1) NLN = 1
                 DELLP2(NL) =  (XSIGP2(NLN)/XSIGP2(NL)) - 1.D0
                 DELLM2(NL) =  (XSIGM2(NLN)/XSIGM2(NL)) - 1.D0                  
             ENDIF
          ENDDO
          
          !	Pair-hopping. Kinetic.
          DO M = 1,2
             DO N = 1,2
                UR_K(M,N) = DCMPLX(0.D0,0.D0)
             ENDDO
          ENDDO
          UR_K(1,1) =   DCMPLX( 1.D0/SQRT(2.D0),0.D0 )
          UR_K(1,2) =   DCMPLX( 1.D0/SQRT(2.D0),0.D0 )
          UR_K(2,1) =   DCMPLX( 1.D0/SQRT(2.D0),0.D0 )
          UR_K(2,2) =   DCMPLX(-1.D0/SQRT(2.D0),0.D0 )
          
          DO M = 1,2
             DO N = 1,2
                URT_K(M,N) = DCONJG(UR_K(N,M))
             ENDDO
          ENDDO
          
          !****** For Hubbard.
          ALPHA_U = dcmplx(acosh(exp(RHUB*dtau/2.0)),0.d0)

          XSIGMA_U(-1) = EXP(ALPHA_U*dcmplx(-1.d0,0.d0))
          XSIGMA_U( 1) = EXP(ALPHA_U*dcmplx(1.d0,0.d0))
          
          RHO = 0.5
          DO NL = -1,1
             IF (NL.NE.0) THEN
                   NLN = NFLIPL(NL)
                   DELTA_U(NL)=(XSIGMA_U(NLN)/XSIGMA_U(NL)) -  DCMPLX(1.D0,0.D0)
             ENDIF
          ENDDO
          
          RETURN
        END SUBROUTINE SALPH
