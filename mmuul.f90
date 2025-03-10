	SUBROUTINE MMUUL(A, NF,NTAU,NFLAG)
	
          !In A  Out A* EXP(D(NF)) * UT(NF)  if NFLAG = 1
          !In A  Out A* U(NF)                if NFLAG = 2

          Use Blockc
          Implicit Real (KIND=8) (A-G,O-Z)
          Implicit Integer (H-N)

          !Arguments:
          COMPLEX (Kind=8), Dimension(:,:) :: A
          Integer :: NF,NTAU,NFLAG

          !	Local
          COMPLEX (Kind=8), Dimension(:),  Allocatable :: V1, V2
          COMPLEX (Kind=8) :: UT(2,2), U(2,2)

          N = Size(A,1)
          ALLOCATE (V1(N), V2(N))
          IF (NFLAG.EQ.3.and.RHUB.GT.zero) THEN
             DO NL = 1,N
                DO J  = 1,Ndim
                   A(NL,J) = A(NL,J) * XSIGMA_U(NSIGL_U(J,NTAU))
                ENDDO
             ENDDO
          ENDIF
          
          !Kinetic.
          NF1 = NF
          
          IF (NFLAG.EQ.2 .AND. RV2.GT. ZERO) THEN ! V2 > 0, (ni - nj)
             DO I = 1,LFAM
                I1 = L_next(I,0  ) 
                I2 = L_next(I,nf1) 
                ! Kenitic
                DO J = 1,N
                    A(J,I1) = DCMPLX(XSIGP2(NAUX_V2(I,nf1,NTAU)),0.D0)*A(J,I1)
                    A(J,I2) = DCMPLX(XSIGM2(NAUX_V2(I,nf1,NTAU)),0.D0)*A(J,I2)   
                ENDDO
             ENDDO
          ENDIF
          
          IF (NFLAG.EQ.1 .AND. RV1.LT. -ZERO) THEN ! V1 < 0, (ni + nj - 1)
             DO I = 1,LFAM
                I1 = L_bonds(I,0  ) 
                I2 = L_bonds(I,nf1) 
                ! Kenitic
                DO J = 1,N
                    A(J,I1) = DCMPLX(XSIGP1(NAUX_V1(I,nf1,NTAU)),0.D0)*A(J,I1)
                    A(J,I2) = DCMPLX(XSIGP1(NAUX_V1(I,nf1,NTAU)),0.D0)*A(J,I2)   
                ENDDO
             ENDDO
          ENDIF

          DEALLOCATE (V1, V2)

          RETURN
    END SUBROUTINE MMUUL
    
