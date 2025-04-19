	SUBROUTINE MMUURM1(A, NF, NTAU,  NFLAG)
          !	In A Out U(NF) * A                if NFLAG = 1
          !	In A Out EXP(D(NF)) * UT(NF) * A  if NFLAG = 2

          Use Blockc
          Implicit Real (KIND=8) (A-G,O-Z)
          Implicit Integer (H-N)
          COMPLEX (KIND=8), Dimension(:,:) :: A
          Integer :: NF, NTAU,  NFLAG

          !Local
          COMPLEX (KIND=8), Dimension(:), Allocatable ::  V1, V2
          COMPLEX (KIND=8) :: UT(2,2), U(2,2)

          N = Size(A,2)

          Allocate (V1(N), V2(N))
             
          !          Kinetic.
          NF1 = NF
          
          IF (NFLAG.EQ.2 .AND. RV2 .GT. ZERO) THEN ! V2 > 0, (ni - nj)
            do no = 1, Norb
               DO I = 1,LFAM
                  I1 = L_next(I,no,0  )
                  I2 = L_next(I,no,nf1)
                  ! Kinetic 
                  DO J = 1,N
                     A(I1,J) =  A(I1,J) / DCMPLX(XSIGP2(NAUX_V2(I,Nf1,NTAU)),0.D0)   
                     A(I2,J) =  A(I2,J) / DCMPLX(XSIGM2(NAUX_V2(I,Nf1,NTAU)),0.D0)
                  ENDDO
               ENDDO
            enddo
          ENDIF
          
          IF (NFLAG.EQ.1 .AND. RV1 .LT. -ZERO) THEN ! V1 < 0, (ni + nj - 1)
             DO I = 1,LFAM
                I1 = L_Bonds(I,0  )
                I2 = L_Bonds(I,nf1)
                ! Kinetic 
                DO J = 1,N
                    A(I1,J) =  A(I1,J) / DCMPLX(XSIGP1(NAUX_V1(I,Nf1,NTAU)),0.D0)   
                    A(I2,J) =  A(I2,J) / DCMPLX(XSIGP1(NAUX_V1(I,Nf1,NTAU)),0.D0)
                ENDDO
             ENDDO
          ENDIF
          
          Deallocate (V1, V2)
          
          Return
    END SUBROUTINE MMUURM1
        
