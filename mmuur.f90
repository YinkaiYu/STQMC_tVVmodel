	SUBROUTINE MMUUR(A, NF, NTAU, NFLAG)

          Use Blockc
          Implicit Real (KIND=8) (A-G,O-Z)
          Implicit Integer (H-N)

          COMPLEX (Kind=8), Dimension(:,:) :: A 
          INTEGER :: NF, NT , NFLAG

          !	Local
          COMPLEX (Kind=8),  Dimension(:), Allocatable ::  V1, V2
          COMPLEX (Kind=8) ::  UT(2,2), U(2,2)

          N = SIZE(A,2)

          ALLOCATE (V1(N), V2(N))
             
          !          Kinetic.
          NF1 = NF
   
          !	WRITE(6,*) 'MMUUL: ',NF, NFLAG, NF1, NN
	
          IF (NFLAG.EQ.2 .AND.  RV2 .GT. ZERO ) THEN ! V2 > 0, (ni - nj)
            do no = 1, Norb
               DO I = 1,LQ
                  I1 = L_next  (I,no,0  )
                  I2 = L_next  (I,no,nf1)
                     
                  !Kinetic
                  DO J = 1,N
                     A(I1,J) = DCMPLX(XSIGP2(NAUX_V2(I,no,nf1,NTAU)),0.D0)*A(I1,J)
                     A(I2,J) = DCMPLX(XSIGM2(NAUX_V2(I,no,nf1,NTAU)),0.D0)*A(I2,J) 
                  ENDDO
               ENDDO
            enddo
          ENDIF
	
          IF (NFLAG.EQ.1 .AND.  RV1 .LT. -ZERO ) THEN ! V1 < 0, (ni + nj - 1)
             DO I = 1,LFAM
                I1 = L_bonds  (I,0  )
                I2 = L_bonds  (I,nf1)
                   
                !Kinetic
                DO J = 1,N
                    A(I1,J) = DCMPLX(XSIGP1(NAUX_V1(I,nf1,NTAU)),0.D0)*A(I1,J)
                    A(I2,J) = DCMPLX(XSIGP1(NAUX_V1(I,nf1,NTAU)),0.D0)*A(I2,J) 
                ENDDO
             ENDDO
          ENDIF
          
          DEALLOCATE (V1, V2)

          RETURN
    END SUBROUTINE MMUUR
    
  
        
