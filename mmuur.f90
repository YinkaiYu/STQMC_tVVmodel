	SUBROUTINE MMUUR(A, NF, NTAU, NFLAG)


          !	In A Out U(NF) * A                if NFLAG = 1
          !	In A Out EXP(D(NF)) * UT(NF) * A  if NFLAG = 2

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

          IF (NFLAG.EQ.3.and.RHUB.GT.zero) THEN
             DO NL = 1,N
                DO I = 1,Ndim
                   A(I,NL) = XSIGMA_U( NSIGL_U(I,NTAU)) * A(I,NL)
                ENDDO
             ENDDO
          ENDIF
	
          IF (NFLAG.EQ.2 .AND.  RJ .GT. ZERO ) THEN
             DO I = 1,LFAM
                I1 = L_bonds  (I,0  )
                I2 = L_bonds  (I,nf1)
                   
                !Kinetic
                DO J = 1,N
                    A(I1,J) = DCMPLX(XSIGP2(NSIGL_K(I,nf1,NTAU)),0.D0)*A(I1,J)
                    A(I2,J) = DCMPLX(XSIGM2(NSIGL_K(I,nf1,NTAU)),0.D0)*A(I2,J) 
                ENDDO
             ENDDO
          ENDIF
          
          DEALLOCATE (V1, V2)

          RETURN
    END SUBROUTINE MMUUR
    
  
        
