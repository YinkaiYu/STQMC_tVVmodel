      SUBROUTINE OBSER(UL,UR,ULRINV,phase)
        Use Blockc
        Use Block_obs
        Use MyMats
        Implicit Real (KIND=8) (A-G,O-Z)
        Implicit Integer (H-N)
!#define DEC
        Interface
           SUBROUTINE CALCGR(UL, UR, ULRINV, GRUP, GRUPC)
             COMPLEX (Kind=8), Dimension(:,:) :: UL, UR, ULRINV
             COMPLEX (Kind=8), Dimension(:,:) :: GRUP, GRUPC
           END SUBROUTINE CALCGR
        end Interface

        Complex(Kind=8) :: Zp, phase
        Complex(Kind=8), Dimension(:,:) :: UL, UR, ULRINV
        Complex(Kind=8), Dimension(:,:), Allocatable :: GRUP, GRUPC
        
        Real (Kind=8) :: mom_x,mom_y,xk_p(2), aimj_p(2)
        Integer :: RX_min, RX_max, RY_min, RY_max
        Complex(Kind=8) :: gk11, gk12

        Allocate(GRUP(Ndim,Ndim), GRUPC(Ndim,Ndim))
        
        ! Green functions.
        CALL CALCGR( UL, UR, ULRINV, GRUP, GRUPC)
        !GRUP (I,J) = <c_i c^+_j >
        !GRUPC (I,J) = <c^+_i c_j >

! QAH & CM
        ! i_sublattice   = 1
        ! i_NNNdirection = 1
        ! j_sublattice   = 1
        ! j_NNNdirection = 1
        ! do ix = 1, nlx
        !     do  iy = 1, nly
        !         do jx = 1, nlx
        !             do jy = 1, nly
        !                 ii_0 = invnlist(ix,iy,i_sublattice,1)
        !                 jj_0 = invnlist(jx,jy,j_sublattice,1)
        !                 ii_n = 
        !                 imj_x = ix-jx
        !                 imj_y = iy-jy
        !             enddo
        !         enddo
        !     enddo
        ! enddo

! old observables
        
        do ii = 1, LQ
            i_0 = L_bonds(ii,0) 
            do nf = 1, Nbond
                i_n = L_bonds(ii,nf)
                kinetic = kinetic - RT1*real(GRUPC(i_0,i_n) + GRUPC(i_n,i_0))*real(phase)
                potential = potential - RV1*real(GRUPC(i_0,i_0)*GRUPC(i_n,i_n) + GRUPC(i_0,i_n)*GRUP(i_n,i_0))*real(phase)
            enddo
        enddo
        
        do ii = 1, Ndim
            density = density + real(GRUPC(ii,ii))*real(phase) / dble(Ndim)
        enddo

        mom_x = 1.0d0/dble(NLx)
        mom_y = 1.0d0/dble(NLy)
        xk_p(1) = mom_x * b1_p(1) + mom_y * b2_p(1)
        xk_p(2) = mom_x * b1_p(2) + mom_y * b2_p(2)

        do ix = 1, nlx
            do  iy = 1, nly
                do jx = 1, nlx
                    do jy = 1, nly
                        imj_x = ix-jx
                        imj_y = iy-jy
                        aimj_p(1) = dble(imj_x)* a1_p(1) + dble(imj_y)* a2_p(1)
                        aimj_p(2) = dble(imj_x)* a1_p(2) + dble(imj_y)* a2_p(2)

                        
                        II = invnlist(ix,iy,1,1)
                        JJ = invnlist(jx,jy,1,1)
                        S0_11 = S0_11 + real( phase * (GRUPC(II,JJ)*GRUP(II,JJ) + (GRUPC(ii,ii)-0.5d0)*(GRUPC(jj,jj)-0.5d0) ) ) / dble(LQ)**2
                        Sk_11 = Sk_11 + real( phase * (GRUPC(II,JJ)*GRUP(II,JJ) + (GRUPC(ii,ii)-0.5d0)*(GRUPC(jj,jj)-0.5d0) ) * exp( cmplx( 0.d0, xk_p(1)*aimj_p(1)+xk_p(2)*aimj_p(2) ) ) ) / dble(LQ)**2

                        II = invnlist(ix,iy,1,1)
                        JJ = invnlist(jx,jy,2,1)
                        S0_12 = S0_12 + real( phase * (GRUPC(II,JJ)*GRUP(II,JJ) + (GRUPC(ii,ii)-0.5d0)*(GRUPC(jj,jj)-0.5d0) ) ) / dble(LQ)**2
                        Sk_12 = Sk_12 + real( phase * (GRUPC(II,JJ)*GRUP(II,JJ) + (GRUPC(ii,ii)-0.5d0)*(GRUPC(jj,jj)-0.5d0) ) * exp( cmplx( 0.d0, xk_p(1)*aimj_p(1)+xk_p(2)*aimj_p(2) ) ) ) / dble(LQ)**2

                        II = invnlist(ix,iy,2,1)
                        JJ = invnlist(jx,jy,1,1)
                        S0_21 = S0_21 + real( phase * (GRUPC(II,JJ)*GRUP(II,JJ) + (GRUPC(ii,ii)-0.5d0)*(GRUPC(jj,jj)-0.5d0) ) ) / dble(LQ)**2
                        Sk_21 = Sk_21 + real( phase * (GRUPC(II,JJ)*GRUP(II,JJ) + (GRUPC(ii,ii)-0.5d0)*(GRUPC(jj,jj)-0.5d0) ) * exp( cmplx( 0.d0, xk_p(1)*aimj_p(1)+xk_p(2)*aimj_p(2) ) ) ) / dble(LQ)**2

                        II = invnlist(ix,iy,2,1)
                        JJ = invnlist(jx,jy,2,1)
                        S0_22 = S0_22 + real( phase * (GRUPC(II,JJ)*GRUP(II,JJ) + (GRUPC(ii,ii)-0.5d0)*(GRUPC(jj,jj)-0.5d0) ) ) / dble(LQ)**2
                        Sk_22 = Sk_22 + real( phase * (GRUPC(II,JJ)*GRUP(II,JJ) + (GRUPC(ii,ii)-0.5d0)*(GRUPC(jj,jj)-0.5d0) ) * exp( cmplx( 0.d0, xk_p(1)*aimj_p(1)+xk_p(2)*aimj_p(2) ) ) ) / dble(LQ)**2

                    enddo
                enddo
            enddo
        enddo

        mom_x = 2.0d0/3.0d0 + 1.0d0*PI/dble(NLx)
        mom_y = 1.0d0/3.0d0 - 1.0d0*PI/dble(NLy)
        xk_p(1) = mom_x * b1_p(1) + mom_y * b2_p(1)
        xk_p(2) = mom_x * b1_p(2) + mom_y * b2_p(2)

        do ix = 1, nlx
            do  iy = 1, nly
                do jx = 1, nlx
                    do jy = 1, nly  
                        i     = invlist(ix,iy)                
                        j     = invlist(jx,jy)
                        iup   = invnlist(ix,iy,1,1)
                        ido   = invnlist(ix,iy,2,1)
                        jup   = invnlist(jx,jy,1,1)
                        jdo   = invnlist(jx,jy,2,1)
                        imj_x = ix-jx ! no npbcx here
                        imj_y = iy-jy ! no npbcy here
                        imj   = invlist(npbcx(imj_x),npbcx(imj_y))
                        ! fermion correlation function (PRL 128, 225701 (2022))
                        aimj_p(1) = dble(imj_x)* a1_p(1) + dble(imj_y)* a2_p(1)
                        aimj_p(2) = dble(imj_x)* a1_p(2) + dble(imj_y)* a2_p(2)
                        fermicor11_deltaq = fermicor11_deltaq + real( phase *  exp( cmplx( 0.d0, xk_p(1)*aimj_p(1)+xk_p(2)*aimj_p(2) ) ) * GRUPC(iup,jup) ) / dble(LQ)
                        fermicor12_deltaq = fermicor12_deltaq + real( phase *  exp( cmplx( 0.d0, xk_p(1)*aimj_p(1)+xk_p(2)*aimj_p(2) ) ) * GRUPC(iup,jdo) ) / dble(LQ)
                        fermicor21_deltaq = fermicor21_deltaq + real( phase *  exp( cmplx( 0.d0, xk_p(1)*aimj_p(1)+xk_p(2)*aimj_p(2) ) ) * GRUPC(ido,jup) ) / dble(LQ)
                        fermicor22_deltaq = fermicor22_deltaq + real( phase *  exp( cmplx( 0.d0, xk_p(1)*aimj_p(1)+xk_p(2)*aimj_p(2) ) ) * GRUPC(ido,jdo) ) / dble(LQ)
                    enddo
                enddo    
            enddo
        enddo


        ! do ix = 1, nlx
        !     do  iy = 1, nly
        !         do no1 = 1,norb
        !             do jx = 1, nlx
        !                 do jy = 1, nly                  
        !                     do no2 = 1,norb
        !                         II = invnlist(ix,iy,no1,1)
        !                         JJ = invnlist(jx,jy,no2,1)
        !                         imj_x = ix-jx
        !                         imj_y = iy-jy
        !                         den(imj_x,imj_y,no1,no2) = den(imj_x,imj_y,no1,no2) + &
        !                         &    phase*(GRUPC(II,JJ)*GRUP(II,JJ)* dble(N_sun) + &
        !                         & (GRUPC(ii,ii)-0.5d0)*(GRUPC(jj,jj)-0.5d0)* dble(N_sun*N_sun)) /cmplx(LQ,0.d0) 
        !                     enddo
        !                 enddo
        !             enddo        
        !         enddo    
        !     enddo
        ! enddo  ! in the SU(N) symmetric model, spin/spinpm are identical.
        
        ! write(6,*) 'matrix of GRUPC'
        ! WRITE(6, '(A7)',advance='no') "Col/Row" 
        ! DO j = 1, Ndim
        !     WRITE(6,'(7X,I5,8X)', advance='no') j
        ! ENDDO
        ! WRITE(6,*) ! \n
        ! DO i = 1, Ndim
        !     WRITE(6, '(I5,2X)',advance='no') i 
        !     DO j = 1, Ndim
        !         WRITE(6,'( "(", F8.4, ",", F8.4, ") " )', advance='no') REAL(GRUPC(i,j)), AIMAG(GRUPC(i,j))
        !     ENDDO
        !     WRITE(6,*) ! \n
        ! ENDDO

        ! write(6,*) 'matrix of den11'
        ! RX_max = NLX-1
        ! RX_min = 1-NLX
        ! RY_max = NLY-1
        ! RY_min = 1-NLY
        ! WRITE(6, '(A7)',advance='no') "Col/Row" 
        ! DO j = RY_min, RY_max
        !     WRITE(6,'(7X,I5,8X)', advance='no') j
        ! ENDDO
        ! WRITE(6,*) ! \n
        ! DO i = RX_min, RX_max
        !     WRITE(6, '(I5,2X)',advance='no') i 
        !     DO j = RY_min, RY_max
        !         WRITE(6,'( "(", F8.4, ",", F8.4, ") " )', advance='no') REAL(den(i,j,1,1))*1000.d0, AIMAG(den(i,j,1,1))*1000.d0
        !     ENDDO
        !     WRITE(6,*) ! \n
        ! ENDDO

        ! write(6,*) 'matrix of den12'
        ! RX_max = NLX-1
        ! RX_min = 1-NLX
        ! RY_max = NLY-1
        ! RY_min = 1-NLY
        ! WRITE(6, '(A7)',advance='no') "Col/Row" 
        ! DO j = RY_min, RY_max
        !     WRITE(6,'(7X,I5,8X)', advance='no') j
        ! ENDDO
        ! WRITE(6,*) ! \n
        ! DO i = RX_min, RX_max
        !     WRITE(6, '(I5,2X)',advance='no') i 
        !     DO j = RY_min, RY_max
        !         WRITE(6,'( "(", F8.4, ",", F8.4, ") " )', advance='no') REAL(den(i,j,1,2))*1000.d0, AIMAG(den(i,j,1,2))*1000.d0
        !     ENDDO
        !     WRITE(6,*) ! \n
        ! ENDDO
       
       Deallocate ( GRUP , GRUPC )
      END SUBROUTINE OBSER


   