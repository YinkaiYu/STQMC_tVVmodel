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
        
        Real (Kind=8) :: aimj_p(2)
        Real (Kind=8) :: mom_x,mom_y,xk_p(2)
        Real (Kind=8) :: mom_shift_x, mom_shift_y, xk_shift_p(2)
        Integer :: RX_min, RX_max, RY_min, RY_max
        Complex(Kind=8) :: gk11, gk12

        Allocate(GRUP(Ndim,Ndim), GRUPC(Ndim,Ndim))
        
        ! Green functions.
        CALL CALCGR( UL, UR, ULRINV, GRUP, GRUPC)
        !GRUP (I,J) = <c_i c^+_j >
        !GRUPC (I,J) = <c^+_i c_j >

        ! temp QAH test

        do ix = 1, NLX
            do iy = 1, NLY
                do i_sublattice = 1, Norb
                    do i_direction = 1, Nnext
                        i = invlist(ix,iy)
                        ii_0 = L_next(i,i_sublattice,0)
                        ii_n = L_next(i,i_sublattice,i_direction)
                        R1 = (-1.0d0)**dble(i_sublattice+i_direction)
                        QAH_temp = QAH_temp + aimag(GRUPC(ii_0,ii_n) - GRUPC(ii_n,ii_0)) * R1 * real(phase) / dble(LQ)
                    enddo
                enddo
            enddo
        enddo

        ! QAH: peak at (0,0)

        mom_x = 0.0d0
        mom_y = 0.0d0
        xk_p(1) = mom_x * b1_p(1) + mom_y * b2_p(1)
        xk_p(2) = mom_x * b1_p(2) + mom_y * b2_p(2)
        mom_shift_x = 1.0d0/dble(NLx)
        mom_shift_y = 1.0d0/dble(NLy)
        xk_shift_p(1) = mom_shift_x * b1_p(1) + mom_shift_y * b2_p(1)
        xk_shift_p(2) = mom_shift_x * b1_p(2) + mom_shift_y * b2_p(2)

        do ix = 1, nlx
            do  iy = 1, nly
                do jx = 1, nlx
                    do jy = 1, nly
                        i = invlist(ix,iy)
                        j = invlist(jx,jy)
                        imj_x = ix-jx
                        imj_y = iy-jy
                        aimj_p(1) = dble(imj_x)* a1_p(1) + dble(imj_y)* a2_p(1)
                        aimj_p(2) = dble(imj_x)* a1_p(2) + dble(imj_y)* a2_p(2)
                        do i_sublattice = 1, Norb
                            do j_sublattice = 1, Norb
                                do i_direction = 1, 2 ! for x or y direction
                                    do j_direction = 1, 2
                                        ii_0 = L_next(i,i_sublattice,0)
                                        ii_n = L_next(i,i_sublattice,i_direction)
                                        jj_0 = L_next(j,j_sublattice,0)
                                        jj_n = L_next(j,j_sublattice,j_direction)
                                        S_QAH(i_sublattice,j_sublattice,i_direction,j_direction) = S_QAH(i_sublattice,j_sublattice,i_direction,j_direction) + real( phase * ( &
                                            &   GRUPC(ii_0,jj_0)*GRUP(ii_n,jj_n) + GRUPC(ii_0,ii_n)*GRUPC(jj_n,jj_0) &
                                            & - GRUPC(ii_0,jj_n)*GRUP(ii_n,jj_0) - GRUPC(ii_0,ii_n)*GRUPC(jj_0,jj_n) &
                                            & - GRUPC(ii_n,jj_0)*GRUP(ii_0,jj_n) - GRUPC(ii_n,ii_0)*GRUPC(jj_n,jj_0) &
                                            & + GRUPC(ii_n,jj_n)*GRUP(ii_0,jj_0) + GRUPC(ii_n,ii_0)*GRUPC(jj_0,jj_n) &
                                        & ) * exp( cmplx( 0.d0, xk_p(1)*aimj_p(1)+xk_p(2)*aimj_p(2) ) ) ) / dble(LQ)**2
                                        S_QAH_shift(i_sublattice,j_sublattice,i_direction,j_direction) = S_QAH_shift(i_sublattice,j_sublattice,i_direction,j_direction) + real( phase * ( &
                                            &   GRUPC(ii_0,jj_0)*GRUP(ii_n,jj_n) + GRUPC(ii_0,ii_n)*GRUPC(jj_n,jj_0) &
                                            & - GRUPC(ii_0,jj_n)*GRUP(ii_n,jj_0) - GRUPC(ii_0,ii_n)*GRUPC(jj_0,jj_n) &
                                            & - GRUPC(ii_n,jj_0)*GRUP(ii_0,jj_n) - GRUPC(ii_n,ii_0)*GRUPC(jj_n,jj_0) &
                                            & + GRUPC(ii_n,jj_n)*GRUP(ii_0,jj_0) + GRUPC(ii_n,ii_0)*GRUPC(jj_0,jj_n) &
                                        & ) * exp( cmplx( 0.d0, xk_shift_p(1)*aimj_p(1)+xk_shift_p(2)*aimj_p(2) ) ) ) / dble(LQ)**2
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        ! CM & VBS: peak at (\pm K), K is the Dirac point

        mom_x = 2.0d0/3.0d0
        mom_y = 1.0d0/3.0d0
        xk_p(1) = mom_x * b1_p(1) + mom_y * b2_p(1)
        xk_p(2) = mom_x * b1_p(2) + mom_y * b2_p(2)
        mom_shift_x = 2.0d0/3.0d0 + 1.0d0/dble(NLX)
        mom_shift_y = 1.0d0/3.0d0 - 1.0d0/dble(NLY)
        xk_shift_p(1) = mom_shift_x * b1_p(1) + mom_shift_y * b2_p(1)
        xk_shift_p(2) = mom_shift_x * b1_p(2) + mom_shift_y * b2_p(2)

        do ix = 1, nlx
            do  iy = 1, nly
                do jx = 1, nlx
                    do jy = 1, nly
                        i = invlist(ix,iy)
                        j = invlist(jx,jy)
                        imj_x = ix-jx
                        imj_y = iy-jy
                        aimj_p(1) = dble(imj_x)* a1_p(1) + dble(imj_y)* a2_p(1)
                        aimj_p(2) = dble(imj_x)* a1_p(2) + dble(imj_y)* a2_p(2)

                        ! CM
                        do i_sublattice = 1, Norb
                            do j_sublattice = 1, Norb
                                ii = invnlist(ix,iy,i_sublattice,1)
                                jj = invnlist(jx,jy,j_sublattice,1)
                                S_CM(i_sublattice,j_sublattice) = S_CM(i_sublattice,j_sublattice) &
                                    & + real( phase * (GRUPC(ii,jj)*GRUP(ii,jj) + (GRUPC(ii,ii)-0.5d0)*(GRUPC(jj,jj)-0.5d0) ) * exp( cmplx( 0.d0, xk_p(1)*aimj_p(1)+xk_p(2)*aimj_p(2) ) ) ) / dble(LQ)**2
                                S_CM_shift(i_sublattice,j_sublattice) = S_CM_shift(i_sublattice,j_sublattice) &
                                    & + real( phase * (GRUPC(ii,jj)*GRUP(ii,jj) + (GRUPC(ii,ii)-0.5d0)*(GRUPC(jj,jj)-0.5d0) ) * exp( cmplx( 0.d0, xk_shift_p(1)*aimj_p(1)+xk_shift_p(2)*aimj_p(2) ) ) ) / dble(LQ)**2
                            enddo
                        enddo
                        
                        ! VBS
                        do i_bond = 1, Nbond
                            do j_bond = 1, Nbond
                                ii_0 = L_bonds(i,0)
                                ii_n = L_bonds(i,i_bond)
                                jj_0 = L_bonds(j,0)
                                jj_n = L_bonds(j,j_bond)
                                S_VBS(i_bond,j_bond) = S_VBS(i_bond,j_bond) + real( phase * ( &
                                    &   GRUPC(ii_0,jj_n)*GRUP(ii_n,jj_0) + GRUPC(ii_0,ii_n)*GRUPC(jj_0,jj_n) &
                                    & + GRUPC(ii_0,jj_0)*GRUP(ii_n,jj_n) + GRUPC(ii_0,ii_n)*GRUPC(jj_n,jj_0) &
                                    & + GRUPC(ii_n,jj_n)*GRUP(ii_0,jj_0) + GRUPC(ii_n,ii_0)*GRUPC(jj_0,jj_n) &
                                    & + GRUPC(ii_n,jj_0)*GRUP(ii_0,jj_n) + GRUPC(ii_n,ii_0)*GRUPC(jj_n,jj_0) &
                                & ) * exp( cmplx( 0.d0, xk_p(1)*aimj_p(1)+xk_p(2)*aimj_p(2) ) ) ) / dble(LQ)**2
                                S_VBS_shift(i_bond,j_bond) = S_VBS_shift(i_bond,j_bond) + real( phase * ( &
                                    &   GRUPC(ii_0,jj_n)*GRUP(ii_n,jj_0) + GRUPC(ii_0,ii_n)*GRUPC(jj_0,jj_n) &
                                    & + GRUPC(ii_0,jj_0)*GRUP(ii_n,jj_n) + GRUPC(ii_0,ii_n)*GRUPC(jj_n,jj_0) &
                                    & + GRUPC(ii_n,jj_n)*GRUP(ii_0,jj_0) + GRUPC(ii_n,ii_0)*GRUPC(jj_0,jj_n) &
                                    & + GRUPC(ii_n,jj_0)*GRUP(ii_0,jj_n) + GRUPC(ii_n,ii_0)*GRUPC(jj_n,jj_0) &
                                & ) * exp( cmplx( 0.d0, xk_shift_p(1)*aimj_p(1)+xk_shift_p(2)*aimj_p(2) ) ) ) / dble(LQ)**2
                            enddo
                        enddo

                    enddo
                enddo
            enddo
        enddo

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

        mom_x = 2.0d0/3.0d0 + 1.0d0/dble(NLx)
        mom_y = 1.0d0/3.0d0 - 1.0d0/dble(NLy)
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

        do ix = 1, nlx
            do iy = 1, nly
                iA = invnlist(ix,iy,1,1)
                iB = invnlist(ix,iy,2,1)
                jx = npbcx(ix+NLX/3)
                jy = npbcy(iy+NLY/3)
                jA = invnlist(jx,jy,1,1)
                jB = invnlist(jx,jy,2,1)
                fermicor11_onethird = fermicor11_onethird + real( GRUPC(iA,jA) + GRUPC(jA,iA) ) / dble(LQ)
                fermicor12_onethird = fermicor12_onethird + real( GRUPC(iA,jB) + GRUPC(jB,iA) ) / dble(LQ)
                fermicor21_onethird = fermicor21_onethird + real( GRUPC(iB,jA) + GRUPC(jA,iB) ) / dble(LQ)
                fermicor22_onethird = fermicor22_onethird + real( GRUPC(iB,jB) + GRUPC(jB,iB) ) / dble(LQ)
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


   