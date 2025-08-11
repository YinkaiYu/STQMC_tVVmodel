SUBROUTINE PREQ(NOBS,Nobs_tot)
     Use Blockc
     Use Block_obs
     Implicit Real (KIND=8) (A-G,O-Z)
     Implicit Integer (H-N)
     COMPLEX (KIND =8) :: Znorm
     COMPLEX (KIND =8) :: Znorm1
     COMPLEX (KIND=8), Dimension(:,:,:,:), Allocatable :: Collect2
     COMPLEX (KIND=8), Dimension(:), Allocatable :: Collect1
     real(kind=8) :: collect3
     complex (kind=8) :: collect4
     real(kind=8), dimension(:,:,:,:), Allocatable :: real_4index
     real(kind=8), dimension(:,:), Allocatable :: real_2index
     Character(32) :: filek,filek1,filek2
     real(kind=8) :: mom_x, mom_y, S_temp
     Integer NOBS_TOT
     Interface
     Subroutine Fourier_Trans(gr,filek)
          Complex (Kind=8), dimension(:,:,:) :: gr
          Character (32) :: filek
     end Subroutine Fourier_Trans
     Subroutine  correlation(gr,filek)
     Complex (Kind=8), dimension(:,:,:) :: gr
          Character (32) :: filek
     end Subroutine correlation   
     Subroutine structurefactor(gr,filek,mom_x,mom_y,no1,no2)
     Complex (Kind=8), dimension(:,:,:,:) :: gr
     Integer :: lp, no1, no2
     Character (32) :: filek
     Real (Kind=8) :: mom_x,mom_y,xk_p(2), aimj_p(2)
     Integer :: RX_min, RX_max, RY_min, RY_max
     end subroutine structurefactor

     end Interface
     !#define DEC
     INCLUDE 'mpif.h'
     INTEGER STATUS(MPI_STATUS_SIZE)
     Integer :: RX_min, RX_max, RY_min, RY_max
     RX_max = NLX-1
     RX_min = 1-NLX
     RY_max = NLY-1
     RY_min = 1-NLY
     CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
     CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
     ZNORM = CMPLX(1.D0,0.D0)/ DCMPLX( DBLE(NOBS), 0.D0 )
     ZNORM1 = CMPLX(1.D0,0.D0)/ DCMPLX( DBLE(NOBS_TOT),0.D0 )

     DEN = ZNORM* DEN
     kinetic = znorm*kinetic; potential = znorm*potential
     density = znorm* density
     QAH_temp = znorm* QAH_temp
     M2 = znorm* M2
     S0_11 = znorm* S0_11
     S0_12 = znorm* S0_12
     S0_21 = znorm* S0_21
     S0_22 = znorm* S0_22
     Sk_11 = znorm* Sk_11
     Sk_12 = znorm* Sk_12
     Sk_21 = znorm* Sk_21
     Sk_22 = znorm* Sk_22
     S_QAH = znorm* S_QAH
     S_CM  = znorm* S_CM
     S_VBS = znorm* S_VBS
     S_QAH_shift = znorm* S_QAH_shift
     S_CM_shift  = znorm* S_CM_shift
     S_VBS_shift = znorm* S_VBS_shift
     fermicor11_deltaq = fermicor11_deltaq / DBLE(NOBS)
     fermicor12_deltaq = fermicor12_deltaq / DBLE(NOBS)
     fermicor21_deltaq = fermicor21_deltaq / DBLE(NOBS)
     fermicor22_deltaq = fermicor22_deltaq / DBLE(NOBS)
     fermicor11_onethird = fermicor11_onethird / DBLE(NOBS)
     fermicor12_onethird = fermicor12_onethird / DBLE(NOBS)
     fermicor21_onethird = fermicor21_onethird / DBLE(NOBS)
     fermicor22_onethird = fermicor22_onethird / DBLE(NOBS)
     phasetot = phasetot/dble(Ncount)

     !Collect.
     N = 1
     Collect3 = 0.0d0
     CALL MPI_REDUCE(kinetic,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     kinetic = Collect3/DBLE(ISIZE)

     N = 1
     Collect3 =0.0d0
     CALL MPI_REDUCE(potential,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     potential = Collect3/DBLE(ISIZE)

     N = 1
     Collect3 =0.0d0
     CALL MPI_REDUCE(M2,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     M2 = Collect3/DBLE(ISIZE)

     N = 1
     Collect3 =0.0d0
     CALL MPI_REDUCE(density,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     density = Collect3/DBLE(ISIZE)

     N = 1
     Collect3 =0.0d0
     CALL MPI_REDUCE(QAH_temp,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     QAH_temp = Collect3/DBLE(ISIZE)

     N = 1
     Collect3 =0.0d0
     CALL MPI_REDUCE(S0_11,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     S0_11 = Collect3/DBLE(ISIZE)
     Collect3 =0.0d0
     CALL MPI_REDUCE(S0_12,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     S0_12 = Collect3/DBLE(ISIZE)
     Collect3 =0.0d0
     CALL MPI_REDUCE(S0_21,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     S0_21 = Collect3/DBLE(ISIZE)
     Collect3 =0.0d0
     CALL MPI_REDUCE(S0_22,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     S0_22 = Collect3/DBLE(ISIZE)

     N = 1
     Collect3 =0.0d0
     CALL MPI_REDUCE(Sk_11,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     Sk_11 = Collect3/DBLE(ISIZE)
     Collect3 =0.0d0
     CALL MPI_REDUCE(Sk_12,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     Sk_12 = Collect3/DBLE(ISIZE)
     Collect3 =0.0d0
     CALL MPI_REDUCE(Sk_21,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     Sk_21 = Collect3/DBLE(ISIZE)
     Collect3 =0.0d0
     CALL MPI_REDUCE(Sk_22,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     Sk_22 = Collect3/DBLE(ISIZE)

     allocate(real_4index(Norb,Norb,2,2))
     N = Norb*Norb*2*2

     real_4index = 0.d0
     CALL MPI_REDUCE(S_QAH,real_4index,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     S_QAH = real_4index/DBLE(ISIZE)
     
     real_4index = 0.d0
     CALL MPI_REDUCE(S_QAH_shift,real_4index,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     S_QAH_shift = real_4index/DBLE(ISIZE)

     allocate(real_2index(Norb,Norb))
     N = Norb*Norb

     real_2index = 0.d0
     CALL MPI_REDUCE(S_CM,real_2index,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     S_CM = real_2index/DBLE(ISIZE)
     
     real_2index = 0.d0
     CALL MPI_REDUCE(S_CM_shift,real_2index,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     S_CM_shift = real_2index/DBLE(ISIZE)

     deallocate(real_2index)
     allocate(real_2index(Nbond,Nbond))
     N = Nbond*Nbond

     real_2index = 0.d0
     CALL MPI_REDUCE(S_VBS,real_2index,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     S_VBS = real_2index/DBLE(ISIZE)
     
     real_2index = 0.d0
     CALL MPI_REDUCE(S_VBS_shift,real_2index,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     S_VBS_shift = real_2index/DBLE(ISIZE)

     Allocate(Collect2(RX_min:RX_max, RY_min:RY_max, norb, norb))
     N = (2*NLX-1)*(2*NLY-1)*Norb*Norb
     Collect2 = CMPLX(0.D0,0.D0)
     CALL MPI_REDUCE(DEN,Collect2,N,MPI_COMPLEX16,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     DEN = Collect2/CMPLX(DBLE(ISIZE),0.D0)

     Collect3 = 0.d0; N=1
     CALL MPI_REDUCE(phasetot,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     phasetot = Collect3/CMPLX(DBLE(ISIZE),0.D0)

     N = 1
     Collect3 =0.0d0
     CALL MPI_REDUCE(fermicor11_deltaq,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     fermicor11_deltaq = Collect3/DBLE(ISIZE)
     Collect3 =0.0d0
     CALL MPI_REDUCE(fermicor12_deltaq,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     fermicor12_deltaq = Collect3/DBLE(ISIZE)
     Collect3 =0.0d0
     CALL MPI_REDUCE(fermicor21_deltaq,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     fermicor21_deltaq = Collect3/DBLE(ISIZE)
     Collect3 =0.0d0
     CALL MPI_REDUCE(fermicor22_deltaq,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     fermicor22_deltaq = Collect3/DBLE(ISIZE)

     N = 1
     Collect3 =0.0d0
     CALL MPI_REDUCE(fermicor11_onethird,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     fermicor11_onethird = Collect3/DBLE(ISIZE)
     Collect3 =0.0d0
     CALL MPI_REDUCE(fermicor12_onethird,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     fermicor12_onethird = Collect3/DBLE(ISIZE)
     Collect3 =0.0d0
     CALL MPI_REDUCE(fermicor21_onethird,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     fermicor21_onethird = Collect3/DBLE(ISIZE)
     Collect3 =0.0d0
     CALL MPI_REDUCE(fermicor22_onethird,Collect3,N,MPI_REAL8,MPI_SUM,&
          & 0,MPI_COMM_WORLD,IERR)
     fermicor22_onethird = Collect3/DBLE(ISIZE)


     IF (IRANK.EQ.0) THEN

          ! mom_x = 0.d0
          ! mom_y = 0.d0
          ! filek = "cdw11"
          ! call structurefactor(den,filek,mom_x,mom_y,1,1)
          ! filek = "cdw12"
          ! call structurefactor(den,filek,mom_x,mom_y,1,2)
          ! filek = "cdw21"
          ! call structurefactor(den,filek,mom_x,mom_y,2,1)
          ! filek = "cdw22"
          ! call structurefactor(den,filek,mom_x,mom_y,2,2)

          !    mom_x = 1.0d0/dble(NLX)
          !    mom_y = 1.0d0/dble(NLY)
          !    filek = "cdw11_dk"
          !    call structurefactor(den,filek,mom_x,mom_y,1,1)
          !    filek = "cdw12_dk"
          !    call structurefactor(den,filek,mom_x,mom_y,1,2)
          ! filek = "cdw21_dk"
          ! call structurefactor(den,filek,mom_x,mom_y,2,1)
          ! filek = "cdw22_dk"
          ! call structurefactor(den,filek,mom_x,mom_y,2,2)


          filek = "phasetot"
          call orderparameter(phasetot,filek)
          filek = "kinetic"
          call orderparameter(kinetic,filek)
          filek = "potential"
          call orderparameter(potential,filek)
          filek = "density"
          call orderparameter(density,filek)
          filek = "QAH_temp"
          call orderparameter(QAH_temp,filek)
          filek = "S0_11"
          call orderparameter(S0_11,filek)
          filek = "S0_12"
          call orderparameter(S0_12,filek)
          filek = "S0_21"
          call orderparameter(S0_21,filek)
          filek = "S0_22"
          call orderparameter(S0_22,filek)
          filek = "Sk_11"
          call orderparameter(Sk_11,filek)
          filek = "Sk_12"
          call orderparameter(Sk_12,filek)
          filek = "Sk_21"
          call orderparameter(Sk_21,filek)
          filek = "Sk_22"
          call orderparameter(Sk_22,filek)
          
          do i_sublattice = 1, Norb
               do j_sublattice = 1, Norb
                    do i_direction = 1, 2 ! for x or y direction
                         do j_direction = 1, 2

                              write(filek, "('S_QAH_','sub',I0,'sub',I0,'drc',I0,'drc',I0)") i_sublattice, j_sublattice, i_direction, j_direction
                              S_temp = S_QAH(i_sublattice,j_sublattice,i_direction,j_direction)
                              call orderparameter(S_temp,filek)

                              write(filek, "('S_QAH_','sub',I0,'sub',I0,'drc',I0,'drc',I0,'_shift')") i_sublattice, j_sublattice, i_direction, j_direction
                              S_temp = S_QAH_shift(i_sublattice,j_sublattice,i_direction,j_direction)
                              call orderparameter(S_temp,filek)

                         enddo
                    enddo
               enddo
          enddo

          do i_sublattice = 1, Norb
               do j_sublattice = 1, Norb

                    write(filek, "('S_CM','sub',I0,'sub',I0)") i_sublattice, j_sublattice
                    S_temp = S_CM(i_sublattice,j_sublattice)
                    call orderparameter(S_temp,filek)

                    write(filek, "('S_CM','sub',I0,'sub',I0,'_shift')") i_sublattice, j_sublattice
                    S_temp = S_CM_shift(i_sublattice,j_sublattice)
                    call orderparameter(S_temp,filek)

               enddo
          enddo

          do i_bond = 1, Nbond
               do j_bond = 1, Nbond

                    write(filek, "('S_VBS','sub',I0,'sub',I0)") i_bond, j_bond
                    S_temp = S_VBS(i_bond,j_bond)
                    call orderparameter(S_temp,filek)

                    write(filek, "('S_VBS','sub',I0,'sub',I0,'_shift')") i_bond, j_bond
                    S_temp = S_VBS_shift(i_bond,j_bond)
                    call orderparameter(S_temp,filek)

               enddo
          enddo
          
          filek = "fermicor11_dk"
          Call orderparameter(fermicor11_deltaq, filek )
          filek = "fermicor12_dk"
          Call orderparameter(fermicor12_deltaq, filek )
          filek = "fermicor21_dk"
          Call orderparameter(fermicor21_deltaq, filek )
          filek = "fermicor22_dk"
          Call orderparameter(fermicor22_deltaq, filek )

          filek = "fermicor11_onethird"
          Call orderparameter(fermicor11_onethird, filek )
          filek = "fermicor12_onethird"
          Call orderparameter(fermicor12_onethird, filek )
          filek = "fermicor21_onethird"
          Call orderparameter(fermicor21_onethird, filek )
          filek = "fermicor22_onethird"
          Call orderparameter(fermicor22_onethird, filek )
          
     ENDIF
     END SUBROUTINE PREQ

     !correlation function         
     Subroutine correlation(gr,filek)
     Use Blockc
     Use Block_obs
     Implicit Real (KIND=8) (A-G,O-Z)
     Implicit Integer (H-N)
     Complex (Kind=8), dimension(:,:,:) :: gr
     Integer :: lp
     Character (32) :: filek
     Real (Kind=8) :: aimj_p(2)
     OPEN (UNIT=20,FILE=filek,STATUS='UNKNOWN', action="write", position="append")
          do nr = 1,LQ
          aimj_p(1) = dble( nlist( nr,1) )
          aimj_p(2) = dble( nlist( nr,2) )
          write(20,*) aimj_p(1), aimj_p(2)
          do no1 = 1,2
               do no2 = 1,2
                    write(20,*) gr(nr,1,1)
               enddo
          enddo
          enddo
     close(20)
     end Subroutine correlation

     !output order parameter
     Subroutine orderparameter(gr1,file)
     Use Blockc
     Use Block_obs
     real(kind=8) :: gr1
     Character (32) :: file
     OPEN (UNIT=20,FILE=file,STATUS='UNKNOWN', action="write", position="append")
          write(20,*) gr1      
     close(20)
     end Subroutine orderparameter

     !structure factor
     Subroutine structurefactor(gr,filek,mom_x,mom_y,no1,no2)
     Use Blockc
     Use Block_obs
     Implicit Integer (H-N)
     Complex (Kind=8), dimension(:,:,:,:) :: gr
     Integer :: lp, no1, no2
     Character (32) :: filek
     Real (Kind=8) :: mom_x,mom_y,xk_p(2), aimj_p(2)
     Integer :: RX_min, RX_max, RY_min, RY_max
     RX_max = NLX-1
     RX_min = 1-NLX
     RY_max = NLY-1
     RY_min = 1-NLY
     gk = cmplx(0.d0,0.d0)
     xk_p(1) = mom_x * b1_p(1) + mom_y * b2_p(1)
     xk_p(2) = mom_x * b1_p(2) + mom_y * b2_p(2)
     do imj_x = RX_min,RX_max
          do imj_y = RY_min,RY_max
               aimj_p(1) = dble(imj_x)* a1_p(1) + dble(imj_y)* a2_p(1)
               aimj_p(2) = dble(imj_x)* a1_p(2) + dble(imj_y)* a2_p(2)
               gk = gk + exp( cmplx( 0.d0, xk_p(1)*aimj_p(1)+xk_p(2)*aimj_p(2) ) ) * gr(imj_x,imj_y,no1,no2)
          enddo
     enddo
     gk = gk/cmplx(LQ,0.d0)
     OPEN (UNIT=20,FILE=filek,STATUS='UNKNOWN', action="write", position="append")
     write(20,*) real(gk)
     close(20)

     !   write(6,*) 'matrix of structurefactor'
     !   write(6,*) filek
     !   WRITE(6, '(A7)',advance='no') "Col/Row" 
     !   DO j = RY_min, RY_max
     !        WRITE(6,'(7X,I5,8X)', advance='no') j
     !   ENDDO
     !   WRITE(6,*) ! \n
     !   DO i = RX_min, RX_max
     !        WRITE(6, '(I5,2X)',advance='no') i 
     !        DO j = RY_min, RY_max
     !             WRITE(6,'( "(", F8.4, ",", F8.4, ") " )', advance='no') REAL(gr(i,j,no1,no1))*1000.d0, AIMAG(gr(i,j,no1,no1))*1000.d0
     !        ENDDO
     !        WRITE(6,*) ! \n
     !   ENDDO

end subroutine structurefactor