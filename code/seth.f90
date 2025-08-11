Subroutine  SetH(HLP2)

   Use Blockc
   Use MyMats
   Implicit Real (KIND=8) (A-G,O-Z)
   Implicit Integer (H-N)
   

   Complex (Kind=8), Dimension(:,:) :: HLP2 
   Complex (Kind=8) :: Z1


!  nearest bond
   IseedHop = 3958195
   HLP2 = CMPLX(0.D0,0.D0)
   do I = 1,LQ
      i_0 = L_Bonds(i,0)
      do nf = 1,Nbond
         i_n = L_Bonds(i,nf)
         random = ranf(IseedHop)
         Z1 = CMPLX(-RT1,0.D0)
         if (Itwist == 1) Z1 = CMPLX(-RT1,0.D0) + TwistX*(random - 0.5d0)
         HLP2(i_0,i_n)  =  Z1
         HLP2(i_n,i_0)  = conjg(Z1)
      enddo
   enddo

   if (Itwist==2) then
      do ii = 1, NDIM
         random = ranf(IseedHop)
         HLP2(ii,ii) = HLP2(ii,ii) + TwistX*(random - 0.5d0)
      enddo
   endif

   if (Itwist==3) then
      ! t2 term twist
      do ix = 1, NLX
         do iy = 1, NLY
            do i_sublattice = 1, Norb
               do i_direction = 1, Nnext
                  i = invlist(ix,iy)
                  ii_0 = L_next(i,i_sublattice,0)
                  ii_n = L_next(i,i_sublattice,i_direction)
                  R2 = - TwistX / dble(LQ) * (-1.0d0)**dble(i_sublattice+i_direction)
                  HLP2(ii_0,ii_n) = CMPLX( 0.d0,  R2)
                  HLP2(ii_n,ii_0) = CMPLX( 0.d0, -R2)
               enddo
            enddo
         enddo
      enddo
   endif

end Subroutine SetH
