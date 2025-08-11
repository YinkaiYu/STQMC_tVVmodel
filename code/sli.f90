    SUBROUTINE SLI
        Use Blockc
        Implicit Real (KIND=8) (A-G,O-Z)
        Implicit Integer (H-N)
     
        LIST = 0; INVLIST = 0; NNLIST = 0

! nlist: Index_site -> (nx,ny)
        NCOUNT = 0
	    do NX = 1,NLX
	        do NY = 1,NLY
	            NCOUNT = NCOUNT + 1
	            LIST(NCOUNT,1) = NX
	            LIST(NCOUNT,2) = NY
	            INVLIST(NX,NY) = NCOUNT
            enddo
        enddo

! nlist: Index_dimen -> (nx,ny,n_orbital, n_spin)
	    NCOUNT = 0	
        do nx = 1, NLX
            do ny = 1, NLY
                do NO = 1, Norb
                    do NS = 1 ,nspin
                        ncount = ncount + 1
                        nlist(NCOUNT,1) = NX
                        nlist(NCOUNT,2) = NY
                        nlist(NCOUNT,3) = NO
                        nlist(NCOUNT,4) = NS
                        invnlist(NX,NY,NO,NS) = NCOUNT
                    enddo
                enddo
            enddo
        enddo
	
!define the nearest neighbor bonds
         do ly = 1, NLY
            do lx = 1, NLX               
                ll = invlist(lx,ly)
                L_bonds(ll,0) = invnlist(lx,ly,1,1)
                L_bonds(ll,1) = invnlist(lx,ly,2,1)
                lln1 = invnlist(npbcx(lx+1),ly,2,1)
                L_bonds(ll,2) = lln1
                lln2 = invnlist(lx,npbcy(ly+1),2,1)
                L_bonds(ll,3) = lln2
            enddo
         enddo
	
!define the next nearest neighbor interaction
        do no = 1, Norb
            do ly = 1, NLY
                do lx = 1, NLX               
                    ll = invlist(lx,ly)
                    L_next(ll,no,0) = invnlist(lx,ly,no,1)
                    L_next(ll,no,1) = invnlist(npbcx(lx+1),ly,no,1)
                    L_next(ll,no,2) = invnlist(lx,npbcy(ly+1),no,1)
                    L_next(ll,no,3) = invnlist(npbcx(lx-1),npbcy(ly+1),no,1)
                enddo
            enddo
        enddo

! define the unite matrix
        ZKRON = DCMPLX(0.D0,0.D0)
	    do I = 1,NDIM
	        ZKRON(I,I) = DCMPLX(1.D0,0.D0)
	    enddo
        
	    RETURN
    END SUBROUTINE SLI
    
    integer FUNCTION Iscalar(vec1,vec2)
        Use Blockc
        Iscalar = vec1(1)*vec2(1) + vec1(2)*vec2(2)       
        return    
    end function Iscalar
    
      
    INTEGER FUNCTION NPBCX(NR)
        use blockc
        NPBCX = NR
        IF (NR.GT.NLX) NPBCX = NR - NLX
        IF (NR.LT.1) NPBCX = NR + NLX
        RETURN
    END FUNCTION NPBCX
    
    INTEGER FUNCTION NPBCY(NR)
        use blockc
        NPBCY = NR
        IF (NR.GT.NLY) NPBCY = NR - NLY
        IF (NR.LT.1) NPBCY = NR + NLY
        RETURN
    END FUNCTION NPBCY
    
     
    INTEGER FUNCTION NPBC_tau(NR)
        use blockc
        NPBC_tau = NR
        IF (NR.GT.LTROT) NPBC_tau = NR - LTROT
        IF (NR.LT.1) NPBC_tau = NR + LTROT
        RETURN  
    END FUNCTION NPBC_tau


    