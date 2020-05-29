      subroutine pbc(dr,b,hb)
      implicit none
      
      double precision dr,b,hb
      
      if(dr.gt.hb) then
        dr = dr - b
      elseif(dr.le.-hb) then
        dr = dr + b
      endif
      
      return
      end
      
      program sc_overl
      implicit none
      
      integer n,i,j,k,nat_sol
      double precision b(3),hb(3),r_hs,r_hs2,dr,rij2
     &                ,r(50000,3)
      
      open(12,file="conf.xyz")
      read(12,*) n,b(1),b(2),b(3),r_hs,nat_sol
      do i=1,3
         hb(i) = b(i)
      enddo
      r_hs2 = r_hs**2
      do i=1,n
           read(12,*) r(i,1),r(i,2),r(i,3)
      enddo
      close(12)
      
      do i=nat_sol+1,n-1
         do j=i+1,n
	    rij2 = 0.0d0
	    do k=1,3
	       dr = r(j,k)-r(i,k)
	       call pbc(dr,b(k),hb(k))
	       rij2 = rij2 + dr**2
	    enddo
	    if(rij2.le.r_hs2) then
	      write(*,'(a,i,i,d,d,a)') "Overlapping pair found! (i,j,rij,r_hs): ("
     &	                 ,i,j,dsqrt(rij2),r_hs,")"
	    endif
	 enddo
      enddo
      
      stop
      end