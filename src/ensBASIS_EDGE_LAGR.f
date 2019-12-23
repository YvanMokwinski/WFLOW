
      subroutine ensBASIS_EDGE_DRLAGR6(trp_,trc_,n_,p_,poff_,c_,boff_)
      implicit none
      integer n_,poff_,boff_
      character trp_
      character trc_
      double precision c_(boff_,*)
      double precision p_(poff_,*)
      double precision r
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
         else
            r = p_(i,1)
         endif
         if (trc_.eq.'T') then          
            c_(i,1) = 1.25d-2*(r*(r*(r*(r*(4.86d2*r-4.05d2)-1.8d2)
     &           +1.35d2)+8.0d0)-4.0d0)
            c_(i,2) = 5.0d-2*(r*(r*(r*(r*(4.05d2
     &           -7.29d2*r)+5.4d2)-2.7d2)-2.7d1)+9.0d0)
            c_(i,3) = 6.25d-2*(r*(r*(r*(r*(1.458d3*r
     &           -4.05d2)-1.404d3)+3.51d2)+2.16d2)-3.6d1)
            c_(i,4) = 5.0d-1*r*(r**2*(2.52d2-2.43d2*r**2)-4.9d1)
            c_(i,5) = 6.25d-2*(r*(r*(r*(r*(1.458d3*r+
     &           4.05d2)-1.404d3)-3.51d2)+2.16d2)+3.6d1)
            c_(i,6) = 5.0d-2*(r*(r*(r*(r*(-7.29d2*r-4
     &           .05d2)+5.4d2)+2.7d2)-2.7d1)-9.0d0)
            c_(i,7) = 1.25d-2*(r*(r*(r*(r*(4.86d2*r+4
     &           .05d2)-1.8d2)-1.35d2)+8.0d0)+4.0d0)
         else  if (trc_.eq.'N') then 
            c_(1,i) = 1.25d-2*(r*(r*(r*(r*(4.86d2*r-4
     &           .05d2)-1.8d2)+1.35d2)+8.0d0)-4.0d0)
            c_(2,i) = 5.0d-2*(r*(r*(r*(r*(4.05d2-7.29
     &           d2*r)+5.4d2)-2.7d2)-2.7d1)+9.0d0)
            c_(3,i) = 6.25d-2*(r*(r*(r*(r*(1.458d3*r-
     &           4.05d2)-1.404d3)+3.51d2)+2.16d2)-3.6d1)
            c_(4,i) = 5.0d-1*r*(r**2*(2.52d2-2.43d2*r**2)-4.9d1)
            c_(5,i) = 6.25d-2*(r*(r*(r*(r*(1.458d3*r+
     &           4.05d2)-1.404d3)-3.51d2)+2.16d2)+3.6d1)
            c_(6,i) = 5.0d-2*(r*(r*(r*(r*(-7.29d2*r-4
     &           .05d2)+5.4d2)+2.7d2)-2.7d1)-9.0d0)
            c_(7,i) = 1.25d-2*(r*(r*(r*(r*(4.86d2*r+4
     &           .05d2)-1.8d2)-1.35d2)+8.0d0)+4.0d0)
         endif
      end do
      end subroutine   



      subroutine ensBASIS_EDGE_LAGR6(trp_,trc_,n_,p_,poff_,c_,boff_)
      implicit none
      integer n_,poff_,boff_
      character trp_
      character trc_
      double precision c_(boff_,*)
      double precision p_(poff_,*)
      double precision r
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
         else
            r = p_(i,1)
         endif
         if (trc_.eq.'T') then          
            c_(i,1) = 1.25d-2*r*(r*(r*(r*(r*(8.1d1*r-
     &           8.1d1)-4.5d1)+4.5d1)+4.0d0)-4.0d0)
            c_(i,2) = 2.5d-2*r*(r*(r*(r*(r*(1.62d2-2.
     &           43d2*r)+2.7d2)-1.8d2)-2.7d1)+1.8d1)
            c_(i,3) = 6.25d-2*r*(r*(r*(r*(r*(2.43d2*r
     &           -8.1d1)-3.51d2)+1.17d2)+1.08d2)-3.6d1)
            c_(i,4) = 2.5d-1*(r**2*(r**2*(1.26d2-8.1d1*r**2)-4.9d
     &           1)+4.0d0)
            c_(i,5) = 6.25d-2*r*(r*(r*(r*(r*(2.43d2*r
     &           +8.1d1)-3.51d2)-1.17d2)+1.08d2)+3.6d1)
            c_(i,6) = 2.5d-2*r*(r*(r*(r*(r*(-2.43d2*r
     &           -1.62d2)+2.7d2)+1.8d2)-2.7d1)-1.8d1)
            c_(i,7) = 1.25d-2*r*(r*(r*(r*(r*(8.1d1*r+
     &           8.1d1)-4.5d1)-4.5d1)+4.0d0)+4.0d0)
         else  if (trc_.eq.'N') then 
            c_(1,i) = 1.25d-2*r*(r*(r*(r*(r*(8.1d1*r-
     &           8.1d1)-4.5d1)+4.5d1)+4.0d0)-4.0d0)
            c_(2,i) = 2.5d-2*r*(r*(r*(r*(r*(1.62d2-2.
     &           43d2*r)+2.7d2)-1.8d2)-2.7d1)+1.8d1)
            c_(3,i) = 6.25d-2*r*(r*(r*(r*(r*(2.43d2*r
     &           -8.1d1)-3.51d2)+1.17d2)+1.08d2)-3.6d1)
            c_(4,i) = 2.5d-1*(r**2*(r**2*(1.26d2-8.1d1*r**2)-4.9d
     &           1)+4.0d0)
            c_(5,i) = 6.25d-2*r*(r*(r*(r*(r*(2.43d2*r
     &           +8.1d1)-3.51d2)-1.17d2)+1.08d2)+3.6d1)
            c_(6,i) = 2.5d-2*r*(r*(r*(r*(r*(-2.43d2*r
     &           -1.62d2)+2.7d2)+1.8d2)-2.7d1)-1.8d1)
            c_(7,i) = 1.25d-2*r*(r*(r*(r*(r*(8.1d1*r+
     &           8.1d1)-4.5d1)-4.5d1)+4.0d0)+4.0d0)
         endif
      end do
      end subroutine   



      subroutine ensBASIS_EDGE_DRLAGR5(trp_,trc_,n_,p_,poff_,c_,boff_)
      implicit none
      integer n_,poff_,boff_
      character trp_
      character trc_
      double precision c_(boff_,*)
      double precision p_(poff_,*)
      double precision r
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
         else
            r = p_(i,1)
         endif
         if (trc_.eq.'T') then          
            c_(i,1) = 1.302083333333333d-3*(r*(r*(r*(2.5d3-
     &           3.125d3*r)+7.5d2)-5.0d2)-9.0d0)
            c_(i,2) = 1.302083333333333d-3*(r*(r*(r*(1.5625
     &           d4*r-7.5d3)-9.75d3)+3.9d3)+1.25d2)
            c_(i,3) = 2.604166666666667d-3*(r*(r*(r*(2.5d3-
     &           1.5625d4*r)+1.275d4)-1.7d3)-1.125d3)
            c_(i,4) = 2.604166666666667d-3*(r*(r*(r*(1.5625
     &           d4*r+2.5d3)-1.275d4)-1.7d3)+1.125d3)
            c_(i,5) = 1.302083333333333d-3*(r*(r*(r*(-1.562
     &           5d4*r-7.5d3)+9.75d3)+3.9d3)-1.25d2)
            c_(i,6) = 1.302083333333333d-3*(r*(r*(r*(3.125d
     &           3*r+2.5d3)-7.5d2)-5.0d2)+9.0d0)
         else  if (trc_.eq.'N') then 
            c_(1,i) = 1.302083333333333d-3*(r*(r*(r*(2.5d3-
     &           3.125d3*r)+7.5d2)-5.0d2)-9.0d0)
            c_(2,i) = 1.302083333333333d-3*(r*(r*(r*(1.5625
     &           d4*r-7.5d3)-9.75d3)+3.9d3)+1.25d2)
            c_(3,i) = 2.604166666666667d-3*(r*(r*(r*(2.5d3-
     &           1.5625d4*r)+1.275d4)-1.7d3)-1.125d3)
            c_(4,i) = 2.604166666666667d-3*(r*(r*(r*(1.5625
     &           d4*r+2.5d3)-1.275d4)-1.7d3)+1.125d3)
            c_(5,i) = 1.302083333333333d-3*(r*(r*(r*(-1.562
     &           5d4*r-7.5d3)+9.75d3)+3.9d3)-1.25d2)
            c_(6,i) = 1.302083333333333d-3*(r*(r*(r*(3.125d
     &           3*r+2.5d3)-7.5d2)-5.0d2)+9.0d0)
         endif
      end do
      end subroutine   



      subroutine ensBASIS_EDGE_LAGR5(trp_,trc_,n_,p_,poff_,c_,boff_)
      implicit none
      integer n_,poff_,boff_
      character trp_
      character trc_
      double precision c_(boff_,*)
      double precision p_(poff_,*)
      double precision r
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
         else
            r = p_(i,1)
         endif
         if (trc_.eq.'T') then          
            c_(i,1) = 1.302083333333333d-3*(r*(r*(r*(
     &           r*(6.25d2-6.25d2*r)+2.5d2)-2.5d2)-9.0d0)+9.0d0)
            c_(i,2) = 1.302083333333333d-3*(r*(r*(r*(
     &           r*(3.125d3*r-1.875d3)-3.25d3)+1.95d3)+1.25d2)-7.5d1)
            c_(i,3) = 2.604166666666667d-3*(r*(r*(r*(
     &           r*(6.25d2-3.125d3*r)+4.25d3)-8.5d2)-1.125d3)+2.25d2)
            c_(i,4) = 2.604166666666667d-3*(r*(r*(r*(
     &           r*(3.125d3*r+6.25d2)-4.25d3)-8.5d2)+1.125d3)+2.25d2)
            c_(i,5) = 1.302083333333333d-3*(r*(r*(r*(
     &           r*(-3.125d3*r-1.875d3)+3.25d3)+1.95d3)-1.25d2)-7.5d1)
            c_(i,6) = 1.302083333333333d-3*(r*(r*(r*(
     &           r*(6.25d2*r+6.25d2)-2.5d2)-2.5d2)+9.0d0)+9.0d0)
         else  if (trc_.eq.'N') then 
            c_(1,i) = 1.302083333333333d-3*(r*(r*(r*(
     &           r*(6.25d2-6.25d2*r)+2.5d2)-2.5d2)-9.0d0)+9.0d0)
            c_(2,i) = 1.302083333333333d-3*(r*(r*(r*(
     &           r*(3.125d3*r-1.875d3)-3.25d3)+1.95d3)+1.25d2)-7.5d1)
            c_(3,i) = 2.604166666666667d-3*(r*(r*(r*(
     &           r*(6.25d2-3.125d3*r)+4.25d3)-8.5d2)-1.125d3)+2.25d2)
            c_(4,i) = 2.604166666666667d-3*(r*(r*(r*(
     &           r*(3.125d3*r+6.25d2)-4.25d3)-8.5d2)+1.125d3)+2.25d2)
            c_(5,i) = 1.302083333333333d-3*(r*(r*(r*(
     &           r*(-3.125d3*r-1.875d3)+3.25d3)+1.95d3)-1.25d2)-7.5d1)
            c_(6,i) = 1.302083333333333d-3*(r*(r*(r*(
     &           r*(6.25d2*r+6.25d2)-2.5d2)-2.5d2)+9.0d0)+9.0d0)
         endif
      end do
      end subroutine   

      subroutine ensBASIS_EDGE_DRLAGR4(trp_,trc_,n_,p_,poff_,c_,boff_)
      implicit none
      integer n_,poff_,boff_
      character trp_
      character trc_
      double precision c_(boff_,*)
      double precision p_(poff_,*)
      double precision r
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
         else
            r = p_(i,1)
         endif
         if (trc_.eq.'T') then          
            c_(i,1) = 1.666666666666667d-1*(r*(r*(1.6d1*r-1.2d1)-
     &           2.0d0)+1.0d0)
            c_(i,2) = 3.333333333333333d-1*(r*(r*(1.2d1-3.2d1*r)+
     &           1.6d1)-4.0d0)
            c_(i,3) = r*(1.6d1*r**2-1.0d1)
            c_(i,4) = 3.333333333333333d-1*(r*(r*(-3.2d1*r-1.2d1)
     &           +1.6d1)+4.0d0)
            c_(i,5) = 1.666666666666667d-1*(r*(r*(1.6d1*r+1.2d1)-
     &           2.0d0)-1.0d0)
         else  if (trc_.eq.'N') then 
            c_(1,i) = 1.666666666666667d-1*(r*(r*(1.6d1*r-1.2d1)-
     &           2.0d0)+1.0d0)
            c_(2,i) = 3.333333333333333d-1*(r*(r*(1.2d1-3.2d1*r)+
     &           1.6d1)-4.0d0)
            c_(3,i) = r*(1.6d1*r**2-1.0d1)
            c_(4,i) = 3.333333333333333d-1*(r*(r*(-3.2d1*r-1.2d1)
     &           +1.6d1)+4.0d0)
            c_(5,i) = 1.666666666666667d-1*(r*(r*(1.6d1*r+1.2d1)-
     &           2.0d0)-1.0d0)
         endif
      end do
      end subroutine   



      subroutine ensBASIS_EDGE_LAGR4(trp_,trc_,n_,p_,poff_,c_,boff_)
      implicit none
      integer n_,poff_,boff_
      character trp_
      character trc_
      double precision c_(boff_,*)
      double precision p_(poff_,*)
      double precision r
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
         else
            r = p_(i,1)
         endif
         if (trc_.eq.'T') then          
            c_(i,1) = 1.666666666666667d-1*r*(r*(r*(4.0d0*r-4.0d0
     &           )-1.0d0)+1.0d0)
            c_(i,2) = 3.333333333333333d-1*r*(r*(r*(4.0d0-8.0d0*r
     &           )+8.0d0)-4.0d0)
            c_(i,3) = r**2*(4.0d0*r**2-5.0d0)+1.0d0
            c_(i,4) = 3.333333333333333d-1*r*(r*(r*(-8.0d0*r-4.0d
     &           0)+8.0d0)+4.0d0)
            c_(i,5) = 1.666666666666667d-1*r*(r*(r*(4.0d0*r+4.0d0
     &           )-1.0d0)-1.0d0)
         else  if (trc_.eq.'N') then 
            c_(1,i) = 1.666666666666667d-1*r*(r*(r*(4.0d0*r-4.0d0
     &           )-1.0d0)+1.0d0)
            c_(2,i) = 3.333333333333333d-1*r*(r*(r*(4.0d0-8.0d0*r
     &           )+8.0d0)-4.0d0)
            c_(3,i) = r**2*(4.0d0*r**2-5.0d0)+1.0d0
            c_(4,i) = 3.333333333333333d-1*r*(r*(r*(-8.0d0*r-4.0d
     &           0)+8.0d0)+4.0d0)
            c_(5,i) = 1.666666666666667d-1*r*(r*(r*(4.0d0*r+4.0d0
     &           )-1.0d0)-1.0d0)
         endif
      end do
      end subroutine   


      subroutine ensBASIS_EDGE_DRLAGR3(trp_,trc_,n_,p_,poff_,c_,boff_)
      implicit none
      integer n_,poff_,boff_
      character trp_
      character trc_
      double precision c_(boff_,*)
      double precision p_(poff_,*)
      double precision r
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
         else
            r = p_(i,1)
         endif
         if (trc_.eq.'T') then          

            c_(i,1) = 5.0d-2*(r*(2.4d1-2.7d1*r)-1.0d0)
            c_(i,2) = 1.25d-1*(r*(2.7d1*r-1.2d1)-9.0d0)
            c_(i,3) = 2.0d-1*(r*(-2.7d1*r-6.0d0)+9.0d0)
            c_(i,4) = 1.25d-1*(r*(2.7d1*r+1.2d1)-5.0d0)
         else  if (trc_.eq.'N') then 

            c_(1,i) = 5.0d-2*(r*(2.4d1-2.7d1*r)-1.0d0)
            c_(2,i) = 1.25d-1*(r*(2.7d1*r-1.2d1)-9.0d0)
            c_(3,i) = 2.0d-1*(r*(-2.7d1*r-6.0d0)+9.0d0)
            c_(4,i) = 1.25d-1*(r*(2.7d1*r+1.2d1)-5.0d0)
         endif
      end do
      end subroutine   



      subroutine ensBASIS_EDGE_LAGR3(trp_,trc_,n_,p_,poff_,c_,boff_)
      implicit none
      integer n_,poff_,boff_
      character trp_
      character trc_
      double precision c_(boff_,*)
      double precision p_(poff_,*)
      double precision r
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
         else
            r = p_(i,1)
         endif
         if (trc_.eq.'T') then          
            c_(i,1) = 5.0d-2*(r*(r*(1.2d1-9.0d0*r)-1.0d0)-2.0d0)
            c_(i,2) = 1.25d-1*(r*(r*(9.0d0*r-6.0d0)-9.0d0)+6.0d0)
     &           
            c_(i,3) = 2.0d-1*(r*(r*(-9.0d0*r-3.0d0)+9.0d0)+3.0d0)
     &           
            c_(i,4) = 1.25d-1*(r*(r*(9.0d0*r+6.0d0)-5.0d0)-2.0d0)
         else  if (trc_.eq.'N') then 
            c_(1,i) = 5.0d-2*(r*(r*(1.2d1-9.0d0*r)-1.0d0)-2.0d0)
            c_(2,i) = 1.25d-1*(r*(r*(9.0d0*r-6.0d0)-9.0d0)+6.0d0)
     &           
            c_(3,i) = 2.0d-1*(r*(r*(-9.0d0*r-3.0d0)+9.0d0)+3.0d0)
     &           
            c_(4,i) = 1.25d-1*(r*(r*(9.0d0*r+6.0d0)-5.0d0)-2.0d0)
         endif
      end do
      end subroutine   







      subroutine ensBASIS_EDGE_DRLAGR2(trp_,trc_,n_,p_,poff_,c_,boff_)
      implicit none
      integer n_,poff_,boff_
      character trp_
      character trc_
      double precision c_(boff_,*)
      double precision p_(poff_,*)
      double precision r
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
         else
            r = p_(i,1)
         endif
         if (trc_.eq.'T') then          
            c_(i,1) = 5.0d-1*(2.0d0*r-1.0d0)
            c_(i,2) = -2.0d0*r
            c_(i,3) = 5.0d-1*(2.0d0*r+1.0d0)
         else  if (trc_.eq.'N') then 
            c_(1,i) = 5.0d-1*(2.0d0*r-1.0d0)
            c_(2,i) = -2.0d0*r
            c_(3,i) = 5.0d-1*(2.0d0*r+1.0d0)
         endif
      end do
      end subroutine   

      subroutine ensBASIS_EDGE_LAGR2(trp_,trc_,n_,p_,poff_,c_,boff_)
      implicit none
      integer n_,poff_,boff_
      character trp_
      character trc_
      double precision c_(boff_,*)
      double precision p_(poff_,*)
      double precision r
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
         else
            r = p_(i,1)
         endif
         if (trc_.eq.'T') then          
            c_(i,1) = 5.0d-1*r*(r-1.0d0)
            c_(i,2) = 1.0d0-r**2
            c_(i,3) = 5.0d-1*r*(r+1.0d0)
         else  if (trc_.eq.'N') then 
            c_(1,i) = 5.0d-1*r*(r-1.0d0)
            c_(2,i) = 1.0d0-r**2
            c_(3,i) = 5.0d-1*r*(r+1.0d0)
         endif
      end do
      end subroutine   



      subroutine ensBASIS_EDGE_DRLAGR1(trp_,trc_,n_,p_,poff_,c_,boff_)
      implicit none
      integer n_,poff_,boff_
      character trp_
      character trc_
      double precision c_(boff_,*)
      double precision p_(poff_,*)
      integer i
      if (trp_.eq.'N') then
         c_(1,1) = p_(1,1)
      endif
      if (trc_.eq.'T') then          
         do i=1,n_
            c_(i,1) = -0.5d0
            c_(i,2) = 0.5d0
         end do
      else  if (trc_.eq.'N') then
         do i=1,n_
            c_(1,i) = -0.5d0
            c_(2,i) = 0.5d0
         end do
      endif
      end subroutine   

      subroutine ensBASIS_EDGE_LAGR1(trp_,trc_,n_,p_,poff_,c_,boff_)
      implicit none
      integer n_,poff_,boff_
      character trp_
      character trc_
      double precision c_(boff_,*)
      double precision p_(poff_,*)
      double precision r
      integer i
      
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
         else 
            r = p_(i,1)
         endif
         if (trc_.eq.'T') then          
            c_(i,1) = (1.0d0-r)*0.5d0
            c_(i,2) = (1.0d0+r)*0.5d0
         else  if (trc_.eq.'N') then 
            c_(1,i) = (1.0d0-r)*0.5d0
            c_(2,i) = (1.0d0+r)*0.5d0
         endif
      end do
      end subroutine   



      subroutine ensBASIS_EDGE_DRLAGR0(trp_,trc_,n_,p_,poff_,c_,boff_)
      implicit none
      integer n_,poff_,boff_
      character trp_
      character trc_
      double precision c_(boff_,*)
      double precision p_(poff_,*)
      integer i
      if (trp_.eq.'N') then
         c_(1,1) = p_(1,1)
      endif
      
      if (trc_.eq.'T') then          
         do i=1,n_
            c_(i,1) = 0.d0
         end do
      else
         do i=1,n_
            c_(1,i) = 0.d0
         end do
      endif
      end subroutine   



      subroutine ensBASIS_EDGE_LAGR0(trp_,trc_,n_,p_,poff_,c_,boff_)
      implicit none
      integer n_,poff_,boff_
      character trp_
      character trc_
      double precision c_(boff_,*)
      double precision p_(poff_,*)
      integer i
      if (trp_.eq.'N') then
         c_(1,1) = p_(1,1)
      endif
      if (trc_.eq.'T') then          
         do i=1,n_
            c_(i,1) = 1.d0
         end do
      else
         do i=1,n_
            c_(1,i) = 1.d0
         end do
      endif
      end subroutine   


