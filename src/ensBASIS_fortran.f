
      subroutine ensBASIS_TRIA_DSLAGR3(trp_,trb_,n_,p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,*)
      double precision p_(poff_,*)
      double precision r,s
      integer i
      character trp_
      character trb_
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         
         if (trb_.eq.'T') then 
         b_(i,1) = 5.0d-1*(s*(-2.7d1*s-5.4d1*r+3.6d1)
     &        +r*(3.6d1-2.7d1*r)-1.1d1)
         b_(i,2) = 0.0d0
         b_(i,3) = 5.0d-1*(s*(2.7d1*s-1.8d1)+2.0d0)
         b_(i,4) = 5.0d-1*(5.4d1*r*s+r*(5.4d1*r-4.5d1))
         b_(i,5) = 5.0d-1*r*(9.0d0-2.7d1*r)
         b_(i,6) = 5.0d-1*r*(2.7d1*r-9.0d0)
         b_(i,7) = 5.0d-1*(5.4d1*r*s-9.0d0*r)
         b_(i,8) = 5.0d-1*(s*(-8.1d1*s-5.4d1*r+7.2d1)
     &        +9.0d0*r-9.0d0)
         b_(i,9) = 5.0d-1*(s*(8.1d1*s+1.08d2*r-9.0d1)
     &        +r*(2.7d1*r-4.5d1)+1.8d1)
         b_(i,10) = r*(2.7d1-2.7d1*r)-5.4d1*r*s         
         else 
         b_(1,i) = 5.0d-1*(s*(-2.7d1*s-5.4d1*r+3.6d1)
     &        +r*(3.6d1-2.7d1*r)-1.1d1)
         b_(2,i) = 0.0d0
         b_(3,i) = 5.0d-1*(s*(2.7d1*s-1.8d1)+2.0d0)
         b_(4,i) = 5.0d-1*(5.4d1*r*s+r*(5.4d1*r-4.5d1))
         b_(5,i) = 5.0d-1*r*(9.0d0-2.7d1*r)
         b_(6,i) = 5.0d-1*r*(2.7d1*r-9.0d0)
         b_(7,i) = 5.0d-1*(5.4d1*r*s-9.0d0*r)
         b_(8,i) = 5.0d-1*(s*(-8.1d1*s-5.4d1*r+7.2d1)
     &        +9.0d0*r-9.0d0)
         b_(9,i) = 5.0d-1*(s*(8.1d1*s+1.08d2*r-9.0d1)
     &        +r*(2.7d1*r-4.5d1)+1.8d1)
         b_(10,i) = r*(2.7d1-2.7d1*r)-5.4d1*r*s         
         endif
      end do
      end subroutine

      subroutine ensBASIS_TRIA_DSLAGRBUBBLE2(trp_,trb_,
     &     n_,p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,*)
      double precision p_(poff_,*)
      character trp_
      character trb_
      double precision r,s
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T') then 
         b_(i,1) = 4.0d0*s+4.0d0*r-3.0d0
         b_(i,2) = 0.0d0
         b_(i,3) = 4.0d0*s-1.0d0
         b_(i,4) = -4.0d0*r
         b_(i,5) = 4.0d0*r
         b_(i,6) = -8.0d0*s-4.0d0*r+4.0d0
         b_(i,7) = r*(2.7d1-2.7d1*r)-5.4d1*r*s         
         else 
         b_(1,i) = 4.0d0*s+4.0d0*r-3.0d0
         b_(2,i) = 0.0d0
         b_(3,i) = 4.0d0*s-1.0d0
         b_(4,i) = -4.0d0*r
         b_(5,i) = 4.0d0*r
         b_(6,i) = -8.0d0*s-4.0d0*r+4.0d0
         b_(7,i) = r*(2.7d1-2.7d1*r)-5.4d1*r*s         
         endif
      end do
      end subroutine

      subroutine ensBASIS_TRIA_DSLAGR2(trp_,trb_,n_,p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      character trp_
      character trb_
      double precision b_(boff_,*)
      double precision p_(poff_,*)
      double precision r,s
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T') then 
         b_(i,1) = 4.0d0*s+4.0d0*r-3.0d0
         b_(i,2) = 0.0d0
         b_(i,3) = 4.0d0*s-1.0d0
         b_(i,4) = -4.0d0*r
         b_(i,5) = 4.0d0*r
         b_(i,6) = -8.0d0*s-4.0d0*r+4.0d0
         else 
         b_(1,i) = 4.0d0*s+4.0d0*r-3.0d0
         b_(2,i) = 0.0d0
         b_(3,i) = 4.0d0*s-1.0d0
         b_(4,i) = -4.0d0*r
         b_(5,i) = 4.0d0*r
         b_(6,i) = -8.0d0*s-4.0d0*r+4.0d0
         endif
      end do
      end subroutine   



      subroutine ensBASIS_TRIA_DSLAGRBUBBLE1(trp_,
     &     trb_,n_,p_,poff_,b_,boff_)
      implicit none
      character trp_
      character trb_
      integer n_,poff_,boff_
      double precision b_(boff_,*)
      double precision p_(poff_,*)
      double precision r,s
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         
         if (trb_.eq.'T') then 
         b_(i,1) = -1.0d0                                                         
         b_(i,2) = 0.0d0                                                         
         b_(i,3) = 1.0d0                                                         
         b_(i,4) = r*(2.7d1-2.7d1*r)-5.4d1*r*s         
         else 
         b_(1,i) = -1.0d0                                                         
         b_(2,i) = 0.0d0                                                         
         b_(3,i) = 1.0d0                                                         
         b_(4,i) = r*(2.7d1-2.7d1*r)-5.4d1*r*s         
         endif
      end do
      end subroutine   



      subroutine ensBASIS_TRIA_DSLAGR1(trp_,trb_,n_,p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      character trp_
      character trb_
      double precision b_(boff_,*)
      double precision p_(poff_,*)

      double precision r,s
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T') then 
         
         b_(i,1) = -1.0d0                                                         
         b_(i,2) = 0.0d0                                                         
         b_(i,3) = 1.0d0                                                         
         else          
         b_(1,i) = -1.0d0                                                         
         b_(2,i) = 0.0d0                                                         
         b_(3,i) = 1.0d0                                                         
         endif
      end do
      end subroutine   


      subroutine ensBASIS_TRIA_DSLAGR0(trp_,trb_,n_,p_,poff_,b_,boff_)
      implicit none
      character trp_
      character trb_
      integer n_,poff_,boff_
      double precision b_(boff_,*)
      double precision p_(poff_,*)
      double precision r,s
      integer i
      do i=1,n_
         if (trb_.eq.'T') then 
         b_(i,1) = 0.0d0                                                         
         else 
         b_(1,i) = 0.0d0                                                         
         endif
      end do
      end subroutine   


      subroutine ensBASIS_TRIA_DRLAGR3(trp_,trb_,n_,p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,*)
      double precision p_(poff_,*)
      double precision r,s
      character trp_
      character trb_
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T') then 
            b_(i,1) = 5.0d-1*(s*(-2.7d1*s-5.4d1*r+3.6d1)
     &           +r*(3.6d1-2.7d1*r)-1.1d1)
            b_(i,2) = 5.0d-1*(r*(2.7d1*r-1.8d1)+2.0d0)
            b_(i,3) = 0.0d0
            b_(i,4) = 5.0d-1*(s*(2.7d1*s+1.08d2*r-4.5d1)
     &           +r*(8.1d1*r-9.0d1)+1.8d1)
            b_(i,5) = 5.0d-1*((9.0d0-5.4d1*r)*s+r*(7.2d1-8.1d1*r)-9.0d0)
            b_(i,6) = 5.0d-1*(5.4d1*r-9.0d0)*s
            b_(i,7) = 5.0d-1*s*(2.7d1*s-9.0d0)
            b_(i,8) = 5.0d-1*s*(9.0d0-2.7d1*s)
            b_(i,9) = 5.0d-1*s*(5.4d1*s+5.4d1*r-4.5d1)
            b_(i,10) = s*(-2.7d1*s-5.4d1*r+2.7d1)
         else if (trb_.eq.'N') then 
            b_(1,i) = 5.0d-1*(s*(-2.7d1*s-5.4d1*r+3.6d1)
     &           +r*(3.6d1-2.7d1*r)-1.1d1)
            b_(2,i) = 5.0d-1*(r*(2.7d1*r-1.8d1)+2.0d0)
            b_(3,i) = 0.0d0
            b_(4,i) = 5.0d-1*(s*(2.7d1*s+1.08d2*r-4.5d1)
     &           +r*(8.1d1*r-9.0d1)+1.8d1)
            b_(5,i) = 5.0d-1*((9.0d0-5.4d1*r)*s+r*(7.2d1-8.1d1*r)-9.0d0)
            b_(6,i) = 5.0d-1*(5.4d1*r-9.0d0)*s
            b_(7,i) = 5.0d-1*s*(2.7d1*s-9.0d0)
            b_(8,i) = 5.0d-1*s*(9.0d0-2.7d1*s)
            b_(9,i) = 5.0d-1*s*(5.4d1*s+5.4d1*r-4.5d1)
            b_(10,i) = s*(-2.7d1*s-5.4d1*r+2.7d1)
         endif
      end do
      end subroutine


      subroutine ensBASIS_TRIA_DRLAGR2(trp_,trb_,n_,p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,*)
      double precision p_(poff_,*)
      double precision r,s
      character trp_
      character trb_
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T') then 
            b_(i,1) = 4.0d0*s+4.0d0*r-3.0d0
            b_(i,2) = 4.0d0*r-1.0d0
            b_(i,3) = 0.0d0
            b_(i,4) = -4.0d0*s-8.0d0*r+4.0d0
            b_(i,5) = 4.0d0*s
            b_(i,6) = -4.0d0*s
         else if (trb_.eq.'N') then 
            b_(1,i) = 4.0d0*s+4.0d0*r-3.0d0
            b_(2,i) = 4.0d0*r-1.0d0
            b_(3,i) = 0.0d0
            b_(4,i) = -4.0d0*s-8.0d0*r+4.0d0
            b_(5,i) = 4.0d0*s
            b_(6,i) = -4.0d0*s
         endif
      end do
      end subroutine   


      subroutine ensBASIS_TRIA_DRLAGR1(trp_,trb_,n_,p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,*)
      double precision p_(poff_,*)

      character trp_
      character trb_
      double precision r,s
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T') then 
            b_(i,1) = -1.0d0                                                         
            b_(i,2) = 1.0d0                                                         
            b_(i,3) = 0.0d0                                                         
         else
            b_(1,i) = -1.0d0                                                         
            b_(2,i) = 1.0d0                                                         
            b_(3,i) = 0.0d0                                                         
         endif         
      end do
      end subroutine   

      subroutine ensBASIS_TRIA_DRLAGRBUBBLE2(trp_,trb_,n_,p_,
     &     poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,*)
      double precision p_(poff_,*)
      character trp_
      character trb_
      double precision r,s
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T') then 
            b_(i,1) = 4.0d0*s+4.0d0*r-3.0d0
            b_(i,2) = 4.0d0*r-1.0d0
            b_(i,3) = 0.0d0
            b_(i,4) = -4.0d0*s-8.0d0*r+4.0d0
            b_(i,5) = 4.0d0*s
            b_(i,6) = -4.0d0*s
            b_(i,7) = s*(-2.7d1*s-5.4d1*r+2.7d1)
         else if (trb_.eq.'N') then 
            b_(1,i) = 4.0d0*s+4.0d0*r-3.0d0
            b_(2,i) = 4.0d0*r-1.0d0
            b_(3,i) = 0.0d0
            b_(4,i) = -4.0d0*s-8.0d0*r+4.0d0
            b_(5,i) = 4.0d0*s
            b_(6,i) = -4.0d0*s
            b_(7,i) = s*(-2.7d1*s-5.4d1*r+2.7d1)
         endif
      end do
      end subroutine   



      subroutine ensBASIS_TRIA_DRLAGRBUBBLE1(trp_,trb_,n_,p_,poff_,
     & b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,*)
      double precision p_(poff_,*)
      character trp_
      character trb_

      double precision r,s
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T') then 
            b_(i,1) = -1.0d0                                                         
            b_(i,2) = 1.0d0                                                         
            b_(i,3) = 0.0d0                                                         
            b_(i,4) = s*(-2.7d1*s-5.4d1*r+2.7d1)
         else if (trb_.eq.'N') then 
            b_(1,i) = -1.0d0                                                         
            b_(2,i) = 1.0d0                                                         
            b_(3,i) = 0.0d0                                                         
            b_(4,i) = s*(-2.7d1*s-5.4d1*r+2.7d1)
         endif
      end do
      end subroutine   


      subroutine ensBASIS_TRIA_DRLAGR0(trp_,trb_,n_,p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,*)
      double precision p_(poff_,*)
      character trp_
      character trb_
      double precision r,s
      integer i
      do i=1,n_
         if (trb_.eq.'T') then 
            b_(i,1) = 0.0d0                                                         
         else
            b_(1,i) = 0.0d0                                                         
         endif
      end do
      end subroutine   



      subroutine ensBASIS_TRIA_LAGR3(trp_,trb_,n_,p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,*)
      double precision p_(poff_,*)
      double precision r,s
      character trp_
      character trb_
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T')then
            b_(i,1) = 5.0d-1*(s*(s*(-9.0d0*s-2.7d1*r+1.8d1)
     &           +r*(3.6d1-2.7d1*r)-1.1d1)+r*(r*(1.8d1-9.0d0*r)-1.1d1)
     &           +2.0d0)
            b_(i,2) = 5.0d-1*r*(r*(9.0d0*r-9.0d0)+2.0d0)
            b_(i,3) = 5.0d-1*s*(s*(9.0d0*s-9.0d0)+2.0d0)
            b_(i,4) = 5.0d-1*(s*(2.7d1*r*s+r*(5.4d1*r-4.5d1))
     &           +r*(r*(2.7d1*r-4.5d1)+1.8d1))
            b_(i,5) = 5.0d-1*(r*(9.0d0-2.7d1*r)*s
     &           +r*(r*(3.6d1-2.7d1*r)-9.0d0))
            b_(i,6) = 5.0d-1*r*(2.7d1*r-9.0d0)*s
            b_(i,7) = 5.0d-1*s*(2.7d1*r*s-9.0d0*r)
            b_(i,8) = 5.0d-1*s*(s*(-2.7d1*s-2.7d1*r+3.6d1)
     &           +9.0d0*r-9.0d0)
            b_(i,9) = 5.0d-1*s*(s*(2.7d1*s+5.4d1*r-4.5d1)
     &           +r*(2.7d1*r-4.5d1)+1.8d1)
            b_(i,10) = s*(r*(2.7d1-2.7d1*r)-2.7d1*r*s)
         else if (trb_.eq.'N')then
            b_(1,i) = 5.0d-1*(s*(s*(-9.0d0*s-2.7d1*r+1.8d1)
     &           +r*(3.6d1-2.7d1*r)-1.1d1)+r*(r*(1.8d1-9.0d0*r)-1.1d1)
     &           +2.0d0)
            b_(2,i) = 5.0d-1*r*(r*(9.0d0*r-9.0d0)+2.0d0)
            b_(3,i) = 5.0d-1*s*(s*(9.0d0*s-9.0d0)+2.0d0)
            b_(4,i) = 5.0d-1*(s*(2.7d1*r*s+r*(5.4d1*r-4.5d1))
     &           +r*(r*(2.7d1*r-4.5d1)+1.8d1))
            b_(5,i) = 5.0d-1*(r*(9.0d0-2.7d1*r)*s
     &           +r*(r*(3.6d1-2.7d1*r)-9.0d0))
            b_(6,i) = 5.0d-1*r*(2.7d1*r-9.0d0)*s
            b_(7,i) = 5.0d-1*s*(2.7d1*r*s-9.0d0*r)
            b_(8,i) = 5.0d-1*s*
     &           (s*(-2.7d1*s-2.7d1*r+3.6d1)+9.0d0*r-9.0d0)
            b_(9,i) = 5.0d-1*s*
     &           (s*(2.7d1*s+5.4d1*r-4.5d1)+r*(2.7d1*r-4.5d1)
     &           +1.8d1)
            b_(10,i) = s*(r*(2.7d1-2.7d1*r)-2.7d1*r*s)
         endif
      end do
      end subroutine



      subroutine ensBASIS_TRIA_LAGR2(trp_,trb_,n_,p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,*)
      double precision p_(poff_,*)
      double precision r,s
      character trp_
      character trb_
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T')then
            b_(i,1) = s*(2.0d0*s+4.0d0*r-3.0d0)+
     &           r*(2.0d0*r-3.0d0)+1.0d0
            b_(i,2) = r*(2.0d0*r-1.0d0)
            b_(i,3) = s*(2.0d0*s-1.0d0)
            b_(i,4) = r*(4.0d0-4.0d0*r)-4.0d0*r*s
            b_(i,5) = 4.0d0*r*s
            b_(i,6) = s*(-4.0d0*s-4.0d0*r+4.0d0)
         else if (trb_.eq.'N')then
            b_(1,i) = s*(2.0d0*s+4.0d0*r-3.0d0)+
     &           r*(2.0d0*r-3.0d0)+1.0d0
            b_(2,i) = r*(2.0d0*r-1.0d0)
            b_(3,i) = s*(2.0d0*s-1.0d0)
            b_(4,i) = r*(4.0d0-4.0d0*r)-4.0d0*r*s
            b_(5,i) = 4.0d0*r*s
            b_(6,i) = s*(-4.0d0*s-4.0d0*r+4.0d0)
         endif
      end do
      end subroutine


      subroutine ensBASIS_TRIA_LAGRBUBBLE2(trp_,trb_,n_,
     &     p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,*)
      double precision p_(poff_,*)
      double precision r,s
      character trp_
      character trb_
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T')then
            b_(i,1) = s*(2.0d0*s+4.0d0*r-3.0d0)+
     &           r*(2.0d0*r-3.0d0)+1.0d0
            b_(i,2) = r*(2.0d0*r-1.0d0)
            b_(i,3) = s*(2.0d0*s-1.0d0)
            b_(i,4) = r*(4.0d0-4.0d0*r)-4.0d0*r*s
            b_(i,5) = 4.0d0*r*s
            b_(i,6) = s*(-4.0d0*s-4.0d0*r+4.0d0)
            b_(i,7) = s*(r*(2.7d1-2.7d1*r)-2.7d1*r*s)
         else if (trb_.eq.'N')then
            b_(1,i) = s*(2.0d0*s+4.0d0*r-3.0d0)+
     &           r*(2.0d0*r-3.0d0)+1.0d0
            b_(2,i) = r*(2.0d0*r-1.0d0)
            b_(3,i) = s*(2.0d0*s-1.0d0)
            b_(4,i) = r*(4.0d0-4.0d0*r)-4.0d0*r*s
            b_(5,i) = 4.0d0*r*s
            b_(6,i) = s*(-4.0d0*s-4.0d0*r+4.0d0)
            b_(7,i) = s*(r*(2.7d1-2.7d1*r)-2.7d1*r*s)
         endif
      end do
      end subroutine








      subroutine ensBASIS_TRIA_LAGR1(trp_,trb_,n_,p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,*)
      double precision p_(poff_,*)
      double precision r,s
      character trp_
      character trb_
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T')then
         b_(i,1) = 1.d0-(r+s)
         b_(i,2) = r
         b_(i,3) = s
         else if (trb_.eq.'N')then
         b_(1,i) = 1.d0-(r+s)
         b_(2,i) = r
         b_(3,i) = s
         endif
      end do
      end subroutine


      subroutine ensBASIS_TRIA_LAGRBUBBLE1(trp_,trb_,n_,
     &     p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,*)
      double precision p_(poff_,*)
      double precision r,s
      character trp_
      character trb_
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T')then
         b_(i,1) = 1.d0-(r+s)
         b_(i,2) = r
         b_(i,3) = s
         b_(i,4) = s*(r*(2.7d1-2.7d1*r)-2.7d1*r*s)         
         else if (trb_.eq.'N')then
         b_(1,i) = 1.d0-(r+s)
         b_(2,i) = r
         b_(3,i) = s
         b_(4,i) = s*(r*(2.7d1-2.7d1*r)-2.7d1*r*s)         
         endif
      end do
      end subroutine


      subroutine ensBASIS_TRIA_LAGR0(trp_,trb_,n_,p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,*)
      double precision p_(poff_,*)
      double precision r,s
      character trp_
      character trb_
      integer i
      do i=1,n_
         if (trb_.eq.'T')then
         b_(i,1) = 1.d0
         else if (trb_.eq.'N')then
         b_(1,i) = 1.d0
         endif
      end do
      end subroutine



      subroutine ensBASIS_TRIA_LAGRBUBBLE0(trp_,trb_,n_,p_,
     &     poff_,b_,boff_)
      implicit none
      integer 		n_,poff_,boff_
      double precision 	b_(boff_,*)
      double precision 	p_(poff_,*)
      double precision 	r,s
      character trp_
      character trb_
      integer i
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T')then
            b_(i,1) = 1.d0
            b_(i,2) = s*(r*(2.7d1-2.7d1*r)-2.7d1*r*s)         
         else if (trb_.eq.'N')then
            b_(1,i) = 1.d0
            b_(2,i) = s*(r*(2.7d1-2.7d1*r)-2.7d1*r*s)         
         endif
      end do
      end subroutine



