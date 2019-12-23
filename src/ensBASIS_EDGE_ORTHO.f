      subroutine ensbasis_edge_ortho6(trp_,trc_,n_,p_,poff_,c_,boff_) 
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
            c_(i,1) = 7.071067811865475d-1
            c_(i,2) = 1.224744871391589d0*r
            c_(i,3) = 3.535533905932738d-1*(6.708203932499369d0*r**2
     &           -2.23606797749979d0)
            c_(i,4) = 3.535533905932738d-1*r*(1.322875655532295d1*r**2
     &           -7.937253933193772d0)
            c_(i,5) = 8.838834764831845d-2*(1.0d0*r**2*(1.05d2*r**2
     &           -9.0d1)+9.0d0)
            c_(i,6) = 8.838834764831845d-2*r*(1.0d0*r**2*(
     &           2.089473617923902d2*r**2-2.32163735324878d2)
     &           +4.9749371855331d1)
            c_(i,7) = 4.419417382415922d-2*(1.0d0*r**2*(1.0d0*r**2*(
     &           8.328823446321815d2*r**2-1.135748651771157d3)
     &           +3.785828839237189d2)-1.802775637731995d1)
         else if (trc_.eq.'N') then  
            c_(1,i) = 7.071067811865475d-1
            c_(2,i) = 1.224744871391589d0*r
            c_(3,i) = 3.535533905932738d-1*(6.708203932499369d0*r**2
     &           -2.23606797749979d0)
            c_(4,i) = 3.535533905932738d-1*r*(1.322875655532295d1*r**2
     &           -7.937253933193772d0)
            c_(5,i) = 8.838834764831845d-2*(1.0d0*r**2*(1.05d2*r**2
     &           -9.0d1)+9.0d0)
            c_(6,i) = 8.838834764831845d-2*r*(1.0d0*r**2*(
     &           2.089473617923902d2*r**2-2.32163735324878d2)
     &           +4.9749371855331d1)
            c_(7,i) = 4.419417382415922d-2*(1.0d0*r**2*(1.0d0*r**2*(
     &           8.328823446321815d2*r**2-1.135748651771157d3)
     &           +3.785828839237189d2)-1.802775637731995d1)
         endif  
      end do 
      end subroutine 
      subroutine ensbasis_edge_drortho6(trp_,trc_,n_,p_,poff_,c_,boff_) 
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
         c_(i,1) = 0.0d0
         c_(i,2) = 1.224744871391589d0
         c_(i,3) = 4.743416490252569d0*r
         c_(i,4) = 3.535533905932738d-1*(3.968626966596886d1*r**2
     &        -7.937253933193772d0)
         c_(i,5) = 8.838834764831845d-2*r*(4.2d2*r**2-1.8d2)
         c_(i,6) = 8.838834764831845d-2*(1.0d0*r**2*(
     &        1.044736808961951d3*r**2
     &        -6.96491205974634d2)+4.9749371855331d1)
         c_(i,7) = 4.419417382415922d-2*r*(1.0d0*r**2*(
     &        4.997294067793089d3*r**2-4.542994607084627d3)
     &        +7.571657678474377d2)
      else if (trc_.eq.'N') then  
         c_(1,i) = 0.0d0
         c_(2,i) = 1.224744871391589d0
         c_(3,i) = 4.743416490252569d0*r
         c_(4,i) = 3.535533905932738d-1*(3.968626966596886d1*r**2
     &        -7.937253933193772d0)
         c_(5,i) = 8.838834764831845d-2*r*(4.2d2*r**2-1.8d2)
         c_(6,i) = 8.838834764831845d-2*(1.0d0*r**2*(
     &        1.044736808961951d3*r**2-6.96491205974634d2)
     &        +4.9749371855331d1)
         c_(7,i) = 4.419417382415922d-2*r*(1.0d0*r**2*(
     &        4.997294067793089d3*r**2-4.542994607084627d3)
     &        +7.571657678474377d2)
      endif  
      end do 
      end subroutine 
      
      
      subroutine ensbasis_edge_ortho5(trp_,trc_,n_,p_,poff_,c_,boff_) 
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
         c_(i,1) = 7.071067811865475d-1
         c_(i,2) = 1.224744871391589d0*r
         c_(i,3) = 3.535533905932738d-1*(6.708203932499369d0*r**2
     &        -2.23606797749979d0)
         c_(i,4) = 3.535533905932738d-1*r*(1.322875655532295d1*r**2
     &        -7.937253933193772d0)
         c_(i,5) = 8.838834764831845d-2*(1.0d0*r**2*(1.05d2*r**
     &        2-9.0d1)+9.0d0)
         c_(i,6) = 8.838834764831845d-2*r*(1.0d0*r**2*(
     &        2.089473617923902d2*r**2-2.32163735324878d2)
     &        +4.9749371855331d1)
      else if (trc_.eq.'N') then  
         c_(1,i) = 7.071067811865475d-1
         c_(2,i) = 1.224744871391589d0*r
         c_(3,i) = 3.535533905932738d-1*(6.708203932499369d0*r**2
     &        -2.23606797749979d0)
         c_(4,i) = 3.535533905932738d-1*r*(1.322875655532295d1*r**2
     &        -7.937253933193772d0)
         c_(5,i) = 8.838834764831845d-2*(1.0d0*r**2*(1.05d2*r**
     &        2-9.0d1)+9.0d0)
         c_(6,i) = 8.838834764831845d-2*r*(1.0d0*r**2*(
     &        2.089473617923902d2*r**2-2.32163735324878d2)
     &        +4.9749371855331d1)
      endif  
      end do 
      end subroutine 
      
      
      subroutine ensbasis_edge_drortho5(trp_,trc_,n_,p_,poff_,c_,boff_) 
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
         c_(i,1) = 0.0d0
         c_(i,2) = 1.224744871391589d0
         c_(i,3) = 4.743416490252569d0*r
         c_(i,4) = 3.535533905932738d-1*(3.968626966596886d1*r**2
     &        -7.937253933193772d0)
         c_(i,5) = 8.838834764831845d-2*r*(4.2d2*r**2-1.8d2)
         c_(i,6) = 8.838834764831845d-2*(1.0d0*r**2*(
     &        1.044736808961951d3*r**2-6.96491205974634d2)
     &        +4.9749371855331d1)
      else if (trc_.eq.'N') then  
         c_(1,i) = 0.0d0
         c_(2,i) = 1.224744871391589d0
         c_(3,i) = 4.743416490252569d0*r
         c_(4,i) = 3.535533905932738d-1*(3.968626966596886d1*r**2
     &        -7.937253933193772d0)
         c_(5,i) = 8.838834764831845d-2*r*(4.2d2*r**2-1.8d2)
         c_(6,i) = 8.838834764831845d-2*(1.0d0*r**2*(
     &        1.044736808961951d3*r**2-6.96491205974634d2)
     &        +4.9749371855331d1)
      endif  
      end do 
      end subroutine 
      
      
      subroutine ensbasis_edge_ortho4(trp_,trc_,n_,p_,poff_,c_,boff_) 
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
         c_(i,1) = 7.071067811865475d-1
         c_(i,2) = 1.224744871391589d0*r
         c_(i,3) = 3.535533905932738d-1*(6.708203932499369d0*r**2
     &        -2.23606797749979d0)
         c_(i,4) = 3.535533905932738d-1*r*(1.322875655532295d1*r**2
     &        -7.937253933193772d0)
         c_(i,5) = 8.838834764831845d-2*(1.0d0*r**2*(1.05d2*r**2
     &        -9.0d1)+9.0d0)
      else if (trc_.eq.'N') then  
         c_(1,i) = 7.071067811865475d-1
         c_(2,i) = 1.224744871391589d0*r
         c_(3,i) = 3.535533905932738d-1*(6.708203932499369d0*r**2
     &        -2.23606797749979d0)
         c_(4,i) = 3.535533905932738d-1*r*(1.322875655532295d1*r**2
     &        -7.937253933193772d0)
         c_(5,i) = 8.838834764831845d-2*(1.0d0*r**2*(1.05d2*r**2
     &        -9.0d1)+9.0d0)
      endif  
      end do 
      end subroutine 
      subroutine ensbasis_edge_drortho4(trp_,trc_,n_,p_,poff_,c_,boff_) 
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
         c_(i,1) = 0.0d0
         c_(i,2) = 1.224744871391589d0
         c_(i,3) = 4.743416490252569d0*r
         c_(i,4) = 3.535533905932738d-1*(3.968626966596886d1*r**2
     &        -7.937253933193772d0)
         c_(i,5) = 8.838834764831845d-2*r*(4.2d2*r**2-1.8d2)
      else if (trc_.eq.'N') then  
         c_(1,i) = 0.0d0
         c_(2,i) = 1.224744871391589d0
         c_(3,i) = 4.743416490252569d0*r
         c_(4,i) = 3.535533905932738d-1*(3.968626966596886d1*r**2
     &        -7.937253933193772d0)
         c_(5,i) = 8.838834764831845d-2*r*(4.2d2*r**2-1.8d2)
      endif  
      end do 
      end subroutine 
      
      
      subroutine ensbasis_edge_ortho3(trp_,trc_,n_,p_,poff_,c_,boff_) 
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
         c_(i,1) = 7.071067811865475d-1
         c_(i,2) = 1.224744871391589d0*r
         c_(i,3) = 3.535533905932738d-1*(6.708203932499369d0*r**2
     &        -2.23606797749979d0)
         c_(i,4) = 3.535533905932738d-1*r*(1.322875655532295d1*r**2
     &        -7.937253933193772d0)
      else if (trc_.eq.'N') then  
         c_(1,i) = 7.071067811865475d-1
         c_(2,i) = 1.224744871391589d0*r
         c_(3,i) = 3.535533905932738d-1*(6.708203932499369d0*r**2
     &        -2.23606797749979d0)
         c_(4,i) = 3.535533905932738d-1*r*(1.322875655532295d1*r**2
     &        -7.937253933193772d0)
      endif  
      end do 
      end subroutine 
      subroutine ensbasis_edge_drortho3(trp_,trc_,n_,p_,poff_,c_,boff_) 
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
         c_(i,1) = 0.0d0
         c_(i,2) = 1.224744871391589d0
         c_(i,3) = 4.743416490252569d0*r
         c_(i,4) = 3.535533905932738d-1*(3.968626966596886d1*r**2
     &        -7.937253933193772d0)
      else if (trc_.eq.'N') then  
         c_(1,i) = 0.0d0
         c_(2,i) = 1.224744871391589d0
         c_(3,i) = 4.743416490252569d0*r
         c_(4,i) = 3.535533905932738d-1*(3.968626966596886d1*r**2
     &        -7.937253933193772d0)
      endif  
      end do 
      end subroutine 
      
      
      
      subroutine ensbasis_edge_ortho2(trp_,trc_,n_,p_,poff_,c_,boff_) 
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
         c_(i,1) = 7.071067811865475d-1
         c_(i,2) = 1.224744871391589d0*r
         c_(i,3) = 3.535533905932738d-1*(6.708203932499369d0*r**2
     &        -2.23606797749979d0)
      else if (trc_.eq.'N') then  
         c_(1,i) = 7.071067811865475d-1
         c_(2,i) = 1.224744871391589d0*r
         c_(3,i) = 3.535533905932738d-1*(6.708203932499369d0*r**2
     &        -2.23606797749979d0)
      endif  
      end do 
      end subroutine 
      subroutine ensbasis_edge_drortho2(trp_,trc_,n_,p_,poff_,c_,boff_) 
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
         c_(i,1) = 0.0d0
         c_(i,2) = 1.224744871391589d0
         c_(i,3) = 4.743416490252569d0*r
      else if (trc_.eq.'N') then  
         c_(1,i) = 0.0d0
         c_(2,i) = 1.224744871391589d0
         c_(3,i) = 4.743416490252569d0*r
      endif  
      end do 
      end subroutine 

      subroutine ensbasis_edge_ortho1(trp_,trc_,n_,p_,poff_,c_,boff_) 
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
         c_(i,1) = 7.071067811865475d-1
         c_(i,2) = 1.224744871391589d0*r
      else if (trc_.eq.'N') then  
         c_(1,i) = 7.071067811865475d-1
         c_(2,i) = 1.224744871391589d0*r
      endif  
      end do 
      end subroutine 
      
      subroutine ensbasis_edge_drortho1(trp_,trc_,n_,p_,poff_,c_,boff_) 
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
         c_(i,1) = 0.0d0
         c_(i,2) = 1.224744871391589d0
      else if (trc_.eq.'N') then  
         c_(1,i) = 0.0d0
         c_(2,i) = 1.224744871391589d0
      endif  
      end do 
      end subroutine 

      subroutine ensbasis_edge_ortho0 (trp_,trc_,n_,p_,poff_,c_,boff_) 
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
         c_(i,1) = 7.071067811865475d-1
      else if (trc_.eq.'N') then  
         c_(1,i) = 7.071067811865475d-1
      endif  
      end do 
      end subroutine 

      subroutine ensbasis_edge_drortho0(trp_,trc_,n_,p_,poff_,c_,boff_) 
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
         c_(i,1) = 0.0d0
      else if (trc_.eq.'N') then  
         c_(1,i) = 0.0d0
      endif  
      end do 
      end subroutine 
