      
      subroutine ensBASIS_TRIA_DRORTHO0(trp_,trb_,n_,p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,n_)
      double precision p_(poff_,n_)
      double precision r,s
      integer i
      character trp_,trb_
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T') then 
            b_(i,1) = 0.d0
         else if (trb_.eq.'N') then 
            b_(1,i) = 0.d0
         endif
      end do
      end subroutine       

      subroutine ensBASIS_TRIA_DRORTHO1(trp_,trb_,n_,p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,n_)
      double precision p_(poff_,n_)
      double precision r,s
      integer i
      character trp_,trb_
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T') then 
         b_(i,1) = 0.0d0                                                         
         b_(i,2) = 6.0d0                                                         
         b_(i,3) = 3.464101615137755d0                                           
         else if (trb_.eq.'N') then 
         b_(1,i) = 0.0d0                                                         
         b_(2,i) = 6.0d0                                                         
         b_(3,i) = 3.464101615137755d0                                           
         endif
      end do
      end subroutine       
      subroutine ensBASIS_TRIA_DRORTHO2(trp_,trb_,n_,p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,n_)
      double precision p_(poff_,n_)
      double precision r,s
      integer i
      character trp_,trb_
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T') then 

       b_(i,1) = 0.0d0                                                         
       b_(i,2) = 6.0d0                                                         
       b_(i,3) = 3.464101615137755d0                                           
       b_(i,4) = 4.898979485566356d1*r-1.959591794226543d1                     
       b_(i,5) = 4.242640687119285d1*s+4.242640687119285d1*r
     1    -2.545584412271571d1                                                             
       b_(i,6) = 3.286335345030997d1*s+1.095445115010332d1*r
     1    -1.095445115010332d1                                                             
         else if (trb_.eq.'N') then 


       b_(1,i) = 0.0d0                                                         
       b_(2,i) = 6.0d0                                                         
       b_(3,i) = 3.464101615137755d0                                           
       b_(4,i) = 4.898979485566356d1*r-1.959591794226543d1                     
       b_(5,i) = 4.242640687119285d1*s+4.242640687119285d1*r
     1    -2.545584412271571d1                                                             
       b_(6,i) = 3.286335345030997d1*s+1.095445115010332d1*r
     1    -1.095445115010332d1                                                             


         endif

      end do
      end subroutine       
      subroutine ensBASIS_TRIA_DRORTHO3(trp_,trb_,n_,p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,n_)
      double precision p_(poff_,n_)


      double precision r,s
      integer i
      character trp_,trb_
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T') then 

         b_(i,1) = 0.0d0                                                         
         b_(i,2) = 6.0d0                                                         
         b_(i,3) = 3.464101615137755d0                                           
         b_(i,4) = 4.898979485566356d1*r-1.959591794226543d1                     
         b_(i,5) = 4.242640687119285d1*s+4.242640687119285d1*r
     1        -2.545584412271571d1                                                             
         b_(i,6) = 3.286335345030997d1*s+1.095445115010332d1*r
     1        -1.095445115010332d1                                                             
         b_(i,7) = r*(2.9698484809835d2*r-2.545584412271571d2)        
     1        +4.242640687119285d1                                                       
         b_(i,8) = 3.333333333333333d-1*((1.234542830362722d3*r        
     1        -3.527265229607776d2)*s+r*(9.259071227720413d2*r        
     2        -9.699979381421385d2)+1.910601999370879d2)                                        
         b_(i,9) = s*(2.656313234541439d2*s+5.312626469082877d2*r        
     1        -3.035786553761644d2)+r*(1.328156617270719d2*r        
     2        -1.897366596101028d2)+5.692099788303083d1                                         
         b_(i,10) = s*(2.244994432064365d2*s+1.795995545651492d2*        
     1        r-1.795995545651492d2)+r*(2.244994432064365d1*r        
     2        -4.48998886412873d1)+2.244994432064365d1                                         
         




      else         if (trb_.eq.'N') then 

         b_(1,i) = 0.0d0                                                         
         b_(2,i) = 6.0d0                                                         
         b_(3,i) = 3.464101615137755d0                                           
         b_(4,i) = 4.898979485566356d1*r-1.959591794226543d1                     
         b_(5,i) = 4.242640687119285d1*s+4.242640687119285d1*r
     1        -2.545584412271571d1                                                             
         b_(6,i) = 3.286335345030997d1*s+1.095445115010332d1*r
     1        -1.095445115010332d1                                                             
         b_(7,i) = r*(2.9698484809835d2*r-2.545584412271571d2)        
     1        +4.242640687119285d1                                                       
         b_(8,i) = 3.333333333333333d-1*((1.234542830362722d3*r        
     1        -3.527265229607776d2)*s+r*(9.259071227720413d2*r        
     2        -9.699979381421385d2)+1.910601999370879d2)                                        
         b_(9,i) = s*(2.656313234541439d2*s+5.312626469082877d2*r        
     1        -3.035786553761644d2)+r*(1.328156617270719d2*r        
     2        -1.897366596101028d2)+5.692099788303083d1                                         
         b_(10,i) = s*(2.244994432064365d2*s+1.795995545651492d2*        
     1        r-1.795995545651492d2)+r*(2.244994432064365d1*r        
     2        -4.48998886412873d1)+2.244994432064365d1                                         
         




      endif




      end do
      end subroutine       
      subroutine ensBASIS_TRIA_DRORTHO4(trp_,trb_,n_,p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,n_)
      double precision p_(poff_,n_)


      double precision r,s
      integer i
      character trp_,trb_
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T') then 


         b_(i,1) = 0.0d0                                                         
         b_(i,2) = 6.0d0                                                         
         b_(i,3) = 3.464101615137755d0                                           
         b_(i,4) = 4.898979485566356d1*r-1.959591794226543d1                     
         b_(i,5) = 4.242640687119285d1*s+4.242640687119285d1*r
     1        -2.545584412271571d1                                                             
         b_(i,6) = 3.286335345030997d1*s+1.095445115010332d1*r
     1        -1.095445115010332d1                                                             
         b_(i,7) = r*(2.9698484809835d2*r-2.545584412271571d2)        
     1        +4.242640687119285d1                                                       
         b_(i,8) = 3.333333333333333d-1*((1.234542830362722d3*r        
     1        -3.527265229607776d2)*s+r*(9.259071227720413d2*r        
     2        -9.699979381421385d2)+1.910601999370879d2)                                        
         b_(i,9) = s*(2.656313234541439d2*s+5.312626469082877d2*r        
     1        -3.035786553761644d2)+r*(1.328156617270719d2*r        
     2        -1.897366596101028d2)+5.692099788303083d1                                         
         b_(i,10) = s*(2.244994432064365d2*s+1.795995545651492d2*        
     1        r-1.795995545651492d2)+r*(2.244994432064365d1*r        
     2        -4.48998886412873d1)+2.244994432064365d1                                         
         b_(i,11) = r*(r*(1.593787940724863d3*r        
     1        -2.125050587633151d3)+7.968939703624316d2)
     &        -7.58946638440411d1                       
         b_(i,12) = (r*(2.760521689826037d3*r        
     1        -1.840347793217358d3)+2.300434741521698d2)*s+r*(r        
     2        *(1.840347793217358d3*r-2.760521689826037d3)
     3        +1.150217370760849d3)-1.204989626511366d2                                                                
         b_(i,13) = 2.5d-1*(s*((1.221880517890354d4*r
     &        -2.71529        
     1        0039756342d3)*s+r*(1.832820776835531d4*r-1.7649385258416        
     2        22d4)+2.884995667241114d3)+r*(r*(4.072935059634513        
     3        d3*r-7.467047609329941d3)+3.903229432149742d3)
     &        -5.091168824543142d2)                                                                 
         b_(i,14) = s*(s*(1.505988047761336d3*s+4.51796414328        
     1        4008d3*r-2.509980079602227d3)+r*(2.710778485970405d3*r-3        
     2        .815169720995385d3)+1.10439123502498d3)+r*(r*(        
     3        3.011976095522672d2*r-7.027944222886235d2)
     &        +5.019960159204453d2)        
     4        -1.003992031840891d2                                                     
         b_(i,15) = 2.0d-1*(s*(s*(6.640783086353596d3*s+8.538        
     1        149682454624d3*r-8.538149682454624d3)+r*(2.8460498941515        
     2        41d3*r-5.692099788303083d3)+2.846049894151541d3)+r*(1.0d        
     3        0*r*(1.897366596101028d2*r-5.692099788303083d2)        
     4        +5.692099788303083d2)-1.897366596101028d2)                                           
         




      else          if (trb_.eq.'N') then 

         b_(1,i) = 0.0d0                                                         
         b_(2,i) = 6.0d0                                                         
         b_(3,i) = 3.464101615137755d0                                           
         b_(4,i) = 4.898979485566356d1*r-1.959591794226543d1                     
         b_(5,i) = 4.242640687119285d1*s+4.242640687119285d1*r
     1        -2.545584412271571d1                                                             
         b_(6,i) = 3.286335345030997d1*s+1.095445115010332d1*r
     1        -1.095445115010332d1                                                             
         b_(7,i) = r*(2.9698484809835d2*r-2.545584412271571d2)        
     1        +4.242640687119285d1                                                       
         b_(8,i) = 3.333333333333333d-1*((1.234542830362722d3*r        
     1        -3.527265229607776d2)*s+r*(9.259071227720413d2*r        
     2        -9.699979381421385d2)+1.910601999370879d2)                                        
         b_(9,i) = s*(2.656313234541439d2*s+5.312626469082877d2*r        
     1        -3.035786553761644d2)+r*(1.328156617270719d2*r        
     2        -1.897366596101028d2)+5.692099788303083d1                                         
         b_(10,i) = s*(2.244994432064365d2*s+1.795995545651492d2*        
     1        r-1.795995545651492d2)+r*(2.244994432064365d1*r        
     2        -4.48998886412873d1)+2.244994432064365d1                                         
         b_(11,i) = r*(r*(1.593787940724863d3*r        
     1        -2.125050587633151d3)+7.968939703624316d2)
     &        -7.58946638440411d1                       
         b_(12,i) = (r*(2.760521689826037d3*r        
     1        -1.840347793217358d3)+2.300434741521698d2)*s+r*(r        
     2        *(1.840347793217358d3*r-2.760521689826037d3)
     3        +1.150217370760849d3)-1.204989626511366d2                                                                
         b_(13,i) = 2.5d-1*(s*((1.221880517890354d4*r
     &        -2.71529        
     1        0039756342d3)*s+r*(1.832820776835531d4*r-1.7649385258416        
     2        22d4)+2.884995667241114d3)+r*(r*(4.072935059634513        
     3        d3*r-7.467047609329941d3)+3.903229432149742d3)
     &        -5.091168824543142d2)                                                                 
         b_(14,i) = s*(s*(1.505988047761336d3*s+4.51796414328        
     1        4008d3*r-2.509980079602227d3)+r*(2.710778485970405d3*r-3        
     2        .815169720995385d3)+1.10439123502498d3)+r*(r*(        
     3        3.011976095522672d2*r-7.027944222886235d2)
     &        +5.019960159204453d2)        
     4        -1.003992031840891d2                                                     
         b_(15,i) = 2.0d-1*(s*(s*(6.640783086353596d3*s+8.538        
     1        149682454624d3*r-8.538149682454624d3)+r*(2.8460498941515        
     2        41d3*r-5.692099788303083d3)+2.846049894151541d3)+r*(1.0d        
     3        0*r*(1.897366596101028d2*r-5.692099788303083d2)        
     4        +5.692099788303083d2)-1.897366596101028d2)                                           
         





      endif




      end do
      end subroutine       
      subroutine ensBASIS_TRIA_DRORTHO5(trp_,trb_,n_,p_,poff_,b_,boff_)
      implicit none
      integer n_,poff_,boff_
      double precision b_(boff_,n_)
      double precision p_(poff_,n_)

      double precision r,s
      integer i
      character trp_,trb_
      do i=1,n_
         if (trp_.eq.'N') then 
            r = p_(1,i)
            s = p_(2,i)
         else if (trp_.eq.'T') then 
            r = p_(i,1)
            s = p_(i,2)
         endif
         if (trb_.eq.'T') then 

       b_(i,1) = 0.0d0                                                         
       b_(i,2) = 6.0d0                                                         
       b_(i,3) = 3.464101615137755d0                                           
       b_(i,4) = 4.898979485566356d1*r-1.959591794226543d1                     
       b_(i,5) = 4.242640687119285d1*s+4.242640687119285d1*r
     1    -2.545584412271571d1                                                             
       b_(i,6) = 3.286335345030997d1*s+1.095445115010332d1*r
     1    -1.095445115010332d1                                                             
       b_(i,7) = r*(2.9698484809835d2*r-2.545584412271571d2)        
     1    +4.242640687119285d1                                                       
       b_(i,8) = 3.333333333333333d-1*((1.234542830362722d3*r        
     1    -3.527265229607776d2)*s+r*(9.259071227720413d2*r        
     2    -9.699979381421385d2)+1.910601999370879d2)                                        
       b_(i,9) = s*(2.656313234541439d2*s+5.312626469082877d2*r        
     1    -3.035786553761644d2)+r*(1.328156617270719d2*r        
     2    -1.897366596101028d2)+5.692099788303083d1                                         
       b_(i,10) = s*(2.244994432064365d2*s+1.795995545651492d2*        
     1    r-1.795995545651492d2)+r*(2.244994432064365d1*r        
     2    -4.48998886412873d1)+2.244994432064365d1                                         
       b_(i,11) = r*(r*(1.593787940724863d3*r        
     1    -2.125050587633151d3)+7.968939703624316d2)-7.58946638440411d1                       
       b_(i,12) = (r*(2.760521689826037d3*r        
     1    -1.840347793217358d3)+2.300434741521698d2)*s+r*(r        
     2    *(1.840347793217358d3*r-2.760521689826037d3)
     3    +1.150217370760849d3)-1.204989626511366d2                                                                
       b_(i,13) = 2.5d-1*(s*((1.221880517890354d4*r
     &      -2.71529        
     1    0039756342d3)*s+r*(1.832820776835531d4*r-1.7649385258416        
     2    22d4)+2.884995667241114d3)+r*(r*(4.072935059634513        
     3    d3*r-7.467047609329941d3)+3.903229432149742d3)-5.0911688245431        
     4    42d2)                                                                 
       b_(i,14) = s*(s*(1.505988047761336d3*s+4.51796414328        
     1    4008d3*r-2.509980079602227d3)+r*(2.710778485970405d3*r-3        
     2    .815169720995385d3)+1.10439123502498d3)+r*(r*(3.01        
     3    1976095522672d2*r-7.027944222886235d2)+5.019960159204453d2)-1.        
     4    003992031840891d2                                                     
       b_(i,15) = 2.0d-1*(s*(s*(6.640783086353596d3*s+8.538        
     1    149682454624d3*r-8.538149682454624d3)+r*(2.8460498941515        
     2    41d3*r-5.692099788303083d3)+2.846049894151541d3)+r*(1.0d        
     3    0*r*(1.897366596101028d2*r-5.692099788303083d2)+5.692099788303        
     4    083d2)-1.897366596101028d2)                                           
       b_(i,16) = 3.333333333333333d-1*(r*(r*(r*(2.40        
     1    0622419290464d4*r-4.364768035073571d4)+2.618860821044143d4)-5.        
     2    819690713431428d3)+3.637306695894642d2)                               
       b_(i,17) = (r*(r*(1.584d4*r-1.728d4)+5.184d3)-        
     1    3.84d2)*s+r*(r*(r*(9.9d3*r-1.944d4)+1.2528d4        
     2    )-2.976d3)+1.98d2                                                     
       b_(i,18) = 6.666666666666667d-2*(s*((r*(3.4508        
     1    28161470809d5*r-1.882269906256805d5)+1.882269906256805d4)*s+1.        
     2    0d0*r*(r*(4.601104215294411d5*r-6.274233020856016d5)+2.2        
     3    58723887508166d5)-1.951983606488538d4)+r*(r*(        
     4    r*(9.585633781863357d4*r-2.161124707183739d5)+1.6103864753530        
     5    44d5)-4.415201014676455d4)+3.369495511200453d3)                       
       b_(i,19) = s*(s*((2.01633330578057d4*r-3.66606        
     1    0555964672d3)*s+r*(4.536749938006282d4*r-4.1243181254602        
     2    56d4)+5.774045375644359d3)+r*(r*(2.419599966936684        
     3    d4*r-4.289290850478666d4)+2.111650880235651d4)-2.4195999669366        
     4    84d3)+r*(r*(r*(2.520416632225712d3*r-6.78221        
     5    2028534643d3)+6.213972642360119d3)-2.162975728019157d3)+2.1079        
     6    84819679686d2                                                         
       b_(i,20) = s*(s*(s*(8.002074730968213d3*s+3.20        
     1    0829892387285d4*r-1.745907214029428d4)+r*(3.086514539087        
     2    74d4*r-4.302414206001092d4)+1.215899666913352d4)+r*(
     3    r*(9.145228263963672d3*r-2.120030188464306d4)+1.496491897739        
     4    51d4)-2.909845356715714d3)+r*(r*(r*(5.715767        
     5    664977295d2*r-1.870614872174388d3)+2.182384017536786d3)-1.0392        
     6    30484541326d3)+1.55884572681199d2                                     
       b_(i,21) = s*(s*(s*(7.238148934637916d3*s+1.28        
     1    6782032824518d4*r-1.286782032824518d4)+r*(7.238148934637        
     2    916d3*r-1.447629786927583d4)+7.238148934637916d3)+r*(1.0        
     3    d0*r*(1.378695035169127d3*r-4.136085105507381d3)+4.13608510550        
     4    7381d3)-1.378695035169127d3)+r*(r*(r*(5.7445        
     5    62646538029d1*r-2.297825058615211d2)+3.446737587922817d2)-2.29        
     6    7825058615211d2)+5.744562646538029d1                                  




         else if (trb_.eq.'N') then 


       b_(1,i) = 0.0d0                                                         
       b_(2,i) = 6.0d0                                                         
       b_(3,i) = 3.464101615137755d0                                           
       b_(4,i) = 4.898979485566356d1*r-1.959591794226543d1                     
       b_(5,i) = 4.242640687119285d1*s+4.242640687119285d1*r
     1    -2.545584412271571d1                                                             
       b_(6,i) = 3.286335345030997d1*s+1.095445115010332d1*r
     1    -1.095445115010332d1                                                             
       b_(7,i) = r*(2.9698484809835d2*r-2.545584412271571d2)        
     1    +4.242640687119285d1                                                       
       b_(8,i) = 3.333333333333333d-1*((1.234542830362722d3*r        
     1    -3.527265229607776d2)*s+r*(9.259071227720413d2*r        
     2    -9.699979381421385d2)+1.910601999370879d2)                                        
       b_(9,i) = s*(2.656313234541439d2*s+5.312626469082877d2*r        
     1    -3.035786553761644d2)+r*(1.328156617270719d2*r        
     2    -1.897366596101028d2)+5.692099788303083d1                                         
       b_(10,i) = s*(2.244994432064365d2*s+1.795995545651492d2*        
     1    r-1.795995545651492d2)+r*(2.244994432064365d1*r        
     2    -4.48998886412873d1)+2.244994432064365d1                                         
       b_(11,i) = r*(r*(1.593787940724863d3*r        
     1    -2.125050587633151d3)+7.968939703624316d2)-7.58946638440411d1                       
       b_(12,i) = (r*(2.760521689826037d3*r        
     1    -1.840347793217358d3)+2.300434741521698d2)*s+r*(r        
     2    *(1.840347793217358d3*r-2.760521689826037d3)
     3    +1.150217370760849d3)-1.204989626511366d2                                                                
       b_(13,i) = 2.5d-1*(s*((1.221880517890354d4*r
     &      -2.71529        
     1    0039756342d3)*s+r*(1.832820776835531d4*r-1.7649385258416        
     2    22d4)+2.884995667241114d3)+r*(r*(4.072935059634513        
     3    d3*r-7.467047609329941d3)+3.903229432149742d3)-5.0911688245431        
     4    42d2)                                                                 
       b_(14,i) = s*(s*(1.505988047761336d3*s+4.51796414328        
     1    4008d3*r-2.509980079602227d3)+r*(2.710778485970405d3*r-3        
     2    .815169720995385d3)+1.10439123502498d3)+r*(r*(3.01        
     3    1976095522672d2*r-7.027944222886235d2)+5.019960159204453d2)-1.        
     4    003992031840891d2                                                     
       b_(15,i) = 2.0d-1*(s*(s*(6.640783086353596d3*s+8.538        
     1    149682454624d3*r-8.538149682454624d3)+r*(2.8460498941515        
     2    41d3*r-5.692099788303083d3)+2.846049894151541d3)+r*(1.0d        
     3    0*r*(1.897366596101028d2*r-5.692099788303083d2)+5.692099788303        
     4    083d2)-1.897366596101028d2)                                           
       b_(16,i) = 3.333333333333333d-1*(r*(r*(r*(2.40        
     1    0622419290464d4*r-4.364768035073571d4)+2.618860821044143d4)-5.        
     2    819690713431428d3)+3.637306695894642d2)                               
       b_(17,i) = (r*(r*(1.584d4*r-1.728d4)+5.184d3)-        
     1    3.84d2)*s+r*(r*(r*(9.9d3*r-1.944d4)+1.2528d4        
     2    )-2.976d3)+1.98d2                                                     
       b_(18,i) = 6.666666666666667d-2*(s*((r*(3.4508        
     1    28161470809d5*r-1.882269906256805d5)+1.882269906256805d4)*s+1.        
     2    0d0*r*(r*(4.601104215294411d5*r-6.274233020856016d5)+2.2        
     3    58723887508166d5)-1.951983606488538d4)+r*(r*(        
     4    r*(9.585633781863357d4*r-2.161124707183739d5)+1.6103864753530        
     5    44d5)-4.415201014676455d4)+3.369495511200453d3)                       
       b_(19,i) = s*(s*((2.01633330578057d4*r-3.66606        
     1    0555964672d3)*s+r*(4.536749938006282d4*r-4.1243181254602        
     2    56d4)+5.774045375644359d3)+r*(r*(2.419599966936684        
     3    d4*r-4.289290850478666d4)+2.111650880235651d4)-2.4195999669366        
     4    84d3)+r*(r*(r*(2.520416632225712d3*r-6.78221        
     5    2028534643d3)+6.213972642360119d3)-2.162975728019157d3)+2.1079        
     6    84819679686d2                                                         
       b_(20,i) = s*(s*(s*(8.002074730968213d3*s+3.20        
     1    0829892387285d4*r-1.745907214029428d4)+r*(3.086514539087        
     2    74d4*r-4.302414206001092d4)+1.215899666913352d4)+r*(
     3    r*(9.145228263963672d3*r-2.120030188464306d4)+1.496491897739        
     4    51d4)-2.909845356715714d3)+r*(r*(r*(5.715767        
     5    664977295d2*r-1.870614872174388d3)+2.182384017536786d3)-1.0392        
     6    30484541326d3)+1.55884572681199d2                                     
       b_(21,i) = s*(s*(s*(7.238148934637916d3*s+1.28        
     1    6782032824518d4*r-1.286782032824518d4)+r*(7.238148934637        
     2    916d3*r-1.447629786927583d4)+7.238148934637916d3)+r*(1.0        
     3    d0*r*(1.378695035169127d3*r-4.136085105507381d3)+4.13608510550        
     4    7381d3)-1.378695035169127d3)+r*(r*(r*(5.7445        
     5    62646538029d1*r-2.297825058615211d2)+3.446737587922817d2)-2.29        
     6    7825058615211d2)+5.744562646538029d1                                  




         endif
      end do
      end subroutine       
      







    
