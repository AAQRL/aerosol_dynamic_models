c--------------------------------------------------------------------c
c     ALL RIGHT RESERVED. USE AND DISTRIBUTE UNDER PERMISSION ONLY.  c
c              (JIANHANC@UCI.EDU, MANDELSH@UCI.EDU)                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Collections of subroutines for calculating U matrices
c---------------------------------------------------------------------
c     Author: Jianhan Chen (jianhanc@uci.edu)
c
c     Purpuse: Given the signal, window range, basis information,
c     construct the U matrices and C vector. (see refences in main)
c
c     Note: Actually only 2-Scale FDM. Uncomment corresponding lines
c     to activate the true multi-scale FDM.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     UMatrix1D_MS0: implementation based on the MS-FDM paper         c
c     Author: Jianhan Chen (jianhanc@uci.edu)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine UMatrix1D_MS0(coef, Nsig, delt, Nb0, Nbc, Nb,
     &     wminl, wmaxl, gsave, U)
      implicit none
      integer i,j,n,k,l,Ml,Nsig,Nb,Nbc,Nb0,M(Nb)
      complex*16 U(Nb,Nb,0:1),coef(0:Nsig),f(Nb),g(Nb),z(Nb),wk(Nb),
     &     gh(Nb,Nb),gsave(Nb),cc,z1,z2
      real*8 delt, wminl, wmaxl, wmin, wmax, dw0, pi
c
      pi=dacos(-1d0)
      Ml=(Nsig-2)/2
      wmin=wminl*delt
      wmax=wmaxl*delt
      Nbc = Nbc + mod(Nbc,2) ! Nb0, Nbc have to be even numbers
      Nb0 = Nb0 + mod(Nb0,2)
      Nb = Nbc + Nb0
c
c     generate the basis
c
      dw0=(wmax-wmin)/Nb0
      do n=1, Nb0
         z(n+Nbc/2)=wmin+(n-0.5)*dw0
         M(n+Nbc/2)=Ml
      enddo
c      do n=Nbc/2, 1, -1  !!! uncomment to generate Multi-scale Basis
c         t=dfloat(2*n-1)/dfloat(Nbc)
c         rhoc=dsin(pi*t/2.0)
c         rhoc=(dsin(pi*t/2.0))**2 !!! replace with any distribution fxn
c         rhoc=dfloat(Nbc)/dfloat(Ml-Nb0)
c         M(n)=Ml*rhoc
c         if(M(n).lt.4) M(n)=4
c         M(Nb-n+1)=M(n)
c         z(n)=z(n+1)-(dw0*Ml)/dfloat(M(n))
c         z(Nb-n+1)=z(Nb-n)+(dw0*Ml)/dfloat(M(n))
c      enddo
      do n=1, Nbc !!! comment this part to disable double-scale basis
         if(n.le.Nbc/2) then
            M(n)=Nbc-1
            z(n)=-pi+(n-0.5)*2*pi/Nbc
         else
            M(n+Nb0)=Nbc-1
            z(n+Nb0)=-pi+(n-0.5)*2*pi/Nbc
         endif
      enddo
c
      do n=1, Nbc/2    ! REMOVE ANY OVERLAP BETWEEN WINDOW AND COARSE BASIS
         if(dreal(z(n)).lt.-pi) then
            k=(z(n)+pi)/(2*pi)
            z(n)=(z(n)-k*2*pi)+2*pi
         endif
         i=1
         do while (dabs(dreal(z(n)-z(i+Nbc/2))).gt.1d-8.and.i.le.Nb0)
            i=i+1
         enddo
         if(dabs(dreal(z(n)-z(i+Nbc/2))).le.1d-8.and.i.le.Nb0)
     &        z(n)=z(n)+5d-2*dw0
c
         if(dreal(z(Nb-n+1)).gt.pi) then
            k=(z(Nb-n+1)-pi)/(2*pi)
            z(Nb-n+1)=(z(Nb-n+1)-k*2*pi)-2*pi
         endif
         i=1
         do while
     &        (dabs(dreal(z(Nb-n+1)-z(i+Nbc/2))).gt.1d-8.and.i.le.Nb0)
            i=i+1
         enddo
         if(dabs(dreal(z(Nb-n+1)-z(i+Nbc/2))).le.1d-8.and.i.le.Nb0)
     &        z(Nb-n+1)=z(Nb-n+1)+5d-2*dw0
      enddo
      do n=1, Nb
         write(11,*) n, dreal(z(n))/delt, M(n)/dfloat(Ml), 1
         z(n)=cdexp((0,1d0)*z(n))
      enddo
      close(11)
c
c     basis generated; calculate the U matrices
c
      call HybridFT(Ml, M, Nb, coef, z, f, g)
c
c     call the matrix elements
c
      do l=0, 1
         if(l.eq.0) then
            do i=1, Nb
               wk(i)=f(i)
            enddo
            call diagME(Nb, M, coef, z, U)
            call rotFTgh(coef,Nb,Nb0,Nbc,Ml,M,z,g,gh)
         else
            do i=1, Nb
               z1=z(i)**M(i)
               U(i,i,l)=(U(i,i,l-1)-coef(l-1)-wk(i)+gh(i,i))/z(i)
     &              +coef(2*M(i)+l)*z1*z1
               wk(i)=wk(i)/z(i)-coef(l)+coef(M(i)+l)*z1
               do j=1, Nb/2
                  z2=z(i)**M(j)
                  gh(j,i)=gh(j,i)/z(i)-coef(M(j)+l)*z2
     &                 +coef(M(i)+M(j)+l)*z1*z2
                  gh(Nb-j+1,i)=gh(j,i)
               enddo
            enddo
         endif
         do i=1, Nb
            do j=1, i-1
               cc=z(j)/z(i)
               U(i,j,l)=coef(l)+(wk(i)-cc*wk(j)+
     &              cc**(-M(i))*gh(i,j)-cc**(M(j)+1)*gh(j,i))
     &              /(1-cc)
               U(j,i,l)=U(i,j,l)
            enddo
         enddo
      enddo
c
      do i=1, Nb                ! generate correct array C
         gsave(i)=f(i)+coef(0)
      enddo
      return
      end
c_____________________________________________________________________
c     END OF "UMatrix1D_MS0" SUBROUTINE
c----------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Name: Umatrix1D_dMS0
c
c     A routine using 'old' MS-FDM to construct U_Nskip & U_Nskip+1 for
c     specified window, and calculate gsave_0. Can be used to do linear
c     phase correction of delay = Nskip * tau by solving:
c
c       U_Nskip+1 Bk = uk U_Nskip Bk
c       dk^(1/2) = gsave' * Bk
c
c     Note: when called with Nskip=0, it will be exactly the same as
c           subroutine Umatrix1D_MS0().
c
c     Author: Jianhan Chen (jianhanc@uci.edu)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Umatrix1D_dMS0(coef,Nsig,Nskip,delt,Nb0,Nbc,Nb,
     &     wminl,wmaxl,gsave,U)
      implicit none
      integer i,j,n,k,l,Ml,Nsig,Nskip,Nb,Nbc,Nb0,M(Nb)
      complex*16 U(Nb,Nb,0:1), coef(0:Nsig+Nskip), gsave(Nb),
     &     wk(Nb), f(Nb), g(Nb), z(Nb), gh(Nb,Nb), cc, z1, z2
      real*8 delt, wminl, wmaxl, wmin, wmax, dw0, pi

      pi=dacos(-1d0)
      Ml=(Nsig-2)/2
      wmin=wminl*delt
      wmax=wmaxl*delt
      Nbc = Nbc + mod(Nbc,2) ! Nb0, Nbc have to be even numbers
      Nb0 = Nb0 + mod(Nb0,2)
      Nb = Nbc + Nb0
c
c     generate the basis
c
      dw0=(wmax-wmin)/Nb0
      do n=1, Nb0
         z(n+Nbc/2)=wmin+(n-0.5)*dw0
         M(n+Nbc/2)=Ml
      enddo
c      do n=Nbc/2, 1, -1  !!! uncomment to generate Multi-scale Basis
c         t=dfloat(2*n-1)/dfloat(Nbc)
c         rhoc=dsin(pi*t/2.0)
c         rhoc=(dsin(pi*t/2.0))**2 !!! replace with any distribution fxn
c         rhoc=dfloat(Nbc)/dfloat(Ml-Nb0)
c         M(n)=Ml*rhoc
c         if(M(n).lt.4) M(n)=4
c         M(Nb-n+1)=M(n)
c         z(n)=z(n+1)-(dw0*Ml)/dfloat(M(n))
c         z(Nb-n+1)=z(Nb-n)+(dw0*Ml)/dfloat(M(n))
c      enddo
      do n=1, Nbc !!! comment this part to disable double-scale basis
         if(n.le.Nbc/2) then
            M(n)=Nbc-1
            z(n)=-pi+(n-0.5)*2*pi/Nbc
         else
            M(n+Nb0)=Nbc-1
            z(n+Nb0)=-pi+(n-0.5)*2*pi/Nbc
         endif
      enddo

      do n=1, Nbc/2             ! check and correct basis positions
         if(dreal(z(n)).lt.-pi) then
            k=(z(n)+pi)/(2*pi)
            z(n)=(z(n)-k*2*pi)+2*pi
         endif
         i=1
         do while (dabs(dreal(z(n)-z(i+Nbc/2))).gt.1d-8.and.i.le.Nb0)
            i=i+1
         enddo
         if(dabs(dreal(z(n)-z(i+Nbc/2))).le.1d-8.and.i.le.Nb0)
     &        z(n)=z(n)+5d-2*dw0

         if(dreal(z(Nb-n+1)).gt.pi) then
            k=(z(Nb-n+1)-pi)/(2*pi)
            z(Nb-n+1)=(z(Nb-n+1)-k*2*pi)-2*pi
         endif
         i=1
         do while
     &        (dabs(dreal(z(Nb-n+1)-z(i+Nbc/2))).gt.1d-8.and.i.le.Nb0)
            i=i+1
         enddo
         if(dabs(dreal(z(Nb-n+1)-z(i+Nbc/2))).le.1d-8.and.i.le.Nb0)
     &        z(Nb-n+1)=z(Nb-n+1)+5d-2*dw0
      enddo
      do n=1, Nb
         write(11,*) n, dreal(z(n))/delt, M(n)/dfloat(Ml), 1
         z(n)=cdexp((0,1d0)*z(n))
      enddo
      close(11)
c
c     basis generated; calculate the U matrices
c
      call HybridFT(Ml, M, Nb, coef(Nskip), z, f, g)
c
c     GENERATE THE MATRICES U_Nskip & U_Nskip+1 (U1, U0 if Nskip=0)
c
      do l=0, 2
         if(l.eq.0) then
            do i=1, Nb
               wk(i)=f(i)
            enddo
            call diagME(Nb, M, coef(Nskip), z, U)
            call rotFTgh(coef(Nskip),Nb,Nb0,Nbc,Ml,M,z,g,gh)
         else
            do i=1, Nb
               z1=z(i)**M(i)
               U(i,i,l)=(U(i,i,l-1)-coef(Nskip+l-1)-wk(i)+gh(i,i))/z(i)
     &              +coef(Nskip+2*M(i)+l)*z1*z1
               wk(i)=wk(i)/z(i)-coef(Nskip+l)+coef(Nskip+M(i)+l)*z1
               do j=1, Nb/2
                  z2=z(i)**M(j)
                  gh(j,i)=gh(j,i)/z(i)-coef(Nskip+M(j)+l)*z2
     &                 +coef(Nskip+M(i)+M(j)+l)*z1*z2
                  gh(Nb-j+1,i)=gh(j,i)
               enddo
            enddo
         endif
c
         do i=1, Nb
            do j=1, i-1
               cc=z(j)/z(i)
               U(i,j,l)=coef(Nskip+l)+(wk(i)-cc*wk(j)+
     &              cc**(-M(i))*gh(i,j)-cc**(M(j)+1)*gh(j,i))
     &              /(1-cc)
               U(j,i,l)=U(i,j,l)
            enddo
         enddo
      enddo
c
c     calculate gsave_Nskip and obtain gsave_0 by back rotating
c
      do i=1, Nb                ! actually gsave_Nskip (vector C)
         gsave(i)=f(i)+coef(Nskip)
      enddo
      do n=Nskip-1,0,-1         ! back-rotate gsave_Nskip to get gsave_0
         do i=1, Nb
            gsave(i)=gsave(i)*z(i)-z(i)**(M(i)+1)*coef(M(i)+n+1)+coef(n)
         enddo
      enddo
      return
      end
c________________________________________________________________________
c     END OF  UMatrix1D_dMS0
c------------------------------------------------------------------------
c     Exact the same purpose as UMatrix1D_dMS0(), but implemented using
c     a more efficient way: instead of using various Mj to multi-scale
c     basis, complex frequency grids are used, Im(zj) <-> Mj.
c
c     Author: Jianhan Chen (jianhanc@uci.edu)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Umatrix1D_dMS1(coef,Nsig,Nskip,delt,Nb0,Nbc,Nb,
     &     wminl,wmaxl,gsave,U)
      implicit none
      integer Nb,i,j,n,M,Nsig,Nskip,l,Nb0,Nbc,Mc
      complex*16 U(Nb,Nb,0:1),coef(0:Nsig),f(Nb),g(Nb),diag(Nb),z(Nb),
     &     Z1,Z2,cc,xi,gsave(Nb),wk(Nb),dk(Nb)
      real*8 delt,wminl,wmaxl,wminl_new,wmaxl_new,pi,z_re,z_im, dw,
     &     r1,rhoc,dw0
      logical MULTI_SCALE
c      MULTI_SCALE=.true.
      MULTI_SCALE=.false.
c
c     NOTE:  set MULTI_SCALE = .false. if only 2-scale needed
c     NOTE:  for multi-scale case, diffierent rhoc(j) distribution
c            can be chosen by inserting new formula for rhoc(j) at line 180
c
      pi=dacos(-1d0)
      xi=(0,1d0)
      M=(Nsig-2)/2
      wminl_new=wminl*delt
      wmaxl_new=wmaxl*delt
      Nbc = Nbc + mod(Nbc,2) ! Nbc have to be even numbers
      Nb = Nbc + Nb0
      dw0=(wmaxl_new-wminl_new)/Nb0
c      z_im=dw0
      z_im=0d0
      do i=1,Nb0
         z_re=wminl_new+(i-0.5d0)*dw0
         z(i+Nbc/2)=dcmplx(z_re,z_im)
      enddo
      if(.not.MULTI_SCALE.and.Nbc.gt.0) then !     two-scale
c        write(*,*) 'generate double-scale basis ...'
         r1=dfloat(Nbc)/dfloat(M+1)
         call Im_Mj(r1,z_im)
         z_im=z_im*delt
         if(Nbc.ne.0) dw=(2d0*pi-(wmaxl_new-wminl_new))/Nbc
         do i=1, Nbc
            if(i.le.Nbc/2) then
               z(i)=dcmplx(-pi+(i-0.5)*2*pi/Nbc,z_im)
            else
               z(i+Nb0)=dcmplx(-pi+(i-0.5)*2*pi/Nbc,z_im)
            endif
         enddo
      else if(Nbc.gt.0) then    ! multi-scale coarse basis
c         write(*,*) 'generate mulit-scale basis ...'
         do n=Nbc/2, 1, -1
            r1=dfloat(2*n-1)/dfloat(Nbc)
c 180        rhoc=dsin(pi*r1/2.0)
            rhoc=(dsin(pi*r1/2.0))**2
c      rhoc=dfloat(Nbc)/dfloat(Ml-Nb0)
c      rhoc=dsin(pi*r1/2.0)**2*(3.0/Mmin)
            Mc=M*rhoc
            if(Mc.lt.4) Mc=4
            rhoc=dfloat(Mc)/dfloat(M)
            call Im_Mj(rhoc,z_im)
            z_im=z_im*delt
            z_re=dreal(z(n+1))-dw0/rhoc
            z(n)=dcmplx(z_re,z_im)
            z_re=dreal(z(Nb-n))+dw0/rhoc
            z(Nb-n+1)=dcmplx(z_re,z_im)
         enddo
c         do n=1, Nbc/2          ! fold back into [-pi,pi]; unnecessary
c            if(dreal(z(n)).lt.-pi) then
c               k=(dreal(z(n))+pi)/(2*pi)
c               z(n)=dcmplx(dreal(z(n))-(k-1)*2*pi,dimag(z(n)))
c            endif
c            if(dreal(z(Nb-n+1)).gt.pi) then
c               k=(z(Nb-n+1)-pi)/(2*pi)
c               z(Nb-n+1)=dcmplx(dreal(z(Nb-n+1))-(k+1)*2*pi,
c     &              dimag(z(Nb-n+1)))
c            endif
c         enddo
      endif
      do n=1, Nb
         write(11,*) n, dreal(z(n))/delt, dimag(z(n))/delt, 1
         z(n)=cdexp(z(n)*xi)
      enddo
      close(11)
c
c      write(*,*) 'calculating the 1D FT series ...'
      call arrays(Nsig,M,coef(Nskip),Nb,z,f,g,diag)
c
c      write(*,*) 'calculating the U matrices ...'
      do l=0,1
         do i=1,Nb
            Z1=z(i)**M
            Z2=Z1**2
            if(l.eq.0) then
               U(i,i,0)=diag(i)
               wk(i)=f(i)
               dk(i)=g(i)
            else
               U(i,i,l)=( U(i,i,l-1)-coef(Nskip+l-1)
     &              - wk(i) + dk(i)*Z1 ) / z(i)
     &              + coef(Nskip+2*M+l) * Z2
               wk(i)=wk(i)/z(i)-coef(Nskip+l)+coef(Nskip+M+l)*Z1
               dk(i)=dk(i)/z(i)-coef(Nskip+M+l)+coef(Nskip+2*M+l)*Z1
            endif
         enddo
         do i=1,Nb
            Z1=z(i)**M
            do j=1,i-1
               Z2=z(j)**M
               cc=z(i)/z(j)
               U(i,j,l)=coef(Nskip+l)
     &              +(wk(j)-cc*Z1*dk(j)
     &              -cc*wk(i)+Z2*dk(i))/(1-cc)
               U(j,i,l)=U(i,j,l)
            enddo
         enddo
      enddo
      do i=1, Nb                ! actually gsave_Nskip
         gsave(i)=f(i)+coef(Nskip)
      enddo
      do n=Nskip-1,0,-1         ! back-rotate gsave_Nskip to get gsave_0
         do i=1, Nb
            gsave(i)=gsave(i)*z(i)-z(i)**(M+1)*coef(M+n+1)+coef(n)
         enddo
      enddo
c      write(*,*) 'UMatrix1D_dMS1 is done'
      return
      end
c__________________________________________________________________________
c     END OF UMatrix_dMS1()
c--------------------------------------------------------------------------
c
c     calculate FT f=sum_1^Mj() and g=sum_(Ml+1)^(Ml+Mj)()
c
      subroutine HybridFT(Ml,M,Nb,coef,z,f,g)
      implicit none
      integer Ml, Nb, M(Nb), i, n
      complex*16 z(Nb), f(Nb), g(Nb), coef(0:100), c1

      do i=1, Nb
         f(i)=(0,0d0)
         g(i)=(0,0d0)
         do n=1, M(i)
            if(mod(n, 100).eq.1) then
               c1=z(i)**n
            else
               c1=c1*z(i)
            endif
            f(i)=f(i)+coef(n)*c1
            g(i)=g(i)+coef(n+Ml)*c1
         enddo
         g(i)=g(i)*z(i)**Ml
      enddo
      return
      end
c
c     'rotate' g(j) to get gh(i,j)=sum_(Mi+1)^(Mi+Mj)()
c
      subroutine rotFTgh(coef,Nb,Nb0,Nbc,Ml,M,z,g,gh)
      implicit none
      integer Ml, Nb, Nb0, Nbc, M(Nb), i, j, n, Mpre
      complex*16 coef(0:100), z(Nbc), g(Nb), gh(Nb, Nb),
     &     c1, c2, gpre

      do i=1, Nb
         gpre=g(i)
         Mpre=Ml
         c1=z(i)**M(i)
         do j=Nb/2, 1, -1
            do n=Mpre, M(j)+1, -1
               c2=z(i)**n
               gpre=gpre-coef(M(i)+n)*c1*c2+coef(n)*c2
            enddo
            gh(j,i)=gpre
            gh(Nb-j+1,i)=gpre
            Mpre=M(j)
         enddo
      enddo
      return
      end
c
c     a subroutine calculating the diagonal matrix elements
c
      subroutine diagME(Nb, M, coef, z, U)
      implicit none
      integer Nb, M(Nb), i, n
      complex*16 coef(0:100), z(Nb), U(Nb,Nb), cc, tmp1, tmp2

      do n=1, Nb
         tmp1=(0,0d0)
         tmp2=(0,0d0)
         do i=1, M(n)
            if(mod(i,100).eq.1) then
               cc=z(n)**i
            else
               cc=cc*z(n)
            endif
            tmp1=tmp1+coef(i)*(i+1)*cc
            tmp2=tmp2+coef(i+M(n))*(M(n)+1-i)*cc
         enddo
         U(n,n)=coef(0)+tmp1+tmp2*cc
      enddo
      return
      end
c
c     modified by J.Chen at 10-21-2001. now calculate:
c
c        f(i) = sum_i=1^M z(i)**i * c(i)
c        g(i) = sum_i=1^M z(i)**i * c(i+M)
c
      subroutine arrays(Nsig,M,coef,Nz,z,f,g,diag)
      implicit none
      integer Nz,i,n,M,Nsig
      complex*16 coef(0:Nsig),f(Nz),g(Nz),diag(Nz),z(Nz),coef1(M),Z1
c      write(6,*) 'construction of 1D arrays for the U-matrices'
      do n=1,M
         coef1(n)=(n+1)*coef(n)
      enddo
      call Slow_FT(M,coef1(1),Nz,z,diag)
      do n=1,M
         coef1(n)=(n+1)*coef(n+M)
      enddo
      call Slow_FT(M,coef1(1),Nz,z,f)
      call Slow_FT(M,coef(M+1),Nz,z,g)
      do i=1,Nz
         Z1=z(i)**M
c        g(i)=g(i)*z(i)**M        ! bad !!!
         diag(i)=diag(i)+coef(0)+((M+2)*g(i)-f(i))*Z1
      enddo
      call Slow_FT(M,coef(1),Nz,z,f)
      return
      end
c
c     More accurate slow FT
      subroutine Slow_FT(M,coef,Nz,z,f)
      implicit none
      integer M,n,Nz,i
      complex*16 coef(M),z(Nz),f(Nz),fff
      do i=1,Nz
         f(i) = (0,0d0)
         fff = coef(M)*z(i)
         do n=M-1,1,-1
            fff = (fff+coef(n))*z(i)
            if(mod(n,100).eq.0.or.n.eq.1) then
               f(i)=f(i)+fff*z(i)**(n-1)
               fff=(0,0d0)
            endif
         enddo
      enddo
      return
      end
c
c     Given t=(Mj+1)/(M+1), calculate z_im which satisfies:
c        sum_{n=0}^M exp(-z_im*n)= Mj+1
c
      subroutine Im_Mj(t, Im)
      implicit none
      real*8 t, x, Im
      if(t.le.0.or.t.gt.1.0) then
         write(*,1999) t
 1999    format('WARNING: t=',E13.7,'; t be between 0 and 1.0')
         t=1d-2
      endif
      if(t.le.0.685) then
         x=0d0
         if(t.gt.1d-2) x=-dexp(-1d0/t)/t
         Im=1d0/t+x-x**2+1.5d0*x**3-8d0*x**4/3d0+125d0*x**5/24d0
      else
         Im=2.82d0*(1d0-t)
      endif
      return
      end
