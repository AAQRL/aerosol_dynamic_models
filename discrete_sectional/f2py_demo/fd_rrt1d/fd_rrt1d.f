c                    ALL RIGHTS RESERVED!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        program FDM_RRT_1D_FinalVersion_CJH
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1. NOTE:    USE OR DISTRIBUTE WITH PERMISSION ONLY.
c
c     2. NAME:    fd_rrt1d.f
c
c     3. AUTHORS: Jianhan Chen (jianhanc@uci.edu)
c                 V. A. Mandelshtam (mandelsh@uci.edu)
c
c     4. DATE:    10-10-2001
c
c     5. PURPOSE: 1D FDM/RRT/DFT
c     
C              FDM = Filter Diagonalization Method
c              RRT = Regularized Resolvent Transform 
c              DFT = Discrete Fourier Transform
C
c     6. GENERAL DESCRIPTION: 
c
c      The program is based on the FDM and RRT algorithms for high resolution 
c      spectral analysis of time domain signals. See references for details.
c      Given a signal c(n) = C(t0+tau*n) n=0, 1, ..., Nsig, calculate the 
c      line list {w_k, d_k} and/or the complex spectrum F(w) by fitting c(n) 
c      to the Auto-Regressive (AR) model,
c
c                 c(n) = sum_k d_k exp(-i n w_k)        (1)
c                 
c      where, w_k is the complex frequency (line position & width)
c            d_k is the complex amplitude 30482(phase & intregral)
c     
c      1) in FDM, the complex spectrum is computed using the formula,
c     
c                 F(w) = sum_k d_k/(w-w_k+i*Gamm)       (2)
c
c         in RRT, it is estimated as,
c
c                 F(w) ~ Cv 1/R(w) Cv'                    (3)
c     
c      2). QZ-algorithm is applied to solve the generalized eigenvalue 
c      problem. ZGESV routine from LAPACK is used to invert matrix R(w).
c
c     7. INPUTS:
C
C      IT IS IMPORTANT THAT YOU UNDERSTAND SOME BASIC IDEAS OF FDM/RRT
C      IN ORDER TO UNDERSTAND THE MEANING OF SOME INPUTS SUCH AS BASIS 
C      SIZE, BASIS DENSITY. IN SUMMARY, FDM/RRT PROCESS THE TIME SIGNAL
C      BY DIVIDING THE WHOLE NYQUIST RANGE INTO ONE OR SEVERAL LOCAL
C      WINDOWS IN THE FREQUENCY DOMAIN. 
c
c      signal: input data file with appropriate first line. The data 
c              file is assumed to contain an ASCII formatted list of
c              c(n\tau). The first line must be as following:
c
c               dim { NSigmax tau }_[1,..,dim] {idat}_[1,...,dim]
c              
c              where dim    : the dimension of the signal (dim=1, 2, ...)
c                    Nsigmax: the maximum signal length
c                    tau    : time step (acquistion inteval)
c                    idat   : data format along each dimension
c               e.g/ 1D example: 1 1024 0.001 1 
c                    2D example: 2 2024 0.001 128 0.0001 1 1
c              
c              Note: 'idat' implies the data tormat. In 1D, it is defined
c              as following:
c                     idat>0 read real signal, idat<0 read complex signal 
c                     |idat| = 1 read c_n
c                     |idat| = 2 read t_n,c_n
c
c      t0    : 1st order (i.e., linear) phase correction (in seconds) 
c              NOTE: t0 is not used in this code.
c      theta : zero order (overall) phase correction
c      method: method used to compute the spectra: FDM, RRT, or DFT
c      Nsig  : number of data points to be used in computing the spectra.
c              It will be automatically reduced to Nsigmax (maximum length)
c              if specified Nsig is greater than Nsigmax.
c      wmino : (see wmaxo)
c      wmaxo : specify the frequency range of the spectra to be computed. 
c              It might be adjusted to fit harmonically with the basis size
c              specified (see Nb0, rho) and to fit in the Nyquist Range.
c
c      threshhold: see par below.
c      par   : file where the linelist will be written (FDM only).
c              Only those with |d(k)| > threshhold will be output.
c      ReSp  : file where Re[F(w)] will be written 
c      ImSp  : file where Im[F(w)] will be written
c      AbsSp : file where Abs[F(w)] will be written
c              NOTE: if file name is specified as 'none' or 'None', no
c                    corresponding output will be written.

c      rho   : basis density (default = 1.1). rho > 1.0 means putting 
c              more basis functions than default. For example, rho=1.1
c              means putting 10% more basis functions then default. 
c              rho should not be less than 1.0. rho = 1.5 should be
c              considered as the maximum. For typical data, rho = 1.1-1.2
c              will usually work the best. if rho < 0, default is used.
c      Nb0   : number of narrow band Fourier basis functions used per 
c              window (see references). Typical value is about 50-300.
c              Large Nb0 means bigger spectral range per window. The result
c              might be more stable but the calculation will be slower for
c              larger Nb0. Except for very long signal (say, Nsig > 10,000),
c              Nb0 = 150 will be a good choice. 
c      Nbc   : number of broad band Fourier basis functions (coarse basis).
c              See second reference for more details. Nbc <= 0 means no
c              coarse basis. Adding coarse basis (Nbc>0) will usually give
c              stabler results and more accurate linelist, especially when
c              the signal is noisy or contains non-localized features such
c              as broad lines, backgrounds. Typical Nbc = 10 ~ 50. A good
c              choice could be Nbc = 20.
c
c      Nsp   : number of points used in plotting the spectra. Defines the
c              digital resolution of output spectra.
c      Gamm  : smoothing parameter. Defines the smallest possible linewidth.
c              Default Gamm = 0.2 * 1/(N*tau). FDM/RRT might give some very
c              narrow lines (spikes) due to noise. A non-zero Gamm will 
c              improve the looking of resulted spectra. if Gamm < 0, default
c              will be used.
c      cheat : multiply all widths by cheat (FDM only). 
c              NOTE: be very careful to use cheat < 1.0. It might lead to 
c                    very misleading results. Use cheat=1.0 only.
c      cheatmore: if is .true. F(w) is computed with Im d_k (FDM only).
c              NOTE: be very careful again. Default: cheatmore = .false.
c
c      ros   : regularization parameter (RRT only). If ros < = 0, no 
c              regularization in RRT. See 3rd reference for more details.
c              In general, bigger ros means more 'regularization', and will
c              lead to more stable result but will lower resolution.
c              Qualitively, ros is related to the amount of noise present
c              in the signal. Noisy data requires bigger regularization.
c              For typical NMR signals with fine S/N (say, 500M magnet,
c              1 mM sample, normal probe, 4-16 scans), ros = 1d-6 should
c              a good guess. If you see artifacts (very easy to tell if 
c              you know what a typical NMR peak will look like), increase
c              ros (say, double the value) and run again. If the spectrum
c              look very clean but the resolution is very poor, it may 
c              indicate too big regularization. Decrease ros and run again.
c              NOTE: Choosing an optimal value for ros is a try-and-error
c                    game. Typically, there will be a flat regime where
c                    the spectra remains similar even changing ros by order 
c                    of a magnitude. Try following:
c
c                              too many artifacts
c                  ros=1d-6---------------  ros = ros*10 = 1d-5 --> cont.
c                            | clean, but low res. 
c                            -------------  ros = ros/10 = 1d-7 --> cont.
c                            | GREAT! 
c                            --------- done!
c
c      *** SAMPLE INPUT FILE: TO PROCESS A SIGNAL STORED IN 'signal.txt'
c
c       'signal.txt'        	                /signal
c       1.57                     		/theta
c       1                               	/ispec
c       1000	                            	/Nsig
c       -10 10 	                  		/wmin wmax
c       'par'	 	                      	/parameters file
c       'ft','none','none' 	          	/ReSp,ImSp,AbsSp files
c       1., 200, -20                          	/rho, Nb0, Nbc
c       3000, 1d-2				/Nsp, Gamm
c       1 F                             	/cheat, cheatmore
c       5d-7                             	/ros
c
c     8. OUTPUTS: can be found in the spectra and parameter files specified
c        in the input file. When multi-scale is used, the position of basis
c        functions can be found in file 'fort.11'.
c
c     9. REFERENCES
c     
c      FDM/RRT:
c       H. Hu et. al., J. Magn. Reson. 134, 76-87 (1998).
c       J. Chen and V. A. Mandelshtam, J. Chem. Phys. 112, 4429-4437 (2000). 
c       J. Chen, et.al., J. Magn. Reson. 147, 129-137 (2000). 
c       V. A. Mandelshtam, Prog. Nuc. Magn. Reson. Spect. 38, 159-196 (2001). 
c      QZ: http://www.netlib.org/toms/535
c      ZGESV/LAPACK: http://www.netlib.org/lapack/index.html
c
c     Questions and comments should be directed to V. A. Mandelshtam.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer Nb,Nb0,Nbc,M,Nsig,idat,Nsig1max,NWin,ispec,Dim,Nz,Nsp,
     &     Nbmax,Nsigmax,Npowermax,iminl,imaxl,imeanl,nw,N0,
     &     i, j, k, n, di
      parameter (Nsigmax=99999,Npowermax=99999,Nbmax=1500)
      real*8 theta,cheat,pi,Omega,Gamm,delt,wmin,wmax,
     &     wmino,wmaxo,rho,dOmega,dW,wminl,wmaxl,wmeanl,
     &     tt,r1,r2,r3,ros, threshhold
      complex*16 wk(Nbmax),dk(Nbmax),phase_corr,Spec(0:Npowermax),
     &     cc,dd,uu,Z1,Z2,Z3,xi,coef(0:Nsigmax),Spec_win(0:Npowermax)
      complex*16 gsave(Nbmax), U(Nbmax*Nbmax*3), heap(4*Nbmax**2)
      character*30 par,signal,AbsSp,ReSp,ImSp,method
      logical cheatmore,checkfile

      write(*,001)
 001  format(/,'  1D FDM/RRT/DFT with Multi-Scale option, Beta version',
     &     /,  '       J. Chen and V. A. Mandelshtam, UC-Irvine',/)
c
c     READ IN THE PARAMTERS
c      
      read(5,*) signal           ! signal file
      read(5,*) theta            ! zero, order phase correction
      read(5,*) method            ! FDM, RRT or DFT
      read(5,*) Nsig             ! signal length to be used
      read(5,*) wmino,wmaxo      ! window size
      read(5,*) par, threshhold  ! linelist output file
      read(5,*) ReSp,ImSp,AbsSp  ! spectra output files
      read(5,*) rho, Nb0, Nbc    ! basis density, basis size, coarse basis
      read(5,*) Nsp, Gamm        ! plotting pnts, smoothing 
      read(5,*) cheat, cheatmore ! used in FDM, see discriptions.
      read(5,*) ros              ! RRT regularization parameter
c
c     CHECK THE PARAMETERS
c
      xi=(0d0,1d0)
      pi=dacos(-1d0)
      open(15,file=signal,status='old')
      write(6,*) 'Read the signal from: ', signal
      read(15,*) Dim, Nsig1max, delt, (i, r1, k=1,Dim-1), idat
      if(Dim.le.0.or.delt.eq.0.or.Nsig1max.le.0) 
     &     stop 'Error: invalid entries in 1st line of the signal!'
      if(Nsig.gt.Nsig1max-1) Nsig=Nsig1max-1
      if(Dim.ne.1) idat = 0     ! multi-D data, idat = 0
      if(Nsig.gt.Nsigmax) stop 'Error: Nsig > Maximum allowed!'
      r2=0d0
      r3=0d0                    ! r3 = sum_n |c(n)|^2
      phase_corr=cdexp(dcmplx(0d0,theta))
      do n=0,Nsig
         if(idat.eq.-1) then
            read(15,*,end=100) r1,r2
         else if(idat.eq.-2) then
            read(15,*,end=100) tt,r1,r2
         else if(idat.eq.1) then
            read(15,*,end=100) r1
         else if(idat.eq.2) then
            read(15,*,end=100) tt,r1
         else
            read(15,*,end=100) r1,r2 ! default format
         endif
         coef(n)=dcmplx(r1,r2)*phase_corr
         r3=r3+cdabs(coef(n))**2
      enddo
      goto 101
 100  Nsig=n-1
 101  M=(Nsig-2)/2
      Nsig=2*M+2
      write(6,9991) 2*M+3
 9991 format(' ', I5, ' c_n are read and to be used')
      ros=ros*r3*(M+1d0)**2     ! actual regularization parameter
c
      wmin=wmino*delt
      wmax=wmaxo*delt
      if(wmin.gt.wmax) then
         r1 = wmin
         wmin = wmax
         wmax = r1
      endif
      if(wmin.gt.pi.or.wmax.lt.-pi)
     &     stop 'Error: invalid spectral range, abort!'
      if(wmin.lt.-pi) wmin=-pi
      if(wmax.gt.pi) wmax=pi
      wmino=wmin/delt
      wmaxo=wmax/delt
      write(6,9992) wmino,wmaxo
 9992 format(' Spectral range to be processed: [',F10.3,', ',F10.3,']')
c      if(Gamm.lt.0) write(6,*) 'Warning: Gamm < 0,',
c     &     ' which means line narrowing instead of smoothing!'
      if(Gamm.le.0) Gamm = 0.4*pi/(Nsig*delt)
      write(*,*) 'Smoothing factor Gamm = ', Gamm
      di = 10*Nsig*(wmax-wmin)/(2*pi)
      Nsp = max(Nsp, di)
      if(Nsp>Npowermax) then
         write(6,*) 'Warning: Nsp is truncated to Npowermax (99999).'
         Nsp=Npowermax
      endif
      dOmega=(wmaxo-wmino)/Nsp
      do j=0,Nsp
         Spec(j)=(0d0,0)
      enddo
c
      write(6,*)
      if(method.eq.'DFT') then
         ispec = -1
         write(6,*) 'Discrete Fourier Transform'
         goto 246               ! DFT
      else if(method.eq.'NDFT') then
         ispec = -2
         write(6,*) 'Discrete Fourier Transform w/o apodization'
         goto 246               ! DFT
      else if(method.eq.'RRT') then
         ispec = 0
         write(6,*) 'Regularized Resolvent Transform (RRT)'
      else 
         ispec = 1
         write(6,*) 'Filter Diagonalization Method (FDM)'
      endif
c
      if(rho.le.0) rho=1.1      ! default basis density
      if(Nb0.lt.2) Nb0=2
      write(6,*) 'Use the frequency grid with rho = ',rho
      Nz=((wmax-wmin)/(2*pi))*(M+1)*rho
      if(Nb0.gt.Nz) Nb0=Nz
      NWin = (2*Nz)/Nb0         ! number of windows
      dW = (wmaxo-wmino)/NWin   ! window size
      Nb0 = dW*delt/pi * (M+1)*rho ! adjusted size of window basis
      Nb0=Nb0+mod(Nb0,2)
      write(6,9993) NWin-1, Nb0
 9993 format(' Divide the spectral range into ',I4,' windows with Nb0 ='
     &     , I4)
      if(Nbc.le.0) then
         Nbc=0
         write(6,*) 'No coarse basis used'
      else
         Nbc=Nbc+mod(Nbc,2)
         write(6,*) 'Add coarse basis with size Nbc = ',Nbc
      endif
      Nb=Nbc+Nb0
      write(6,9994) Nb
 9994 format(' The actual basis size per window is Nb = ',I4)
      if(Nb.gt.Nbmax) stop 'Error: Nb > Nbmax, abort! Reduce the size'
c
      if(ispec.eq.0) then
         write(6,*) 'Regularization paramter q^2 = ', ros         
      else if(ispec.gt.0) then
         if(cheat.ne.1d0) write(6,*) 'using cheat = ', cheat
         if(cheatmore) write(6,*) 'cheatmore = .true.,',
     &        ' |dk| is used to construct F(w)'
      endif
c_______________________________________________________________
c
      if(checkfile(par).and.ispec.gt.0) then
         threshhold = max(0d0, threshhold)
         open(7,file=par)
         write(7,8880) threshhold
 8880    format('# Parameter files generated by fd_rrt1d.f',/,
     &        '# Only lines with  |d(k)| > ',F10.7, ' are recorded',/
     &        '# w: complex frequencies, d: complex amplitudes',/,/,
     &        '#',7x,'Re w',18x,'Im w',15x,'Re d',14x,'Im d')
         write(6,8881) threshhold, par
 8881    format(' {wk,dk} for lines with |d(k)| > ',F8.6, 
     &        ' will be writen to : ',A14)
      endif
      write(6,*)
      wmeanl=wmino
      do nw=1, NWin-1
         wminl=wmeanl
         wmeanl=wminl+dW
         wmaxl=wmeanl+dW
         write(6,9995) nw, wminl, wmaxl
 9995    format('  window No.',I4,': [wmin,wmax]=[',F10.3,',',F10.3,']')
         iminl=(wminl-wmino)/dOmega
         imeanl=(wmeanl-wmino)/dOmega
         imaxl=(wmaxl-wmino)/dOmega
         do i=iminl, imaxl
            Spec_win(i) = dcmplx(0,0)
         enddo
c
         call UMatrix1D_MS0(coef, Nsig, delt, Nb0, Nbc, Nb,
     &        wminl, wmaxl, gsave, U)
c         call Umatrix1D_dMS0(coef,Nsig, 0, delt, Nb0, Nbc, Nb,
c     &        wminl, wmaxl, gsave, U)
c         call Umatrix1D_dMS1(coef,Nsig, 0, delt, Nb0, Nbc, Nb,
c     &        wminl, wmaxl, gsave, U)
c
         if(ispec.eq.0) then    ! RRT
            call ROS2K_1D(Nb,U,imaxl-iminl+1,wmino+iminl*dOmega,
     &         dOmega,delt,gsave,ros,heap,Spec_win(iminl),Gamm)
         else                   ! FDM
            call FD_1D(Nb, U, gsave, delt, wk, dk)
            do k=1,Nb           ! F(w) ~ sum_k {...}
               dd=dk(k)
               cc=wk(k)
               if(dabs(dimag(cc)).lt.Gamm) then
                  cc=dcmplx(dreal(cc),-Gamm)
               else
                  cc=dcmplx(dreal(cc),cheat*dimag(cc))
               endif                  
               do i=iminl, imaxl
                  Omega=wmino+i*dOmega
                  if(cheatmore) then
                     Spec_win(i)=Spec_win(i)-xi*delt*dd*
     &                    dreal(1/(1-cdexp(xi*delt*(Omega-cc))))
                  else
                     Spec_win(i)=Spec_win(i)-xi*delt*dd/
     &                    (1-cdexp(xi*delt*(Omega-cc)))
                  endif
               enddo
               if(checkfile(par)) then
                  if(dreal(wk(k)).ge.wminl.and.dreal(wk(k)).le.wmaxl
     &                 .and.cdabs(dk(k)).ge.threshhold) 
     &                 write(7,13) dreal(wk(k)),dimag(wk(k)),
     &                 dreal(dk(k)),dimag(dk(k))
               endif
 13            format(E22.16,3E18.10)
            enddo
         endif
c
c     add up the indivual contribution of each window
c
         if(iminl.eq.0) then
            do i=iminl, imeanl
               Spec(i)=Spec(i)+Spec_win(i)
            enddo
         else
            di=imeanl-iminl
            do i=iminl, imeanl
               Spec(i)=Spec(i)+Spec_win(i)*
     &              dsin(0.5*pi*dfloat(i-iminl)/float(di))**2
            enddo
         endif
         if(imaxl.ge.Nsp-1) then
            do i=imeanl+1, imaxl
               Spec(i)=Spec(i)+Spec_win(i)
            enddo
         else
            di=imaxl-imeanl
            do i=imeanl+1, imaxl
               Spec(i)=Spec(i)+Spec_win(i)*
     &              dsin(0.5*pi*dfloat(imaxl-i)/dfloat(di))**2
            enddo
         endif
         if(checkfile(par).and.ispec.gt.0) write(7,*)
 177     continue
      enddo
c      
 246  if(ispec.lt.0) then       ! DFT
         if(ispec.eq.-1) then
            N0=Nsig/10
            do n=0,Nsig
               if (n.gt.N0) then
                  r1 = dfloat(n-N0)/dfloat(Nsig+1-N0)
                  coef(n)=coef(n)*dexp(r1)*(1D0-r1)
               endif
            enddo
         endif
         do i=0,Nsp
            Omega=wmino+i*dOmega
            uu=cdexp((0,1d0)*Omega*delt)
            Z1=-delt*xi         !*cdexp(dcmplx(0d0,t0*Omega))
            Spec(i)=Spec(i)+coef(0)*Z1/2
            do n=1,Nsig
               Z1=Z1*uu
               Spec(i)=Spec(i)+coef(n)*Z1
            enddo
         enddo
      endif
c
      if(ispec.ne.-1) then      ! first point correction
         cc=-coef(0)*delt*xi/2.0
         do i=0, Nsp
            Spec(i)=Spec(i)-cc
         enddo
      endif
c
      write(6,*)                ! output spectra
      if(checkfile(ReSp)) then
         write(6,*) 'write Re[I(w)] to file: ', ReSp
         open(14,file=ReSp)
         do i=0,Nsp
            write(14,986) wmino+i*dOmega,dreal(Spec(i))
         enddo
         close(14)
      endif
      if(checkfile(ImSp)) then
         write(6,*) 'write Im[I(w)] to file: ', ImSp
         open(14,file=ImSp)
         do i=0,Nsp
            write(14,986) wmino+i*dOmega,dimag(Spec(i))
         enddo
         close(14)
      endif
      if(checkfile(AbsSp)) then
         write(6,*) 'write Abs[I(w)] to file: ', AbsSp
         open(14,file=AbsSp)
         do i=0,Nsp
            write(14,986) wmino+i*dOmega,cdabs(Spec(i))
         enddo
         close(14)
      endif
 986  format(E14.7,E12.4)
      write(6,*)
      write(6,*) 'End of program: successful excution'
      stop
      end
c      
      logical function checkfile(file)
      character*30 file
      checkfile=.true.
      if(file.eq.'none'.or.file.eq.'NONE') checkfile=.false.
      return
      end
c
c     FDM subroutine: give U and gsave, calculate dk, wk
c
      subroutine FD_1D(Nb,U,gsave,delt,wk,dk)
      implicit none
      integer Nb, k, j
      real*8 delt
      complex*16 U(Nb,Nb,0:1), gsave(Nb), dk(Nb), wk(Nb), ZZ(Nb,Nb)
c      
      call CQZ(Nb,U(1,1,1),U,ZZ,wk) ! solve GEP
      do k=1,Nb                 ! get the frequency in the correct units
         wk(k)=(0,1d0)*cdlog(wk(k))/delt    
      enddo
      do k=1,Nb                 ! compute the amplitudes
         dk(k)=(0,0d0)
         do j=1,Nb
            dk(k)=dk(k) + ZZ(j,k)*gsave(j)
         enddo
         dk(k) = dk(k)**2
      enddo
      call cpiksrt(Nb,wk,dk)    ! sort according to fequency
      return
      end
c
c     The frequencies are sorted in increasing order.
c
      subroutine cpiksrt(Nb,wr,wi)
      implicit real*8(a-h,o-z)
      complex*16 wr(Nb),wi(Nb),war,wai
      do j=2,Nb
         war=wr(j)
         wai=wi(j)
         do i=j-1,1,-1
            if(dreal(wr(i)).le.dreal(war))go to 10
            wr(i+1)=wr(i)
            wi(i+1)=wi(i)
         enddo
         i=0
 10      continue
         wr(i+1)=war
         wi(i+1)=wai
      enddo
      return
      end
c
c     calculate spectra according to Regularized Resolvent Transform
c     NOTE: first point correction is done in the main program
c
      subroutine ROS2K_1D(Nb, U, Nsp, wmin, dOmega, delt, 
     &     gsave, ros, Ull, spec, gamm)
      implicit none
      integer Nb, Nsp, i, j, n1, IPIV(Nb), INFO
      real*8 wmin, wmax, delt, ros, dOmega, Omega, gamm
      complex*16 U(Nb,Nb,0:1),gsave(Nb),B(Nb,Nb),C(Nb),A(Nb,Nb),
     &     ALPHA,BETA,cc,zz,xi,Ull(Nb,Nb,4),spec(Nsp),tmp(Nb)
c
      xi=(0d0,1d0)
      ALPHA=dcmplx(1d0,0d0)     ! needed in calling BLAS routines
      BETA=dcmplx(0d0,0d0)      ! needed in calling BLAS routines
c
c     calculate U00 = U0 * U0+,  ZGEMM is a BLAS level 3 routine
c
      call ZGEMM ( 'N', 'C', Nb, Nb, Nb, ALPHA, U(1,1,0), Nb, 
     &     U(1,1,0), Nb, BETA, Ull(1,1,1), Nb ) 
c
      if(ros.lt.0) then         ! just in case ...
         write(6,001) 
 001     format('WARNING: ros = ',E9.3,' should be non-negative, ',
     &        'thus set to zero')
         ros=0d0
      endif
c
c     calculate U0l=U0*Ul+, Ul0=U1*U0+, Ull=Ul*Ul+
c
      call ZGEMM ( 'N', 'C', Nb, Nb, Nb, ALPHA, U, Nb, U(1,1,1), Nb,
     &     BETA, Ull(1,1,2), Nb ) 
      call ZGEMM ( 'N', 'C', Nb, Nb, Nb, ALPHA, U(1,1,1), Nb, U, Nb,
     &     BETA, Ull(1,1,3), Nb ) 
      call ZGEMM ( 'N', 'C', Nb, Nb, Nb, ALPHA, U(1,1,1), Nb, 
     &     U(1,1,1), Nb, BETA, Ull(1,1,4), Nb ) 
c
      do n1=1, Nsp
         Omega=wmin+(n1-1)*dOmega
         zz=cdexp(-xi*(Omega+xi*gamm)*delt)
c
c     A = U - z*S ( where S is the overlap matrix )
c
         do i=1, Nb
            C(i)=gsave(i)
            do j=1, Nb
               A(i,j)=zz*U(i,j,0)-U(i,j,1)
            enddo
         enddo
c     
c     calculate B = (U-z*S) * (U-z*S)+ + ros, ros = q**2
c     
         cc=dcmplx(dreal(zz),-dimag(zz))
         do i=1, Nb
            do j=1, Nb
               B(i,j) = Ull(i,j,1)*zz*cc
     &              - Ull(i,j,2)*zz - Ull(i,j,3)*cc + Ull(i,j,4)
            enddo
         enddo
         if(ros.gt.0) then
            do i=1, Nb
               B(i,i) = B(i,i) + ros
            enddo
         endif
c
c     solve B * X = C, with result X stored in C;
c                      Using ZGESV or ZSYSV Lapack routines
c     
         call ZGESV( Nb, 1, B, Nb, IPIV, C, Nb, INFO )
         if(INFO.ne.0) then
            write(6,*) 'WARNING: ZGESV::INFO = ',INFO
            write(6,*) 'Error: ZGESV failed, unable to continue'
            stop
         endif
c     
c     calculate I(z) = C' * (U-z*S)+ * X ( X is in C here)
c            
         call ZGEMV ( 'C', Nb, Nb, ALPHA, A, Nb, C, 1,
     &        BETA, tmp, 1 )    !  ZGEMV is a BLAS level 2 routine
         cc=(0,0d0)
         do i=1, Nb
            cc=cc+gsave(i)*tmp(i)
         enddo
         spec(n1) = -cc*delt*zz*xi
      enddo
      return
      end
