c       To run these module you need to copy the two following
c       data files in the same directory: chisquare.dat and ftest.dat
c       Then change the variable name21 and name17 line 350 and 370 of
c       this file. Lastly compile it, executing the following line:
c       f2py -c -m mtm2 mtm2.f
c       Specify in your Pythonpath the path to the new mtm2 module built
c       and then import mtm2 and mtm_fct to enjoy the new mtm functions!
c
c
        subroutine mtm_mann(spec_raw,spec_resh,spec_smoo,
     $   spec_conf,recon_sig,time_serie,nscan,dt,npi,nwin,f1,f2,
     $   inorm,ispec,iresh,ithresh,inoise,ilog,ismooth,
     $  isignal,irecon,nsignals,iplotresh,iplotftest,iplotsmoo,
     $  iplotraw,iplotconf,icon)
c
c     driver for MTM spectral analysis/signal reconstruction code
c
c     (c) Michael E. Mann 12/96 modified November 2nd 2009 to be
c       operated as a python function
c
c
      parameter (maxlen=50000,zero=0.0,maxsignal=100)
      parameter (nlim=32768)
        integer nscan,npi,nwin,inorm,ispec,iresh,ithresh
        integer inoise,ilog,ismooth,isignal,irecon,nsignals,icon
        integer iplotresh,iplotftest,iplotsmoo,iplotraw,iplotconf
        real time_serie(nscan)
        real spec_raw(nlim,3),spec_resh(nlim,3)
        real spec_smoo(nlim,2),spec_conf(nlim,5)
        real recon_sig(nlim,3)
        real dt,f1,f2
cf2py   intent(in) time_serie
cf2py   intent(in) dt
cf2py   intent(in) f1
cf2py   intent(in) f2
cf2py   intent(in) nscan
cf2py   intent(hide) nscan
cf2py   intent(in) npi
cf2py   intent(in) nwin
cf2py   intent(in) inorm,
cf2py   intent(in) ispec
cf2py   intent(in) iresh
cf2py   intent(in) ithresh
cf2py   intent(in) inoise
cf2py   intent(in) ilog
cf2py   intent(in) ismooth
cf2py   intent(in) isignal
cf2py   intent(in) irecon
cf2py   intent(in) nsignals
cf2py   intent(in) iplotresh
cf2py   intent(in) iplotftest
cf2py   intent(in) iplotsmoo
cf2py   intent(in) iplotraw
cf2py   intent(in) iplotconf
cf2py   intent(in) icon
cf2py   intent(out) spec_raw
cf2py   intent(out) spec_smoo
cf2py   intent(out) spec_conf
cf2py   intent(out) spec_resh
cf2py   intent(out) recon_sig
        real a(maxlen),rcon(maxlen),signal(maxlen)
      real fsignal(maxsignal),confsignal(maxsignal)
      integer irec(maxsignal)
      real anpi,fthresh,f0
      integer nf0,nbnd
      real fconf(6),conf(4)
      real specraw0(nlim),specmed0(nlim),specresh0(nlim),
     $   harmonic0(nlim),specback0(nlim),ftest(nlim)
      real whiteraw0,whiterob0,rhoraw0,rhorob0,tauraw0,taurob0
c
c     read in data
c
        do i=1,nscan
         a(i) = time_serie(i)
      end do
        do i=nscan+1,maxlen
         a(i)=zero
      end do
        demn=0
      do i=1,nscan
         demn = demn+a(i)
      end do
      demn = demn/float(nscan)
c
c     calculate Rayleigh and Nyquist frequencies
c
      fray=1.0/(nscan*dt)
      fny = 0.5/dt
c
c
c     SET DEFAULTS
c
c     indicate that no reconstructions/plotting possible until
c     spectrum has been calculated...
c
c
c     resolution/variance tradeoff
c
      anpi=float(npi)
      fray=1.0/(nscan*dt)
      fny = 0.5/dt
      bndwdth = 2.0*anpi*fray
c
c     spectrum
c
      frange = f2-f1
      if (ithresh.eq.0) fthresh=50.0
      if (ithresh.eq.1) fthresh=90.0
        if (ithresh.eq.2) fthresh=95.0
        if (ithresh.eq.3) fthresh=99.0
        if (ithresh.eq.4) fthresh=99.5
        if (ithresh.eq.5) fthresh=99.9
c
c     null hypothesis
c
      fsmooth = frange/5.0
      if (fsmooth.lt.bndwdth) fsmooth=bndwdth
c
c     reconstruction
c
      do i=1,maxsignal
         irec(i)=0
         fsignal(i)=0.0
         confsignal(i)=0.0
      end do
c
c
c
       call spec(a,nwin,npi,dt,nscan,f1,f2,
     $     inoise,ismooth,ilog,fsmooth,
     $     isignal,ispec,inorm,iresh,ithresh,
     $     nf0,df,specraw0,specmed0,specresh0,harmonic0,
     $     specback0,ftest,
     $     conf,fconf,nbnd,
     $     whiteraw0,whiterob0,rhoraw0,rhorob0,tauraw0,taurob0)

       call findpeaks(flow,fhigh,isignal,
     $      nf0,df,nbnd,specraw0,specmed0,specresh0,harmonic0,
     $      specback0,ftest,
     $      conf,fconf,
     $      fsignal,confsignal,nsignals)

        call display(nf0,f1,df,
     $        spec_raw,spec_resh,spec_smoo,spec_conf,
     $     specraw0,specmed0,specresh0,harmonic0,
     $     specback0,ftest,conf,fconf,
     $     whiteraw0,whiterob0,rhoraw0,rhorob0,
     $     tauraw0,taurob0,
     $     iplotresh,iplotftest,iplotsmoo,
     $     iplotraw,iplotconf)

      if (irecon.eq.1) then
        do i=1,maxlen
           signal(i)=0.0
           rcon(i)=0.0
        end do
        do j=1,nsignals
         if (irec(j).eq.1) then
           f0 = fsignal(j)
           call recon(a,nwin,npi,dt,nscan,bndwdth,f0,icon,rcon)
           do i=1,nscan
              signal(i)=signal(i)+rcon(i)
           end do
         endif
        end do
          do i=1,nscan
            recon_sig(i,1)=float(i-1)*dt
                recon_sig(i,2)=signal(i)
                recon_sig(i,3)=a(i)-demn
        end do
        endif
c
      end
c
c
c
c Subroutine Spec
c =============================
      subroutine spec(a0,nwin0,ip,dt0,nscan0,
     $   flow,fhigh,
     $   inoise0,ismooth,ilog0,fsmooth,
     $   isignal,ispec,inorm,iresh,iper,
     $   nf0,df,specraw0,specmed0,specresh0,harmonic0,specback0,
     $   ftest,conf,fconf,nbnd,
     $   whiteraw0,whiterob0,rhoraw0,rhorob0,tauraw0,taurob0)
c
c     (c) Michael Mann
c
c     MTM procedure of Mann and Lees (1996) adapted fromm original MTM code
c     of J. Park, to perform multiple null-hypothesis testing, both harmonic
c     and quasiperiodic signal detection, and robust background noise
c     estimation.
c
      integer maxwin,maxlen,maxdof,nlim,nmed
      real big,small,tol
      parameter (maxlen=50000)
      parameter (maxwin=8,maxdof=2*maxwin,nlim=32768)
      parameter (nmed=2000)
      parameter (big=1.0e+36,small=1.0e-12,tol=1e-4)
c
        character*2000 name21, name17
      real dt0
      real a0(maxlen)
      integer ip,nscan0,nwin0,ilog0
      integer inoise,inoise0
      integer nbnd,nbnd0
      integer idone(nlim)
      real rho,rho0
      real MSE
      real*8 avar,el
      real val(nmed)
      real demean(nlim),specmed
      real specraw0(nlim)
      real specresh0(nlim),harmonic0(nlim),specback0(nlim),
     $      specmed0(nlim)
      real rednoise(nlim),rednoise0(nlim)
      real ftest(nlim)
      real base(nlim)
      integer iharm(nlim)
      real fconf(6)
      real conf(4)
      real adum,sumpeak
      complex zed
        real thresh
c
c     note that following limits are hard-wired into the code:
c       max # tapers = 8
c       max # frequency paris = 16400
c       max # points in dataseries = 32800
c
c
      common/taperz/ta(16400,8),tai(16400,8),dcf(16400,8),
     $   amu(16400,8)
      common/staple/tad(32800),tad1(32800)
      common/stap2/tas(32800,16)
      common/misc/npad
      common/npiprol/anpi
      common/data/a(maxlen)
      common/work/b(32800)
      common/work2/fspec(16400)
      common/BLK1/white,specmed(nlim),fnynew,ddf,f0
      common/BLK2/ilog,nfreq
c
      dimension dum(35)
      real *8 b1(32800),b2(32800)
      dimension reim(2),el(10)
      equivalence (zed,reim),(iy,dum)
c
c     define some constants
c
      pi=3.14159265358979
      radian=180./pi
      tiny = 1e-6
c
c     asign parameters from main
c
      nscan = nscan0
      do i=1,maxlen
         a(i)=a0(i)
      end do
      npi = ip
      anpi=float(ip)
      nwin=nwin0
      dt = dt0
      nscan = nscan0
      ilog = ilog0
c
c     zero pad to first power of 2 > # data points
c     constrain minimum padding to 1024 points...
c
      npts = nscan
      npad=npts-1
      ij=0
c
c     set zero padding to first power of 2 greater than
c     number of data points
c
c     impose minimum value of 1024
c
1000  npad=npad/2
      ij=ij+1
      if(npad.ne.0) go to 1000
c
      npad=max(2**ij,1024)
c
c     determine frequency limits, associated DFT points, etc.
c
      fmin=flow
      fmax=fhigh
      fnynew = fhigh
      ddf=1.0/(npad*dt)
      df = ddf
      eps=.1*ddf
      n1=(fmin+eps)/ddf
      n2=((fmax+eps)/ddf)
      fmin=n1*ddf
      fmax=n2*ddf
      nf=n2-n1+1
      nf0 = nf
      nfreq=nf
      f0 = fmin
      nmin=n1*2+1
      nmax=n2*2+2
c
c     determine Rayleigh, Nyquist frequencies, bandwidths, etc.
c
      fny = 0.5/dt
      fray=1.0/(npts*dt)
      bndwdth = 2.0*anpi*fray
      halfwdth = bndwdth/2.0
      nbnd = bndwdth/ddf
      nbnd0 = fray/ddf
      timefreq = 2.0*anpi
c
c     demean the series
c
      demn = 0.0
      do i=1,nscan
         demn=demn+a(i)
      end do
      demn = demn/float(nscan)
      do i=1,nscan
         demean(i)=a(i)-demn
      end do
c
c
c     determine raw noise parameters
c
c
      inoise = inoise0
      if (inoise.gt.0) then
         rho0 = 0.0
      else
c
c      determine raw lag 1 autocorrelation coefficient
c
       var = 0.0
       do i=1,nscan
         var = var + demean(i)**2
       end do
       var = var/float(nscan)
       sd = sqrt(var)
       c1 = 0.0
       icount = 0
       do i=2,nscan
         icount = icount + 1
         c1 = c1 + demean(i-1)*demean(i)
       end do
       c1 = c1/float(icount)
       rho0 = c1/var
      endif
c
c     determine chi-square values for confidence level
c     determination of spectrum
c
c      open (unit=21,file='/home/geovault-00/brel/Pythonlib/mtm_f2py/chisquare.dat',status='old')
      name21='/home/filipe/Dropbox/pymodules/python-ssa-mtm/mtm/chisquar
     +e.dat'
        open (unit=21,file=name21,status='old')
         do i=2,maxdof
            idofs = 2*nwin
            if (i.eq.idofs) then
               read (21,*) idum,conf(1),conf(2),conf(3),conf(4)
            else
               read (21,*)
            endif
         end do
      close (unit=21)
c
      iftest = 1
      jper = 3
c
c     read in f-test data if appropriate
c
      if (nwin.eq.1) iftest = 0
      if (iftest.eq.1) then
c       open (unit=17,file='/home/geovault-00/brel/Pythonlib/mtm_f2py/ftest.dat',status='unknown')
      name17='/home/filipe/Dropbox/pymodules/python-ssa-mtm/mtm/ftest.da
     +t'
        open (unit=17,file=name17,status='unknown')
       do i=2,14,2
         if (i.eq.2*nwin-2) then
            read (17,*) idum,fconf(1),fconf(2),fconf(3),fconf(4),
     $            fconf(5),fconf(6)
         else
            read (17,*)
         endif
       end do
       close (unit=17)
       thresh = fconf(iper)
      endif
c
c     determine amplitude theshold for reshaping
c
      afact = 0.0
      if ((iresh.eq.1).and.(isignal.ne.2)) then
         afact = conf(1)
         if (iper.eq.2) afact = conf(2)
         if (iper.eq.3) afact = conf(3)
         if (iper.gt.3) afact = conf(4)
      endif
c
c     create output files
c
c     DETERMINE SPECTRUM
c
c     calculate the eigentapers
c
      call taper2(npts,nwin,el)
c
c     normalization:
c
c     mult by dt if we wish absolute spectral estimate
c         e.g. analysis of time-limited signal
c     divide by npts if we wish amplitude spectrum per unit time
c
      if (inorm.eq.2) then  !  amp spec (time limited)
        anrm=1.0/dt
      else
        if(inorm.eq.1) then  ! amp spect per unit time
          anrm=float(npts)
        else
          anrm=1.
        endif
      endif
c
c     perform convolution of the series with each datataper
c
      do iwin=1,nwin
        do i=1,nscan
          b1(i)=demean(i)*tas(i,iwin)
          b2(i)=0.0
        end do
        j=nscan+1
        do i=j,npad
          b1(i)=0.0
          b2(i)=0.0
        end do
        call fft2(b1,b2,npad)
        sum=0.
        j = 0
        do i=1,npad/2
          j = j+1
          ta(j,iwin)=b1(i)/anrm
          tai(j,iwin)=b2(i)/anrm
          b(i) = b1(i)
          b(i+1) = b2(i)
          sumi=b1(i)*b1(i)+b2(i)*b2(i)
          sum=sum+sumi
        end do
      end do
c
      if (ispec.eq.0) then
c
c       calculate "high-resolution" spectrum
c
        call hires(ta,tai,el,nwin,nf-1,amu)
c
c       determine psd as squared amplitude spectrum
c       dissallow zero values
c
        do i=1,nf-1
          amu(i,1)=abs(amu(i,1))**2
        end do
      endif
c
      if (ispec.eq.1) then
c
c       calculate adaptively weighted spectrum
c
        avar=0.d0
        do i=1,nf-1
          avar=avar+demean(i)*demean(i)
        end do
c
c       avar is a factor that scales the bias factor in the adaptive
c       weighting scheme.
c
        if (inorm.eq.1) then ! amp spect per unit time
          avar=avar/(npts*npts)
        elseif(inorm.eq.2) then  ! absolute amp spect
          avar=avar*dt*dt
        endif
c
        call adwait(ta,tai,dcf,el,nwin,nf-1,amu,amu(1,2),avar)
c
c       normalized psd, disallow negative values
c
        do i=1,nf-1
          amu(i,1)=abs(amu(i,1))
        end do
      endif
c
c     initialize raw, reshaped, and harmonic spectra as the psd estimated
c     above, and the f-test as null
c
      do i=1,nf-1
         iharm(i)=0
         specraw0(i)=amu(i,1)
         specresh0(i)=amu(i,1)
         harmonic0(i)=amu(i,1)
         idone(i)=0
      end do
c
      if (iftest.eq.1) then
c
c
c       perform f-test for phase coherence if indicated
c
        call regre(ta,tai,nf-1,nwin,amu)
c
c
        do i=1,nf-1
           ftest(i)=amu(i,3)
        end do
      endif
c
c     determine the average (white) power level of the raw spectrum
c
      white0 = 0.0
      do i=1,nf-1
         white0 = white0+specraw0(i)
      end do
      white0 = white0/float(nf-1)
c
c     set default noise parameters as raw noise parameters
c
      rho = rho0
      white = white0
c
      if (ismooth.eq.1) then
c
c      determine median smoothed spectrum and associated
c     "robust" average (white) power level
c
       white = 0.0
       nsmooth = fsmooth/ddf
       do j=1,nf-1
         if1 = j-nsmooth/2
         if2 = j+nsmooth/2
         if (if1.lt.1) then
           if1=1
c           if2=nsmooth+1
         endif
         if (if2.gt.nf-1) then
           if2=nf-1
c           if1=nf-nsmooth-1
         endif
         nblk = if2-if1+1
         do i=1,nblk
            val(i)=specraw0(if1+i-1)
         end do
c
c        sort spectrum in this block
c
         kmid = (nblk+1)/2
         do kk1=1,nblk
            do kk2=kk1+1,nblk
               if (val(kk2).lt.val(kk1)) then
                  adum = val(kk2)
                  val(kk2)=val(kk1)
                  val(kk1)=adum
               endif
            end do
         end do
         specmed(j)=val(kmid)
         specmed0(j)=specmed(j)
         white = white + specmed(j)
       end do
       white = white/float(nf-1)
c
       if (inoise.ne.2) then
c
c         Unless user has selected the "local white" assumption,
c         we attempt to fit a parametric (white or red) noise
c         background robustly
c
c         determine best fit spectrum of the form
c
c         rednoise = white*(1.0-rho**2)/
c                     (1.0-2.0*rho*cos(arg)+rho**2)
c
c         to the median-smoother.
c
c         note that the "white noise" is  a trivial case
c
c        "white" is the robust white noise power level
c         estimated above from the median smoothed
c         spectrum
c
c         do a global search of the interval [0,1)
c         to find the optimal rho as determined by minimum MSE
c
          amin = big
          do rho1=0.0,0.999,0.001
           amiss = MSE(rho1)
           if (amiss.lt.amin) then
            rho=rho1
            amin = amiss
           endif
          end do
       endif
      endif
c
c     determine raw and robust estimates of the decorrelation
c     timescale of the noise
c
      tau = 1e+16
      tau0 = 1e+16
      if (rho.gt.0.0) tau = -dt/log(rho)
      if (rho0.gt.0.0) tau0 = -dt/log(rho0)
c
c     determine the noise background
c
      do i=1,nf-1
         ff = fmin+(i-1)*ddf
         freq_norm = ff/fnynew
         arg = freq_norm*pi
         rednoise0(i) = white0*(1.0-rho0**2)/
     $      (1.0-2.0*rho0*cos(arg)+rho0**2)
         rednoise(i) = white*(1.0-rho**2)/
     $      (1.0-2.0*rho*cos(arg)+rho**2)
c
c        if the "harmonic signal" option is selected,
c        the estimated null "base" spectrum is zero,
c        otherwise it is the noise background determined above
c
         if (isignal.eq.2) then
            base(i)=0.0
         else
           if (inoise.lt.2) then
             if (ismooth.eq.1) then
                base(i) = rednoise(i)
             else
                base(i) = rednoise0(i)
             endif
           else
             base(i) = specmed(i)*conf(1)
           endif
           specback0(i)=base(i)
         endif
      end do
c
c     now perform reshaping procedure if indicated
c
      if (iresh.eq.1) then
c
c       do a poor mans reshaping -detect harmonic peaks and then
c       interpolate the continuous spectrum across the effected
c       bandwidth  -only reshape if harmonic peak is greater than
c       the significance level in terms of overall power that was
c       indicated for the F-test harmonic detection procedure
c
c       note that for harmonic signal assumption, base=0 and all
c       signficant f-test peaks will be reshaped
c
c       note: reshaping is only done at frequencies outside
c       the secular band
c
        do i=nbnd,nf-1
           if ((ftest(i).gt.thresh).
     $         and.(specraw0(i).gt.afact*base(i))) iharm(i)=1
        end do
c
        do i=nbnd,nf-1
c
c           determine frequency points at boarder of the Rayleigh bandwidth
c           and the full spectral bandwidth
c
c           if frequency falls in a band within reshaping was already
c           performed, we skip
c
            if (idone(i).ne.1) then
c
            ipre0 = i-nbnd0/2-1
            iaft0 = i+nbnd0/2+1
            ipre = i-nbnd/2
            iaft = i+nbnd/2
            if (ipre.lt.nbnd+1) ipre=nbnd+1
            if (iaft.gt.nf-2) iaft=nf-2
            if (ipre0.lt.nbnd) ipre0=nbnd
            if (iaft0.gt.nf-1) iaft0=nf-1
c
c           reshaped spectrum is estimating by assuming that the
c           spectrum is continuous across the Rayleigh spectral
c           bandwidth within which a periodic signal is detected.
c           this gap is linearly interpolated across that bandwidth
c
c           "harmonic0" is the sum of the reshaped "continuous background"
c           and the detected spectral line. The line component is estimated
c           by assuming that the total power within the reshaped region
c           is represented by the sum of the continuous background and
c           a narrow peak with the Rayleigh resolution
c
            sumpeak = 0.0
            if (iharm(i).eq.1) then
             do j=ipre,iaft
              specresh0(j)=0.5*(specraw0(ipre-1)+specraw0(iaft+1))
              harmonic0(j)=specresh0(j)
             end do
             do j=ipre,iaft
               idone(j)=1
               sumpeak = sumpeak+specraw0(j)
             end do
             do j=ipre0,iaft0
                harmonic0(j)=sumpeak/(iaft0-ipre0+1)
             end do
            endif

            endif

        end do
      endif
c
c     for periodic signal assumption, the estimated noise background
c     is simply the reshaped (ie, estimated continuous) spectrum
c
      if (isignal.eq.2) then
         do i=1,nf-1
            specback0(i)=specresh0(i)
         end do
      endif
c
      whiteraw0=white0
      whiterob0=white
      rhoraw0=rho0
      rhorob0=rho
      tauraw0 = tau0
      taurob0 = tau
c
8888  continue
c
      return
      end
c
      real function MSE(rho)
c
      parameter (nlim=32768)
      COMMON /BLK1/white,specmed(nlim),fnynew,ddf,f0
      COMMON /BLK2/ilog,nfreq
      real rho
      real pie,dff,freq_norm,ff,arg,small,rednoise,val1,val2
      pie = 3.1415926
      small = 1e-12
      dff = 0.0
      do j=1,nfreq
         ff = f0+(j-1)*ddf
         freq_norm = ff/fnynew
         arg = freq_norm*pie
         rednoise = white*(1.0-rho**2)/
     $      (1.0-2.0*rho*cos(arg)+rho**2)
         if (ilog.eq.0) then
           dff = dff + (specmed(j)-rednoise)**2
         else
           val1 = abs(specmed(j))
           val2 = abs(rednoise)
           if (val1.lt.small) val1=small
           if (val2.lt.small) val2=small
           dff = dff +
     $     (log(val1)-log(val2))**2
         endif
      end do
      MSE = dff
      return
      end
c
      subroutine hires(ta,tai,el,nwin,nf,ares)
      real*8 el
      dimension ta(16400,1),tai(16400,1),el(1),ares(1)
      do j=1,nf
        ares(j)=0.
      end do
      do i=1,nwin
        a=1./(el(i)*nwin)
        do j=1,nf
          ares(j)=ares(j)+a*(ta(j,i)*ta(j,i)+tai(j,i)*tai(j,i))
        end do
      end do
      do j=1,nf
        ares(j)=sqrt(ares(j))
      end do
      return
      end
c
c
      subroutine adwait(ta,tai,dcf,el,nwin,nf,ares,degf,avar)
c
c  this version uses Thomson's algorithm for calculating
c  the adaptive spectrum estimate
c
      real*8 avar,spw,as,das,tol,el,a1,bias,scale,ax,fn,fx
      dimension ta(16400,1),tai(16400,1),el(1),ares(1),degf(1)
      dimension spw(10),bias(10),dcf(16400,1)
c
c  set tolerance for iterative scheme exit
c
      tol=3.d-4
      jitter=0
      scale=avar
c
c  we scale the bias by the total variance of the frequency transform
c  from zero freq to the nyquist
c  in this application we scale the eigenspectra by the bias in order to avoid
c  possible floating point overflow
c
      do 200 i=1,nwin
  200 bias(i)=(1.d0-el(i))
      do 100 j=1,nf
        do 150 i=1,nwin
  150   spw(i)=(ta(j,i)*ta(j,i)+tai(j,i)*tai(j,i))/scale
c
c  first guess is the average of the two lowest-order eigenspectral estimates
c
        as=(spw(1)+spw(2))/2.d0
        do 300 k=1,20
c
c  find coefficients
c
          fn=0.d0
          fx=0.d0
          do 350 i=1,nwin
            a1=dsqrt(el(i))*as/(el(i)*as+bias(i))
            a1=a1*a1
            fn=fn+a1*spw(i)
            fx=fx+a1
  350     continue
          ax=fn/fx
          das=dabs(ax-as)
          if(das/as.lt.tol) go to 400
  300   as=ax
c
c  flag if iteration does not converge
c
      jitter=jitter+1
  400 continue
      ares(j)=as*scale
c
c  calculate degrees of freedom
c
      df=0.
      do 450 i=1,nwin
      dcf(j,i)=dsqrt(el(i))*as/(el(i)*as+bias(i))
  450 df=df+dcf(j,i)*dcf(j,i)
c
c  we normalize degrees of freedom by the weight of the first eigenspectrum
c  this way we never have fewer than two degrees of freedom
c
  100 degf(j)=df*2./(dcf(j,1)*dcf(j,1))
      return
      end
c
      subroutine regre(sr,si,nf,nwin,amu)
c
      real b
      real *8 junk
      dimension sr(16400,1),si(16400,1),amu(16400,1)
      common/tapsum/b(10),junk(10)
c
c     "b" is the DFT of Slepian eigentapers at zero frequency
c     "sr" and "si" are the eigenspectra
c     "amu" contains line frequency estimates and f-test parameter
c
      sum=0.
      do i=1,nwin
        sum=sum+b(i)*b(i)
      end do
      do i=1,nf
        amu(i,1)=0.
        amu(i,2)=0.
        do j=1,nwin
          amu(i,1)=amu(i,1)+sr(i,j)*b(j)
          amu(i,2)=amu(i,2)+si(i,j)*b(j)
        end do
        amu(i,1)=amu(i,1)/sum
        amu(i,2)=amu(i,2)/sum
        sum2=0.
        do j=1,nwin
          sumr=sr(i,j)-amu(i,1)*b(j)
          sumi=si(i,j)-amu(i,2)*b(j)
          sum2=sum2+sumr*sumr+sumi*sumi
        end do
        amu(i,3)=(nwin-1)*(amu(i,2)**2+amu(i,1)**2)*sum/sum2
      end do
      return
      end
c
c
      subroutine taper2(n,nwin,el)
c
c     generate slepian tapers
c     "ta" is a real*4 array
c
c     written by J. Park
c
      real*8 el,a,z,pi,ww,cs,ai,an,eps,rlu,rlb
      real*8 dfac,drat,gamma,bh,ell
      real dump
      common/npiprol/anpi
      common/tapsum/tapsum(10),ell(10)
      common/misc/npad
      common/work/ip(32800)
      common/taperz/z(65536),dump(393728)
      common/stap2/ta(32800,16)
      dimension a(32800,8),el(10)
      data pi/3.14159265358979d0/
      equivalence (a(1,1),ta(1,1))
      an=dfloat(n)
      ww=dble(anpi)/an
      cs=dcos(2.d0*pi*ww)
c
c     note:
c
c initialize matrix for eispack subroutine
c this matrix is not the bandwidth retention factor matrix A (for minimizing
c spectral leakage) described in various multitaper papers
c --e.g. Thomson (1982) [proc ieee], Park et al (1987) [jgr].
c the bandwidth matrix x returns a cluster of eigenvectors
c (which are the slepian tapers) with eigenvalues very close to unity.
c since the spacing of eigenvalues can be comparable to machine precision,
c numerical instability can occur for time-bandwidth product n .ge. 6
c (that is, 6pi-prolate tapers)
c also, the bandwidth retention matrix A is toeplitz, but full, and it is not
c feasible to calculate tapers explicitly for long time series (in an earlier
c code i interpolated an m-point taper from the eigenvectors of a 128x128
c bandwidth retention matrix). the following is cribbed from a 1978 paper
c by Slepian (Prolate spheroidal wave functions, Fourier analysis
c and uncertainty, Bell System Tech Journal, v57, pp1371-1430, 1978)
c which gives a three-term recursion that is satisfied by the
c slepian tapers. the three-term recursion can be manipulated to show that the
c slepian tapers are the eigenvectors of a tridiagonal matrix a, with
c well-spaced eigenvalues that are (practically speaking) unrelated to the
c eigenvalues of x. using this matrix a we obtain both numerical stability
c and speed (tridiagonal matrices can be decomposed rapidly enough to
c calculate m-point tapers on the fly).
c
      do i=0,n-1
        ai=dfloat(i)
        a(i+1,1)=-cs*((an-1.d0)/2.d0-ai)**2
        a(i+1,2)=-ai*(an-ai)/2.d0
c
c       next statement is eispack routine tridib, see its documentation
c
        a(i+1,3)=a(i+1,2)**2
      end do
      eps=1.e-13
      m11=1
      call tridib(n,eps,a(1,1),a(1,2),a(1,3),rlb,rlu,m11,nwin,el,ip,
     x       ierr,a(1,4),a(1,5))
      call tinvit(n,n,a(1,1),a(1,2),a(1,3),nwin,el,ip,z,ierr,
     x            a(1,4),a(1,5),a(1,6),a(1,7),a(1,8))
c
c  note:
c
c  we calculate the eigenvalues of the dirichlet-kernel problem
c  i.e. the bandwidth retention factors
c  from slepian 1978 asymptotic formula, gotten from thomson 1982 eq 2.5
c  supplemented by the asymptotic formula for k near 2n from Slepian (1978)
c  more precise values of these parameters, perhaps useful in adaptive
c  spectral estimation, can be calculated explicitly using the
c  rayleigh-quotient formulas in Thomson (1982) and Park et al (1987)
c
      dfac=an*pi*ww
      drat=8.d0*dfac
      dfac=4.d0*dsqrt(pi*dfac)*dexp(-2.d0*dfac)
      do k=1,nwin
        el(k)=1.d0-dfac
        dfac=dfac*drat/k  ! is this correct formula? yes,but fails as k -> 2n
      end do
      gamma=dlog(8.d0*an*dsin(2.d0*pi*ww))+0.5772156649d0
      do k=1,nwin
        bh=-2.d0*pi*(an*ww-dfloat(k-1)/2.d0-.25d0)/gamma
        ell(k)=1.d0/(1.d0+dexp(pi*bh))
      end do
      do i=1,nwin
        el(i)=dmax1(ell(i),el(i))
      end do
c
c  normalize the eigentapers to preserve power for a white process
c  i.e. they have rms value unity
c  "tapsum" is the average of the eigentaper, should be near zero for
c  antisymmetric tapers
c
      do k=1,nwin
        kk=(k-1)*n
        tapsum(k)=0.
        tapsq=0.
        do i=1,n
          aa=z(kk+i)
          ta(i,k)=aa
          tapsum(k)=tapsum(k)+aa
          tapsq=tapsq+aa*aa
        end do
        aa=sqrt(tapsq/n)
        tapsum(k)=tapsum(k)/aa
        do i=1,n
          ta(i,k)=ta(i,k)/aa
        end do
      end do
c
c  the real FFT will preserve amplitudes with zeropadding
c  for example, a(i)=1.,i=1,100 will transform at zero freq
c  to b(f=0)=100 no matter how much zero padding is done
c  therefore we need not doctor the taper normalization,
c  but wait until the fft to force the desired units
c
      return
      end
c
c
c
c Subroutine findpeaks
c ===========================
      subroutine findpeaks(flow,fhigh,isignal,
     $   nf0,df,nbnd,specraw0,specmed0,specresh0,harmonic0,
     $   specback0,ftest,
     $   conf,fconf,
     $   fsignal,confsignal,nsignals)
c
      parameter (maxsignal=40,nlarge=300)
      parameter (nlim=32768)
      real fsignal(maxsignal),confsignal(maxsignal)
      real fsig(nlarge),confsig(nlarge),sigrat(nlarge)
      integer null(nlim)
      integer nbnd,isignal
      real fconf(6),conf(4),sig(3)
      real ratio,ratmax,ratio0,dummy
      real specraw0(nlim),specmed0(nlim),specresh0(nlim),
     $   harmonic0(nlim),specback0(nlim),ftest(nlim)
c      real whiteraw0,whiterob0,rhoraw0,rhorob0,tauraw0,taurob0
c
c     determine central frequency of all signals above the 90%
c     level (or the maxsignal most significant signals)
c
      nsignals = 0
      do i=1,nf0
         null(i)=0
      end do
      sig(1)=99.0
      sig(2)=95.0
      sig(3)=90.0
c
c     loop search in order of decreasing signficance:99,95,90% conf
c
      do k=1,3
       if (isignal.lt.2) then
        thresh = conf(5-k)
       else
        thresh = fconf(5-k)
       endif
       do i=1,nf0-1
         ff = flow+float(i-1)*df
c
         if (isignal.lt.2) then
            ratio0 = specraw0(i)/specback0(i)
         else
            ratio0 = ftest(i)
         endif
         if ((ratio0.gt.thresh).and.(null(i).eq.0)) then
c
           nsignals=nsignals+1
c
c          if max signal capacity exceeded, quit
c
           if (nsignals.gt.nlarge) then
              nsignals=nsignals-1
              goto 888
           endif
c
c          peak identified -- search for "center" of peak
c          within the spectral bandwidth
c
           ratmax = ratio0
           imax = i
           fmax = ff
           do j=max(1,i-nbnd),min(i+nbnd,nf0-1)
              f0 = flow+float(j-1)*df
              if (isignal.lt.2) then
                ratio = specraw0(j)/specback0(j)
              else
                ratio = ftest(j)
              endif
              if (ratio.gt.ratmax) then
                imax = j
                ratmax=ratio
                fmax = f0
              endif
           end do
c
c          center (significance-wise) of peak identified
c
c          if the center of the peak now falls within one
c          bandwidth of a previously recognized peak, we
c          must throw it out
c
           if (null(imax).eq.1) then
              nsignals = nsignals-1
              goto 444
           endif
c
c          otherwise, we happily add to the list
c
           fsig(nsignals)=fmax
           confsig(nsignals)=sig(k)
           sigrat(nsignals)=ratmax
c
c          now block out the bandwidth from future peak identification
c
           do j=max(1,imax-nbnd),min(imax+nbnd,nf0-1)
              null(j)=1
           end do
c
         endif
444      continue
       end do
      end do
888   continue
c
c     sort signals in terms of significance
c
      do i=1,nsignals
         do j=i+1,nsignals
            if (sigrat(j).gt.sigrat(i)) then
              dummy = sigrat(j)
              sigrat(j)=sigrat(i)
              sigrat(i) = dummy
              dummy = fsig(j)
              fsig(j)=fsig(i)
              fsig(i) = dummy
              dummy =  confsig(j)
              confsig(j)=confsig(i)
              confsig(i) = dummy
            endif
         end do
      end do
c
c     now store only the first "maxsignal" of these
c
      nsignals=min(maxsignal,nsignals)
c
c     sort the remaining peaks by increasing frequency
c
      do i=1,nsignals
         do j=i+1,nsignals
            if (fsig(j).lt.fsig(i)) then
              dummy = fsig(i)
              fsig(i)=fsig(j)
              fsig(j)=dummy
              dummy = confsig(i)
              confsig(i)=confsig(j)
              confsig(j)=dummy
            endif
         end do
      end do
      do i=1,nsignals
        fsignal(i)=fsig(i)
        confsignal(i)=confsig(i)
      end do
      return
      end
c
c
c
c Subroutine display
c =======================
      subroutine display(nf0,f1,df,
     $        spec_raw,spec_resh,spec_smoo,spec_conf,
     $     specraw0,specmed0,specresh0,harmonic0,
     $     specback0,ftest,conf,fconf,
     $     whiteraw0,whiterob0,rhoraw0,rhorob0,
     $     tauraw0,taurob0,
     $     iplotresh,iplotftest,iplotsmoo,
     $     iplotraw,iplotconf)
c
      parameter (nlim=32768)
      integer iplotraw,iplotresh,iplotftest,iplotsmoo,
     $    iplotconf
      real fconf(6),conf(4)
        real spec_raw(nlim,3),spec_resh(nlim,3)
        real spec_smoo(nlim,2),spec_conf(nlim,5)
      real specraw0(nlim),specmed0(nlim),specresh0(nlim),
     $   harmonic0(nlim),specback0(nlim),ftest(nlim)
      real whiteraw0,whiterob0,rhoraw0,rhorob0,tauraw0,taurob0
      real f1,df
      integer nf0
c
c     Writting output arrays
c
        do i=1,nf0-1
            ff = f1+float(i-1)*df
         if (iplotraw.eq.1) then
                spec_raw(i,1) = ff
                spec_raw(i,2) = specraw0(i)
                spec_raw(i,3) = ftest(i)
            endif
         if (iplotresh.eq.1) then
                spec_resh(i,1) = ff
                spec_resh(i,2) = specresh0(i)
                spec_resh(i,3) = harmonic0(i)
            endif
         if (iplotsmoo.eq.1) then
                spec_smoo(i,1) = ff
                spec_smoo(i,2) = specmed0(i)
            endif
         if (iplotconf.eq.1) then
                spec_conf(i,1) = ff
                spec_conf(i,2) = specback0(i)*conf(1)
            spec_conf(i,3) = specback0(i)*conf(2)
                spec_conf(i,4) = specback0(i)*conf(3)
                spec_conf(i,5) = specback0(i)*conf(4)
            endif
      end do
c
c       Writting mtm_spec.inf file
c
c      open (unit=27,file='mtm_spec_inf',status='unknown')
c      write (27,*) 'raw determination:'
c      write (27,*) '    white noise variance: ',whiteraw0
c      write (27,*) '    rho: ',rhoraw0
c      write (27,*) '    tau: ',tauraw0
c      write (27,*) 'robust determination:'
c      write (27,*) '    white noise variance: ',whiterob0
c      write (27,*) '    rho: ',rhorob0
c      write (27,*) '    tau: ',taurob0
c      close (unit=27)
c
8888  continue
c
652   format (f7.4,e32.8)
653   format (f7.4,2e30.8)
654   format (f7.4,1x,4e17.3)
655   format (f7.4,1x,5e14.2)
c
      return
      end
c
c
c
c Subroutine recon
c =====================
      subroutine recon(a0,nwin0,npi,dt0,nscan,bndwdth,ff0,icon,rcon)
c
c     perform time-reconstruction of signal based on inversion
c     of multitaper envelope.
c
c     modified by M.E. Mann 12/96 based on original code of J. Park
c
c  program is limited to 16 windows
c  program is limited to 4100 frequency-ftest estimation pairs
c
      parameter (maxlen=50000,maxsignal=20)
      real a0(maxlen),dt0
      integer nwin0,nscan,npi
      real*8 ell,env0,cs,sn,ex
      real bndwdth,f0,ff0
      character*4 chead(158)
      complex zed
      common/npiprol/anpi
      common/header/iah(158)
      common/data/a(maxlen)
      common/nnaamm/iwflag
      common/work/b(8194)
      common/pltopt/npad,nst,mpts,nmin,nmax,t(3),fmin,fmax,nf,
     x                                ddt,nwin,inc,tt(3)
      common/taperz/ta(4100,16),tai(4100,16)
      common/staple/tad(4000)
      common/stap2/tas(8192,16),ex(500)
      common/stap4/cs(8200),sn(8200)
      common/envel1/ell(20)
      dimension dum(35),ahead(158),env0(2,1000)
      dimension reim(2),time(2)
      real *8 b1(8194),b2(8194)
      real env(3,500,2),aa(8200,3),acent(8200),rcon(maxlen)
      equivalence (zed,reim),(iy,dum),(iah,ahead),(iah,chead),
     x          (tad,env0)
      iwflag=0
      pi=3.14159265358979
c
      do i=1,maxlen
         a(i)=a0(i)
      end do
      time(1)=0.
      dt = dt0
      time(2)=dt/1.
      fzero=0.
      fny=0.5/dt
      npts=nscan
      nwin = nwin0
      anpi = float(npi)
c
c     set frequency range
c
      f1 = 0.0
      f2 = fny
c
      npad=2
      do while (npad.lt.nscan)
        npad=npad*2
      end do
      npad=8192
      ddf=1./(npad*dt)
      eps=.1*ddf
      nn1=(f1+eps)/ddf
      fmin=nn1*ddf
      nn2=(f2+eps)/ddf
      fmax=nn2*ddf
      nf=nn2-nn1+1
      nmin=nn1*2+1
      nmax=nn2*2+2
      t(1)=fmin
      t(2)=ddf
      t(3)=ddf
      ftfrq=ddf
  102 format(i6,g12.4,i6)
      n1=1
      n2=nscan
 2660 call taper1(npts,nwin,ell)
c
c     demean series
c
      demn=0.
      do i=n1,n2
         demn=demn+a(i)
      end do
      demn=demn/(n2-n1+1)
      do i=n1,n2
         acent(i)=a(i)-demn
      end do
c
c  loop over # of windows
c  normalization anrm:  divide by npts as we wish spectrum per unit time
      anrm=npts
      icont=0
      do iwin=1,nwin
        do i=n1,n2
          b1(i-n1+1)=acent(i)*tas(i+1-n1,iwin)
          b2(i-n1+1)=0.0
        end do
        j=npts+1
        do i=j,npad
          b1(i)=0.
          b2(i)=0.
        end do
        call fft2(b1,b2,npad)
        j = 0
        do i=1,npad/2
          j = j+1
          b(i)=b1(i)
          b(i+1)=b2(i)
          ta(i,iwin)=b1(i)
          tai(i,iwin)=b2(i)
        end do
      end do
c
c     decimate the tapers
c
1660  inc=(npts-1)/500+1
      j=0
      do i=1,npts,inc
        j=j+1
        do k=1,nwin
          tas(j,k)=tas(i,k)
        end do
      end do
      mpts=j
      ddt=dt*inc
c
c     need to take into account bandwith retention factor
c     (leads to double counting of secular variance)
c
c     if specified carrier frequency is within one bandwidth
c     of f=0, then specify its carrier as zero to avoid
c     difficulty in estimating bandwidth retention factor
c     for 0<f<bndwdth
c
      f0=ff0
      if (f0.lt.bndwdth) f0=0.0
c
c     reconstruction
c
      anorm = 1.0
      if (f0.lt.eps) anorm=0.5
c
c
      jmax = 1
      if (icon.eq.0) jmax=3
      do i=1,jmax
        itype = icon-1
        if (icon.eq.0) itype=i-1
        call envel(f0,itype)
        do j=1,mpts
           env(i,j,1)=env0(1,j)
           env(i,j,2)=env0(2,j)
        end do
      end do
c
      call cossin(npts,cs,sn,1.d0,0.d0,dble(f0),dble(dt))
c
      finc=float(inc)
      do jj=1,jmax
        do i=1,mpts
           env0(1,i)=env(jj,i,1)
           env0(2,i)=env(jj,i,2)
        end do
c
c       linear interpolate the envelope from mpts (<500) to npts
c
        do i=1,mpts-1
          ii=(i-1)*inc
          do j=1,inc
            a1=(env0(1,i)*(inc+1-j)+env0(1,i+1)*(j-1))/finc
            a2=(env0(2,i)*(inc+1-j)+env0(2,i+1)*(j-1))/finc
            a(n1-1+nscan+ii+j)=a1*cs(ii+j)+a2*sn(ii+j)
            aa(n1-1+ii+j,jj)=anorm*a(n1-1+nscan+ii+j)
          end do
        end do
c
c       linear extrapolate envelope one point
c
        ii=(mpts-1)*inc
        e1=2.*env0(1,mpts)-env0(1,mpts-1)
        e2=2.*env0(2,mpts)-env0(2,mpts-1)
        do j=1,inc
          a1=(env0(1,mpts)*(inc+1-j)+e1*(j-1))/finc
          a2=(env0(2,mpts)*(inc+1-j)+e2*(j-1))/finc
          a(n1-1+nscan+ii+j)=a1*cs(ii+j)+a2*sn(ii+j)
          aa(n1-1+ii+j,jj)=anorm*a(n1-1+nscan+ii+j)
        end do
      end do
c
c     default constraint = pre-specified
c
      w1keep = 1.0
      w2keep = 0.0
      w3keep = 0.0
c
c     if indicated, determine weights on the 3 possible
c     reconstructions based on minimum-misfit criterion
c
      abest = 1.0e+36
      step = 0.02
      if (icon.eq.0) then
        do w1=0,1,step
          do w2=0,1,step
            w3 = 1.0-w1-w2
            if ((w3.ge.0.0).and.(w3.le.1.0)) then
              amiss = 0.0
              do i=1,nscan
                signal = w1*aa(i,1)+w2*aa(i,2)+w3*aa(i,3)
                amiss = amiss + (signal-acent(i))**2
              end do
              if (amiss.lt.abest) then
                 abest=amiss
                 w1keep = w1
                 w2keep = w2
                 w3keep = w3
              endif
            endif
           end do
        end do
        do i=1,nscan
          rcon(i)=w1keep*aa(i,1)+w2keep*aa(i,2)+w3keep*aa(i,3)
        end do
      else
        do i=1,nscan
           rcon(i)=aa(i,1)
        end do
        do i=nscan+1,maxlen
           rcon(i)=0.0
        end do
      endif
      return
      end
c
c
      subroutine taper1(n,nwin,el)
c
c     generate Slepian tapers
c     ta is a real*4 array
c
c     written by J. Park
c
      real*8 el,a,z,pi,ww,cs,ai,an,eps,rlu,rlb,ex
      real*8 dfac,drat,gamma,bh,ell
      common/nnaamm/iwflag
      common/npiprol/anpi
      common/tapsum/tapsum(20),ell(20)
      common/work/ip(8194)
      common/taperz/z(65600)  ! we use this common block for scratch space
      common/stap2/ta(8192,16),ex(500)
      dimension a(8192,8),el(10)
      data iwflag/0/,pi/3.14159265358979d0/
      equivalence (a(1,1),ta(1,1))
      an=dfloat(n)
      ww=dble(anpi)/an
      cs=dcos(2.d0*pi*ww)
      do i=0,n-1
        ai=dfloat(i)
        a(i+1,1)=-cs*((an-1.d0)/2.d0-ai)**2
        a(i+1,2)=-ai*(an-ai)/2.d0
        a(i+1,3)=a(i+1,2)**2        ! required by eispack routine
      end do
      eps=1.e-13
      m11=1
      call tridib(n,eps,a(1,1),a(1,2),a(1,3),rlb,rlu,m11,nwin,el,ip,
     x       ierr,a(1,4),a(1,5))
      call tinvit(n,n,a(1,1),a(1,2),a(1,3),nwin,el,ip,z,ierr,
     x            a(1,4),a(1,5),a(1,6),a(1,7),a(1,8))
c
c  we calculate the eigenvalues of the dirichlet-kernel problem
c  i.e. the bandwidth retention factors
c  from slepian 1978 asymptotic formula, gotten from thomson 1982 eq 2.5
c  supplemented by the asymptotic formula for k near 2n from slepian 1978
c
      dfac=an*pi*ww
      drat=8.d0*dfac
      dfac=4.d0*dsqrt(pi*dfac)*dexp(-2.d0*dfac)
      do k=1,nwin
        el(k)=1.d0-dfac
        dfac=dfac*drat/k  ! note that this doesn't hold as k -> 2n
      end do
      gamma=dlog(8.d0*an*dsin(2.d0*pi*ww))+0.5772156649d0
      do k=1,nwin
        bh=-2.d0*pi*(an*ww-dfloat(k-1)/2.d0-.25d0)/gamma
        ell(k)=1.d0/(1.d0+dexp(pi*bh))
      end do
      do i=1,nwin
        el(i)=dmax1(ell(i),el(i))
      end do
c
c   normalize the eigentapers to preserve power for a white process
c   i.e. they have rms value unity
c
      do k=1,nwin
        kk=(k-1)*n
        tapsum(k)=0.
        tapsq=0.
        do i=1,n
          aa=z(kk+i)
          ta(i,k)=aa
          tapsum(k)=tapsum(k)+aa
          tapsq=tapsq+aa*aa
        end do
        aa=sqrt(tapsq/n)
        tapsum(k)=tapsum(k)/aa
        do i=1,n
          ta(i,k)=ta(i,k)/aa
        end do
      end do
  101 format(80a)
c
c     The real DFT preserves amplitudes with zeropadding
c     therefore we need not doctor the taper normalization,
c      but wait until the fft to force the desired units
c
      iwflag=1
      return
      end
c
      subroutine envel(f0,inorm)
c
      real*8 cs,sn,env0,amp,ell,dex,ex
        real qc
      complex*16 g,d,env,dum,za,amp0,amp1,qrsave
      common/pltopt/npad,nst,npts,nmin,nmax,t(3),fmin,fmax,nf,
     x                                       dt,nwin,inc,tt(3)
      common/envel1/ell(20)
      common/taperz/ta(4100,16),tai(4100,16)
      common/staple/env(1000)
      common/stap2/tas(8192,16),ex(500)
      common/stap3/g(150000)
      common/stap4/cs(8200),sn(8200)
      common/stap5/d(100),za(100),qrsave(100),dum2(2000),dum(8200)
      dimension env0(2,1000),amp(2)
      equivalence (env0,env),(amp0,amp)
      tt(1)=0.
      tt(2)=dt
      pi=3.14159265358979
      rad=180./pi
      eps=1.e-4
      jsmoo = inorm
      iboost=0
      fmin=t(1)
      fmax=t(1)+(nf-1)*t(2)
      ff0 = f0
      qc=0.
      dex=dble(qc*pi/npts)
      call boost(iboost,npts,ex,dex)
      call fillup(f0)
c  dot product with decaying envelope for fixed element
c  then integrate the kernels twice
      call premult(iboost,jsmoo)
      call zqrdc(g,npts,npts,nwin,qrsave,ip,dumm,0)
      call zbacktr(nwin,npts,g,d,dum,2)
c  mult by za to get a scalar
      amp0=dcmplx(0.d0,0.d0)
      do k=1,nwin
        amp0=amp0+conjg(za(k))*dum(k)  ! note: conjugated
      end do
c  next we calculate double backtransform of za
c  and the triangular matrix r
      call zbacktr(nwin,npts,g,za,dum,2)
c  mult by za to get a scalar
      amp1=dcmplx(0.d0,0.d0)
      do k=1,nwin
        amp1=amp1+conjg(za(k))*dum(k)  ! note: conjugated
      end do
      amp0=amp0/amp1
c  subtract from eigenspectra
      sum1=0.
      do k=1,nwin
        sum1=sum1+(cdabs(d(k)))**2
      end do
      do k=1,nwin
        d(k)=d(k)-za(k)*amp0
      end do
      sum2=0.
      do k=1,nwin
        sum2=sum2+(cdabs(d(k)))**2
      end do
      rat=sum2/sum1
c  solve for m-tilde
c  we backtransform once, then multiply by q
      do i=1,npts
        env(i)=dcmplx(0.,0.)
      end do
      call zbacktr(nwin,npts,g,d,env,1)
      call zqrsl(g,npts,npts,nwin,qrsave,env,env,
     x dum,dum,dum,dum,10000,info)
      call postmult(jsmoo,npts,env,amp0,cs)
      npts2=npts*2
      return
      end
c
      subroutine zbacktr(n,npts,g,d,dum,ick)
      complex*16 g(npts,1),dum(1),d(1)
c
c  g is assumed to be the qr decomposition of a matrix
c  we only use the upper triangle of g
c  ick.eq.1 => only perform backtransform with lower triangle
c  the lower triangular backtransform:
c
      do k=1,n
        dum(k)=d(k)
      end do
      do i=1,n
        i1=i-1
        if(i.gt.1) then
          do j=1,i1
            dum(i)=dum(i)-conjg(g(j,i))*dum(j)
c            dum(i)=dum(i)-(g(j,i))*dum(j)
          end do
        endif
        dum(i)=dum(i)/conjg(g(i,i))
c        dum(i)=dum(i)/(g(i,i))
      end do
      if(ick.eq.1) return
c  the upper triangular backtransform:
      do i=n,1,-1
        ip1=i+1
        if(i.lt.n) then
          do j=ip1,n
            dum(i)=dum(i)-g(i,j)*dum(j)
c            dum(i)=dum(i)-conjg(g(i,j))*dum(j)
          end do
        endif
        dum(i)=dum(i)/g(i,i)
c        dum(i)=dum(i)/conjg(g(i,i))
      end do
      return
      end
c
      subroutine boost(iboost,npts,ex,dex)
      implicit real*8 (a-h,o-z)
      dimension ex(1)
      if(iboost.eq.1) then
c  boost the envelope by multiplying tapers by decreasing exponential
        dex=dexp(-dex)
        ex(1)=1.d0
        do i=2,npts
          ex(i)=ex(i-1)*dex
        end do
      else  ! dont boost the envelope
        do i=1,npts
          ex(i)=1.d0
        end do
      endif
      return
      end
c
      subroutine fillup(f0)
c
c  to fill data vector and kernel matrix
c
      real*8 cs,sn,ex
      complex*16 g,d,junk1,junk3
      real *4 junk2
      common/pltopt/npad,nst,npts,nmin,nmax,t(3),fmin,fmax,nf,
     x                                dt,nwin,inc,tt(3)
      common/taperz/ta(4100,16),tai(4100,16)
      common/stap2/tas(8192,16),ex(500)
      common/stap3/g(150000)
      common/stap4/cs(8200),sn(8200)
      common/stap5/d(100),junk1(200),junk2(2000),junk3(8200)
      n1=(f0-fmin)/t(2)+1
      f1=(fmin+(n1-1)*t(2))
      df1=f0-f1
      if(df1*2.0.gt.t(2)) then
        n1=n1+1
        f1=f1+t(2)
        df1=f1-f0
      else
        df1=-df1    ! since formula uses f1-f0
      endif
c  ************* careful! make certain that freq and dt are reciprocal units!
c  ************* if dt is sec, df1 is mhz, need to scale by 1000.
c  ************* same if dt is kyr and df1 is cyc/myr
      call cossin(npts,cs,sn,1.d0,0.d0,dble(df1),dble(dt))
c      call cossin(npts,cs,sn,1.d0,0.d0,dble(df1),dble(dt/1000.))
c  assemble data vector
c  multiply by 2.0/inc to obtain proper normalization
c  2.0 from <cos(t)**2>=1/2
c  and inc from decimation of tapers in the data kernels
      do k=1,nwin
        d(k)=dcmplx(ta(n1,k),tai(n1,k))*2.d0/inc
      end do
c  mult the tapers by sinusoid to get kernels
c  note that we use the complex conjugate
      do k=1,nwin
        ii=(k-1)*npts
        jj=(k-1)*npts
        j=0
        do i=1,npts
          j=j+1
          g(jj+j)=ex(i)*dcmplx(tas(i,k)*cs(j),-tas(i,k)*sn(j))
        end do
      end do
      return
      end
c
      subroutine premult(iboost,jsmoo)
      real*8 cs,sn,cs0
      complex*16 g,d,za,junk1,junk3
      real *4 junk2
      common/pltopt/npad,nst,npts,nmin,nmax,t(3),fmin,fmax,nf,
     x                                dt,nwin,inc,tt(3)
      common/stap3/g(150000)
      common/stap4/cs(8200),sn(8200)
      common/stap5/d(100),za(100),junk1(100),
     x    junk2(2000),junk3(8200)
      pi=3.14159265358979
      eps=1.e-4
c  dot product with decaying envelope for fixed element
      if(iboost.eq.0) then
        cs0=qc*pi/(npts-1)
        cs0=dexp(-cs0)
        cs(1)=1.d0
        do i=2,npts
          cs(i)=cs(i-1)*cs0
        end do
      else   ! if we have boosted the kernels, the fixed element is constant
        do i=1,npts
          cs(i)=1.d0
        end do
      endif
      do k=1,nwin
        za(k)=dcmplx(0.,0.)
        kk=(k-1)*npts
        do i=1,npts
          za(k)=za(k)+cs(i)*g(kk+i)
        end do
        za(k)=conjg(za(k))  ! must conjugate to adhere to conventions
      end do
c  then integrate the kernels twice
      if(jsmoo.gt.0) then
        do ksmoo=1,jsmoo
          do k=1,nwin
            ii=(k-1)*npts
            nptsm=npts-1
            do i=nptsm,1,-1
              g(ii+i)=g(ii+i)+g(ii+i+1)
            end do
c  little fudge factor to ensure penalty of first derivative
            g(ii+1)=g(ii+1)/eps
          end do
        end do
      endif
      return
      end
c
      subroutine postmult(jsmoo,npts,env,amp0,cs)
      real*8 cs(1)
      complex*16 env(1),amp0
      eps=1.e-4
      if(jsmoo.gt.0) then
c  post normalize to the envelope fluctuation
        do ksmoo=1,jsmoo
          env(1)=env(1)/eps
          do i=2,npts
            env(i)=env(i)+env(i-1)
          end do
        end do
      endif
c  add the decay envelope back in
      do i=1,npts
        env(i)=env(i)+amp0*cs(i)
      end do
      return
      end
c
      subroutine gfill(npts,nwin,ib,jb,df,ddt,g,ex,tas)
c  will fill ib,jb block of kernel matrix g with tapers
c  g has structure
c
c               |g11 g12|
c          g  = |g21 g22|
c
c  where each block is nptsxnwin, gij is ith carrier freq and jth fft freq
c  here ib is i jb is j
c
      implicit real*8 (a-h,o-z)
c      real*4 tas,tt(3),dum2
      real*4 tas,tt(3)
      real*8 junk1
      complex*16 g
      common/stap4/cs(500),sn(500),junk1(15400)
      dimension g(npts,2,nwin,2),ex(1),tas(8192,1)
c
      call cossin(npts,cs,sn,1.d0,0.d0,df,ddt)
c  mult the tapers by sinusoid to get kernels
c  note that we use the complex conjugate
      tt(1)=0.
      tt(2)=ddt
      do k=1,nwin
        do i=1,npts
          g(i,ib,k,jb)=ex(i)*dcmplx(tas(i,k)*cs(i),-tas(i,k)*sn(i))
        end do
      end do
      return
      end
c
c
c
c Subroutines miscsubs
c =======================
c
c     calculate cos function by recursion relations for efficiency
c
      SUBROUTINE COSSIN(N,C,S,C0,S0,F,DT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(1),S(1)
      DATA PI/3.14159265358979D0/
      C(1)=C0
      S(1)=S0
      CS=DCOS(2.D0*PI*F*DT)
      SN=DSIN(2.D0*PI*F*DT)
      DO 100 I=2,N
      C(I)=C(I-1)*CS-S(I-1)*SN
  100 S(I)=C(I-1)*SN+S(I-1)*CS
      RETURN
      END
c
      SUBROUTINE FFT2(ar,ai,N)
c
c     Jeff Park's fft
c
c     2 real input arrays rather than complex
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension ar(1),ai(1)
      mex=dlog(dble(float(n)))/.693147d0
      nv2=N/2
      nm1=N-1
      j=1
      do 7 i=1,nm1
      if(i .ge. j) go to 5
      tr=ar(j)
      ti=ai(j)
      ar(j)=ar(i)
      ai(j)=ai(i)
      ar(i)=tr
      ai(i)=ti
   5  k=nv2
   6  if(k .ge. j) go to 7
      j=j-k
      k=k/2
      go to 6
   7  j=j+k
      pi=3.14159265358979d0
      do 20 l=1,mex
      le=2**l
      le1=le/2
      ur=1.0
      ui=0.
      wr=dcos(pi/le1 )
      wi=dsin (pi/le1)
      do 20 j=1,le1
      do 10 i=j,N,le
      ip=i+le1
      tr=ar(ip)*ur - ai(ip)*ui
      ti=ai(ip)*ur + ar(ip)*ui
      ar(ip)=ar(i)-tr
      ai(ip)=ai(i) - ti
      ar(i)=ar(i)+tr
      ai(i)=ai(i)+ti
  10  continue
      utemp=ur
      ur=ur*wr - ui*wi
      ui=ui*wr + utemp*wi
  20  continue
      return
      end
c
c     subroutines TRIDIB, TINVIT from DEISPK
c
c
      SUBROUTINE TRIDIB(N,EPS1,D,E,E2,LB,UB,M11,M,W,IND,IERR,RV4,RV5)
      INTEGER I,J,K,L,M,N,P,Q,R,S,II,M1,M2,M11,M22,TAG,IERR,ISTURM
      DOUBLE PRECISION D(N),E(N),E2(N),W(M),RV4(N),RV5(N)
      DOUBLE PRECISION U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,TST1,TST2,EPSLON
      INTEGER IND(M)
      IERR = 0
      TAG = 0
      XU = D(1)
      X0 = D(1)
      U = 0.0D0
      DO 40 I = 1, N
         X1 = U
         U = 0.0D0
         IF (I .NE. N) U = DABS(E(I+1))
         XU = DMIN1(D(I)-(X1+U),XU)
         X0 = DMAX1(D(I)+(X1+U),X0)
         IF (I .EQ. 1) GO TO 20
         TST1 = DABS(D(I)) + DABS(D(I-1))
         TST2 = TST1 + DABS(E(I))
         IF (TST2 .GT. TST1) GO TO 40
   20    E2(I) = 0.0D0
   40 CONTINUE
      X1 = N
      X1 = X1 * EPSLON(DMAX1(DABS(XU),DABS(X0)))
      XU = XU - X1
      T1 = XU
      X0 = X0 + X1
      T2 = X0
      P = 1
      Q = N
      M1 = M11 - 1
      IF (M1 .EQ. 0) GO TO 75
      ISTURM = 1
   50 V = X1
      X1 = XU + (X0 - XU) * 0.5D0
      IF (X1 .EQ. V) GO TO 980
      GO TO 320
   60 IF (S - M1) 65, 73, 70
   65 XU = X1
      GO TO 50
   70 X0 = X1
      GO TO 50
   73 XU = X1
      T1 = X1
   75 M22 = M1 + M
      IF (M22 .EQ. N) GO TO 90
      X0 = T2
      ISTURM = 2
      GO TO 50
   80 IF (S - M22) 65, 85, 70
   85 T2 = X1
   90 Q = 0
      R = 0
  100 IF (R .EQ. M) GO TO 1001
      TAG = TAG + 1
      P = Q + 1
      XU = D(P)
      X0 = D(P)
      U = 0.0D0
      DO 120 Q = P, N
         X1 = U
         U = 0.0D0
         V = 0.0D0
         IF (Q .EQ. N) GO TO 110
         U = DABS(E(Q+1))
         V = E2(Q+1)
  110    XU = DMIN1(D(Q)-(X1+U),XU)
         X0 = DMAX1(D(Q)+(X1+U),X0)
         IF (V .EQ. 0.0D0) GO TO 140
  120 CONTINUE
  140 X1 = EPSLON(DMAX1(DABS(XU),DABS(X0)))
      IF (EPS1 .LE. 0.0D0) EPS1 = -X1
      IF (P .NE. Q) GO TO 180
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 900
  180 X1 = X1 * (Q - P + 1)
      LB = DMAX1(T1,XU-X1)
      UB = DMIN1(T2,X0+X1)
      X1 = LB
      ISTURM = 3
      GO TO 320
  200 M1 = S + 1
      X1 = UB
      ISTURM = 4
      GO TO 320
  220 M2 = S
      IF (M1 .GT. M2) GO TO 940
      X0 = UB
      ISTURM = 5
      DO 240 I = M1, M2
         RV5(I) = UB
         RV4(I) = LB
  240 CONTINUE
      K = M2
  250    XU = LB
         DO 260 II = M1, K
            I = M1 + K - II
            IF (XU .GE. RV4(I)) GO TO 260
            XU = RV4(I)
            GO TO 280
  260    CONTINUE
  280    IF (X0 .GT. RV5(K)) X0 = RV5(K)
  300    X1 = (XU + X0) * 0.5D0
         IF ((X0 - XU) .LE. DABS(EPS1)) GO TO 420
         TST1 = 2.0D0 * (DABS(XU) + DABS(X0))
         TST2 = TST1 + (X0 - XU)
         IF (TST2 .EQ. TST1) GO TO 420
  320    S = P - 1
         U = 1.0D0
         DO 340 I = P, Q
            IF (U .NE. 0.0D0) GO TO 325
            V = DABS(E(I)) / EPSLON(1.0D0)
            IF (E2(I) .EQ. 0.0D0) V = 0.0D0
            GO TO 330
  325       V = E2(I) / U
  330       U = D(I) - X1 - V
            IF (U .LT. 0.0D0) S = S + 1
  340    CONTINUE
         GO TO (60,80,200,220,360), ISTURM
  360    IF (S .GE. K) GO TO 400
         XU = X1
         IF (S .GE. M1) GO TO 380
         RV4(M1) = X1
         GO TO 300
  380    RV4(S+1) = X1
         IF (RV5(S) .GT. X1) RV5(S) = X1
         GO TO 300
  400    X0 = X1
         GO TO 300
  420    RV5(K) = X1
      K = K - 1
      IF (K .GE. M1) GO TO 250
  900 S = R
      R = R + M2 - M1 + 1
      J = 1
      K = M1
      DO 920 L = 1, R
         IF (J .GT. S) GO TO 910
         IF (K .GT. M2) GO TO 940
         IF (RV5(K) .GE. W(L)) GO TO 915
         DO 905 II = J, S
            I = L + S - II
            W(I+1) = W(I)
            IND(I+1) = IND(I)
  905    CONTINUE
  910    W(L) = RV5(K)
         IND(L) = TAG
         K = K + 1
         GO TO 920
  915    J = J + 1
  920 CONTINUE
  940 IF (Q .LT. N) GO TO 100
      GO TO 1001
  980 IERR = 3 * N + ISTURM
 1001 LB = T1
      UB = T2
      RETURN
      END
c
      SUBROUTINE TINVIT(NM,N,D,E,E2,M,W,IND,Z,
     X                  IERR,RV1,RV2,RV3,RV4,RV6)
      INTEGER I,J,M,N,P,Q,R,S,II,IP,JJ,NM,ITS,TAG,IERR,GROUP
      DOUBLE PRECISION D(N),E(N),E2(N),W(M),Z(NM,M),
     X       RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)
      DOUBLE PRECISION U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,NORM,ORDER,EPSLON,
     X       PYTHAG
      INTEGER IND(M)
      IERR = 0
      IF (M .EQ. 0) GO TO 1001
      TAG = 0
      ORDER = 1.0D0 - E2(1)
      Q = 0
  100 P = Q + 1
      DO 120 Q = P, N
         IF (Q .EQ. N) GO TO 140
         IF (E2(Q+1) .EQ. 0.0D0) GO TO 140
  120 CONTINUE
  140 TAG = TAG + 1
      S = 0
      DO 920 R = 1, M
         IF (IND(R) .NE. TAG) GO TO 920
         ITS = 1
         X1 = W(R)
         IF (S .NE. 0) GO TO 510
         XU = 1.0D0
         IF (P .NE. Q) GO TO 490
         RV6(P) = 1.0D0
         GO TO 870
  490    NORM = DABS(D(P))
         IP = P + 1
         DO 500 I = IP, Q
  500    NORM = DMAX1(NORM, DABS(D(I))+DABS(E(I)))
         EPS2 = 1.0D-3 * NORM
         EPS3 = EPSLON(NORM)
         UK = Q - P + 1
         EPS4 = UK * EPS3
         UK = EPS4 / DSQRT(UK)
         S = P
  505    GROUP = 0
         GO TO 520
  510    IF (DABS(X1-X0) .GE. EPS2) GO TO 505
         GROUP = GROUP + 1
         IF (ORDER * (X1 - X0) .LE. 0.0D0) X1 = X0 + ORDER * EPS3
  520    V = 0.0D0
         DO 580 I = P, Q
            RV6(I) = UK
            IF (I .EQ. P) GO TO 560
            IF (DABS(E(I)) .LT. DABS(U)) GO TO 540
            XU = U / E(I)
            RV4(I) = XU
            RV1(I-1) = E(I)
            RV2(I-1) = D(I) - X1
            RV3(I-1) = 0.0D0
            IF (I .NE. Q) RV3(I-1) = E(I+1)
            U = V - XU * RV2(I-1)
            V = -XU * RV3(I-1)
            GO TO 580
  540       XU = E(I) / U
            RV4(I) = XU
            RV1(I-1) = U
            RV2(I-1) = V
            RV3(I-1) = 0.0D0
  560       U = D(I) - X1 - XU * V
            IF (I .NE. Q) V = E(I+1)
  580    CONTINUE
         IF (U .EQ. 0.0D0) U = EPS3
         RV1(Q) = U
         RV2(Q) = 0.0D0
         RV3(Q) = 0.0D0
  600    DO 620 II = P, Q
            I = P + Q - II
            RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)
            V = U
            U = RV6(I)
  620    CONTINUE
         IF (GROUP .EQ. 0) GO TO 700
         J = R
         DO 680 JJ = 1, GROUP
  630       J = J - 1
            IF (IND(J) .NE. TAG) GO TO 630
            XU = 0.0D0
            DO 640 I = P, Q
  640       XU = XU + RV6(I) * Z(I,J)
            DO 660 I = P, Q
  660       RV6(I) = RV6(I) - XU * Z(I,J)
  680    CONTINUE
  700    NORM = 0.0D0
         DO 720 I = P, Q
  720    NORM = NORM + DABS(RV6(I))
         IF (NORM .GE. 1.0D0) GO TO 840
         IF (ITS .EQ. 5) GO TO 830
         IF (NORM .NE. 0.0D0) GO TO 740
         RV6(S) = EPS4
         S = S + 1
         IF (S .GT. Q) S = P
         GO TO 780
  740    XU = EPS4 / NORM
         DO 760 I = P, Q
  760    RV6(I) = RV6(I) * XU
  780    DO 820 I = IP, Q
            U = RV6(I)
            IF (RV1(I-1) .NE. E(I)) GO TO 800
            U = RV6(I-1)
            RV6(I-1) = RV6(I)
  800       RV6(I) = U - RV4(I) * RV6(I-1)
  820    CONTINUE
         ITS = ITS + 1
         GO TO 600
  830    IERR = -R
         XU = 0.0D0
         GO TO 870
  840    U = 0.0D0
         DO 860 I = P, Q
  860    U = PYTHAG(U,RV6(I))
         XU = 1.0D0 / U
  870    DO 880 I = 1, N
  880    Z(I,R) = 0.0D0
         DO 900 I = P, Q
  900    Z(I,R) = RV6(I) * XU
         X0 = X1
  920 CONTINUE
      IF (Q .LT. N) GO TO 100
 1001 RETURN
      END
      DOUBLE PRECISION FUNCTION EPSLON (X)
      DOUBLE PRECISION X
      DOUBLE PRECISION A,B,C,EPS
      A = 4.0D0/3.0D0
   10 B = A - 1.0D0
      C = B + B + B
      EPS = DABS(C-1.0D0)
      IF (EPS .EQ. 0.0D0) GO TO 10
      EPSLON = EPS*DABS(X)
      RETURN
      END
      DOUBLE PRECISION FUNCTION PYTHAG(A,B)
      DOUBLE PRECISION A,B
      DOUBLE PRECISION P,R,S,T,U
      P = DMAX1(DABS(A),DABS(B))
      IF (P .EQ. 0.0D0) GO TO 20
      R = (DMIN1(DABS(A),DABS(B))/P)**2
   10 CONTINUE
         T = 4.0D0 + R
         IF (T .EQ. 4.0D0) GO TO 20
         S = R/T
         U = 1.0D0 + 2.0D0*S
         P = U*P
         R = (S/U)**2 * R
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
c================================================================
c   routines from zlinpack that are used in envelope estimation
c================================================================
c
      SUBROUTINE ZQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,JOB)
      INTEGER LDX,N,P,JOB
      INTEGER JPVT(1)
      COMPLEX*16 X(LDX,1),QRAUX(1),WORK(1)
      INTEGER J,JP,L,LP1,LUP,MAXJ,PL,PU
      DOUBLE PRECISION MAXNRM,DZNRM2,TT
      COMPLEX*16 ZDOTC,NRMXL,T
      LOGICAL NEGJ,SWAPJ
      COMPLEX*16 CSIGN,ZDUM,ZDUM1,ZDUM2
      DOUBLE PRECISION CABS1
      DOUBLE PRECISION DREAL,DIMAG
      COMPLEX*16 ZDUMR,ZDUMI
      DREAL(ZDUMR) = ZDUMR
      DIMAG(ZDUMI) = (0.0D0,-1.0D0)*ZDUMI
      CSIGN(ZDUM1,ZDUM2) = CDABS(ZDUM1)*(ZDUM2/CDABS(ZDUM2))
      CABS1(ZDUM) = DABS(DREAL(ZDUM)) + DABS(DIMAG(ZDUM))
      PL = 1
      PU = 0
      IF (JOB .EQ. 0) GO TO 60
         DO 20 J = 1, P
            SWAPJ = JPVT(J) .GT. 0
            NEGJ = JPVT(J) .LT. 0
            JPVT(J) = J
            IF (NEGJ) JPVT(J) = -J
            IF (.NOT.SWAPJ) GO TO 10
               IF (J .NE. PL) CALL ZSWAP(N,X(1,PL),1,X(1,J),1)
               JPVT(J) = JPVT(PL)
               JPVT(PL) = J
               PL = PL + 1
   10       CONTINUE
   20    CONTINUE
         PU = P
         DO 50 JJ = 1, P
            J = P - JJ + 1
            IF (JPVT(J) .GE. 0) GO TO 40
               JPVT(J) = -JPVT(J)
               IF (J .EQ. PU) GO TO 30
                  CALL ZSWAP(N,X(1,PU),1,X(1,J),1)
                  JP = JPVT(PU)
                  JPVT(PU) = JPVT(J)
                  JPVT(J) = JP
   30          CONTINUE
               PU = PU - 1
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
      IF (PU .LT. PL) GO TO 80
      DO 70 J = PL, PU
         QRAUX(J) = DCMPLX(DZNRM2(N,X(1,J),1),0.0D0)
         WORK(J) = QRAUX(J)
   70 CONTINUE
   80 CONTINUE
      LUP = MIN0(N,P)
      DO 200 L = 1, LUP
         IF (L .LT. PL .OR. L .GE. PU) GO TO 120
            MAXNRM = 0.0D0
            MAXJ = L
            DO 100 J = L, PU
               IF (DREAL(QRAUX(J)) .LE. MAXNRM) GO TO 90
                  MAXNRM = DREAL(QRAUX(J))
                  MAXJ = J
   90          CONTINUE
  100       CONTINUE
            IF (MAXJ .EQ. L) GO TO 110
               CALL ZSWAP(N,X(1,L),1,X(1,MAXJ),1)
               QRAUX(MAXJ) = QRAUX(L)
               WORK(MAXJ) = WORK(L)
               JP = JPVT(MAXJ)
               JPVT(MAXJ) = JPVT(L)
               JPVT(L) = JP
  110       CONTINUE
  120    CONTINUE
         QRAUX(L) = (0.0D0,0.0D0)
         IF (L .EQ. N) GO TO 190
            NRMXL = DCMPLX(DZNRM2(N-L+1,X(L,L),1),0.0D0)
            IF (CABS1(NRMXL) .EQ. 0.0D0) GO TO 180
               IF (CABS1(X(L,L)) .NE. 0.0D0)
     *            NRMXL = CSIGN(NRMXL,X(L,L))
               CALL ZSCAL(N-L+1,(1.0D0,0.0D0)/NRMXL,X(L,L),1)
               X(L,L) = (1.0D0,0.0D0) + X(L,L)
               LP1 = L + 1
               IF (P .LT. LP1) GO TO 170
               DO 160 J = LP1, P
                  T = -ZDOTC(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
                  CALL ZAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
                  IF (J .LT. PL .OR. J .GT. PU) GO TO 150
                  IF (CABS1(QRAUX(J)) .EQ. 0.0D0) GO TO 150
                     TT = 1.0D0 - (CDABS(X(L,J))/DREAL(QRAUX(J)))**2
                     TT = DMAX1(TT,0.0D0)
                     T = DCMPLX(TT,0.0D0)
                     TT = 1.0D0
     *                    + 0.05D0*TT
     *                      *(DREAL(QRAUX(J))/DREAL(WORK(J)))**2
                     IF (TT .EQ. 1.0D0) GO TO 130
                        QRAUX(J) = QRAUX(J)*CDSQRT(T)
                     GO TO 140
  130                CONTINUE
                        QRAUX(J) = DCMPLX(DZNRM2(N-L,X(L+1,J),1),0.0D0)
                        WORK(J) = QRAUX(J)
  140                CONTINUE
  150             CONTINUE
  160          CONTINUE
  170          CONTINUE
               QRAUX(L) = X(L,L)
               X(L,L) = -NRMXL
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
      RETURN
      END
c
      SUBROUTINE ZQRSL(X,LDX,N,K,QRAUX,Y,QY,QTY,B,RSD,XB,JOB,INFO)
      INTEGER LDX,N,K,JOB,INFO
      COMPLEX*16 X(LDX,1),QRAUX(1),Y(1),QY(1),QTY(1),B(1),RSD(1),XB(1)
      INTEGER I,J,JJ,JU,KP1
      COMPLEX*16 ZDOTC,T,TEMP
      LOGICAL CB,CQY,CQTY,CR,CXB
      COMPLEX*16 ZDUM
      DOUBLE PRECISION CABS1
      DOUBLE PRECISION DREAL,DIMAG
      COMPLEX*16 ZDUMR,ZDUMI
      DREAL(ZDUMR) = ZDUMR
      DIMAG(ZDUMI) = (0.0D0,-1.0D0)*ZDUMI
      CABS1(ZDUM) = DABS(DREAL(ZDUM)) + DABS(DIMAG(ZDUM))
      INFO = 0
      CQY = JOB/10000 .NE. 0
      CQTY = MOD(JOB,10000) .NE. 0
      CB = MOD(JOB,1000)/100 .NE. 0
      CR = MOD(JOB,100)/10 .NE. 0
      CXB = MOD(JOB,10) .NE. 0
      JU = MIN0(K,N-1)
      IF (JU .NE. 0) GO TO 40
         IF (CQY) QY(1) = Y(1)
         IF (CQTY) QTY(1) = Y(1)
         IF (CXB) XB(1) = Y(1)
         IF (.NOT.CB) GO TO 30
            IF (CABS1(X(1,1)) .NE. 0.0D0) GO TO 10
               INFO = 1
            GO TO 20
   10       CONTINUE
               B(1) = Y(1)/X(1,1)
   20       CONTINUE
   30    CONTINUE
         IF (CR) RSD(1) = (0.0D0,0.0D0)
      GO TO 250
   40 CONTINUE
         IF (CQY) CALL ZCOPY(N,Y,1,QY,1)
         IF (CQTY) CALL ZCOPY(N,Y,1,QTY,1)
         IF (.NOT.CQY) GO TO 70
            DO 60 JJ = 1, JU
               J = JU - JJ + 1
               IF (CABS1(QRAUX(J)) .EQ. 0.0D0) GO TO 50
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  T = -ZDOTC(N-J+1,X(J,J),1,QY(J),1)/X(J,J)
                  CALL ZAXPY(N-J+1,T,X(J,J),1,QY(J),1)
                  X(J,J) = TEMP
   50          CONTINUE
   60       CONTINUE
   70    CONTINUE
         IF (.NOT.CQTY) GO TO 100
            DO 90 J = 1, JU
               IF (CABS1(QRAUX(J)) .EQ. 0.0D0) GO TO 80
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  T = -ZDOTC(N-J+1,X(J,J),1,QTY(J),1)/X(J,J)
                  CALL ZAXPY(N-J+1,T,X(J,J),1,QTY(J),1)
                  X(J,J) = TEMP
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
         IF (CB) CALL ZCOPY(K,QTY,1,B,1)
         KP1 = K + 1
         IF (CXB) CALL ZCOPY(K,QTY,1,XB,1)
         IF (CR .AND. K .LT. N) CALL ZCOPY(N-K,QTY(KP1),1,RSD(KP1),1)
         IF (.NOT.CXB .OR. KP1 .GT. N) GO TO 120
            DO 110 I = KP1, N
               XB(I) = (0.0D0,0.0D0)
  110       CONTINUE
  120    CONTINUE
         IF (.NOT.CR) GO TO 140
            DO 130 I = 1, K
               RSD(I) = (0.0D0,0.0D0)
  130       CONTINUE
  140    CONTINUE
         IF (.NOT.CB) GO TO 190
            DO 170 JJ = 1, K
               J = K - JJ + 1
               IF (CABS1(X(J,J)) .NE. 0.0D0) GO TO 150
                  INFO = J
                  GO TO 180
  150          CONTINUE
               B(J) = B(J)/X(J,J)
               IF (J .EQ. 1) GO TO 160
                  T = -B(J)
                  CALL ZAXPY(J-1,T,X(1,J),1,B,1)
  160          CONTINUE
  170       CONTINUE
  180       CONTINUE
  190    CONTINUE
         IF (.NOT.CR .AND. .NOT.CXB) GO TO 240
            DO 230 JJ = 1, JU
               J = JU - JJ + 1
               IF (CABS1(QRAUX(J)) .EQ. 0.0D0) GO TO 220
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  IF (.NOT.CR) GO TO 200
                     T = -ZDOTC(N-J+1,X(J,J),1,RSD(J),1)/X(J,J)
                     CALL ZAXPY(N-J+1,T,X(J,J),1,RSD(J),1)
  200             CONTINUE
                  IF (.NOT.CXB) GO TO 210
                     T = -ZDOTC(N-J+1,X(J,J),1,XB(J),1)/X(J,J)
                     CALL ZAXPY(N-J+1,T,X(J,J),1,XB(J),1)
  210             CONTINUE
                  X(J,J) = TEMP
  220          CONTINUE
  230       CONTINUE
  240    CONTINUE
  250 CONTINUE
      RETURN
      END
c
      SUBROUTINE ZAXPY(N,ZA,ZX,INCX,ZY,INCY)
      COMPLEX*16 ZX(1),ZY(1),ZA
      DOUBLE PRECISION DCABS1
      IF(N.LE.0)RETURN
      IF (DCABS1(ZA) .EQ. 0.0D0) RETURN
      IF (INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        ZY(IY) = ZY(IY) + ZA*ZX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
   20 DO 30 I = 1,N
        ZY(I) = ZY(I) + ZA*ZX(I)
   30 CONTINUE
      RETURN
      END
c
      SUBROUTINE  ZCOPY(N,ZX,INCX,ZY,INCY)
      COMPLEX*16 ZX(1),ZY(1)
      INTEGER I,INCX,INCY,IX,IY,N
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        ZY(IY) = ZX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
   20 DO 30 I = 1,N
        ZY(I) = ZX(I)
   30 CONTINUE
      RETURN
      END
c
      SUBROUTINE  ZSCAL(N,ZA,ZX,INCX)
      COMPLEX*16 ZA,ZX(1)
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
      IX = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        ZX(IX) = ZA*ZX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
   20 DO 30 I = 1,N
        ZX(I) = ZA*ZX(I)
   30 CONTINUE
      RETURN
      END
c
      SUBROUTINE  ZSWAP (N,ZX,INCX,ZY,INCY)
      COMPLEX*16 ZX(1),ZY(1),ZTEMP
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        ZTEMP = ZX(IX)
        ZX(IX) = ZY(IY)
        ZY(IY) = ZTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
   20 DO 30 I = 1,N
        ZTEMP = ZX(I)
        ZX(I) = ZY(I)
        ZY(I) = ZTEMP
   30 CONTINUE
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION DZNRM2( N, ZX, INCX)
      LOGICAL IMAG, SCALE
      INTEGER          NEXT
      DOUBLE PRECISION CUTLO, CUTHI, HITEST, SUM, XMAX, ABSX, ZERO, ONE
      COMPLEX*16      ZX(1)
      DOUBLE PRECISION DREAL,DIMAG
      COMPLEX*16 ZDUMR,ZDUMI
      DREAL(ZDUMR) = ZDUMR
      DIMAG(ZDUMI) = (0.0D0,-1.0D0)*ZDUMI
      DATA         ZERO, ONE /0.0D0, 1.0D0/
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
      IF(N .GT. 0) GO TO 10
         DZNRM2  = ZERO
         GO TO 300
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
      DO 210 I=1,NN,INCX
         ABSX = DABS(DREAL(ZX(I)))
         IMAG = .FALSE.
         GO TO NEXT,(30, 50, 70, 90, 110)
   30 IF( ABSX .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      SCALE = .FALSE.
   50 IF( ABSX .EQ. ZERO) GO TO 200
      IF( ABSX .GT. CUTLO) GO TO 85
      ASSIGN 70 TO NEXT
      GO TO 105
  100 ASSIGN 110 TO NEXT
      SUM = (SUM / ABSX) / ABSX
  105 SCALE = .TRUE.
      XMAX = ABSX
      GO TO 115
   70 IF( ABSX .GT. CUTLO ) GO TO 75
  110 IF( ABSX .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / ABSX)**2
         XMAX = ABSX
         GO TO 200
  115 SUM = SUM + (ABSX/XMAX)**2
      GO TO 200
   75 SUM = (SUM * XMAX) * XMAX
   85 ASSIGN 90 TO NEXT
      SCALE = .FALSE.
      HITEST = CUTHI/FLOAT( N )
   90 IF(ABSX .GE. HITEST) GO TO 100
         SUM = SUM + ABSX**2
  200 CONTINUE
      IF(IMAG) GO TO 210
         ABSX = DABS(DIMAG(ZX(I)))
         IMAG = .TRUE.
      GO TO NEXT,(  50, 70, 90, 110 )
  210 CONTINUE
      DZNRM2 = DSQRT(SUM)
      IF(SCALE) DZNRM2 = DZNRM2 * XMAX
  300 CONTINUE
      RETURN
      END
c
c      COMPLEX FUNCTION ZDOTC*16(N,ZX,INCX,ZY,INCY)
      COMPLEX*16 FUNCTION ZDOTC(N,ZX,INCX,ZY,INCY)
      COMPLEX*16 ZX(1),ZY(1),ZTEMP
      ZTEMP = (0.0D0,0.0D0)
      ZDOTC = (0.0D0,0.0D0)
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        ZTEMP = ZTEMP + DCONJG(ZX(IX))*ZY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      ZDOTC = ZTEMP
      RETURN
   20 DO 30 I = 1,N
        ZTEMP = ZTEMP + DCONJG(ZX(I))*ZY(I)
   30 CONTINUE
      ZDOTC = ZTEMP
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION DCABS1(Z)
      COMPLEX*16 Z,ZZ
      DOUBLE PRECISION T(2)
      EQUIVALENCE (ZZ,T(1))
      ZZ = Z
      DCABS1 = DABS(T(1)) + DABS(T(2))
      RETURN
      END
