c                                                                       exm01910
      subroutine modps(m,n)                                             exm01920
c                                                                       exm01930
c     read in radial model parameters : radius r, velocity v, damping q exm01940
c     compute parameters a and b of v=a*(r/rc)**b law, rc=reference r   exm01950
c     print model parameters if flag is set (ind.gt.0)                  exm01960
c     m specifies input file, n output file                             exm01970
c     input : rc,ind,r(i),v(i,j),q(i,j); j=1,2 for p,s waves            exm01980
c     output: nlay,a(i,j),b(i,j)                                        exm01990
c     nlay=number of layers, maximum 99                                 exm02000
c     r,v,q,rc and nlay are put in common block /model/                 exm02010
c     a and b are double precision and are put in common block /dmod/   exm02020
c                                                                       exm02030
c     subroutines called : none                                         exm02040
c     language fortran iv                                               exm02050
c     run on ibm                                                        exm02060
c     reference doornbos : asphericity and ellipticity corrections,     exm02070
c     in seismological algorithms (1988)                                exm02080
c                                                                       exm02090
c     include 'model.com'
      common/model/r(120),v(120,2),q(120,2),rc,nlay           
      common/dmod/a(120,2),b(120,2)                                    
      double precision  a,b      
c
      double precision  rlg,vlg                                         exm02130
      dimension  vi(2),vvi(2)                                           exm02140
      logical lind                                                      exm02150
c     data  eps/1.e-06/,nmax/99/,qmax/10000./                           exm02160
      data  eps/1.e-06/,nmax/120/,qmax/10000./    
      open(m,file=
     +  '/home/samhaug/Doornbos_raytracing/rmod.dat',status='old')
c                                     modified   by S. K.  1995 Sep. 9
c     read reference radius rc, flag to print model (ind.gt.0)          exm02170
      read(m,1)  rc,ind                                                 exm02180
    1 format(f10.3,i5)                                                  exm02190
      lind=(ind.gt.0)                                                   exm02200
c     read model, compute a and b in v=a*(r/rc)**b                      exm02210
      if(lind)  write(n,2)                                              exm02220
    2 format('1model parameters'/'0number',5x,'r',8x,'vp',7x,'a',12x,'b exm02230
     1',11x,'q',7x,'vs',7x,'a',12x,'b',11x,'q')                         exm02240
      do 4  i=1,nmax                                                    exm02250
      read(m,3,end=5)  r(i),v(i,1),q(i,1),v(i,2),q(i,2)                 exm02260
c   3 format(5f10.2)                                                    exm02270
    3 format(f12.2,f8.3,f12.2,f8.3,f12.2)   
c                          modified by Kneshima Satoshi  1995 Sep. 9 
    4 continue                                                          exm02280
      i=i+1                                                             exm02290
    5 nlay=i-1                                                          exm02300
      i=1                                                               exm02310
      do 6  j=1,2                                                       exm02320
      if(q(i,j).lt.eps)  q(i,j)=qmax                                    exm02330
      vvi(j)=v(i,j)                                                     exm02340
      if(vvi(j).lt.eps)  vvi(j)=qmax                                    exm02350
      b(i,j)=0.                                                         exm02360
    6 a(i,j)=v(i,j)                                                     exm02370
      if(lind)  write(n,7)  i,r(i),(v(i,j),a(i,j),b(i,j),q(i,j),j=1,2)  exm02380
    7 format(1x,i4,f11.2,2(f9.3,2d13.5,f8.1))                           exm02390
      rri=r(i)                                                          exm02400
      do 11  i=2,nlay                                                   exm02410
      ri=r(i)                                                           exm02420
      do 10  j=1,2                                                      exm02430
      if(q(i,j).lt.eps)  q(i,j)=qmax                                    exm02440
      vi(j)=v(i,j)                                                      exm02450
      if(vi(j).lt.eps)  vi(j)=qmax                                      exm02460
      if(ri-rri.gt.eps)  go to 9                                        exm02470
      if(ri.lt.rri)  write(n,8)                                         exm02480
    8 format(' error in model')                                         exm02490
      b(i,j)=0.                                                         exm02500
      a(i,j)=v(i,j)                                                     exm02510
      go to 10                                                          exm02520
    9 rlg=ri/rri                                                        exm02530
      vlg=vi(j)/vvi(j)                                                  exm02540
      b(i,j)=dlog(vlg)/dlog(rlg)                                        exm02550
      a(i,j)=v(i,j)*((rc/r(i))**b(i,j))                                 exm02560
   10 vvi(j)=vi(j)                                                      exm02570
      if(lind)  write(n,7)i,ri,(v(i,j),a(i,j),b(i,j),q(i,j),j=1,2)      exm02580
   11 rri=ri                                                            exm02590
      if(lind)  write(n,12)  rc                                         exm02600
   12 format('0normalizing reference radius',f10.3)                     exm02610
      return                                                            exm02620
      end                                                               exm02630
c                                                                       exm02640
      subroutine select(ips,rmax,vrs,nls,n1,n2)                         exm02650
c                                                                       exm02660
c     calculates vrs and nls for given rmax,n1 is somewhere above nls   exm02670
c     n2 is search step.  vrs=vp if ips=1, vrs=vs if ips.gt.1           exm02680
c     uses parameters in common blocks /model/ and /dmod/               exm02690
c     modification of existing routine of unknown origin                exm02700
c                                                                       exm02710
c     subroutines called : none                                         exm02720
c     language fortran iv                                               exm02730
c     run on ibm                                                        exm02740
c     reference doornbos : asphericity and ellipticity corrections,     exm02750
c     in seismological algorithms (1988)                 
c                                                         
C        n2:  search step to be input            2008/11/20
C
C     common/model/r(120),v(120,2),q(120,2),rc,nlay    
C     common/dmod/a(120,2),b(120,2)   
C     double precision  a,b          
      include 'model.com'
      data  eps/1.e-06/                        
C     write(6,*) 'SELECT: ', ips,rmax,vrs,nls,n1,n2
      i=1                                                               exm02820
      if(ips.gt.1)  i=2                                                 exm02830
      nl=n1                                                             exm02840
      nld=n2                                                            exm02850
    1 if(r(nl)-rmax.gt.eps)  go to 2                                    exm02860
      if(rmax-r(nl).lt.eps)  go to 7                                    exm02870
      nl=nl+1                                                           exm02880
      go to 1                                                           exm02890
    2 if(rmax.gt.r(1))  go to 4                                         exm02900
      vrs=v(1,i)                                                        exm02910
      nls=1                                                             exm02920
      go to 10                                                          exm02930
    3 nld=max0(1,nld/2)                                                 exm02940
C     write(6,*) 'SELECT: ', rmax, r(nl), r(nl+1)
      if(rmax.ge.r(nl))  go to 5                                        exm02950
    4 nl=nl-nld                                                         exm02960
      go to 3                                                           exm02970
    5 if(rmax.lt.r(nl+1))  go to 6                                      exm02980
      nl=nl+nld                                                         exm02990
      go to 3                                                           exm03000
    6 if(rmax-r(nl).gt.eps)  go to 9                                    exm03010
    7 if(nl.lt.2)  go to 8                                              exm03020
      if(rmax-r(nl-1).lt.eps)  nl=nl-1                                  exm03030
    8 vrs=v(nl,i)                                                       exm03040
      nls=nl                                                            exm03050
      go to 10                                                          exm03060
    9 nl=nl+1                                                           exm03070
      vrs=a(nl,i)*((rmax/rc)**b(nl,i))                                  exm03080
      nls=nl                                                            exm03090
   10 return                                                            exm03100
      end                                                               exm03110
c
c  
      subroutine setvln(nlay)
c
c   Loop to find layer number (nrm) and velocity (vrm) at upper endpoint
c           of each ray pieces (segments)
c                                                 KANESHIMA Satoshi 11/23/1996
      include 'raygm.com'
      character*2 ps(2),sp(6)
      data  ps(1),ps(2)/' p',' s'/
c
      nld=nlay/2
      do i=1,nrseg
         if(itign(i).gt.1)  itign(i)=1
         if(itign(i).lt.1)  itign(i)=-1
         tign(i)=itign(i)
         rm(i,3)=imult(i)
         rmax=rm(i,1)
         ips=isp(i)
         if(ips.lt.1)  ips=1
         if(ips.gt.2)  ips=2
         sp(i)=ps(ips)
         rm(i,4)=ips
c                    write(6,*)'SETVLN nld:   ',nld,'  nlay:   ', nlay
         call select(ips,rmax,vrs,nls,nlay,nld)
         nrm(i)=nls
         vrm(i)=vrs
      enddo
      return
      end
