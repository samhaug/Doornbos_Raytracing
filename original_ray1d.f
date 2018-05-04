c     RAY1D    For spherical earth model. Modified from EXASP.f 
c     file contains main program, plus subroutine trefps(modified)    
c                                                                    
c     main program                                                  
c                                                                  
c     radial model is read in subroutine modps (from separate file)   
c        model parameters are put in common block /model/            
c     travel times are computed for rays with given ray parameter p
c      the number of rays between these points is nphase               
c      the precision of integration along a ray is controlled by the  
c      parameter stab. usually, stab.le.1.                           
c     each ray is specified by -                                    
c      p=initial ray parameter (in s/deg)                          
c      nmod=nr of ray pieces                                      
c     each ray piece is specified by -                           
c      (rm(i,j),j=1,2),imult(i),isp(i),itign(i),i=1,nmod        
c      rm(i,1)=rmax, rm(i,2)=rmin, i.e. ray stays between rmax and rmin 
c      imult (i)=nr of legs between rmax and rmin, isp(i)=1/2=p/s      
c      itign(i) specifies sense of direction of the ray               
c      itign(i).ge.1 if ray goes down from its upper endpoint (set to 1)
c      itign(i).lt.1 if ray goes up to its upper endpoint (set to -1)  
c     n.b : ray pieces must be input in the order followed by the ray  
c                                                                     
c     subroutines called : modps,trefps
c     language fortran iv                                           
c     run on ibm                                                    
c     reference doornbos : asphericity and ellipticity corrections, 
c     in seismological algorithms (1988)                            
c                                                                  
c     include 'model.com'
      common/model/r(120),v(120,2),q(120,2),rc,nlay           
      common/dmod/a(120,2),b(120,2)                                    
      double precision  a,b      
      double precision  p,t,dddp,qamp                                 
      double precision  delta
      dimension  rm(14,4),nrm(14),vrm(14),tign(14)
      dimension itign(14),imult(14),isp(14) 
c     integer  in/5/,ut/6/                                             
      integer  in, ut
      character*2 ps(2),sp(14)                                         
      data  ps(1),ps(2)/' p',' s'/                                   
      data  inmod/1/,inex/7/,delta/0.0/, in/5/, ut/6/
      data  pmin/0./,eps/.1e-06/,copi/57.29578/,nit/9/              
      data  nphmax/150/,nmodm/14/,pi/3.141593/,twopi/6.283185/       
      open(8,file='ray1d.out',status='old')
c                  
      call modps(inmod,ut)                                         
c                   
  100 read(in,*,end=200)  nphase,stab  
      if(nphase.lt.1)  nphase=1                        
      if(nphase.gt.nphmax)  nphase=nphmax             
      if(stab.lt.eps)  stab=1.                       
c     loop over phases                              
      do 300  iphase=1,nphase                      
c  
      read(in,*)  p,nmod
c                                                 
   20 deltag=delta*copi                          
c                                               
      if(nmod.gt.nmodm)  nmod=nmodm            
      read(in,*)  ((rm(i,j),j=1,2),imult(i),isp(i),itign(i),i=1,nmod)   
c                    
      nld=nlay/2                             
      pmax=0.                               
C
      do 13  i=1,nmod                      
C
        if(itign(i).gt.1)  itign(i)=1       
        if(itign(i).lt.1)  itign(i)=-1     
        tign(i)=itign(i)                  
c                     
        rm(i,3)=imult(i)                  
        rmax=rm(i,1)                     
        ips=isp(i)                      
        if(ips.lt.1)  ips=1            
        if(ips.gt.2)  ips=2                                 
        sp(i)=ps(ips)                                      
        rm(i,4)=ips                                       
        call select(ips,rmax,vrs,nls,nlay,nld)           
        if(pmax.lt.rmax/vrs)  pmax=rmax/vrs             
        nrm(i)=nls                                     
        vrm(i)=vrs                                    
   13 enddo                                         
      p=p*copi                                     
      gps=eps                                     
      ni=nit                                     
c                               
      call trefps(nmod,rm,nrm,vrm,delta,p,pmax,pmin,t,dddp,qamp,rlow,  
     + gps,ni,itign)  
C
c               
      p0=p/copi                                                    
      write(ut,*) 't=', t, ' p=', p, ' pmax=',pmax, ' pmin=', pmin
      write(ut,*) 't*=', -qamp
      wdelta = delta*copi
      write(ut,*) 'delta=', wdelta
c     
      p=p/copi 
      close (8)
c     call rpath(delta)
      write(ut,21)  p,t,wdelta
      write(ut,*) ''
   21 format('p= ',f7.4,'   t= ',f9.3,'  delta= ',f10.3) 
c                   
  300 continue   
  200 continue   
      end       
c              
c             
      subroutine trefps(nmod,rm,nrm,vrm,delta,p,pmax,pmin,t,dddp,q,rlow,
     +ep,ni,itign)
C
c               Input itign(i) included      KANESHIMA S. July 13, 1996
c     computes p,t,dddp,q,rlow - ray specified by nmod,rm,nrm,vrm,delta 
c     on input,p=initial ray parameter, within limits pmax and pmin    
c      ep=absolute error tolerance for delta,ni=max. nr. of iterations  
c     on return,p=final ray parameter,ep=delta(given)-delta(computed),  
c      ni=actual nr. of iterations                                      
c                                                                       
c     nmod=nr of ray pieces                                             
c     each ray piece is specified by -                                  
c      rm(i,j),j=1,4, nrm(i), vrm(i)                                    
c      rm(i,1)=rmax, rm(i,2)=rmin, rm(i,3)=nr of legs,rm(i,4)=1/2=p/s   
c      nrm(i)=layer nr at upper endpoint of ray piece                   
c      vrm(i)=velocity at upper endpoint of ray piece                   
c     t=travel time                                                     
c     dddp=d(delta)/d(p)=gradient of delta versus ray parameter function
c     q=integrated damping effect (equivalent to -t*)                   
c     rlow=minimum level of the ray (turning point, or lowest reflector)
c     t in seconds, delta in radians, p in s/rad.                       
c     in argument list are p,t,dddp,q double precision                 
c                                                                       
c     subroutines called : rayps                                        
c     language fortran iv                                               
c     run on ibm                                                        
c     reference doornbos : asphericity and ellipticity corrections,     
c     in seismological algorithms (1988)                                
c                                                                       
      dimension  rm(14,4),nrm(14),vrm(14),itign(14)     
      double precision  p,t,dddp,ad,p1,pshift,t1,del,d1ddp,dt           
      double precision  ddel,ddddp,d1d,a1d,dd,vk,q,q1amp,dqamp,delta    
      data  eps/1.e-06/                                                 
      n=0                                                               
      ad=100000.                                                        
      p1=p                                                              
      pshift=0.                                                         
    1 p1=p1+pshift                                                     
c 
c     write(6,*) 'p1=', p1, ' pshift=', pshift, ' delta=', delta
      if(p1.gt.pmax)  p1=pmax                                          
      if(p1.lt.pmin)  p1=pmin                                          
      t1=0.                                                            
      del=0.                                                           
      d1ddp=0.                                                         
      q1amp=0.                                                         
      r1low=100000.                                                    
c                   Loop for ray segments (or ray pieces)
      do 2  i=1,nmod                                                  
       rmax=rm(i,1)    
       rmin=rm(i,2) 
       xmult=rm(i,3) 
       imult = rm(i,3)
       nps=1  
       if(rm(i,4).gt.1.5)  nps=2 
       vk=vrm(i)   
       nk=nrm(i)    
       iud = itign(i)
       write(6,*)'nps=',nps,' rmax=',rmax,' rmin=',rmin,' vk=', vk
c      write(6,*) 'nk=', nk,' p1=',p1 
       call rayps(nps,rmax,rmin,vk,nk,p1,dt,ddel,ddddp,dqamp,rll,iud,
     1 imult,i) 
c               Input iud imult i included    KANESHIMA S. July 14, 1996
c      write(6,*)' delta=',delta,' del=',del
       t1=t1+xmult*dt   
       del=del+xmult*ddel  
       delta = del
       d1ddp=d1ddp+xmult*ddddp 
       q1amp=q1amp+xmult*dqamp  
    2 r1low=amin1(r1low,rll)                                            exa02830
      d1d=delta-del                                                     exa02840
      a1d=dabs(d1d)                                                     exa02850
      n=n+1                                                             exa02860
      if(a1d.lt.ad)  go to 3                                            exa02870
      pshift=pshift/2.                                                  exa02880
      if(n-ni) 1,4,4                                                    exa02890
    3 ad=a1d                                                            exa02900
      dd=d1d                                                            exa02910
      t=t1                                                              exa02920
      dddp=d1ddp                                                        exa02930
      q=q1amp                                                           exa02940
      rlow=r1low                                                        exa02950
      pshift=0.                                                         exa02960
      if(dabs(dddp).gt.eps)  pshift=dd/dddp                             exa02970
      p=p1                                                              exa02980
      ni=n                                                              exa03000
    4 ep=dd                                                             exa03010
      return                                                            exa03030
      end                                                               exa03040
