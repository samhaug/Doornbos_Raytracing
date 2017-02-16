c                                                                       exr02700
      subroutine rayps(nps,rmax,rmin,vk,nk,pp,t,delta,dddp,qamp,rlow,
     1iud,imult,nmod)
c        Input iud =tign(i) imult =imult(i) nmod included KANESHIMA S.Jul.14,1996
c     computes t,delta,dddp,qamp,rlow for ray parameter pp              exr02720
c     between levels rmax and rmin                                      exr02730
c     nps=1/2 for p/s                                                   exr02740
c     vk=velocity at rmax                                               exr02750
c     nk=layer number at rmax                                           exr02760
c     pp in s/radian, delta in radian                                   exr02770
c     model parameters in common blocks /model/,/dmod/                  exr02780
c     output - t,delta,dddp,qamp,rlow                                   exr02790
c     dddp=d(delta)/d(p)=gradient of delta versus ray parameter functionexr02800
c     qamp=integrated damping effect (equivalent to -t*)                exr02810
c     rlow=minimum level of ray piece (turning point, or rmin)          exr02820
c     in argument are pp,t,delta,dddp,qamp double precision             exr02830
c     modification of existing routine of unknown origin                exr02840
c                                                                       exr02850
c     subroutines called : none                                         exr02860
c     language fortran iv                                               exr02870
c     run on ibm                                                        exr02880
c     reference doornbos : asphericity and ellipticity corrections,     exr02890
c     in seismological algorithms (1988)                                exr02900
c                                                                       exr02910
c     common/model/r(99),v(99,2),q(99,2),rc,nlay                        exr02920
c     common/dmod/a(99,2),b(99,2)                                       exr02930
      include 'model.com'
c     double precision  a,b                                             exr02940
      double precision  pp,pe,t,delta,dddp,qamp                         exr02950
      double precision  eta,e,f,g,e0,f0,g0,vm,dt,fc                     exr02960
      double precision  vk
      logical  loop                                                     exr02970
      data  eps/1.e-06/,copi/57.29578/
      open(8,file='ray1d.out',status='old')
      loop=.true.                                                       exr02990
      delta=0.                                                          exr03000
      t=0.                                                              exr03010
      dddp=0.                                                           exr03020
      qamp=0.                                                           exr03030
      rlow=rmin                                                         exr03040
      pe=eps                                                            exr03050
      pe=dmax1(pp,pe)                                                   exr03060
      nl=nk                                                             exr03070
      eta=rmax/vk                                                       exr03080
 3000 format('ray parameter p (s/radian)= ',f8.3,' (s/deg)= ',f6.3) 
      write(8,2000) nmod, delta*copi, rmax, t, nl, nps, iud
    1 rlow=r(nl)                                                        exr03090
      if(eta-pp.lt.eps)  go to 9                                        exr03100
      f=dsqrt(eta*eta-pp*pp)                                            exr03110
      e=datan(f/pe)                                                     exr03120
      g=1./f                                                            exr03130
    2 n1=nl                                                             exr03140
      nl=nl-1                                                           exr03150
      if(nl.lt.1)  go to 4                                              exr03160
      if(r(nl)-rmin.lt.eps)  go to 5                                    exr03170
      eta=r(nl)/v(nl,nps)                                               exr03180
      if(r(n1)-r(nl).lt.eps)  go to 1                                   exr03190
      if(eta-pp.gt.eps)  go to 7                                        exr03200
    3 fc=1./(1.-b(n1,nps))                                              exr03210
      rlow=rc*((pp*a(n1,nps)/rc)**fc)                                   exr03220
      e0=0.                                                             exr03230
      f0=0.                                                             exr03240
      g0=0.                                                             exr03250
      loop=.false.                                                      exr03260
      go to 8                                                           exr03270
    4 if(rmin.lt.eps)  go to 3                                          exr03280
      vm=v(1,nps)                                                       exr03290
      go to 6                                                           exr03300
    5 vm=v(nl,nps)                                                      exr03310
      if(rmin-r(nl).gt.eps)  vm=a(n1,nps)*((rmin/rc)**b(n1,nps))        exr03320
    6 eta=rmin/vm                                                       exr03330
      if(eta-pp.lt.eps)  go to 3                                        exr03340
      rlow=rmin                                                         exr03350
      loop=.false.                                                      exr03360
    7 f0=dsqrt(eta*eta-pp*pp)                                           exr03370
      e0=datan(f0/pe)                                                   exr03380
      g0=1./f0                                                          exr03390
      fc=1./(1.-b(n1,nps))                                              exr03400
    8 delta=delta+(e-e0)*fc                                             exr03410
      dt=(f-f0)*fc                                                      exr03420
      t=t+dt                                                            exr03430
      write(8,2000) nmod, delta*copi, r(nl), t, nl, nps, iud
 1000 format(a3,f6.2,a3,f9.3,a3,f9.1,a4,i8) 
 2000 format(i8,f7.2,f11.1,f13.3,3i8) 
      dddp=dddp-(g-g0)*fc                                               exr03440
      qamp=qamp-dt/q(n1,nps)                                            exr03450
      e=e0                                                              exr03460
      f=f0                                                              exr03470
      g=g0                                                              exr03480
      if(loop)  go to 2                                                 exr03490
      if(pp.lt.eps)  delta=0.                                           exr03500
    9 return                                                            exr03510
      end                                                               exr03520
c                                                                       exr03530
