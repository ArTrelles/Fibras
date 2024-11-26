      program PopSynth
      implicit none
      integer MNSc,MN2D,MN3D,MN,MNiso
      parameter (MN=5000000)
      !parameter (MN=361115) ! *** Max number of particles in input file
      !parameter (MN=361115)
      parameter (MNSc=4096) !Number of grid points in extinction map Schegel
      parameter (MN3D=1024) !Number of grid points in extinction map 3D
      parameter (MN2D=1024) !Number of grid points in extinction map 2D
      parameter (MNiso=34)  !Number of isochrones in isochromes file
c from  IsoNo.txt (name of isochrone as function of metallicity)      
      integer i,niso(0:MNiso-1),Ns,j,k,jindex(0:6200),npr
     &     ,IsoNo(0:MNiso-1),l,q,inputdata,GenSel,per3
     &     ,IsoStart(0:MNiso-1),ageindex,np,m,per2,lbegin,lend,ngen
      integer dummyint,nparts
      double precision MiMax,Afac,z,PAge,vgen,igen,avgen,aigen,z2
      double precision bgen,abgen,bvgen,exb
      double precision frac,Vigen,dummy,dmod,mag2,rd,circ,z1,mtemp,mp
      double precision Age(0:MN-1),logL_Lo(0:MN-1),logG(0:MN-1),zsh,
     &     Mi(0:MN-1),Mc(0:MN-1),logT(0:MN-1),
     &     C_O(0:MN-1),VMU(0:MN-1),
     &     VMB(0:MN-1),VMV(0:MN-1),VMR(0:MN-1),VMI(0:MN-1),
     &     VMJ(0:MN-1),VMH(0:MN-1),VMK(0:MN-1),
     &     mbol(0:MN-1)

      double precision Teff, surface_gravity
      double precision vl_p,vb_p,vlo_p,vbo_p,vrots,mtn,Nsf,z3,z4
      double precision rshift, FeH, alphaFe
      double precision M_hec,period,pmode,logMdot,intIMF,vcir
      double precision mtotalneed
      double precision d_p,Milb,mag,mgen,sumgen,R,tss,month
      double precision x_p,y_p,z_p,vx_p,vy_p,vz_p,m_p,ex,ex31,ex148
      double precision x_pn,y_pn,z_pn,vx_pn,vy_pn,vz_pn
      integer il,ib,ie1,ie2,ie3,ie4,ie5,ie6,ie7,ie8
      integer il3,ib3,id3,jflag
      double precision dl,db,wl1,wl2,wb1,wb2,wd1,wd2,exS,ex3
      double precision dl3,db3,dd3,CAfactor
      double precision leS(0:MNSc-1),beS(0:MNSc-1),le3(0:MN3D-1),
     &     be3(0:MN3D-1),de3(0:99),datS(0:16777215),dat3(0:104857599)
      double precision l_p,b_p,a(6),ao(6),ae(6),p(2),po(2),pe(2), 
     &     ap(4),apo(4),ape(4)
      double precision M_PI,t0,uh,vh,wh,makms,mdkms,makmso,mdkmso
      double precision uhn,vhn,whn
      parameter (M_PI=4.0d0*datan(1.0d0)) !pi value
      double precision alpha,delta,g,gb,gl,GRVS
      double precision Uo,Vo,Wo,xo_p,yo_p,zo_p,vxo_p,vyo_p,vzo_p,
     &     vro_p,vroto_p,rop,od_p,rp,vr_p,vrot_p
      double precision lo_p,bo_p
      double precision d_pnew,fsample
      INTEGER  Ngrid_i,Ngrid_f,Np_2 ! grid y particulas
      INTEGER  it,l_2,il2,jl2,kl2 ! iteradores
      INTEGER  max_cols,col ! columnas
      INTEGER  lim_inf,lim_sup
      PARAMETER (lim_inf=-200,lim_sup=200)
      real*8 delta_cell
      real*8 nip
      INTEGER ip,jp,kp
      !valores auxiliares para calcular rms
      real*8 v1,v2,v3,v4,v5,v6
      !PARAMETER(Np_2=249307,max_cols=13,Ngrid_i=-50,Ngrid_f=50)
      PARAMETER(Np_2=861787,max_cols=14,Ngrid_i=-200,Ngrid_f=200)
      real*8 xp,yp,zp,Resf
      real*8 vxp,vyp,vzp
      real*8 xms,yms,zms
      real*8 vxms,vyms,vzms
      real*8 nb_ind
      double precision gr_n,r_n

c     real*8 VXs,VYs,VZs 
      INTEGER seed
      REAL (Kind=8), DIMENSION (:,:,:), ALLOCATABLE :: grid,nb
      REAL(Kind=8), DIMENSION (:,:,:), ALLOCATABLE :: Bx,By,Bz
      REAL(Kind=8), DIMENSION (:,:,:), ALLOCATABLE :: Bvx,Bvy,Bvz
      REAL(Kind=8), DIMENSION (:,:,:), ALLOCATABLE :: xms2,yms2,zms2
      REAL(Kind=8), DIMENSION (:,:,:), ALLOCATABLE :: vxms2,vyms2,vzms2
      REAL (Kind=8), DIMENSION(:,:), ALLOCATABLE :: coordinates
      real GAUSS,RANDd ! Generador gaussiano
c Variables definidas por Arturo para recuperar gr y r
      double precision r7,gr7,z7
      real*8 ran2
c   Distancia al sol
      real*8  rpn
      real*8 prio2,priot  
      integer prio1,iden

      character*40 FN
      character*25 ONAME
      character*1 check
      character*1 hash
      include 'parameters.h'
      include 'Mockinioptions.h'

c      *** For simpler file ***
      open(30,file='ini/input_200.dat',
     &   status='unknown')
  
      ALLOCATE(coordinates(Np_2,max_cols))

     
      !OPEN(15, FILE='nuevo_archivo_n.txt',status='unknown')

      DO it = 1,Np_2
        READ(30,*) (coordinates(it,col),col=1,max_cols)
      END DO

      delta_cell = 1.
      nip=1.0000000E+00
      v1=0.000000E+0
      v2=0.000000E+0
      v3=0.000000E+0
      v4 =0.000000E+00
      v5 =0.000000E+00
      v6 =0.000000E+00

      fsample=0.1d0


     

      REWIND(30)

c *** GEnrating MOck DAta From Sims GEMODAFS version 1.0 31/08/16 ***

c *** For ASCII input file ini/input.dat ***
     
c *** TEMP particle ***
      open(20,FILE='ini/UBV/IsoNo.txt', STATUS='OLD')
      do i=0,MNiso-1
         read(20,'(I10)') IsoNo(i) ! Read isochrones metallicity*1e6
      enddo
      close(20)
     
      per2=0
      per3=0
      ngen=0
      GenSel=0
c *** Read extinction arrays ***
      if (EXTINCT.eq.1) then
         if (Schlegel.eq.1) then !Schlegel extinction array
c table from http://irsa.ipac.caltech.edu/applications/DUST/
            open(51,file='ini/ScEx1.bin',status='old'
     &           ,form='unformatted')
            write(6,*) 'Reading Extinction arrays'
            do i=0,16777215
               read(51) datS(i)
            enddo

c Define londitud, latitud bins (gal coords) for extinction
            dl=(2.0d0*M_PI)/4095.0d0 
            db=M_PI/4095.0d0   
            write(6,*) 'Schlegel extinction array Read'
         else
            open(51,file='ini/2dEx1.bin',status='old'
     &           ,form='unformatted') !2D analytic extinction
            write(6,*) 'Reading Extinction arrays'
            do i=0,1048575
               read(51) datS(i)
            enddo
             write(6,*) '2D analytical extinction array Read'

c Define londitud, latitud bins (gal coords) for extinction
             dl=(2.0d0*M_PI)/1023.0d0 
             db=M_PI/1023.0d0   
         endif
         close(51)
         open(51,file='ini/3dEx1.bin',status='old',form='unformatted')
         do i=0,104857599
            read(51) dat3(i)
         enddo
         close(51)
         write(6,*) '3D factor array Read'
         write(6,*) 'Extinction running'

c Define longitud, latitud and distance bins (gal coords) for extinction
c     0<l<2*pi 
         dl3=(2.0d0*M_PI)/1023.0d0  
c     -pi/2<b<pi/2, with +pi/2 for +ve
         db3=M_PI/1023.0d0    
         dd3=3.47712123/99.0d0 
      else
         write(6,*) '*** No Extinction Added ***'
      endif
c *** Read in (Padova) Isochrone Files ***   
      j=0
      k=1
      jindex(0)=0
      do i=0,MNiso-1
         write(FN,'(a18,i5.5,a4)')'ini/UBV/output_0.0',IsoNo(i),'.dat'
         open(20,FILE=FN,STATUS='OLD')
c *** 10 headers read and discarded ***
         do l=1,9
            read(20,'(A150)')
         enddo
         hash='#'
c *** Read data ***
         niso(i)=0
         write(ONAME,'(a12,i5.5,a4)')'ini/UBV/V-VI',IsoNo(i),'.dat'
         do l=0,30000
            read(20,*,end=600) check
            if (index(check,hash).ne.1) then
               backspace(20)
               read(20,*) Age(j),Mi(j),Mc(j),logL_Lo(j),logT(j)
     &              ,logG(j),mbol(j),VMU(j),VMB(j),VMV(j),VMR(j)
     &              ,VMI(j),VMJ(j),VMH(j),VMK(j),C_O(j),M_hec
     &              ,period,pmode,logMdot,intIMF
               j=j+1
               niso(i)=niso(i)+1
            else
               read(20,*,end=600) check
               jindex(k)=j
               k=k+1
            endif
         enddo 
  600    close(20)
         IsoStart(i)=j
      enddo
      write(6,*) 'Isochrones Read'
 1940 format(17(1pE13.5),I10,3(1pE13.5))
c *** Open files to read & write ***
c *** open output file (binary or ASCII) ***
      if (Binary.eq.0) then
       open(10,file='output/Generated_testzB3.dat',status='unknown')
c     &     form='unformatted')
c         open(15,file='ini/input.dat',status='unknown')
      else
      open(10,file='output/GeneratedStarslebig.bin',status='unknown'
     &     ,form='unformatted')
      endif


c      write(10, '(A)') 'col1 col2 col3 col4 col5 col6 col7 col8 col9 col
c     &10 col11 col12 col13 col14 col15 col16 col17 col18 col19 col20 col
c     &21 col22 col23 col24 col25 col26 col27 col28 col29 col30 col31 col
c     &32 col33 col34 col35 col36 col37 col38 col39 col40 col41'

      !write(10,*) 'TARGETID     ','RA   ','DEC  ',
      !&      'PRIORITY    ','SUBPRIORITY','TYPE_ID'
c *** Open input file ***
c *** For example file ***
c      open(30,file='ini/lbselsl90d5.dat',status='unknown')
c 2030 format(16(1pE13.5),I10)
!c *** For simpler file ***
!      open(30,file='ini/input.dat',status='unknown')
c 2030 format(9(1pE13.5))
      np=50000000
      !np=361115 !max number of particles
      !np=361115
      sumgen=0.0d0
      mtotalneed=0.0d0
      nparts=0
c   16 format(9(1pE13.5)) 



      seed=123
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do q=0,np-1
         mtemp=0.0d0
c *** For example file ***
c         read(30,2030,end=603) dummy,dummy,dummy,x_p,y_p,z_p,dummy,
c     &        dummy,vx_p,vy_p,vz_p,dummy,dummy,z,tss,m_p,dummyint
c *** x,y,z(kpc),vx,vy,vz(km/s),metallicity(z),age(Gyr),mass(solar)
         !read(30,*,end=603) x_p,y_p,z_p,vx_p,vy_p,vz_p,m_p,tss,
         !+zsh,z3,z4,vcir,z1
c        data.columns = ['x', 'y', 'z', 'vx', 'vy', 'vz', 'm_p', 'Fe/H',
c        'alpha/Fe', 'tss', 'rshift', 'z3', 'z4','vcirc']
         
         read(30,*,end=603) x_p,y_p,z_p,vx_p,vy_p,vz_p,m_p,FeH,
     &   alphaFe,tss,rshift,z3,z4,vcir
         z=z4*z3
c         write(15,16) x_p,y_p,z_p,vx_p,vy_p,vz_p,z,tss,m_p
         nparts=nparts+1
         x_p=((x_p)/100.0d0)
         y_p=((y_p)/100.0d0)
         z_p=((z_p)/100.0d0) !U100 kpc
         ! En esa coordenada se estimara la densidad de EF

         vx_p=(vx_p+usun)/VUKMS
         vy_p=(vy_p+vlsr+vsun)/VUKMS
         vz_p=(vz_p+wsun)/VUKMS
         z=z*1e6               ! Metallicity in isochrone filename is in
c z*10^6

c *** if z is lower than the one in the first isochrone file, assigns the
c lower metallicity in files, if higher, the highest ***

         if (z.lt.100.0d0) then
            z=100.0d0
         else if (z.gt.30000.0d0) then
            z=30000.0d0
         endif 

c *** definitions of usefull variables and changes of units ***

         rp=dsqrt((x_p+xsun)**2+(y_p+ysun)**2)
         vr_p=((vx_p*(x_p+xsun)+vy_p*(y_p+ysun)+vz_p*(z_p+zsun))/rp)*
     +VUKMS
         vrot_p=(-1.0d0*((x_p+xsun)*vy_p-(y_p+ysun)*vx_p)/rp)*VUKMS
         rp=dsqrt((x_p+xsun)**2+(y_p+ysun)**2+(z_p+zsun)**2)
         d_p=dsqrt(x_p**2+y_p**2+z_p**2)
c         dmod=5.0d0*(dlog10(d_p*100000.0d0)-1.0d0)! distance term for 
c        V-Mv=5(log(r)-1)
         dmod=5.0d0*(dlog10(d_p*100000.0d0)-1.0d0)
         PAge=log10(tss*1000000000.0d0) !log(age) in yrs

c *** Cut in maximum radius

         if(rp.lt.150.0/100.) then

c *** Find approproate isochrone segment ***
c *** Random age & metallicity if unknown ***
c     z=(ran1(idum)*0.9d0+0.1d0)*19000.0d0
c     PAge=log10(ran1(idum)*12000000000.0d0)
         ageindex=nint((PAge-6.6d0)/0.02d0)

c limits in max and min allowed ages (depend on isochrones)
         if (ageindex.lt.0) then 
            ageindex=0
         else if (ageindex.gt.181) then
            ageindex=181
         endif

c *** extinction (CTIO) ***
         if (EXTINCT.eq.1) then
c *** Make grids l,b,d ***
            if (Schlegel.eq.1) then
               do i=0,MNSc-1
                  leS(i)=i*dl
                  beS(i)=i*db-M_PI/2.0d0
               enddo
            else
            do i=0,MN2D-1
                  leS(i)=i*dl
                  beS(i)=i*db-M_PI/2.0d0
            enddo
            endif
            do i=0,MN3D-1
               le3(i)=i*dl3
               be3(i)=i*db3-M_PI/2.0d0
            enddo
            do i=0,99
               de3(i)=i*dd3-2.0d0
            enddo


c ***  


c *** l,b,d from cartesian input ***
            l_p=-datan2(y_p,-x_p)+M_PI
            b_p=dasin(z_p/d_p)
            vrots=(-1.0d0*(x_p*vy_p-y_p*vx_p)/d_p)*VUKMS
            vl_p=vrots/(4.74d0*(d_p*1.0e5))
            vb_p=vz_p*VUKMS/(4.74d0*(d_p*1.0e5))

c *** Find il,ib,id (position in the extinction map grid)
c     for each particle in the input ***
            il=int(l_p/dl)
            ib=int((b_p+(M_PI/2.0d0))/db)
            il3=int(l_p/dl3)
            ib3=int((b_p+(M_PI/2.0d0))/db3)
            id3=int((dlog10(d_p*100.0d0)+2.0d0)/dd3)
            if (float(id3).gt.99.0) then
               id3=99
            else if (float(id3).lt.0.0) then
               id3=0
            endif
c *** Find weighting (extinction interpolated between contiguous
c     extinction bins to the real position of particle in 
c     extinction Schegel or 2D bin) ***
 c *** Assigning isochrones to particle depending on metals
c *** total age bins = 181
c All isochrones has been stored in same array, now select the entire
c isochrone for the proper metals (lbegin to lend)
         MiMax=0.0d0
         do i=1,MNiso-1
            if (z.le.dble(IsoNo(i))) then
             
c *** Calc MiMax: maximum stellar mass to take in the isochrone 
c                 depending on the age of the population (sim particle) ***
c WARNING: We are taking only the initial mass of the star, not accounting
c for mass loss in stellar evolution (Mc(l)=mass at present).
                  if (Mi(l).gt.MiMax) then
                     MiMax=Mi(l)
                  endif
               enddo
c *** Calc constants for the IMF (imftype=1 is Salpeter)
               if(imftype.eq.1) then
               Afac=-m_p*0.35d0/(MiMax**(-0.35d0)-Msmin**(-0.35d0))
               endif
c *** Calc Milb that is the minimum star mass limit for stars that will contribute
c to the total magnitud (not interested into fainter stars that
c     will not be seen.
c We take advantage that in the isochrone library star contribution is 
c ordered by mass from low to high.
               Milb=10.0d0**12.0d0 !A very high number just in case
                !Interpolation between
                     Milb=(Mi(l)-Mi(l-1))*frac+Mi(l-1)!present and previous to
c the limit, mass.
                     goto 667
                  endif

c If star particle is too far to be observed (The more massive 
c stars when splitting the stellar population have Vmag>limit)
c selected min mass will be higher than selected max mass, then
c no stars are generated.

               enddo
            
c *** Calc number to generate *** (new part)
  667          continue
             if(imftype.eq.1) then
               Nsf=((Afac/(-1.35d0))*(MiMax**(-1.35d0)-Milb
     &              **(-1.35d0)))
               Ns=int(Nsf)
c select randomly if we choose one more or one less star (Nsf may not be
c an integer)
               if (ran2(idum).lt.(Nsf-Ns)) then
                  Ns=Ns+1
               endif

c Control of the total stellar mass generated if Nstars=Ns
               mtn=(Afac/(-0.35d0))
     &              *(MiMax**(-0.35d0)-Milb**(-0.35d0))
c Control for Ns=0 (if error on m_p or particle too far, too faint)
               if (Ns.lt.0.0d0) then
                  Ns=0
                  mtn=0.0d0
               endif
              endif
               mtotalneed=mtotalneed+mtn
               if (Ns.eq.0) then
                  goto 601
               endif
c *** Generate stars ***


                 rpn = sqrt((x_pn-xsun*100.0d0)**2.
     &          +(y_pn-ysun*100.0d0)**2.
     &           +(z_pn-zsun*100.0d0)**2.)     
               

ccc               if(r7.le.20.and.r7.ge.12.and.
ccc     &          gr7.gt.(-0.7).and.gr7.lt.2.5)then
c                 hasta este punto ao(3) viene en arcsec       
                  call CriteriaTar(idum,rpn,r7,gr7,z7,ao,ae,apo,
     &                  po,priot,iden)
                   !POr consistencia de unidades es mejor regresarlas a
                   ! arcsec
                       ao(3) = ao(3)/1.d3
                       ae(3) = ae(3)/1.d3
                       !write(6,*) ao(4),'ao4'
                       ao(4) = ao(4)/1.d3
                       ao(5) = ao(5)/1.d3
                   
c                  call  UVWH_TO_CARTE(Uo,Vo,Wo,VXs,VYs,VZs)


!G_m = data[:,20]

!G_mag = data[:,22]
color_GBRP=data[:,23]
!G_error = data[:,24]
c_err = data[:,25]

!parallax_error  ae(3)
!paraobs  ao(3)    
!snr abs(ae(3)/ao(3)).lt.0.3
!sncerr = np.abs(color_GBRP / c_err)

!valid_parallax = (sncerr > 10) & (np.abs(G_error) < 0.5 ) &
!(parallax_error < 0.5) & (np.abs(parallax_error / paraobs) < 0.2)



ccc                    if(iden.ne.8.and.iden.ne.0)then
                     if(pe(1).lt.0.5)then
                     if(r7.ge.16.and.r7.le.19)then
                     if(gr7.lt.2.5.and.gr7.gt.-1)then
                     if(ae(3)*1.d3.lt.0.5)then
                     if(abs(ae(3)/ao(3)).lt.0.1)then 
                     if(abs(po(2)/pe(2)).gt.10)then
                     if(iden.gt.0)then 
                                    
                              if (Binary.eq.0) then
c *** Ascii ***   
                        !write(6,*) x_pn, xo_p      
c                       if(vxo_p.gt.50)then 
c                        write(6,*) xo_p,vxo_p 
c                       endif
                        
                          !if(iden.eq.8)then
                          !   write(6,*) iden,m
                          ! endif
c                                     write(10,*) 
 
                        

                      write(10, '(6(F20.16," "),6(F20.16," "),2(F20.16," 
     &"),2(F20.16," "),2(F20.16," "),2(F20.16," "),F20.16," ",F20.16," "
     &,F20.16," ",F25.16," ",F25.16," ",F20.16," ",F20.16," ",F20.16," "
     &,F20.16," ",F25.16," ",F25.16," ",F20.16," ",F20.16," ",I10,F20.16
     &," ",F20.16," ",F20.16," ",F20.16," ",F20.16," ",F20.16," "F20.16)
     &')(ao(j),j=1,6),(ae(j),j=1,6),(po(j),j=1,2),(pe(j),j=1,2),(apo(j),
     &j=1,2),(ape(j),j=1,2),xo_p,yo_p,zo_p,vx_pn,vy_pn,vz_pn,vxo_p,vyo_p
     &,vzo_p,mgen,z,PAge,priot,iden,r7,gr7,z7,lo_p,bo_p,r_n,gr_n
c                                    write(10,*)
c     &                             ,vgen,d_p,x_p,y_p,z_p,vx_p,vy_p,vz_p
c     &                             ,mgen 
                        
c *** Binary ***
                             else
                                  write(10) 
     &(ao(j),j=1,6),(ae(j),j=1,6),(po(j),j=1,2),(pe(j),j=1,2),(apo(j),
     &j=1,2),(ape(j),j=1,2),xo_p,yo_p,zo_p,vx_pn,vy_pn,vz_pn,vxo_p,vyo_p
     &,vzo_p,mgen,z,PAge,priot,iden,r7,gr7,z7,lo_p,bo_p,r_n,gr_n

c                                write(10) (a(j),j=1,6),(ao(j),j=1,6)
c     &                             ,(ae(j),j=1,6),(p(j),j=1,2)
c     &                             ,(po(j),j=1,2),(pe(j),j=1,2)
c     &                             ,vgen,d_p,x_p,y_p,z_p,vx_p,vy_p,vz_p
c     &                             ,mgen,GRVS,VIgen,ex31,avgen
c     &                             ,avgen-aigen,z,PAge,od_p,xo_p,yo_p
c     &                             ,zo_p,vxo_p,vyo_p,vzo_p,vro_p,vroto_p
c     &                             ,rp,vr_p,vrot_p,rop,makms,mdkms
c     &                             ,makmso,mdkmso,vl_p,vb_p,vlo_p,vbo_p
c     &                             ,l_p,b_p
                              endif
                             GenSel=GenSel+1
                         endif              
                       endif
                      endif
                    endif
                  endif
                 endif
              endif
                     !       endif
                          !endif
                          !endif
c *** END CUSTOM SELECTIONS ***
                           d_p=d_p/LUKPC
                        endif
c                 este endif termina la condicion para fsample
                        endif
                        goto 602
                     endif
                  enddo
  602             sumgen=sumgen+mgen
                  mtemp=mtemp+mgen
                  ngen=ngen+1
               enddo
               goto 601
            endif
         enddo
c *** Percentage done ***
  601    if (int((float(q)/float(nparts))*100.0d0).gt.per2) then
            per2=int((float(q)/float(nparts))*100.0d0)
            if (inputdata.lt.2) then
               write(6,*) per2,'% complete'
            endif
         endif
       endif
      enddo 
  102 format(63(1pE13.5))
      close(11)

      DEALLOCATE(coordinates)
      DEALLOCATE(nb)
      DEALLOCATE(Bx,By,Bz)
      DEALLOCATE(Bvx,Bvy,Bvz)

      !valores sigma
      DEALLOCATE(xms2,yms2,zms2)
      DEALLOCATE(vxms2,vyms2,vzms2)

c *** Finish ***
  603 write(6,*) '100 % complete'
      write(6,*) 'You generated',sumgen,' solar masses of stars!'
      write(6,*) 'You needed',mtotalneed,' solar masses of stars!'
      if (sumgen.gt.mtotalneed) then
         write(6,*) 'You have',sumgen-mtotalneed,' excess'
      else
         write(6,*) 'You have',mtotalneed-sumgen,' deficit'
      endif
  
      write(6,*) 'You generated',ngen,' stars',GenSel,' selected'
      close(30)
      close(10)


      END


  
