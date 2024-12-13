      program nn
!                                                                       
!     NEAREST NEIGBOUR DISTANCES                                        
!         (USING LAPW STRUCTURE FILE)                                   
!     WRITTEN BY K.SCHWARZ AND P.BLAHA                                  
!     OCTOBER 1988, QTP, U.OF FLORIDA, GAINESVILLE                      
!                                                                       
      use struk
      use variable_fields
      IMPLICIT REAL*8 (A-H,O-Z)
!
!     Constant parameter definition
!
!c      INCLUDE 'param.inc'
!      PARAMETER (NATO= 64)                                              
      PARAMETER (NNN=100000)                                              
!      PARAMETER (NDIF= 99)                                              
      PARAMETER (mshell= 100)    
!     NATO    NUMBER OF  NONEQUIVALENT ATOMS            
!     NNN     NUMBER OF NEIGHBOURS                                      
!     NDIF    NUMBER OF ALL ATOMS IN UNITCELL                           
!     MSHELL  Max. Nr. of SHELLS
!                                              
!      parameter (dlimit=1.d-4)                                          
!      LOGICAL           ovlap(NATO),used(nato*ndif)
      LOGICAL           ORTHO,error
      LOGICAL           SHOWALL
!      CHARACTER*4       LATTIC                                          
!      CHARACTER*10      NAME                         
      CHARACTER*10      KNAME
!      CHARACTER*10      namen(ndif)                          
      CHARACTER*11      STATUS,FORM                                     
      CHARACTER*79      TITLE                                           
      CHARACTER*80      FNAME,fname1                                           
!                                                                       
!-----------------------------------------------------------------------
!      real*8  POS(3,NDIF),V(NATO)   
!      integer IATNR(NATO),MULT(NATO)                          
      real*8  A(3),PIA(3),VI   
!                                                                       
!                                                                       
!      dimension zz(nato), icnt(ndif,mshell,48),shdist(ndif,mshell)
!      dimension zzorig(nato)
!      dimension shellz(ndif,mshell,48),iz(ndif,mshell),ishell(ndif*nato)
!      dimension ityp(ndif*nato),imult(ndif*nato),POSN(3,NDIF)
      DIMENSION DIF(3),XX(3),PP(3),P(3),DISTS(NNN),NR(NNN),PNN(3,NNN)   
      DIMENSION NNAT(NNN),help(3)                                            
!      DIMENSION         R0(NATO),R0N(NDIF),DH(NATO),JRJ(NATO),JRJN(NDIF)
!      dimension         rmtt(NATO),rmttn(NDIF),zzn(NDIF),zzo(ndif) 
!      COMMON /ROTMAT/   ROTLOC(3,3,NATO),ROTIJ(3,3,NDIF),TAUIJ(3,NDIF)  

      LOGICAL, allocatable :: ovlap(:),used(:)
      real*8,  allocatable :: rmt(:),zz(:),shdist(:,:)
      real*8,  allocatable :: zzorig(:),shellz(:,:,:)
      real*8,  allocatable :: V(:),POSN(:,:)
      real*8,  allocatable :: R0(:),R0N(:),DH(:)  
      real*8,  allocatable :: rmtt(:),rmttn(:),zzn(:),zzo(:) 
      real*8,  allocatable :: ROTLOC(:,:,:),ROTIJ(:,:,:),TAUIJ(:,:)  
      integer,  allocatable :: JRJ(:),JRJN(:)
      integer,  allocatable :: IATNR(:),MULT(:),icnt(:,:,:) 
      integer,  allocatable :: iz(:,:),ishell(:),ityp(:),imult(:)
      CHARACTER*10, allocatable :: namen(:) 

! MAGINT insertion
      integer NNAT_EQ(NNN)
! END MAGINT insertion


!      COMMON /GENER  /  BR2(3,3)                                        
!      COMMON /CHAR/     LATTIC,NAME(NATO)                               
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      logical there
!     This should be how much larger (in 1D) the DFT distances are relative to true
      DFTERR = 1.01D0        !This is about right for PBE -- could look in case.in0?
!
!     See if the user has calibrated the lattice parameter scaling difference
      inquire(file='.latcalib',exist=there)
      if(there)then
        open(unit=99,file='.latcalib')
!       This should be how much larger (in 1D) the DFT distances are relative to true
        read(99,*)DFTERR
        close(unit=99)
      endif
      SHOWALL=.TRUE.
      write(*,*) 'specify nn-bondlength factor: (usually=2) [and optionally dlimit, dstmax (about 1.d-5, 20)]'
      dlimit=0.d0
      dstmax=-1
      read(*,'(a)') fname
      read(fname,*,end=12345) dfac,dlimit,dd
      dstmax=dd
12345 continue 
      if(dfac .lt. 0.)then
        SHOWALL=.FALSE.
        dfac=abs(dfac)
      endif
      if(dlimit.lt.1.d-7) dlimit=1.d-5
!      write(*,*) ' specify dlimit for equivalency (default=1.d-4)'
!      read(*,*,end=1234) dlimit
!      write(*,*) ' specify number of shells for equivalency '
      ishellmax=99
!      read(*,*,end=1234) ishellmax
 1234 continue
      iarg=iargc()
      call getarg(iarg,fname)
!      call getarg(2,fname)
!      if(fname.eq.'      ') call getarg(1,fname)
      OPEN(1,FILE=fname,STATUS='OLD',ERR=8000)
!      OPEN(1,FILE='//zeus/usr/lapw/def/nn.def',STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL
      OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM, &
      ERR=8002)
      if(iunit.eq.20) fname1=trim(fname)
      GOTO 8003
 8000 WRITE(*,*) ' ERROR IN OPENING NN.DEF !!!!'
      STOP 'NN.DEF'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS, &
      '  FORM:',FORM
      STOP 'OPEN FAILED'
 8001 CONTINUE
!       CALL XUFLOW(0)                                                  
      RTB2=SQRT(3.)/2.                                                  
!                                                                       
!                                                                       
!.....START READING FILE STRUCT, writing new struct (21)                   
      READ(20,1510) TITLE                                               
      write(66,1510) TITLE                                               
      write(21,1510) TITLE                                               
      READ(20,1511) LATTIC,NAT,title                                          
      write(66,1511) LATTIC,NAT,title                                          
!    allocate nato-arrays
      allocate (  name(nat) )
      allocate ( ovlap(nat) )
      allocate ( rmt(nat),zz(nat) )
      allocate ( zzorig(nat),V(nat),R0(nat),DH(nat),rmtt(nat) )
      allocate ( ROTLOC(3,3,nat),JRJ(nat),IATNR(nat),MULT(nat) ) 
! allocate pos to max value (48*nat), reduce later
      allocate ( POS(3,48*nat*2) )
!      IF(NAT.GT.NATO) STOP 'NATO TOO SMALL'                             
!     READ IN LATTICE CONSTANTS                                         
      read(20,17) A(1),A(2),A(3),alpha,beta,gamma
      if(alpha.eq.0.d0) alpha=90.d0                                       
      if(beta .eq.0.d0) beta =90.d0                                       
      if(gamma.eq.0.d0) gamma=90.d0                                       
      write(66,17) A(1),A(2),A(3),alpha,beta,gamma
 17   FORMAT(6F10.6)                                                    
!      dstmax=min(40.d0,3.99*min(a(1),a(2),a(3)))     
!      if(dstmax .lt. 0)dstmax=max(20.d0,0.61d0*max(a(1),a(2),a(3)))  
      if(dstmax .lt. 0)dstmax=max(20.d0,1.11d0*max(a(1),a(2),a(3)))  
      if(nat .eq. 1) then
        dstmax=max(20.d0,a(1),a(2),a(3))*1.1
      else
!        dstmax=min(40.d0,dstmax)
      endif
      print*,'DSTMAX:',dstmax                  
      iix=5
      iiy=5
      iiz=5
      if(lattic(1:1).eq.'P'.or.lattic(1:1).eq.'H') then
        iix=max(1,nint(dstmax/a(1)))+1
        iiy=max(1,nint(dstmax/a(2)))+1
        iiz=max(1,nint(dstmax/a(3)))+1
      endif
177   continue
      if(nat .gt. 1)then
       if( (iix-1)*a(1) .gt. dstmax*2.d0)then
        iix=max(1,iix-1)
        goto 177
       else if( (iiy-1)*a(2) .gt. dstmax*2.d0)then
        iiy=max(1,iiy-1)
        goto 177
       else if( (iiz-1)*a(3) .gt. dstmax*2.d0)then
        iiz=max(1,iiz-1)
 !        goto 177
       endif
      endif
      if(lattic(1:3).eq.'CXY') then     !fix for orthorombic CXY with very different a,b
        iix=max(iix,iiy)
        iiy=iix
      endif
      print*,'iix,iiy,iiz',iix,iiy,iiz,a(1)*iix,a(2)*iiy,a(3)*iiz

!     INDEX COUNTS ALL ATOMS IN THE UNIT CELL, JATOM COUNTS ONLY THE    
!     NONEQUIVALENT ATOMS                                               
      INDEX=0                                                           
      DO 50 JATOM = 1,NAT                                               
         INDEX=INDEX+1                                                  
!         IF(INDEX.GT.NDIF) STOP 'NDIF TOO SMALL'                        
         READ(20,1012) IATNR(JATOM),( POS(J,INDEX),J=1,3 ),MULT(JATOM)  
         write(66,1012) IATNR(JATOM),( POS(J,INDEX),J=1,3 ),MULT(JATOM)  
         IF (MULT(JATOM).EQ.0) THEN                                     
            write(66,1020) JATOM,INDEX,MULT(JATOM)                       
            STOP ' NNN: MULT EQ 0'                                      
         ENDIF                                                          
            DO 55 M=1,MULT(JATOM)-1                                     
               INDEX=INDEX+1                                            
!               IF(INDEX.GT.NDIF) STOP 'NDIF TOO SMALL'                  
               READ(20,1011) IATNR(JATOM),( POS(J,INDEX),J=1,3)         
               write(66,1011) IATNR(JATOM),( POS(J,INDEX),J=1,3)         
 55         CONTINUE                                                    
         READ(20,1050) NAME(JATOM),JRJ(JATOM),R0(JATOM),RMTT(jatom), &
                       zz(jatom)  
          zzorig(jatom)=zz(jatom)      
          if(name(jatom)(3:3).ne.' ') then
            write(66,*) 'NAMED ATOM: ',name(jatom), 'Z changed to ', &
                         'IATNR+999 to determine equivalency'
            write(*,*) 'NAMED ATOM: ',name(jatom), 'Z changed to ', &
                         'IATNR+999 to determine equivalency'
            zz(jatom)=999+jatom
          endif
         write(66,1049) NAME(JATOM),JRJ(JATOM),R0(JATOM),RMTT(jatom), &
                   zzorig(jatom)        
         if((jrj(jatom)/2)*2.eq.jrj(jatom)) then
          write(*,*) 'WARNING: JRJ of atom',jatom,' is even:',jrj(jatom)
          write(*,*) 'CHANGE it to ODD number !!!!'
          write(66,*) 'WARNING: JRJ of atom',jatom,' is even:', &
                                             jrj(jatom)
          write(66,*) 'CHANGE it to ODD number !!!!'
         endif
         DH(JATOM)=LOG(RMTT(jatom)/R0(JATOM)) / (JRJ(JATOM)-1)                 
         RMT(JATOM)=R0(JATOM)*EXP( DH(JATOM)*(JRJ(JATOM)-1) )           
      READ(20,1051) ((ROTLOC(I1,I2,JATOM),I1=1,3),I2=1,3)               
      write(66,1051) ((ROTLOC(I1,I2,JATOM),I1=1,3),I2=1,3)               
 50   CONTINUE                                                          
!    allocate ndif
      call reduce_alloc_r8_2d(3,48*nat,3,index)
      allocate ( shdist(index,mshell),shellz(index,mshell,48) )
      allocate ( namen(index),POSN(3,index),R0N(index),rmttn(index) )
      allocate ( zzn(index),zzo(index),ROTIJ(3,3,index),TAUIJ(3,index) )
      allocate ( JRJN(index),icnt(index,mshell,48),iz(index,mshell) )   
!    allocate nato*ndif
      allocate ( used(nat*index),ishell(nat*index),ityp(nat*index) )
      allocate ( imult(nat*index) )
!                                                                       
!                                                                       
!.....SET UP LATTICE, AND IF REQUIRED ROTATION MATRICES                 
      CALL DIRLAT (nat,ortho,alpha,beta,gamma)                                
!                                           
      write(66,*) ' '
      write(66,*) 'Bond-Valence Sums are calculated for current lattice parameters'
      write(66,*) 'and rescaled ones by 1 %. (You can put scaling into  .latcalib)' 
      pi=4.d0*atan(1.d0)
      cosgam=cos(gamma/180.d0*pi)
      singam=sin(gamma/180.d0*pi)           
      INDEX=0                                                           
      DO 200 JATOM=1,NAT                                                
        do 149 i=1,nat
 149       ovlap(i)=.true.
      DO 190 M=1,MULT(JATOM)                                            
      INDEX=INDEX+1                                                     
      DO 150 J=1,3                                                      
  150 XX(J)=POS(J,INDEX)                                                
      if(SHOWALL.or.(M.EQ.1))write(66,1)JATOM,M,NAME(JATOM),XX(1),XX(2),XX(3)                   
    1 FORMAT(/,' ATOM:',I3,2X,'EQUIV.',I3,2X,A10,' AT',3F10.5)          
      NC=0             
          DO 180 I1=-iix,iix                                     
          DO 180 I2=-iiy,iiy                                                
          DO 180 I3=-iiz,iiz                                                
          IF(ortho) THEN                                   
            P(1)=I1*BR2(1,1)+I2*BR2(2,1)+I3*BR2(3,1)                    
            P(2)=I1*BR2(1,2)+I2*BR2(2,2)+I3*BR2(3,2)                    
            P(3)=I1*BR2(1,3)+I2*BR2(2,3)+I3*BR2(3,3)                    
            ELSE                                                        
            P(1)=I1                                                     
            P(2)=I2                                                     
            P(3)=I3                                                     
            IF(LATTIC(1:3).eq.'CXZ') THEN
              P(1)=I1*0.5d0+i3*0.5d0
              P(2)=I2
              P(3)=-I1*0.5d0+i3*0.5d0
            END IF
          ENDIF                                                         
          K=0                                                           
      DO 120 JAT=1,NAT                                                  
      DO 110 MM=1,MULT(JAT)                                             
      K=K+1                                                             
      DIST=0.                                                           
      DO 100 L=1,3                                                      
      PP(L)=POS(L,K)+P(L)
  100 DIF(L)=XX(L)-PP(L)                                                
      IF (.not.ortho) THEN                                       
        help(1)=dif(1)  
        help(2)=dif(2)  
        help(3)=dif(3)  
      if(lattic(1:1).eq.'R') then
        dif(1)=help(1)*BR2(1,1)+help(2)*BR2(2,1)+help(3)*BR2(3,1)             
        dif(2)=help(1)*BR2(1,2)+help(2)*BR2(2,2)+help(3)*BR2(3,2)             
        dif(3)=help(1)*BR2(1,3)+help(2)*BR2(2,3)+help(3)*BR2(3,3)           
      elseif(lattic(1:3).eq.'CXZ') then
        dif(1)=help(1)*singam            
        dif(2)=(help(1)*cosgam*a(1)+help(2)*a(2))/a(2)             
        dif(3)=help(3)           
      else
        dif(1)=(help(1)*BR2(1,1)*a(1)+help(2)*BR2(2,1)*a(2)+ &
                help(3)*BR2(3,1)*a(3))/a(1)             
        dif(2)=(help(1)*BR2(1,2)*a(1)+help(2)*BR2(2,2)*a(2)+ &
                help(3)*BR2(3,2)*a(3))/a(2)             
        dif(3)=(help(1)*BR2(1,3)*a(1)+help(2)*BR2(2,3)*a(2)+ &
                help(3)*BR2(3,3)*a(3))/a(3)           
      endif
      ENDIF                                                             
      DO 103 L=1,3                                                      
  103 DIST=DIST+DIF(L)*DIF(L)*A(L)*A(L)                                 
      DIST=SQRT(DIST)                                                   
      IF(DIST.GT.dstmax) GO TO 110                                        
      IF(DIST.LT..001) GO TO 110                                        
      NC=NC+1         
      if(nc.gt.nnn) stop ' nnn too small'                                     
      DISTS(NC)=DIST                                                    
      NNAT(NC)=JAT                                                      
! MAGINT insertion
      NNAT_EQ(NC)=MM
! END MAGINT insertion
      DO 105 L=1,3                                                      
  105 PNN(L,NC)=PP(L)                                                   
!     write(66,2)JAT,NAME(JAT),PP(1),PP(2),PP(3),DIST                    
!   2 FORMAT(' TO ATOM:',I2,2X,A10,' AT',3F8.4,                         
!    * ' IS ',F10.5,' A.U.')                                            
  110 CONTINUE                                                          
  120 CONTINUE                                                          
  180 CONTINUE                                                          
      CALL ORD2(DISTS,NR,NC)                                            
      N1=1                                                              
      N2=NR(N1)                                                         
      N3=NNAT(N2)             
      SUMRAD=RMT(JATOM)+RMT(N3)                                         
      IF(M.EQ.1) THEN                                                   
      IF(SUMRAD.LE.DISTS(N1))THEN                                       
      WRITE(*,7)JATOM,NAME(JATOM),N3,NAME(N3)                           
    7 FORMAT(/,3X,2(' ATOM',I3,2X,A10))                                 
      WRITE(*,5)JATOM,RMT(JATOM),N3,RMT(N3),SUMRAD,DISTS(1)             
      write(66,5)JATOM,RMT(JATOM),N3,RMT(N3),SUMRAD,DISTS(1)             
    5 FORMAT(' RMT(',I3,')=',F7.5,' AND RMT(',I3,')=',F7.5,/,            &
       ' SUMS TO',F8.5,2X,'LT.  NN-DIST=',F8.5)                         
      ELSE                                                              
      ovlap(n3)=.false.                        
      write(66,4) JATOM,RMT(JATOM),N3,RMT(N3), SUMRAD, DISTS(N1)         
    4 FORMAT(/,'   ERROR !!!!!!!!!!!!!!!',                               &
      /,' RMT(',I3,')=',F7.5,' AND RMT(',I3,')=',F7.5,/,                 &
       ' SUMS TO',F8.5,' GT NNN-DIST=',F8.5)                            
      WRITE(*,4) JATOM,RMT(JATOM),N3,RMT(N3), SUMRAD, DISTS(N1)         
      ENDIF                                                             
      ENDIF                                                             
!.....determination of "equal" atoms
!
      olddis=0.d0
      ishell(index)=0
!
      bva=0.0
      bvaerr=0.0
      DO 185 N1=1,NC                                                    
      N2=NR(N1)                                                         
      N3=NNAT(N2)                                                       
!
!     Include BVA
      i1=nint(zzorig(jatom))
      i2=nint(zzorig(N3))
      call bvan(i1,i2,scale,dbva)
!      write(88,88)i1,i2,scale,dbva,dists(n1)*0.529177,bva
88    format(2i3,4F10.5)
      if(scale .gt. 0)then
                val=exp( (dbva - dists(n1)*0.529177) /scale)
                bva=bva+val
                valerr=exp( (dbva - dists(n1)*0.529177/DFTERR) /scale)
                bvaerr=bvaerr+valerr
      endif 
!
      SUMRAD=RMT(JATOM)+RMT(N3)                                        
      if(dists(n1).lt.dfac*dists(1)) then
! MAGINT insertion
       if(SHOWALL.or.(M.EQ.1))write(66,3)N3,NAME(N3),(PNN(L,N2),L=1,3),&
                   DISTS(N1),DISTS(N1)*0.529177,NNAT_EQ(N2)
!                   DISTS(N1),DISTS(N1)*0.529177,NNAT_EQ(N2)
! END MAGINT insertion
       IF(ovlap(n3).and.SUMRAD.GE.DISTS(N1)) THEN
        ovlap(n3)=.false.                        
        write(66,4) JATOM,RMT(JATOM),N3,RMT(N3),SUMRAD,DISTS(N1)
        WRITE(*,4) JATOM,RMT(JATOM),N3,RMT(N3), SUMRAD, DISTS(N1)
       end if         
      end if         
!
      if((dists(n1)-olddis).gt.dlimit) then
!.....new shell
        ishell(index)=ishell(index)+1
        iz(index,ishell(index))=1
        olddis=dists(n1)  
        if(ishell(index).ge.mshell) goto 187
        icnt(index,ishell(index),iz(index,ishell(index)))=1
        shellz(index,ishell(index),iz(index,ishell(index)))=zz(n3)
        shdist(index,ishell(index))=olddis
      else
!.....old shell
        do i=1,iz(index,ishell(index))
        if(zz(n3).eq.shellz(index,ishell(index),i)) then
          icnt(index,ishell(index),i)=icnt(index,ishell(index),i)+1
          goto 186
        endif
        enddo
        iz(index,ishell(index))=iz(index,ishell(index))+1
        shellz(index,ishell(index),iz(index,ishell(index)))=zz(n3)
        icnt(index,ishell(index),iz(index,ishell(index)))=1          
      endif
 186  continue
!
  185 CONTINUE
  187  CONTINUE
     if( (SHOWALL.or.(M.EQ.1)).and.(bva .gt. 0.001D0) .and. (M.eq.1) )then
        write(66,666)JATOM,M,NAME(JATOM),bva,bvaerr
666     format('Atom ',i3,' equiv ',i2,' ',a,' Bond-Valence Sum ',2F8.2)
!       write(*,*)'Bond-valence sum ',bva
     endif

!
!.....limit shells to some maximum
      write(67,*) '---------------- ATOM',index
      do i=1,ishell(index)
      if(shdist(index,i).gt.dstmax-dstmax/3.d0) then
         ishell(index)=i
         goto 190
      endif
!c
      write(67,*) ' SHELL:',i,shdist(index,i)
        do j=1,iz(index,i)
        write(67,*) icnt(index,i,j),shellz(index,i,j)
        enddo
!c
      enddo
!
  190 CONTINUE                                                          
  200 CONTINUE                                                          
!                      
      write(66,552)
 552  format(//,'SHELL STRUCTURE for all ATOMS:',/, &
      'ATOM  | DISTANCE   #of NEIGHBORS   Z |')
      INDEX=0                                                           
      DO  JATOM=1,NAT                
      DO  M=1,MULT(JATOM)                                            
        INDEX=INDEX+1
        zzo(index)=zz(jatom)                                                  
        write(66,553)index,((shdist(index,i),icnt(index,i,j), &
                     shellz(index,i,j),j=1,iz(index,i)),i=1,4)
 553    format(i3,1x,8('|',f6.3,i3,f6.1))
      enddo
      enddo
      write(66,*)
!
      inat=0
      do i=1,index
      ityp(i)=i
      imult(i)=0
      used(i)=.false.  
      enddo      
      write(66,554)
 554  format(/,'LISTING of INDEX of all EQUIVALENT ATOMS:')
      ityp1=0    
      ij=0          
      i=0                           
      do 500 i0=1,nat
      do 500 i00=1,mult(i0)
      i=i+1
        if(used(i)) goto 500
        write(66,*)
        do 501 j=i,index
          if(zz(i0).ne.zzo(j)) goto 501
!     compare atom with index i with shells of other atoms (Fuhr)
          if (ishell(i).ne.ishell(j)) then
             write(66,559) i,j,ishell(i),ishell(j)
 559         format(' atom:',i4,' and ATOM:',i4,' differ in number of shells',i4,'ne',i4)
             goto 501
          endif
!     compare atom with index i with all other atoms
          do 510 i1=1,ishell(i)-1
            if(abs(shdist(i,i1)-shdist(j,i1)).gt.dlimit.and.shdist(i,i1).lt.2.d000*max(a(1),a(2),a(3))) then
!            if(abs(shdist(i,i1)-shdist(j,i1)).gt.dlimit) then
              write(66,550) i,j,i1,shdist(i,i1),shdist(j,i1)
 550          format(' atom:',i4,' and ATOM:',i4,' differ in shell',i3,2f10.5)
              goto 501
            endif
            do 511 i2=1,iz(i,i1)
             do j2=1,iz(j,i1)
             if((icnt(i,i1,i2).eq.icnt(j,i1,j2)).and. &
             (shellz(i,i1,i2).eq.shellz(j,i1,j2))) goto 511
             if(shdist(i,i1).gt.2.d0000*max(a(1),a(2),a(3))) goto 511
             enddo
              write(66,551) i,j,i1,icnt(i,i1,i2),icnt(j,i1,j2-1),shellz(i,i1,i2),shellz(j,i1,j2-1),shdist(i,i1)
 551          format(' atom:',i4,' and ATOM:',i4,' differ in shell',i3,2i4,2f7.1,f10.5)
             goto 501
 511        continue
 510      continue
          write(66,555) i,j
 555      format(' ATOM:',i4,' and ATOM:',i4,' are equivalent')
          if(i.eq.j) then
            ityp1=ityp1+1
            namen(ityp1)=name(i0)
            JRJN(ityp1)=jrj(i0)
            R0N(ityp1)=r0(i0)
            RMTTN(ityp1)=rmtt(i0)
            zzn(ityp1)=zzorig(i0)
          endif
          ityp(j)=ityp1
          if(inat.lt.ityp(j)) inat=ityp(j)
          imult(ityp1)=imult(ityp1)+1
          used(j)=.true.
          ij=ij+1
          posn(1,ij)=pos(1,j)
          posn(2,ij)=pos(2,j)
          posn(3,ij)=pos(3,j)
 501    continue
 500  continue
!
      write(66,*)      
      error=.false.
      INDEX=0                                                           
      DO  JATOM=1,NAT                
      write(66,556) jatom,mult(jatom),imult(jatom)  
 556  format(/,' ATOM KIND:',i4,'  OLD and NEW MULTIPLICITY:  ',2i4)
      if(mult(jatom).ne.imult(jatom)) then
           error=.true.
           write(66,*)'WARNING: MULT not equal. The new multiplicity is', &
           ' different from the old one'  
           write(6,*)'WARNING: Mult not equal. PLEASE CHECK outputnn-file'
      end if 
      DO  M=1,MULT(JATOM)                                            
      INDEX=INDEX+1                                                     
 557  format(5x,'ATOM INDEX:',i4,'  OLD and NEW ATOM KIND:',2i4)
      if(jatom.ne.ityp(index)) then
           error=.true.
           write(66,557) index,jatom,ityp(index)
           write(66,*) 'WARNING: ITYP not equal. The new type is', &
           ' different from the old one'  
           write(6,*)'WARNING: ityp not equal. PLEASE CHECK outputnn-file'
      endif
      enddo
      enddo
!
      if(.not.error) stop 'NN ENDS' 
      write(66,*)
      write(66,*) 'NEW LIST of EQUIVALENT POSITIONS written to', &
                  ' case.struct_nn'   
!.....write new struct file
!
      write(21,1511) LATTIC,inat,title                                  
      write(21,17) A(1),A(2),A(3),alpha,beta,gamma
      index=0
      do jatom=1,inat
      INDEX=INDEX+1                                                     
         write(21,1012) -JATOM,( POSN(J,INDEX),J=1,3 ),iMULT(JATOM),8  
         write(66,*)
         write(66,1011) -JATOM,( POSN(J,INDEX),J=1,3 ),NAMEN(JATOM), &
                        zzn(jatom) 
            DO M=1,iMULT(JATOM)-1                                     
               INDEX=INDEX+1                                            
               write(21,1011) -JATOM,( POSN(J,INDEX),J=1,3)         
               write(66,1011) -JATOM,( POSN(J,INDEX),J=1,3)         
            ENDDO                                                    
         if(jatom.lt.10) then
         write(21,1149) NAMEN(JATOM)(1:2),jatom,JRJN(JATOM),R0N(JATOM),RMTTN(jatom),zzn(jatom)       
         else if(jatom.lt.100) then
         write(21,1148) NAMEN(JATOM)(1:2),jatom,JRJN(JATOM),R0N(JATOM),RMTTN(jatom),zzn(jatom)       
         else
         write(21,1147) NAMEN(JATOM)(1:2),jatom,JRJN(JATOM),R0N(JATOM),RMTTN(jatom),zzn(jatom)       
         endif
         write(21,1512) 1.0,0.0,0.0
         write(21,1512) 0.0,1.0,0.0
         write(21,1512) 0.0,0.0,1.0
      enddo
         write(21,*) 0
!                 
      print*,' '
      print*," NN created a new ",trim(fname1),"_nn file"   
                                                     
      STOP 'NN created a new CASE.STRUCT_NN FILE'                          
!
! MAGINT insertion
!    3 FORMAT(' ATOM:',I3,2X,A10,'AT',3F10.6, &
!       ' IS',F9.5,' A.U.',f10.5,' ANG',i6)  
    3 FORMAT(' ATOM:',I3,2X,A10,'AT',3F8.4,                              &
       ' IS',F9.5,' A.U.',f10.5,' ANG',i6)               
! END MAGINT insertion
 1011 FORMAT(4X,I4,4X,F10.8,3X,F10.8,3X,F10.8,5x,a10,f5.1)         
 1012 FORMAT(4X,I4,4X,F10.8,3X,F10.8,3X,F10.8,/,15X,I2,17X,I2)     
 1020 FORMAT(///,3X,'ERROR IN EBND  : MULT(JATOM)=0 ...',                &
             /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)           
 1049 FORMAT(A10,' NPT=',I5,'  R0=',F10.8,' RMT=',F10.4,'   Z:',f5.1)      
 1149 FORMAT(A2,i1,7x,' NPT=',I5,'  R0=',F10.8,' RMT=',F10.4,'   Z:',f5.1)  
 1148 FORMAT(A2,i2,6x,' NPT=',I5,'  R0=',F10.8,' RMT=',F10.4,'   Z:',f5.1)  
 1147 FORMAT(A2,i3,5x,' NPT=',I5,'  R0=',F10.8,' RMT=',F10.4,'   Z:',f5.1)  
 1050 FORMAT(A10,5X,I5,5X,F10.8,5X,F10.4,5x,f5.1)                           
 1051 FORMAT(20X,3F10.7)                                                
 1510 FORMAT(A79)                                                       
 1511 FORMAT(A4,23X,I3,/,13X,A4)                                        
 1512 FORMAT(20x,3f10.7)                                                       
      END                                                               
