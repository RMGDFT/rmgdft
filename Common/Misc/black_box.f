c************************** SVN Revision Information **************************
c **    $Id: black_box.f 1004 2008-08-09 18:37:21Z froze $    **
c******************************************************************************
 
*-----------------------------------------------------------------------
      subroutine symmetry(ibrav,s,nsym,irg,irt,ftau,nat,tau,ityp,nks,
     &                    xk,wk,celldm,nr1,nr2,nr3, wflag)
*-----------------------------------------------------------------------

* the subroutine generates the crystal cell in real space given the bravais
* lattice index ibrav. Look for the right symmetry operations, defines
* the true point group and check the consistency of the special points
* used. From the pwscf library for the real space multigrid program.
*(Last revised by MBN, apr. 1996).
* In the following is the meaning of the ibrav parameter and the physical
* quantities requested:
*-----------------------------------------------------------------------
*
*     ibrav is the structure index:
*
*     ibrav        structure
*       1          cubic P (sc)
*       2          cubic F (fcc)
*       3          cubic I (bcc)
*       4          Hexagonal and Trigonal P
*       5          Trigonal R
*       6          Tetragonal P (st)
*      61          Tetragonal (001) SL
*      62          Orthorhombic (110) SL
*       7          Tetragonal I (bct)
*       8          Orthorhombic P
*      12          Monoclinic P
*      14          Triclinic P
*
*     The remaining bravais lattices are not programmed - use ibrav=8
*     for body-centered orthorhombic, etc.
*
*   indpg is the point-group index, defined as follows:
*
* indpg   group     indpg    group     indpg   group      indpg  group
*  1    1  (C1 )      9    3m  (C3v)    17   4/mmm(D4h)    25   222(D2 )
*  2   <1> (Ci )     10   <3>m (D3d)    18   6    (C6 )    26   mm2(C2v)
*  3    2  (C2 )     11    4   (C4 )    19   <6>  (C3h)    27   mmm(D2h)
*  4    m  (C1h)     12    <4> (S4 )    20   6/m  (C6h)    28   23 (T  )
*  5    2/m(C2h)     13    4/m (C4h)    21   622  (D6 )    29   m3 (Th )
*  6    3  (C3 )     14    422 (D4 )    22   6mm  (C6v)    30   432(O  )
*  7    <3>(C3i)     15    4mm (C4v)    23   <6>m2(D3h)    31 <4>3m(Td )
*  8    32 (D3 )     16   <4>2m(D2d)    24   6/mmm(D6h)    32   m3m(Oh )
*
*-----------------------------------------------------------------------
*
*  point group bravais lattice    ibrav  celldm(2)-celldm(6)
*.......................................................................
*     432,<4>3m,m3m     sc          1     not used
*.......................................................................
*     23,m3             sc          1         "
*.......................................................................
*     432,<4>3m,m3m    fcc          2         "
*.......................................................................
*     23,m3            fcc          2         "
*.......................................................................
*     432,<4>3m,m3m    bcc          3         "
*     23,m3             sc          1         "
*.......................................................................
*     432,<4>3m,m3m    fcc          2         "
*.......................................................................
*     23,m3            fcc          2         "
*.......................................................................
*     432,<4>3m,m3m    bcc          3         "
*.......................................................................
*     23,m3            bcc          3         "
*.......................................................................
*     622,6mm,
*     <6>m2,6/mmm      hex(p)       4      celldm(3)=c/a
*.......................................................................
*     6,<6>,6/m,       hex(p)
*     32,3m,<3>m      trig(p)       4         "
*.......................................................................
*     3,<3>           trig(p)       4         "
*.......................................................................
*     32,3m,<3>m      trig(r)       5     celldm(4)=cos(aalpha)
*.......................................................................
*     3,<3>           trig(r)       5         "
*.......................................................................
*     422,4mm,
*     <4>2m,4/mmm     tetr(p)       6      celldm(3)=c/a
*.......................................................................
*     4,<4>,4/m       tetr(p)       6         "
*.......................................................................
*     422,4mm,
*     <4>2m,4/mmm     tetr(i)       7         "
*.......................................................................
*     4,<4>,4/m       tetr(i)       7         "
*.......................................................................
*     222,mm2,mmm     orth(p)       8     above + celldm(2)=b/a
*.......................................................................
*     2,m,2/m         mcln(p)      12     above + celldm(4)=cos(ab)
*.......................................................................
*     1,<1>           tcln(p)      14     celldm(2)= b/a
*                                         celldm(3)= c/a
*                                         celldm(4)= cos(bc)
*                                         celldm(5)= cos(ac)
*                                         celldm(6)= cos(ab)
*-----------------------------------------------------------------------
*
*   The special axis is the z-axis, one basal-plane vector is along x,
*   and the other basal-plane vector is at angle beta for monoclinic
*   (beta is not actually used), at 120 degrees for trigonal and hexagonal(p)
*   groups, and at 90 degrees for remaining groups, excepted fcc, bcc,
*   tetragonal(i), for which the crystallographic vectors are as follows:
*
*   fcc bravais lattice.
*   ====================
*
*   a1=(a/2)(-1,0,1), a2=(a/2)(0,1,1), a3=(a/2)(-1,1,0).
*
*   bcc bravais lattice.
*   ====================
*
*   a1=(a/2)(1,1,1), a2=(a/2)(-1,1,1), a3=(a/2)(-1,-1,1).
*
*   tetragonal (i) bravais lattices.
*   ================================
*   a1=(a/2,a/2,c/2), a2=(a/2,-a/2,c/2), a3=(-a/2,-a/2,c/2).
*
*   trigonal(r) groups.
*   ===================
*
*   for these groups, the z-axis is chosen as the 3-fold axis, but the
*   crystallographic vectors form a three-fold star around the z-axis,
*   and the primitive cell is a simple rhombohedron. if c is the cosine
*   of the angle between any pair of crystallographic vectors, and if
*   tx=sqrt((1-c)/2), ty=sqrt((1-c)/6), tz=sqrt((1+2c)/3), the crystal-
*   lographic vectors are:
*
*         a1=a(0,2ty,tz),  a2=a(tx,-ty,tz),  a3=a(-tx,-ty,tz).
*
*-----------------------------------------------------------------------

      implicit none
 
      integer npk, nrot, nrotot,  nat, nks, nsym, nr1, nr2, nr3
      integer natp, nksp, wflag
      parameter(natp=600,nksp=2000)

      integer s(3,3,48), irt(48,natp), irg(48),krg(48),ityp(natp),
     & ftau(3,48), nks0, na,nb, i,j,k,ncos,ik,jk,ic,jc, is,ir,jr,
     & n,isym,irot,n6,nf,m,ibrav,iwork(6)
      real*8 at(3,3),bg(3,3),tau(3,natp),rtau(3,48,natp),xau(3,natp), 
     &       rau(3,natp),xks(3,48),xkg(3),w(48),work(6),sm1(3,3),
     &       xk(3,nksp),wk(nksp),ft(3),ft1, ft2, ft3, one, sw, ddsum,
     &       hxgrid,hygrid,hzgrid
      integer invs(3,3,48)
      real*8 celldm(6),omega
      logical invsym
      external checksym, error, ddsum, coset, matinv

      if(nat .ge. natp) then
        print *,'symmetry  NAT  .GT.   NATP', nat, natp
        stop
      endif

      if(nks .ge. nksp) then
        print *,'symmetry  NKS  .GT.   NKSP', nks, nksp
        stop
      endif


c generate direct space vectors 
      call latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega,1)
c reciprocal space vectors
      call recips_f(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))

      if(wflag .eq. 0)then
        write(6,'(5x,''crystal axis:  ''/
     &          3(''('',3f7.4,'')  '') )') ((at(i,j), i=1,3), j=1,3)
        write(6,'(5x,''reciprocal axis:  ''/
     &          3(''('',3f7.4,'')  '') )') ((bg(i,j), i=1,3), j=1,3)
        print*,' cell volume ', omega
      endif
      if (ibrav.eq.4 .or. ibrav.eq.5) then
         call hexsym(at,s,nrot)
      else if(ibrav.gt.1 .or. ibrav.le.14) then
         call cubicsym(at,s,nrot)
      else
         call error('setup','wrong ibrav',1)
      end if

      call sgama (nrot,nat,s,at,bg,tau,ityp,nsym,nr1,nr2,nr3,
     &            irg,irt,ftau,rtau,nks,nks,xk,wk,xau,rau,
     &            invsym, wflag)

      if(wflag .eq. 0)then
        write(6,'(5x,''symmetry operations (in crystal axis):  '')')
      endif
      do n6 = 0,nsym/6
         nf = min(nsym-6*n6,6)
         if(wflag .eq. 0)then
           write(6,'( )')
           do m=1,3
              write (6,'(6(3i3,2x))')  ((s(m,j,n6*6+n), j=1,3), n=1,nf)
           end do
         endif
      end do
      do 10 isym = 1,nsym
         irot = irg(isym)
         if (ftau(1,irot).ne.0 .or. ftau(2,irot).ne.0 .or.
     +   ftau(3,irot).ne.0) then
           if(wflag .eq. 0)then
              write (6,
     +'('' the following symmetry operations have'',
     +  '' fractionary translations:'')')
           endif
            go to 20
         end if
   10 continue
   20 do 30 isym = 1,nsym
         irot = irg(isym)
         if (ftau(1,irot).ne.0 .or. ftau(2,irot).ne.0 .or.
     +   ftau(3,irot).ne.0) then
            ft1 = at(1,1)*ftau(1,irot)/nr1 + at(1,2)*ftau(2,irot)/nr2 +
     +      at(1,3)*ftau(3,irot)/nr3
            ft2 = at(2,1)*ftau(1,irot)/nr1 + at(2,2)*ftau(2,irot)/nr2 +
     +      at(2,3)*ftau(3,irot)/nr3
            ft3 = at(3,1)*ftau(1,irot)/nr1 + at(3,2)*ftau(2,irot)/nr2 +
     +      at(3,3)*ftau(3,irot)/nr3
            if(wflag .eq. 0)then
              write (6,'(5x,''op. # '',i2,5x,''tau = '',3f11.7)') irot,
     +        ft1,ft2,ft3
            endif
         end if
   30 continue

      if(wflag .eq. 0)then
        write(6,'(/i5,5x,''special points and weights:  '')') nks
        do n=1,nks
           write(6,'(5x,'' '',3f12.8,''  '',1x,f12.8)')
     &           (xk(i,n), i=1,3),wk(n)
        end do
      endif
      return
      end

*-----------------------------------------------------------------------
      subroutine sgama (nrot,nat,s,at,bg,tau,ityp,nsym,nr1,nr2,nr3,
     &                  irg,irt,ftau,rtau,npk,nks,xk,wk,xau,rau,
     &                  invsym, wflag)
*-----------------------------------------------------------------------
*
*     given a point group, this routine finds the subgroup which is
*     the point group of the crystal under consideration
*     non symmorphic groups allowed, ftau is the fractionary translation
*
c     nrot     number of symmetry operations (ouput of cubic/hexsym)
c     nat      number of atoms 
c     s        point group operations
c     at       real space lattice vectors
c     bg       reciprocal space lattice space vectors
c     tau      atomic positions
c     ityp     kind of atom
c     nsym     order of the true point group
c     nr1,2,3  real space mesh
c     irg(48)  if the i-th symm. op. belongs to the true point group of the
c              crystal => n=n+1 , irg(n) = i      
c     irt(48,nat) if the i-th symm. op. moves the nb atom on the na atom
c              and ityp(na)=ityp(nb) => irt(i,na) = nb 
c     ftau     non primitive translations
c     rtau     coordinates of the rotated atoms
c     npk      maximum number of special k points
c     nks      number of true k points given the actual symmetry
c     xk       k-point coordinates
c     wk       k-point weights
c     xau      coordinates of the atom in the crystalline basis
c     rau      coordinates of the rotated atom in the crystalline basis
c     invsym   inverse of the point group operations

      implicit none
      logical satm, latm, sym(48), invsym
      integer npk, nrot, nrotot, nat, nks, nsym, nr1, nr2, nr3
      integer wflag
      integer n, s(3,3,48), irt(48,nat), irg(48), krg(48), ityp(nat),
     & ftau(3,48), nks0, na,nb, i,j,k, ncos, ik,jk, ic,jc, is,ir,jr,
     & iwork(6)
      real*8 at(3,3), bg(3,3), tau(3,nat), rtau(3,48,nat), xau(3,nat), 
     &       rau(3,nat), xks(3,48), xkg(3), w(48), work(6), sm1(3,3),
     &       xk(3,npk), wk(npk), ft(3), ft1, ft2, ft3, one, sw, ddsum,
     &       ascal
      integer invs(3,3,48)
      external checksym, error, ddsum, coset, matinv
*
* inversion symmetry is added to the candidate symmetries
*
      do ir=1, nrot
         do j=1, 3
            do i=1, 3
               s(i,j,ir+nrot) = - s(i,j,ir)
            end do
         end do
      end do
      nrotot=2*nrot
*
* xau are the atomic coordinate in the basis of the primitive translation
* vectors (crystal basis)
*                              
      do 21 na=1, nat
         do 20 k=1,3
            xau(k,na) = bg(1,k)*tau(1,na) + bg(2,k)*tau(2,na) +
     &                  bg(3,k)*tau(3,na)
  20     continue
  21  continue
*
      do 190 ir=1, nrotot
*
* rau are the crystal coordinates of the rotated atom
*
         do 170 na=1, nat
            do 180 k=1, 3
                rau(k,na) = s(1,k,ir)*xau(1,na) +
     &                      s(2,k,ir)*xau(2,na) +
     &                      s(3,k,ir)*xau(3,na)
 180        continue
 170     continue
*
*      first attempt: no fractionary translation
*
         do 171 i=1,3
            ftau(i,ir) = 0
            ft  (i)      = 0.0
 171     continue
*
         call checksym(ir,nat,ityp,xau,rau,ft,sym,irt)
*
         if(.not.sym(ir)) then
            do 160 na=1, nat
               do 150 nb = 1, nat
                  if(ityp(nb).eq.ityp(na)) then
*
*      second attempt: check all possible fractionary translations
*
                   ft(1)=rau(1,na)-xau(1,nb)-nint(rau(1,na)-xau(1,nb))
                   ft(2)=rau(2,na)-xau(2,nb)-nint(rau(2,na)-xau(2,nb))
                   ft(3)=rau(3,na)-xau(3,nb)-nint(rau(3,na)-xau(3,nb))
                   call checksym(ir,nat,ityp,xau,rau,ft,sym,irt)
                   if(sym(ir)) then
*
* convert ft to fast-fourier transform
* coordinates for later use in symmetrizations
*
                        ft1 = ft(1)*nr1
                        ft2 = ft(2)*nr2
                        ft3 = ft(3)*nr3
                        if(abs(ft1-nint(ft1))/nr1 .gt. 1.0e-5 .or.
     &                     abs(ft2-nint(ft2))/nr2 .gt. 1.0e-5 .or.
     &                     abs(ft3-nint(ft3))/nr3 .gt. 1.0e-5) then
                           if(wflag .eq. 0)then
                           write(6,'(5x,''warning: symmetry operation'',
     &                      '' # '',i2,'' not allowed.   fractionary '',
     &                      ''translation:''/5x,3f11.7,''  in crystal'',
     &                      '' coordinates'')')       ir, ft
                           endif
                           sym(ir) = .false.
                        end if
                        ftau(1,ir) = nint(ft1)
                        ftau(2,ir) = nint(ft2)
                        ftau(3,ir) = nint(ft3)
                        go to 190
                     end if
                  end if
 150           continue
 160        continue
         end if
 190  continue   
*
      nsym = 0
      do 200 ir=1, nrotot
*
* find the inverse matrix of the symmetry operation for k vectors rotation
*
         do i = 1,3
            do j = 1,3
               sm1(i,j) = float(s(i,j,ir))
            end do
         end do
         call matinv(sm1,3,3,iwork)
         do i = 1,3
            do j = 1,3
               invs(i,j,ir) = nint(sm1(i,j))
            end do
         end do
*
         if(sym(ir)) then
            nsym = nsym +1
            krg(nsym)=ir
            irg(nsym)=ir
         end if
 200  continue
*
      call coset(nsym,nrotot,krg,s)
*
      ncos = nrotot/nsym
      nks0 = nks
      do 390 jk=1, nks0
*
* k vectors are projected along the crystal axis bg for the reciprocal lattice
*
         do 300 k=1,3
            xkg(k) = at(1,k)*xk(1,jk) + at(2,k)*xk(2,jk) +
     &               at(3,k)*xk(3,jk)
 300     continue
*
         do 350 ir=1,nrotot
            jr = krg(ir)
*
            do 310 k=1, 3
               xks(k,ir)=invs(k,1,jr)*xkg(1)+invs(k,2,jr)*xkg(2)
     +                  +invs(k,3,jr)*xkg(3)
 310        continue
 350     continue
*
         do 380 ic = 1, ncos
            ir = (ic-1)*nsym + 1
            latm = .false.
*
*   latm = .true. if the present k-vector is equivalent to some previous one
*
            do 377 jc =1, ic-1
               do 376 is = 1, nsym
*
*   satm = .true. if the present symmetry operation makes the ir and ik
*   k-vectors equivalent (nb: or equivalent to minus each other)
*
                  ik = (jc-1)*nsym + is
                  satm =       abs(  xks(1,ir)-xks(1,ik) -
     -                          nint(xks(1,ir)-xks(1,ik)) ) .lt. 1.0e-5
     &                   .and. abs(  xks(2,ir)-xks(2,ik) -  
     -                          nint(xks(2,ir)-xks(2,ik)) ) .lt. 1.0e-5
     &                   .and. abs(  xks(3,ir)-xks(3,ik) -
     -                          nint(xks(3,ir)-xks(3,ik)) ) .lt. 1.0e-5
     &             .or.        abs(  xks(1,ir)+xks(1,ik) - 
     -                          nint(xks(1,ir)+xks(1,ik)) ) .lt. 1.0e-5
     &                   .and. abs(  xks(2,ir)+xks(2,ik) -
     -                          nint(xks(2,ir)+xks(2,ik)) ) .lt. 1.0e-5
     &                   .and. abs(  xks(3,ir)+xks(3,ik) -
     -                          nint(xks(3,ir)+xks(3,ik)) ) .lt. 1.0e-5
                  latm = latm .or. satm
                  if( satm .and. w(jc).ne.0.0 ) then
                     w(jc) = w(jc) + 1.0
                     go to 378
                  end if
 376           continue
 377        continue
*
 378        continue
            if(latm) then
               w(ic) = 0.0
            else
               w(ic) = 1.0
            end if
 380     continue
*
cwgs         sw = wk(jk)/ddsum(ncos,w,1)
cwgs         wk(jk) = sw * w(1)
cwgs         do 385 ic = 2,ncos
cwgs            ir = (ic-1)*nsym + 1
cwgs            if(w(ic).ne.0.0) then
cwgs               nks = nks + 1
               if(nks.gt.npk) call error
     &         ('sgama','too many k-points',nks)
cwgs               wk(nks) = sw * w(ic)
cwgs               do 320 k=1,3
cwgs                xk(k,nks) = bg(k,1)*xks(1,ir) + bg(k,2)*xks(2,ir)
cwgs     +                      + bg(k,3)*xks(3,ir)
cwgs 320           continue
cwgs            end if
cwgs 385     continue
 390  continue
*

      one=ddsum(nks,wk,1)
      if(one.gt.0.0) then
        do ir = 1, nks
         wk(ir) = wk(ir) / one
        enddo  
      endif
      invsym=sym(nrot+1)
*
* we copy the symm. op. on the first 24 positions of s, and the last 24
* are those containing the inversion symmetry
*
      is = 0
      do ir=1, nsym
  10     is=is+1
         if(.not.sym(is)) go to 10
         if(is.ne.ir) then
            do na=1, nat
               irt(ir,na) = irt(is,na)
            end do
            do i=1, 3
               ftau(i,ir)=ftau(i,is)
               do j=1, 3
                  s(i,j,ir)=s(i,j,is)
               end do
             end do
             sym(ir)=.true.
         end if
      end do
*
      if(invsym) then
         if( ftau(1,nrot+1).ne.0 .or. ftau(2,nrot+1).ne.0 .or.
     &       ftau(3,nrot+1).ne.0 ) then
            if(wflag .eq. 0)then
              write(6,
     &          '(5x,''inversion symmetry+fractionary translation'')')
            endif
            invsym = .false.
         else
            if(wflag .eq. 0)then
              write(6,'(5x,''inversion symmetry'')')
            endif
         end if
* si sottrae l'inversione dalle operazioni di simmetria
         nsym=nsym/2
      else
         if(wflag .eq. 0)then
           write(6,'(5x,''no inversion symmetry'')')
         endif
      end if
* si ricopiano le inverse in s(i,j,ir+24)
      do ir=1, nsym
         do j=1, 3 
            do i=1, 3
               s(i,j,ir+24)=invs(i,j,irg(ir) )
             end do
         end do
         irg(ir)=ir
      end do
*
      ascal = 1.0/nsym
      do ir =1,nks
        wk(ir) = wk(ir) * ascal
      enddo
*
      return
      end

*-----------------------------------------------------------------------
      subroutine checksym(ir,nat,ityp,xau,rau,ft,sym,irt)
*-----------------------------------------------------------------------
*
      implicit none
      integer nat, irt(48,nat), ityp(*), na, nb, ir
      real*8 xau(3,nat), rau(3,nat), ft(3)
      logical sym(48), eqvect
      external eqvect
*
      do 10 na = 1, nat
         do 20 nb = 1, nat
            sym(ir) = ityp(na).eq.ityp(nb) .and.
     +            eqvect(rau(1,na),xau(1,nb),ft)
            if(sym(ir)) then
*
* the rotated atom does coincide with one of the like atoms
* keep track of which atom the rotated atom coincides with
*
               irt(ir,na)=nb
               go to 10
            end if
 20      continue
*
* the rotated atom does not coincide with any of the like atoms
* s(ir) + ft is not a symmetry operation
*
         return
 10   continue
*
* s(ir) + ft is a symmetry operation
*
      return
      end
*
      subroutine cubicsym(at,is,nrot)
* Provides symmetry operations (excepted Inversion) for all cubic and
* lower-symmetry (excepted Hexagonal and Trigonal) bravais lattices
      implicit none
      real*8 at(3,3), s(3,3,24), overlap(3,3), rat(3), rot(3,3), value
      integer is(3,3,48), iwork(3), j, m, k, n, nrot
      external matinv
      data s/ 
     +           1, 0, 0,  0, 1, 0,  0, 0, 1   ,
     +          -1, 0, 0,  0,-1, 0,  0, 0, 1   ,
     +          -1, 0, 0,  0, 1, 0,  0, 0,-1   ,
     +           1, 0, 0,  0,-1, 0,  0, 0,-1   ,
     +           0, 1, 0,  1, 0, 0,  0, 0,-1   ,
     +           0,-1, 0, -1, 0, 0,  0, 0,-1   ,
     +           0,-1, 0,  1, 0, 0,  0, 0, 1   ,
     +           0, 1, 0, -1, 0, 0,  0, 0, 1   ,
     +           0, 0, 1,  0,-1, 0,  1, 0, 0   ,
     +           0, 0,-1,  0,-1, 0, -1, 0, 0   ,
     +           0, 0,-1,  0, 1, 0,  1, 0, 0   ,
     +           0, 0, 1,  0, 1, 0, -1, 0, 0   ,
     +          -1, 0, 0,  0, 0, 1,  0, 1, 0   ,
     +          -1, 0, 0,  0, 0,-1,  0,-1, 0   ,
     +           1, 0, 0,  0, 0,-1,  0, 1, 0   ,
     +           1, 0, 0,  0, 0, 1,  0,-1, 0   ,
     +           0, 0, 1,  1, 0, 0,  0, 1, 0   ,
     +           0, 0,-1, -1, 0, 0,  0, 1, 0   ,
     +           0, 0,-1,  1, 0, 0,  0,-1, 0   ,
     +           0, 0, 1, -1, 0, 0,  0,-1, 0   ,
     +           0, 1, 0,  0, 0, 1,  1, 0, 0   ,
     +           0,-1, 0,  0, 0,-1,  1, 0, 0   ,
     +           0,-1, 0,  0, 0, 1, -1, 0, 0   ,
     +           0, 1, 0,  0, 0,-1, -1, 0, 0   /


      do j = 1,3
         do k = 1,3
            overlap(k,j) = at(1,k)*at(1,j) + 
     +                     at(2,k)*at(2,j) + 
     +                     at(3,k)*at(3,j) 
         end do
      end do

      call matinv(overlap,3,3,iwork)
      nrot = 1
      do n = 1,24
         do j = 1,3
            do m = 1,3
               rat(m) = s(m,1,n)*at(1,j) + 
     +                  s(m,2,n)*at(2,j) +
     +                  s(m,3,n)*at(3,j)
            end do
            do k = 1,3
               rot(k,j) = at(1,k)*rat(1) + 
     +                    at(2,k)*rat(2) + 
     +                    at(3,k)*rat(3) 
            end do
         end do
         do j = 1,3
            do k = 1,3
               value = overlap(j,1)*rot(1,k) + 
     +                 overlap(j,2)*rot(2,k) +
     +                 overlap(j,3)*rot(3,k)
               if (abs(float(nint(value))-value).gt.1.0e-8) then
* if a noninteger is obtained, this implies that this operation 
* is not a symmetry operation for the given lattice
                  go to 10
               end if
               is(k,j,nrot) = nint(value)
            end do
         end do
         nrot = nrot + 1
  10     continue
      end do
      nrot=nrot-1

      return
      end

c
c-----------------------------------------------------------------------
      subroutine coset(nsym,nrot,irg,s)
c-----------------------------------------------------------------------
c
      implicit none
      integer nsym,ir,is,nc,nrot
      integer nwsym, s(3,3,48), irg(nsym,*)
      logical flag(48)
c
      do 10 ir=1, nrot
         flag(ir) = .true.
  10  continue
      do 11 is=1, nsym
         flag(irg(is,1)) = .false.
  11  continue
c
      do 35 nc =2, nrot/nsym
         do 33 ir = 1, nrot
            if(flag(ir)) then
               do 31 is = 1, nsym
                  irg(is,nc) = nwsym(irg(is,1),ir,s,nrot)
                  flag(irg(is,nc)) = .false.
  31           continue
               go to 35
            end if
  33     continue
  35  continue
c
      return
      end
*-----------------------------------------------------------------------
      logical function eqvect(x,y,f)
*-----------------------------------------------------------------------
      real*8 x(3), y(3), f(3)
      eqvect =abs(x(1)-y(1)-f(1)-nint(x(1)-y(1)-f(1))).lt.1.0e-5 .and.
     &        abs(x(2)-y(2)-f(2)-nint(x(2)-y(2)-f(2))).lt.1.0e-5 .and.
     &        abs(x(3)-y(3)-f(3)-nint(x(3)-y(3)-f(3))).lt.1.0e-5
      return
      end
c
c----------------------------------------------------------------------
      subroutine error(routin,messag,ierr)
c----------------------------------------------------------------------
c
      character*(*) routin, messag
      if(ierr.eq.0) return
      write(6,*) ' '
      write(6,'(1x,79(''!''))')
      write(6,'(5x,''from '',a,'' : error #'',i10)') routin,ierr
      write(6,'(5x,a)') messag
      write(6,'(1x,79(''!''))')
      if(ierr.gt.0) then
         write(6,'(''     stopping ...'')')
         stop
      else
         write(6,*) ' '
         return
      end if
      end

c
c-----------------------------------------------------------------------
      integer function nwsym(ir,is,s,nrot)
c-----------------------------------------------------------------------
c
      implicit none
      integer s(3,3,48), ss(3,3)
      integer i,j,jr,is,ir,nrot
      do 11 i=1,3
         do 10 j=1,3
            ss(i,j) = s(i,1,is)*s(1,j,ir) + s(i,2,is)*s(2,j,ir) +
     &                s(i,3,is)*s(3,j,ir)
  10     continue
  11  continue
c
      do 22 jr=1, nrot
         do 21 i=1,3
            do 20 j=1,3
               if( ss(i,j) .ne. s(i,j,jr) ) go to 22
  20        continue
  21     continue
         nwsym=jr
         return
  22  continue
c
      return
      end

	function ddsum(n,vect,inc)

	implicit none
        integer n, inc, i
	real*8 vect(n),sum,ddsum


	if (n.lt.0) return
        if (inc.lt.0) return

	sum=0
	
	do 10 i=1,n,inc
	sum=sum+vect(i)
 10	continue
	ddsum=sum
	return
	end
c
c---------------------------------------------------------------------
      subroutine recips_f(a1,a2,a3,b1,b2,b3)
c---------------------------------------------------------------------
c
c   generates the reciprocal lattice vectors b1,b2,b3 given the real
c   space vectors a1,a2,a3.  the b's are units of 2pi/a.
      implicit none
      integer i,j,k,iperm,ir,l
      real*8 s,den,a1(3),a2(3),a3(3),b1(3),b2(3),b3(3)
      den=0
      i=1
      j=2
      k=3
      s=1.0
    1 do 2 iperm=1,3
      den=den+s*a1(i)*a2(j)*a3(k)
      l=i
      i=j
      j=k
      k=l
    2 continue
      i=2
      j=1
      k=3
      s=-s
      if(s.lt.0.0) go to 1
      i=1
      j=2
      k=3
      do 5 ir=1,3
      b1(ir)=(a2(j)*a3(k)-a2(k)*a3(j)) / den
      b2(ir)=(a3(j)*a1(k)-a3(k)*a1(j)) / den
      b3(ir)=(a1(j)*a2(k)-a1(k)*a2(j)) / den
      l=i
      i=j
      j=k
      k=l
    5 continue
      return
      end
c
c         *****     end of dmixp package     *****
c
c
c*********************************************************************
c***          direct- and reciprocal-lattice section               ***
c*********************************************************************
c
c---------------------------------------------------------------------
      subroutine hexsym(at,is,nrot)
* Provides symmetry operations for Hexagonal and Trigonal lattices
* (excepted inversion). The c axis is assumed to be along the z axis
      implicit none
      real*8 at(3,3), s(3,3,12), overlap(3,3), rat(3), rot(3,3), value
     ,      ,sin3, cos3, msin3, mcos3
* sin3 = sin(pi/3) and so on
      integer is(3,3,48), iwork(3), j, m, k, n, nrot
      parameter ( sin3 = 0.866025403784438597d0, cos3 = 0.5, 
     &           msin3 =-0.866025403784438597d0,mcos3 =-0.5  )
      external matinv
      data s/ 
     +           1, 0, 0,  0, 1, 0,  0, 0, 1   ,
     +          -1, 0, 0,  0,-1, 0,  0, 0, 1   ,
     +          -1, 0, 0,  0, 1, 0,  0, 0,-1   ,
     +           1, 0, 0,  0,-1, 0,  0, 0,-1   ,
     +    cos3, sin3, 0, msin3, cos3, 0, 0, 0, 1,
     +    cos3,msin3, 0,  sin3, cos3, 0, 0, 0, 1,
     +   mcos3, sin3, 0, msin3,mcos3, 0, 0, 0, 1,
     +   mcos3,msin3, 0,  sin3,mcos3, 0, 0, 0, 1,
     +    cos3,msin3, 0, msin3,mcos3, 0, 0, 0,-1,
     +    cos3, sin3, 0,  sin3,mcos3, 0, 0, 0,-1,
     +   mcos3,msin3, 0, msin3, cos3, 0, 0, 0,-1,
     +   mcos3, sin3, 0,  sin3, cos3, 0, 0, 0,-1  /


      do j = 1,3
         do k = 1,3
            overlap(k,j) = at(1,k)*at(1,j) + 
     +                     at(2,k)*at(2,j) + 
     +                     at(3,k)*at(3,j) 
         end do
      end do

      call matinv(overlap,3,3,iwork)

      nrot = 1
      do n = 1,12
         do j = 1,3
            do m = 1,3
               rat(m) = s(m,1,n)*at(1,j) + 
     +                  s(m,2,n)*at(2,j) +
     +                  s(m,3,n)*at(3,j)
            end do
            do k = 1,3
               rot(k,j) = at(1,k)*rat(1) + 
     +                    at(2,k)*rat(2) + 
     +                    at(3,k)*rat(3) 
            end do
         end do
         do j = 1,3
            do k = 1,3
               value = overlap(j,1)*rot(1,k) + 
     +                 overlap(j,2)*rot(2,k) +
     +                 overlap(j,3)*rot(3,k)
               if (abs(float(nint(value))-value).gt.1.0e-8) then
* if a noninteger is obtained, this implies that this operation 
* is not a symmetry operation for the given lattice
                  go to 10
               end if
               is(k,j,nrot) = nint(value)
            end do
         end do
         nrot = nrot+1
  10     continue
      end do
      nrot = nrot-1

      return
      end
*
      SUBROUTINE MATINV(A,LDA,N,IWORK)
       implicit none
       integer lda, n, info
       real*8 a(LDA,*),work1(256)
       integer iwork(*)
       integer niwork(256)
       if(n .ge. 256) then
         CALL ERROR('MATINV','N GT 256',N)
       endif
       call dgetrf(n, n, a, lda, niwork, info)
       CALL ERROR('MATINV','INFO NE 0',INFO)
       call dgetri(n, a, lda, niwork, work1, n, info)

      RETURN
      END
*
*-----------------------------------------------------------------------
      subroutine symrho(rho,nr1,nr2,nr3,nsym,s,irg,ftau)
*-----------------------------------------------------------------------
*
*     symmetrize the charge density.
*
CDIR$ inline
      implicit none
      integer nr1, nr2, nr3, nsym, irg(48), s(3,3,48),
     +        ftau(3,48), irot, isym, i, j, k, ri, rj, rk
      real*8 rho(nr1,nr2,nr3), sum
      external ruotaijk
*
      if (nsym.eq.1) return
*
      do 30 k = 1,nr3
         do 20 j = 1,nr2
            do 10 i = 1,nr1
               rho(i,j,k) = -rho(i,j,k)
   10       continue
   20    continue
   30 continue
*

      do 80 k = 1,nr3
         do 70 j = 1,nr2
            do 60 i = 1,nr1
*
               if (rho(i,j,k).lt.0) then
                  sum = 0.0
                  do 40 isym = 1,nsym
                     irot = irg(isym)
                     call ruotaijk(s(1,1,irot),ftau(1,irot),i,j,k,nr1,
     +               nr2,nr3,ri,rj,rk)

                     sum = sum - rho(ri,rj,rk)
   40             continue


*
*     sum contains the symmetrized charge density at point r.
*     now fill the star of r with this sum.
*
                  do 50 isym = 1,nsym
                     irot = irg(isym)
                     call ruotaijk(s(1,1,irot),ftau(1,irot),i,j,k,nr1,
     +               nr2,nr3,ri,rj,rk)
                     rho(ri,rj,rk) = sum
   50             continue
               end if
*
   60       continue
   70    continue
   80 continue
*

      return
      end
*
      subroutine ruotaijk(s,ftau,i,j,k,nr1,nr2,nr3,ri,rj,rk)
      implicit none
      integer s(3,3), ftau(3), i, j, k, nr1, nr2, nr3, ri, rj, rk
*
      ri=s(1,1)*(i-1)+s(2,1)*(j-1)+s(3,1)*(k-1) - ftau(1)
      ri=mod(ri,nr1)+1
      if(ri.lt.1) ri=ri+nr1
      rj=s(1,2)*(i-1)+s(2,2)*(j-1)+s(3,2)*(k-1) - ftau(2)
      rj=mod(rj,nr2)+1
      if(rj.lt.1) rj=rj+nr2
      rk=s(1,3)*(i-1)+s(2,3)*(j-1)+s(3,3)*(k-1) - ftau(3)
      rk=mod(rk,nr3)+1
      if(rk.lt.1) rk=rk+nr3
*
      return
      end

c-----------------------------------------------------------------------
      subroutine fsymforces (force,s,irg,irt,nat,ibrav,nsym,
     &                       celldm,nr1,nr2,nr3)
c-----------------------------------------------------------------------
 
* This subroutine symmetrize the force vector for every atom.
* Last revised by MBN - may 1996

      integer s(3,3,48), irt(48,nat), irg(48)
      real*8 force(3,nat),at(3,3),bg(3,3),celldm(6),omega
cEmil      intrinsic my_pe
      integer fmy_pe
 

* generate direct space vectors
      call latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega,1)
* reciprocal space vectors
      call recips_f(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
 
*
* transform to crystal axis...
*
      do na=1, nat
         call trnvect(force(1,na),at,bg,-1)
      end do
*
* ...symmetrize...
*
      call symvect(nat,force,nsym,s,irg,irt)
*
* ... and transform back to cartesian axis
*
      do na=1, nat
         call trnvect(force(1,na),at,bg, 1)
      end do
      return
      end
*

      subroutine trnvect(vect,at,bg,iflag)
*-----------------------------------------------------------------------
*
* transforms a vector from crystal to cartesian axis (iflag.gt.0) and viceversa
*
      implicit none
      integer iflag, i, j
      real*8 vect(3), work(3), at(3,3), bg(3,3)
*
      if(iflag.gt.0) then
*
* forward transformation
*
         do 10 i = 1,3
            work(i)=vect(i)*sqrt(at(1,i)**2+at(2,i)**2+at(3,i)**2)
  10     continue
*
         do 20 i = 1,3
            vect(i) = 0.0
            do 30 j = 1,3
               vect(i) = vect(i) + work(j)*bg(i,j)
 30         continue
 20      continue
      else
*
* backward transformation
*
         do 110 i = 1,3
            work(i) = 0.0
            do 120 j = 1,3
               work(i) = work(i) + vect(j)*at(j,i)
 120        continue
 110     continue
*
         do 130 i=1,3
            vect(i)=work(i)/sqrt(at(1,i)**2+at(2,i)**2+at(3,i)**2)
 130     continue
      end if
*
      return
      end
           
*
*-----------------------------------------------------------------------
      subroutine symvect(nat,vect,nsym,s,irg,irt)
*-----------------------------------------------------------------------
*
* symmetrizes a vector in the crystal axis basis
*
      implicit none
      integer nat,nsym,na,nar,isym,irot,l,irg(48),irt(48,nat),
     &        s(3,3,48)
      integer idx, ii
      real*8 rnsym, vect(3,nat), work(3,600)
*
      if(nsym.eq.1) return
      if(nat .ge. 600)then
        print *, 'Error in symvect NAT .ge. 600'
        stop
      endif
*
      do idx=1,nat
        work(1,idx) = 0.0
        work(2,idx) = 0.0
        work(3,idx) = 0.0
      enddo
      do 100 na=1, nat
         do 200 isym = 1,nsym
            irot = irg(isym)
            nar = irt(irot,na)
            do 190 l = 1,3
               work(l,na) = work(l,na) + s(l,1,irot)*vect(1,nar) +
     +         s(l,2,irot)*vect(2,nar) + s(l,3,irot)*vect(3,nar)
  190       continue
  200    continue
  100 continue
      rnsym = 1.0/nsym
      do idx=1,nat
        work(1,idx) = work(1,idx) * rnsym
        work(2,idx) = work(2,idx) * rnsym
        work(3,idx) = work(3,idx) * rnsym
        vect(1,idx) = work(1,idx)
        vect(2,idx) = work(2,idx)
        vect(3,idx) = work(3,idx)
      enddo
*
      return
      end
           
