c234567890

      integer i,j,k,l,tmpline,tmpint,tmpint1,tmpint2,ttt,tta,ttb,ttu
      integer countat,begin,endd,countchgrp,pdbnCGbeads
      integer atom1,atom2,atom3,atom4
      real tmpdist,tmpdistx,tmpdisty,tmpdistz
 
c DATA for the pdb file
      integer pdbnres,pdbnatom,pdbatnumb(10000),pdbresnum(10000)
      real ax(10000),ay(10000),az(10000)
      character*4 tmp
      character*3 pdbatname(10000),pdbresname(10000)
c CG FF data
      character*3 resname(20)
      integer nbeads(20),nbond(20),nangle(20),ndihe(20),nconst(20)
      character*3 tmpchar,atomname(20,7)
      character*3 beadtype(20,7)
      integer chargrp(20,7)
      real charbead(20,7)
      integer bondatom1(20,10),bondatom2(20,10),bondtype(20,10)
      real bondist(20,10),bondk(20,10)
      integer constatom1(20,10),constatom2(20,10),consttype(20,10)
      real constdist(20,10)
      integer angleatom1(20,10),angleatom2(20,10),angleatom3(20,10)
      integer angletype(20,10)
      real angle(20,10),anglek(20,10)
      integer diheatom1(20,10),diheatom2(20,10),diheatom3(20,10)
      integer diheatom4(20,10),dihetype(20,10)
      real dihe(20,10),dihek(20,10)
c DATA for topology (top)
      character*3 topatname(10000),tpresname(1000),CGname(1000)
      real atcharge(1000,7)
      integer tprestype(1000),tpnbeads(1000),beginres(1000)
      character*3 tpbeadtype(1000,7),tpatomname(1000,7)
      integer tpatnumb(1000,7),tpchrgrp(1000,7)
      real tpatcharge(1000,7)
      real*8 vector1(3),vector2(3),angleCACACA,angleCACACB,angleCBCACA
c AA FF data
      integer tmp2,natomresAA(20),pdbnatomAA,pdbnresAA
      integer natomresCG(20)
      character*3 nameresAA(20)
      character*120 tmpbigchar
c Stuff for the COMs
      integer countatcom
      integer countatCG,countatAA,countatSC
      real tmpcomx,tmpcomy,tmpcomz
      real tmpcom1x,tmpcom1y,tmpcom1z
      real tmpcom2x,tmpcom2y,tmpcom2z
      real tmpcom3x,tmpcom3y,tmpcom3z
      real tmpcom4x,tmpcom4y,tmpcom4z
      real comBBx(1000),comBBy(1000),comBBz(1000)
      real comSCx(1000,4),comSCy(1000,4),comSCz(1000,4)
c Paramters
      parameter     (pi=3.141592654)


C--------------------------------------------------------------
C This program is ment to read a pdb file and generate a new
C pdb file corresponding to the topology of FFCG-2.1.
C The main idea is to use the pdb structure to contruct the
C beads of the CG model.
C All comments to x.periole@rug.nl; NYC-CSICUNY-03-05-2007
C--------------------------------------------------------------
c  Important Notes for Correct Functioning:
c  0- OF COURSE the use of this program is at your own risks
c  1- the pdb file should be cleaned of hydrogens, if any,
c     HEADER/REMARK etc, and Nterm and Cterm atoms.
c  3- the first line of the pdb file must contain the number 
c     of residues and the number of atoms.
c  4- The pdb file should not contain a CHAIN-ID (between the 
c     resname and the resnumber) 
C--------------------------------------------------------------

C  Define stuff for the CG model

C  Integrate the info of the CG FF
c  Note : the resname are red in cg.dat and then they are 
c  refered as the index as the order of appearance
C 

      open(unit=11,name="cg-2.1.dat",status='old')
c     read(11,*) tmpline

      do i=1,20
      read(11,222) resname(i),nbeads(i),nbond(i),nangle(i),ndihe(i)
     .,nconst(i)
c     write(6,222) resname(i),nbeads(i),nbond(i),nangle(i),ndihe(i)
c    .,nconst(i)
222   format(A3,4X,I1,3X,I1,3X,I1,3x,I1,3x,I1)
      do j=1,nbeads(i)
      read(11,*) tmpint1,beadtype(i,j),tmpint2,tmpchar,atomname(i,j)
     .,chargrp(i,j),charbead(i,j)
c     write(6,*) tmpint1,beadtype(i,j),tmpint2,tmpchar,atomname(i,j)
c    .,chargrp(i,j),charbead(i,j)
      enddo
      do j=1,nbond(i)
      read(11,*) bondatom1(i,j),bondatom2(i,j),bondtype(i,j)
     .,bondist(i,j),bondk(i,j)
c     write(6,*) bondatom1(i,j),bondatom2(i,j),bondtype(i,j)
c    ,,bondist(i,j),bondk(i,j)
      enddo
      do j=1,nconst(i)
      read(11,*) constatom1(i,j),constatom2(i,j),consttype(i,j)
     .,constdist(i,j)
c     write(6,*) constatom1(i,j),constatom2(i,j),consttype(i,j)
c    .,constdist(i,j)
      enddo
      do j=1,nangle(i)
      read(11,*) angleatom1(i,j),angleatom2(i,j),angleatom3(i,j)
     .,angletype(i,j),angle(i,j),anglek(i,j)
c     write(6,*) angleatom1(i,j),angleatom2(i,j),angleatom3(i,j)
c    .,angletype(i,j),angle(i,j),anglek(i,j)
      enddo
      do j=1,ndihe(i)
      read(11,*) diheatom1(i,j),diheatom2(i,j),diheatom3(i,j)
     .,diheatom4(i,j),dihetype(i,j),dihe(i,j),dihek(i,j)
c     write(6,*) diheatom1(i,j),diheatom2(i,j),diheatom3(i,j)
c    .,diheatom4(i,j),dihetype(i,j),dihe(i,j),dihek(i,j)
      enddo
      enddo

C  Read the info concerning residues

      open(unit=22,name="AA-2.1.dat",status='old')

      read (22,*) tmpbigchar
      read (22,*) tmpbigchar
      read (22,*) tmpbigchar
      read (22,*) tmp2
      do i = 1, tmp2
      read (22,*) nameresAA(i),natomresAA(i),natomresCG(i)
c     write (6,*) nameresAA(i),natomresAA(i),natomresCG(i)
      enddo
     
C  Read the pdb file -

c first read the number of residues and atoms

      read (5,*)  pdbnresAA,pdbnatomAA
      if (pdbnresAA .le. 0 .or. pdbnatomAA .le. 0 ) then 
      write(6,*) "The PDB file you give as an input is not formatted
     .correctly"
      write(6,*) "Tip: The first line should be the number of residues 
     .and atoms."
      stop
      endif

      do i = 1,pdbnatomAA
c     read (5,*) tmp,pdbatnumb(i),pdbatname(i),pdbresname(i),pdbresnum(i),ax(i),
c    .ay(i),az(i)
      read (5,662) pdbatnumb(i),pdbatname(i),pdbresname(i),pdbresnum(i),
     .ax(i),ay(i),az(i)
c     write (6,662) pdbatnumb(i),pdbatname(i),pdbresname(i),pdbresnum(i)
c    .,ax(i),ay(i),az(i)
      enddo

C  getting the info for each residue

c  here the counter of beads for the CG
      pdbnCGbeads=0
      countat = 1
      do i = 1,pdbnresAA 

c define the residue type

       beginres(i)=countat
       tpresname(i)=pdbresname(beginres(i))

        do j =1,20
        if (tpresname(i) .eq. resname(j)) then
        tprestype(i) = j
c       write(6,*) i,tpresname(i),resname(j),tprestype(i),countat
        endif
        enddo

        if ((tprestype(i) .lt. 1) .and. (tprestype(i) .gt. 20)) then
        write (6,*) "Something is wrong with the attribution of the
     .residue index corresponding to its name in the pdb file!"
        stop
        endif

c here adding the number of beads for the current residue

        pdbnCGbeads=pdbnCGbeads+nbeads(tprestype(i))

c initialize current COM for CA

        tmpcomx = 0.0
        tmpcomy = 0.0
        tmpcomz = 0.0
        countatcom = 0.0

c get the COM  for BB

        do k=1,natomresAA(tprestype(i))
c        write(6,663) pdbatname(beginres(i)+k-1)
         if((pdbatname(beginres(i)+k-1) .eq. "CA ")) then
c        write(6,663) pdbatname(beginres(i)+k-1)
         tmpcomx = ax(beginres(i)+k-1)
         tmpcomy = ay(beginres(i)+k-1)
         tmpcomz = az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
         countat = countat + 1
        enddo
663   format(A4)

c calculate the COM for CA ; note the unusual way to calculate(old stuff)

         if(countatcom .eq. 1) then 
         comBBx(i) = tmpcomx
         comBBy(i) = tmpcomy
         comBBz(i) = tmpcomz
         endif

C Sidechains
c initialize the COM for the SC

        tmpcomx = 0.0
        tmpcomy = 0.0
        tmpcomz = 0.0
        tmpcom1x = 0.0
        tmpcom1y = 0.0
        tmpcom1z = 0.0
        tmpcom2x = 0.0
        tmpcom2y = 0.0
        tmpcom2z = 0.0
        tmpcom3x = 0.0
        tmpcom3y = 0.0
        tmpcom3z = 0.0
        tmpcom4x = 0.0
        tmpcom4y = 0.0
        tmpcom4z = 0.0
        countatcom = 0.0
        countattmp = 0.0

c GLY; ALA no COM for SC

c CYS
       if(tprestype(i) .eq. 3) then
        do k=1,natomresAA(tprestype(i))
         if(pdbatname(beginres(i)+k-1) .eq. "CB "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "SG "  ) then 
         tmpcomx = tmpcomx + ax(beginres(i)+k-1)
         tmpcomy = tmpcomy + ay(beginres(i)+k-1)
         tmpcomz = tmpcomz + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
        enddo

c calculate the COM for SC

         if(countatcom .eq. 2) then
         comSCx(i,1) = tmpcomx/2.0
         comSCy(i,1) = tmpcomy/2.0
         comSCz(i,1) = tmpcomz/2.0
         else
         write(6,*) "something wrong in the Cys ",i
         stop
         endif
         countattmp = countattmp + 1
       endif

c VAL
       if(tprestype(i) .eq. 4) then
        do k=1,natomresAA(tprestype(i))
         if(pdbatname(beginres(i)+k-1) .eq. "CB "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CG1"  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CG2"  ) then
         tmpcomx = tmpcomx + ax(beginres(i)+k-1)
         tmpcomy = tmpcomy + ay(beginres(i)+k-1)
         tmpcomz = tmpcomz + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
         countattmp = countattmp + 1
        enddo

c calculate the COM for SC

         if(countatcom .eq. 3) then
         comSCx(i,1) = tmpcomx/3.0
         comSCy(i,1) = tmpcomy/3.0
         comSCz(i,1) = tmpcomz/3.0
         else
         write(6,*) "something wrong in the Val ",i
         stop
         endif
       endif

c LEU 

       if(tprestype(i) .eq. 5) then
        do k=1,natomresAA(tprestype(i))
         if(pdbatname(beginres(i)+k-1) .eq. "CB "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CG "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CD1"  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CD2"  ) then 
         tmpcomx = tmpcomx + ax(beginres(i)+k-1)
         tmpcomy = tmpcomy + ay(beginres(i)+k-1)
         tmpcomz = tmpcomz + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
         countattmp = countattmp + 1
        enddo

c calculate the COM for SC

         if(countatcom .eq. 4) then
         comSCx(i,1) = tmpcomx/4.0
         comSCy(i,1) = tmpcomy/4.0
         comSCz(i,1) = tmpcomz/4.0
         else
         write(6,*) "something wrong in the Leu ",i
         stop
         endif
       endif

c ILE
       if(tprestype(i) .eq. 6) then
        do k=1,natomresAA(tprestype(i))
c        write(6,*) pdbatname(beginres(i)+k-1)
         if(pdbatname(beginres(i)+k-1) .eq. "CB "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CG1"  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CG2"  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CD1"  ) then 
c        write(6,*) pdbatname(beginres(i)+k-1)
         tmpcomx = tmpcomx + ax(beginres(i)+k-1)
         tmpcomy = tmpcomy + ay(beginres(i)+k-1)
         tmpcomz = tmpcomz + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
         countattmp = countattmp + 1
        enddo

c calculate the COM for SC

         if(countatcom .eq. 4) then
         comSCx(i,1) = tmpcomx/4.0
         comSCy(i,1) = tmpcomy/4.0
         comSCz(i,1) = tmpcomz/4.0
         else
         write(6,*) "something wrong in the Ile ",i
         stop
         endif
       endif

c MET
       if(tprestype(i) .eq. 7) then
        do k=1,natomresAA(tprestype(i))
c        write(6,*) pdbatname(beginres(i)+k-1)
c        write(6,*) tmpcomx,tmpcomy,tmpcomz
         if(pdbatname(beginres(i)+k-1) .eq. "CB "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CG "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "SD "  .or.  
     .      pdbatname(beginres(i)+k-1) .eq. "CE "  ) then
         tmpcomx = tmpcomx + ax(beginres(i)+k-1)
         tmpcomy = tmpcomy + ay(beginres(i)+k-1)
         tmpcomz = tmpcomz + az(beginres(i)+k-1)
c        write(6,*) tmpcomx,tmpcomy,tmpcomz
         countatcom = countatcom+1
         endif
         countattmp = countattmp + 1
        enddo

c calculate the COM for SC

         if(countatcom .eq. 4) then
         comSCx(i,1) = tmpcomx/4.0
         comSCy(i,1) = tmpcomy/4.0
         comSCz(i,1) = tmpcomz/4.0
         else
         write(6,*) "something wrong in the Met ",i
         stop
         endif
       endif

c PRO
       if(tprestype(i) .eq. 8) then
        do k=1,natomresAA(tprestype(i))
         if(pdbatname(beginres(i)+k-1) .eq. "CB "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CG "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CD "  ) then
         tmpcomx = tmpcomx + ax(beginres(i)+k-1)
         tmpcomy = tmpcomy + ay(beginres(i)+k-1)
         tmpcomz = tmpcomz + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
         countattmp = countattmp + 1
        enddo

c calculate the COM for SC

         if(countatcom .eq. 3) then
         comSCx(i,1) = tmpcomx/3.0
         comSCy(i,1) = tmpcomy/3.0
         comSCz(i,1) = tmpcomz/3.0
         else
         write(6,*) "something wrong in the Pro ",i
         stop
         endif
       endif

c ASN
       if(tprestype(i) .eq. 9) then
        do k=1,natomresAA(tprestype(i))
         if(pdbatname(beginres(i)+k-1) .eq. "CB "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CG "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "OD1"  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "ND2"  ) then
         tmpcomx = tmpcomx + ax(beginres(i)+k-1)
         tmpcomy = tmpcomy + ay(beginres(i)+k-1)
         tmpcomz = tmpcomz + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
         countattmp = countattmp + 1
        enddo

c calculate the COM for SC

         if(countatcom .eq. 4) then
         comSCx(i,1) = tmpcomx/4.0
         comSCy(i,1) = tmpcomy/4.0
         comSCz(i,1) = tmpcomz/4.0
         else
         write(6,*) "something wrong in the Asn ",i
         stop
         endif
       endif

c GLN
       if(tprestype(i) .eq. 10) then
        do k=1,natomresAA(tprestype(i))
         if(pdbatname(beginres(i)+k-1) .eq. "CB "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CG "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CD "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "OE1"  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "NE2"  ) then
         tmpcomx = tmpcomx + ax(beginres(i)+k-1)
         tmpcomy = tmpcomy + ay(beginres(i)+k-1)
         tmpcomz = tmpcomz + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
         countattmp = countattmp + 1
        enddo

c calculate the COM for SC

         if(countatcom .eq. 5) then
         comSCx(i,1) = tmpcomx/5.0
         comSCy(i,1) = tmpcomy/5.0
         comSCz(i,1) = tmpcomz/5.0
         else
         write(6,*) "something wrong in the Gln ",i
         stop
         endif
       endif

c LYS ; here two beads! 
       if(tprestype(i) .eq. 11) then
        do k=1,natomresAA(tprestype(i))
c COM1
         if(pdbatname(beginres(i)+k-1) .eq. "CB "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CG "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CD "  ) then
         tmpcom1x = tmpcom1x + ax(beginres(i)+k-1)
         tmpcom1y = tmpcom1y + ay(beginres(i)+k-1)
         tmpcom1z = tmpcom1z + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
c COM2  
         if(pdbatname(beginres(i)+k-1) .eq. "CE "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "NZ "  ) then
         tmpcom2x = tmpcom2x + ax(beginres(i)+k-1)
         tmpcom2y = tmpcom2y + ay(beginres(i)+k-1)
         tmpcom2z = tmpcom2z + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
         countattmp = countattmp + 1
        enddo

c calculate the COM for SC

         if(countatcom .eq. 5) then
         comSCx(i,1) = tmpcom1x/3.0
         comSCy(i,1) = tmpcom1y/3.0
         comSCz(i,1) = tmpcom1z/3.0
         comSCx(i,2) = tmpcom2x/2.0
         comSCy(i,2) = tmpcom2y/2.0
         comSCz(i,2) = tmpcom2z/2.0
         else
         write(6,*) "something wrong in the Lys ",i
         stop
         endif
       endif

c ASP
       if(tprestype(i) .eq. 12) then
        do k=1,natomresAA(tprestype(i))
         if(pdbatname(beginres(i)+k-1) .eq. "CB "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CG "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "OD1"  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "OD2"  ) then
         tmpcomx = tmpcomx + ax(beginres(i)+k-1)
         tmpcomy = tmpcomy + ay(beginres(i)+k-1)
         tmpcomz = tmpcomz + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
         countattmp = countattmp + 1
        enddo

c calculate the COM for SC

         if(countatcom .eq. 4) then
         comSCx(i,1) = tmpcomx/4.0
         comSCy(i,1) = tmpcomy/4.0
         comSCz(i,1) = tmpcomz/4.0
         else
         write(6,*) "something wrong in the Asp ",i
         stop
         endif
       endif

c GLU
       if(tprestype(i) .eq. 13) then
        do k=1,natomresAA(tprestype(i))
         if(pdbatname(beginres(i)+k-1) .eq. "CB "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CG "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CD "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "OE1"  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "OE2"  ) then
         tmpcomx = tmpcomx + ax(beginres(i)+k-1)
         tmpcomy = tmpcomy + ay(beginres(i)+k-1)
         tmpcomz = tmpcomz + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
         countattmp = countattmp + 1
        enddo

c calculate the COM for SC

         if(countatcom .eq. 5) then
         comSCx(i,1) = tmpcomx/5.0
         comSCy(i,1) = tmpcomy/5.0
         comSCz(i,1) = tmpcomz/5.0
         else
         write(6,*) "something wrong in the Glu ",i
         stop
         endif
       endif

c THR
       if(tprestype(i) .eq. 14) then
        do k=1,natomresAA(tprestype(i))
         if(pdbatname(beginres(i)+k-1) .eq. "CB "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "OG1"  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CG2"  ) then
         tmpcomx = tmpcomx + ax(beginres(i)+k-1)
         tmpcomy = tmpcomy + ay(beginres(i)+k-1)
         tmpcomz = tmpcomz + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
         countattmp = countattmp + 1
        enddo

c calculate the COM for SC

         if(countatcom .eq. 3) then
         comSCx(i,1) = tmpcomx/3.0
         comSCy(i,1) = tmpcomy/3.0
         comSCz(i,1) = tmpcomz/3.0
         else
         write(6,*) "something wrong in the Thr ",i
         stop
         endif
       endif

c SER
       if(tprestype(i) .eq. 15) then
        do k=1,natomresAA(tprestype(i))
         if(pdbatname(beginres(i)+k-1) .eq. "CB "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "OG "  ) then
         tmpcomx = tmpcomx + ax(beginres(i)+k-1)
         tmpcomy = tmpcomy + ay(beginres(i)+k-1)
         tmpcomz = tmpcomz + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
         countattmp = countattmp + 1
        enddo

c calculate the COM for SC

         if(countatcom .eq. 2) then
         comSCx(i,1) = tmpcomx/2.0
         comSCy(i,1) = tmpcomy/2.0
         comSCz(i,1) = tmpcomz/2.0
         else
         write(6,*) "something wrong in the Ser ",i
         stop
         endif
       endif

c ARG ; here two beads!
       if(tprestype(i) .eq. 16) then
        do k=1,natomresAA(tprestype(i))
c COM1
         if(pdbatname(beginres(i)+k-1) .eq. "CB "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CG "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CD "  ) then
         tmpcom1x = tmpcom1x + ax(beginres(i)+k-1)
         tmpcom1y = tmpcom1y + ay(beginres(i)+k-1)
         tmpcom1z = tmpcom1z + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
c COM2
         if(pdbatname(beginres(i)+k-1) .eq. "NE "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CZ "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "NH1"  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "NH2"  ) then
         tmpcom2x = tmpcom2x + ax(beginres(i)+k-1)
         tmpcom2y = tmpcom2y + ay(beginres(i)+k-1)
         tmpcom2z = tmpcom2z + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
         countattmp = countattmp + 1
        enddo

c calculate the COM for SC

         if(countatcom .eq. 7) then
         comSCx(i,1) = tmpcom1x/3.0
         comSCy(i,1) = tmpcom1y/3.0
         comSCz(i,1) = tmpcom1z/3.0
         comSCx(i,2) = tmpcom2x/4.0
         comSCy(i,2) = tmpcom2y/4.0
         comSCz(i,2) = tmpcom2z/4.0
         else
         write(6,*) "something wrong in the Arg ",i
         stop
         endif
       endif

c HIS
       if(tprestype(i) .eq. 17) then
        do k=1,natomresAA(tprestype(i))
         if(pdbatname(beginres(i)+k-1) .eq. "CB "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "CG "  ) then
         tmpcom1x = tmpcom1x + ax(beginres(i)+k-1)
         tmpcom1y = tmpcom1y + ay(beginres(i)+k-1)
         tmpcom1z = tmpcom1z + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif

         if(pdbatname(beginres(i)+k-1) .eq. "ND1"       ) then
         tmpcom2x = tmpcom2x + ax(beginres(i)+k-1)
         tmpcom2y = tmpcom2y + ay(beginres(i)+k-1)
         tmpcom2z = tmpcom2z + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif

         if(pdbatname(beginres(i)+k-1) .eq. "NE2"       ) then
         tmpcom3x = tmpcom3x + ax(beginres(i)+k-1)
         tmpcom3y = tmpcom3y + ay(beginres(i)+k-1)
         tmpcom3z = tmpcom3z + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
         countattmp = countattmp + 1
        enddo

c calculate the COM for SC

         if(countatcom .eq. 4) then
         comSCx(i,1) = tmpcom1x/2.0
         comSCy(i,1) = tmpcom1y/2.0
         comSCz(i,1) = tmpcom1z/2.0
         comSCx(i,2) = tmpcom2x/1.0
         comSCy(i,2) = tmpcom2y/1.0
         comSCz(i,2) = tmpcom2z/1.0
         comSCx(i,3) = tmpcom3x/1.0
         comSCy(i,3) = tmpcom3y/1.0
         comSCz(i,3) = tmpcom3z/1.0
         else
         write(6,*) "something wrong in the His ",i
         stop
         endif
       endif

c PHE
       if(tprestype(i) .eq. 18) then
        do k=1,natomresAA(tprestype(i))
         if(pdbatname(beginres(i)+k-1) .eq. "CD1"       ) then
         tmpcom1x = tmpcom1x + ax(beginres(i)+k-1)
         tmpcom1y = tmpcom1y + ay(beginres(i)+k-1)
         tmpcom1z = tmpcom1z + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif

         if(pdbatname(beginres(i)+k-1) .eq. "CD2"       ) then
         tmpcom2x = tmpcom2x + ax(beginres(i)+k-1)
         tmpcom2y = tmpcom2y + ay(beginres(i)+k-1)
         tmpcom2z = tmpcom2z + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif

         if(pdbatname(beginres(i)+k-1) .eq. "CZ "       ) then
         tmpcom3x = tmpcom3x + ax(beginres(i)+k-1)
         tmpcom3y = tmpcom3y + ay(beginres(i)+k-1)
         tmpcom3z = tmpcom3z + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
         countattmp = countattmp + 1
        enddo

c calculate the COM for SC

         if(countatcom .eq. 3) then
         comSCx(i,1) = tmpcom1x/1.0
         comSCy(i,1) = tmpcom1y/1.0
         comSCz(i,1) = tmpcom1z/1.0
         comSCx(i,2) = tmpcom2x/1.0
         comSCy(i,2) = tmpcom2y/1.0
         comSCz(i,2) = tmpcom2z/1.0
         comSCx(i,3) = tmpcom3x/1.0
         comSCy(i,3) = tmpcom3y/1.0
         comSCz(i,3) = tmpcom3z/1.0
         else
         write(6,*) "something wrong in the Phe ",i
         stop
         endif
       endif

c TYR
       if(tprestype(i) .eq. 19) then
        do k=1,natomresAA(tprestype(i))
         if(pdbatname(beginres(i)+k-1) .eq. "CD1"       ) then
         tmpcom1x = tmpcom1x + ax(beginres(i)+k-1)
         tmpcom1y = tmpcom1y + ay(beginres(i)+k-1)
         tmpcom1z = tmpcom1z + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif

         if(pdbatname(beginres(i)+k-1) .eq. "CD2"       ) then
         tmpcom2x = tmpcom2x + ax(beginres(i)+k-1)
         tmpcom2y = tmpcom2y + ay(beginres(i)+k-1)
         tmpcom2z = tmpcom2z + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif

         if(pdbatname(beginres(i)+k-1) .eq. "CZ "  .or.
     .      pdbatname(beginres(i)+k-1) .eq. "OH "  ) then
         tmpcom3x = tmpcom3x + ax(beginres(i)+k-1)
         tmpcom3y = tmpcom3y + ay(beginres(i)+k-1)
         tmpcom3z = tmpcom3z + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
         countattmp = countattmp + 1
        enddo

c calculate the COM for SC

         if(countatcom .eq. 4) then
         comSCx(i,1) = tmpcom1x/1.0
         comSCy(i,1) = tmpcom1y/1.0
         comSCz(i,1) = tmpcom1z/1.0
         comSCx(i,2) = tmpcom2x/1.0
         comSCy(i,2) = tmpcom2y/1.0
         comSCz(i,2) = tmpcom2z/1.0
         comSCx(i,3) = tmpcom3x/2.0
         comSCy(i,3) = tmpcom3y/2.0
         comSCz(i,3) = tmpcom3z/2.0
         else
         write(6,*) "something wrong in the Tyr ",i
         stop
         endif
       endif

c TRP
       if(tprestype(i) .eq. 20) then
        do k=1,natomresAA(tprestype(i))
         if(pdbatname(beginres(i)+k-1) .eq. "CG "       ) then
         tmpcom1x = tmpcom1x + ax(beginres(i)+k-1)
         tmpcom1y = tmpcom1y + ay(beginres(i)+k-1)
         tmpcom1z = tmpcom1z + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif

         if(pdbatname(beginres(i)+k-1) .eq. "NE1"       ) then
         tmpcom2x = tmpcom2x + ax(beginres(i)+k-1)
         tmpcom2y = tmpcom2y + ay(beginres(i)+k-1)
         tmpcom2z = tmpcom2z + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif

         if(pdbatname(beginres(i)+k-1) .eq. "CZ2"       ) then
         tmpcom3x = tmpcom3x + ax(beginres(i)+k-1)
         tmpcom3y = tmpcom3y + ay(beginres(i)+k-1)
         tmpcom3z = tmpcom3z + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif

         if(pdbatname(beginres(i)+k-1) .eq. "CE3"       ) then
         tmpcom4x = tmpcom4x + ax(beginres(i)+k-1)
         tmpcom4y = tmpcom4y + ay(beginres(i)+k-1)
         tmpcom4z = tmpcom4z + az(beginres(i)+k-1)
         countatcom = countatcom+1
         endif
         countattmp = countattmp + 1
        enddo

c calculate the COM for SC

         if(countatcom .eq. 4) then
         comSCx(i,1) = tmpcom1x/1.0
         comSCy(i,1) = tmpcom1y/1.0
         comSCz(i,1) = tmpcom1z/1.0
         comSCx(i,2) = tmpcom2x/1.0
         comSCy(i,2) = tmpcom2y/1.0
         comSCz(i,2) = tmpcom2z/1.0
         comSCx(i,3) = tmpcom3x/1.0
         comSCy(i,3) = tmpcom3y/1.0
         comSCz(i,3) = tmpcom3z/1.0
         comSCx(i,4) = tmpcom4x/1.0
         comSCy(i,4) = tmpcom4y/1.0
         comSCz(i,4) = tmpcom4z/1.0
         else
         write(6,*) "something wrong in the Trp ",i
         stop
         endif
       endif

      enddo

C Here we write the pdb file according to the FFCG-2.1

c First the number of residues and the number of beads are written in a HEADER!

      write(6,109) pdbnresAA,pdbnCGbeads
109   format('HEADER number of residues: ',I4," number of beads: ",I4)
      countatCG = 1 
      countatAA = 1 
      do i = 1,pdbnresAA
       ttt = countatCG
       tta = countatAA
       ttu = beginres(i)
c      write(6,662) countatCG,pdbatname(tta),pdbresname(ttu),
       write(6,662) countatCG," CA ",pdbresname(ttu),
     .pdbresnum(ttu),comBBx(i),comBBy(i),comBBz(i)
       countatSC = 1
       countatCG = countatCG + 1
        do j = 1, (natomresCG(tprestype(i))-1)
        ttb = countatSC
        write(6,662) countatCG,"SC ",pdbresname(ttu),
     .pdbresnum(ttu),comSCx(i,ttb),comSCy(i,ttb),comSCz(i,ttb)
        countatSC = countatSC + 1
        countatCG = countatCG + 1
        enddo
c      countat = countat + natomresCG(tprestype(i))
       countatAA = countatAA + natomresAA(tprestype(i))
      enddo

662   format('ATOM',I7,1X,A4,1X,A3,2X,I4,4X,3F8.3)

      end
