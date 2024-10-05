module globalvar
implicit double precision (A-H,O-Z) 
character (len =80) :: rec
! automatic meshing
integer :: icode ,imesh , idmesh 
! geomet
integer :: nnode ,nele , ntype
integer , dimension (:) ,allocatable , save :: imat , itype
integer , dimension (: ,:) ,allocatable , save :: in
real(4) ,dimension (: ,:) ,allocatable , save :: coord , ctype 
! scode 
integer :: ndof
integer , dimension (: ,:) ,allocatable , save :: idof 
! loads
real(4) ,dimension (:) ,allocatable , save :: vloads , therm 
! assemb
real(4) ,dimension (: ,:) ,allocatable , save :: vk 
real(4) ,dimension (6 ,6) :: st
! solve
real(4) ,dimension (:) ,allocatable , save :: vdisp 
REAL(4), DIMENSION(:,:), ALLOCATABLE, SAVE :: VKK
! stress2d
integer :: imean
real(4) :: totalvolume
real(4) ,dimension (:) ,allocatable , save :: sumweight , sigvm
real(4) ,dimension (: ,:) ,allocatable , save :: signod , epsnod
real(4) ,dimension (4) :: eps ,eps0 ,sig
real(4) ,dimension (3) :: sigp , epsp
! joints
real(4) ,dimension (:) ,allocatable , save :: cincl !not used
! evolve
integer :: check ,nrk ,nelek ,nrmin ,iconv , itermax
real(4) :: volumetot ,volk , efactor ,rr ,rmv ,alphan ,v0 ,rr0 ,er ,rrmax ,rmvmax
real(4) ,dimension (:) ,allocatable , save :: ievol 
end module globalvar


program cst 
use globalvar
! title 
character (len =80) :: title 
! opening of the files
open (11 , file =' input.txt ',status ='unknown ')
open (12 , file ='output.txt ',status ='unknown ')
open (9, file ='mkk .txt ',status ='unknown ')
! reading title 
read (11 ,*) title
write (12 , '(1x,a30)') title 
! call the subroutines
call evolve ! reading data
call mesh
call geomet
call scode
call loads
call allocating
call assemb
call joints
! call solve
if( icode /=2) then
call stress2d
if( check ==0) call plot2d
end if
! closing  
close (11)
close (9)
! Use of evolutaionary procedure 
evol :if(check >0) then
if( icode ==2) exit evol 
call evolve
end if evol
! closing last files 
close (12)
end program cst



! subroutine mesh 
subroutine mesh 
use globalvar
integer ::i,j,nnodex ,nnodey ,ind ,ind1 ,ind2 , nsub (2) 
real(4) ::dx ,dy ,x0 ,y0 , dimen (2)
! reading / writing the input data 
read (11 ,*) rec
read (11 ,*) imesh ! type of mesh 
if( imesh ==0) then
write (12 , '(/ ,1x ," Manual Mesh ") ') 
do i=1 ,6
read (11 ,*) rec 
end do
return
else if ( imesh ==1) then
write (12 , '(/ ,1x ," Automatic Mesh ") ') 
end if
read (11 ,*) rec
read (11 ,*) idmesh  ! shape of the mesh
read (11 ,*) rec
read (11 ,*) x0 ,y0
read (11 ,*) rec
read (11 ,*) ( dimen (i),i=1 ,2) ,( nsub (i),i=1 ,2) 
! evaluating number of nodes , elements
dx= dimen (1)/ nsub (1) 
dy= dimen (2)/ nsub (2)
nnodex = nsub (1) +1
nnodey = nsub (2) +1
nnode = nnodex * nnodey 
nele =2* nsub (1) * nsub (2)
! coordinates 
allocate ( coord (nnode ,2))
coord =0. 
ind =0
do i=1, nnodex 
do j=1, nnodey
ind =ind +1
coordx =x0+dx *(i -1)
coordy =y0+dy *(j -1)
coord (ind ,1)= coordx
coord (ind ,2)= coordy
end do 
end do
! incidence
nblockx = nsub (1) /2 
nblocky = nsub (2) /2 
ntriangblock =8
ind =1
ind1 =1
ind2 = ind1
allocate (in(nele ,3)) 
in (: ,:) =0
if( idmesh ==1) then 
do i=1, nblocky
do j=1, nblockx
in(ind ,1) = ind1
in(ind ,2) = ind1 + nnodey
in(ind ,3) = ind1 + nnodey +1
in(ind +1 ,1)= ind1 + nnodey
in(ind +1 ,2)= ind1 +2* nnodey
in(ind +1 ,3)= ind1 + nnodey +1
in(ind +2 ,1)= ind1 +2* nnodey
in(ind +2 ,2)= ind1 +2* nnodey +1
in(ind +2 ,3)= ind1 + nnodey +1
in(ind +3 ,1)= ind1 +2* nnodey +1
in(ind +3 ,2)= ind1 +2* nnodey +2
in(ind +3 ,3)= ind1 + nnodey +1
in(ind +4 ,1)= ind1 +2* nnodey +2
in(ind +4 ,2)= ind1 + nnodey +2
in(ind +4 ,3)= ind1 + nnodey +1
in(ind +5 ,1)= ind1 + nnodey +2
in(ind +5 ,2)= ind1 +2
in(ind +5 ,3)= ind1 + nnodey +1
in(ind +6 ,1)= ind1 +2
in(ind +6 ,2)= ind1 +1
in(ind +6 ,3)= ind1 + nnodey +1 
in(ind +7 ,1)= ind1 +1
in(ind +7 ,2)= ind1
in(ind +7 ,3)= ind1 + nnodey +1 
ind =ind + ntriangblock
ind1 = ind2 +(2* j)* nnodey
end do
ind2 = ind2 +2 
ind1 = ind2
end do
else if( idmesh ==2) then 
do i=1, nblocky
do j=1, nblockx
in(ind ,1) = ind1
in(ind ,2) = ind1 + nnodey
in(ind ,3) = ind1 +1
in(ind +1 ,1)= ind1 + nnodey
in(ind +1 ,2)= ind1 + nnodey +1
in(ind +1 ,3)= ind1 +1
in(ind +2 ,1)= ind1 + nnodey
in(ind +2 ,2)= ind1 +2* nnodey +1
in(ind +2 ,3)= ind1 + nnodey +1
in(ind +3 ,1)= ind1 + nnodey
in(ind +3 ,2)= ind1 + nnodey *2
in(ind +3 ,3)= ind1 + nnodey *2+1
in(ind +4 ,1)= ind1 + nnodey *2+1
in(ind +4 ,2)= ind1 + nnodey *2+2
in(ind +4 ,3)= ind1 + nnodey +2
in(ind +5 ,1)= ind1 + nnodey +1
in(ind +5 ,2)= ind1 + nnodey *2+1
in(ind +5 ,3)= ind1 + nnodey +2
in(ind +6 ,1)= ind1 +1
in(ind +6 ,2)= ind1 + nnodey +1
in(ind +6 ,3)= ind1 + nnodey +2
in(ind +7 ,1)= ind1 +1
in(ind +7 ,2)= ind1 + nnodey +2
in(ind +7 ,3)= ind1 +2 
ind =ind + ntriangblock 
ind1 = ind2 +(2* j)* nnodey
end do
ind2 = ind2 +2 
ind1 = ind2
end do
else if( idmesh ==3) then 
do i=1, nblocky
do j=1, nblockx
in(ind ,1) = ind1
in(ind ,2) = ind1 + nnodey
in(ind ,3) = ind1 + nnodey +1
in(ind +1 ,1)= ind1 + nnodey
in(ind +1 ,2)= ind1 +2* nnodey +1
in(ind +1 ,3)= ind1 + nnodey +1
in(ind +2 ,1)= ind1 + nnodey
in(ind +2 ,2)= ind1 + nnodey *2
in(ind +2 ,3)= ind1 + nnodey *2+1
in(ind +3 ,1)= ind1 +2* nnodey +1
in(ind +3 ,2)= ind1 +2* nnodey +2
in(ind +3 ,3)= ind1 + nnodey +1
in(ind +4 ,1)= ind1 +2* nnodey +2
in(ind +4 ,2)= ind1 + nnodey +2
in(ind +4 ,3)= ind1 + nnodey +1
in(ind +5 ,1)= ind1 +1
in(ind +5 ,2)= ind1 + nnodey +1
in(ind +5 ,3)= ind1 + nnodey +2
in(ind +6 ,1)= ind1 +1
in(ind +6 ,2)= ind1 + nnodey +2
in(ind +6 ,3)= ind1 +2
in(ind +7 ,1)= ind1 +1
in(ind +7 ,2)= ind1
in(ind +7 ,3)= ind1 + nnodey +1 
ind =ind + ntriangblock
ind1 = ind2 +(2* j)* nnodey
end do
ind2 = ind2 +2 
ind1 = ind2
end do
else if( idmesh ==4) then 
do i=1, nblocky
do j=1, nblockx
in(ind ,1) = ind1
in(ind ,2) = ind1 + nnodey
in(ind ,3) = ind1 +1
in(ind +1 ,1)= ind1 + nnodey
in(ind +1 ,2)= ind1 + nnodey +1
in(ind +1 ,3)= ind1 +1
in(ind +2 ,1)= ind1 + nnodey
in(ind +2 ,2)= ind1 + nnodey *2
in(ind +2 ,3)= ind1 + nnodey +1
in(ind +3 ,1)= ind1 + nnodey *2
in(ind +3 ,2)= ind1 + nnodey *2+1
in(ind +3 ,3)= ind1 + nnodey +1
in(ind +4 ,1)= ind1 + nnodey *2+1
in(ind +4 ,2)= ind1 + nnodey *2+2
in(ind +4 ,3)= ind1 + nnodey +2
in(ind +5 ,1)= ind1 + nnodey +1
in(ind +5 ,2)= ind1 + nnodey *2+1
in(ind +5 ,3)= ind1 + nnodey +2
in(ind +6 ,1)= ind1 + nnodey +2
in(ind +6 ,2)= ind1 +2
in(ind +6 ,3)= ind1 + nnodey +1
in(ind +7 ,1)= ind1 +2
in(ind +7 ,2)= ind1 +1
in(ind +7 ,3)= ind1 + nnodey +1 
ind =ind + ntriangblock
ind1 = ind2 +(2* j)* nnodey
end do
ind2 = ind2 +2 
ind1 = ind2
end do 
end if
return
end subroutine mesh




! subroutine geomet 
subroutine geomet 
use globalvar
integer :: idnode ,idtype ,idele ,i,j, idcode
! reding type of problem 
read (11 ,*) rec
read (11 ,*) icode
if ( icode == 0) then
write (12 , '(1x ," PLANE STRESS PROBLEM ") ') 
else if ( icode == 1) then
write (12 , '(1x ," PLANE STRAIN PROBLEM ") ') 
else
write (12 , '(1x ," AXISYMMETRIC PROBLEM ") ') 
end if
! reading the type of weight functions (0= arithmetic average , 1= mean on volume , 2= energy based weight )
read (11 ,*) rec
read (11 ,*) imean
!manual meshing
if( imesh ==0) then
! reding writing nodes 
read (11 ,*) rec
read (11 ,*) nnode 
if(nnode <2) then
write (12 , '(/ ,1x ," incorrect number of nodes ") ') 
stop
return 
end if
end if
write (12 , '(1x ," NUMBER OF NODES =" ,1x,i3)') nnode 
write (12 ,101) 'NODE ','COORDINATE X','COORDINATE Y' 
! automatic meshing
if( imesh ==0) allocate ( coord (nnode ,2) )
do i=1, nnode
if( imesh ==0) then
read (11 ,*) idnode , coord (idnode ,:) 
write (12 ,102) idnode , coord (idnode ,:)
else
write (12 ,102)i ,( coord (i,j),j=1 ,2) 
endif
end do
! reading materials 
read (11 ,*) rec 
read (11 ,*) ntype 
if(ntype <1) then
write (12 , '(/ ,1x ," incorrect number of material type ") ') 
stop
return 
end if
allocate ( ctype (ntype ,5))
write (12 , '(/ ,1x ," NUMBER OF MATERIALS =" ,1x,i3)') ntype
write (12 ,103) 'MAT ','YOUNG ','POISSON ','THERMAL ','WEIGHT ','THICKNESS'
do i=1, ntype
read (11 ,*) idtype , ctype (idtype ,:) 
write (12 ,104) idtype , ctype (idtype ,:) 
! check 
idcode =0
if (( ctype ( idtype ,1) <=0)) idcode =1
if (( ctype ( idtype ,2) <0).or .( ctype (idtype ,2) >0.5)) idcode =2 
if( ctype ( idtype ,3) <0) idcode =3
if( ctype ( idtype ,4) <0) idcode =4 
if( idcode >0) then
write (12 , '(/ ,1x ," error in material type ",i2 ,i1)')idtype , idcode 
stop
return 
end if
end do
if( imesh ==0) then
! R/W connectivity lattrix
read (11 ,*) rec
read (11 ,*) nele 
if(nele <1) then
write (12 , '(/ ,1x ," DATA ERROR specify ") ') 
stop
return 
end if
allocate (in(nele ,3)) 
end if
allocate ( itype ( nele )) 
itype (:) =1
if( imesh ==0) then
write (12 , '(/ ,1x ," NUMBER OF ELEMENTS =" ,1x,i3)') nele 
write (12 ,105) 'ELEM ','NODE1 ','NODE2 ','NODE3 ','TYPE ' 
do i=1, nele
read (11 ,*) idele ,( in(idele ,j),j=1 ,3) ,itype ( idele ) 
write (12 ,106) idele ,( in(idele ,j),j=1 ,3) ,itype ( idele )
end do 
endif

! format list
101 FORMAT (1X,A5 ,5X,A12 ,5X,A12)
102 FORMAT (1X,I5 ,5X,F12 .3 ,5X,F12 .3)
103 FORMAT (1X,A5 ,5x,a12 ,4(5x,a10))
104 FORMAT (1X,I5 ,5x,f12 .1 ,4(5x,f10 .3))
105 FORMAT (1X ,5(a5 ,5x))
106 FORMAT (1X ,5(i5 ,5x))
return
end subroutine geomet



! subroutine scode 
Subroutine scode 
Use globalvar
Integer ::i,j,nres ,nlink ,idnode , jdir 
Integer :: nmast , nslave
Allocate ( idof (nnode ,2)) 
! initialization idof
idof (: ,:) =0
! reading / writing restraints 
Read (11 ,*) rec
Read (11 ,*) nres
Write (12 , '(/ ,1x ," NUMBER OF RESTRAINTS =" ,1x,i3)') nres 
Write (12 ,101) 'node ','direction '
Do i=1, nres
Read (11 ,*) idnode , jdir 
Write (12 ,102) idnode , jdir 
Idof (idnode , jdir )=-1
End do
! links
Read (11 ,*) rec 
Read (11 ,*) nlink 
If(nlink <0) then
Write (12 , '(1x ," NUMBER OF LINKS IS INCORRECT ") ') 
Stop
Return 
End if
Write (12 , '(/ ,1x ," NUMBER OF LINKS =" ,1x,i3)') nlink 
If(nlink >0) then
Write (12 ,103) 'MASTER ','SLAVE ','DIRECTION ' 
Do i=1, nlink
Read (11 ,*) nmast ,nslave , jdir 
Write (12 ,104) nmast ,nslave , jdir 
Idof (nslave , jdir )= nmast
End do 
End if
!dof
Ndof =0
Write (12 , '(/ ,1x ," DOF MATRIX ") ') 
Write (12 , '(4 A8)')'NODE ','UX ','UY ' 
Ext :do i=1, nnode
Int :do j=1 ,2
If( idof (i,j) <0) cycle int ! restraint 
If( idof (i,j) ==0) then !dof
ndof = ndof +1 
Idof (i,j)= ndof
Else ! link
Idnode = idof (i,j)
Idof (i,j)= idof (idnode ,j)
End if 
End do int
! writing dof matrix
Write (12 , '(3(4x,i4))')i ,( idof (i ,:)) 
End do ext
write (* ,*) 'end scode '
Return !end of the subroutine goes back to main program 
! format list
101 FORMAT (1X,A4 ,T12 ,A9)
102 FORMAT (1X,I4 ,T18 ,I3)
103 FORMAT (1X,A6 ,T16 ,A5 ,T26 ,A9)
104 FORMAT (1X,I4 ,T18 ,I3 ,T32 ,I3)
End subroutine scode



! subroutine loads 
subroutine loads 
use globalvar
integer ::i,nnl ,idnode ,jdir ,jdof ,nel ,idele , npre
real(4) ::b(3) ,c(3)
write (12 , '(/ ," NODAL AND ELEMENT LOADS ") ') 
allocate ( vloads ( ndof ))
vloads (:) =0.
! reading / writing of the nodal loads 
read (11 ,*) rec
read (11 ,*) nnl 
If(nnl <0) then
Write (12 , '(1x ," incorrect number of nodal loads  ") ') 
Stop
Return 
End if
write (12 , '(/ ,1x ," NUMBER OF NODAL LOADS =" ,1x,i3)')nnl 
If(nnl >0) then
Write (12 ,101) 'NODE ','DIR ','FORCE ' 
Do i=1, nnl
Read (11 ,*) idnode ,jdir , pforce
If (( idnode <1) .or .( idnode > nnode )) then
write (12 , '(1x ," nodal load error") ') 
stop
return 
End if
If (( jdir <1) .or .( jdir >3) ) then
write (12 , '(1x ," the direction of nodal load is incorrect ") ') 
stop
return 
End if
Write (12 ,102) idnode ,jdir , pforce 
jdof = idof (idnode , jdir )
vloads ( jdof )= vloads ( jdof )+ pforce
End do
Write (12 , '(/ ,1x ," VECTOR OF NODAL LOADS ") ') 
Write (12 ,103) 'DOF ','FORCE '
Do i=1, ndof
Write (12 ,104)i, vloads (i) 
End do
End if
! Surface Loads 
read (11 ,*) rec 
read (11 ,*) npre 
if(npre <0) then
write (12 , '(/ ,1x ," Number of surface loads incorrect ") ') 
stop
return 
end if
write (12 , '(/ ,1x ," NUMBER OF SURFACE LOADS =" ,1x,i3)') npre 
if(npre >0) then
write (12 ,105) 'N1 ','N2 ','PX1 ','PX2 ','PY1 ','PY2 ' 
do i=1, npre
read (11 ,*) n1 ,n2 ,px1 ,px2 ,py1 ,py2 
write (12 ,106) n1 ,n2 ,px1 ,px2 ,py1 ,py2 
dx = coord (n2 ,1) - coord (n1 ,1)
dy = coord (n2 ,2) - coord (n1 ,2) 
al = sqrt (dx **2. + dy **2.)
fx1 = (2.* px1 +px2 )*al /6. 
fx2 = (2.* px2 +px1 )*al /6. 
fy1 = (2.* py1 +py2 )*al /6. 
fy2 = (2.* py2 +py1 )*al /6.
j1x = idof (n1 ,1) 
j2x = idof (n2 ,1) 
j1y = idof (n1 ,2) 
j2y = idof (n2 ,2) 
if(j1x >0) then
vloads ( j1x) = vloads (j1x )+fx1
end if
if(j2x >0) then
vloads ( j2x) = vloads (j2x )+fx2
end if
if(j1y >0) then
vloads ( j1y) = vloads (j1y )+fy1
end if
if(j2y >0) then
vloads ( j2y) = vloads (j2y )+fy2
end if 
end do
end if
! Volume Loads
do i = 1, nele
it = itype (i) 
weight = ctype (it ,4) 
thick = ctype (it ,5) 
n1 = in(i ,1)
n2 = in(i ,2)
n3 = in(i ,3)
x1 = coord (n1 ,1)
x2 = coord (n2 ,1)
x3 = coord (n3 ,1)
y1 = coord (n1 ,2)
y2 = coord (n2 ,2)
y3 = coord (n3 ,2)
a1 = x2*y3 -x3*y2
a2 = x3*y1 -x1*y3
a3 = x1*y2 -x2*y1 
area = (a1+a2+a3)/2.
force = -weight * thick * area /3. 
do j = 1,3
jdof = idof (in(i,j) ,2) 
if (jdof >0) then
vloads ( jdof ) = vloads ( jdof )+ force
end if 
end do
end do
! Thermal Loads 
read (11 ,*) rec 
read (11 ,*) ntherm 
if( ntherm <0) then
write (12 , '(/ ,1x ," Number of thermal loads incorrect ") ') 
stop
return 
end if
write (12 , '(/ ,1x ," NUMBER OF THERMAL LOADS =" ,1x,i3)') ntherm 
allocate ( therm ( nele ))
therm (:) =0.
if( ntherm >0) then
do i=1, ntherm
read (11 ,*) idele , ther
therm ( idele ) = ther
end do
do i = 1, nele
n1 = in(i ,1)
n2 = in(i ,2)
n3 = in(i ,3)
x1 = coord (n1 ,1) 
x2 = coord (n2 ,1) 
x3 = coord (n3 ,1) 
y1 = coord (n1 ,2) 
y2 = coord (n2 ,2) 
y3 = coord (n3 ,2) 
b(1) = y2 -y3
b(2) = y3 -y1
b(3) = y1 -y2
c(1) = -(x2 -x3)
c(2) = -(x3 -x1)
c(3) = -(x1 -x2)
it = itype (i)
ee = ctype (it ,1) 
pois = ctype (it ,2) 
alpha = ctype (it ,3) 
if ( icode ==0) then
ck = 1.- pois
else if ( icode ==1) then
ck = (1.+ pois ) *(1. -2.* pois )
end if
coeff = ee* alpha * therm (i)* thick /(2.* ck) 
do j = 1,3
fx = coeff *b(j)
fy = coeff *c(j)
jx = idof (in(i,j) ,1)
jy = idof (in(i,j) ,2) 
vloads (jx) = vloads (jx)+fx 
vloads (jy) = vloads (jy)+fy
end do 
end do
end if
write (* ,*) 'end loads ' 
return
! format list
101 FORMAT (1X,A4 ,3X,A4 ,3X,A20)
102 FORMAT (1X,I4 ,3X,I4 ,3X,F20 .10)
103 FORMAT (1X,A4 ,10X,A12)
104 FORMAT (1X,I4 ,10X,F12 .3)
105 FORMAT (1X ,2(a7 ,5X) ,4(a11 ,5X))
106 FORMAT (1X ,2(I7 ,5X) ,4( F11 .3 ,5X))
End subroutine loads 



! subroutine allocating
subroutine allocating 
use globalvar
allocate (vk(ndof , ndof ))
allocate ( vdisp ( ndof ))
allocate ( signod (nnode ,3))
allocate ( epsnod (nnode ,3))
allocate ( sumweight ( nnode ))
allocate ( sigvm ( nele ))
allocate ( ievol ( nele )) 
ievol (:) =1
allocate (vkk( nsky ))
write (* ,*) 'end allocating ' 
return
end subroutine allocating 



! subroutine assemb
subroutine assemb 
use globalvar 
integer ::i,j, ncode (6) 
vk (: ,:) =0.
do ne =1, nele 
n1=in(ne ,1) 
n2=in(ne ,2) 
n3=in(ne ,3)
ncode (1:2) = idof (n1 ,1:2)
ncode (3:4) = idof (n2 ,1:2)
ncode (5:6) = idof (n3 ,1:2) 
if ( icode /=2) then
call mkk2d (ne) 
else
call mkk3d (ne) 
end if
ext : do i=1 ,6
ic= ncode (i)
if(ic <0) cycle ext 
int : do j=1 ,6
jc= ncode (j)
if(jc <0) cycle int
vk(ic ,jc)=vk(ic ,jc)+st(i,j)
end do int 
end do ext
end do 
return
end subroutine assemb 


! subroutine mkk2d
subroutine mkk2d (ne) 
use globalvar 
integer ::n1 ,n2 ,n3 ,it
real(4) ::aa ,ee ,pois ,tmp (3 ,6) ,b(3 ,6) ,d(3 ,3) ,x1 ,x2 ,x3 ,y1 ,y2 ,y3 ,a1 ,a2 ,a3 ,b1 ,b2 ,b3 ,c1 ,c2 ,c3
! geometric properties of the element 
n1=in(ne ,1) ! from the incident matrix 
n2=in(ne ,2)
n3=in(ne ,3)
x1= coord (n1 ,1)  ! node coordinates 
x2= coord (n2 ,1)
x3= coord (n3 ,1) 
y1= coord (n1 ,2) 
y2= coord (n2 ,2)
y3= coord (n3 ,2)
! inner compatibility matrix
a1=x2*y3 -x3*y2 
a2=x3*y1 -x1*y3 
a3=x1*y2 -x2*y1 
aa =( a1+a2+a3) /2. 
b1=y2 -y3
b2=y3 -y1 
b3=y1 -y2 
c1 =-(x2 -x3) 
c2 =-(x3 -x1) 
c3 =-(x1 -x2)
! computation of B
b(: ,:) =0. ! initialization of the compatibility matrix
b(1 ,1)=b1 !for constanst strain triangle
b(2 ,2)=c1
b(3 ,1)=c1
b(3 ,2)=b1
b(1 ,3)=b2
b(2 ,4)=c2
b(3 ,3)=c2
b(3 ,4)=b2
b(1 ,5)=b3
b(2 ,6)=c3
b(3 ,5)=c3
b(3 ,6)=b3
b(: ,:)=b(: ,:) /(2.* aa)
! save in tmeporary file 
write (9 ,*)b(: ,:)
! characteristics of the element 
it= itype (ne)
thick = ctype (it ,5) ! thickness of cst 
ee= ctype (it ,1)  !Young 's modulus 
pois = ctype (it ,2)  ! Poisson 's ration 
! accounting for evolution
if( ievol (ne) ==0) then
ee= efactor *ee
end if
! computation of D
d(: ,:) =0.
if( icode ==1) then
ee=ee /(1. - pois **2.) 
pois = pois /(1. - pois )
end if
d(1 ,1)=ee /(1. - pois **2.)
d(1 ,2)=ee* pois /(1. - pois **2.) 
d(2 ,1)=d(1 ,2)
d(2 ,2)=d(1 ,1)
d(3 ,3)=ee /(2.*(1.+ pois ))
! save in temporary file 
write (9 ,*)d(: ,:)
! initialization
tmp (: ,:) =0. 
st (: ,:) =0.
! computation of K 
tmp = matmul (d,b)
st= matmul ( transpose (b),tmp ) 
st= thick *aa*st
! save in temporary file 
write (9 ,*) st (: ,:)
! print garbage
if( check ==0) then  ! check for printing of evolved structure
write (99 , '(/ ,1x ," ELEMENT NUMBER : ",i5)')ne
write (99 , '(/ ,1x ,"b VALUES ") ')
write (99 , '(3(2x,f12 .2))')b1 ,b2 ,b3
write (99 , '(/ ,1x ,"c VALUES ") ')
write (99 , '(3(2x,f12 .2))')c1 ,c2 ,c3
write (99 , '(/ ,1x ,"D MATRIX ") ') 
do i=1 ,3
write (99 , '(3(2x,f12 .2))')d(i ,:) 
end do
! print stiffness matrix
write (99 , '(/ ,1x ," STIFFNESS MATRIX ") ') 
do i=1 ,6
write (99 , '(6(2x,f12 .2))')st(i ,:) 
end do
end if 
return
end subroutine mkk2d 

! subroutine mkk3d
subroutine mkk3d (ne)  !for the axisymmetric problem 
use globalvar
integer ::n1 ,n2 ,n3 ,it
real(4) ::aa ,ee ,pois ,tmp (4 ,6) ,b(4 ,6) ,d(4 ,4) ,x1 ,x2 ,x3 ,y1 ,y2 ,y3 ,a1 ,a2 ,a3 ,b1 ,b2 ,b3 ,c1 ,c2 ,c3 ,zm ,fn1 ,fn2 ,fn3 , coeffd
! geometric properties of the element 
n1=in(ne ,1) ! from the incident matrix 
n2=in(ne ,2)
n3=in(ne ,3)
x1= coord (n1 ,1)  ! node coordinates 
x2= coord (n2 ,1)
x3= coord (n3 ,1) 
y1= coord (n1 ,2) 
y2= coord (n2 ,2) 
y3= coord (n3 ,2)
! inner compatibility matrix
a1=x2*y3 -x3*y2 
a2=x3*y1 -x1*y3 
a3=x1*y2 -x2*y1 
aa =( a1+a2+a3) /2. 
b1=y2 -y3
b2=y3 -y1 
b3=y1 -y2 
c1 =-(x2 -x3) 
c2 =-(x3 -x1) 
c3 =-(x1 -x2)
! characteristics of the center of gravity 
thick =( x1+x2+x3)/3. ! mean radius
zm =( y1+y2+y3) /3.  ! mean height 
fn1 =a1+b1* thick +c1*zm
fn2 =a2+b2* thick +c2*zm 
fm3 =a3+b3* thick +c3*zm
! computation of b
b(: ,:) =0. 
b(1 ,1)=b1
b(1 ,3)=b2
b(1 ,5)=b3
b(2 ,2)=c1
b(2 ,4)=c2
b(2 ,6)=c3
b(3 ,1)=fn1/ thick
b(3 ,3)=fn2/ thick
b(3 ,5)=fn3/ thick
b(4 ,1)=c1
b(4 ,2)=b1
b(4 ,3)=c2
b(4 ,4)=b2
b(4 ,5)=c3
b(4 ,6)=b3
b(: ,:)=b(: ,:) /(2.* aa)
! save in tmeporary file 
write (9 ,*)b(: ,:)
! characteristics of the element 
it= itype (ne)
ee= ctype (it ,1)  !Young 's modulus 
pois = ctype (it ,2)  ! Poisson 's ration 
! computation of D
d(: ,:) =0.
coeffd =ee *(1. - pois ) /((1.+ pois ) *(1. -2* pois ))
d(1 ,1)= coeffd
d(1 ,2)= coeffd * pois /(1. - pois )
d(1 ,3)=d(1 ,2)
d(2 ,1)=d(1 ,2)
d(2 ,2)=d(1 ,1)
d(2 ,3)=d(1 ,3)
d(3 ,1)=d(1 ,3)
d(3 ,2)=d(2 ,3)
d(3 ,3)=d(1 ,1)
d(4 ,4)= coeffd *(1 -2.* pois ) /(2.*(1 - pois ))
! save in temporary file 
write (9 ,*)d(: ,:)
! initialization
tmp (: ,:) =0. 
st (: ,:) =0.
! computation of K 
tmp = matmul (d,b)
st= matmul ( transpose (b),tmp ) 
st= thick *aa*st
! save in temporary file 
write (9 ,*) st (: ,:) 
return
end subroutine mkk3d




! subroutine joints 
subroutine joints 
use globalvar
integer ::i, jdir
integer :: nelr ,nell ,ntdis ,nrdis , nincl
integer :: idnode ,n1 ,n2 
real(4) :: stiff ,adis , vkjj
write (12 , '(/ ," ELASTIC RESTRAINTS /LINKS - IMPOSED DISPLACEMENTS ") ') 
! elastic restraints
read (11 ,*) rec 
read (11 ,*) nelr 
if(nelr <0) then
write (12 , '(1x ," NUMBER OF ELASTIC RESTRAINTS INCORRECT ") ') 
stop
return 
end if
write (12 , '(1x ," NUMBER OF ELASTIC RESTRAINTS =" ,1x,i3)') nelr 
if(nelr >0) then
write (12 ,101) 'node ','direction ','stiffness ' 
do i = 1, nelr
read (11 ,*) idnode ,jdir , stiff 
write (12 ,102) idnode ,jdir , stiff 
jdof = idof (idnode , jdir )
vk(jdof , jdof ) = vk(jdof , jdof )+ stiff
end do 
end if
! elastic links 
read (11 ,*) rec 
read (11 ,*) nell 
if(nell <0) then
write (12 , '(1x ," number of elastic links incorrect ") ') 
stop
return 
end if
write (12 , '(1x ," NUMBER OF ELASTIC LINKS =" ,1x,i3)') nell 
if(nell >0) then
write (12 ,103) 'node1 ','node2 ','direction ','stiffness ' 
do i = 1, nell
read (11 ,*) n1 ,n2 ,jdir , stiff 
write (12 ,104) n1 ,n2 ,jdir , stiff 
jdof1 = idof (n1 , jdir )
jdof2 = idof (n2 , jdir )
vk(jdof1 , jdof1 ) = vk(jdof1 , jdof1 )+ stiff 
vk(jdof2 , jdof2 ) = vk(jdof2 , jdof2 )+ stiff 
vk(jdof1 , jdof2 ) = vk(jdof1 , jdof2 )-stiff 
vk(jdof2 , jdof1 ) = vk(jdof1 , jdof2 )
end do 
end if
! imposed absolute displacements 
read (11 ,*) rec
read (11 ,*) ntdis 
if(ntdis <0) then
write (12 , '(1x ," number of imposed absolute displacements incorrect ") ')
stop 
return
end if
write (12 , '(1x ," NUMBER OF IMPOSED ABSOLUTE DISPLACEMENTS =" ,1x,i3 )') ntdis
if(ntdis >0) then
write (12 ,105) 'node ','direction ','displacement ' 
do i = 1, ntdis
read (11 ,*) idnode ,jdir , adis 
write (12 ,106) idnode ,jdir , adis 
jdof = idof (idnode , jdir )
vkjj = vk(jdof , jdof ) 
vk(jdof ,1: ndof ) = 0.
vloads (1: ndof ) = vloads (1: ndof )-vk (1: ndof , jdof )* adis 
vk (1: ndof , jdof ) = 0.
vk(jdof , jdof ) = vkjj 
vloads ( jdof ) = vkjj * adis
end do 
end if
! imposed relative displacements 
read (11 ,*) rec
read (11 ,*) nrdis 
if(nrdis <0) then
write (12 , '(1x ," number of imposed relative displacements incorrect ") ')
stop 
return
end if
write (12 , '(1x ," NUMBER OF IMPOSED RELATIVE DISPLACEMENTS =" ,1x,i3)' ) nrdis
if(nrdis >0) then
write (12 ,107) 'node1 ','node2 ','direction ',' displacement ' 
do i = 1, nrdis
read (11 ,*) n1 ,n2 ,jdir , adis 
write (12 ,108) n1 ,n2 ,jdir , adis 
jdof1 = idof (n1 , jdir )
jdof2 = idof (n2 , jdir )
if (( jdof1 == -1).or .( jdof2 == -1)) then
write (12 , '(" displacements are applied to a constrained node !") ')
stop 
return
end if
vkjj = vk(jdof1 , jdof1 )+vk(jdof2 , jdof2 ) 
do k = 1, ndof
vk(jdof1 ,k) = vk(jdof1 ,k)+vk(jdof2 ,k)
end do
vloads ( jdof1 ) = vloads ( jdof1 )+ vloads ( jdof2 ) 
do k =1, ndof
vk(jdof2 ,k) = 0.
vk(k, jdof1 ) = vk(k, jdof1 )+vk(k, jdof2 ) 
vloads (k) = vloads (k)-vk(k, jdof2 )* adis 
vk(k, jdof2 ) = 0.
end do
vk(jdof1 , jdof1 ) = vk(jdof1 , jdof1 )+ vkjj 
vk(jdof1 , jdof2 ) = -vkjj
vk(jdof2 , jdof1 ) = -vkjj 
vk(jdof2 , jdof2 ) = vkjj
vloads ( jdof1 ) = vloads ( jdof1 )-vkjj * adis 
vloads ( jdof2 ) = vkjj * adis
end do 
end if

! application of restraints / imposed displacements 
if ( ityp ==0) then
! elastic restraints
vk(i,i)=vk(i,i)+ stiff
else
! imposed displacements
vkjj = vk(i,i) 
do k=1, ndof
vk(i,k)=0.
vloads (k)= vloads (k)-vk(k,i)* adis 
vk(k,i)=0.
end do
vk(i,i)= vkjj 
vloads (i)= vkjj * adis
end if 
end do
end if 
return
! format list
101 format (2x,a4 ,t15 ,a9 ,t30 ,a9)
102 format (3x,i3 ,t21 ,i3 ,t30 ,f9 .4)
103 format (1x,a5 ,t19 ,a5 ,t30 ,a9 ,t45 ,a9)
104 format (1x,i5 ,t19 ,i5 ,t30 ,i9 ,t45 ,f9 .4)
105 format (2x,a4 ,t15 ,a9 ,t27 ,a12)
106 format (2x,i4 ,t15 ,i9 ,t27 ,f12 .4)
107 format (1x,a5 ,t19 ,a5 ,t30 ,a9 ,t42 ,a12)
108 format (1x,i5 ,t19 ,i5 ,t30 ,i9 ,t42 ,f12 .4) 
end subroutine joints




! subroutine solve 
SUBROUTINE SOLVE
USE GLOBALVAR
INTEGER :: I,J,K,JDOF
REAL(4) :: DNOD(2), DMX(3), DMN(3), PVAL, TEMP(NDOF, NDOF) !INITIALIZATION OF THE VECTOR OF THE CONSTANT TERMS
VDISP = VLOADS 
!GAUSSIAN ELIMINATION METHOD DO I=1,NDOF
DO K=I+1,NDOF
VK(I,K) = VK(I,K)/VK(I,I)
END DO
VDISP(I) = VDISP(I)/VK(I,I) VK(I,I) = 1.
DO J=I+1,NDOF
DO K = I+1,NDOF
VK(J,K) = VK(J,K)-VK(J,I)*VK(I,K) END DO
VDISP(J) = VDISP(J)-VK(J,I)*VDISP(I) VK(J,I) = 0.
END DO 
END DO
!BACKWARD SUBSTITUTION
DO I=NDOF,1,-1
DO K=I+1,NDOF
VDISP(I) = VDISP(I)-VK(I,K)*VDISP(K)
END DO
END DO
!writing nodal displacements
if(check==0)then !check for printing of evolved structure
write(12,’(/,/,"SUMMARY OF THE RESULTS")’)
write(12,’(/,1x,"NODAL DISPLACEMENTS")’)
write(12,102)’NODE’,’Ux’,’Uy’
end if
dnod=0. !initialization vector of nodal displacements
dmx=-1.e6 !initialization of highest value
dmn=1.e6 !inizialization of lowest value
do i=1,nnode do j=1,2
if(idof(i,j)<=0)then !indentifying restraints 
dnod(j)=0. 
else
jdof=idof(i,j)
dnod(j)=vdisp(jdof) 
 end if
pval=dmx(j) 
dmx(j)=max(pval,dnod(j)) 
pval=dmn(j) 
dmn(j)=min(pval,dnod(j)) 
end do
if(check==0)then 
write(12,103)i,(dnod(j),j=1,2) end if
end do
if(check==0)then 
write(12,104)’MAX’,(dmx(j),j=1,2) !max
write(12,104)’MIN’,(dmn(j),j=1,2) !min
write(12,’(/,1x,"DISPLACEMENT AT DEGREE OF FREEDOM")’) 
write(12,105)’DOF’,’DISPLACEMENT’
do i=1,ndof
write(12,106)i,vdisp(i) 
end do
end if
return
!format list
102 format(2x,a5,t15,2(a20,2x))
103 format(2x,i4,t15,2(f15.4,2x)) 
104 format(2x,a5,t15,2(f15.4,2x)) 
105 format(2x,a5,t15,a20)
106 format(2x,i4,t15,f15.4)
end subroutine solve





! subroutine stress2d 
subroutine stress2d 
use globalvar 
integer ::i,j,k, ncode (6)
real(4) :: pigr =3.14159265359 , qs (6) ,b(3 ,6) ,d(3 ,3) ,ed1el ( nele ),tmp (6) ,theta
! rewind 
rewind (9)
! writing the type of weight functions (0= arithmetic average , 1= mean on volume , 2= energy based weight )
if( check ==0) then  ! check for printing of evolved structure 
write (12 , '(/ ,1x ," NODAL STRESS / STRAINS BASE ON ") ')
if( imean ==0) then
write (12 , '(1x ," ARITHMETIC MEAN ") ') 
else if( imean ==1) then
write (12 , '(1x ," MEAN ON VOLUME ") ') 
else if( imean ==2) then
write (12 , '(1x ," MEAN ON THE ENERGY ") ') 
else
write (12 , '(1x ," INCORRECT WEIGHT PARAMETER SPECIFIED ") ') 
end if
end if
! initialization
edtot =0. 
edtotw =0. 
volumetot =0. 
edptot =0.
sigvm (:) =0. !von mises stresses
! evaluation of the deformation energy from the external work 
edtotw = dot_product (vdisp , vloads )
edtotw = edtotw /2.
! deformation and stresses of each element 
do ne =1, nele
read (9 ,*)b
read (9 ,*)d 
read (9 ,*) st 
n1=in(ne ,1) 
n2=in(ne ,2) 
n3=in(ne ,3)
ncode (1:2) = idof (n1 ,:)
ncode (3:4) = idof (n2 ,:)
ncode (5:6) = idof (n3 ,:) 
do i=1 ,6
nndof = ncode (i) 
if( nndof == -1) then
qs(i)=0.
else
qs(i)= vdisp ( nndof )
end if 
end do
do i=1 ,3 
eps (i)=0. 
eps0 (i)=0. 
do j=1 ,6
eps (i)=eps (i)+b(i,j)*qs(j)
end do 
end do
it= itype (ne)
ee= ctype (it ,1) 
pois = ctype (it ,2) 
alpha = ctype (it ,3) 
if( icode ==0) then
eps0 (1) = alpha * therm (ne) 
eps0 (2) = eps0 (1)
else
eps0 (1) =(1. - pois )* alpha * therm (ne) 
eps0 (2) = eps0 (1)
end if
do i=1 ,3 
sig (i)=0. 
do j=1 ,3
sig (i)=sig (i)+d(i,j)*( eps(j)-eps0 (j))
end do 
end do
! dependent stress / strain component 
if( icode ==0) then
sig (4) =0.
eps (4) =- pois /(1. - pois )*( eps (1)+eps (2)) +(1.+ pois )/(1. - pois )* alpha * therm (ne)
else
eps (4) =0.
sig (4) = pois *( sig (1) +sig (2))-ee* alpha * therm (ne) 
end if
! principal strain / stresses
c1 =( sig (1)+ sig (2)) /2.
c2= dsqrt ((( sig (1) -sig (2) ) /2.) **2.+ sig (3) **2.)
sigp (1) =c1+c2
sigp (2) =c1 -c2
sigvm (ne)= dsqrt ( sigp (1) **2.+ sigp (2) **2. - sigp (1)* sigp (2))  !Von Mises stresses
c1 =( eps (1)+ eps (2)) /2.
c2= dsqrt ((( eps (1) -eps (2) ) /2.) **2.+ eps (3) **2.)
epsp (1) =c1+c2 
epsp (2) =c1 -c2
if( sig (1) /= sig (2)) then
theta = datan (2.* sig (3) /( sig (1) -sig (2)))/2.
else
theta = pigr /4. 
if(sig (3) <0) then
theta =- theta
end if 
end if
if( sig (1) <sig (2)) then 
theta = theta + pigr /2.
else
if(sig (3) <0) then 
theta = theta + pigr
end if 
end if
theta = theta *180./ pigr
! write (12 , '(/ ,1x ," DEFORMATION IN ELEMENT N:" ,1x,i5) ')ne
! write (12 ,101) 'EPSX ','EPSY ',' GAMMAXY ','EPSZ ','EPS1 ','EPS2 ','THETA'
! write (12 ,102) (eps(i),i=1 ,4) ,( epsp (i),i=1 ,2) ,theta

! write (12 ,101) 'SIGX ','SIGY ','TAUXY ','SIGZ ','SIG1 ','SIG2 ','THETA ' 
! write (12 ,102) (sig(i),i=1 ,4) ,( sigp (i),i=1 ,2) ,theta
! weights
it= itype (ne) 
thick = ctype (it ,5)
x1= coord (n1 ,1)  ! node coordinates 
x2= coord (n2 ,1)
x3= coord (n3 ,1) 
y1= coord (n1 ,2) 
y2= coord (n2 ,2) 
y3= coord (n3 ,2) 
a1=x2*y3 -x3*y2 
a2=x3*y1 -x1*y3 
a3=x1*y2 -x2*y1 
aa =( a1+a2+a3) /2. 
volume =aa* thick
volumetot = volumetot + volume
! start of evolve stuff 
if( ievol (ne) ==1) then
volk = volk + volume
end if
!end of evolve stuff 
ed1el (ne)=0.
tmp = matmul (st ,qs)
ed1el (ne)= dot_product (qs ,tmp)/2. 
ed1elp = ed1el (ne) *100./ edtotw 
edptot = edptot + ed1elp
! write (12 , '(/ ,1x ," ENERGY IN THE ELEMENT N.:" ,1x,i5) ')ne 
! write (12 ,103) 'VOLUME ','ENERGY ',' ENERGY %'
! write (12 ,104) volume , ed1el (ne),ed1elp 
! nodal deformation / stress
if( imean ==0) then  ! arithmetic mean
wmed =1.
else if( imean ==1) then
wmed = volume ! volume mean
else
wmed = ed1el (ne)  ! energy mean 
end if
do i=1 ,3
node =in(ne ,i)
sumweight ( node )= sumweight ( node )+ wmed 
do j=1 ,3
epsnod (node ,j)= epsnod (node ,j)+eps (j)* wmed 
signod (node ,j)= signod (node ,j)+sig (j)* wmed
end do 
end do
end do
do i=1, nnode 
do j=1 ,3
epsnod (i,j)= epsnod (i,j)/ sumweight (i) 
signod (i,j)= signod (i,j)/ sumweight (i)
end do 
end do
! writing many informations
if( check ==0) then  ! check for printing of evolved structure
! WRITING OF NODAL DEFORMATIONS
write (12 , '(/ ,1x ," MEAN NODAL DEFORMATIONS ") ') 
write (12 ,105) 'NODE ','EPSX ','EPSY ','GAMMAXY ' 
do k=1, nnode
write (12 ,106)k ,( epsnod (k,i),i=1 ,3) 
end do
! WRITING OF NODAL STRESSES
write (12 , '(/ ,1x ," MEAN NODAL STRESSES ") ') 
write (12 ,105) 'NODE ','SIGX ','SIGY ','TAUXY ' 
do k=1, nnode
write (12 ,107)k ,( signod (k,i),i=1 ,3) 
end do
! WRITING OF ENERGY STORED IN THE STRUCTURE
write (12 , '(/ ,1x ," ENERGY STORED IN THE STRUCTURE ") ')
write (12 ,108) 'TOTAL VOLUME ','TOTAL STRAIN ENERGY ','TOTAL ENERGY %'
write (12 ,109) volumetot ,edtot , edptot
! WRITING OF WORK DONE BY EXTERNAL FORCES
write (12 , '(/ ,1x ," WORK DONE BY THE EXTERNAL FORCE :" ,1x,f20 .10) ')edtotw
end if 
return
! format list
101 format (2x ,7(a7 ,5x))
102 format (2x ,6( f10 .3 ,2x),f7 .2)
103 format (2x ,3(a7 ,5x))
104 format (2x ,3( f7 .2 ,5x))
105 format (2x,a4 ,5x ,3( a20 ,5x))
106 format (2x,i4 ,5x ,3( f20 .10 ,5x))
107 format (2x,i4 ,5x ,3( f20 .5 ,5x))
108 format (2x,a20 ,5x,a20 ,5x,a20)
109 format (2x,f20 .10 ,5x,f20 .10 ,5x,f20 .10)
end subroutine stress2d




! subroutine evolve
subroutine evolve 
use globalvar
! defining variables 
integer ::kk ,i,j
! reading parameters from input 
if( check ==0) then
read (11 ,*) rec 
read (11 ,*) check 
if(check >0) then
read (11 ,*) efactor ,rr0 ,er ,alphan ,rrmax ,rmvmax , itermax
write (12 ,101) 'EFACTOR ','RR0 ','ER ','ALPHAN ','RRMAX ','RMVMAX ',' MAXITER '
write (12 ,102) efactor ,rr0 ,er ,alphan ,rrmax ,rmvmax , itermax 
end if
return 
end if
! initializing variables 
if( check ==1) then
v0= volumetot  ! initial volume
volk =v0 ! volume at k loop
nelek = nele  ! elements at k loop
rr= rr0  ! rejection ratio at k loop 
rmv =0.
iconv =0 ! convergence indicator 
end if
! main evolution loop 
itloop :do kk =1, itermax
write (*, '(" entering loop number : ",i3)')kk 
! element removal
nrk =0
sigmalim =rr* maxval ( sigvm ) 
ext :do ne =1, nele
if( ievol (ne) ==0) cycle ext ! element already eliminated 
if( sigvm (ne) >= sigmalim ) cycle ext  ! element to keep 
! write (*,'(i3 ,5x,f20 .7) ')ne , sigvm (ne)
ievol (ne)=0 
nrk =nrk +1
end do ext 
nelek =nelek -nrk
nrmin =int( alphan * nelek )
rmv =1. - volk /v0  ! removed volume
write (*, '(" removed elements = ",i3 ," | number of total elements = ",i4 ,"| updated volume = ",f14 .10 ," | rr = ",f5 .3 ," | rmv = ",f5 .3) ')nrk ,nelek ,volk ,rr ,rmv
! solution of the new structure
volk =0.
open (9, file ='mkk .txt ',status ='unknown ') 
call assemb
call solsky 
! call solve 
call stress2d
if(nrk <= nrmin )rr=rr+er  
iconv =1
write (*, '(/ ,1x ," MAXIMUM RR REACHED ") ') 
end if
if(rmv > rmvmax ) then  
iconv =2
write (*, '(/ ,1x ," MAXIMUM RMV REACHED ") ') 
end if
if (( iconv ==1) .or .( iconv ==2) ) then 
call plot2d
close (9) 
exit itloop
end if 
close (9)
end do itloop
! writing of nodal deformations
write (12 , '(/ ,1x ," MEAN NODAL DEFORMATIONS ") ') 
write (12 ,105) 'NODE ','EPSX ','EPSY ','GAMMAXY ' 
do k=1, nnode
write (12 ,106)k ,( epsnod (k,i),i=1 ,3) 
end do
! writing of nodal stresses
write (12 , '(/ ,1x ," MEAN NODAL STRESSES ") ') 
write (12 ,105) 'NODE ','SIGX ','SIGY ','TAUXY ' 
do k=1, nnode
write (12 ,107)k ,( signod (k,i),i=1 ,3) 
end do
! writing of energy stored in the structure
write (12 , '(/ ,1x ," ENERGY STORED IN THE STRUCTURE ") ')
write (12 ,108) 'TOTAL VOLUME ','TOTAL STRAIN ENERGY ','TOTAL ENERGY %'
write (12 ,109) volumetot ,edtot , edptot
! writing of work done by external forces
write (12 , '(/ ,1x ," WORK DONE BY THE EXTERNAL FORCE :" ,1x,f20 .10) ')edtotw
write (12 , '(" NEW SET OF ELEMENTS ") ') 
! writing element data
j=0
write (12 , '(a8 ,2x,a20)')'ELEMENT ','SIGVM ' 
do i=1, nele
if( ievol (i) ==1) then
write (12 , '(i8 ,2x,f20 .10) ')i, sigvm (i) 
j=j+1
end if 
end do
write (12 , '(" TOTAL NUMBER OF ELEMENTS =" ,1x,i3)')j 
! format list
101 format (1x,a10 ,2x ,6(a7 ,2x))
102 format (1x,e10 .2 ,2x ,5( f7 .3 ,2x),i7)
105 format (2x,a4 ,5x ,3( a20 ,5x))
106 format (2x,i4 ,5x ,3( f20 .10 ,5x))
107 format (2x,i4 ,5x ,3( f20 .5 ,5x))
108 format (2x,a20 ,5x,a20 ,5x,a20)
109 format (2x,f20 .10 ,5x,f20 .10 ,5x,f20 .10)
return
end subroutine evolve




! subroutine plot2d 
subroutine plot2d 
use globalvar
integer :: i,j, ncode (6) ,jdof
real(4) :: pigr =3.14159265359 , temp , tempp
real(4) :: qs (6) ,ndisp (2) ,b(3 ,6) ,d(3 ,3) ,ed1el ( nele ),tmp (6) 
! rewind of the temporary file
rewind (9)
! opening file " plot .dat" (for post - processing ) 
open (33 , file ='plot .dat ',status ='unknown ')
! type of analysis 
! 0 = plain stress 
! 1 = plain strain 
! 2 = axisymmetric 
write (33 ,*) icode
! number of nodes and coordinates 
write (33 ,*) nnode
do i=1, nnode
write (33 ,*) ( coord (i,j),j=1 ,2) 
end do
! number of elements and incidence 
write (33 ,*) nele
do i=1, nele
write (33 ,*) (in(i,j),j=1 ,3) 
end do
! strains and stresses in the elements 
do ne =1, nele
read (9 ,*) b
read (9 ,*) d
read (9 ,*) st 
n1 = in(ne ,1) 
n2 = in(ne ,2) 
n3 = in(ne ,3)
ncode (1:2) = idof (n1 ,:)
ncode (3:4) = idof (n2 ,:)
ncode (5:6) = idof (n3 ,:) 
do i=1 ,6
nndof = ncode (i) 
if ( nndof == -1) then
qs(i) = 0.
else
qs(i) = vdisp ( nndof )
end if
end do
do i=1 ,3
eps (i) = 0. 
eps0 (i) = 0. 
do j=1 ,6
eps (i) = eps (i)+b(i,j)*qs(j)
end do 
end do
it = itype (ne)
ee = ctype (it ,1) 
pois = ctype (it ,2) 
alpha = ctype (it ,3) 
if ( icode ==0) then
eps0 (1) = alpha * therm (ne) 
eps0 (2) = eps0 (1)
else
eps0 (1) = (1. - pois )* alpha * therm (ne) 
eps0 (2) = eps0 (1)
end if 
do i=1 ,3
sig (i) = 0. 
do j=1 ,3
sig (i) = sig (i)+d(i,j)*( eps(j)-eps0 (j))
end do 
end do
! dependent stresses / deformations 
if ( icode ==0) then
sig (4) = 0.
eps (4) = -pois /(1. - pois )*( eps (1)+eps (2)) +(1.+ pois )/(1. - pois )* alpha * therm (ne)
else
eps (4) = 0.
sig (4) = pois *( sig (1) +sig (2))-ee* alpha * therm (ne) 
end if
! principal stresses / deformations
c1 = ( sig (1)+ sig (2)) /2.
c2 = sqrt ((( sig (1) -sig (2) ) /2.) **2.+ sig (3) **2.) 
sigp (1) = c1+c2
sigp (2) = c1 -c2
c1 = ( eps (1)+ eps (2)) /2.
c2 = sqrt ((( eps (1) -eps (2) ) /2.) **2.+ eps (3) **2.) 
epsp (1) = c1+c2
epsp (2) = c1 -c2
if( sig (1) /= sig (2)) then
theta = atan (2.* sig (3) /( sig (1) -sig (2)))/2. 
else
theta = pigr /4. 
if (sig (3) <0) then
theta = -theta
end if 
end if
if( sig (1) <sig (2)) then 
theta = theta + pigr /2.
else
if (sig (3) <0) then 
theta = theta + pigr
end if 
end if
theta = theta *180./ pigr 
if(check >0) then
temp =sig (1) 
tempp =eps (1)
sig (1) = sigvm (ne)  ! displaying Von Mises stresses 
eps (1) =- ievol (ne) ! displaying active elements
endif
write (33 ,*) (eps(j),j=1 ,4) ,( epsp (j),j=1 ,2) 
write (33 ,*) (sig(j),j=1 ,4) ,( sigp (j),j=1 ,2) 
if(check >0) then
sig (1) = temp 
eps (1) = tempp
endif 
end do
! nodal displacements 
do i=1, nnode
do j=1 ,2
if ( idof (i,j) <=0) then 
ndisp (j)=0.
else
jdof = idof (i,j) 
ndisp (j) = vdisp ( jdof )
endif 
end do
write (33 ,*) ( ndisp (j),j=1 ,2) 
end do
edtot = 0. 
edtotw = 0. 
edptot = 0.
! evaluation of the external work 
edtotw = dot_product (vdisp , vloads ) 
edtotw = edtotw /2.
! evaluation of internal energy 
rewind (9)
do ne =1, nele
read (9 ,*) b
read (9 ,*) d
read (9 ,*) st 
n1 = in(ne ,1) 
n2 = in(ne ,2) 
n3 = in(ne ,3)
ncode (1:2) = idof (n1 ,:)
ncode (3:4) = idof (n2 ,:)
ncode (5:6) = idof (n3 ,:) 
do i=1 ,6
nndof = ncode (i) 
if ( nndof == -1) then
qs(i) = 0.
else
qs(i) = vdisp ( nndof )
endif 
end do
ed1el (ne) = 0.
tmp = matmul (st ,qs)
ed1el (ne) = dot_product (qs ,tmp) 
ed1el (ne) = ed1el (ne)/2.
ed1elp = ed1el (ne) *100./ edtotw 
edtot = edtot + ed1el (ne) 
edptot = edptot + ed1elp
write (33 ,*) ed1elp 
end do
write (33 ,*) edtot
! closing file plot .dat 
close (33)
return
end subroutine plot2d
