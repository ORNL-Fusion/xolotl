
!*********************************************************************************
MODULE mod1
integer :: loopsize, nm0, mobtruncate, tcut, nx, nc, mobmax, nmob, nsm, nmloop, j, i, &
	kk, mattop, matbott, matleft, matright, cap12nnz, cap34nnz, virtd34nnz, thickness, fluxtype, &
	iblankleft, iblankright

integer*8 :: maxnnz

logical :: maxnnzresized

real*8 :: bias, D00, D01, D0v, eb_min, r0, Em_loop, flux, rate, micros, tempr, &
        ceq(2), efi, efv, pi, Emi, Emv, dis_recomb, Em_loop2, iblankx
real*8, allocatable :: Ebin(:), Ebvn(:), mobil3(:), tt0(:), arraytemp(:), &
        mq(:), D0(:), Eb(:,:), r(:), x(:), dx(:), dx2(:), d(:), k(:,:), dis(:,:), zz1&
        (:), zz2(:,:), zz3(:,:), sum_array(:), arraytemp2(:), mattemp(:,:), &
	mattemp2(:,:), u(:,:), imp(:,:), casrhocoeff(:), casrho(:), casyield(:)

integer, allocatable :: mobil12(:,:), jt(:), hx(:), vy(:), mx(:), my(:),  mj(:), smx(:),&
        smy(:), smj(:), smj_in_mj(:), jt_in_mj(:,:), jt_in_smj(:,:), allowx(:), allowy(:),&
        allowj(:), allowj_nx(:,:), smj_in_allowj(:), old_allowj_in_new(:), capture(:,:),&
        capture1(:), capture3a(:), cap12(:), cap34(:), cap12ja(:), cap34ia(:), dissociation&
        (:,:), dissociation1(:,:), dissociation2(:,:), dissociation3a(:,:), dissociation3b&
        (:,:), dissociation4(:,:), dissoc12(:),dissoc34(:), virtualdis(:,:), virtualdis3a&
        (:), virtualdis4(:), virtd34(:), virtd34ia(:), intarraytemp(:), dbl_smj(:)&
        , coef1(:), coef2(:,:), jt_in_allowj(:), intarraytemp2(:), jt_in_allowj_mat1(:,:),&
	intmattemp(:,:)

logical, allocatable :: mask(:), maskmat(:,:), jt_in_allowj_mat2(:,:)

character*2 ef_system
character*2 eb_system

END MODULE mod1   

!********************************************************************************
MODULE modsolverinterface

integer :: ntspan, ntspan0, neq, old_neq, mm, ntt0, nh, nv, ni, nallow, old_nh, old_nv, old_ni, &
    old_nallow, num_threads
integer :: hr, minu, sec, hundsec   ! for timing purpose
integer*8 :: dfdy_nnz, npds
integer, allocatable :: dfdy_in(:), dfdy_jn(:), dfdy_ia(:), dfdy_ja(:), he_bound(:), &
	v_bound(:), i_bound(:)
real*8 ::  t0, rtol, atol, ILUdroptol
real*8, allocatable :: tspan(:), y0(:), dfdy_val(:), dfdy_csr(:)
character*1 set_htry
real*8::  htry

END MODULE modsolverinterface

!********************************************************************************
MODULE modfilehandles

integer fh2 /502/, fh3 /503/, fh4 /504/, fh5 /505/, fh6 /506/

END MODULE modfilehandles

!*********************************************************************************
MODULE mod2

contains

!************ define an indexing function *************
function indx(z1,z2,nz)
! ** z1: helium number; z2: V number; **
use mod1
use modsolverinterface
implicit none
integer :: z1(nz), z2(nz), indx(nz), nz
indx = z1*(ni+nv+1)+z2+ni+1
return
end function indx

function febi(n,nn)
implicit none
integer n(nn), nn
real*8 febi(nn)
!febi=4.33d0-5.76d0*(dble(n)**(2.d0/3.d0)-(dble(n)-1.d0)**(2.d0/3.d0))   ! for iron
febi=1.5d0*(4.33d0-5.76d0*(dble(n)**(2.d0/3.d0)-(dble(n)-1.d0)**(2.d0/3.d0)))
return
end function febi

function febv(n,nn)
integer n(nn), nn
real*8 febv(nn)
!febv=1.73d0-2.59d0*(dble(n)**(2.d0/3.d0)-(dble(n)-1.d0)**(2.d0/3.d0))   ! for iron
febv=3.d0-4.5625d0*(dble(n)**(2.d0/3.d0)-(dble(n)-1.d0)**(2.d0/3.d0))
return
end function febv

!************ a tool for finding array a in array b ***********
!%%%%%%%%%%% warning : this is only valid for cases w/o repeated elements in array b %%%%%%%%%
function xufind(na,a,nb,b)
! ** find a in b, return positions of a elements in b **
implicit none
integer :: na, nb, a(na), b(nb), xufind(na), i, j
do i=1,na
    xufind(i)=0
    do j=1,nb
        if (b(j)==a(i)) then
            xufind(i)=j
            exit
        end if
    end do
!    if (xufind(i)==0) print*,'Warning: a(',i,') was not found in b'
end do 
return
end function xufind

!************ a polynomial function for calculating depth-dependent defect generation probability ***********
function xupoly9(coeff,xdum,nxdum)
implicit none
integer :: nxdum, expn(10), aa
real*8 :: coeff(10), xdum(nxdum), xupoly9(nxdum)
expn=(/(aa,aa=9,0,-1)/)
xupoly9=0.d0
do aa=1,10
	xupoly9=xupoly9+coeff(aa)*xdum**expn(aa)
end do
return
end function xupoly9

END MODULE mod2    

!********************************************************************************
PROGRAM main
use mod1
use mod2
use modsolverinterface
use modfilehandles
use omp_lib
implicit none
real*8 tstart, tend, dclock
real*8, allocatable :: t_txt(:), y_txt(:,:)
logical logi(1)
external dclock, xusolver
pi=3.1415926535d0

!global tprevious k d dis u nh nv ni nallow old_nh old_nv old_ni old_nallow he_bound i_bound v_bound he_subbound i_subbound v_subbound ALLOWJ ABSTOL
print*,'set num_threads: '
read*, num_threads
print*,'set droptol for ILU linear tool (default: 0.01): '
read*, ILUdroptol
print*,'set foil thickness (60,84 or 108, integer, in nm): '
read*, thickness
print*,'set temperature: '
read*, tempr
print*,'set ion fluxtype (1: 1.56d-3 <=> 5d-4 dpa/sec, 2: 1.56d-4, 3: 1.56d-5): '
read*, fluxtype
print*,'set Em for I: '
read*, Emi
print*,'set Em for V: '
read*, Emv
print*,'set Em for defining intermediate loops:'
read*, Em_loop
print*,'set Em for large loops starting from I21:'
read*, Em_loop2
print*,'enter V-I recombination radius: '
read*, dis_recomb
print*,'set surface clearance depth for interstitials: '
read*, iblankx

call omp_set_num_threads(num_threads)

ef_system='AB'    !'AB' or 'MD', formation energy system
D00= 2.d11
D01=  D00
D0v= 2.d11

!%%%%%% input simple mobile species %%%%%%%
nsm=2        ! number of single mobile species, i.e., I, V, and/or He, either 2 or 3
allocate(smx(nsm),smy(nsm),smj(nsm),smj_in_mj(nsm))
smx=(/0,0/)
smy=(/-1,1/)

allocate(dbl_smj(nsm))

nm0= 21       ! number of small mobile species to be included in the mobil array; intersitital loops, if mobile, will be included later, separately.
allocate(mobil12(nm0,2),mobil3(nm0))
mobil12(:,1)=0       ! He-number
mobil12(:,2)= (/(kk,kk=-20,-1),(kk,kk=1,1)/)      ! V (or I) -number
mobil3(:11)= Em_loop-(Em_loop-Emi)/dble(11)*dble((/(kk,kk=1,11)/))  
mobil3(12:20)=Emi
mobil3(21:)=Emv      ! E_m migration energies of the small mobile species


!%%%%%% input time grids %%%%%%%%%
!ntspan=size((/0.d0,10.d0**(dble((/(kk,kk=-300,600,2)/))/100.d0)/))        ! check the length of requested tspan
if (fluxtype==1) then
	ntspan=size((/(kk,kk=0,10),(kk,kk=12,30,2),(kk,kk=60,3000,30)/))
	allocate(tspan(ntspan))
	!tspan=(/0.d0,10.d0**(dble((/(kk,kk=-300,600,2)/))/100.d0)/)
	tspan=dble((/(kk,kk=0,10),(kk,kk=12,30,2),(kk,kk=60,3000,30)/))
else if (fluxtype==2) then
	ntspan=size((/(kk,kk=0,10),(kk,kk=12,30,2),(kk,kk=35,300,5),(kk,kk=600,3000,300)/))
	allocate(tspan(ntspan))
	tspan=dble((/(kk,kk=0,10),(kk,kk=12,30,2),(kk,kk=35,300,5),(kk,kk=600,3000,300)/))
else if (fluxtype==3) then
	ntspan=size((/(kk,kk=0,10),(kk,kk=15,100,5),(kk,kk=150,3000,50),(kk,kk=6000,30000,3000)/))
	allocate(tspan(ntspan))
	tspan=dble((/(kk,kk=0,10),(kk,kk=15,100,5),(kk,kk=150,3000,50),(kk,kk=6000,30000,3000)/))
else
	print*,'fluxtype option not set correctly'
	stop
end if

ntspan0=ntspan

!%%%%%% input other parameters %%%%%%
!tempr= 80.d0     ! temperature in degree C
tempr=tempr+273.d0
nh=  0  !nh0; 
print*, 'input nv: '
read*, nv
print*,'input ni: '
read*, ni
print*,'input ntt0 (number of previously completed data points): '
read*, ntt0
print*,'enter smallest loop size: '
read*, loopsize
print*,'enter largest loop size (mobtruncate): '
read*, mobtruncate

tstart=dclock()

if (ntt0==0) then
	old_nh=nh
	old_nv=nv
	old_ni=ni
	set_htry='n'
else 
	print*,'enter old_nh:'
	read*, old_nh
	print*,'enter old_nv:'
	read*, old_nv
	print*,'enter old_ni:'
	read*, old_ni
	print*,'enter old_nallow:'
	read*, old_nallow
	print*,'set initial step size ? (''y'' for ''yes''; ''n'' for ''no''):'
	read*, set_htry
	if (set_htry=='y') then
		print*,'enter initial step size now: '
		read*, htry
	end if
end if

eb_min= -6.d0      ! cut-off value of binding energy, a negative value
bias=  1.2d0
r0=  0.d0   ! base for trapping radius
atol= 1.d-10   ! absolute tolerance
rtol= 1.d-3   ! relative tolerance
nx=19
allocate(x(nx),dx(nx-1),dx2(nx-2), arraytemp(13))
arraytemp=dble((/(kk,kk=0,thickness,thickness/12)/))   !12 is the unit thickness scale for the TEM fringes used in the Argonne experiments
x((/1,(kk,kk=5,15),19/))=arraytemp
x(2:4)=x(1)+x(5)/4.d0*dble((/(kk,kk=1,3)/))
x(16:18)=x(15)+x(5)/4.d0*dble((/(kk,kk=1,3)/))
deallocate(arraytemp)
print*,'x(1:6): ', x(1:6)
print*
print*,'x(14:19): ', x(14:)

dx=x(2:)-x(:nx-1)
dx2=(dx(2:)+dx(:nx-2))/2.d0
maxnnz=3e8   ! ** 1e6 is an initial estimate for the maximum no. of non-zeros *******

do kk=1,nx
	if (x(kk)<=iblankx) then 
		continue
	else
		exit
	end if
end do
iblankleft=kk-1

do kk=nx,1,-1
	if (x(nx)-x(kk)<=iblankx) then
		continue
	else 
		exit
	end if
end do
iblankright=kk+1
print*,'iblankleft= ',iblankleft
print*,'iblankright= ', iblankright


! open files for read/write
open(fh2,file='tdata',form='unformatted')
open(fh3,file='ydata',form='unformatted')
open(fh4,file='parameters.txt')
open(fh5,file='x_allowx_allowy',form='unformatted')
open(fh6,file='nt_completed.txt')


!rate=1.56d-3         ! in ions/ nm^2 /sec, 1.56d-3 corresponds to 5.d-4 dpa/sec.
if (fluxtype==1) then
	rate=1.56d-3
else if (fluxtype==2) then
	rate=1.56d-4
else if (fluxtype==3) then
	rate=1.56d-5
end if

allocate(imp(29,nx),casyield(29),casrho(nx),casrhocoeff(10))  ! cascade imp. for I4~I, and V~V2
	casyield=(/1.0314d-01, 0.d0, 0.d0, 0.d0, 2.0476d-01, 0.d0, 0.d0, 0.d0, 2.7050d-01, 0.d0, 0.d0, 6.0033d-01, 7.7354d-01, 1.8047d+00, 2.6454d+00, 4.0223d+00, 5.7084d+00, 8.7514d+00, 1.8810d+01, 1.2529d+02, 1.6741d+02, 2.4475d+01, 7.4273d+00, 4.4693d+00, 2.5481d+00, 0.d0, 0.d0, 0.d0, 1.2808d+00/)
do i=1,29
imp(i,:)=rate*casyield(i)/108.d0   !/x(nx)
end do

allocate(Ebin(int(5e6)),Ebvn(int(5e6)))
Ebin=febi((/(mm,mm=1,int(5e6))/),int(5e6))
Ebvn=febv((/(mm,mm=1,int(5e6))/),int(5e6))    


do while (.true.)

    ! open previous files
    
    if (nh<old_nh .or. nv<old_nv .or. ni<old_ni) then
        print*,'ERROR: New phase space does not cover old phase space completely at ntt0= ', ntt0
        stop
    end if
    if (ntt0 /= 0) then
	if (allocated(tt0)) deallocate(tt0)        
	allocate(tt0(ntt0))
	do mm=1,ntt0
	        read(fh2), tt0(mm)
	end do
        if (tt0(ntt0)==tspan(ntspan)) exit
        if (old_nh==nh .and. old_nv==nv .and. old_ni==ni .and. tt0(ntt0)/=tspan(1)) &
            print*,'New tspan does not agree with the end of previous session, automatic truncation will be performed'
        do tcut=1,ntspan
            if (tspan(tcut)==tt0(ntt0)) exit
        end do
        
        allocate(arraytemp(ntspan-tcut+1))
        arraytemp=tspan(tcut:)
        deallocate(tspan)
        allocate(tspan(ntspan-tcut+1))
        tspan=arraytemp
        ntspan=ntspan-tcut+1
        deallocate(arraytemp)
    end if
    t0=tspan(1)
    
    ! calculate different lengths
    nc=(nh+1)*(nv+ni+1)

    ! set total cluster composition index
    if (allocated (jt)) deallocate(jt,hx,vy)
    allocate(jt(nc),hx(nc),vy(nc))
    jt=(/(mm,mm=1,nc)/)
    hx=floor(dble(jt-1)/(ni+nv+1))  ! hx=He No. ;
    vy=jt-hx*(ni+nv+1)-ni-1      ! vy=V No.; if negative means I No.

    ! set mobile species
    mobmax=min(ni,mobtruncate)
    nmloop=max(mobmax-loopsize+1,0)
    nmob=nmloop+nm0
    if (allocated(mx)) deallocate(mx,my,mj,mq,D0)
    allocate(mx(nmob),my(nmob),mj(nmob),mq(nmob),D0(nmob))
    mx=(/(0,mm=1,nmloop), mobil12(:,1)/)
    my=(/(mm,mm=-mobmax,-loopsize), mobil12(:,2)/)
    mq=(/(Em_loop2,mm=1,nmloop), mobil3/)
    D0=(/(D01*dble((/(mm,mm=mobmax,loopsize,-1)/))**(-0.7d0)), (/(D00/dble(mm),mm=20,1,-1)/), (D0v,kk=1,1)/)
print*,'D0(nmloop+1:)= ', D0(nmloop+1:)
    mj=indx(mx,my,nmob)
    
    ! identify single mobile species, i.e., I, V and He, for the purpose of
    ! analyzing dissociation.
    smj=indx(smx,smy,nsm)
    smj_in_mj=xufind(nsm,smj,nmob,mj)
    
    if (allocated(jt_in_mj)) deallocate(jt_in_mj,jt_in_smj,jt_in_allowj)
    allocate(jt_in_mj(nc,2),jt_in_smj(nc,2),jt_in_allowj(nc),jt_in_allowj_mat1(-ni:nv,0:nh),&
	jt_in_allowj_mat2(-ni:nv,0:nh))
    jt_in_mj=0
    jt_in_smj=0
    jt_in_allowj=0

    jt_in_mj(mj,1)=1
    jt_in_mj(mj,2)=(/(mm,mm=1,nmob)/)
    jt_in_smj(smj,1)=1
    jt_in_smj(smj,2)=(/(mm,mm=1,nsm)/)

    ! define binding energies
    if (allocated(Eb)) deallocate(Eb)
    allocate(Eb(nc,nsm))
    Eb(1:ni-1,1)=Ebin(ni:2:-1)  ! type 1 is I-releasing
    Eb(ni+3:ni+nv+1,2)=Ebvn(2:nv) ! type 2 is V-releasing

    !%%%%%%% set forbidden and allowed index %%%%%%%%
    if (allocated(mask)) deallocate(mask)
    allocate(mask(nc))
    
    !$OMP PARALLEL PRIVATE(mm), DEFAULT(shared)
    !$OMP DO
    do mm=1,nc
	mask(mm)=.false.
	if (hx(mm)==0 .and. vy(mm)==0) then
		mask(mm)=.true.
	else if (hx(mm)>0 .and. vy(mm)<0) then
		mask(mm)=.true.
	else if (Eb(mm,1)<eb_min) then
		mask(mm)=.true.
	else if (Eb(mm,2)<eb_min) then
		mask(mm)=.true.		
	end if		
    end do
    !$OMP END DO
    !$OMP END PARALLEL 
    nallow = nc-count(mask)

    neq=nallow*nx
    if (allocated(allowx)) deallocate(allowx,allowy,allowj,allowj_nx,smj_in_allowj)
    allocate(allowx(nallow),allowy(nallow),allowj(nallow),allowj_nx(neq,2),smj_in_allowj(nsm))
    allowx = pack(hx,.not.(mask))
    allowy = pack(vy,.not.(mask))
    allowj = pack(jt,.not.(mask))
    allowj_nx(:,1)=reshape(spread(allowx,2,nx),(/neq/))
    allowj_nx(:,2)=reshape(spread(allowy,2,nx),(/neq/))
    smj_in_allowj=xufind(nsm,smj,nallow,allowj)
    jt_in_allowj(allowj)=(/(mm,mm=1,nallow)/)
    jt_in_allowj_mat1=reshape(jt_in_allowj,(/ni+nv+1,nh+1/))
    jt_in_allowj_mat2=jt_in_allowj_mat1/=0

    if (old_nh/=nh .or. old_nv/=nv .or. old_ni/=ni) then
        if (allocated(old_allowj_in_new)) deallocate(old_allowj_in_new)
        allocate(old_allowj_in_new(old_nallow))
        if (allocated(mask)) deallocate(mask)
        allocate(mask(nallow))
	!$OMP PARALLEL PRIVATE(mm),DEFAULT(shared)
	!$OMP DO
	do mm=1,nallow
		mask(mm)=.true.
		if (allowx(mm)>old_nh) then
			mask(mm)=.false.		
		else if (allowy(mm)<-old_ni) then
			mask(mm)=.false.		
		else if (allowy(mm)>old_nv) then
			mask(mm)=.false.		
		end if
	end do
	!$OMP END DO
	!$OMP END PARALLEL
        old_allowj_in_new=pack((/(mm,mm=1,nallow)/),mask)
    end if

    if (allocated(mask)) deallocate(mask)
    if (allocated(he_bound)) deallocate(he_bound,v_bound,i_bound)
    allocate(mask(neq))
    mask=allowj_nx(:,1)==nh
    allocate(he_bound(count(mask)))
    he_bound=pack((/(mm,mm=1,neq)/),mask)
    mask=allowj_nx(:,2)==nv
    allocate(v_bound(count(mask)))
    v_bound=pack((/(mm,mm=1,neq)/),mask)
    mask=allowj_nx(:,2)==-ni
    allocate(i_bound(count(mask)))
    i_bound=pack((/(mm,mm=1,neq)/),mask)    

    ! identify captures of mobile species and index the captures for later
    ! calculation of diffusivities in the main loop.

call gettim(hr,minu,sec,hundsec)
print*, 'calculating capture at ', hr,' : ', minu, ' : ', sec, ' : ', hundsec

    if (allocated(capture)) deallocate(capture,capture1,capture3a,cap12,cap34,&
	cap12ja,cap34ia)
    allocate(capture(nc,nmob),intmattemp(nc,nmob),cap12ja(nmob+1),cap34ia(nc+1),cap12&
	(nmob),cap34(nc))
    capture =0; cap12 = 0; cap34 =0; 

!$OMP PARALLEL PRIVATE(i,j,maskmat,matleft,matright,mattop,matbott,mask),DEFAULT(SHARED)
if (allocated(maskmat)) deallocate(maskmat)
allocate(maskmat(-ni:nv,0:nh))
!$OMP Do
    do j=1,nmob
	maskmat=.false.
	matleft=0
	matright=nh-mx(j)
	mattop=max(-ni,-ni-my(j))
	matbott=min(nv,nv-my(j))
	maskmat(mattop:matbott,matleft:matright)=jt_in_allowj_mat2(mattop:matbott,matleft:&
		matright) .and. jt_in_allowj_mat2(mattop+my(j):matbott+my(j),matleft+mx(j):&
		matright+mx(j))
	if (mx(j)==0) then
		if (-my(j)>=-ni .and. -my(j)<=nv) then
			if (jt_in_allowj_mat2(-my(j),0)==.true.) then
				maskmat(-my(j),0)=.true.
			end if
		end if
	end if

        cap12(j) = count(maskmat)
        intmattemp(1:cap12(j),j) = allowj(pack(jt_in_allowj_mat1,maskmat))
        capture(intmattemp(1:cap12(j),j),j) = indx(hx(intmattemp(1:cap12(j),j))+mx(j),vy(intmattemp(&
           1:cap12(j),j))+my(j),cap12(j))
    end do
!$OMP END Do
!$OMP BARRIER
!$OMP SINGLE
    cap12nnz=sum(cap12)
    if (allocated(capture1)) deallocate(capture1)
    allocate(capture1(cap12nnz))
    cap12ja(1)=1
    do j=2,nmob+1
	cap12ja(j)=cap12ja(j-1)+cap12(j-1)
    end do
!$OMP END SINGLE
!$OMP BARRIER
!$OMP DO
    do j=1,nmob
	capture1(cap12ja(j):cap12ja(j+1)-1)=intmattemp(1:cap12(j),j)
    end do
!$OMP END DO
!$OMP BARRIER

    if (allocated(mask)) deallocate(mask)
    allocate(mask(nmob))
!$OMP Do
    do i=1,nallow
        j = allowj(i)
        mask = capture(j,:) /= 0
        cap34(j) = count(mask)
        intmattemp(j,1:cap34(j)) = pack((/(mm,mm=1,nmob)/),mask)
    end do
!$OMP END Do
!$OMP BARRIER
!$OMP SINGLE
    cap34nnz=sum(cap34)
    allocate(capture3a(cap34nnz))
    cap34ia(1)=1
    do j=2,nc+1
	cap34ia(j)=cap34ia(j-1)+cap34(j-1)
    end do
!$OMP END SINGLE
!$OMP BARRIER
!$OMP DO
    do j=1,nc
	capture3a(cap34ia(j):cap34ia(j+1)-1)=intmattemp(j,1:cap34(j))
    end do
!$OMP END DO
!$OMP END PARALLEL
call gettim(hr,minu,sec,hundsec)
print*, 'done calculating capture at ', hr,' : ', minu, ' : ', sec, ' : ', hundsec


    ! identify dissociations of single mobile species for later calculation of
    ! dissoc. rate in the main loop.

call gettim(hr,minu,sec,hundsec)
print*, 'calculating dissociation at ', hr,' : ', minu, ' : ', sec, ' : ', hundsec

    if (allocated(dissociation)) deallocate(dissociation,dissociation1,dissociation2,&
        dissoc12,dissociation3a,dissociation3b,dissociation4,dissoc34)
    allocate(dissociation(nc,nsm),dissociation1(nc,nsm),dissociation2(nc,nsm),dissoc12(nsm),&
        dissociation3a(nc,nsm),dissociation3b(nc,nsm),dissociation4(nc,nsm),dissoc34(nc))
    dissociation =0; dissociation1 =0; dissociation2 =0; dissoc12 =0;
    dissociation3a =0; dissociation3b =0; dissociation4 =0; dissoc34 =0;
if (allocated(maskmat)) deallocate(maskmat)
allocate(maskmat(-ni:nv,0:nh))    
    do j=1,nsm
	! general conditions
	maskmat=.false.
	matleft=smx(j)
	matright=nh
	mattop=max(-ni,-ni+smy(j))
	matbott=min(nv,nv+smy(j))
	maskmat(mattop:matbott,matleft:matright)=jt_in_allowj_mat2(mattop:matbott,matleft:&
		matright) .and. jt_in_allowj_mat2(mattop-smy(j):matbott-smy(j),matleft-smx(j):&
		matright-smx(j))

	! special conditions
	if (j==1) then
		maskmat(1:nv,0) =.false. ! forbid I-release from V_n clusters.
	else if (j==2) then
		maskmat(-ni:-1,0) =.false. ! forbid V-release from I_n clusters	
	end if

        dissoc12(j) = count(maskmat)
        dissociation1(1:dissoc12(j),j) = allowj(pack(jt_in_allowj_mat1,maskmat))
        dissociation(dissociation1(1:dissoc12(j),j),j) = indx(hx(dissociation1(1:dissoc12(j),j)&
            )-smx(j),vy(dissociation1(1:dissoc12(j),j))-smy(j),dissoc12(j))
        dissociation2(1:dissoc12(j),j) = dissociation(dissociation1(1:dissoc12(j),j),j)
    end do
    if (allocated(mask)) deallocate(mask)
    allocate(mask(nsm))
    do i=1,nallow
            j = allowj(i)
            mask = dissociation(j,:) /= 0
            dissoc34(j) = count(mask) 
            dissociation3a(j,1:dissoc34(j)) = pack((/(mm,mm=1,nsm)/),mask)
            dissociation3b(j,1:dissoc34(j)) = pack(smj,mask)
            dissociation4(j,1:dissoc34(j)) = dissociation(j,dissociation3a(j,1:dissoc34(j)))
    end do

call gettim(hr,minu,sec,hundsec)
print*, 'done calculating dissociation at ', hr,' : ', minu, ' : ', sec, ' : ', hundsec

    ! identify virtual dissociations which will be used to facilitate later
    ! analysis on the formation of clusters by capturing which is in the inverse direction.

call gettim(hr,minu,sec,hundsec)
print*, 'calculating virtual dis at ', hr,' : ', minu, ' : ', sec, ' : ', hundsec

    if (allocated(virtualdis3a)) deallocate(virtualdis3a,virtualdis4,&
        virtd34,virtd34ia)
    allocate(virtualdis(nc,nmob),virtd34(nc),virtd34ia(nc+1))
    virtualdis =0; virtd34 =0;
!$OMP PARALLEL PRIVATE(j,maskmat,matleft,matright,mattop,matbott,mm,intarraytemp,mask,i), &
!$OMP  DEFAULT(SHARED)
if (allocated(maskmat)) deallocate(maskmat)
allocate(maskmat(-ni:nv,0:nh))

!$OMP do
    do j=1,nmob
	maskmat=.false.
	matleft=mx(j)
	matright=nh
	mattop=max(-ni,-ni+my(j))
	matbott=min(nv,nv+my(j))
	maskmat(mattop:matbott,matleft:matright)=jt_in_allowj_mat2(mattop:matbott,matleft:&
		matright) .and. jt_in_allowj_mat2(mattop-my(j):matbott-my(j),matleft-mx(j):&
		matright-mx(j))
        mm = count(maskmat)
        allocate(intarraytemp(mm))
        intarraytemp = allowj(pack(jt_in_allowj_mat1,maskmat))  
        virtualdis(intarraytemp,j) = indx(hx(intarraytemp)-mx(j),vy(intarraytemp)-my(j),mm)
        deallocate(intarraytemp)
    end do
!$OMP end do
!$OMP barrier
    if (allocated(mask)) deallocate(mask)
    allocate(mask(nmob))
!$OMP do
    do i=1,nallow
            j = allowj(i)
            mask = virtualdis(j,:) /= 0
            virtd34(j) = count(mask)
            intmattemp(j,1:virtd34(j)) = pack((/(mm,mm=1,nmob)/),mask)
    end do
!$OMP end do
!$OMP BARRIER
!$OMP SINGLE
    virtd34nnz=sum(virtd34)
    allocate(virtualdis3a(virtd34nnz),virtualdis4(virtd34nnz))
    virtd34ia(1)=1
    do i=2,nc+1
	virtd34ia(i)=virtd34ia(i-1)+virtd34(i-1)
    end do
!$OMP END SINGLE
!$OMP BARRIER
!$OMP do
    do j=1,nc
	if (virtd34(j)/=0) then
		virtualdis3a(virtd34ia(j):virtd34ia(j+1)-1)=intmattemp(j,1:virtd34(j))
		virtualdis4(virtd34ia(j):virtd34ia(j+1)-1)=virtualdis(j,intmattemp(j,1:virtd34(j)))
	end if
    end do
!$OMP end do
!$OMP end parallel

deallocate(virtualdis,intmattemp)

call gettim(hr,minu,sec,hundsec)
print*, 'done calculating virtual dis at ', hr,' : ', minu, ' : ', sec, ' : ', hundsec

    micros=0.d0  ! 2.5e-4;

    ! define formation energies for I and V
    select case (ef_system)             !######
    case ('AB')
        efi = 6.5d0; efv = 3.d0
    case ('MD')
        efi = 4.88d0; efv = 1.71d0
    case default 
        stop 'Ef-system not recognized; restart the program.'
    end select


    ! calculate unrelaxed interface radii (= capture radii - r00) and volume;
    if (allocated(r)) deallocate(r)
    allocate(r(nc))
        r = abs(dble(vy))**(1.d0/3.)*0.139d0+r0    !######
	r(1:ni-1)=(dble((/(kk,kk=ni,2,-1)/))/64.1d0/.273d0/pi)**.5d0

    ! process initial conditions
    if (ntt0==0) then
        ! calculate initial conditions
        if (allocated(y0)) deallocate(y0)
        allocate(y0(neq))
	if (allocated(mattemp)) deallocate(mattemp)
	allocate(mattemp(nallow,nx))
	mattemp=0.d0
	mattemp(smj_in_allowj(2),2:nx-1)=64.1d0*exp(-efv/8.617d-5/tempr)
	mattemp(smj_in_allowj(1),iblankleft+1:iblankright-1)=64.1d0*exp(-efi/8.617d-5/tempr)
	y0=reshape(mattemp,(/neq/))
	deallocate(mattemp)

        ! save first set of data
        write(fh2), tspan(1)
        ntt0=1
        write(fh3), y0
        old_nallow=nallow
        write(fh4,'(9(a15))') 'nh','nv','ni','nallow','old_nh','old_nv','old_ni','old_nallow','nx'
        write(fh4,'(9(I15))') nh, nv, ni,nallow, old_nh, old_nv, old_ni, old_nallow, nx
        rewind(fh4)
        write(fh5), x, allowx, allowy
        rewind(fh5)
        write(fh6,*), ntt0
        rewind(fh6)
    else if (old_nh==nh .and. old_nv==nv .and. old_ni==ni) then
        ! load initial conditions
        if (allocated(y0)) deallocate(y0)
        allocate(y0(neq))
	rewind(fh3)
	do mm=1,ntt0
	        read(fh3), y0
	end do
    else
        ! tim=clock; !sprintf('Loading initial conditions at %g; %g; %g;', tim(4:6))     %// Test
	old_neq=old_nallow*nx
        allocate(mattemp(old_neq,ntt0))
	rewind(fh3)
	do mm=1,ntt0
	        read(fh3), mattemp(:,mm)
	end do
        rewind(fh3)

        allocate(mattemp2(neq,ntt0))
        mattemp2=0.d0
        do j=1,nx
            mattemp2(old_allowj_in_new+(j-1)*nallow,:)=mattemp((/(mm,mm=1,old_nallow)/)+(j-1)*old_nallow,:)
        end do
        deallocate(mattemp)
        
	do mm=1,ntt0
	        write(fh3), mattemp2(:,mm)
	end do

        write(fh4,'(9(a15))') 'nh','nv','ni','nallow','old_nh','old_nv','old_ni','old_nallow','nx'
        write(fh4,'(9(I15))') nh, nv, ni,nallow, old_nh, old_nv, old_ni, old_nallow, nx
        rewind(fh4)
        write(fh5), x, allowx, allowy
        rewind(fh5)
        
        if (allocated(y0)) deallocate(y0)
        allocate(y0(neq))
        y0=mattemp2(:,ntt0)
        deallocate(mattemp2)
    end if

    dbl_smj=indx(2*smx,2*smy,nsm)     
    if (allocated(zz2)) deallocate(zz2,zz3,sum_array,coef1,coef2)
    allocate(zz2(nc,nx),zz3(nc,nx),sum_array(nc),coef1(virtd34nnz),coef2(nc,nsm))
    zz2=0.d0; zz3=0.d0; sum_array=0.d0; coef1=0; coef2=0;

    if (allocated(mask)) deallocate(mask)    
    !$OMP PARALLEL PRIVATE(kk, j, mask, mm), DEFAULT(shared)
    !$OMP DO
    do kk=1,nallow
        j=allowj(kk)
        coef1(virtd34ia(j):virtd34ia(j+1)-1)=1
        if (allocated(mask)) deallocate(mask)
        allocate(mask(virtd34(j)))
        mask=jt_in_mj(virtualdis4(virtd34ia(j):virtd34ia(j+1)-1),1) .and. virtualdis4(virtd34ia(j):&
		virtd34ia(j+1)-1)/=mj(virtualdis3a(virtd34ia(j):virtd34ia(j+1)-1))
        coef1(virtd34ia(j)-1+pack((/(mm,mm=1,virtd34(j))/),mask))=2
        coef2(j,1:dissoc34(j))=1
        deallocate(mask)
        allocate(mask(dissoc34(j)))
        mask=jt_in_smj(dissociation4(j,1:dissoc34(j)),1) .and. dissociation4&
            (j,1:dissoc34(j)) /= dissociation3b(j,1:dissoc34(j))
        coef2(j,pack((/(mm,mm=1,dissoc34(j))/),mask))=2
    end do
    !$OMP END DO
    !$OMP END PARALLEL 

    !!! set global k d dis u
    if (allocated(d)) deallocate(d,k,dis,u)
    allocate(d(nc),k(nc,nmob),dis(nc,nsm),u(nc,nx))
    d=0.d0; k=0.d0; dis=0.d0; 

    ! calculate all the rates
    ceq(1)=64.1d0*exp(-efi/8.617d-5/tempr)  ! equilibrium conc. of smj(1), i.e., I
    ceq(2)=64.1d0*exp(-efv/8.617d-5/tempr)   ! equilibrium conc. of smj(2), i.e., V

    d(mj)=D0*exp(-mq/8.617d-5/tempr)   ! diffusivities of mobile species

    do kk=1,nmob
        allocate(intarraytemp(cap12(kk)))
        intarraytemp=capture1(cap12ja(kk):cap12ja(kk+1)-1)
        
        k(intarraytemp,kk)=4.d0*pi*(r(intarraytemp)+r(mj(kk)))*(d(intarraytemp)+d(mj(kk))) 
        deallocate(intarraytemp)
        if (mx(kk)==0) then          !******************
            if (my(kk)<0) k(1:ni,kk)=k(1:ni,kk)*bias       !******************
        end if       !******************
    end do
	k(smj(1),smj_in_mj(2))=	4.d0*pi*dis_recomb*(d(smj(1))+d(smj(2)))
	k(smj(2),smj_in_mj(1))= k(smj(1),smj_in_mj(2))

    do kk=1, nsm
        allocate(intarraytemp(dissoc12(kk)))
        intarraytemp=dissociation1(1:dissoc12(kk),kk)
        dis(intarraytemp,kk)=k(dissociation(intarraytemp,kk),smj_in_mj(kk))*64.1d0*exp(-Eb&
            (intarraytemp,kk)/8.617d-5/tempr)
        deallocate(intarraytemp)
    end do
   
    call xusolver
    if (ntt0==ntspan0) exit
end do    ! end the outermost loop

tend=dclock()
print*,'total time used before writing .txt: ',tend-tstart, ' seconds'

deallocate(dfdy_ia, dfdy_ja, dfdy_csr,capture, capture1, capture3a, dissociation,&
        dissociation1, dissociation2, dissociation3a, dissociation3b, dissociation4,  &
	virtualdis3a, virtualdis4, coef1)

END PROGRAM main

!********************* ode setup **************************************
SUBROUTINE xuode(ndum,tdum,ydum,ypdum)

use mod1
use modsolverinterface
use omp_lib
implicit none

integer ::  sn, zz, ndum
real *8 :: tdum, ydum(ndum), ypdum(ndum), up(nc,nx), upx(nc,nx), dudx(nc,nx-1)

call omp_set_num_threads(num_threads)
u = 0.d0; up =0.d0; upx=0.d0; dudx=0.d0
ypdum=0.d0

u(allowj,:)=reshape(ydum,(/nallow,nx/))

up(ni-19:ni,2:nx-1)=up(ni-19:ni,2:nx-1)+imp(1:20,2:nx-1)
up(ni+2:ni+10,2:nx-1)=up(ni+2:ni+10,2:nx-1)+imp(21:29,2:nx-1)
dudx(mj,:)=(u(mj,2:)-u(mj,:nx-1))/spread(dx,1,nmob)
upx(mj,2:nx-1)=(dudx(mj,2:)-dudx(mj,:nx-2))/spread(dx2,1,nmob)
up(mj,2:nx-1)=up(mj,2:nx-1)+spread(d(mj),2,nx-2)*upx(mj,2:nx-1)

if (allocated(zz1)) deallocate(zz1)
!$OMP PARALLEL PRIVATE(kk,sn,zz,zz1), DEFAULT(shared)
!$OMP DO SCHEDULE (static)
do kk=1, nallow
	sn = allowj(kk)

	!!%% formation rate by all possible trappings
	allocate(zz1(virtd34(sn)))
	do zz = 1, virtd34(sn)
		zz1(zz)=k(virtualdis4(virtd34ia(sn)-1+zz),virtualdis3a(virtd34ia(sn)-1+zz))/coef1(virtd34ia(sn)-1+zz)
	end do
	up(sn,2:nx-1) = up(sn,2:nx-1)+matmul(zz1,u(mj(virtualdis3a(virtd34ia(sn):virtd34ia(sn+1)-1)),2:nx-1)*u &
			(virtualdis4(virtd34ia(sn):virtd34ia(sn+1)-1),2:nx-1))
	deallocate(zz1)
    ! annihilation rate by type I-III dissociations
    sum_array(sn)=sum(dis(sn,dissociation3a(sn,1:dissoc34(sn)))/coef2(sn,1:dissoc34(sn)))
    up(sn,2:nx-1)=up(sn,2:nx-1)-sum_array(sn)*u(sn,2:nx-1) 

   ! annihilation by all possible trappings
    if (.not.(jt_in_mj(sn,1))) then
        zz2(sn,2:nx-1)=matmul(k(sn,capture3a(cap34ia(sn):cap34ia(sn+1)-1)),u(mj(capture3a(cap34ia(sn):&
		cap34ia(sn+1)-1)),2:nx-1))
        up(sn,2:nx-1)=up(sn,2:nx-1)-zz2(sn,2:nx-1)*u(sn,2:nx-1)

    else 
	zz=jt_in_mj(sn,2) 
        zz3(sn,2:nx-1)=matmul(k(capture1(cap12ja(zz):cap12ja(zz+1)-1),zz),u(capture1(cap12ja(zz):cap12ja&
		(zz+1)-1),2:nx-1))
        up(sn,2:nx-1)=up(sn,2:nx-1)-zz3(sn,2:nx-1)*u(sn,2:nx-1)-k(sn,zz)*u(sn,2:nx-1)**2
        if (sn==smj(1) .or. sn==smj(2)) up(sn,2:nx-1)=up(sn,2:nx-1)+k(smj(1),smj_in_mj(2))*ceq(1)&
            *ceq(2)
        if (hx(sn)==0) then
            if (vy(sn)<0) then
                up(sn,2:nx-1)=up(sn,2:nx-1)-micros*bias*d(sn)*u(sn,2:nx-1)
            else
                up(sn,2:nx-1)=up(sn,2:nx-1)-micros*d(sn)*u(sn,2:nx-1)
                if (vy(sn)==1) up(sn,2:nx-1)=up(sn,2:nx-1)+micros*d(sn)*ceq(2)
            end if
        end if
    end if


    ! formation by type I-III dissociations
    if (.not.(jt_in_smj(sn,1))) then
        do zz=1,nsm
            if (capture(sn,smj_in_mj(zz))/=0) then
                up(sn,2:nx-1)=up(sn,2:nx-1)+dis(capture(sn,smj_in_mj(zz)),zz)*u(capture(sn,&
                    smj_in_mj(zz)),2:nx-1)
            end if
        end do
    else 
	zz=jt_in_smj(sn,2)
        up(sn,2:nx-1)=up(sn,2:nx-1)+matmul(dis(dissociation1(1:dissoc12(zz),zz),zz),u(dissociation1&
            (1:dissoc12(zz),zz),2:nx-1))+dis(dbl_smj(zz),zz)*u(dbl_smj(zz),2:nx-1) 
    end if

end do
!$OMP END DO
!$OMP END PARALLEL

up(1:ni,1:iblankleft)=0.d0
up(1:ni,iblankright:nx)=0.d0

ypdum=reshape(up(allowj,:),(/ndum/))
END SUBROUTINE xuode

!---------------------------------------------------------------
SUBROUTINE xujac(ndum,tdum,ydum)

use mod1
use modsolverinterface
use omp_lib
implicit none
integer ::  sn, zz, ndum, lit, xx, xk, lociter, i0, i1, rownnz(neq), locmin, locmax
integer*8 :: cnt, xu1, loccnt, locnzmax
real *8 :: tdum, ydum(ndum)
integer, allocatable :: locin(:),locjn(:)
real*8, allocatable:: locval(:)
logical :: locerror, gerror, first
call omp_set_num_threads(num_threads)
call gettim(hr,minu,sec,hundsec)
print*, 'starting xujac at ', hr,' : ', minu, ' : ', sec, ' : ', hundsec
if (allocated(dfdy_in)) deallocate(dfdy_in,dfdy_jn,dfdy_val)
if (allocated(dfdy_ia)) deallocate(dfdy_ia,dfdy_ja,dfdy_csr)
allocate(dfdy_ia(neq+1))

if (allocated(arraytemp)) deallocate(arraytemp)
if (allocated(zz1)) deallocate(zz1)
dfdy_ia(1)=1
cnt=0
locnzmax=maxnnz

u=0.d0
u(allowj,:)=reshape(ydum,(/nallow,nx/))

!$OMP PARALLEL PRIVATE(first,lociter,loccnt,locerror,locin,locjn,locval,xk,i0,kk,xu1,xx,sn,&
!$OMP	arraytemp,zz,i1,locmin,locmax,lit,zz1), DEFAULT(shared)

allocate(arraytemp(nc))
211 arraytemp=0.d0
first=.true.
lociter=0
loccnt=0
locerror=.false.
gerror=.false.
allocate(locin(locnzmax),locjn(locnzmax),locval(locnzmax))

!$OMP DO SCHEDULE (STATIC)
do xk=1,neq
    if (gerror) cycle      ! in case of local space shortage, quickly finish the loop
    if (first) i0= xk
    kk=mod(xk-1,nallow)+1   ! index of the cluster in allowj
    xu1=xk-kk
    xx=xu1/nallow+1           ! index of x-grid

        sn=allowj(kk)
	lit=0
	if (xx==1 .or. xx==nx .or. ((xx<=iblankleft .or. xx>=iblankright) .and. kk<=ni)) then   
		lit=lit+1
		loccnt=loccnt+1
		if (loccnt>locnzmax) then
			locerror=.true.
			gerror=.true.
			goto 555
		else 
			locin(loccnt)=xu1+kk
			locjn(loccnt)=xu1+kk
			locval(loccnt)= 0.d0
			goto 555
		end if
	end if
	
        !!! formation by all possible trappings
	allocate(zz1(virtd34(sn)))
	do zz = 1, virtd34(sn)
		zz1(zz)=k(virtualdis4(virtd34ia(sn)-1+zz),virtualdis3a(virtd34ia(sn)-1+zz))/coef1(virtd34ia(sn)-1+zz)
	end do
        arraytemp(mj(virtualdis3a(virtd34ia(sn):virtd34ia(sn+1)-1))) = arraytemp(mj(virtualdis3a(virtd34ia&
		(sn):virtd34ia(sn+1)-1)))+zz1*u(virtualdis4(virtd34ia(sn):virtd34ia(sn+1)-1),xx)
        locmin=minval(mj(virtualdis3a(virtd34ia(sn):virtd34ia(sn+1)-1)))
        locmax=maxval(mj(virtualdis3a(virtd34ia(sn):virtd34ia(sn+1)-1)))
        arraytemp(virtualdis4(virtd34ia(sn):virtd34ia(sn+1)-1)) = arraytemp(virtualdis4(virtd34ia(sn):&
		virtd34ia(sn+1)-1))+zz1*u(mj(virtualdis3a(virtd34ia(sn):virtd34ia(sn+1)-1)),xx)
	deallocate(zz1)
        locmin=min(locmin,minval(virtualdis4(virtd34ia(sn):virtd34ia(sn+1)-1)))
        locmax=max(locmax,maxval(virtualdis4(virtd34ia(sn):virtd34ia(sn+1)-1)))        
        !!! annihilation by type I-III dissociations
	sum_array(sn)=sum(dis(sn,dissociation3a(sn,1:dissoc34(sn)))/coef2(sn,1:dissoc34(sn)))
        arraytemp(sn) = arraytemp(sn)-sum_array(sn)
	locmin=min(locmin,sn)
	locmax=max(locmax,sn)

        !!! annihilation by all possible trappings
        if (.not.(jt_in_mj(sn,1))) then
	     zz2(sn,1:nx)=matmul(k(sn,capture3a(cap34ia(sn):cap34ia(sn+1)-1)),u(mj(capture3a(cap34ia(sn):&
		cap34ia(sn+1)-1)),1:nx))
             arraytemp(sn) = arraytemp(sn)-zz2(sn,xx)
             arraytemp(mj(capture3a(cap34ia(sn):cap34ia(sn+1)-1))) = arraytemp(mj(capture3a(cap34ia(sn):&
		cap34ia(sn+1)-1)))-k(sn,capture3a(cap34ia(sn):cap34ia(sn+1)-1))*u(sn,xx)
	     locmin=min(locmin,minval(mj(capture3a(cap34ia(sn):cap34ia(sn+1)-1))))
	     locmax=max(locmax,maxval(mj(capture3a(cap34ia(sn):cap34ia(sn+1)-1))))
        else 
	    zz=jt_in_mj(sn,2) 
	    zz3(sn,1:nx)=matmul(k(capture1(cap12ja(zz):cap12ja(zz+1)-1),zz),u(capture1(cap12ja(zz):cap12ja&
			(zz+1)-1),1:nx))
            arraytemp(sn) = arraytemp(sn)-zz3(sn,xx)-2*k(sn,zz)*u(sn,xx)
            arraytemp(capture1(cap12ja(zz):cap12ja(zz+1)-1)) = arraytemp(capture1(cap12ja(zz):cap12ja(zz+1)&
			-1))-k(capture1(cap12ja(zz):cap12ja(zz+1)-1),zz)*u(sn,xx)
            if (hx(sn)==0) then
                if (vy(sn)<0) then
                    arraytemp(sn) = arraytemp(sn)-micros*bias*d(sn)
                else
                    arraytemp(sn) = arraytemp(sn)-micros*d(sn)
                end if
            end if    
	    locmin=min(locmin,minval(capture1(cap12ja(zz):cap12ja(zz+1)-1)))
	    locmax=max(locmax,maxval(capture1(cap12ja(zz):cap12ja(zz+1)-1)))
        end if

        ! formation by type I-III dissociations
        if (.not.(jt_in_smj(sn,1))) then
            do zz=1,nsm
                if (capture(sn,smj_in_mj(zz))/=0) then
                    arraytemp(capture(sn,smj_in_mj(zz))) = arraytemp(capture(sn,smj_in_mj(zz)))+&
                        dis(capture(sn,smj_in_mj(zz)),zz)
		    locmin=min(locmin,capture(sn,smj_in_mj(zz)))
		    locmax=max(locmax,capture(sn,smj_in_mj(zz)))
                end if
            end do
        else 
	    zz=jt_in_smj(sn,2)
            arraytemp(dissociation1(1:dissoc12(zz),zz)) = arraytemp(dissociation1(1:dissoc12(zz),zz))&
                +dis(dissociation1(1:dissoc12(zz),zz),zz)
            arraytemp(dbl_smj(zz)) = arraytemp(dbl_smj(zz))+dis(dbl_smj(zz),zz)
	    locmin=min(locmin,minval(dissociation1(1:dissoc12(zz),zz))) 
	    locmax=max(locmax,maxval(dissociation1(1:dissoc12(zz),zz))) 
	    locmin=min(locmin,dbl_smj(zz))
	    locmax=max(locmax,dbl_smj(zz))
        end if
        
            
        ! if there are spatial grids, diffusion terms should be included, here for sn at xx and xx-1; record 
	! the diffusion term for sn at xx-1 first, then the reaction terms at xx, then diffusion term for sn at xx+1
	if (jt_in_mj(sn,1)) then
		if (xx>1 .and. xx<nx) then
			arraytemp(sn)=arraytemp(sn)-d(sn)*2*(1/(x(xx+1)-x(xx))+1/(x(xx)-x(xx-1)))/(x(xx+1)-x(xx-1))
			lit=lit+1
			loccnt=loccnt+1
			if (loccnt>locnzmax) then
				locerror=.true.
				gerror=.true.
				goto 555
			else 
				locin(loccnt)=xu1+kk
				locjn(loccnt)=xu1+kk-nallow
				locval(loccnt)= d(sn)*2*(1/(x(xx)-x(xx-1)))/(x(xx+1)-x(xx-1))
			end if
		end if
	end if     
        
        ! now record non-zeros due to reaction terms at xx spatial point

	do zz=locmin,locmax
		if (arraytemp(zz)/=0.d0 .or. zz==sn) then
			lit=lit+1
			loccnt=loccnt+1
			if (loccnt>locnzmax) then
				locerror=.true.
				gerror=.true.
				goto 555
			else 
				locin(loccnt)=xu1+kk
				locjn(loccnt)=xu1+jt_in_allowj(zz)
				locval(loccnt)=arraytemp(zz)
			end if
			arraytemp(zz)=0.d0
		end if
	end do

	if (jt_in_mj(sn,1)) then
		if (xx>1 .and. xx<nx) then
			lit=lit+1
			loccnt=loccnt+1
			if (loccnt>locnzmax) then
				locerror=.true.
				gerror=.true.
				goto 555
			else 
				locin(loccnt)=xu1+kk
				locjn(loccnt)=xu1+kk+nallow
				locval(loccnt)= d(sn)*2*(1/(x(xx+1)-x(xx)))/(x(xx+1)-x(xx-1))
			end if
		end if
	end if     

555	rownnz(xk)=lit
	first=.false.
	lociter=lociter+1
end do
!$OMP END DO
if (locerror) gerror=.true.
!$OMP BARRIER

if (gerror) then
	deallocate(locin,locjn,locval)
	!$OMP SINGLE
	locnzmax=locnzmax*1.5
	print*,'locnzmax is increased by 50%'
	!$OMP END SINGLE
	go to 211
else
	i1=i0+lociter-1
	!$OMP BARRIER
	!$OMP MASTER
	do zz=2,neq+1
		dfdy_ia(zz)=dfdy_ia(zz-1)+rownnz(zz-1)
	end do
	cnt=dfdy_ia(neq+1)-1
	allocate(dfdy_ja(cnt),dfdy_csr(cnt))
	if (npds==0) allocate(dfdy_in(cnt))	
	!$OMP END MASTER
	!$OMP BARRIER
	if (loccnt>0) then
		if (npds==0) dfdy_in(dfdy_ia(i0):dfdy_ia(i1+1)-1)=locin(1:loccnt)
		dfdy_ja(dfdy_ia(i0):dfdy_ia(i1+1)-1)=locjn(1:loccnt)
		dfdy_csr(dfdy_ia(i0):dfdy_ia(i1+1)-1)=locval(1:loccnt)
	end if
end if

!$OMP END PARALLEL

if (dfdy_nnz /= cnt) then
	print*,'dfdy_nnz has been changed : ', cnt
	dfdy_nnz=cnt
end if

if (npds==0) then
	allocate(dfdy_jn(dfdy_nnz),dfdy_val(dfdy_nnz))
	dfdy_jn=dfdy_ja
	dfdy_val=dfdy_csr
end if

call gettim(hr,minu,sec,hundsec)
print*, 'exiting xujac at ', hr,' : ', minu, ' : ', sec, ' : ', hundsec

END SUBROUTINE xujac


!********************************************************************************
SUBROUTINE xuoutput(v1,n1,v2)
use modfilehandles
use modsolverinterface
implicit none
integer n1, i
real*8  v1(n1), v2(neq,n1)

do i=1,n1
write(fh2),v1(i)
write(fh3),v2(:,i)
end do

ntt0=ntt0+n1
write(fh6,*), ntt0
rewind(fh6)
END SUBROUTINE xuoutput 

!********************************************************************************
!********************************************************************************
SUBROUTINE xusolver
#include 'ilupack_fortran.h'
#include 'ilupackmacros.h'
use omp_lib
use modsolverinterface
implicit none
external:: xuode, xujac
integer :: nsteps, nfailed, nfevals, ndecomps, nsolves, next, &
            maxk, max1tok(5), kI(5,5), kJ(5,5), pardiso_type, &
            maxit, solverk, klast, nconhk, job(8), info, &
            error, iparm(64),iparm0(64), iter, nout_new, kopt, jj, pardisolver
integer*8  pt(64), idum            
integer, allocatable :: solverkcap(:), speye_ia(:), speye_ja(:), Miter_ia(:), Miter_ja(:), &
            out_kI(:,:), ib(:), jb(:)
real*8 :: htspan, tfinal, t, threshold, hmax, G(5), alpha(5), &
            invGa(5), erconst(5), difU(5,5), hmin, rh, absh, h, tdel, abshlast,&
            hinvGak, ddum, difRU(5,5), tnew, minnrm, newnrm, errit, rate, oldnrm,&
            err, hopt, errkm1, hkm1, temp, errkp1, hkp1, dparm(64)
real*8, allocatable :: y(:), ynew(:), f0(:), yp(:), wt(:), f1(:), &
	    dfdt(:), dif(:,:), speye_csr(:), Miter_csr(:), Miterb(:), &
	    tout(:), yout(:), psi(:), mattemp(:,:), pred(:), difkp1(:), invwt(:), &
	    rhs(:),  del(:), tout_new(:), yout_new(:,:), out_s(:), rhs1(:), del1(:)
logical :: bdf, Jcurrent, havrate, done, at_hmin, nofailed, gotynew, tooslow, &
            NNreset_dif, Mcurrent, mask1, mask2, mask3
logical, allocatable :: mask0(:)

integer  ILUPACKFACTOR, ILUPACKSOLVER, ILUPACKNNZ
external ILUPACKINIT, ILUPACKSOL, ILUPACKDELETE,&
	ILUPACKFACTOR, ILUPACKSOLVER,&
	ILUPACKNNZ, ILUPACKINFO
!    ILUPACK external parameters
integer   matching, ilumaxit, lfil, lfilS, nrestart, ierr
character ordering*20
REAL*8     droptol, droptolS, condest, restol, elbow
!    variables that cover the and pass the C-pointers
INTEGER*8 param,PREC

call omp_set_num_threads(num_threads)
if (allocated(y)) &
	deallocate(y,f0,yp,wt,f1,dfdt,pred,difkp1,invwt,rhs,del,ynew,psi,mask0)

allocate(y(neq),f0(neq),yp(neq),wt(neq),f1(neq),dfdt(neq),pred(neq),difkp1(neq),&
    invwt(neq),rhs(neq),del(neq),ynew(neq),psi(neq),mask0(neq),rhs1(neq),del1(neq))

nsteps   = 0 
nfailed  = 0 
nfevals  = 0 
npds     = 0 
ndecomps = 0 
nsolves  = 0 
info = 0
htspan = abs(tspan(2) - tspan(1)) 
next = 2        
tfinal = tspan(ntspan) 

call xuode(neq,t0,y0,f0)

if (rtol < 100.d0 * epsilon(1.d0) ) then
  rtol = 100.d0 * epsilon(1.d0)
  print*,'Warning-- RelTolIncrease: RelTol has been increased to ',rtol
end if

threshold = atol / rtol 
hmax=0.1d0*(tfinal-t0) 


t=t0
y=y0

bdf=.false.
maxk=5    ! order of approximation, the higher order, the more precise, up to 5

maxit=40   !  no. of newton iterations to be taken

G = (/1.d0, 3.d0/2, 11.d0/6, 25.d0/12, 137.d0/60/) 
if (bdf) then
    alpha = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0/) 
else
    alpha = (/-37.d0/200, -1.d0/9, -0.0823d0, -0.0415d0, 0.d0/) 
end if
invGa = 1.d0/ (G * (1.d0 - alpha)) 
erconst = alpha * G + (1.d0/ (/2.,3.,4.,5.,6./)) 

difU(1,1:5)=(/-1,-2,-3,-4,-5/)
difU(2,1:5)=(/0,  1,  3,  6,  10/)
difU(3,1:5)=(/0,  0, -1, -4, -10/)
difU(4,1:5)=(/0,  0,  0,  1,   5/)
difU(5,1:5)=(/0,  0,  0,  0,  -1/)

max1tok=(/1,2,3,4,5/)
kJ=spread((/1,2,3,4,5/),1,5)
kI=spread((/1,2,3,4,5/),2,5)

yp=f0

call xujac(neq,t,y) !output sparse Jac to dfdy_in,dfdy_jn,dfdy_val in mod1
npds=npds+1

Jcurrent = .true.

hmin =1.d-300
!hmin = 16.d0*epsilon(1.d0)*abs(t)

if (set_htry=='n') then

	!$OMP PARALLEL PRIVATE(mm), DEFAULT(shared)
	!$OMP DO
	do mm=1,neq
	wt(mm)=max(abs(y(mm)), threshold)
	end do 
	!$OMP END DO
	!$OMP END PARALLEL 

	rh = 1.25d0 * maxval(abs(yp / wt)) / sqrt(rtol) 

	absh = min(hmax, htspan) 
	if (absh * rh > 1.d0) then
	    absh = 1.d0 / rh 
	end if
	absh = max(absh, hmin) 

	h=absh
	tdel=min(sqrt(epsilon(1.d0))*max(abs(t),abs(t+h)),absh)
	call xuode(neq,t+tdel,y,f1)

	nfevals = nfevals + 1 
	dfdt = (f1 - f0) / tdel 
	call mkl_dcoomv('n', neq, neq, 1.d0, 'GF', dfdy_val, dfdy_in, dfdy_jn, dfdy_nnz, &
	    yp, 1.d0, dfdt)
	deallocate(dfdy_val,dfdy_in,dfdy_jn)
	rh = 1.25d0 * sqrt(0.5d0 * maxval(abs( dfdt / wt)) / rtol) 

	absh = min(hmax, htspan) 
	if (absh * rh > 1.d0) then
	    absh = 1.d0 / rh 
	end if
	absh = max(absh, hmin) 

else
	absh = min(hmax, max(hmin, htry))
end if

h=absh

solverk=1
if (allocated(solverkcap)) deallocate(solverkcap)
allocate(solverkcap(solverk))
solverkcap=(/(mm,mm=1,solverk)/)
klast=solverk
abshlast=absh

if (allocated(dif)) deallocate(dif)
allocate(dif(neq,maxk+2))
dif=0.d0
dif(:,1)=h*yp

hinvGak = h * invGa(solverk) 
nconhk = 0  

if (allocated(speye_ia)) deallocate(speye_ia,speye_ja,speye_csr)
allocate(speye_ia(neq+1),speye_ja(neq),speye_csr(neq))
speye_ia=(/(mm,mm=1,neq+1)/)
speye_ja=(/(mm,mm=1,neq)/)
speye_csr=(/(1.d0,mm=1,neq)/)

if (allocated(Miter_csr)) deallocate(Miter_csr)
if (allocated(ib)) deallocate(ib,jb,Miterb)
allocate(Miter_csr(dfdy_nnz))
allocate(ib(neq+1),jb(dfdy_nnz),Miterb(dfdy_nnz))
call mkl_dcsradd('n',0,idum,neq,neq,speye_csr,speye_ja,speye_ia,-hinvGak, &
dfdy_csr,dfdy_ja,dfdy_ia,Miter_csr,dfdy_ja,dfdy_ia,dfdy_nnz,info)
!dfdy_csr,dfdy_ja,dfdy_ia,Miter_csr,Miter_ja,Miter_ia,dfdy_nnz,info)
ib=dfdy_ia
jb=dfdy_ja
Miterb=Miter_csr

if (info /= 0) then
    print*,'The following ERROR was detected by mkl_dcsradd (1st): ', info
    STOP
end if
   
!initialize ILUPACK
     call ILUPACKINIT(neq,ib,jb,Miterb,matching,ordering,droptol, &
	droptolS, condest, restol, ilumaxit, elbow, lfil, lfilS, nrestart)
	ordering='metisn'
	print*,'ordering is : ', ordering

!    maximum weight matching
!    default value is different from zero, matching turned on
    matching=0

!    multilevel orderings
	
!    'amd' (default) Approximate Minimum Degree
!    'mmd'           Minimum Degree            
!    'rcm'           Reverse Cuthill-McKee
!    'metisn'        Metis multilevel nested dissection by nodes
!    'metise'        Metis multilevel nested dissection by edges
!    'pq'            ddPQ strategy by Saad

!    threshold for ILU, default: 1e-2
      droptol= ILUdroptol

!    threshold for the approximate Schur complements, default: 0.1*droptol
      droptolS=0.1d0*droptol

!    norm bound for the inverse factors L^{-1}, U^{-1}, default: 1e2
      condest=1.d12

!    relative error for the backward error (SPD case: relative energy
!    norm) used during the iterative solver, default: sqrt(eps)
      restol=1.d-4

!    maximum number of iterative steps, default: 500
      ilumaxit=500

call gettim(hr,minu,sec,hundsec)
print*,'calling ILUPACKfactor, location 1, at ', hr,' : ', minu, ' : ', sec, ' : ', hundsec
ierr=ILUPACKFACTOR(param,PREC,neq,ib,jb,Miterb,matching,ordering,droptol, &
droptolS, condest, restol, ilumaxit, elbow, lfil, lfilS, nrestart)
      if (ierr.eq.-1) then
         write (6,'(A)') 'Error. input matrix may be wrong.'
      elseif (ierr.eq.-2) then
         write (6,'(A)') 'matrix L overflow, increase elbow and retry'
      elseif (ierr.eq.-3) then
         write (6,'(A)') 'matrix U overflow, increase elbow and retry'
      elseif (ierr.eq.-4) then
         write (6,'(A)') 'Illegal value for lfil'
      elseif (ierr.eq.-5) then
         write (6,'(A)') 'zero row encountered'
      elseif (ierr.eq.-6) then
         write (6,'(A)') 'zero column encountered'
      elseif (ierr.eq.-7) then
         write (6,'(A)') 'buffers are too small'
      elseif (ierr.ne.0) then
         write (6,'(A,I3)') 'zero pivot encountered at step number',ierr
      endif
      if (ierr.ne.0) then
	print*,'ILUPACKFACTOR, location 1, found error, stopping...'
	stop
      end if

call gettim(hr,minu,sec,hundsec)
print*,'done calling ILUPACKfactor, location 1, at ', hr,' : ', minu, ' : ', sec, ' : ', hundsec
    
ndecomps = ndecomps + 1 
havrate = .false. 


!The main loop

done=.false.
at_hmin = .false. 

do while (.not.(done))

    hmin = 1.d-300
    absh = min(hmax, max(hmin, absh))
    if (absh == hmin) then
        if (at_hmin) then
            absh = abshlast   ! required by stepsize recovery
        endif
        at_hmin = .true.
    else
        at_hmin = .false. 
    end if
    h = absh 

    if (1.1d0*absh >= abs(tfinal - t)) then
    h = tfinal - t 
    absh = abs(h) 
    done = .true. 
    end if

    if ((absh /= abshlast) .or. (solverk /= klast)) then
 
        difRU=(dble(kI) - 1.d0 - dble(kJ)*(absh/abshlast)) / dble(kI)
        do mm=2,5
        difRU(mm,:)=difRU(mm-1,:)*difRU(mm,:)
        enddo
	if (allocated(mattemp)) deallocate(mattemp)
        allocate(mattemp(5,5))
        call dgemm('n','n',5,5,5,1.d0,difRU,5,difU,5,0.d0,mattemp,5)
        difRU = mattemp 
        deallocate(mattemp)

	if (allocated(mattemp)) deallocate(mattemp)
        allocate(mattemp(neq,solverk))
        call dgemm('n','n',neq,solverk,solverk,1.d0,dif(:,solverkcap),neq, &
            difRU(solverkcap,solverkcap),solverk,0.d0,mattemp,neq)
        dif(:,solverkcap) = mattemp
        deallocate(mattemp)

        hinvGak = h * invGa(solverk) 
        nconhk = 0 

	!reevaluate Miter using new dfdy or new hinvGak

        call mkl_dcsradd('n',0,idum,neq,neq,speye_csr,speye_ja,speye_ia,-hinvGak, &
    dfdy_csr,dfdy_ja,dfdy_ia,Miter_csr,dfdy_ja,dfdy_ia,dfdy_nnz,info)
	ib=dfdy_ia
	jb=dfdy_ja
	Miterb=Miter_csr
        if (info /= 0) then
        print*,'The following ERROR was detected by mkl_dcsradd (2nd): ', info
        STOP
        end if

	ierr=ILUPACKFACTOR(param,PREC,neq,ib,jb,Miterb,matching,ordering,droptol, &
	droptolS, condest, restol, ilumaxit, elbow, lfil, lfilS, nrestart)
        if (ierr.eq.-1) then
           write (6,'(A)') 'Error. input matrix may be wrong.'
        elseif (ierr.eq.-2) then
           write (6,'(A)') 'matrix L overflow, increase elbow and retry'
        elseif (ierr.eq.-3) then
           write (6,'(A)') 'matrix U overflow, increase elbow and retry'
        elseif (ierr.eq.-4) then
           write (6,'(A)') 'Illegal value for lfil'
        elseif (ierr.eq.-5) then
           write (6,'(A)') 'zero row encountered'
        elseif (ierr.eq.-6) then
           write (6,'(A)') 'zero column encountered'
        elseif (ierr.eq.-7) then
           write (6,'(A)') 'buffers are too small'
        elseif (ierr.ne.0) then
           write (6,'(A,I3)') 'zero pivot encountered at step number',ierr
        endif
        if (ierr.ne.0) then
  		print*,'ILUPACKFACTOR, location 2, found error, stopping...'
		stop
        end if

        ndecomps = ndecomps + 1 
        havrate = .false. 
    end if    

    ! LOOP FOR ADVANCING ONE STEP.
    nofailed = .true.                       ! no failed attempts
    do while (.true.)                            ! Evaluate the formula.    
        gotynew=.false.
        do while (.not.(gotynew))

            ! Compute the constant terms in the equation for ynew.
            call dgemv('n',neq,solverk,1.d0,dif(:,solverkcap),neq,G(solverkcap) * invGa(solverk),1,0.d0,psi,1)

            ! Predict a solution at t+h.
            tnew = t + h 
            if (done) then
                tnew = tfinal    ! Hit end point exactly.
            end if
            h = tnew - t       ! Purify h.
            pred = y + sum(dif(:,solverkcap),2) 
            ynew = pred 

            ! The difference, difkp1, between pred and the final accepted
            ! ynew is equal to the backward difference of ynew of order
            ! k+1. Initialize to zero for the iteration to compute ynew.
            difkp1 = 0.d0 
	    !$OMP PARALLEL PRIVATE (mm), DEFAULT(shared)
	    !$OMP DO
            do mm=1,neq
            invwt(mm) = 1.d0 / max(max(abs(y(mm)),abs(ynew(mm))),threshold)
            end do
	    !$OMP END DO
	    !$OMP END PARALLEL 
            minnrm = 100.d0*epsilon(1.d0)*maxval(abs(ynew * invwt)) 

            ! Iterate with simplified Newton method.
            tooslow = .false. 
            do iter = 1,maxit
                call xuode(neq,tnew,ynew,rhs)
                rhs = hinvGak*rhs -  (psi+difkp1)
		rhs=rhs+1.d-10
		rhs1=1.d-10
		call ILUPACKSOL (param,PREC,rhs,del,neq)
		call ILUPACKSOL (param,PREC,rhs1,del1,neq)
		del=del-del1
                newnrm = maxval(abs(del * invwt))
                difkp1 = difkp1 + del 
                ynew = pred + difkp1 

                if (newnrm <= minnrm) then
                    gotynew = .true.
                    exit
                else if (iter == 1) then
                    if (havrate) then
                        errit = newnrm * rate / (1.d0 - rate)
                        if (errit <= 0.05d0*rtol) then       ! More stringent when using old rate.
                            gotynew = .true.
                            exit
                        end if
                    else
                        rate = 0.d0
                    end if
                else if (newnrm > 0.9d0*oldnrm) then
                    tooslow = .true.
                    exit
                else
                    rate = max(0.9d0*rate, newnrm / oldnrm)
                    havrate = .true.
                    errit = newnrm * rate / (1.d0 - rate)
                    if (errit <= 0.5d0*rtol) then
                        gotynew = .true.
                        exit
                    else if (iter == maxit) then
                        tooslow = .true.
                        exit
                    else if (0.5d0*rtol < errit*rate**(maxit-iter)) then
                        tooslow = .true.
                        exit
                    end if
                end if

                oldnrm = newnrm
            end do                               ! end of Newton loop
            nfevals = nfevals + iter
            nsolves = nsolves + iter

            if (tooslow) then

                    nfailed = nfailed + 1
                    ! Speed up the iteration by forming new linearization or reducing h.
                    if ((.not.(Jcurrent)) .or. (.not.(Mcurrent))) then
                        if (.not.(Jcurrent)) then
                            call xujac(neq,t,y)   ! update dfdy sparse matrix
			deallocate(Miter_csr,ib,jb,Miterb)
			allocate(Miter_csr(dfdy_nnz),ib(neq+1),jb(dfdy_nnz),Miterb(dfdy_nnz))
                            npds = npds + 1
                            Jcurrent = .true.
			    ib=dfdy_ia
			    jb=dfdy_ja
                        end if
                    else if (absh <= hmin) then
                        print*,'ODE solver:IntegrationTolNotMet: ','Failure at t= ', t
                        print*,'Unable to meet integration tolerances without reducing &
                            the step size below the smallest value allowed: ',hmin,'at time t.'
                        print*,'absh now is: ', absh
                        STOP
                    else
                        abshlast = absh
                        absh = max(0.3d0 * absh, hmin)
                        h = absh
                        done = .false.

                        difRU=(dble(kI) - 1.d0 - dble(kJ)*(absh/abshlast)) / dble(kI)
                        do mm=2,5
                        difRU(mm,:)=difRU(mm-1,:)*difRU(mm,:)
                        enddo
			if (allocated(mattemp)) deallocate(mattemp)
                        allocate(mattemp(5,5))
                        call dgemm('n','n',5,5,5,1.d0,difRU,5,difU,5,0.d0,mattemp,5)
                        difRU = mattemp 
                        deallocate(mattemp)

			if (allocated(mattemp)) deallocate(mattemp)
                        allocate(mattemp(neq,solverk))
                        call dgemm('n','n',neq,solverk,solverk,1.d0,dif(:,solverkcap),neq, &
                            difRU(solverkcap,solverkcap),solverk,0.d0,mattemp,neq)
                        dif(:,solverkcap) = mattemp
                        deallocate(mattemp)

                        hinvGak = h * invGa(solverk)
                        nconhk = 0
                    end if


                    !reevaluate Miter using new dfdy or new hinvGak
                    call mkl_dcsradd('n',0,idum,neq,neq,speye_csr,speye_ja,speye_ia,-hinvGak, &
                    dfdy_csr,dfdy_ja,dfdy_ia,Miter_csr,dfdy_ja,dfdy_ia,dfdy_nnz,info)
			ib=dfdy_ia
			jb=dfdy_ja
			Miterb=Miter_csr

                    if (info /= 0) then
                        print*, 'The following ERROR was detected by mkl_dcsradd (3rd): ', info
                        STOP
                    end if

		ierr=ILUPACKFACTOR(param,PREC,&
			neq,ib,jb,Miterb,matching,ordering,droptol, &
			droptolS, condest, restol, ilumaxit, elbow, lfil, lfilS, nrestart)
		if (ierr.eq.-1) then
		   write (6,'(A)') 'Error. input matrix may be wrong.'
		elseif (ierr.eq.-2) then
		   write (6,'(A)') 'matrix L overflow, increase elbow and retry'
		elseif (ierr.eq.-3) then
		   write (6,'(A)') 'matrix U overflow, increase elbow and retry'
		elseif (ierr.eq.-4) then
		   write (6,'(A)') 'Illegal value for lfil'
		elseif (ierr.eq.-5) then
		   write (6,'(A)') 'zero row encountered'
		elseif (ierr.eq.-6) then
		   write (6,'(A)') 'zero column encountered'
		elseif (ierr.eq.-7) then
		   write (6,'(A)') 'buffers are too small'
		elseif (ierr.ne.0) then
		   write (6,'(A,I3)') 'zero pivot encountered at step number',ierr
		endif
		if (ierr.ne.0) then
	  		print*,'ILUPACKFACTOR, location 3 (tooslow), found error, stopping...'
			stop
		end if

                    ndecomps = ndecomps + 1 
                    havrate = .false. 
            end if
        end do     ! end of while loop for getting ynew
        err = maxval(abs(difkp1 * invwt)) * erconst(solverk)

        if (err > rtol) then                       ! Failed step
            nfailed = nfailed + 1
            if (absh <= hmin) then
                print*,'ODE solver:IntegrationTolNotMet: ','Failure at t= ', t
                print*,'Unable to meet integration tolerances without reducing &
                    the step size below the smallest value allowed: ',hmin,'at time t.'
                print*,'absh now is: ', absh
                STOP
            end if

            abshlast = absh
            if (nofailed) then
                nofailed = .false.
                hopt = absh * max(0.1d0, 0.833d0*(rtol/err)**(1.d0/(solverk+1.d0))) ! 1/1.2
                if (solverk > 1) then
                    errkm1 = maxval(abs((dif(:,solverk) + difkp1) * invwt)) * erconst(solverk-1)

                    hkm1 = absh * max(0.1d0, 0.769d0*(rtol/errkm1)**(1.d0/solverk)) ! 1/1.3
                    if (hkm1 > hopt) then
                        hopt = min(absh,hkm1)      ! don't allow step size increase
                        solverk = solverk - 1
                        if (allocated(solverkcap)) deallocate(solverkcap)
                        allocate(solverkcap(solverk))
                        solverkcap=(/(mm,mm=1,solverk)/)
                    end if
                end if
                absh = max(hmin, hopt)
            else
                absh = max(hmin, 0.5d0 * absh)
            end if
            h = absh
            if (absh < abshlast) then
                done = .false.
            end if

            difRU=(dble(kI) - 1.d0 - dble(kJ)*(absh/abshlast)) / dble(kI)
            do mm=2,5
            difRU(mm,:)=difRU(mm-1,:)*difRU(mm,:)
            enddo
	    if (allocated(mattemp)) deallocate(mattemp)
            allocate(mattemp(5,5))
            call dgemm('n','n',5,5,5,1.d0,difRU,5,difU,5,0.d0,mattemp,5)
            difRU = mattemp 
            deallocate(mattemp)

	    if (allocated(mattemp)) deallocate(mattemp)
            allocate(mattemp(neq,solverk))
            call dgemm('n','n',neq,solverk,solverk,1.d0,dif(:,solverkcap),neq, &
                difRU(solverkcap,solverkcap),solverk,0.d0,mattemp,neq)
            dif(:,solverkcap) = mattemp
            deallocate(mattemp)

            hinvGak = h * invGa(solverk)
            nconhk = 0

            !reevaluate Miter using new dfdy or new hinvGak
            call mkl_dcsradd('n',0,idum,neq,neq,speye_csr,speye_ja,speye_ia,-hinvGak, &
            dfdy_csr,dfdy_ja,dfdy_ia,Miter_csr,dfdy_ja,dfdy_ia,dfdy_nnz,info)
		ib=dfdy_ia
		jb=dfdy_ja
		Miterb=Miter_csr

            if (info /= 0) then
                print*, 'The following ERROR was detected by mkl_dcsradd (4th): ', info
                STOP
            end if
		ierr=ILUPACKFACTOR(param,PREC,&
			neq,ib,jb,Miterb,matching,ordering,droptol, &
			droptolS, condest, restol, ilumaxit, elbow, lfil, lfilS, nrestart)
		if (ierr.eq.-1) then
		   write (6,'(A)') 'Error. input matrix may be wrong.'
		elseif (ierr.eq.-2) then
		   write (6,'(A)') 'matrix L overflow, increase elbow and retry'
		elseif (ierr.eq.-3) then
		   write (6,'(A)') 'matrix U overflow, increase elbow and retry'
		elseif (ierr.eq.-4) then
		   write (6,'(A)') 'Illegal value for lfil'
		elseif (ierr.eq.-5) then
		   write (6,'(A)') 'zero row encountered'
		elseif (ierr.eq.-6) then
		   write (6,'(A)') 'zero column encountered'
		elseif (ierr.eq.-7) then
		   write (6,'(A)') 'buffers are too small'
		elseif (ierr.ne.0) then
		   write (6,'(A,I3)') 'zero pivot encountered at step number',ierr
		endif
		if (ierr.ne.0) then
	  		print*,'ILUPACKFACTOR, location 4(failed step), found error, stopping...'
			stop
		end if

            ndecomps = ndecomps + 1
            havrate = .false.
        else                                ! Successful step
            exit

        end if
    end do   ! while true

!print*, 'done while true'

    nsteps = nsteps + 1

    dif(:,solverk+2) = difkp1 - dif(:,solverk+1)
    dif(:,solverk+1) = difkp1
    do mm = solverk,1,-1
        dif(:,mm) = dif(:,mm) + dif(:,mm+1)
    end do

    NNreset_dif = .false.

    nout_new =  0
    do while (next <= ntspan)
        if ((tnew - tspan(next)) < 0) then
            exit
        end if
        nout_new = nout_new + 1 
        next = next + 1 
    end do
print*,'tnew now is: ', tnew
    if (nout_new > 0) then

	if (allocated(tout_new)) deallocate(tout_new,out_s,yout_new)
        allocate(tout_new(nout_new),out_s(nout_new),yout_new(neq,nout_new))

        tout_new = tspan(next-nout_new : next-1) 
        out_s=(tout_new-tnew)/h
        yout_new=spread(ynew,2,nout_new)

        if (solverk == 1) then
            call dgemm('n','n',neq,nout_new,1,1.d0,dif(:,1),neq,reshape(out_s,(/1,nout_new/)),1,1.d0,yout_new,neq)
        else                    ! cumprod collapses vectors
	    if (allocated(out_kI)) deallocate(out_kI)
	    if (allocated(mattemp)) deallocate(mattemp)
            allocate(out_kI(solverk,nout_new),mattemp(solverk,nout_new))
            out_kI=spread(solverkcap,2,nout_new)
            mattemp=(spread(out_s,1,solverk)+dble(out_kI)-1.d0)/dble(out_kI)
            do mm=2,solverk
                mattemp(mm,:)=mattemp(mm-1,:)*mattemp(mm,:)
            end do
            call dgemm('n','n',neq,nout_new,solverk,1.d0,dif(:,solverkcap),neq,mattemp,solverk,1.d0,yout_new,neq)
            deallocate(out_kI,mattemp)
        end if

        call xuoutput(tout_new,nout_new,yout_new)
        print*, 'now tout_new is : ', tout_new(nout_new)	
        deallocate(tout_new,out_s,yout_new)
    end if
	

    if (done) then
	print*,'done, exiting...'
        exit
    end if

    klast = solverk
    abshlast = absh 
    nconhk = min(nconhk+1,maxk+2) 
    if (nconhk >= solverk + 2) then
        temp = 1.2d0*(err/rtol)**(1.d0/(solverk+1.d0)) 
        if (temp > 0.1d0) then
            hopt = absh / temp 
        else
            hopt = 10.d0*absh 
        end if
        kopt = solverk 
        if (solverk > 1) then
                errkm1 = maxval(abs(dif(:,solverk) * invwt)) * erconst(solverk-1) 
            
            temp = 1.3d0*(errkm1/rtol)**(1.d0/solverk) 
            if (temp > 0.1d0) then
                hkm1 = absh / temp 
            else
                hkm1 = 10.d0*absh 
            end if
            if (hkm1 > hopt) then
                hopt = hkm1 
                kopt = solverk - 1 
            end if
        end if
        if (solverk < maxk) then
                errkp1 = maxval(abs(dif(:,solverk+2) * invwt)) * erconst(solverk+1) 

            temp = 1.4d0*(errkp1/rtol)**(1.d0/(solverk+2.d0)) 
            if (temp > 0.1) then
                hkp1 = absh / temp 
            else
                hkp1 = 10.d0*absh 
            end if
            if (hkp1 > hopt) then
                hopt = hkp1 
                kopt = solverk + 1 
            end if
        end if
        if (hopt > absh) then
            absh = hopt 
            if (solverk /= kopt) then
                solverk = kopt 
                if (allocated(solverkcap)) deallocate(solverkcap)
                allocate(solverkcap(solverk))
                solverkcap = (/(mm,mm=1,solverk)/) 
            end if
        end if
    end if

     ! Advance the integration one step.
    t = tnew 
    y = ynew 
    Jcurrent = .false.
    Mcurrent = .true.                    ! Constant mass matrix I or M.

end do ! while ~done

call ILUPACKDELETE(param,PREC)

111 print*, 'tnew now is: ', tnew
print*,'nsteps: ',nsteps, '		nfailed: ',nfailed
print*,'nfevals: ',nfevals,'			npds: ',npds
print*,'ndecomps: ',ndecomps,'			nsolves: ',nsolves

END SUBROUTINE xusolver
        
!********************************************************************************
!********************************************************************************
