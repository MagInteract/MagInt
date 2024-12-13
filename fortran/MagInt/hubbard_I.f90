!----------------------------------------------------------------------------------------!
! MagInteract library
! Written by  : Leonid V. Pourovskii (CPHT Ecole Polytechnique) 2012-2022
! Email: leonid@cpht.polytechnique.fr
! ---------------------------------------------------------------------------------------!
!"" Implements the force-theorem Hubbard-I (FT-HI) approach to intersite
!   exchange interactions in correlated insulators.
!   The formalism is given in L. V. Pourovskii Phys. Rev. B 94, 115117 (2016)
!   This Python-3 version is based on TRIQS library """
!
! ---------------------------------------------------------------------------------------!
! Hubbard-I version of MagInteract.
! Includes routins for reading and manipulating atomic eigenstates as well as
! Hubbard-I solver modified for MagInteract as described in the main gf_HI_fullU_INT
! subroutine
! ---------------------------------------------------------------------------------------!
MODULE hubbard_I_data
  !use, intrinsic :: iso_c_binding
  REAL(8) :: B
  type :: occup
     LOGICAL :: ifdiag, run
     INTEGER :: ndeg
     INTEGER :: n
     !COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: Hn
     !REAL(8), DIMENSION(:), ALLOCATABLE :: En
     !INTEGER, DIMENSION(:), ALLOCATABLE  :: st_n
     !INTEGER, DIMENSION(:), ALLOCATABLE  :: narr
     COMPLEX(8), pointer :: Hn(:,:)
     REAL(8), pointer :: En(:)
     INTEGER, pointer  :: st_n(:)
     INTEGER, pointer  :: narr(:)
     end type occup
END MODULE hubbard_I_data


MODULE states_GN
  REAL(KIND=8), PARAMETER :: tol_mel=1d-6
  INTEGER, allocatable  :: st_gn(:)
  INTEGER, allocatable  :: narr_g(:)
  real(8), allocatable :: occ(:)
  integer, allocatable :: arr(:)
  INTERFACE OP
       MODULE PROCEDURE OP_1el, OP_2el
  ENDINTERFACE OP
!

 CONTAINS

 FUNCTION c_dag_c_dag_cc(vec,m,m1,m2,m3)
  ! calculates c^dag_m c^dag_m1 c_m2 c_m3 |vec>
  !USE states_GN
  IMPLICIT NONE
  COMPLEX(KIND=8), INTENT(IN) :: vec(:)
  INTEGER,  INTENT(IN) :: m,m1,m2,m3
  COMPLEX(KIND=8) :: c_dag_c_dag_cc(size(vec))
 !
  INTEGER :: N, nso
  INTEGER :: k, l, ks, ls, p, ps, q, qs, k1, ind, i
  REAL(KIND=8) :: fsign, fsign1
 ! calc c_m|vec>
  N=size(vec)
  nso=size(arr)
  c_dag_c_dag_cc=0d0
  DO i=1,N
     occ = merge(1.d0,0.d0,btest(st_gn(i),arr(0:nso-1)))
     k=ibclr(st_gn(i),m2-1)
     IF(k==st_gn(i)) CYCLE
     occ(m2)=0
     ks=SUM(occ(1:m2-1))
     l=ibclr(k,m3-1)
     IF(l==k) CYCLE
     occ(m3)=0
     ls=SUM(occ(1:m3-1))
     p=ibset(l,m1-1)
     IF(p==l) CYCLE
     occ(m1)=1
     ps=SUM(occ(1:m1-1))
     q=ibset(p,m-1)
     IF(q==p) CYCLE
     qs=SUM(occ(1:m-1))
     fsign=(-1d0)**(ks+ls+ps+qs)
     ind=narr_g(q)
     if (ind<1.OR.ind>N) THEN
         WRITE(*,*)'ERROR in c_dag_c_dag_cc', ind,p,narr_g(p)
         STOP
     ENDIF
     c_dag_c_dag_cc(ind)=c_dag_c_dag_cc(ind)+vec(i)*fsign
  ENDDO
  RETURN
 END FUNCTION c_dag_c_dag_cc

 FUNCTION c_dag_c(vec,m,m1)
  ! calculates c^dag_m c_m1 |vec>
  !USE states_GN
  IMPLICIT NONE
  COMPLEX(KIND=8), INTENT(IN) :: vec(:)
  INTEGER,  INTENT(IN) :: m,m1
  COMPLEX(KIND=8) :: c_dag_c(size(vec))
 !
  INTEGER :: N, nso
  INTEGER :: k, l, ls, p, ps, k1, ind
  REAL(KIND=8) :: fsign, fsign1
 ! calc c_m|vec>
  N=size(vec)
  nso=size(arr)
  c_dag_c=0d0
  DO k=1,N
     occ = merge(1.d0,0.d0,btest(st_gn(k),arr(0:nso-1)))
     l=ibclr(st_gn(k),m1-1)
     IF(l==st_gn(k)) CYCLE
     occ(m1)=0
     ls=SUM(occ(1:m1-1))
     fsign=(-1d0)**ls
     p=ibset(l,m-1)
     IF(p==l) CYCLE
     ps=SUM(occ(1:m-1))
     fsign1=(-1d0)**ps
     ind=narr_g(p)
     if (ind<1.OR.ind>N) THEN
         WRITE(*,*)'ERROR in c_dag_c', ind,p,narr_g(p)
         WRITE(*,*)'N= ',N
         STOP
     ENDIF
     c_dag_c(ind)=c_dag_c(ind)+vec(k)*fsign*fsign1
  ENDDO
  RETURN
 END FUNCTION c_dag_c

 FUNCTION OP_2el(vec,mat,nso,n)
 ! calculates OP_el|vec>, where OP_2el is
 ! Sum_{m,m1,m2,m3} mat_m_m1_m3_m2 c^dag_m c^dag_m1 c_m2 c_m3
 !USE states_GN
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nso, n
  COMPLEX(KIND=8), INTENT(IN) :: vec(n)
  COMPLEX(KIND=8), INTENT(IN) :: mat(nso,nso,nso,nso)
  COMPLEX(KIND=8) :: OP_2el(n)
  REAL(KIND=8) :: dd
  INTEGER :: m, m1, m2, m3
  OP_2el=0d0
  !nso=size(arr)
  DO m=1,nso
     DO m1=1,nso
        DO m2=1,nso
           DO m3=1,nso
             dd = SQRT(REAL(mat(m,m1,m2,m3)*CONJG(mat(m,m1,m2,m3)),KIND=8))
             IF (dd > tol_mel) THEN
               OP_2el=OP_2el+c_dag_c_dag_cc(vec,m,m1,m2,m3)*mat(m,m1,m2,m3)
             ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  RETURN
 END FUNCTION OP_2el
!
 FUNCTION OP_1el(vec,mat,nso,n)
  ! calculates OP|vec>, where OP is Sum_{m,m1} mat_m_m1 c^dag_m c_m1
  !USE states_GN
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nso, n
  COMPLEX(KIND=8), INTENT(IN) :: vec(n)
  COMPLEX(KIND=8), INTENT(IN) :: mat(nso,nso)
  COMPLEX(KIND=8) :: OP_1el(n)
  INTEGER :: m, m1
  OP_1el=0d0
  DO m=1,nso
     DO m1=1,nso
        IF (ABS(mat(m,m1)) > tol_mel) THEN
               OP_1el=OP_1el+c_dag_c(vec,m,m1)*mat(m,m1)
        ENDIF
     ENDDO
  ENDDO
  RETURN
 END FUNCTION OP_1el

 FUNCTION ovl(vec,vec1)
  ! calculates <vec|vec1>
  !USE states_GN
  IMPLICIT NONE
  COMPLEX(KIND=8), INTENT(IN) :: vec(:), vec1(:)
  COMPLEX(KIND=8) :: ovl
  INTEGER :: i, N
  N=size(vec)
  ovl=0d0
  DO i=1,N
     ovl=ovl+CONJG(vec(i))*vec1(i)

  ENDDO
  RETURN
 END FUNCTION ovl

 SUBROUTINE setup_aux_arrays(Nat,nso,n)
  !USE states_GN
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Nat, nso
  INTEGER, INTENT(OUT) :: n
  INTEGER :: i, k, nstate
  REAL(KIND=8), EXTERNAL :: factor

  n=NINT(factor(nso)/factor(nso-Nat)/factor(Nat))
  ! Fill up array of states
  nstate=2**nso
  ALLOCATE(st_gn(n),arr(0:nso-1),occ(nso),narr_g(0:nstate-1))
  FORALL( i = 0:nso-1 ) arr(i) = i
  k=0
  narr_g=-1
  DO i=0,nstate-1
    occ = merge( 1.d0, 0.d0, btest( i, arr(0:nso-1) ) )
    IF(SUM(occ)==Nat) THEN
       k=k+1
       st_gn(k)=i
       narr_g(i)=k
    ENDIF
  ENDDO
  RETURN
  END SUBROUTINE setup_aux_arrays


 SUBROUTINE free_aux_arrays
  !USE states_GN
  IMPLICIT NONE
  DEALLOCATE(st_gn,arr,occ,narr_g)
  RETURN
 END SUBROUTINE free_aux_arrays

END MODULE states_GN


SUBROUTINE OP_exp_val(OP_mat,Nat,expvals,nvec,nso)
  USE hubbard_I_data
  USE states_GN
  IMPLICIT NONE
  INTEGER,  INTENT(IN) :: Nat, nso, nvec
  COMPLEX(KIND=8), INTENT(IN) :: OP_mat(nso,nso)
  COMPLEX(KIND=8), INTENT(OUT) :: expvals(nvec)

  INTEGER :: i, ivec, n, k, nstate
  COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: vecs
  COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: tmp_vec


  CALL setup_aux_arrays(Nat,nso,n)
  ALLOCATE(tmp_vec(n),vecs(n,nvec))
  OPEN(111,file='STATES',status='old')
  DO i=1,nvec
     DO k=1,n
        READ(111,*)vecs(k,i)
     ENDDO
     tmp_vec=OP(vecs(:,i),OP_mat,nso,n)
     expvals(i)=ovl(vecs(:,i),tmp_vec)
  ENDDO
  CLOSE(111)
  DEALLOCATE(vecs,tmp_vec)
  CALL free_aux_arrays
  RETURN
END SUBROUTINE OP_exp_val

SUBROUTINE OP_exp_val_invec(resvec,expval,OP_mat,invec,Nat,nso,n)
!
! applies OP_mat to the vector invec
! returns the expectation nalue<invec|OP|invec>
! and normalized OP|invec>
!
  USE hubbard_I_data
  USE states_GN
  IMPLICIT NONE
  INTEGER,  INTENT(IN) :: Nat,nso,n
  COMPLEX(KIND=8), INTENT(IN) :: OP_mat(nso,nso)
  COMPLEX(KIND=8), INTENT(IN) :: invec(n)
  COMPLEX(KIND=8), INTENT(out) :: expval
  COMPLEX(KIND=8), INTENT(out) :: resvec(n)

  COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: tmp_vec
  COMPLEX(KIND=8) :: norm

  INTEGER :: nn


  CALL setup_aux_arrays(Nat,nso,nn)
  IF (nn.NE.n) THEN
      WRITE(*,*)'WRONG dimension nso,Nat,n ',nso,Nat,n,' in OP_on_state'
      STOP
  ENDIF
  resvec=OP(invec,OP_mat,nso,nn)
  expval=ovl(invec,resvec)
  norm=ovl(resvec,resvec)
  resvec=resvec/SQRT(norm)
  CALL free_aux_arrays
  RETURN
END SUBROUTINE OP_exp_val_invec


SUBROUTINE prn_vectors(Nat,nso,n,ST_num,tresh)
! print fock states and coeeficients contributing to
! ST_num first vectors in STATES
! above treshold tresh
USE states_GN
IMPLICIT NONE
INTEGER,  INTENT(IN) :: Nat, nso, n, ST_num
REAL(KIND=8),  INTENT(IN) :: tresh

INTEGER :: i, k, inum, nn
COMPLEX(KIND=8):: coff

CALL setup_aux_arrays(Nat,nso,nn)
IF (nn.NE.n) THEN
    WRITE(*,*)'WRONG dimension nso,Nat,n ',nso,Nat,n,' in OPs_on_state'
    STOP
ENDIF
! Read vector number ST_num
OPEN(111,file='STATES',status='old')
OPEN(112,file='STATES_main_Fock_components')
DO i=1,ST_num
   WRITE(112,'(/,a,i5)')'VECTOR : ',i
   DO k=1,n
      READ(111,*)coff
      IF(ABS(coff) > tresh) THEN
         inum=st_gn(k)
         occ = merge( 1.d0, 0.d0, btest( inum, arr(0:nso-1) ) )
         WRITE(112,'(2F10.4,5x,14I1)')REAL(coff,KIND=8),AIMAG(coff),NINT(occ)
      ENDIF
   ENDDO
ENDDO
CLOSE(111)
CLOSE(112)
CALL free_aux_arrays
RETURN
END SUBROUTINE prn_vectors

SUBROUTINE prn_input_vector(Nat,nso,n,vec,tresh)
! print fock states and coeeficients contributing to
! ST_num first vectors in STATES
! above treshold tresh
USE states_GN
IMPLICIT NONE
INTEGER,  INTENT(IN) :: Nat, nso, n
REAL(KIND=8),  INTENT(IN) :: tresh
COMPLEX(KIND=8), INTENT(IN) :: vec(n)

INTEGER :: i, k, inum, nn
COMPLEX(KIND=8):: coff

CALL setup_aux_arrays(Nat,nso,nn)
IF (nn.NE.n) THEN
    WRITE(*,*)'WRONG dimension nso,Nat,n ',nso,Nat,n,' in OPs_on_state'
    STOP
ENDIF
OPEN(112,file='Vector_Fock_components')
DO k=1,n
   IF(ABS(vec(k)) > tresh) THEN
      inum=st_gn(k)
      occ = merge( 1.d0, 0.d0, btest( inum, arr(0:nso-1) ) )
      WRITE(112,'(2F10.4,5x,14I1)')REAL(vec(k),KIND=8),AIMAG(vec(k)),NINT(occ)
   ENDIF
ENDDO
CLOSE(112)
CALL free_aux_arrays
RETURN
END SUBROUTINE prn_input_vector

SUBROUTINE get_vec(n,iv,vec)
IMPLICIT NONE
!
! Read and return vector number iv from file STATES
!
INTEGER,  INTENT(IN) :: n, iv
COMPLEX(KIND=8), INTENT(OUT) :: vec(n)
INTEGER :: k, i
INTEGER :: IO_STATUS

!write(*,*)' Vector lenght in Fortan:'
!write(*,*)n

OPEN(111,file='STATES',status='old')
DO i=1,n
   DO k=1,n
      READ(111,*,IOSTAT=IO_STATUS) vec(k)
      IF (IO_STATUS /= 0) THEN
          PRINT *, 'Error reading from file'
      END IF
   ENDDO
   IF (i==iv) EXIT
ENDDO
CLOSE(111)
RETURN
END SUBROUTINE get_vec


SUBROUTINE mat_el_vecs(Nat,nso,n,ST_l,ST_r,val,OP_mat)
!
! computes matrix element <ST_l|OP_mat|ST_r>,
! where OP_mat is a one-electron operator matrix
!
USE states_GN
IMPLICIT NONE
INTEGER,  INTENT(IN) :: Nat, nso, n
COMPLEX(KIND=8), INTENT(IN) :: ST_l(n), ST_r(n)
COMPLEX(KIND=8), INTENT(IN) :: OP_mat(nso,nso)
COMPLEX(KIND=8), INTENT(out):: val

INTEGER :: nn
COMPLEX(KIND=8) :: ST_res(n)

CALL setup_aux_arrays(Nat,nso,nn)
IF (nn.NE.n) THEN
    WRITE(*,*)'WRONG dimension nso,Nat,n ',nso,Nat,n,' in OPs_on_state'
    STOP
ENDIF

!WRITE(*,'(/,a)')'OP mat'
!DO i=1,nso
!   write(*,'(10F9.4)')(REAL(OP_mat(i,j)),j=1,nso)
!ENDDO

ST_res=OP(ST_r,OP_mat,nso,nn)
val=ovl(ST_l,ST_res)
CALL free_aux_arrays
RETURN
END SUBROUTINE mat_el_vecs


SUBROUTINE mat_el_vecs_2el(Nat,nso,n,ST_l,ST_r,val,OP_mat)
!
! computes matrix element <ST_l|OP_mat|ST_r>,
! where OP_mat is a two-electron operator matrix
!
USE states_GN
IMPLICIT NONE
INTEGER,  INTENT(IN) :: Nat, nso, n
COMPLEX(KIND=8), INTENT(IN) :: ST_l(n), ST_r(n)
COMPLEX(KIND=8), INTENT(IN) :: OP_mat(nso,nso,nso,nso)
COMPLEX(KIND=8), INTENT(out):: val

INTEGER :: nn
COMPLEX(KIND=8) :: ST_res(n)

CALL setup_aux_arrays(Nat,nso,nn)
IF (nn.NE.n) THEN
    WRITE(*,*)'WRONG dimension nso,Nat,n ',nso,Nat,n,' in OPs_on_state'
    STOP
ENDIF

ST_res=OP(ST_r,OP_mat,nso,nn)
val=ovl(ST_l,ST_res)
CALL free_aux_arrays
RETURN
END SUBROUTINE mat_el_vecs_2el

SUBROUTINE OP_on_state(Nat,nso,n,ST_num,ST_res,norm,OP_mat)
!
! Applies operator OP to the state number ST_num in the file states
! Returns normalized result + its norm= <st_init OP*| OP st_init>
!
USE states_GN
IMPLICIT NONE
INTEGER,  INTENT(IN) :: Nat, nso, n, ST_num
COMPLEX(KIND=8), INTENT(IN) :: OP_mat(nso,nso)
COMPLEX(KIND=8), INTENT(OUT) :: ST_res(n)
REAL(KIND=8), INTENT(OUT) :: norm

INTEGER :: nn, i, k
COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: ST_init

CALL setup_aux_arrays(Nat,nso,nn)
IF (nn.NE.n) THEN
    WRITE(*,*)'WRONG dimension nso,Nat,n ',nso,Nat,n,' in OPs_on_state'
    STOP
ENDIF
ALLOCATE(ST_init(n))

! Read vector number ST_num
OPEN(111,file='STATES',status='old')
DO i=1,ST_num
   DO k=1,n
      READ(111,*)ST_init(k)
   ENDDO
ENDDO
CLOSE(111)

ST_res=OP(ST_init,OP_mat,nso,nn)
norm=REAL(ovl(ST_res,ST_res),KIND=8)
ST_res=ST_res/SQRT(norm)
CALL free_aux_arrays
DEALLOCATE(ST_init)
RETURN
END SUBROUTINE OP_on_state


SUBROUTINE OPs_on_state(ST_res,ST_init,OP_mat,Nat,nso,n,Ntimes)
!
! Applies Ntimes times operator OP to the state ST_init
!
USE states_GN
IMPLICIT NONE
INTEGER,  INTENT(IN) :: Nat, nso, n, Ntimes
COMPLEX(KIND=8), INTENT(IN) :: ST_init(n)
COMPLEX(KIND=8), INTENT(IN) :: OP_mat(nso,nso)
COMPLEX(KIND=8), DIMENSION(n,Ntimes+1), INTENT(OUT) :: ST_res

INTEGER :: i, k, isize, nstate, nn
REAL(KIND=8), EXTERNAL :: factor


CALL setup_aux_arrays(Nat,nso,nn)
IF (nn.NE.n) THEN
    WRITE(*,*)'WRONG dimension nso,Nat,n ',nso,Nat,n,' in OPs_on_state'
    STOP
ENDIF
!APPLY the operator
ST_res(:,1)=ST_init
DO i=1,Ntimes
   ST_res(:,i+1)=OP(ST_res(:,i),OP_mat,nso,nn)
ENDDO
CALL free_aux_arrays
RETURN
END SUBROUTINE OPs_on_state


SUBROUTINE fock_space_rmat(fs_rmat,RM,Nat,nso,n)
!
! calculates rotation matrix for vectors in Fock space for occupancy Nat
! from 1-electron rotation matrix RM
!
USE states_GN
IMPLICIT NONE
INTEGER,  INTENT(IN) :: Nat, nso, n
COMPLEX(KIND=8), INTENT(IN) :: RM(nso,nso)
COMPLEX(KIND=8), INTENT(OUT):: fs_rmat(n,n)
!COMPLEX(KIND=8), DIMENSION(n,n), INTENT(OUT):: fs_rmat
!
INTEGER(KIND=4), PARAMETER :: max_st=1000000
INTEGER :: i, k, m, m1, ind, ind1, n_st, j, n_st_in, n_st_out, nn, ind_in
INTEGER(KIND=8) :: p
INTEGER :: crops(nso,nso),croin(nso),indmax(nso)
COMPLEX(KIND=8) :: coff(nso,nso)
INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: st_arr_in,st_arr_out
COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: st_coff_in,st_coff_out

CALL setup_aux_arrays(Nat,nso,nn)
IF (nn.NE.n) THEN
    WRITE(*,*)'WRONG dimension nso,Nat,n ',nso,Nat,n,' in OPs_on_state'
    STOP
ENDIF

ALLOCATE(st_arr_in(max_st),st_arr_out(max_st),st_coff_in(max_st),st_coff_out(max_st))

fs_rmat=0d0

DO i=1,n
    occ = merge(1.d0,0.d0,btest(st_gn(i),arr(0:nso-1)))
    crops=0
    coff=0d0
    indmax=0
    ind=0
    DO m=1,nso
       p=ibclr(st_gn(i),m-1)
       IF (p==st_gn(i)) CYCLE
       ind=ind+1
       ind1=0
       DO m1=1,nso
          IF (abs(RM(m,m1)) > tol_mel) THEN
              ind1=ind1+1
              coff(ind,ind1)=RM(m,m1)
              crops(ind,ind1)=m1
          ENDIF
       ENDDO
       indmax(ind)=ind1
    ENDDO
    st_arr_in=0
    st_coff_in=0d0
    ! initial vacuum vector
    n_st_in=1
    st_coff_in(1)=1.0
    ind_in=0
    ! apply superpositions of creation operators to create expansion of vector i
    ! into rotated Fock space
    DO k=1,ind
       DO j=1,indmax(k)
           m=crops(k,j)
           CALL c_dag_to_set(st_arr_in,st_coff_in,st_arr_out,st_coff_out,n_st_in,n_st_out,ind_in,m,nso,max_st)
           st_coff_out(ind_in+1:n_st_out)=st_coff_out(ind_in+1:n_st_out)*coff(k,j)
           ind_in=n_st_out
       ENDDO
       st_coff_in=0d0
       st_arr_in=0
       st_coff_in(1:n_st_out)=st_coff_out(1:n_st_out)
       st_arr_in(1:n_st_out)=st_arr_out(1:n_st_out)
       n_st_in=n_st_out
       ind_in=0
    ENDDO
    DO k=1,n_st_in
       ind=narr_g(st_arr_in(k))
       IF (ind<1.OR.ind>n) THEN
         WRITE(*,*)'ERROR in fock_space_rmat',i,k,st_arr_out(k),ind,narr_g(st_arr_out(k))
         STOP
       ENDIF
       fs_rmat(i,ind)=fs_rmat(i,ind)+st_coff_in(k)
    ENDDO
ENDDO

DEALLOCATE(st_arr_in,st_arr_out,st_coff_in,st_coff_out)
CALL free_aux_arrays

RETURN
END SUBROUTINE fock_space_rmat

SUBROUTINE c_dag_to_set(vecs_in,coff_in,vecs_out,coff_out,n_in,n_out,ind_in,m,nso,max_st)
!
! applies creation operator c^+_m to each vector in the set vecs_in
! vectors are supposed to be multiplied by coefficients stored in coff_in
!
! the resulting vectors are stored in vecs_out and their prefactors in coff_out
! starting from index ind_in+1

  USE states_GN
  IMPLICIT NONE
  INTEGER(KIND=4):: max_st
  INTEGER,  INTENT(IN) :: n_in, ind_in, m, nso
  INTEGER(KIND=8),  INTENT(IN) :: vecs_in(max_st)
  COMPLEX(KIND=8), INTENT(IN) :: coff_in(max_st)
  INTEGER(KIND=8),  INTENT(INOUT) :: vecs_out(max_st)
  COMPLEX(KIND=8), INTENT(INOUT) :: coff_out(max_st)
  INTEGER,  INTENT(OUT):: n_out

  INTEGER :: i, j, k, l, ps, ind
  INTEGER(KIND=8) :: p
  REAL(KIND=8) :: fsign

  ind=0
  !WRITE(*,*)'n_in, ind_in, m, nso, max_st ',n_in, ind_in, m, nso, max_st
  DO i=1,n_in
     occ = merge(1.d0,0.d0,btest(vecs_in(i),arr(0:nso-1)))
     !WRITE(*,*)'occ: ',occ
     p=ibset(vecs_in(i),m-1)
     IF(p==vecs_in(i)) CYCLE
     ind=ind+1
     ps=SUM(occ(1:m-1))
     fsign=(-1d0)**(ps)
     vecs_out(ind_in+ind)=p
     coff_out(ind_in+ind)=coff_in(i)*fsign
     !occ = merge(1.d0,0.d0,btest(p,arr(0:nso-1)))
     !WRITE(*,*)'occ-new: ',occ
  ENDDO
  n_out=ind+ind_in
  !WRITE(*,*)'c_dag_to_set ',vecs_out(1:n_out)
END SUBROUTINE c_dag_to_set


SUBROUTINE gf_HI_fullU_INT(e0f,U,ummss,zmsb,nlm,Iwmax,nmom,ns,atocc,atmag,temp,verbosity,N_lev,Z,calc_off_diag, &
                       remove_CF,LadBS,LadOP,CalcOvl,Nbas,Mnat,StBas,GF0,Tail0,GF,Tail, OvlMat)

! Computes atomic GF with the full 4-index U 8.10.2007
! Finite temperature version
! (by L.V. Pourovskii)
!
!  Last change 09.12.2016
!
! Included tail calculation for the triqs package
! M. Aichhorn 10-2009
!
! GF0 /output/ - the GF for T=0, Tail0 is its tail
! GF /output/ - atomic Green's function
! e0f /input/ - atomic level position e0f_mm' C*_m C_m'
! U  /input/ - full 4-index U (orbitals and spins)
! ummss  /input/ - full 2-index U (orbitals and spins)
! zmsb(Iwmax)/input/ - complex energy mesh
! nlm, ns /input/ - orbital and spin degeneracy
! atocc, atmag /output/ - occupancy and magnetic moment of the atom
! temp /input/ - temperature
! verbosity - 0: no text output, 1: basics, 2: all
! N_lev: the degeneracy of the lowest-energy CF multiplet
! Z /output/ : statistical sun
! remove_CF/input/ suppress CF splitting between the levels of the lowest-energy CF multiplet
!
! New 30.01.2012
! Calculates atomic GF corresponding to each atomic level of the lowest-energy CF multiplet taken as
! the ground state. The output GF is an array of all this GF. GF0 is the usual GF computed for T=0
! NEW 10.07.2012
! calc_off_diag: calculate off-diagonal elements of GF=Four.Trans[-<A|T[f_m(tau)f^{dag}]|B>
! GF/Tail are matrices [N_lev,N_lev] of atomic levels

! New 22.12.2014
! LadBS/input/ if True regenerate N_lev states starting from state 1 by applying LadOP
! LadOP/input/ matrix of the corresponding ladder operator in the e0f basis
! CalcOvl/input/ if True calculate overlap of N_lev first states with N_lev states of "standard basis"
! StBas/input(CalcOvl=True)/output(LadBS=True)/ standard basis of states
! NBas/input/ number of states in the standard basis
! OvlMat/output/ [Nlev,Nlev] overlap matrix


  USE hubbard_I_data
  USE states_GN
  IMPLICIT NONE

! Input/output variables
  INTEGER, INTENT(in) :: nlm, ns, Iwmax, nmom, Nbas, Mnat
  COMPLEX(KIND=8), INTENT(in) :: e0f(nlm*ns,nlm*ns)
  COMPLEX(KIND=8), INTENT(in) :: zmsb(Iwmax)
  REAL(KIND=8), INTENT(in) :: ummss(nlm*ns,nlm*ns)
  REAL(KIND=8), INTENT(in) :: U(nlm*ns,nlm*ns,nlm*ns,nlm*ns)
  REAL(KIND=8), INTENT(in) :: temp
  integer, intent(in) :: verbosity
  integer, intent(in) :: N_lev
  LOGICAL, intent(in) :: calc_off_diag, remove_CF, LadBS, CalcOvl
  COMPLEX(KIND=8), INTENT(in) :: LadOP(nlm*ns,nlm*ns)
  COMPLEX(KIND=8), INTENT(inout) :: StBas(Mnat,Nbas)


  COMPLEX(KIND=8), INTENT(out) :: GF(N_lev,N_lev,nlm*ns,nlm*ns,Iwmax)
  Complex(kind=8), intent(out) :: Tail(N_lev,N_lev,nmom,nlm*ns,nlm*ns)
  COMPLEX(KIND=8), INTENT(out) :: GF0(nlm*ns,nlm*ns,Iwmax)
  Complex(kind=8), intent(out) :: Tail0(nmom,nlm*ns,nlm*ns)
  Complex(kind=8), intent(out) :: OvlMat(Nbas,N_lev)
!  Complex(kind=8), intent(out), dimension(Nbas,N_lev) :: OvlMat

  REAL(KIND=8), INTENT(out) :: atocc, atmag
  REAL(KIND=8), INTENT(out) :: Z
! Local
  TYPE(occup), DIMENSION(0:nlm*ns) :: N_occ
!  REAL(KIND=8) :: U(nlm*ns,nlm*ns,nlm*ns,nlm*ns)
!  COMPLEX(KIND=8) :: U(nlm*ns,nlm*ns,nlm*ns,nlm*ns)
  COMPLEX(KIND=8) :: zener
  COMPLEX(KIND=8) :: GF_lev(nlm*ns,nlm*ns,Iwmax)
  Complex(kind=8) :: Tail_lev(nmom,nlm*ns,nlm*ns)
  real(8) :: nomin, denom, E_B, tresh, maxexp, norm, ge
  real(8) :: fsign, Eground, Zterm
  real(8) :: atorb, atmag_lev
  real(8), allocatable :: E_A(:), docc(:), ener(:)
  !real(8), allocatable :: E_A(:), ummss(:,:), occ(:), docc(:), ener(:)
  integer, allocatable :: nground(:)
  INTEGER, PARAMETER :: numexp=650
  integer :: i, j, m, m1, is, is1, iom, ls
  integer :: k, l, ideg, ie, i1, k1, Nat, NN
  integer :: iloc, kl
  INTEGER :: nso, nstate
  REAL(KIND=8), EXTERNAL :: factor
  LOGICAL :: Efirst, Left, Right
  REAL(KIND=8), PARAMETER :: tol=1d-7
!
  COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: vec_tmp
  COMPLEX(KIND=8) :: Lad_mel
  REAL(KIND=8) :: phase
!

  IF (LadBS.AND.CalcOvl) THEN
      STOP 'LadBS and CalcOvl cannot be TRUE simulteneously'
  ENDIF

  WRITE(*,*)"!!!! gf_HI_fullU_INT is starting!"

  nso    = 2 * nlm
  nstate = 2**nso

! Diagonalize H_at with the 2-index U to estimate the ground state
! occupancy


  ALLOCATE( occ(nso), arr(0:nso-1) )
  FORALL( i = 0:nso-1 ) arr(i) = i

  allocate( docc(nso), ener(nso) )
  DO i=1,nso
     ener(i)     = REAL(e0f(i,i))
  ENDDO

!  ummss(1:nlm,1:nlm)         = ujmn
!  ummss(nlm+1:nso,nlm+1:nso) = ujmn
!  ummss(1:nlm,nlm+1:nso)     = umn
!  ummss(nlm+1:nso,1:nlm)     = umn

  allocate( E_A(0:nstate-1),nground(0:nstate-1) )
!
!    Initialize energy state E_A, A={ n_i sigma } and calculate Z
!
  do i = 0,nstate - 1
     occ = merge( 1.d0, 0.d0, btest( i, arr(0:nso-1) ) )
     E_A(i) = dot_product( ener, occ )
     docc = matmul( occ, ummss )
     E_A(i) = E_A(i) + dot_product( docc, occ ) / 2.d0
  enddo
  ge=MINVAL(E_A)
  nground=0
  atocc=0d0
  atmag=0d0
  DO i=0,nstate-1
     IF(ABS(ge-E_A(i)) < 1d-9) THEN
        nground(i)=1
        occ = merge( 1.d0, 0.d0, btest( i, arr(0:nso-1) ) )
        atocc=atocc+SUM(occ)
        atmag=atmag+SUM(occ(1:nlm))-SUM(occ(nlm+1:nso))
     ENDIF
  ENDDO
  norm=SUM(nground)
  atocc=atocc/norm
  atmag=atmag/norm
  if (verbosity>0) write(*,'(/,a,f12.5)')'Atomic occupancy with 2-ind U :',atocc
  if (verbosity>0) write(*,'(/,a,f13.7)')'Ground state energy with 2-ind U :',ge

!  Set up and diagonalize H_at matrices for N-1, N, N+1

!  CALL vertex4ind(Ur,U,nlm,ns)
!  DO i=1,nso
!     DO j=1,nso
!        DO k=1,nso
!           DO l=1,nso
!              IF (ABS(U(i,j,k,l)) > 1d-5) THEN
!                  write(*,*) i,j,k,l,REAL(U(i,j,k,l),KIND=8)
!              ENDIF
!           ENDDO
!        ENDDO
!     ENDDO
!  ENDDO

  Nat=NINT(atocc)
  iloc=0
  DO i=0,nso
     N_occ(i)%ifdiag=.FALSE.
     N_occ(i)%run=.FALSE.
     N_occ(i)%ndeg=0
  ENDDO
  IF(Nat==nso) THEN
     N_occ(nso-1:nso)%run=.TRUE.
  ELSEIF(Nat==0) THEN
     N_occ(0:1)%run=.TRUE.
  ELSE
     N_occ(Nat-1:Nat+1)%run=.TRUE.
  ENDIF

  DO WHILE(iloc==0)

     DO i=0,nso
        if (verbosity>0) write(*,'(a,I7)')'==> Starting N = ',i
        IF(.NOT.N_occ(i)%run.OR.N_occ(i)%ifdiag) CYCLE
        N_occ(i)%ifdiag=.TRUE.
        N_occ(i)%n=NINT(factor(nso)/factor(nso-i)/factor(i))
        NN=N_occ(i)%n
        ALLOCATE(N_occ(i)%Hn(NN,NN),N_occ(i)%En(NN))
        ALLOCATE(N_occ(i)%st_n(NN),N_occ(i)%narr(0:nstate-1))
        CALL diagH(N_occ(i)%Hn,N_occ(i)%En,e0f,U,N_occ(i)%st_n,arr,N_occ(i)%narr,nso,nstate,i,N_occ(i)%n,verbosity)
        if (verbosity>1) WRITE(*,'(a,I7,a,F14.7)')'The lowest energy for N= ',i,' is ',N_occ(i)%En(1)
        if (verbosity>1) write(*,'(a,I7,a,/)')'i = ',i,' done! <=='

     ENDDO

     Efirst=.TRUE.
     DO i=0,nso
        IF(N_occ(i)%ifdiag.AND.Efirst) THEN
           Eground=N_occ(i)%En(1)
           Nat=i
           Efirst=.FALSE.
        ELSEIF(N_occ(i)%ifdiag) then
           if(N_occ(i)%En(1)<Eground) THEN
              Eground=N_occ(i)%En(1)
              Nat=i
           endif
        ENDIF
     ENDDO

     IF((Nat.NE.0).and.(.NOT.N_occ(Nat-1)%ifdiag)) THEN
        N_occ(Nat-1)%run=.TRUE.
     ELSEIF((Nat.NE.nso) .and.(.NOT.N_occ(Nat+1)%ifdiag)) THEN
        N_occ(Nat+1)%run=.TRUE.
     ELSE
        iloc=1
     ENDIF

  ENDDO

  IF(remove_CF) THEN
     if (verbosity>0) write(*,'(/,a)')'CF splitting for GSM is removed !'
     N_occ(Nat)%En(1:N_lev)=Eground
  ENDIF

!  CALL calc_Z_occ(nso,nat,N_occ,arr,temp,zerotemp,Z,atocc,atmag,Eground)

  atocc=0d0
  atmag=0d0

  Eground=N_occ(Nat)%En(1)
  Z=0d0
  DO k=0,nso
     IF(.NOT.N_occ(k)%ifdiag) CYCLE
     DO i=1,N_occ(k)%n
        Zterm=EXP((Eground-N_occ(k)%En(i))/temp)
        IF(Zterm < tol) EXIT
        Z=Z+Zterm
        N_occ(k)%ndeg=i
        DO l=1,N_occ(k)%n
           occ = merge(1.d0,0.d0,btest(N_occ(k)%st_n(l),arr(0:nso-1)))
           atocc=atocc+SUM(occ)*N_occ(k)%Hn(l,i)*CONJG(N_occ(k)%Hn(l,i))*Zterm
           atmag=atmag+(SUM(occ(1:nlm))-SUM(occ(nlm+1:nso)))*N_occ(k)%Hn(l,i)*CONJG(N_occ(k)%Hn(l,i))*Zterm
        ENDDO
     ENDDO
  ENDDO

  atocc=atocc/Z
  atmag=atmag/Z


  if (verbosity>0) WRITE(*,'(/,a,i2,a,f13.6)') &
       &'The ground state has occupancy ',Nat, &
       &' and energy ',Eground
  if (verbosity>0) WRITE(*,'(a,i5,a)') &
       &'Transitions from  ',N_occ(Nat)%ndeg, &
       &' atomic states are included in GF'

  if (verbosity>0) WRITE(*,'(a,f13.6)')'Z = ',Z

  if (verbosity>0) write(*,'(/,a,f12.5)')'Atomic occupancy  :',atocc
  if (verbosity>0) write(*,'(/,a,f12.5)')'Atomic mag. mom.  :',atmag

  OPEN(450,file='ATOMIC_LEVELS')
  OPEN(320,file='STATES')
  WRITE(450,'(a)')'    #      E(eV)                OrbMM       SpinMM      TotM        TotMM'
  DO i=1,N_occ(Nat)%n
    atorb=0d0
    atmag_lev=0d0
    DO k=1,N_occ(Nat)%n
       WRITE(320,* ) N_occ(Nat)%Hn(k,i) ! real(N_occ(Nat)%Hn(k,i)), aimag(N_occ(Nat)%Hn(k,i))
    ENDDO
    DO k=1,N_occ(Nat)%n
       occ = merge(1.d0,0.d0,btest(N_occ(Nat)%st_n(k),arr(0:nso-1)))
       atmag_lev=atmag_lev+(SUM(occ(1:nlm))-SUM(occ(nlm+1:nso)))* &
 &         N_occ(Nat)%Hn(k,i)*CONJG(N_occ(Nat)%Hn(k,i))
       DO kl=1,nlm
        m=kl-(nlm+1)/2
        atorb=atorb+m*(occ(kl)+occ(nlm+kl))* &
 &         N_occ(Nat)%Hn(k,i)*CONJG(N_occ(Nat)%Hn(k,i))
      ENDDO
     ENDDO
     WRITE(450,'(i5,F16.9,5x,4f12.5)')i,N_occ(Nat)%En(i),atorb,atmag_lev, &
     &        atorb+atmag_lev/2d0, atorb+atmag_lev
  ENDDO
  CLOSE(450)
  CLOSE(320)

  ! Calculate the basis on N_lev first states from LadOP
  IF (LadBS) THEN
      ALLOCATE(st_gn(N_occ(Nat)%n),narr_g(0:nstate-1))
      narr_g=N_occ(Nat)%narr
      st_gn=N_occ(Nat)%st_n
      ALLOCATE(vec_tmp(N_occ(Nat)%n))
      if (verbosity>0) WRITE(*,'(/,a)')'Build N_lev states from operator LadOP'
      DO i=1,N_lev-1
         vec_tmp=OP(N_occ(Nat)%Hn(:,i),LadOP,nso,N_occ(Nat)%n)
         Lad_mel=SQRT(ovl(vec_tmp,vec_tmp))
         if (verbosity>0) WRITE(*,'(a,i2,a,i2,a,F9.5,a,F9.5)')'<',i+1,'|LadOP|',i,'> , VAL = ',ABS(Lad_mel),&
               ' PHASE = ',atan2(aimag(Lad_mel),real(Lad_mel))
         N_occ(Nat)%Hn(:,i+1)=vec_tmp/ABS(Lad_mel)
      ENDDO
      DO i=1,N_lev
         StBas(:,i)=N_occ(Nat)%Hn(:,i)
      ENDDO
      DEALLOCATE(st_gn,narr_g,vec_tmp)

  ELSE IF(CalcOvl) THEN
  ! Calculate overlap with the standard basis
      i=size(StBas,2)
      IF (i.NE.Nbas) STOP 'Wrong number of states in StMat'
      i=size(StBas,1)
      IF (i.NE.N_occ(Nat)%n) STOP 'Wrong vector length in Stmat'
      OvlMat=0d0
      DO i=1,Nbas
         DO k=1,N_lev
            OvlMat(i,k)=ovl(StBas(1:N_occ(Nat)%n,i),N_occ(Nat)%Hn(1:N_occ(Nat)%n,k))
         ENDDO
      ENDDO
   ENDIF

  ! END DEBUG

  ! Compute the Green's function

  if (verbosity>0) WRITE(*,'(/,a)')'Start GF calculations'

  GF0=(0d0,0d0)
  Tail0=(0d0,0d0)
  LEFT=.FALSE.
  RIGHT=.FALSE.
  k=Nat; l=Nat
  IF(Nat>0) THEN
     k=Nat-1
     LEFT=.TRUE.
  ENDIF
  IF(Nat<nso) THEN
     l=Nat+1
     RIGHT=.TRUE.
  ENDIF
  IF(LEFT.OR.RIGHT) then
     CALL add_to_GF_N(GF0,Tail0,arr,nso,nmom,Nat,Iwmax,l-k+1,N_occ(k:l),zmsb,Z,Eground,temp,LEFT,RIGHT)
  endif
  IF (verbosity>1) THEN
   WRITE(*,*)'Normalization of atomic GF:'
   DO m=1,nso
     WRITE(*,*)Tail0(1,m,m)
   ENDDO
  ENDIF

  GF0=(0d0,0d0)
  Tail0=(0d0,0d0)
  DO i=0,nso
     LEFT=.FALSE.
     RIGHT=.FALSE.
     k=i; l=i
     IF(i>0.AND.N_occ(i)%ndeg>0) THEN
        k=i-1
        LEFT=.TRUE.
     ENDIF
! Fix 15.11.2011
!    IF(i<nso.AND.N_occ(i)%ndeg>0.AND.N_occ(i+1)%ndeg==0) THEN
     IF(i<nso.AND.N_occ(i)%ndeg>0) THEN
        l=i+1
        RIGHT=.TRUE.
     ENDIF

     IF(LEFT.OR.RIGHT) then
        CALL add_to_GF_N(GF0,Tail0,arr,nso,nmom,Nat,Iwmax,l-k+1,N_occ(k:l),zmsb,Z,Eground,temp,LEFT,RIGHT)
     endif
  ENDDO
  IF (verbosity>1) THEN
      WRITE(*,*)'Normalization of atomic GF:'
      DO m=1,nso
        WRITE(*,*)Tail0(1,m,m)
      ENDDO
  ENDIF


! Compute the Green's function for all levels in the GSM

  if (verbosity>0) WRITE(*,'(/,a)')'Start GF calculations for each level of GSM'

  GF=(0d0,0d0)
  Tail=(0d0,0d0)
  LEFT=.FALSE.
  RIGHT=.FALSE.
  k=Nat; l=Nat
  IF(Nat>0) THEN
     k=Nat-1
     LEFT=.TRUE.
  ENDIF
  IF(Nat<nso) THEN
     l=Nat+1
     RIGHT=.TRUE.
  ENDIF
  N_occ(Nat)%ndeg=N_lev
  DO i=1,N_lev
     DO j=1,N_lev
        IF (calc_off_diag) THEN
           i1=j
        ELSE
           i1=i
        ENDIF
!        WRITE(*,*) GF_lev(0,0,0:4)
!        WRITE(*,*) '   '
        CALL add_to_GF_matrix(GF_lev,Tail_lev,arr,nso,nmom,i,i1,Iwmax,l-k+1,N_occ(k:l),zmsb,temp,LEFT,RIGHT)
        GF(i,i1,:,:,:)=GF_lev
        Tail(i,i1,:,:,:)=Tail_lev
        IF (verbosity>1.AND.i==i1) THEN
         WRITE(*,*)'Normalization of atomic GF: atomic level :',i
         DO m=1,nso
           WRITE(*,*)Tail(i,i,1,m,m)
         ENDDO
       ENDIF
       IF (.NOT.calc_off_diag) EXIT
     ENDDO
  ENDDO


  deallocate( occ, E_A, docc, ener, arr, nground )
  DO i=0,nso
     IF(N_occ(i)%ifdiag) DEALLOCATE(N_occ(i)%Hn,N_occ(i)%En,N_occ(i)%st_n,N_occ(i)%narr)
  ENDDO

  RETURN
END SUBROUTINE gf_HI_fullU_INT

REAL(KIND=8) function factor(N)

! factorial of N

  IMPLICIT NONE
  INTEGER, INTENT(in) :: N
  INTEGER :: i

  factor=1d0
  IF(N==0) RETURN
  DO i=1,N
     factor=factor*i
  ENDDO
  RETURN
END function factor

SUBROUTINE diagH(H,E,e0f,U,st,arr,narr,nso,nstate,Nat,n,verbosity)

! Initilize and diagonalize Hat for occupancy Nat

  IMPLICIT NONE
  INTEGER, INTENT(in) :: nso, nstate, Nat, n, verbosity
  INTEGER, INTENT(in) :: arr(nso)
  COMPLEX(KIND=8), INTENT(in) :: e0f(nso,nso)
  REAL(KIND=8), INTENT(in) :: U(nso,nso,nso,nso)
!  COMPLEX(KIND=8), INTENT(in) :: U(nso,nso,nso,nso)
  ! Output
  COMPLEX(KIND=8), INTENT(out) :: H(n,n)
  INTEGER, INTENT(out) :: st(n)
  REAL(KIND=8), INTENT(out) :: E(n)
  ! Locals
  INTEGER :: occ(nso), narr(0:nstate-1)
  INTEGER :: i, k, l, p, q, m, m1, m2, m3, m4, INFO
  INTEGER :: ks, ls, ps, qs
  REAL(KIND=8) :: fsign
  REAL(KIND=8), DIMENSION(3*n-2) :: rwork
  COMPLEX(KIND=8), DIMENSION(n*n) :: WORK

  k=0
  narr=-1
  DO i=0,nstate-1
     occ = merge( 1.d0, 0.d0, btest( i, arr(1:nso) ) )
     IF(SUM(occ)==Nat) THEN
        k=k+1
        st(k)=i
        narr(i)=k
     ENDIF
  ENDDO
  IF(k.NE.n) THEN
     WRITE(*,*)'Error in diagH: array size is different from number of states'
     STOP
  ENDIF
  H=(0d0,0d0)
  if (verbosity>0) WRITE(*,'(a,i7)')'Set up the Hamiltonian for N = ',Nat
  ! Set up one-particle term
  DO m=1,nso
     DO m1=1,nso
        IF(ABS(e0f(m,m1)) < 1d-9) CYCLE
        DO i=1,n
           occ = NINT(merge( 1.d0, 0.d0, btest( st(i), arr(1:nso) ) ))
           k=ibclr(st(i),m1-1)
           IF(k==st(i)) CYCLE
           occ(m1)=0
           ks=SUM(occ(1:m1-1))
           l=ibset(k,m-1)
           IF(l==k) CYCLE
           ls=SUM(occ(1:m-1))
           fsign=(-1d0)**(ks+ls)
           H(narr(l),i)=H(narr(l),i)+e0f(m,m1)*fsign
        ENDDO
     ENDDO
  ENDDO
  ! Add interaction term
  DO m=1,nso
     DO m1=1,nso
        DO m2=1,nso
           DO m3=1,nso
              IF(ABS(U(m,m1,m2,m3)) < 1d-9) CYCLE
              DO i=1,n
                 occ = NINT(merge( 1.d0, 0.d0, btest( st(i), arr(1:nso) ) ))
                 k=ibclr(st(i),m2-1)
                 IF(k==st(i)) CYCLE
                 occ(m2)=0
                 ks=SUM(occ(1:m2-1))
                 l=ibclr(k,m3-1)
                 IF(l==k) CYCLE
                 occ(m3)=0
                 ls=SUM(occ(1:m3-1))
                 p=ibset(l,m1-1)
                 IF(p==l) CYCLE
                 occ(m1)=1
                 ps=SUM(occ(1:m1-1))
                 q=ibset(p,m-1)
                 IF(q==p) CYCLE
                 qs=SUM(occ(1:m-1))
                 fsign=(-1d0)**(ks+ls+ps+qs)
                 H(narr(q),i)=H(narr(q),i)+0.5d0*U(m,m1,m2,m3)*fsign
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

! Diagonalize H

  if (verbosity>1) WRITE(*,'(a,i7)')'The Hamiltonian is set up'
  CALL ZHEEV('V','U',n,H,n,E,WORK,n*n,RWORK,INFO)
  IF(INFO.NE.0) THEN
     WRITE(*,*)'diagH : error in the matrix diagonalization'
     STOP
  ENDIF
  if (verbosity>1) WRITE(*,'(a)')'The Hamiltonian is diagonalized'

  RETURN
END SUBROUTINE diagH


SUBROUTINE add_to_GF_N(GF,Tail,arr,nso,nmom,Nat,Iwmax,num,N_occ,zmsb,Z,Eground,temp,LEFT,RIGHT)

! Adds to the GF the contribution associated to transitions from the states with
! occupancy N

  USE hubbard_I_data
  IMPLICIT NONE

  INTEGER, INTENT(in) :: nso, Nat, Iwmax, num,nmom
  LOGICAL, INTENT(in) :: LEFT, RIGHT
  REAL(KIND=8), INTENT(in) :: Z, Eground, temp
  COMPLEX(KIND=8), INTENT(in) :: zmsb(Iwmax)
  INTEGER, INTENT(in) :: arr(0:nso-1)
  TYPE(occup), INTENT(in) :: N_occ(100)

  COMPLEX(KIND=8), INTENT(inout) :: GF(nso,nso,Iwmax)
  complex(kind=8), intent(inout) :: Tail(nmom,nso,nso)

  COMPLEX(KIND=8) :: Gstore(nso,nso,Iwmax)
  complex(kind=8) :: Tailstore(nmom,nso,nso)
  COMPLEX(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: anmat,crmat
  COMPLEX(KIND=8) :: zener
  INTEGER :: i, k, k1, num1, m, m1, l, ls, ie
  INTEGER :: occ(nso)
  REAL(KIND=8) :: fsign, ecoff

  Gstore=GF
  GF=(0d0,0d0)
  Tailstore=Tail
  Tail=(0d0,0d0)

  IF(LEFT) THEN
     ! Matrix elements |<N-1|d_m|N>|
     ALLOCATE(anmat(nso,N_occ(1)%n,N_occ(2)%ndeg))
     anmat=0d0
     DO i=1,N_occ(2)%ndeg
        DO k=1,N_occ(2)%n
           DO m=1,nso
              occ = merge(1.d0,0.d0,btest(N_occ(2)%st_n(k),arr(0:nso-1)))
              l=ibclr(N_occ(2)%st_n(k),m-1)
              IF(l==N_occ(2)%st_n(k)) CYCLE
              ls=SUM(occ(1:m-1))
              fsign=(-1d0)**ls
              DO k1=1,N_occ(1)%n
                 anmat(m,k1,i)=anmat(m,k1,i)+CONJG(N_occ(2)%Hn(k,i))*N_occ(1)%Hn(N_occ(1)%narr(l),k1)*fsign
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     ! Compute contribution to GF
     DO i=1,N_occ(2)%ndeg              ! sum over contributing ground states
        DO k=1,N_occ(1)%n              ! sum over excited states
           IF(k >N_occ(1)%ndeg) THEN
              ecoff=EXP((Eground-N_occ(2)%En(i))/temp)
           ELSE
              ecoff=EXP((Eground-N_occ(1)%En(k))/temp)+EXP((Eground-N_occ(2)%En(i))/temp)
           ENDIF
           DO ie=1,Iwmax
              zener=1d0/(zmsb(ie)-N_occ(2)%En(i)+N_occ(1)%En(k))*ecoff
              DO m=1,nso
                 DO m1=1,nso
                    GF(m,m1,ie)=GF(m,m1,ie)+anmat(m1,k,i)*CONJG(anmat(m,k,i))*zener
                 ENDDO
              ENDDO
           ENDDO

           ! Tails:
           do ie=1,nmom   ! number of moments
              do m=1,nso
                 do m1=1,nso
                    Tail(ie,m,m1) = Tail(ie,m,m1) + ecoff *anmat(m1,k,i)*CONJG(anmat(m,k,i)) * &
                         &(N_occ(2)%En(i)-N_occ(1)%En(k))**(ie-1)
                 enddo
              enddo
           enddo

        ENDDO
     ENDDO
     DEALLOCATE(anmat)
  ENDIF

  IF(RIGHT) THEN
     ! Matrix elements |<N+1|d_m^+|N>|
     num1=num-1
     ALLOCATE(crmat(nso,N_occ(num)%n,N_occ(num1)%ndeg))
     crmat=0d0
     DO i=1,N_occ(num1)%ndeg
        DO k=1,N_occ(num1)%n
           DO m=1,nso
              occ = merge(1.d0,0.d0,btest(N_occ(num1)%st_n(k),arr(0:nso-1)))
              l=ibset(N_occ(num1)%st_n(k),m-1)
              IF(l==N_occ(num1)%st_n(k)) CYCLE
              ls=SUM(occ(1:m-1))
              fsign=(-1d0)**ls
              DO k1=1,N_occ(num)%n
                 crmat(m,k1,i)=crmat(m,k1,i)+CONJG(N_occ(num1)%Hn(k,i))*N_occ(num)%Hn(N_occ(num)%narr(l),k1)*fsign
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     ! Compute contribution to GF
     DO i=1,N_occ(num1)%ndeg
     !   DO k=1,N_occ(num)%n
     ! Fix 15.11.2011: only states above ndeg are included
        DO k=N_occ(num)%ndeg+1,N_occ(num)%n
           IF(k >N_occ(num)%ndeg) THEN
              ecoff=EXP((Eground-N_occ(num1)%En(i))/temp)
           ELSE
              ecoff=EXP((Eground-N_occ(num)%En(k))/temp)+EXP((Eground-N_occ(num1)%En(i))/temp)
           ENDIF
           DO ie=1,Iwmax
              zener=1d0/(zmsb(ie)-N_occ(num)%En(k)+N_occ(num1)%En(i))*ecoff
              DO m=1,nso
                 DO m1=1,nso
                    GF(m,m1,ie)=GF(m,m1,ie)+crmat(m,k,i)*CONJG(crmat(m1,k,i))*zener
                 ENDDO
              ENDDO
           ENDDO

           ! Tails:
           do ie=1,nmom   ! number of moments
              do m=1,nso
                 do m1=1,nso
                    Tail(ie,m,m1) = Tail(ie,m,m1) + ecoff * crmat(m,k,i)*CONJG(crmat(m1,k,i)) * &
                         &(N_occ(num)%En(k)-N_occ(num1)%En(i))**(ie-1)
                 enddo
              enddo
           enddo

        ENDDO
     ENDDO
     DEALLOCATE(crmat)
  ENDIF
  GF=Gstore+GF/Z
  Tail = Tailstore + Tail/Z
  RETURN
END SUBROUTINE add_to_GF_N



SUBROUTINE add_to_GF_matrix(GF,Tail,arr,nso,nmom,i_lev,j_lev,Iwmax,num,N_occ,zmsb,temp,LEFT,RIGHT)

! Evaluates GSM GF-matrix element [i_lev,j_lev]  using eq. 7 from PRB 94, 115117
!

  USE hubbard_I_data

  IMPLICIT NONE

  INTEGER, INTENT(in) :: nso, i_lev, j_lev,Iwmax, num, nmom
  LOGICAL, INTENT(in) :: LEFT, RIGHT
  COMPLEX(KIND=8), INTENT(in) :: zmsb(Iwmax)
  INTEGER, INTENT(in) :: arr(0:nso-1)
  REAL(KIND=8), INTENT(in) :: temp
  TYPE(occup), INTENT(in) :: N_occ(100) ! MODIFIED!

  COMPLEX(KIND=8), INTENT(out) :: GF(nso,nso,Iwmax)
  complex(kind=8), intent(out) :: Tail(nmom,nso,nso)

  COMPLEX(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: anmat,crmat
  COMPLEX(KIND=8) :: zener
  INTEGER :: i, j, ist, k, k1, num1, m, m1, l, ls, ie
  INTEGER :: occ(nso), istate(2)
  REAL(KIND=8) :: fsign, ecoff, de, fexp, Eav

  GF=(0d0,0d0)
  Tail=(0d0,0d0)
  istate(1)=i_lev
  istate(2)=j_lev
! Define E_av as average between those i_lev and j_lev to enforce propre behavior of G-matrix upon conjugation
  Eav=0.5*(N_occ(2)%En(istate(2))+N_occ(num-1)%En(istate(1)))
!
  IF(LEFT) THEN
     ! Matrix elements |<N-1|d_m|N>|
     ALLOCATE(anmat(nso,N_occ(1)%n,N_occ(2)%ndeg))
     anmat=0d0
     DO ist=1,2
       i=istate(ist)
       IF(ist == 2.AND.istate(1)==istate(2)) EXIT
       DO k=1,N_occ(2)%n
        DO m=1,nso
           occ = merge(1.d0,0.d0,btest(N_occ(2)%st_n(k),arr(0:nso-1)))
           l=ibclr(N_occ(2)%st_n(k),m-1)
           IF(l==N_occ(2)%st_n(k)) CYCLE
           ls=SUM(occ(1:m-1))
           fsign=(-1d0)**ls
           DO k1=1,N_occ(1)%n
              anmat(m,k1,i)=anmat(m,k1,i)+CONJG(N_occ(2)%Hn(k,i))*N_occ(1)%Hn(N_occ(1)%narr(l),k1)*fsign
           ENDDO
        ENDDO
       ENDDO
     ENDDO
     ! Compute contribution to GF
     DO k=1,N_occ(1)%n              ! sum over excited states
        ! substitute E of istate(2) by average
        !de=N_occ(1)%En(k)-N_occ(2)%En(istate(2))
        de=N_occ(1)%En(k)-Eav
        fexp=EXP(-de/temp)+1d0
        !fexp=1.0
        DO ie=1,Iwmax
           !zener=1d0/(zmsb(ie)-N_occ(2)%En(istate(2))+N_occ(1)%En(k))
           zener=fexp/(zmsb(ie)+de)
           DO m=1,nso
                 DO m1=1,nso
                    GF(m,m1,ie)=GF(m,m1,ie)+anmat(m1,k,istate(1))*CONJG(anmat(m,k,istate(2)))*zener
                 ENDDO
           ENDDO
        ENDDO

           ! Tails:
        do ie=1,nmom   ! number of moments
           do m=1,nso
              do m1=1,nso
                 Tail(ie,m,m1) = Tail(ie,m,m1) + anmat(m1,k,istate(1))*CONJG(anmat(m,k,istate(2))) * &
                      &(Eav-N_occ(1)%En(k))**(ie-1)*fexp
                      !&(N_occ(2)%En(istate(2))-N_occ(1)%En(k))**(ie-1)*fexp
              enddo
           enddo
        enddo

     ENDDO
     DEALLOCATE(anmat)
  ENDIF

  IF(RIGHT) THEN
     ! Matrix elements |<N+1|d_m^+|N>|
     num1=num-1
     ALLOCATE(crmat(nso,N_occ(num)%n,N_occ(num1)%ndeg))
     crmat=0d0
     DO ist =1,2
       i=istate(ist)
       IF(ist == 2.AND.istate(1)==istate(2)) EXIT
       DO k=1,N_occ(num1)%n
        DO m=1,nso
           occ = merge(1.d0,0.d0,btest(N_occ(num1)%st_n(k),arr(0:nso-1)))
           l=ibset(N_occ(num1)%st_n(k),m-1)
           IF(l==N_occ(num1)%st_n(k)) CYCLE
           ls=SUM(occ(1:m-1))
           fsign=(-1d0)**ls
           DO k1=1,N_occ(num)%n
              crmat(m,k1,i)=crmat(m,k1,i)+CONJG(N_occ(num1)%Hn(k,i))*N_occ(num)%Hn(N_occ(num)%narr(l),k1)*fsign
           ENDDO
        ENDDO
       ENDDO
     ENDDO
     ! Compute contribution to GF
     DO k=1,N_occ(num)%n
        ! substitute E of istate(2) by average
        !de=N_occ(num)%En(k)-N_occ(num1)%En(istate(1))
        de=N_occ(num)%En(k)-Eav
        !fexp=EXP(-de/temp)+1d0
        fexp=1.0
        DO ie=1,Iwmax
           !zener=1d0/(zmsb(ie)-N_occ(num)%En(k)+N_occ(num1)%En(istate(1)))
           zener=fexp/(zmsb(ie)-de)
           DO m=1,nso
              DO m1=1,nso
                 GF(m,m1,ie)=GF(m,m1,ie)+crmat(m,k,istate(1))*CONJG(crmat(m1,k,istate(2)))*zener
              ENDDO
           ENDDO
        ENDDO

           ! Tails:
        do ie=1,nmom   ! number of moments
           do m=1,nso
              do m1=1,nso
                 Tail(ie,m,m1) = Tail(ie,m,m1) +  crmat(m,k,istate(1))*CONJG(crmat(m1,k,istate(2))) * &
                      &(N_occ(num)%En(k)-Eav)**(ie-1)*fexp
                      !&(N_occ(num)%En(k)-N_occ(num1)%En(istate(1)))**(ie-1)*fexp
              enddo
           enddo
        enddo

     ENDDO
     DEALLOCATE(crmat)
  ENDIF
  RETURN
END SUBROUTINE add_to_GF_matrix


SUBROUTINE gf_HI_fullU(GF,Tail,e0f,U,ummss,zmsb,nlm,Iwmax,nmom,ns,atocc,atmag,temp,verbosity, &
                       remove_split,n_lev)

!
! Computes atomic GF with the full 4-index U 8.10.2007
! Finite temperature version
! (by L.V. Pourovskii)
!
!  Last change 23.10.2007
!
! Included tail calculation for the triqs package
! M. Aichhorn 10-2009
!
! GF /output/ - atomic Green's function
! e0f /input/ - atomic level position e0f_mm' C*_m C_m'
! U  /input/ - full 4-index U (orbitals and spins)
! ummss  /input/ - full 2-index U (orbitals and spins)
! zmsb(Iwmax)/input/ - COMPLEX energy mesh
! nlm, ns /input/ - orbital and spin degeneracy
! atocc, atmag /output/ - occupancy and magnetic moment of the atom
! temp /input/ - temperature
! verbosity/input/ - 0: no text output, 1: basics, 2: all
! remove_split/input/ - True: remove splitting between n_lev first levels of
!                       GS occupancy
! n_lev/input/ - the number of levels for which the splitting is removed, see
!                  above
  USE hubbard_I_data

  IMPLICIT NONE

! Input/output variables
  INTEGER, INTENT(in) :: nlm, ns, Iwmax, nmom
  COMPLEX(8), INTENT(in) :: e0f(nlm*ns,nlm*ns)
  COMPLEX(8), INTENT(in) :: zmsb(Iwmax)
!  REAL(8), INTENT(in) :: umn(nlm,nlm), ujmn(nlm,nlm)
!  REAL(8), INTENT(in) :: ur(nlm,nlm,nlm,nlm)
  REAL(8), INTENT(in) :: ummss(nlm*ns,nlm*ns)
  REAL(8), INTENT(in) :: U(nlm*ns,nlm*ns,nlm*ns,nlm*ns)
  REAL(8), INTENT(in) :: temp
  INTEGER, INTENT(in) :: verbosity, n_lev
  LOGICAL, INTENT(in) :: remove_split
  COMPLEX(8), INTENT(out) :: GF(nlm*ns,nlm*ns,Iwmax)
  COMPLEX(8), intent(out) :: Tail(nmom,nlm*ns,nlm*ns)
  REAL(8), INTENT(out) :: atocc, atmag
! Local
  TYPE(occup), DIMENSION(0:nlm*ns) :: N_occ
  integer, allocatable :: arr(:)
!  REAL(8) :: U(nlm*ns,nlm*ns,nlm*ns,nlm*ns)
  COMPLEX(8) :: zener
  real(8) :: Z, nomin, denom, E_B, tresh, maxexp, norm, ge, atorb
  real(8) :: fsign, Eground, Zterm
  real(8), allocatable :: E_A(:), occ(:), docc(:), ener(:)
  integer, allocatable :: nground(:)
  INTEGER, PARAMETER :: numexp=650
  integer :: i, j, m, m1, is, is1, iom, ls
  integer :: k, kl, l, ideg, ie, i1, k1, Nat, NN
  integer :: iloc
  INTEGER :: nso, nstate
  REAL(8), EXTERNAL :: factor
  REAL(8), PARAMETER :: tol=1d-7
  LOGICAL :: Efirst, Left, Right

!  print*, U

  WRITE(*,'()')

  nso    = 2 * nlm
  nstate = 2**nso

! Diagonalize H_at with the 2-index U to estimate the ground state
! occupancy


  ALLOCATE( occ(nso), arr(0:nso-1) )
  FORALL( i = 0:nso-1 ) arr(i) = i

  allocate( docc(nso), ener(nso) )
  DO i=1,nso
     ener(i)     = REAL(e0f(i,i))
  ENDDO

!  ummss(1:nlm,1:nlm)         = ujmn
!  ummss(nlm+1:nso,nlm+1:nso) = ujmn
!  ummss(1:nlm,nlm+1:nso)     = umn
!  ummss(nlm+1:nso,1:nlm)     = umn

  allocate( E_A(0:nstate-1),nground(0:nstate-1) )
!
!    Initialize energy state E_A, A={ n_i sigma } and calculate Z
!
  do i = 0,nstate - 1
     occ = merge( 1.d0, 0.d0, btest( i, arr(0:nso-1) ) )

     E_A(i) = dot_product( ener, occ )
     docc = matmul( occ, ummss )
     E_A(i) = E_A(i) + dot_product( docc, occ ) / 2.d0
  enddo
  ge=MINVAL(E_A)
  nground=0
  atocc=0d0
  atmag=0d0
  DO i=0,nstate-1
     IF(ABS(ge-E_A(i)) < 1d-9) THEN
        nground(i)=1
        occ = merge( 1.d0, 0.d0, btest( i, arr(0:nso-1) ) )
        atocc=atocc+SUM(occ)
        atmag=atmag+SUM(occ(1:nlm))-SUM(occ(nlm+1:nso))
     ENDIF
  ENDDO
  norm=SUM(nground)
  atocc=atocc/norm
  atmag=atmag/norm
  if (verbosity>0) write(*,'(/,a,f12.5)')'Atomic occupancy with 2-ind U :',atocc
  if (verbosity>0) write(*,'(/,a,f13.7)')'Ground state energy with 2-ind U :',ge

!  Set up and diagonalize H_at matrices for N-1, N, N+1

!  CALL vertex4ind(Ur,U,nlm,ns)

  Nat=NINT(atocc)
  iloc=0
  DO i=0,nso
     N_occ(i)%ifdiag=.FALSE.
     N_occ(i)%run=.FALSE.
     N_occ(i)%ndeg=0
  ENDDO
  IF(Nat==nso) THEN
     N_occ(nso-1:nso)%run=.TRUE.
  ELSEIF(Nat==0) THEN
     N_occ(0:1)%run=.TRUE.
  ELSE
     N_occ(Nat-1:Nat+1)%run=.TRUE.
  ENDIF

  DO WHILE(iloc==0)

     DO i=0,nso
        if (verbosity>0) write(*,'(a,I7)')'==> Starting N = ',i
        IF(.NOT.N_occ(i)%run.OR.N_occ(i)%ifdiag) CYCLE
        N_occ(i)%ifdiag=.TRUE.
        N_occ(i)%n=NINT(factor(nso)/factor(nso-i)/factor(i))
        NN=N_occ(i)%n
        ALLOCATE(N_occ(i)%Hn(NN,NN),N_occ(i)%En(NN))
        ALLOCATE(N_occ(i)%st_n(NN),N_occ(i)%narr(0:nstate-1))
        CALL diagH(N_occ(i)%Hn,N_occ(i)%En,e0f,U,N_occ(i)%st_n,arr,N_occ(i)%narr,nso,nstate,i,N_occ(i)%n,verbosity)
        if (verbosity>1) WRITE(*,'(a,I7,a,F14.7)')'The lowest energy for N= ',i,' is ',N_occ(i)%En(1)
        if (verbosity>1) write(*,'(a,I7,a,/)')'i = ',i,' done! <=='

     ENDDO

     Efirst=.TRUE.
     DO i=0,nso
        IF(N_occ(i)%ifdiag.AND.Efirst) THEN
           Eground=N_occ(i)%En(1)
           Nat=i
           Efirst=.FALSE.
        ELSEIF(N_occ(i)%ifdiag) then
           if(N_occ(i)%En(1)<Eground) THEN
              Eground=N_occ(i)%En(1)
              Nat=i
           endif
        ENDIF
     ENDDO

     IF((Nat.NE.0).and.(.NOT.N_occ(Nat-1)%ifdiag)) THEN
        N_occ(Nat-1)%run=.TRUE.
     ELSEIF((Nat.NE.nso) .and.(.NOT.N_occ(Nat+1)%ifdiag)) THEN
        N_occ(Nat+1)%run=.TRUE.
     ELSE
        iloc=1
     ENDIF

  ENDDO

  IF (remove_split) THEN
      ! remove splitting between first n_lev levels
      !
      OPEN(450,file='gs_multiplet.dat') !write these levels before removing splitting
      WRITE(450,'(a)') &
       '    #     E       M_orb       M_spin      J_tot       M_tot'
      DO i=1,n_lev
         atorb=0d0
         atmag=0d0
         DO k=1,N_occ(Nat)%n
            occ = merge(1.d0,0.d0,btest(N_occ(Nat)%st_n(k),arr(0:nso-1)))
            atmag=atmag+(SUM(occ(1:nlm))-SUM(occ(nlm+1:nso)))*N_occ(Nat)%Hn(k,i)*CONJG(N_occ(Nat)%Hn(k,i))
            DO kl=1,nlm
               m=kl-(nlm+1)/2
               atorb=atorb+m*(occ(kl)+occ(nlm+kl))*N_occ(Nat)%Hn(k,i)*CONJG(N_occ(Nat)%Hn(k,i))
            ENDDO
         ENDDO
         !WRITE(450,*)i,N_occ(Nat)%En(i)
         WRITE(450,'(i5,F15.8,4f9.5)')i,N_occ(Nat)%En(i),atorb,atmag,atorb+atmag/2d0,atorb+atmag
      ENDDO
      CLOSE(450)
      Eground=N_occ(Nat)%En(1)
      N_occ(Nat)%En(1:n_lev)=Eground
  ENDIF
  !END remove splitting

  atocc=0d0
  atmag=0d0
  Eground=N_occ(Nat)%En(1)
  Z=0d0
  DO k=0,nso
     IF(.NOT.N_occ(k)%ifdiag) CYCLE
     DO i=1,N_occ(k)%n
        Zterm=EXP((Eground-N_occ(k)%En(i))/temp)
        IF(Zterm < tol) EXIT
        Z=Z+Zterm
        N_occ(k)%ndeg=i
        DO l=1,N_occ(k)%n
           occ = merge(1.d0,0.d0,btest(N_occ(k)%st_n(l),arr(0:nso-1)))
           atocc=atocc+SUM(occ)*N_occ(k)%Hn(l,i)*CONJG(N_occ(k)%Hn(l,i))*Zterm
           atmag=atmag+(SUM(occ(1:nlm))-SUM(occ(nlm+1:nso)))*N_occ(k)%Hn(l,i)*CONJG(N_occ(k)%Hn(l,i))*Zterm
        ENDDO
     ENDDO
  ENDDO

  atocc=atocc/Z
  atmag=atmag/Z


  if (verbosity>0) WRITE(*,'(/,a,i2,a,f13.6)') &
       &'The ground state has occupancy ',Nat, &
       &' and energy ',Eground
  if (verbosity>0) WRITE(*,'(a,i5,a)') &
       &'Transitions from  ',N_occ(Nat)%ndeg, &
       &' atomic states are included in GF'

  if (verbosity>0) WRITE(*,'(a,f13.6)')'Z = ',Z

  if (verbosity>0) write(*,'(/,a,f12.5)')'Atomic occupancy  :',atocc
  if (verbosity>0) write(*,'(/,a,f12.5)')'Atomic mag. mom.  :',atmag

  OPEN(450,file='atomic_levels.dat')
  OPEN(320,file='states.dat')
  WRITE(450,'(a)')'    #     E       M_orb       M_spin      J_tot       M_tot'
  DO i=1,N_occ(Nat)%n
    atorb=0d0
    atmag=0d0
    WRITE(320,*)N_occ(Nat)%Hn(1:N_occ(Nat)%n,i)
    DO k=1,N_occ(Nat)%n
       occ = merge(1.d0,0.d0,btest(N_occ(Nat)%st_n(k),arr(0:nso-1)))
       atmag=atmag+(SUM(occ(1:nlm))-SUM(occ(nlm+1:nso)))*N_occ(Nat)%Hn(k,i)*CONJG(N_occ(Nat)%Hn(k,i))
       DO kl=1,nlm
          m=kl-(nlm+1)/2
          atorb=atorb+m*(occ(kl)+occ(nlm+kl))*N_occ(Nat)%Hn(k,i)*CONJG(N_occ(Nat)%Hn(k,i))
       ENDDO
    ENDDO
    !WRITE(450,*)i,N_occ(Nat)%En(i)
    WRITE(450,'(i5,F15.8,4f9.5)')i,N_occ(Nat)%En(i),atorb,atmag,atorb+atmag/2d0, atorb+atmag
  ENDDO
  CLOSE(450)
  CLOSE(320)

! Compute the Green's function

  if (verbosity>0) WRITE(*,'(/,a)')'Start GF calculations'
  ! Diagonalize H for additional occupancies if needed for GF
  DO k=0,nso
     IF(N_occ(k)%ndeg > 0) THEN
        IF(k>0.AND..NOT.N_occ(k-1)%ifdiag) N_occ(k-1)%run=.TRUE.
        IF(k<nso.AND..NOT.N_occ(k+1)%ifdiag) N_occ(k+1)%run=.TRUE.
     ENDIF
  ENDDO

  DO i=0,nso
     IF(.NOT.N_occ(i)%ifdiag.AND.N_occ(i)%run) THEN
        N_occ(i)%run=.TRUE.
        N_occ(i)%n=NINT(factor(nso)/factor(nso-i)/factor(i))
        NN=N_occ(i)%n
        ALLOCATE(N_occ(i)%Hn(NN,NN),N_occ(i)%En(NN))
        ALLOCATE(N_occ(i)%st_n(NN),N_occ(i)%narr(0:nstate-1))
        CALL diagH(N_occ(i)%Hn,N_occ(i)%En,e0f,U,N_occ(i)%st_n,arr,N_occ(i)%narr,nso,nstate,i,N_occ(i)%n,verbosity)
        if (verbosity>1) WRITE(*,'(/,a,I7,a,F14.7)')'The lowest energy for N= ',i,' is ',N_occ(i)%En(1)
     ENDIF
  ENDDO

  GF=(0d0,0d0)
  Tail=(0d0,0d0)
  DO i=0,nso
     LEFT=.FALSE.
     RIGHT=.FALSE.
     k=i; l=i
     IF(i>0.AND.N_occ(i)%ndeg>0) THEN
        k=i-1
        LEFT=.TRUE.
     ENDIF
! Fix 15.11.2011
!    IF(i<nso.AND.N_occ(i)%ndeg>0.AND.N_occ(i+1)%ndeg==0) THEN
     IF(i<nso.AND.N_occ(i)%ndeg>0) THEN
        l=i+1
        RIGHT=.TRUE.
     ENDIF

     IF(LEFT.OR.RIGHT) then
        CALL add_to_GF_N(GF,Tail,arr,nso,nmom,Nat,Iwmax,l-k+1,N_occ(k:l),zmsb,Z,Eground,temp,LEFT,RIGHT)
     endif
  ENDDO
  IF (verbosity>1) THEN
      WRITE(*,*)'Normalization of atomic GF:'
      DO m=1,nso
        WRITE(*,*)Tail(1,m,m)
      ENDDO
  ENDIF

  deallocate( occ, E_A, docc, ener, arr, nground )
  DO i=0,nso
     IF(N_occ(i)%ifdiag) DEALLOCATE(N_occ(i)%Hn,N_occ(i)%En,N_occ(i)%st_n,N_occ(i)%narr)
  ENDDO

  RETURN
END SUBROUTINE gf_HI_fullU


