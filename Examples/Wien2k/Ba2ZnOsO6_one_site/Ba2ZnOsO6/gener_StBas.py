from h5 import *
from MagInt.HubbardI_interact import *
from hubbardI.hubbard_I import mat_el_vecs, get_vec
#from MagInt.HubbardI_interact import get_vec, mat_el_vecs, mat_el_vecs_2el, vertex4ind, ops_on_state, prn_vectors
from MagInt.utils import *
import numpy
import math
from MagInt.Multipolar import *
import sys


beta=40
scale=1.0
U_int= 3.20
J_hund =0.50
#
l=2
N=2
spin=1.0
nlm=2*l+1
nlms=2*nlm
vlen=int(factorial(nlms)/factorial(N)/factorial(nlms-N))
read_eal=True
#read_eal=True
useSO=True
split=0.02
#split=0.00001

tol_coff=0.000001

#effective J basis
J=2.0
L=1
N_lev=int(2*J+1)
N_bas=int((2*L+1)*(2*spin+1))
N_mult=int(abs(L+spin)-abs(L-spin))+1
J_mult=[J-i for i in range(N_mult)] 
mult_deg=[int(J_mult[m]*2+1) for m in range(N_mult)]
Nlm=2*L+1
Nlms=Nlm*2

print(N_lev,N_bas,N_mult,J_mult,mult_deg)

#read_eal=True
#=

# CF fit, used if read_eal=False
CF_Lcoff=[[4,0],[4,4]]

#CF_Ham_file='CF_Hamiltonian_Non_SP_CFP.dat'
#spin_pol_L=False
#eal_fname='eal_site_0.dat'
eal_fname='eal_last.dat'
spin_pol_L=False



def gener_J_mat(j):
    Ndim=int(2*j)+1
    Jz=numpy.zeros((Ndim,Ndim),complex)
    Jp=numpy.zeros((Ndim,Ndim),complex)
    Jm=numpy.zeros((Ndim,Ndim),complex)
    for i in range(Ndim):
        #m=j-i
        m=i-j
        Jz[i,i]=m
        if i>0:
           Jp[i,i-1]=math.sqrt((j-m+1.0)*(j+m))
           Jm[i-1,i]=math.sqrt((j-m+1.0)*(j+m))
    return Jz,Jp,Jm

eal={}
# read eal
nso=nlms
if read_eal:
    eal_tmp=numpy.loadtxt(eal_fname)
    eal['ud']=numpy.zeros((nso,nso),complex)
    for i in range(nso):
        for j in range(nso):
            eal['ud'][i,j]=eal_tmp[i,2*j]+1j*eal_tmp[i,2*j+1]
            #eal['ud'][i,j]=eal_tmp[i,j]



Ndim=nlms
tt=numpy.loadtxt('real_d_harms')
transmat=numpy.zeros((Ndim,Ndim),complex)
for i in range(Ndim):
    for j in range(Ndim):
        transmat[i,j]=tt[i,2*j]+1j*tt[i,2*j+1]

tmat_cub=transmat[0:3,0:5]

eal_cub=numpy.dot(transmat,numpy.dot(eal['ud'],transmat.conjugate().transpose()))

# DEBUG : set eg-t2g coupling to zero
print_arr(eal_cub,log='\neal_cub before removing t2g-eg coupling')
for i in range(Ndim):
    for j in range(Ndim):
        if (i>5 or j>5) and i!=j: eal_cub[i,j]=0.0
eal['ud']=numpy.dot(transmat.conjugate().transpose(),numpy.dot(eal_cub,transmat))
print_arr(eal_cub,log='\neal_cub')


Sz=numpy.zeros((Ndim,Ndim),complex)
Sp=numpy.zeros((Ndim,Ndim),complex)
Sm=numpy.zeros((Ndim,Ndim),complex)


# generate L,S,J matrices
for i in range(nlm):
    Sp[i,i+nlm]=1.0
    Sm[i+nlm,i]=1.0
    Sz[i,i]=0.5
    Sz[i+nlm,i+nlm]=-0.5

Lz=numpy.zeros((Ndim,Ndim),complex)
Lp=numpy.zeros((Ndim,Ndim),complex)
Lm=numpy.zeros((Ndim,Ndim),complex)

psL_bas=numpy.zeros((3,5),complex)
sq12=1.0/numpy.sqrt(2.0)
psL_bas[0,3]=-1.0
psL_bas[1,0]=-sq12
psL_bas[1,4]=sq12
psL_bas[2,1]=1.0

Tz,Tp,Tm=gener_J_mat(1)


Tz=numpy.dot(psL_bas.conjugate().transpose(),numpy.dot(Tz,psL_bas))
Tp=numpy.dot(psL_bas.conjugate().transpose(),numpy.dot(Tp,psL_bas))
Tm=numpy.dot(psL_bas.conjugate().transpose(),numpy.dot(Tm,psL_bas))


lz=numpy.kron(numpy.identity(2),Tz)
lp=numpy.kron(numpy.identity(2),Tp)
lm=numpy.kron(numpy.identity(2),Tm)

Jz=Sz+lz
Jp=Sp+lp
Jm=Sm+lm
Jx=0.5*(Jp+Jm)
Jy=0.5j*(Jm-Jp)

Sx=0.5*(Sp+Sm)
Sy=0.5j*(Sm-Sp)

#Lx=0.5*(Lp+Lm)
#Ly=0.5j*(Lm-Lp)

# Physical moment
llz,llp,llm=gener_J_mat(2)
Lz=numpy.kron(numpy.identity(2),llz)
Lp=numpy.kron(numpy.identity(2),llp)
Lm=numpy.kron(numpy.identity(2),llm)
Lx=0.5*(Lp+Lm)
Ly=0.5j*(Lm-Lp)


M={}
M['x']=Lx+2.0*Sx
M['y']=Ly+2.0*Sy
M['z']=Lz+2.0*Sz


print_arr(M['z'],log='\nMz')

eal['ud']+=split*M['z']


S = HubbardI_interact(beta = beta, l = 2, n_lev=N_lev,U_int=U_int*scale,J_hund=J_hund*scale,verbosity=2, use_spin_orbit=True)
S.Nmoments=4
#
#
S.set_ud_levels( eal = eal )
S.run_HI()



states=[]

vecs_bas = np.zeros((N_bas, vlen), dtype=np.cdouble, order='F')
vec = np.zeros(vlen, dtype=np.cdouble, order='F')
for i in range(N_lev):
    vecs_bas[i, :] = my_get_vec(n=vlen, iv=i, vec=vec)
    ovl=mat_el_vecs(nat=N,nso=nlms,n=vlen,st_l=vecs_bas[i,:],st_r=vecs_bas[i,:],op_mat=M['z'])
    print('%s <Mz> = %12.6f '%(i,ovl.real))
for i in range(N_lev):
    ovl=mat_el_vecs(nat=N,nso=nlms,n=vlen,st_l=vecs_bas[i,:],st_r=vecs_bas[i,:],op_mat=Lz)
    print('%s  <Lz> = %12.6f '%(i,ovl.real))
for i in range(N_lev):
    ovl=mat_el_vecs(nat=N,nso=nlms,n=vlen,st_l=vecs_bas[i,:],st_r=vecs_bas[i,:],op_mat=Sz)
    print('%s <Sz> = %12.6f '%(i,ovl.real))
print('\n\n Jz mel:')
for i in range(N_lev):
    ovl=mat_el_vecs(nat=N,nso=nlms,n=vlen,st_l=vecs_bas[i,:],st_r=vecs_bas[i,:],op_mat=Jz)
    print('%s <Jz> = %12.6f %12.6f'%(i,ovl.real,ovl.imag))
    states.append('|%s>'%(int(round(ovl.real))))
    print(states[-1])
print('\n\n pseudo-Lz mel:')
for i in range(N_lev):
    ovl=mat_el_vecs(nat=N,nso=nlms,n=vlen,st_l=vecs_bas[i,:],st_r=vecs_bas[i,:],op_mat=lz)
    print('%s <lz> = %12.6f %12.6f'%(i,ovl.real,ovl.imag))
print('\n\n J+ mel:')
for i in range(N_lev-1):
    ovl=mat_el_vecs(nat=N,nso=nlms,n=vlen,st_l=vecs_bas[i+1,:],st_r=vecs_bas[i,:],op_mat=Jp)
    vecs_bas[i+1,:]*=numpy.sign(ovl.real)
    ovl=mat_el_vecs(nat=N,nso=nlms,n=vlen,st_l=vecs_bas[i+1,:],st_r=vecs_bas[i,:],op_mat=Jp)
    print('%s <J+> = %12.6f %12.6f'%(i,ovl.real,ovl.imag))
# read symm eal
eal_tmp=numpy.loadtxt(eal_fname)
for i in range(nlms):
    for j in range(nlms):
        eal['ud'][i,j]=eal_tmp[i,2*j]+1j*eal_tmp[i,2*j+1]


S = HubbardI_interact(beta = beta, l = 2, n_lev=N_lev,U_int=U_int,J_hund=J_hund,verbosity=2, use_spin_orbit=True)
S.Nmoments=4
#
#
S.set_ud_levels( eal = eal )
S.run_HI()

at_lev=numpy.loadtxt('ATOMIC_LEVELS')

E=(at_lev[:,1]-at_lev[0,1])*1000

vecs_lev = np.zeros((N_lev, vlen), dtype=np.cdouble, order='F')
vec = np.zeros(vlen, dtype=np.cdouble, order='F')
for i in range(N_lev):
    vecs_lev[i, :] = my_get_vec(n=vlen, iv=i, vec=vec)

I=numpy.identity(nlms,complex)

#tmat to S2 pseudospin
tmat_S2=numpy.zeros((N_lev,N_lev))

print('\n\n CF ENERGIES and STATES:\n')
for i in range(N_lev):
    exp_coff=[]
    for j in range(N_bas):
        coff=mat_el_vecs(nat=N,nso=nlms,n=vlen,st_l=vecs_bas[j,:],st_r=vecs_lev[i,:],op_mat=I)
        if abs(coff*coff/N/N) > 1e-5: exp_coff.append([j,coff/N])
    st=''
    norm=0.0
    cf_val=[exp_coff[k][1] for k in range(len(exp_coff))]
    iarr=[k[0] for k in sorted(enumerate(cf_val), key=lambda x:abs(x[1]),reverse=True)]
    fsign=numpy.sign(exp_coff[iarr[0]][1])
    vecs_lev[i,:]*=fsign
    for icf in iarr:
        exp_coff[icf][1]*=fsign
        cf=exp_coff[icf]
        if abs(cf[1]) > tol_coff:
            if  abs(cf[1].imag) < tol_coff:
                st+="%+9.6f%s"%(cf[1].real,states[cf[0]])
            else:
                st+="[%+9.6f%+9.6f*I]%s"%(cf[1].real,cf[1].imag,states[cf[0]])
            norm+=cf[1].real*cf[1].real
            # coefficients in tmat_S2 are rounded either to 1 or to 1/sqrt(2)
            tmat_S2[i,cf[0]]=round(cf[1].real) if abs(abs(cf[1].real)-1.0) < 5e-2 else numpy.sign(cf[1].real)*sq12 
    #print '\n%s E= %6.0f  $%s$   Norm=%7.4f'%(i+1,E[i],st,math.sqrt(norm))
    print(' %6.1f & $%s$  \\\\ \n'%(E[i],st))


print_arr(tmat_S2,log='\ntmat_S2')
# trnsformation CF-> S2
vecs_stbas=numpy.dot(tmat_S2.transpose(),vecs_lev)

H_S2=numpy.zeros((N_lev,N_lev))
for i in range(N_lev):
    H_S2+=E[i]*numpy.outer(tmat_S2[i,:],tmat_S2[i,:])

print_arr(H_S2,log='\n\n# remnant CF in pseudoJ=2 basis:',decd=6,prn_zero_imag=True)

# check that it has right Jz values
print('\n\n J mel of final pseudoJ=2 basis:')
Jmat={}
for dir in ['x','y','z']: Jmat[dir]=numpy.zeros((N_lev,N_lev),complex)
for i in range(N_lev):
    for j in range(N_lev):
        Jmat['x'][i,j]=mat_el_vecs(nat=N,nso=nlms,n=vlen,st_l=vecs_stbas[i,:],st_r=vecs_stbas[j,:],op_mat=Jx)
        Jmat['y'][i,j]=mat_el_vecs(nat=N,nso=nlms,n=vlen,st_l=vecs_stbas[i,:],st_r=vecs_stbas[j,:],op_mat=Jy)
        Jmat['z'][i,j]=mat_el_vecs(nat=N,nso=nlms,n=vlen,st_l=vecs_stbas[i,:],st_r=vecs_stbas[j,:],op_mat=Jz)

print('\n\n J+ mel of final pseudoJ=2 basis:')
for i in range(N_lev-1):
    print(mat_el_vecs(nat=N,nso=nlms,n=vlen,st_l=vecs_stbas[i+1,:],st_r=vecs_stbas[i,:],op_mat=Jp))

for dir in ['x','y','z']:
    print_arr(Jmat[dir],log='\n\nJmat[%s]'%dir,decd=5)

# Calculate Mag. mom. matrices in this basis
Mmat={}
for dir in ['x','y','z']:
    Mmat[dir]=numpy.zeros((N_lev,N_lev),complex)
    for i in range(N_lev):
        for j in range(N_lev):
            Mmat[dir][i,j]=mat_el_vecs(nat=N,nso=nlms,n=vlen,st_l=vecs_stbas[i,:],st_r=vecs_stbas[j,:],op_mat=M[dir])
    print_arr(Mmat[dir],log='# M%s'%dir,decd=6,prn_zero_imag=True)

print("\nPhysical spin and orbital moment:")
Lops={'x':Lx,'y':Ly,'z':Lz}
Sops={'x':Sx,'y':Sy,'z':Sz}
for dir in ['x','y','z']:
    Mmat[dir]=numpy.zeros((N_lev,N_lev),complex)
    for i in range(N_lev):
        for j in range(N_lev):
            Mmat[dir][i,j]=mat_el_vecs(nat=N,nso=nlms,n=vlen,st_l=vecs_stbas[i,:],st_r=vecs_stbas[j,:],op_mat=Sops[dir])
    print_arr(2.0*Mmat[dir],log='# Mspin %s'%dir,decd=6,prn_zero_imag=True)

    Mmat[dir]=numpy.zeros((N_lev,N_lev),complex)
    for i in range(N_lev):
        for j in range(N_lev):
            Mmat[dir][i,j]=mat_el_vecs(nat=N,nso=nlms,n=vlen,st_l=vecs_stbas[i,:],st_r=vecs_stbas[j,:],op_mat=Lops[dir])
    print_arr(Mmat[dir],log='# Morb %s'%dir,decd=6,prn_zero_imag=True)

ar = HDFArchive('Standard_Basis.h5', 'a')
ar['Ion0'] = numpy.transpose(vecs_stbas[0:N_lev, :])
# ar['MMats' + str(i_sh)] = Mmat
# ar['H_S2' + str(i_sh)] = H_S2
del ar

