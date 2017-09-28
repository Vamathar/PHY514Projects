import numpy as np
from numba import jit

@jit
def initialize_positions_and_velocities(rx,ry,vx,vy, Nx, Ny,L):
  dx=L/Nx;
  dy=L/Ny;
  np.random.seed(0)
  for i in range(Nx):
    for j in range(Ny):
      rx[i*Ny+j]=dx*(i+0.5)
      ry[i*Ny+j]=dy*(j+0.5)
      
      u=np.random.random() #This is box muller
      v=np.random.random()
      vx[i*Ny+j]=np.sqrt(-2*np.log(u))*np.cos(2.*np.pi*v)
      vy[i*Ny+j]=np.sqrt(-2*np.log(u))*np.sin(2.*np.pi*v)
  #subtract net velocity to avoid global drift
  vxav=sum(vx)/vx.size
  vyav=sum(vy)/vx.size
  vx-=vxav
  vy-=vyav

@jit
def force(rsq):
    rsqinv=1./rsq
    return 24.*(rsqinv**4-2.*rsqinv**7)

@jit
def potential(rsq):
  rsqinv=1./rsq
  r6inv=rsqinv*rsqinv*rsqinv
  return -4*r6inv*(1-r6inv)

@jit
def compute_kinetic_energy(vx,vy):
  return 0.5*sum(vx*vx+vy*vy)

@jit
def compute_potential_energy(rx,ry,rcut,L):
  rcutsq=rcut*rcut
  rcutv=potential(rcutsq) #shift the potential to avoid jump at rc
  Epot=0. 
  for i in range(rx.size):
    for j in range(i):
      dx=rx[i]-rx[j]
      dy=ry[i]-ry[j]
      #minimum image convention
      if(dx > L/2.): dx=dx-L
      if(dx <-L/2.): dx=dx+L
      if(dy > L/2.): dy=dy-L
      if(dy <-L/2.): dy=dy+L
      #print dx,dy
      #compute the distance
      rsq=dx*dx+dy*dy
      if(rsq < rcutsq):
        Epot+=potential(rsq)-rcutv
  return Epot 
      
@jit
def compute_forces(rx,ry,dV_drx, dV_dry, N, L, rcut):
  rcutsq=rcut*rcut
  for i in range(N):
    for j in range(i):
      dx=rx[i]-rx[j] ; 
      dy=ry[i]-ry[j] ; 
      #minimum image convention
      if(dx > L/2.): dx=dx-L
      if(dx <-L/2.): dx=dx+L
      if(dy > L/2.): dy=dy-L
      if(dy <-L/2.): dy=dy+L
      #compute the distance
      rsq=dx*dx+dy*dy
      #check if we are < the cutoff radius
      if(rsq < rcutsq):
        #here is the call of the force calculation
        dV_dr=force(rsq)
        
        #here the force is being added to the particle. Note the additional dx
        dV_drx[i]+=dx*dV_dr
        dV_drx[j]-=dx*dV_dr
        dV_dry[i]+=dy*dV_dr
        dV_dry[j]-=dy*dV_dr
 
@jit
def euler(rx,ry,vx,vy,dV_drx,dV_dry):
  deltat=0.001
  #update the positions
  rx+=deltat*vx
  ry+=deltat*vy
  
  #update the velocities
  vx-=deltat*dV_drx
  vy-=deltat*dV_dry

@jit
def verlet(rx,ry,vx,vy,dV_drx,dV_dry,N,L,rcut):
  deltat=0.001
  rx+=(deltat*vx-0.5*deltat*deltat*dV_drx)
  ry+=(deltat*vy-0.5*deltat*deltat*dV_dry)
  rebox(rx,ry,L)
  dV_drx_old = dV_drx.copy()
  dV_dry_old = dV_dry.copy()
  dV_drx*=0.
  dV_dry*=0.
  compute_forces(rx,ry,dV_drx,dV_dry,N,L,rcut)
  vx-=0.5*deltat*(dV_drx_old+dV_drx)
  vy-=0.5*deltat*(dV_dry_old+dV_dry)
 
  #put back into box:
@jit
def rebox(rx,ry,L):
  for i in range(rx.size):
    if rx[i] > L:
      rx[i]=rx[i]-L
    if rx[i] < 0:
      rx[i]=rx[i]+L
    if ry[i] > L:
      ry[i]=ry[i]-L
    if ry[i] < 0:
      ry[i]=ry[i]+L

def print_result(rxlog,rylog,vxlog,vylog):
  fr=open("positions.dat",'w')
  fv=open("velocities.dat",'w')

  for j in range(rxlog.shape[1]):
    for i in range(rxlog.shape[0]):
      fr.write(str(rxlog[i,j])+" "+str(rylog[i,j])+'\n')
      fv.write(str(vxlog[i,j])+" "+str(vylog[i,j])+'\n')
    fr.write('\n')
    fv.write('\n')

@jit
def diracdelta(r,dr,rtest):
  if (r <= rtest) and (rtest < r+dr):
    return 1.
  else:
    return 0.

@jit
def avgcorr(r,dr,rx,ry,L):
  N = len(rx)
  total=0.
  denominator = 0
  for i in range(N):
    for j in range(N):
      if i != j:
        denominator+=1
        dx=rx[i]-rx[j]
        dy=ry[i]-ry[j]
        rtest=np.sqrt(dx*dx+dy*dy)
        total+=diracdelta(r,dr,rtest)
  return total/denominator

@jit
def correlation(r,dr,rx,ry,L):
  density=len(rx)/(L*L)
  coeff = 1./(density*2.*np.pi*r*dr)
  return coeff*avgcorr(r,dr,rx,ry,L)

@jit
def run(Nx,Ny,L,rcut,Nstep):
  N = Nx*Ny
  vx=np.zeros(N)
  vy=np.zeros(N)
  rx=np.zeros(N)
  ry=np.zeros(N)

  rxlog=np.zeros([Nstep,N])
  rylog=np.zeros([Nstep,N])
  vxlog=np.zeros([Nstep,N])
  vylog=np.zeros([Nstep,N])
  elog=np.zeros([Nstep,3])

  initialize_positions_and_velocities(rx,ry,vx,vy,Nx,Ny,L)
  for i in range(Nstep):
    dV_drx=np.zeros(N)
    dV_dry=np.zeros(N)
    compute_forces(rx,ry,dV_drx,dV_dry, N, L, rcut)

    #propagate using forward Euler
    #euler(rx,ry,vx,vy,dV_drx,dV_dry)
    verlet(rx,ry,vx,vy,dV_drx,dV_dry,N,L,rcut)
    #make sure we're still in the box
    rebox(rx,ry,L)

    #keep track for printing
    rxlog[i]=rx
    rylog[i]=ry
    vxlog[i]=vx
    vylog[i]=vy
 
    #get some observables
    Epot=compute_potential_energy(rx,ry,rcut,L)
    Ekin=compute_kinetic_energy(vx,vy)
    elog[i]=[Ekin,Epot,Ekin+Epot]
  #print result
  #print_result(rxlog,rylog,vxlog,vylog)
  nr = 1000
  dr = L/nr
  rlist = np.linspace(dr,np.sqrt(2)*L,nr)
  glist = np.zeros(len(rlist))
  i = 0
  for r in rlist:
    for time in range(Nstep):
      glist[i]+=correlation(r,dr,rxlog[time],rylog[time],L)
    glist[i]=glist[i]/Nstep
    i+=1
  corrlist=np.array([rlist,glist])
  np.savetxt("xpositions.dat",rxlog,delimiter=",")
  np.savetxt("ypositions.dat",rylog,delimiter=",")
  np.savetxt("xvelocities.dat",vxlog,delimiter=",")
  np.savetxt("yvelocities.dat",vylog,delimiter=",")
  np.savetxt("energies.dat",elog,delimiter=",")
  np.savetxt("correlation.dat",corrlist,delimiter=",")

run(2,2,5,2.5,100)

