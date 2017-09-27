import numpy as np


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

def force(rsq):
  #
  #
  #
  #  Implement the force equation here. Note that the force is minus the derivative of the potential. What is passed into the function is the square of the interparticle distance. What is used on line 79 is dx times the result of this function
  #
  #
  #

def potential(rsq):
  rsqinv=1./rsq
  r6inv=rsqinv*rsqinv*rsqinv
  return -4*r6inv*(1-r6inv)

def compute_kinetic_energy(vx,vy):
  return 0.5*sum(vx*vx+vy*vy)

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
 
def euler(rx,ry,vx,vy,dV_drx,dV_dry):
  deltat=0.001
  #update the positions
  rx+=deltat*vx
  ry+=deltat*vy
  
  #update the velocities
  vx-=deltat*dV_drx
  vy-=deltat*dV_dry

 
  #put back into box:
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

def main():

  #simulation parameters
  Nx=5; Ny=5; N=Nx*Ny #set particles onto a grid initially
  L=5
  Nstep=10000
  rcut=2.5 # a usual choice for the cutoff radius


  vx=np.zeros(N)
  vy=np.zeros(N)
  rx=np.zeros(N)
  ry=np.zeros(N)

  rxlog=np.zeros([Nstep,N])
  rylog=np.zeros([Nstep,N])
  vxlog=np.zeros([Nstep,N])
  vylog=np.zeros([Nstep,N])

  initialize_positions_and_velocities(rx,ry,vx,vy,Nx,Ny,L)
  for i in range(Nstep):
    dV_drx=np.zeros(N)
    dV_dry=np.zeros(N)
    compute_forces(rx,ry,dV_drx,dV_dry, N, L, rcut)

    #propagate using forward Euler
    euler(rx,ry,vx,vy,dV_drx,dV_dry)
    #
    #
    #
    # Replace the FW Euler with a velocity verlet
    #
    #
    #

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
    print i,Epot,Ekin,Epot+Ekin
  #print result
  print_result(rxlog,rylog,vxlog,vylog)


if __name__ == "__main__":
    main()
