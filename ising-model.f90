program ising_model
implicit none
real*8::T,J_ising,E,M,r,mag,Ei,Ef,dE,u,h
real*8::avg_M,avg_E,avg_M_N,avg_E_N,avg_M_2,avg_E_2,chi,cv,avg_M_4,avg_E_4,bc
integer::L,niter,seed,N,i,j,k,a1,a2,a3,a4,a5,a6,time,mm,nn,ll,Temp,N_eqm,N_stat
integer,dimension(:,:,:),allocatable::spin

print*, 'Enter the size of the lattice'
read*,L

print*, 'Desired number of iterations'
read*,niter

!print*, 'Desired Temp'
!read*,T 

J_ising=1.0d0

seed=27863

allocate(spin(L,L,L))
E=0.0d0
M=0.0d0
N=L*L*L

N_eqm=10000
N_stat=10

call random_seed

!initializing lattice
open(10,file='initial_lattice.dat')
do i=1,L
 do j=1,L
  do k=1,L
   call random_number(r)
!  spin(i,j,k)=1
   if(r<0.5)then
    spin(i,j,k)=-1
   else
    spin(i,j,k)=1
   end if
   write(10,*) float(spin(i,j,k))
  end do
 end do 
end do
close(10)

!caculations
do i=1,L
 do j=1,L
  do k=1,l
   a1=i+1;a2=i-1;a3=j+1;a4=j-1;a5=k+1;a6=k-1 !neighbours
   !setting boundary conditions
   if(i==L)a1=1
   if(i==1)a2=L
   if(j==L)a3=1
   if(j==1)a4=L
   if(k==L)a5=1
   if(k==1)a6=L
   
   M=M+spin(i,j,k)
   E=E-J_ising*float(spin(i,j,k)*(spin(a1,j,k)+spin(a2,j,k)+spin(i,a3,k)+spin(i,a4,k)+spin(i,j,a5)+spin(i,j,a6)))

  end do
 end do
end do

mag=M/(float(N))
E=E*0.5d0 !Check Energy variation 
print*,'Initial energy E, E per spin =',E ,E/float(N)
print*,'Initial magnetization M, M per spin =',M ,mag

!Evolution
open(11,file='ising_model.dat')

do Temp=470,380,-2  !Temp Loop
 T=float(Temp)/100.0d0 
 avg_M=0.0d0; avg_E=0.0d0
 avg_M_N=0.0d0; avg_E_N=0.0d0 !<E>,<M>
 avg_M_2=0.0d0; avg_E_2=0.0d0 !<E^2>,<M^2>
 avg_M_4=0.0d0; avg_E_4=0.0d0 !<E^4>,<M^4>

 do time=1,niter
  do mm=1,L
   do nn=1,L
    do ll=1,L
     call random_number(r);  i=int(r*float(L))+1
     call random_number(r);  j=int(r*float(L))+1
     call random_number(r);  k=int(r*float(L))+1
    
     a1=i+1;a2=i-1;a3=j+1;a4=j-1;a5=k+1;a6=k-1 !neighbours
     !setting boundary conditions
     if(i==L)a1=1
     if(i==1)a2=L
     if(j==L)a3=1
     if(j==1)a4=L
     if(k==L)a5=1
     if(k==1)a6=L
    
     Ei=-J_ising*float(spin(i,j,k)*(spin(a1,j,k)+spin(a2,j,k)+spin(i,a3,k)+spin(i,a4,k)+spin(i,j,a5)+spin(i,j,a6)))
   
     spin(i,j,k)=-spin(i,j,k)

     Ef=-J_ising*float(spin(i,j,k)*(spin(a1,j,k)+spin(a2,j,k)+spin(i,a3,k)+spin(i,a4,k)+spin(i,j,a5)+spin(i,j,a6)))

     dE=Ef-Ei

     if(dE<=0.0d0)then
            E=E+dE
            M=M+(2.0d0*float(spin(i,j,k)))
     else
            u=exp(-dE/(T))
            call random_number(h)
            if(h<u)then
                    E=E+dE
                    M=M+(2.0d0*float(spin(i,j,k)))
            else
                    spin(i,j,k)=-spin(i,j,k)
            end if
     end if       

    end do
   end do
  end do
  !write(11,*)time,M/float(N),E/float(N)
 !end do

 if(time .gt. N_eqm)then
 ! if(mod(time,n_stat) .eq. 0)then
   mag=abs(M)/float(N)
   avg_M=avg_M + mag; avg_E=avg_E + E/float(N)

   avg_M_N= avg_M_N + abs(M); avg_E_N = avg_E_N + E
   avg_M_2 = avg_M_2 + (M*M); avg_E_2 = avg_E_2 + (E*E)
   avg_M_4 = avg_M_4 + (M*M*M*M); avg_E_4 = avg_E_4 + (E*E*E*E) 
 !  end if
 end if
 end do

 avg_M=avg_M/float(niter - N_eqm); avg_E =avg_E/float(niter - N_eqm)
 avg_M_N=avg_M_N/float(niter - N_eqm); avg_E_N =avg_E_N/float(niter - N_eqm)
 avg_M_2=avg_M_2/float(niter-N_eqm); avg_E_2 = avg_E_2/float(niter-N_eqm)
 avg_M_4=avg_M_4/float(niter-N_eqm); 

 cv=(avg_E_2-avg_E_N*avg_E_N)/(T*T); !cv=cv/float(niter-n_eqm)
 chi=(avg_M_2-avg_M_N*avg_M_N)/(T); !chi=chi/float(niter-n_eqm)
 bc=1-(avg_M_4/(3.0d0*(avg_M_2)*(avg_M_2)))
 write(11,*)T,avg_M,avg_E,cv,chi,bc
end do 
close(11)
deallocate(spin)
end program ising_model
