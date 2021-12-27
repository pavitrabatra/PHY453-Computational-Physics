program variational_principle

        implicit none
        integer,parameter:: dp=selected_real_kind(14,200)
        real(dp)::V0,a,b
        integer::n,pw
        real(dp) :: x, dx, modulus, p
        real(dp),allocatable::E(:),H(:,:),k(:)
        real(dp),allocatable::work(:)
        real(dp),parameter:: pi=3.14159265358979_dp
        complex(dp)::psi
        integer::i,j,lwork,info,t

        write(*,"('Give V0')")
        read(*,*) V0

        write(*,"('Give b')")
        read(*,*) b

        write(*,"('Give a')")
        read(*,*) a 

        write(*,"('Give n')")
        read(*,*) n

        pw=2*n+1

        allocate (E(pw), H(pw,pw), k(pw), work(3*pw))
       
        H(:,:)=0.0_dp
        k(1)=0.0_dp
                
        do i=2,pw-1,2
        k(i)=(i/2)*2.0_dp*pi/a
        k(i+1)=-(i/2)*2.0_dp*pi/a
        end do

          do i=1,pw
          do j=1,pw
          if (i==j) then
                  H(i,j)=k(i)**2 - (V0/a)*b
          else
                  H(i,j)=-(V0/a) * sin((k(j)-k(i))*b/2.0_dp) /(k(j)-k(i))*2.0_dp
          end if
          end do
          end do
        
          lwork =3*pw
          call dsyev ('V','U', pw, H, pw, E, work, lwork, info)
          if(info /= 0) stop 'Diagonalisation failed'

          open(69,file='energy-20181155.dat',action="write")
          dx=0.01_dp
          modulus=0.0_dp
          t= nint(a/2.0_dp/dx)
          do i=-t,t
          x=dx*i
          psi=0.0_dp
           do j=1,pw
           psi=psi+H(j,1)*exp((0.0,1.0)*k(j)*x)/sqrt(a)
           end do
          p=psi*conjg(psi)
          modulus =modulus+p*dx
          end do
          write(69,*) E(1)
          deallocate(E,H,work,k)
          close(69)
end program variational_principle






        
