program dirichlet
        implicit none

        integer, parameter :: lx=34, ly=34
        real*8 :: old_temp(1:lx,1:ly),temp(1:lx,1:ly)
        integer :: i,j,ii,jj,kk
        real*8 :: dx,dy,bound_temp_1,bound_temp_2,increment_temp,prefactor
        integer :: test, counter

        bound_temp_1=3.7d0
        bound_temp_2=0.4d0
        increment_temp=0.1d0
        
        old_temp=0.0d0
        do i=1,ly
        old_temp(1,i) = bound_temp_1
        old_temp(lx,i) = bound_temp_2
        end do

        do i=2,lx-1
        old_temp(i,1) = old_temp(1,1) - dfloat(i-1)*increment_temp
        old_temp(i,ly) = old_temp(1,ly) - dfloat(i-1)*increment_temp
        end do

        temp=old_temp

        dx=0.05d0;dy=0.05d0 !Will get cancelled since equal so wont really be there in the expression
        test=0
        counter=0
        prefactor=(0.5d0*dx*dx*dy*dy)/(dx*dx + dy*dy)

        do
        counter=counter+1
        test=0

        do jj=2,ly-1
         do ii=2,lx-1
          temp(ii,jj)=0.25*(old_temp(ii-1,jj)+old_temp(ii+1,jj)+old_temp(ii,jj-1)+old_temp(ii,jj+1))
          end do
        end do

        do jj=2,ly-1
         do ii=2,lx-1
          if((abs(temp(jj,ii)-old_temp(jj,ii)).gt.0.0001d0)) test=1
          end do
        end do

        if(test.eq.0) exit
        old_temp=temp
        end do

        write(*,*) 'counter', counter

        open(13,file='dirichlet-20181155.dat')
        do ii=1,lx
         do jj=1,ly
          write(13,*) ii,jj,temp(ii,jj)
          end do
        end do
end program dirichlet
