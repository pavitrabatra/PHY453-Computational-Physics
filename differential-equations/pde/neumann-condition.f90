program neumann 
        implicit none

        integer, parameter :: lx=34, ly=34
        real*8 :: old_temp(1:lx,1:ly),temp(1:lx,1:ly)
        integer :: i,j,ii,jj,kk
        real*8 :: dx,dy,shift_temp,prefactor
        integer :: test, counter
        real*8 :: A(ly),B(ly),C(lx),D(lx)


        A=-70.0d0; B=-40.0d0; C=20.0d0; D=-10.0d0
        
        old_temp=0.0d0

        dx=1.0d0;dy=1.0d0
        test=0
        counter=0
        prefactor=(0.5d0*dx*dx*dy*dy)/(dx*dx + dy*dy)

        do
        counter=counter+1
        test=0

        !boundaries without corners
        do j=2,ly-1
          temp(1,j)=0.25d0*(2.0d0*old_temp(2,j)-2.0d0*dx*A(j)+old_temp(1,j-1)+old_temp(1,j+1))
          temp(lx,j)=0.25d0*(2.0d0*old_temp(lx-1,j)+2.0d0*dx*B(j)+old_temp(lx,j-1)+old_temp(lx,j+1))
        end do

        do i=2,lx-1
          temp(i,1)=0.25d0*(old_temp(i+1,1)+old_temp(i-1,1)+2.0d0*old_temp(i,2)-2.0d0*dx*C(i))
          temp(i,ly)=0.25d0*(old_temp(i+1,ly)+old_temp(i-1,ly)+2.0d0*old_temp(i,ly-1)+2.0d0*dx*D(i))
        end do

        !corners
        temp(1,1)=0.5d0*(old_temp(1,2)-dx*C(1)+old_temp(2,1)-dx*A(1))
        temp(1,ly)=0.5d0*(old_temp(1,ly-1)+dx*D(1)+old_temp(2,ly)-dx*A(ly))
        temp(lx,1)=0.5d0*(old_temp(lx-1,1)+dx*B(1)+old_temp(lx,2)-dx*C(lx))
        temp(lx,ly)=0.5d0*(old_temp(lx-1,ly)+dx*B(ly)+old_temp(lx,ly-1)+dx*D(lx))

        !updating
        do jj=2,ly-1
         do ii=2,lx-1
          temp(ii,jj)=0.25*(old_temp(ii-1,jj)+old_temp(ii+1,jj)+old_temp(ii,jj-1)+old_temp(ii,jj+1))
          end do
        end do

        
        !checking
        do jj=1,ly
         do ii=1,lx
          if((abs(temp(jj,ii)-old_temp(jj,ii))).gt.0.00001d0) test=1
          end do
        end do

        if(test.eq.0) exit
        
        shift_temp=2000-temp(1,1)
        temp=temp+shift_temp


        old_temp=temp
        end do

        write(*,*) 'counter', counter

        open(13,file='neumann-20181155.dat')
        do ii=1,lx
         do jj=1,ly
          write(13,*) ii,jj,temp(ii,jj)
          end do
        end do
end program neumann
