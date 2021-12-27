program numerov
 implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  integer :: grid, i, icl
  integer :: nodes, hnodes, ncross, j, niter
  real(dp) :: xmax, dx, c1, modulus, c2, jump, c3
  real(dp) :: Emax, Emin, e
  real(dp), allocatable :: x(:), y(:), p(:), V(:), f(:)

  write(*,"('Enter xmax')") 
  read (*,*) xmax
  write(*,"('Enter number of grid points')")
  read (*,*) grid
  
  
  allocate ( x(0:grid), y(0:grid), p(0:grid), V(0:grid), f(0:grid) )
  dx =  xmax/grid 
  c1=dx*dx/12.0_dp

  do i = 0, grid
     x(i) = float(i) * dx
     V(i) = 0.5_dp * x(i)*x(i)
  end do

  open (20,file='numerov-20181155.dat')

  search_loop: do 
     write(*,"('Enter number of nodes(-1 to stop)')") 
     read (*,*) nodes

     if (nodes < 0) then
        close(20)
        deallocate ( f, V, p, y, x )
        stop 
     end if

     Emax=maxval (V(:))
     Emin=minval (V(:))

     write(*,"('Trial energy (0=search with bisection)')") 
     read (*,*) e

     if ( e == 0.0_dp ) then
        e = 0.5_dp * (Emin + Emax)
        niter = 1000
     else
             niter = 1
     endif
     
     iterate: do j = 1, niter
        f(0)=c1*(2.0_dp*(V(0)-e))
        icl=-1
        do i=1,grid
           f(i)=c1*2.0_dp*(V(i)-e)
           if ( f(i) == 0.0_dp) f(i)=1.d-20
           if ( f(i) /= sign(f(i),f(i-1)) ) icl=i
        end do
        
        if (icl >= grid-2) then
           deallocate ( f, V, p, y, x )
           print *, 'last change of sign too far'
           stop 1
        else if (icl < 1) then
           deallocate ( f, V, p, y, x )
           print *, 'no classical turning point'
           stop 1
        end if
        
        f = 1.0_dp - f
        y = 0.0_dp

        hnodes = nodes/2

        if (2*hnodes == nodes) then
           y(0) = 1.0_dp
           y(1) = 0.5_dp*(12.0_dp-10.0_dp*f(0))*y(0)/f(1)
        else
           y(0) = 0.0_dp
           y(1) = dx
        end if

        ncross=0
        do i =1,icl-1
           y(i+1)=((12.0_dp-10.0_dp*f(i))*y(i)-f(i-1)*y(i-1))/f(i+1)
           if ( y(i) /= sign(y(i),y(i+1)) ) ncross=ncross+1
        end do
        c3 = y(icl)

        if (2*hnodes == nodes) then
           ncross = 2*ncross
        else
           ncross = 2*ncross+1
        end if

        if ( niter > 1 ) then
           if (ncross /= nodes) then
              if ( j == 1) &
                   print '("Bisection Energy Nodes Discontinuity")'
              print *, j, e, ncross
              if (ncross > nodes) then
                 Emax = e
              else 
                 Emin = e
              end if
              e = 0.5_dp * (Emax+Emin)
              cycle
           end if
        else
           print *, e, ncross, nodes
        end if

        y(grid) = dx
        y(grid-1) = (12.0_dp-10.0_dp*f(grid))*y(grid)/f(grid-1)
        modulus = 1.0d100
        do i = grid-1,icl+1,-1
           y(i-1)=((12.0_dp-10.0_dp*f(i))*y(i)-f(i+1)*y(i+1))/f(i-1)
       
           if ( abs(y(i-1)) > modulus ) then
              y(i-1:grid) = y(i-1:grid) / modulus
           endif
        end do

        c3 = c3/y(icl)
        y(icl:) = y(icl:)*c3

        modulus = (2.0_dp*dot_product (y, y) - y(0)*y(0)) * dx 
        y = y / sqrt(modulus)

        if ( niter > 1 ) then
           jump = ( y(icl+1) + y(icl-1) - (14.0_dp-12.0_dp*f(icl))*y(icl) ) / dx
           print *, j, e, nodes, jump
           if (jump*y(icl) > 0.0_dp) then
              Emax = e
           else
              Emin = e
           endif
           e = 0.5_dp * (Emax+Emin)
           if ( Emax-Emin < 1.d-10) exit iterate
        endif
        
     end do iterate

     modulus = 0.0_dp
     p(icl:) = 0.0_dp
     do i=0,icl
        c2 = (e - x(i)**2/2.0_dp)
        if ( c2 > 0.0_dp) then
           p(i) = 1.0_dp/sqrt(c2)
        else
           p(i) = 0.0_dp
        end if
        modulus = modulus + 2.0_dp*dx*p(i)
     enddo

     modulus = modulus - dx*p(0)
     p(:icl-1) = p(:icl-1)/modulus

     !x<0
     do i=grid,1,-1
        if ( abs(y(i)) < 1.0D-50 ) y(i) = 0.0_dp
        write (20,*) -x(i), (-1)**nodes*y(i), y(i)*y(i), p(i), V(i)
     enddo

     !x>0
     do i=0,grid
        write (20,*) x(i), y(i), y(i)*y(i), p(i), V(i)
     enddo
     
  end do search_loop
  print*,e
end program numerov
