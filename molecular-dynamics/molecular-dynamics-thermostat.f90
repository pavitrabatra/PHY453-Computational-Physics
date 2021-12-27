program md
        implicit none
        real*8::V,KE,TE
        real*8,parameter::dt=0.005
        integer,parameter::n=2000
        integer::i,j,nn,niter,jj,kk
        real*8::avg_vx,avg_vy,avg_vz
        real*8,parameter::rc=2.5d0,rs=2.0d0*rc
        real*8::p,q,m,dt2by2,dtby2,p2,q2,m2 
        real*8::dx,dy,dz,pot
        real*8::r,dist,force,s,px,py,pz
        integer,dimension(n-1)::no_ngb
        integer,dimension(n-1,n)::n_list
        real*8,dimension(3*n)::pos,vel,F,old_F
        
        open(unit=21,file='result2.dat')
        open(unit=22,file='initial2.dat')

        write(*,"('Number of iterations')")
        read(*,*) niter

        dt2by2=0.5d0*dt*dt
        dtby2=0.5d0*dt

        call random_seed

        !!Initializing
        
        !Position
        do i=1,n-1
        pos(3*i-2)= int(i/100) 
        pos(3*i-1)= int((i-pos(3*i-2)*100)/10)
        pos(3*i)=i-pos(3*i-2)*100-pos(3*i-1)*10
        end do
        pos(3*n-2)= 20 
        pos(3*n-1)= 20
        pos(3*n)= 20
        write(22,*) pos

        !Neighbour List
        call nbr(pos,no_ngb,n_list)
        
        !Force
        call calc_F(pos,F,V)

        !Velocity
        do i=1,n
        call random_number(p)
        call random_number(q)
        call random_number(m)
        vel(3*i-2)=sqrt(12.0d0)*(p-0.5)
        vel(3*i-1)=sqrt(12.0d0)*(q-0.5)
        vel(3*i)=sqrt(12.0d0)*(m-0.5)
        end do

        avg_vx=0.0d0
        avg_vy=0.0d0
        avg_vz=0.0d0
        
        do i=1,n
        avg_vx=avg_vx+vel(3*i-2)
        avg_vy=avg_vy+vel(3*i-1)
        avg_vz=avg_vz+vel(3*i)
        end do

        avg_vx=avg_vx/dfloat(n)
        avg_vy=avg_vy/dfloat(n)
        avg_vz=avg_vz/dfloat(n)

        do i=1,n
        vel(3*i-2)=avg_vx-vel(3*i-2)
        vel(3*i-1)=avg_vy-vel(3*i-1)
        vel(3*i)=avg_vz-vel(3*i)
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!
        
        !!Updating
        
        do nn=1,niter
        
        !Position
        do i=1,n
         pos(3*i-2)=pos(3*i-2)+vel(3*i-2)*dt+0.5*F(3*i-2)*dt*dt
         pos(3*i-1)=pos(3*i-1)+vel(3*i-1)*dt+0.5*F(3*i-1)*dt*dt
         pos(3*i)=pos(3*i)+vel(3*i)*dt+0.5*F(3*i)*dt*dt
        end do
        
        !Neighbour List
        if(mod(nn,100) .eq. 0) call nbr(pos,no_ngb,n_list)

        !Force
        old_F = F
        call calc_F(pos,F,V)

        !Velocity
        do i=1,n 
         vel(3*i-2)=vel(3*i-2)+(dt*0.5d0*(old_F(3*i-2)+F(3*i-2)))
         vel(3*i-1)=vel(3*i-1)+(dt*0.5d0*(old_F(3*i-1)+F(3*i-1)))
         vel(3*i)=vel(3*i)+(dt*0.5d0*(old_F(3*i)+F(3*i)))
        end do
         
        !Kinetic Energy
        KE=0.0d0
        
        do i=1,n
         KE=KE+vel(3*i-2)*vel(3*i-2)
         KE=KE+vel(3*i-1)*vel(3*i-1)
         KE=KE+vel(3*i)*vel(3*i)
        end do
        
        KE=0.50d0*KE
        
        !Thermostat
        s=sqrt((3/2)/(KE/n))
        if (mod(nn,100)==0) then
        vel=s*vel
        KE=1.50d0*n
        end if
        
        !Momentum
        px=0.0d0
        py=0.0d0
        pz=0.0d0
        do i=1,n
        px=px+vel(3*i-2)
        py=py+vel(3*i-1)
        pz=pz+vel(3*i)
        end do

        !Total Energy
        TE=0.0d0
        TE=(KE+V)/n
        
        write(21,*) nn,TE,px,py,pz
        print*, nn
        end do

end program md

real*8 function force(r) result(res)
        implicit none
        real*8,parameter::rc=2.5d0,sigma=1.0d0
        real*8,parameter::sigma6=sigma**6,sigma12=sigma**12
        real*8,parameter::fc=4.0d0*((12.0d0*sigma12/(rc**13))-(6.0d0*sigma6/(rc**7)))
        real*8::r

        res=4.0d0*((12.0d0*sigma12/(r**13))-(6.0d0*sigma6/(r**7)))-fc
end function


real*8 function pot(r) result(res)
        implicit none
        real*8,parameter::rc=2.5d0,sigma=1.0d0
        real*8,parameter::sigma6=sigma**6,sigma12=sigma**12
        real*8,parameter::fc=4.0d0*((12.0d0*sigma12/(rc**13))-(6.0d0*sigma6/(rc**7)))
        real*8,parameter::vc=fc*rc+4.0d0*(((sigma/rc)**12)-((sigma/rc)**6))
        real*8::r
        
        res=4.0d0*(((sigma/r)**12)-((sigma/r)**6))+fc*r-vc
end function

real*8 function dist(dx,dy,dz) result(res)
        implicit none
        real*8::dx,dy,dz
        integer::lx=20,ly=20,lz=20
        real*8::llx,lly,llz
        llx=dfloat(lx)
        lly=dfloat(ly)
        llz=dfloat(lz)
         if (abs(dx) .gt. lx/2.0d0) dx=(llx-abs(dx))*(-dx/abs(dx))
         if (abs(dy) .gt. ly/2.0d0) dy=(lly-abs(dy))*(-dy/abs(dy))
         if (abs(dz) .gt. lz/2.0d0) dz=(llz-abs(dz))*(-dz/abs(dz))
         res =sqrt(dx*dx+dy*dy+dz*dz)
end function

subroutine calc_F(pos,F,V)
        implicit none 
        integer::i,jj,kk
        integer,parameter::n=2000
        real*8::dx,dy,dz
        real*8::r,dist
        real*8,dimension(3*n),intent(out)::F 
        real*8,parameter::rc=2.5d0
        real*8,dimension(3*n),intent(in)::pos
        real*8::force,pot
        integer::aa
        integer,dimension(n-1,n)::n_list
        integer,dimension(n-1)::no_ngb
        real*8,intent(out)::V

        call nbr(pos,no_ngb,n_list)

        V=0.0d0
        F=0
        do jj=1,n-1
         do kk=1,no_ngb(jj)
         
         aa=n_list(jj,kk)

         dx=pos(3*jj-2)-pos(3*aa-2)
         dy=pos(3*jj-1)-pos(3*aa-1)
         dz=pos(3*jj)-pos(3*aa)
         r=dist(dx,dy,dz)

         if(r<rc) then
                 F(3*jj-2)=F(3*jj-2)+force(r)*(dx/r)
                 F(3*jj-1)=F(3*jj-1)+force(r)*(dy/r)
                 F(3*jj)=F(3*jj)+force(r)*(dz/r)

                 F(3*aa-2)=F(3*aa-2)+force(r)*(-dx/r)
                 F(3*aa-1)=F(3*aa-1)+force(r)*(-dy/r)
                 F(3*aa)=F(3*aa)+force(r)*(-dz/r)
                 V=V+pot(r)
         end if
         end do
        end do
end subroutine

subroutine nbr(pos,no_ngb,n_list)
        implicit none
        integer::jj,kk
        integer,parameter::n=2000
        real*8::dx,dy,dz,r,dist
        real*8,parameter::rc=2.5d0,rs=2.0d0*rc
        real*8,dimension(3*n),intent(in)::pos
        integer,dimension(n-1),intent(out)::no_ngb
        integer,dimension(n-1,n),intent(out)::n_list

        no_ngb=0
        n_list=0
        do jj=1,n-1
         do kk=jj+1,n
         dx=pos(3*kk-2)-pos(3*jj-2)
         dy=pos(3*kk-1)-pos(3*jj-1)
         dz=pos(3*kk)-pos(3*jj)
         r=dist(dx,dy,dz)

         if (r<rs) then
         no_ngb(jj)=no_ngb(jj)+1
         n_list(jj,no_ngb(jj))=kk
         end if
         end do
        end do
end subroutine
