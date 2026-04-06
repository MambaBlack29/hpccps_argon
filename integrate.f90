subroutine integrate(t, EQMDStep, TotAtom, Mass, Box, Temp, Rcut, Sig, Eps, AtomLabel, TimeStep, r, v, Force, KE, PE)
    use general, only: dp
    use conversions
    implicit none
    integer :: i, j, TotAtom, Step, EQMDStep
    real(kind=dp) :: Box2, t, TimeStep2, ScaleTemp
    real(kind=dp), intent(in) :: Box, Mass, Temp, Rcut, Eps, Sig
    real(kind=dp) :: sumv(3), sumv2, en, Tins, TimeStep
    character(len=5) :: AtomLabel(TotAtom)

    real(kind=dp), intent(out) :: KE, PE
    real(kind=dp), intent(inout) :: r(TotAtom, 3), v(TotAtom, 3), Force(TotAtom, 3)

    Step = int(t/TimeStep)
    TimeStep2 = TimeStep*TimeStep
    sumv = 0.d0
    do i = 1, TotAtom
        r(i, :) = r(i, :) + TimeStep*v(i, :) + 0.5d0*TimeStep2*(Force(i, :)/Mass)      ! r(t+dt)
        !applying scaling 
        if (Step < EQMDStep) then
            call compute_temp(TotAtom, Mass, v, Tins)
            ScaleTemp = dsqrt(Temp/Tins)
        else
            ScaleTemp = 1.d0
        endif
        v(i, :) = v(i, :)*ScaleTemp + 0.5d0*TimeStep*(Force(i, :)/Mass)                      !  v(t + dt/2 )
    enddo

!  if(Step==EQMDSTEP) then         ! shifting momentum to zero  
!    sumv(1)=sum(v(:,1)) 
!    sumv(2)=sum(v(:,2)) 
!    sumv(3)=sum(v(:,3)) 
!    sumv=sumv/real(TotAtom) 
!    do i=1,TotAtom
!      v(i,:)=v(i,:)-sumv 
!    enddo
!  endif

    call force_calc(TotAtom, Box, Rcut, r, Sig, Eps, Force, PE)

    sumv = 0.d0
    sumv2 = 0.d0
    do i = 1, TotAtom
        v(i, :) = v(i, :) + 0.5d0*TimeStep*(Force(i, :)/Mass)     ! v(t+dt)   here force is updated force at t+dt
        sumv = sumv + v(i, :)
        sumv2 = sumv2 + dot_product(v(i, :), v(i, :))
    end do
    sumv2 = Mass*sumv2
    Tins = (sumv2)/real(3.d0*TotAtom - 3.d0)
    KE = sumv2/real(2.d0)

 ! checking whether coordinates are outside the box 

    do i = 1, TotAtom
        do j = 1, 3
            if (r(i, j) > Box) then
                r(i, j) = r(i, j) - Box
            elseif (r(i, j) < 0.d0) then
                r(i, j) = r(i, j) + Box
            endif
        enddo
    enddo

! write(*,*) "VELOCITIES"  
! write(206,"(I8,6F20.12)") Step,Mass*sumv 
! write(*,*)  "Temperature", Tins*TempConv


! writing Trajectory  
! open(unit=201,file='md.traj',action='write') 
! write(201,"(I5)") TotAtom 
! write(201,*) Step+1 
! write(5000,*) 
! write(*,*) "STEP = ",Step+1
! do i=1,TotAtom
!   write(201,"(a5,2x,3F12.4)") AtomLabel(i),r(i,:)*LengthConv
! enddo
  return
end subroutine integrate

