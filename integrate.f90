subroutine integrate(t, EQMDStep, TotAtom, Mass, Box, Temp, Rcut, Sig, Eps, AtomLabel, TimeStep, r, v, Force, KE, PE, vlist, nvlist, rlist, rskin)
    use general, only: dp
    use conversions
    implicit none
    integer :: i, j, TotAtom, Step, EQMDStep
    real(kind=dp) :: t, TimeStep2, ScaleTemp
    real(kind=dp), intent(in) :: Box, Mass, Temp, Rcut, Eps, Sig, rskin
    real(kind=dp) :: sumv(3), sumv2, Tins, TimeStep
    character(len=5) :: AtomLabel(TotAtom)

    real(kind=dp), intent(out) :: KE, PE
    real(kind=dp), intent(inout) :: r(TotAtom, 3), v(TotAtom, 3), Force(TotAtom, 3)

    ! verlet list variables
    integer, intent(inout) :: vlist(TotAtom, 500), nvlist(TotAtom)
    real(kind=dp), intent(inout) :: rlist(TotAtom, 3)

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

    if ( maxval( sum( ( (r - rlist) - Box*anint((r - rlist)/Box) )**2, dim=2 ) ) > (0.5d0*rskin)**2 ) then
            call new_verlet(TotAtom, Box, Rcut, r, vlist, nvlist, rlist, rskin)
    end if

    call force_calc(TotAtom, Box, Rcut, r, Sig, Eps, Force, PE, vlist, nvlist)

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

  return
end subroutine integrate

