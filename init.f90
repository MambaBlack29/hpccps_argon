subroutine initialize(TotAtom,CoorFileName,Temp,Mass,Box,r,v,AtomLabel)
    use general, only: dp
    use conversions
! use ifport
    implicit none

    integer :: CheckAtom
    character(len=300), intent(in) :: CoorFileName
    integer, intent(in) :: TotAtom
    real(kind=8), intent(in) :: Mass
    real(kind=dp), intent(inout) :: r(TotAtom, 3)
    real(kind=dp), intent(out) :: v(TotAtom, 3)

    character(len=5), intent(out) :: AtomLabel(TotAtom)
    integer :: i, j

    real(kind=dp) :: sumv(3), sumv2, ScaleTemp, Box, Temp

    open(unit=2, file=CoorFileName, action='read')
    read(2, *) CheckAtom
    read(2, *)
    if (CheckAtom /= TotAtom) then
        write(*, *) "Check Input Files --- Total No. of Atoms"
        stop
    endif

! reading coordinates  and stored in x array 

    do i = 1, TotAtom
        read(2, *) AtomLabel(i), r(i, :)
        r(i, :) = r(i, :)/LengthConv               ! in au
    enddo

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


!assign vocties 

    sumv = 0.d0
    do i = 1, TotAtom
        v(i, :) = (/rand() - 0.5d0, rand() - 0.5d0, rand() - 0.5d0/)
        sumv = sumv + v(i, :)
    enddo

    sumv = sumv/real(TotAtom)

! zeroing the linear momentum  and scaling 

  sumv2 = 0.d0
  do i = 1, TotAtom
    v(i, :) = (v(i, :) - sumv)
    sumv2 = sumv2 + dot_product(v(i, :), v(i, :))
  enddo

  sumv2 = sumv2*Mass

  ScaleTemp = dsqrt((3.d0*TotAtom - 3.d0)*Temp/(sumv2))
  v = v*ScaleTemp             ! velocity scaling

 !recalculate temperature  ---  for checking 

    sumv2 = 0.d0
    sumv = 0.d0
    do i = 1, TotAtom
        sumv2 = sumv2 + dot_product(v(i, :), v(i, :))
        sumv = sumv + v(i, :)
    enddo

    write(5000, "(a,4F12.3)") "Temperature is set to ", Mass*sumv2*TempConv/real(3.d0*TotAtom - 3.d0)

! write(*,*) "VELOCITIES"  
! write(*,*) "INIT VELOCITIES"
! write(*,"(3E20.12)") (v(i,:),i=1,2) 

  return
end subroutine initialize

