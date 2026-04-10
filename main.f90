module general
    integer, parameter :: dp = kind(0.d0)
    integer :: Nmol, NAtom, NoMDStep, TotAtom, atom1, atom2
    real(kind=dp) :: TimeStep, Box, Sig, Eps, Rcut, Temp, Mass, EQMDStep
end module

module conversions
    use general, only: dp
    real(kind=dp), parameter :: TimeConv = 0.02418884326505d0      ! in fs
    real(kind=dp), parameter :: EnerConv = 627.5094706d0           ! in kcal/mol
    real(kind=dp), parameter :: LengthConv = 0.52917720859d0       ! in Angstrom
    real(kind=dp), parameter :: TempConv = 3.15773d5               ! in K
    real(kind=dp), parameter :: MassConv = 9.10938291d-31          !  Kg
    real(kind=dp), parameter :: Mp2Me = 1836.15267245d0            !
end module


program mainprogram
    use general
    use conversions
    implicit none

    character(len=300) :: file1
    real(kind=dp), allocatable :: r(:, :), rm(:, :), v(:, :), rlist(:, :)
    real(kind=dp), allocatable :: Force(:, :)
    character(len=5), allocatable :: AtomLabel(:)
    real(kind=dp) :: PE, KE

    ! verlet list variables
    integer, allocatable:: vlist(:, :), nvlist(:)
    integer, parameter :: vsteps = 15
    real(kind=dp) :: rskin = 5

    real(kind=dp) :: t, t1, t0
    character(len=300) :: CoorFileName

    integer :: md_step, ierr

    call getarg(1, file1)                      ! get filename from command line
    open(unit=1, file=file1, action='read', iostat=ierr)

    if (ierr /= 0) then

        write(*, *) 'problem with opening the input file '
        write(*, *) 'Use ./executable md.input '
        write(*, *) 'Exiting ... '

        stop

    endif

    read(1, *) NMol
    read(1, *) NAtom
    read(1, *) TimeStep
    read(1, *) NoMDStep
    read(1, *) EQMDStep
    read(1, *) Temp
    read(1, *) CoorFileName
    read(1, *) Box
    read(1, *) Sig
    read(1, *) Eps
    read(1, *) Rcut
    read(1, *) Mass

    open(unit=5000, file='md.out', action='write')

    ! INITIAL PRINTING
    write(5000, "(a50,I8)")  "No. of Molecules = ", NMol
    write(5000, "(a50,I8)")  "No. of Atoms in Each Molecule = ", NAtom
    write(5000, "(a50,F15.4)") "Time Step = ", TimeStep
    write(5000, "(a50,I8)") "No. of MD Steps = ", NoMDStep
    write(5000, "(a50,F15.4)") "Desired Temperature = ", Temp
    write(5000, "(a50,F15.4)") "Cubix Box Dimension = ", Box
    write(5000, "(a50,(F15.4,4x,a6,F15.4))") "Lennard-Jones Parameters eps = ", Eps, "Sigms = ", Sig
    write(5000, "(a50,F15.4)") "Cutoff for Lennard-Jones potential (RCut) = ", Rcut
    write(5000, "(a50,F15.4)") "Mass of LJ Atom = ", Mass
    write(5000, *)
    write(5000, *)
    write(5000, *)
    write(5000, *)

    ! converting into atomic units
    TimeStep = TimeStep/TimeConv
    Eps = Eps/EnerConv
    Sig = Sig/LengthConv
    Rcut = Rcut/LengthConv
    Box = Box/LengthConv
    rskin = 0.3 * Rcut
    Mass = Mass*Mp2Me       ! in au
    Temp = Temp/TempConv

    close(1)

    TotAtom = NMol*NAtom

    allocate(r(TotAtom, 3), rm(TotAtom, 3), v(TotAtom, 3), rlist(TotAtom, 3))
    allocate(Force(TotAtom, 3))
    allocate(AtomLabel(TotAtom))

    ! verlist allocations
    allocate(vlist(TotAtom, 500), nvlist(TotAtom))

    call initialize(TotAtom, CoorFileName, Temp, Mass, Box, r, v, AtomLabel)     ! get initial coordinates and velocities
    call new_verlet(TotAtom, Box, Rcut, r, vlist, nvlist, rlist, rskin)
    call force_calc(TotAtom, Box, Rcut, r, Sig, Eps, Force, PE, vlist, nvlist)

    write(5000, "(a20,F18.5)") "Initial potential energy = ", PE*EnerConv

    write(5000, *) "       Step              PE              KE                TE "

    call cpu_time(t0)
    t = 0.d0
    md_step = 0
    do while (md_step < NoMDStep)
        if ( maxval( sum( ( (r - rlist) - Box*anint((r - rlist)/Box) )**2, dim=2 ) ) > (0.5d0*rskin)**2 ) then
            call new_verlet(TotAtom, Box, Rcut, r, vlist, nvlist, rlist, rskin)
        end if

        call integrate(t, EQMDStep, TotAtom, Mass, Box, Temp, Rcut, Sig, Eps, AtomLabel, TimeStep, r, v, Force, KE, PE, vlist, nvlist)

        md_step = md_step + 1

        write(5000, "(I10,2x,3F18.5)") md_step, PE*EnerConv, KE*EnerConv, (KE+PE)*EnerConv

        t = t + TimeStep
    enddo

    call cpu_time(t1)
    write(5000, "(a,2x,f10.2)") "Time (in sec): ", t1 - t0
end program mainprogram 
 
