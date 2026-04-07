subroutine new_verlet(TotAtom,Box,Rcut,r, vlist, nvlist)
    use general, only: dp, atom1, atom2
    implicit none
    integer, intent(in) :: TotAtom
    real(kind=dp), intent(in) :: Box, Rcut
    real(kind=dp), intent(in) :: r(TotAtom, 3)
    integer, intent(out) :: vlist(TotAtom, 200), nvlist(TotAtom)
    real(kind=dp) :: r2, dr(3), R2cut

    R2cut = Rcut*Rcut

    nvlist = 0
    vlist = 0
    do atom1 = 1, TotAtom - 1
        do atom2 = atom1 + 1, TotAtom
            dr = r(atom1, :) - r(atom2, :)
            dr = dr - Box*anint(dr/Box)
            r2 = dot_product(dr, dr)
            if (r2 <= R2cut) then
                nvlist(atom1) = nvlist(atom1) + 1
                nvlist(atom2) = nvlist(atom2) + 1
                vlist(atom1, nvlist(atom1)) = atom2
                vlist(atom2, nvlist(atom2)) = atom1
            endif
        enddo
    enddo
end subroutine new_verlet
