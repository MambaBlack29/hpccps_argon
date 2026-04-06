subroutine compute_temp(TotAtom,Mass,v,Tins)
  use general, only: dp
  use conversions 
  implicit none
  integer :: TotAtom,i
  real(kind=dp),intent(in) :: v(TotAtom,3),Mass
  real(kind=dp),intent(out) :: Tins
  real(kind=dp) :: sumv2 


  sumv2=0.d0
  do i=1,TotAtom
    sumv2=sumv2+dot_product(v(i,:),v(i,:))
  enddo 
  
  sumv2=Mass*sumv2
  Tins=sumv2/real(3.d0*TotAtom-3.d0)

end subroutine 
