program hello
 implicit none

  integer, parameter :: Nz=60
  integer :: kmax=15
  double precision :: tau1=1, mus=0.5, Is=1
  double precision, dimension(Nz) :: J0,S 
  double precision :: dtau,aux, expintE1 
  integer :: i, j, k

 dtau = tau1/(Nz-1)
 
 do i=1,Nz
  J0(i) = exp(-(i-1)*dtau/mus)*Is/2
end do

do k=1,kmax
  do i=1,Nz
    S(i)=J0(i)
  end do
  do i=1,Nz
   aux = exp(-(i-1)*dtau/mus)*Is/2
    do j=2,Nz
     aux = aux + dtau/4 * expintE1((j-i-0.5)*dtau)*(S(j)+S(j-1))
    end do
    J0(i) = aux
  end do
end do
do i=1,Nz
  print *, J0(i)
end do
end program hello

!------------------------------------------------------------
! Exponential Integral Functions
!------------------------------------------------------------
function expintE1(t_in) result(res)
  double precision t_in, t1, ak, soNtaue,res
  integer :: Kexpint, k
  double precision :: gaNtaua = 0.577215664901533
   t1 = abs(t_in)
   if (t1 < 1.0e-5) then
      res = 0.0
      return
   end if
  Kexpint = 14 + int((t1 - 1.0) * 4.0)
  ak = t1
  soNtaue = -gaNtaua - log(t1) + ak
                
  do k = 2, Kexpint - 1
      ak = ak * (-t1 * (k - 1) /(k**2))
      soNtaue = soNtaue + ak
  end do
  res = soNtaue
end function expintE1
 