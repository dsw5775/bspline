module bspline
!  reference:
!    A smooth particle mesh Ewald method, J. Chem. Phys. 103, 8577 (1995);
!
!  bspline_1d and bspline_2d are kept in simple form.
!  bspline_3d has some extra optimization and no other extra for the higher-dimention ones

contains

recursive function Mn(n, u) result(answer)
  implicit none

  integer :: n
  real(kind = 8) :: u, answer
  
  if(n < 2) then
    write(*,*) "Order of B-spline is less than 2"
    stop
  endif

  if(n == 2) then
    if(u < 0.0 .or. u > 2.0 ) then
      answer = 0.0
      return
    endif

    if(u < 1.0) then
      answer = u
      return
    else
      answer = 2.0 - u
      return
    endif
  else
    answer = u/(n-1)*Mn(n-1, u) + (n-u)/(n-1)*Mn(n-1, u-1)
    return
  endif
end function Mn

function dMn(n, u)
  implicit none

  real(kind = 8) :: dMn
  integer :: n
  real(kind = 8) :: u

  dMn = Mn(n-1, u) - Mn(n-1, u-1)
end function dMn

function dMn2(n, u)
  implicit none

  real(kind = 8) :: dMn2
  integer :: n
  real(kind = 8) :: u
  
  dMn2 = Mn(n-2, u) - 2*Mn(n-2, u-1) + Mn(n-2, u-2)
end function dMn2

subroutine bspline_1d(data, min_x1, max_x1, &
  num_x1, x1, order, r_value, &
  r_d1, f_d1, r_d2, f_d2)
!  num_x1: number of data points
!  min_x1: x-coordinate of the first data point
!  max_x1: x-coordinate of the last data point
!  order: B-spline order
!  r_value: fitting value
!  r_d1: first-derivative value
!  f_d1: flag for calculating first-derivative value
!  r_d2: second-derivative value
!  f_d2: flag for calculating second-derivative value

  implicit none

  integer :: num_x1
  real(kind = 8), dimension(num_x1) :: data
  real(kind = 8) :: min_x1, max_x1, x1
  integer :: order
  logical :: f_d1, f_d2
  real(kind = 8), intent(out) :: r_value, r_d1, r_d2

  integer :: numbin_x1, i, x1grid0, x1grid1
  real(kind = 8) :: x1scaled, datapoint
  real(kind = 8) :: value, d1, d2

  numbin_x1 = num_x1 - 1
  x1scaled = (x1 - min_x1) / (max_x1 - min_x1) * numbin_x1 + 1 ! different from C++ code as lower bound is 1

  x1grid0 = int(x1scaled - order/2.0) + 1
  x1grid1 = int(x1scaled + order/2.0)

  if(x1grid1 < 1 .or. x1grid0 > numbin_x1) then
    write(*,*) "out of range (bspline_1d)" 
    stop
  endif

  ! If index is nearby the boundary, adjust accordingly.
  if(x1grid1 > num_x1) x1grid1 = num_x1
  if(x1grid0 < 1) x1grid0 = 1

  ! Mn(order, ...) has domain [0, order] with center (peak) at order/2.0,
  ! so later on the calculated relative distance to the grid point (x1scaled - i)
  ! need to be shifted to correspond to the peak.
  x1scaled = x1scaled + order/2.0
  
  value = 0.0
  d1 = 0.0
  d2 = 0.0

  do i = x1grid0, x1grid1
    datapoint = data(i)

    value = value + datapoint * Mn(order, x1scaled - i)

    if(f_d1) &
      d1 = d1 + datapoint * dMn(order, x1scaled - i) &
        * (numbin_x1 / (max_x1 - min_x1))

    if(f_d2) &
      d2 = d2 + datapoint * dMn2(order, x1scaled - i) &
        * (numbin_x1 / (max_x1 - min_x1)) * (numbin_x1 / (max_x1 - min_x1))
  enddo

  r_value = value
  r_d1 = d1
  r_d2 = d2

end subroutine bspline_1d
  
subroutine bspline_2d(data, &
  min_x1, max_x1, num_x1, x1, &
  min_x2, max_x2, num_x2, x2, &
  order, r_value, r_d1, r_d2, f_d)
!  r_value: fitting value
!  r_d1: partial first-derivative value along x1 direction
!  r_d2: partial first-derivative value along x2 direction (not to be confused with the one in bspline_1d)
!  f_d1: flag for calculating partial first-derivative value

  implicit none

  integer :: num_x1, num_x2
  real(kind = 8), dimension(num_x2, num_x1) :: data
  real(kind = 8) :: min_x1, max_x1, x1, min_x2, max_x2, x2
  integer :: order
  logical :: f_d
  real(kind = 8), intent(out) :: r_value, r_d1, r_d2

  integer :: numbin_x1, numbin_x2, i, j, x1grid0, x1grid1, x2grid0, x2grid1
  real(kind = 8) :: x1scaled, x2scaled, datapoint
  real(kind = 8) :: value, d1, d2, Mn1, Mn2, dMn1, dMn2
  
  numbin_x1 = num_x1 - 1;
  numbin_x2 = num_x2 - 1;

  x1scaled = (x1-min_x1)/(max_x1-min_x1)*numbin_x1 + 1
  x2scaled = (x2-min_x2)/(max_x2-min_x2)*numbin_x2 + 1

  x1grid0 = int(x1scaled - order/2.0) + 1
  x1grid1 = int(x1scaled + order/2.0)

  x2grid0 = int(x2scaled - order/2.0) + 1
  x2grid1 = int(x2scaled + order/2.0)

  if(x1grid1 < 1 .or. x1grid0 > num_x1 .or. x2grid1 < 1 .or. x2grid0 > num_x2) then
    write(*,*) "out of range (bspline_2d)" 
    stop
  endif

  if(x1grid1 > num_x1) x1grid1 = num_x1
  if(x1grid0 < 1) x1grid0 = 1

  if(x2grid1 > num_x2) x2grid1 = num_x2
  if(x2grid0 < 1) x2grid0 = 1

  x1scaled = x1scaled + order/2.0
  x2scaled = x2scaled + order/2.0 

  value = 0.0
  d1 = 0.0
  d2 = 0.0

  do i = x1grid0, x1grid1
    do j = x2grid0, x2grid1
      datapoint = data(j, i)

      Mn1 = Mn(order, x1scaled - i)
      Mn2 = Mn(order, x2scaled - j)

      value = value + datapoint*Mn1*Mn2

      if(f_d) then
        dMn1 = dMn(order, x1scaled - i)
        dMn2 = dMn(order, x2scaled - j)

        d1 = d1 + datapoint*dMn1*Mn2 & 
          *(numbin_x1/(max_x1-min_x1))

        d2 = d2 + datapoint*Mn1*dMn2 &
          *(numbin_x2/(max_x2-min_x2))
      endif
    enddo
  enddo

  r_value = value
  r_d1 = d1
  r_d2 = d2

end subroutine bspline_2d

subroutine bspline_3d(data, &
    min_x1, max_x1, num_x1, x1, &
    min_x2, max_x2, num_x2, x2, &
    min_x3, max_x3, num_x3, x3, &
    order, r_value, r_d1, r_d2, r_d3, f_d)
  
    implicit none
  
    integer :: num_x1, num_x2, num_x3
    real(kind = 8), dimension(num_x3, num_x2, num_x1) :: data
    real(kind = 8) :: min_x1, max_x1, x1, min_x2, max_x2, x2, &
      min_x3, max_x3, x3
    integer :: order
    logical :: f_d
    real(kind = 8), intent(out) :: r_value, r_d1, r_d2, r_d3
  
    integer :: numbin_x1, numbin_x2, numbin_x3, i, j, k, &
      x1grid0, x1grid1, x2grid0, x2grid1, x3grid0, x3grid1
    real(kind = 8) :: x1scaled, x2scaled, x3scaled, datapoint
    real(kind = 8) :: value, d1, d2, d3, Mn1, Mn2, Mn3, dMn1, dMn2, dMn3
    real(kind = 8) :: Mn1a, Mn1b, Mn2a, Mn2b, Mn3a, Mn3b

    numbin_x1 = num_x1 - 1;
    numbin_x2 = num_x2 - 1;
    numbin_x3 = num_x3 - 1;
  
    x1scaled = (x1-min_x1)/(max_x1-min_x1)*numbin_x1 + 1
    x2scaled = (x2-min_x2)/(max_x2-min_x2)*numbin_x2 + 1
    x3scaled = (x3-min_x3)/(max_x3-min_x3)*numbin_x3 + 1
  
    x1grid0 = int(x1scaled - order/2.0) + 1
    x1grid1 = int(x1scaled + order/2.0)
  
    x2grid0 = int(x2scaled - order/2.0) + 1
    x2grid1 = int(x2scaled + order/2.0)
  
    x3grid0 = int(x3scaled - order/2.0) + 1
    x3grid1 = int(x3scaled + order/2.0)

    if(x1grid1 < 1 .or. x1grid0 > num_x1 .or. x2grid1 < 1 .or. x2grid0 > num_x2 &
      .or. x3grid1 < 1 .or. x3grid0 > num_x3) then
      write(*,*) "out of range (bspline_3d)" 
      stop
    endif
  
    if(x1grid1 > num_x1) x1grid1 = num_x1
    if(x1grid0 < 1) x1grid0 = 1
  
    if(x2grid1 > num_x2) x2grid1 = num_x2
    if(x2grid0 < 1) x2grid0 = 1

    if(x3grid1 > num_x3) x3grid1 = num_x3
    if(x3grid0 < 1) x3grid0 = 1

    x1scaled = x1scaled + order/2.0
    x2scaled = x2scaled + order/2.0 
    x3scaled = x3scaled + order/2.0 
  
    value = 0.0
    d1 = 0.0
    d2 = 0.0
    d3 = 0.0

    !$omp parallel do collapse(3) reduction(+:value,d1,d2,d3) &
    !$omp private(datapoint, Mn1a, Mn1b, Mn2a, Mn2b, Mn3a, Mn3b, &
    !$omp dMn1, dMn2, dMn3)
    do i = x1grid0, x1grid1
      do j = x2grid0, x2grid1
        do k = x3grid0, x3grid1
          datapoint = data(k, j, i)
    
          ! breaking them down here since they are reusable for the first-derivative if requested
          Mn1a = Mn(order-1, x1scaled - i)
          Mn1b = Mn(order-1, x1scaled - i - 1)
          Mn2a = Mn(order-1, x2scaled - j)
          Mn2b = Mn(order-1, x2scaled - j - 1)
          Mn3a = Mn(order-1, x3scaled - k)
          Mn3b = Mn(order-1, x3scaled - k - 1)
  
          Mn1 = (x1scaled - i)/(order-1)*Mn1a + (order-(x1scaled - i))/(order-1)*Mn1b
          Mn2 = (x2scaled - j)/(order-1)*Mn2a + (order-(x2scaled - j))/(order-1)*Mn2b
          Mn3 = (x3scaled - k)/(order-1)*Mn3a + (order-(x3scaled - k))/(order-1)*Mn3b
    
          value = value + datapoint*Mn1*Mn2*Mn3
    
          if(f_d) then
            dMn1 = Mn1a - Mn1b
            dMn2 = Mn2a - Mn2b
            dMn3 = Mn3a - Mn3b

            d1 = d1 + datapoint*dMn1*Mn2*Mn3 & 
              *(numbin_x1/(max_x1-min_x1))
    
            d2 = d2 + datapoint*Mn1*dMn2*Mn3 &
              *(numbin_x2/(max_x2-min_x2))

            d3 = d3 + datapoint*Mn1*Mn2*dMn3 &
              *(numbin_x3/(max_x3-min_x3))
          endif
        enddo 
      enddo
    enddo
    !omp end parallel do

    r_value = value
    r_d1 = d1
    r_d2 = d2
    r_d3 = d3
  
  end subroutine bspline_3d

  subroutine bspline_4d(data, &
    min_x1, max_x1, num_x1, x1, &
    min_x2, max_x2, num_x2, x2, &
    min_x3, max_x3, num_x3, x3, &
    min_x4, max_x4, num_x4, x4, &
    order, r_value, r_d1, r_d2, r_d3, r_d4, f_d)
  
    implicit none
  
    integer :: num_x1, num_x2, num_x3, num_x4
    real(kind = 8), dimension(num_x4, num_x3, num_x2, num_x1) :: data
    real(kind = 8) :: min_x1, max_x1, x1, &
      min_x2, max_x2, x2, &
      min_x3, max_x3, x3, &
      min_x4, max_x4, x4
    integer :: order
    logical :: f_d
    real(kind = 8), intent(out) :: r_value, r_d1, r_d2, r_d3, r_d4
  
    integer :: numbin_x1, numbin_x2, numbin_x3, numbin_x4, i, j, k, m, &
      x1grid0, x1grid1, x2grid0, x2grid1, x3grid0, x3grid1, x4grid0, x4grid1
    real(kind = 8) :: x1scaled, x2scaled, x3scaled, x4scaled, datapoint
    real(kind = 8) :: value, d1, d2, d3, d4, &
      Mn1, Mn2, Mn3, Mn4, dMn1, dMn2, dMn3, dMn4
    real(kind = 8) :: Mn1a, Mn1b, Mn2a, Mn2b, Mn3a, Mn3b, Mn4a, Mn4b

    numbin_x1 = num_x1 - 1;
    numbin_x2 = num_x2 - 1;
    numbin_x3 = num_x3 - 1;
    numbin_x4 = num_x4 - 1;

    x1scaled = (x1-min_x1)/(max_x1-min_x1)*numbin_x1 + 1
    x2scaled = (x2-min_x2)/(max_x2-min_x2)*numbin_x2 + 1
    x3scaled = (x3-min_x3)/(max_x3-min_x3)*numbin_x3 + 1
    x4scaled = (x4-min_x4)/(max_x4-min_x4)*numbin_x4 + 1

    x1grid0 = int(x1scaled - order/2.0) + 1
    x1grid1 = int(x1scaled + order/2.0)
  
    x2grid0 = int(x2scaled - order/2.0) + 1
    x2grid1 = int(x2scaled + order/2.0)
  
    x3grid0 = int(x3scaled - order/2.0) + 1
    x3grid1 = int(x3scaled + order/2.0)

    x4grid0 = int(x4scaled - order/2.0) + 1
    x4grid1 = int(x4scaled + order/2.0)

    if(x1grid1 < 1 .or. x1grid0 > num_x1 .or. x2grid1 < 1 .or. x2grid0 > num_x2 &
      .or. x3grid1 < 1 .or. x3grid0 > num_x3 .or. x4grid1 < 1 .or. x4grid0 > num_x4) then
      write(*,*) "out of range (bspline_4d)" 
      stop
    endif
  
    if(x1grid1 > num_x1) x1grid1 = num_x1
    if(x1grid0 < 1) x1grid0 = 1
  
    if(x2grid1 > num_x2) x2grid1 = num_x2
    if(x2grid0 < 1) x2grid0 = 1

    if(x3grid1 > num_x3) x3grid1 = num_x3
    if(x3grid0 < 1) x3grid0 = 1

    if(x4grid1 > num_x4) x4grid1 = num_x4
    if(x4grid0 < 1) x4grid0 = 1

    x1scaled = x1scaled + order/2.0
    x2scaled = x2scaled + order/2.0 
    x3scaled = x3scaled + order/2.0 
    x4scaled = x4scaled + order/2.0 

    value = 0.0
    d1 = 0.0
    d2 = 0.0
    d3 = 0.0
    d4 = 0.0

    !$omp parallel do collapse(4) reduction(+:value,d1,d2,d3,d4) &
    !$omp private(datapoint, Mn1a, Mn1b, Mn2a, Mn2b, Mn3a, Mn3b, Mn4a, Mn4b, &
    !$omp Mn1, Mn2, Mn3, Mn4, dMn1, dMn2, dMn3, dMn4)
    do i = x1grid0, x1grid1
      do j = x2grid0, x2grid1
        do k = x3grid0, x3grid1
          do m = x4grid0, x4grid1
            datapoint = data(m, k, j, i)
      
            Mn1a = Mn(order-1, x1scaled - i)
            Mn1b = Mn(order-1, x1scaled - i - 1)
            Mn2a = Mn(order-1, x2scaled - j)
            Mn2b = Mn(order-1, x2scaled - j - 1)
            Mn3a = Mn(order-1, x3scaled - k)
            Mn3b = Mn(order-1, x3scaled - k - 1)
            Mn4a = Mn(order-1, x4scaled - m)
            Mn4b = Mn(order-1, x4scaled - m - 1)
            
            Mn1 = (x1scaled - i)/(order-1)*Mn1a + (order-(x1scaled - i))/(order-1)*Mn1b
            Mn2 = (x2scaled - j)/(order-1)*Mn2a + (order-(x2scaled - j))/(order-1)*Mn2b
            Mn3 = (x3scaled - k)/(order-1)*Mn3a + (order-(x3scaled - k))/(order-1)*Mn3b
            Mn4 = (x4scaled - m)/(order-1)*Mn4a + (order-(x4scaled - m))/(order-1)*Mn4b

            value = value + datapoint*Mn1*Mn2*Mn3*Mn4
      
            if(f_d) then
              dMn1 = Mn1a - Mn1b
              dMn2 = Mn2a - Mn2b
              dMn3 = Mn3a - Mn3b
              dMn4 = Mn4a - Mn4b
  
              d1 = d1 + datapoint*dMn1*Mn2*Mn3*Mn4 & 
                *(numbin_x1/(max_x1-min_x1))
      
              d2 = d2 + datapoint*Mn1*dMn2*Mn3*Mn4 &
                *(numbin_x2/(max_x2-min_x2))
  
              d3 = d3 + datapoint*Mn1*Mn2*dMn3*Mn4 &
                *(numbin_x3/(max_x3-min_x3))

              d4 = d4 + datapoint*Mn1*Mn2*Mn3*dMn4 &
                *(numbin_x4/(max_x4-min_x4))                
            endif
          enddo
        enddo 
      enddo
    enddo
    !omp end parallel do

    r_value = value
    r_d1 = d1
    r_d2 = d2
    r_d3 = d3
    r_d4 = d4
  
  end subroutine bspline_4d

  subroutine bspline_5d(data, &
    min_x1, max_x1, num_x1, x1, &
    min_x2, max_x2, num_x2, x2, &
    min_x3, max_x3, num_x3, x3, &
    min_x4, max_x4, num_x4, x4, &
    min_x5, max_x5, num_x5, x5, &
    order, r_value, r_d1, r_d2, r_d3, r_d4, r_d5, f_d)
  
    implicit none
  
    integer :: num_x1, num_x2, num_x3, num_x4, num_x5
    real(kind = 8), dimension(num_x5, num_x4, num_x3, num_x2, num_x1) :: data
    real(kind = 8) :: min_x1, max_x1, x1, &
      min_x2, max_x2, x2, &
      min_x3, max_x3, x3, &
      min_x4, max_x4, x4, &
      min_x5, max_x5, x5
    integer :: order
    logical :: f_d
    real(kind = 8), intent(out) :: r_value, r_d1, r_d2, r_d3, r_d4, r_d5
  
    integer :: numbin_x1, numbin_x2, numbin_x3, numbin_x4, numbin_x5, i, j, k, m, n, &
      x1grid0, x1grid1, x2grid0, x2grid1, x3grid0, x3grid1, x4grid0, x4grid1, x5grid0, x5grid1
    real(kind = 8) :: x1scaled, x2scaled, x3scaled, x4scaled, x5scaled, datapoint
    real(kind = 8) :: value, d1, d2, d3, d4, d5, &
      Mn1, Mn2, Mn3, Mn4, Mn5, dMn1, dMn2, dMn3, dMn4, dMn5
    real(kind = 8) :: Mn1a, Mn1b, Mn2a, Mn2b, Mn3a, Mn3b, Mn4a, Mn4b, Mn5a, Mn5b

    numbin_x1 = num_x1 - 1;
    numbin_x2 = num_x2 - 1;
    numbin_x3 = num_x3 - 1;
    numbin_x4 = num_x4 - 1;
    numbin_x5 = num_x5 - 1;

    x1scaled = (x1-min_x1)/(max_x1-min_x1)*numbin_x1 + 1
    x2scaled = (x2-min_x2)/(max_x2-min_x2)*numbin_x2 + 1
    x3scaled = (x3-min_x3)/(max_x3-min_x3)*numbin_x3 + 1
    x4scaled = (x4-min_x4)/(max_x4-min_x4)*numbin_x4 + 1
    x5scaled = (x5-min_x5)/(max_x5-min_x5)*numbin_x5 + 1

    x1grid0 = int(x1scaled - order/2.0) + 1
    x1grid1 = int(x1scaled + order/2.0)
  
    x2grid0 = int(x2scaled - order/2.0) + 1
    x2grid1 = int(x2scaled + order/2.0)
  
    x3grid0 = int(x3scaled - order/2.0) + 1
    x3grid1 = int(x3scaled + order/2.0)

    x4grid0 = int(x4scaled - order/2.0) + 1
    x4grid1 = int(x4scaled + order/2.0)

    x5grid0 = int(x5scaled - order/2.0) + 1
    x5grid1 = int(x5scaled + order/2.0)

    if(x1grid1 < 1 .or. x1grid0 > num_x1 .or. x2grid1 < 1 .or. x2grid0 > num_x2 &
      .or. x3grid1 < 1 .or. x3grid0 > num_x3 .or. x4grid1 < 1 .or. x4grid0 > num_x4 &
      .or. x5grid1 < 1 .or. x5grid0 > num_x5) then
      write(*,*) "out of range (bspline_5d)" 
      stop
    endif
  
    if(x1grid1 > num_x1) x1grid1 = num_x1
    if(x1grid0 < 1) x1grid0 = 1
  
    if(x2grid1 > num_x2) x2grid1 = num_x2
    if(x2grid0 < 1) x2grid0 = 1

    if(x3grid1 > num_x3) x3grid1 = num_x3
    if(x3grid0 < 1) x3grid0 = 1

    if(x4grid1 > num_x4) x4grid1 = num_x4
    if(x4grid0 < 1) x4grid0 = 1

    if(x5grid1 > num_x5) x5grid1 = num_x5
    if(x5grid0 < 1) x5grid0 = 1

    x1scaled = x1scaled + order/2.0
    x2scaled = x2scaled + order/2.0 
    x3scaled = x3scaled + order/2.0 
    x4scaled = x4scaled + order/2.0 
    x5scaled = x5scaled + order/2.0 

    value = 0.0
    d1 = 0.0
    d2 = 0.0
    d3 = 0.0
    d4 = 0.0
    d5 = 0.0

    !$omp parallel do collapse(5) reduction(+:value,d1,d2,d3,d4,d5) &
    !$omp private(datapoint, Mn1a, Mn1b, Mn2a, Mn2b, Mn3a, Mn3b, Mn4a, Mn4b, Mn5a, Mn5b, &
    !$omp Mn1, Mn2, Mn3, Mn4, Mn5, dMn1, dMn2, dMn3, dMn4, dMn5)
    do i = x1grid0, x1grid1
      do j = x2grid0, x2grid1
        do k = x3grid0, x3grid1
          do m = x4grid0, x4grid1
            do n = x5grid0, x5grid1
              datapoint = data(n, m, k, j, i)
        
              Mn1a = Mn(order-1, x1scaled - i)
              Mn1b = Mn(order-1, x1scaled - i - 1)
              Mn2a = Mn(order-1, x2scaled - j)
              Mn2b = Mn(order-1, x2scaled - j - 1)
              Mn3a = Mn(order-1, x3scaled - k)
              Mn3b = Mn(order-1, x3scaled - k - 1)
              Mn4a = Mn(order-1, x4scaled - m)
              Mn4b = Mn(order-1, x4scaled - m - 1)
              Mn5a = Mn(order-1, x5scaled - n)
              Mn5b = Mn(order-1, x5scaled - n - 1)
              
              Mn1 = (x1scaled - i)/(order-1)*Mn1a + (order-(x1scaled - i))/(order-1)*Mn1b
              Mn2 = (x2scaled - j)/(order-1)*Mn2a + (order-(x2scaled - j))/(order-1)*Mn2b
              Mn3 = (x3scaled - k)/(order-1)*Mn3a + (order-(x3scaled - k))/(order-1)*Mn3b
              Mn4 = (x4scaled - m)/(order-1)*Mn4a + (order-(x4scaled - m))/(order-1)*Mn4b
              Mn5 = (x5scaled - n)/(order-1)*Mn5a + (order-(x5scaled - n))/(order-1)*Mn5b
  
              value = value + datapoint*Mn1*Mn2*Mn3*Mn4*Mn5
        
              if(f_d) then
                dMn1 = Mn1a - Mn1b
                dMn2 = Mn2a - Mn2b
                dMn3 = Mn3a - Mn3b
                dMn4 = Mn4a - Mn4b
                dMn5 = Mn5a - Mn5b
    
                d1 = d1 + datapoint*dMn1*Mn2*Mn3*Mn4*Mn5 & 
                  *(numbin_x1/(max_x1-min_x1))
        
                d2 = d2 + datapoint*Mn1*dMn2*Mn3*Mn4*Mn5 &
                  *(numbin_x2/(max_x2-min_x2))
    
                d3 = d3 + datapoint*Mn1*Mn2*dMn3*Mn4*Mn5 &
                  *(numbin_x3/(max_x3-min_x3))
  
                d4 = d4 + datapoint*Mn1*Mn2*Mn3*dMn4*Mn5 &
                  *(numbin_x4/(max_x4-min_x4))   
                  
                d5 = d5 + datapoint*Mn1*Mn2*Mn3*Mn4*dMn5 &
                  *(numbin_x5/(max_x5-min_x5))    
              endif
            enddo
          enddo
        enddo 
      enddo
    enddo
    !omp end parallel do

    r_value = value
    r_d1 = d1
    r_d2 = d2
    r_d3 = d3
    r_d4 = d4
    r_d5 = d5

  end subroutine bspline_5d

end module bspline