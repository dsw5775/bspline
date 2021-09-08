program testbspline
! testing bspline_1d to bspline_5d with order 6
! you can change order with command line parameter, e.g. ./bspline 4

use bspline

implicit none

integer :: num_x1, num_x2, num_x3, num_x4, num_x5
real(kind = 8), allocatable, dimension(:) :: data1
real(kind = 8), allocatable, dimension(:,:) :: data2
real(kind = 8), allocatable, dimension(:,:,:) :: data3
real(kind = 8), allocatable, dimension(:,:,:,:) :: data4
real(kind = 8), allocatable, dimension(:,:,:,:,:) :: data5
real(kind = 8) :: min_x1, max_x1, x1, &
  min_x2, max_x2, x2, &
  min_x3, max_x3, x3, &
  min_x4, max_x4, x4, &
  min_x5, max_x5, x5
integer :: order

integer :: numbin_x1, numbin_x2, numbin_x3, numbin_x4, numbin_x5, i, j, k, m, n
real(kind = 8) :: width_x1, width_x2, width_x3, width_x4, width_x5
real(kind = 8) :: value, d1, d2, d3, d4, d5

real :: start, end
character(len = 2) :: arg
order = 6
call get_command_argument(1, arg)
if(len_trim(arg) /= 0) then
  read (arg,'(I2)') order
endif

num_x1 = 100+1
min_x1 = 0.0d0
max_x1 = 1.0d0
width_x1 = (max_x1-min_x1) / (num_x1-1)
x1 = 0.5d0

! test bspline_1d for data generated from function x^2
allocate(data1(num_x1))
do i = 1, num_x1
  data1(i) = (width_x1*(i-1)) * (width_x1*(i-1))
enddo

call cpu_time(start)
call bspline_1d(data1, min_x1, max_x1, &
  num_x1, x1, order, value, &
  d1, .true., d2, .true.)
call cpu_time(end)

write(*,*) "bspline_1d test:"
write(*,*) "Calculated:"
write(*,'(3F20.6)') value, d1, d2
write(*,*) "Expected:"
write(*,'(3F20.6)') x1*x1, 2*x1, 2.0
write(*,'(A6F10.6A/)') "Time:", end - start, " sec"

! test bspline_2d for data generated from function x1^2*x2
num_x1 = 100
min_x1 = 0.005d0
max_x1 = 0.995d0
width_x1 = (max_x1-min_x1) / (num_x1-1)
x1 = 0.3d0
num_x2 = 24+1
min_x2 = -10000.0d0
max_x2 = 10000.0d0
width_x2 = (max_x2-min_x2) / (num_x2-1)
x2 = 4000.0d0

allocate(data2(num_x2, num_x1))
do i = 1, num_x1
  do j = 1, num_x2
    data2(j, i) = (min_x1+width_x1*(i-1) ) * (min_x1+width_x1*(i-1) ) * &
    (min_x2+width_x2*(j-1) )
  enddo
enddo

call cpu_time(start)
call bspline_2d(data2, min_x1, max_x1, num_x1, x1, &
  min_x2, max_x2, num_x2, x2, &
  order, value, d1, d2, .true.)
call cpu_time(end)

write(*,*) "bspline_2d test:"
write(*,*) "Calculated:"
write(*,'(3F20.6)') value, d1, d2
write(*,*) "Expected:"
write(*,'(3F20.6)') x1*x1*x2, 2*x1*x2, x1*x1
write(*,'(A6F10.6A/)') "Time:", end - start, " sec"

! test bspline_3d for data generated from function x1^2*x2*x3
num_x3 = 24+1
min_x3 = -0.01d0
max_x3 = 0.01d0
width_x3 = (max_x3-min_x3) / (num_x3-1)
x3 = 0.005d0

allocate(data3(num_x3, num_x2, num_x1))
do i = 1, num_x1
  do j = 1, num_x2
    do k = 1, num_x3
      data3(k, j, i) = (min_x1+width_x1*(i-1) ) * (min_x1+width_x1*(i-1) ) * &
      (min_x2+width_x2*(j-1) ) * (min_x3+width_x3*(k-1) )
    enddo
  enddo
enddo

call cpu_time(start)
call bspline_3d(data3, min_x1, max_x1, num_x1, x1, &
  min_x2, max_x2, num_x2, x2, &
  min_x3, max_x3, num_x3, x3, &
  order, value, d1, d2, d3, .true.)
call cpu_time(end)

write(*,*) "bspline_3d test:"
write(*,*) "Calculated:"
write(*,'(4F20.6)') value, d1, d2, d3
write(*,*) "Expected:"
write(*,'(4F20.6)') x1*x1*x2*x3, 2*x1*x2*x3, x1*x1*x3, x1*x1*x2
write(*,'(A6F10.6A/)') "Time:", end - start, " sec"

! test bspline_4d for data generated from function x1^2*x2*x3*x4
num_x4 = 24+1
min_x4 = -10000.0d0
max_x4 = 10000.0d0
width_x4 = (max_x4-min_x4) / (num_x4-1)
x4 = 6000.0d0

allocate(data4(num_x4, num_x3, num_x2, num_x1))
do i = 1, num_x1
  do j = 1, num_x2
    do k = 1, num_x3
      do m = 1, num_x4
        data4(m, k, j, i) = (min_x1+width_x1*(i-1) ) * (min_x1+width_x1*(i-1) ) * &
        (min_x2+width_x2*(j-1) ) * (min_x3+width_x3*(k-1) ) * (min_x4+width_x4*(m-1) )
      enddo
    enddo
  enddo
enddo

call cpu_time(start)
call bspline_4d(data4, min_x1, max_x1, num_x1, x1, &
  min_x2, max_x2, num_x2, x2, &
  min_x3, max_x3, num_x3, x3, &
  min_x4, max_x4, num_x4, x4, &
  order, value, d1, d2, d3, d4, .true.)
call cpu_time(end)

write(*,*) "bspline_4d test:"
write(*,*) "Calculated:"
write(*,'(5F20.6)') value, d1, d2, d3, d4
write(*,*) "Expected:"
write(*,'(5F20.6)') x1*x1*x2*x3*x4, 2*x1*x2*x3*x4, x1*x1*x3*x4, x1*x1*x2*x4, x1*x1*x2*x3
write(*,'(A6F10.6A/)') "Time:", end - start, " sec"

! test bspline_5d for data generated from function x1^2*x2*x3*x4*x5
num_x5 = 24+1
min_x5 = -0.01d0
max_x5 = 0.01d0
width_x5 = (max_x5-min_x5) / (num_x5-1)
x5 = 0.007d0

allocate(data5(num_x5, num_x4, num_x3, num_x2, num_x1))
do i = 1, num_x1
  do j = 1, num_x2
    do k = 1, num_x3
      do m = 1, num_x4
        do n = 1, num_x5
          data5(n, m, k, j, i) = (min_x1+width_x1*(i-1) ) * (min_x1+width_x1*(i-1) ) * &
          (min_x2+width_x2*(j-1) ) * (min_x3+width_x3*(k-1) ) * (min_x4+width_x4*(m-1) ) &
          * (min_x5+width_x5*(n-1) )
        enddo
      enddo
    enddo
  enddo
enddo

call cpu_time(start)
call bspline_5d(data5, min_x1, max_x1, num_x1, x1, &
  min_x2, max_x2, num_x2, x2, &
  min_x3, max_x3, num_x3, x3, &
  min_x4, max_x4, num_x4, x4, &
  min_x5, max_x5, num_x5, x5, &
  order, value, d1, d2, d3, d4, d5, .true.)
call cpu_time(end)

write(*,*) "bspline_5d test:"
write(*,*) "Calculated:"
write(*,'(6F20.6)') value, d1, d2, d3, d4, d5
write(*,*) "Expected:"
write(*,'(6F20.6)') x1*x1*x2*x3*x4*x5, 2*x1*x2*x3*x4*x5, x1*x1*x3*x4*x5, &
  x1*x1*x2*x4*x5, x1*x1*x2*x3*x5, x1*x1*x2*x3*x4
write(*,'(A6F10.6A/)') "Time:", end - start, " sec"

deallocate(data1, data2, data3, data4)

end program testbspline