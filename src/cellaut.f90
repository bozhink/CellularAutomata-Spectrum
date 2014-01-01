! file: cellaut.f90
! author: Bozhin Karaivanov
! date: 28.07.2013
!
! This simulation uses periodic boundaries conditions

module cellular_automata

   public :: cellaut
   public :: decode
   
contains

   subroutine cellauth (T, N, A, code)
   ! T = length in time
   ! N = size of the matrix
   ! A = CA output matrix
   ! code = vector of size 8, containing 0 and 1
      implicit none
      integer, parameter :: bits = 8
      integer, intent (in) :: T, N
      integer, dimension(bits), intent(in) :: code
      integer, dimension(T, N), intent(in out) :: A
      integer :: i, j, i1
      
      do i = 2, T
         i1 = i-1
         A(i,1) = decode(code, A(i1,N), A(i1,1), A(i1,2))
         do j = 2, N-1
            A(i,j) = decode(code, A(i1,j-1), A(i1,j), A(i1,j+1))
         enddo
         A(i,N) = decode(code, A(i1,N-1), A(i1,N), A(i1,1))
      enddo
   end subroutine cellauth
   
   
   function decode (code, cell0, cell1, cell2) result (res)
      implicit none
      integer, parameter :: bits=8
      integer, dimension(bits), intent(in) :: code
      integer, intent(in) :: cell0, cell1, cell2
      integer :: res
      select case (cell0)
         case (0)
            select case (cell1)
               case (0)
                  select case (cell2)
                     case (0)
                        res = code(8)
                     case (1)
                        res = code(7)
                     case default
                        res = -1
                  end select
               case (1)
                  select case (cell2)
                     case (0)
                        res = code(6)
                     case (1)
                        res = code(5)
                     case default
                        res = -1
                  end select
               case default
                  res = -1
            end select
         case (1)
            select case (cell1)
               case (0)
                  select case (cell2)
                     case (0)
                        res = code(4)
                     case (1)
                        res = code(3)
                     case default
                        res = -1
                  end select
               case (1)
                  select case (cell2)
                     case (0)
                        res = code(2)
                     case (1)
                        res = code(1)
                     case default
                        res = -1
                  end select
               case default
                  res = -1
            end select
         case default
            res = -1
      end select
   end function decode


end module cellular_automata
