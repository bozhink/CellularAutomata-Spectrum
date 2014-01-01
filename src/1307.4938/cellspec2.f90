program cellular_automata_spectrum
   use cellular_automata
   use spectrum
   implicit none
   integer, parameter :: T = 1024, N = 1000, bits=8, TN=T*N
   integer, parameter :: number_of_iterations = 30
   double complex, dimension (T, N) :: Sf
   double precision, dimension (T) :: S
   integer, dimension (T, N) :: A
   integer, dimension (bits) :: code
   integer :: i, j, iclock, tic, tic1, toc
   integer :: code_number, codenumber
   character(len=33) :: arg
   data A /TN*0/
   !data code /0,0,0,1,1,1,1,0/ ! rule 30
   !data code /0,1,1,0,1,1,1,0/ ! rule 110
   data code /0,0,1,1,1,0,1,1/
   
   do codenumber = 251, 251
      print*,codenumber
      call number_code (codenumber, code)
      print*,code_number (code)
   
      write (unit=arg, fmt="('ca_code_',i3.3,'_t_',i6.6,'_n_',i6.6,'.dat')") codenumber, T, N
      open (unit=10, file=arg, status='unknown', action='write', position='append', iostat=i)
      if (i/=0) then
         stop 'Error opening output file!'
      endif
      ! Initialization
      call system_clock(iclock)
      call srand(iclock)
      tic = iclock
      toc = tic
      write (unit=10, fmt="('# Cellular automata spectrum')")
      write (unit=10, fmt="('# Number of lattice positions (chain length):',1x,i5)") N
      write (unit=10, fmt="('# Length in time:',1x,i5)") T
      write (unit=10, fmt="('# Number of generated spectra:',1x,i5)") number_of_iterations
      write (unit=10, fmt="('# Code number:',1x,i6)") codenumber
      write (unit=10, fmt="('#')")
      
      call spectrum_init(T)
      ! Iterations
      do j = 1, number_of_iterations
         tic1 = toc
         write (unit=10, fmt="(30('#'))")
         write (unit=10, fmt="('# Iteration #',i2)") j
         write (unit=10, fmt="('# channel   amplitude')")
         write (unit=10, fmt="('#',1x,28('='))")
         do i = 1, N
            A(1,i) = int(2*rand())
         end do
         call cellauth (T, N, A, code)
         call fourier2 (T, N, A, Sf)
         call spec(T, N, Sf, S)
         do i = 1, N
            write (unit=10, fmt="(i9,1x,e18.10)") i-1, S(i)
         end do
         call system_clock(toc)
         print "('Time of iteration #',i3,':',1x,f10.3,'s.')", j, (toc-tic1)/1000.0
      end do
      close(unit=10)
      call system_clock(toc)
      print "(36('='))"
      print "('Time:',1x,f10.3,'s.')", (toc-tic)/1000.0
   end do
   call spectrum_finalize()
end program cellular_automata_spectrum

function code_number(code) result (res)
   implicit none
   integer, dimension(8), intent(in) :: code
   integer :: res
   integer :: i
   res = 2*code(1)
   do i = 2, 7
      res = 2*(res + code(i))
   end do
   res = res + code(8)
end function code_number

subroutine number_code(number, code)
   implicit none
   integer, intent(in) :: number
   integer, dimension(8), intent(out) :: code
   integer :: idivv, imodd, i
   code = (/0,0,0,0,0,0,0,0/)
   idivv = number
   do i = 8,1,-1
      imodd = idivv
      idivv = idivv / 2
      code(i) = imodd - 2*idivv
   end do
end subroutine number_code
