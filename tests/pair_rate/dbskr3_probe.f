      program dbskr3_probe
      implicit none

      integer n, i
      parameter (n = 600)
      double precision chi_min, chi_max, log_chi, chi, x, y
      double precision dbskr3
      external dbskr3

      chi_min = 1.d-4
      chi_max = 1.d1

      print *, 'chi,x,dbskr3_x_1'
      do i = 0, n - 1
        log_chi = log10(chi_min)
     +            + (log10(chi_max) - log10(chi_min)) * i / (n - 1)
        chi = 10.d0 ** log_chi
        x = 2.d0 / (3.d0 * chi)
        y = dbskr3(x, 1)
        write (*,'(1pe26.16e3,",",1pe26.16e3,",",1pe26.16e3)')
     +       chi, x, y
      enddo

      end
