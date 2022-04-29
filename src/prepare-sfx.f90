program main
   use sfx
   implicit none

   integer,parameter :: n_sim = 1000
   type(crystal) cryst(n_sim)
   type(panel) geom(64)
   real(8) xs, ys
   integer n, m

   call read_panels('cspad_lysozyme_cxidb_refined.geom',geom)

   open(20,file='lysozyme-mosflm-1000.stream',status='old')
   open(30,file="SPOTS.TXT",status='replace')
   open(40,file="wavelengths.lst",status='replace')
   do n = 1, n_sim
      call read_one_image(20,cryst(n))
      write(40,*) cryst(n)%lambda
      do m = 1, cryst(n)%num_peaks
         call panel2cart(cryst(n)%obs(m)%fs,cryst(n)%obs(m)%ss,cryst(n),geom,xs,ys)
         write(30,'(I5,2F9.2)') n, xs, ys
      end do
   end do
   close(20)
   close(30)
   close(40)

   stop

end program
