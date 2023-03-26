!=====================================================================================
! (c) IHEP 2021. All rights reserved.
! Author: Geng Zhi
! Institute of High Energy Physics, Chinese Academy of Science (IHEP, CAS). 
! If you have any problem with the program, please contact author with the following 
! email: gengz@ihep.ac.cn
!======================================================================================

program main
  use indexing
  use mpi
  implicit none

  character(200) line
  character(20) val, key

  logical enable_rescut, enable_uc_scan, enable_gridscan, enable_refcut

  integer is(3), js(3), flag, eof, nrefl, nrefl_low, ncentx, ncenty, cent_step
  integer num, interval, m, n, id, tmp, split_nz(2), nz_start
  integer nx, ny, nz, i, npsi, nphi, ntheta, n1, n2, n3, n4, n5, niter
  integer node, np, ierr, nfound, refcut, status

  real(8) :: om(3,3), rescut, res, tmp1, tmp2, cell_range, cell_step, local_angle_range, local_angle_step
  real(8) :: x(10000), y(10000), sub_x(10000), sub_y(10000), tmp_x(10000),tmp_y(10000), score1, score2, score_best
  real(8) :: delpsi, delphi, deltheta, psi, phi, theta
  real(8) :: qx, qy, qz, h, k, l, xo, yo
  real(8) :: rot(3,3), UB_B(3,3), om_best(3,3)
  real(8) :: best(5), psi_best, phi_best, theta_best, psi_init, phi_init, theta_init
  real(8) :: centx_best, centy_best, distance_best, ua_best, ub_best, uc_best
  real(8), allocatable :: best_all(:,:)
  real time1, time2


  type(fgsl_multimin_fminimizer) :: min_fslv
  type(fgsl_multimin_function) :: func
  type(data),target :: datas
  type(fgsl_vector) :: x_in, step_size
  real(8), allocatable, target :: xv(:), stepv(:)
  real(8), pointer :: xptr(:)
  type(c_ptr) :: ptr

  time1 = MPI_Wtime()

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)


  allocate(best_all(5, np))

  !read in configure file
  open(30, file='param.config', status='old')

  do while(.true.)

      read(30,'(a)', iostat=eof) line
      if(eof .lt. 0) exit
      if(line(1:1).eq.'#' .or. len_trim(adjustl(line)).eq.0) cycle

      call extract_values(line, key, val)

      if(key .eq. 'Space_group') then
          read(val, *) spg_num
      else if(key .eq. 'Crystal_system') then
          crystal_system = val
      else if(key .eq. 'Unit_cell_a') then
          read(val, *) uc(1)
      else if(key .eq. 'Unit_cell_b') then
          read(val, *) uc(2)
      else if(key .eq. 'Unit_cell_c') then
          read(val, *) uc(3)
      else if(key .eq. 'Unit_cell_alpha') then
          read(val, *) uc(4)
      else if(key .eq. 'Unit_cell_beta') then
          read(val, *) uc(5)
      else if(key .eq. 'Unit_cell_gama') then
          read(val, *) uc(6)
      else if(key .eq. 'Wavelength') then
          read(val, *) lambda
      else if(key .eq. 'Distance') then
          read(val, *) distance
      else if(key .eq. 'Pixel_size') then
          read(val, *) pixsize
      else if(key .eq. 'XCenter') then
          read(val, *) centx
      else if(key .eq. 'YCenter') then
          read(val, *) centy
      else if(key .eq. 'Enable_rescut') then
          read(val, *) tmp
          if(tmp .eq. 0) then
              enable_rescut = .false.
          else
              enable_rescut = .true.
          end if
      else if(key .eq. 'Rescut') then
          read(val, *) rescut
      else if(key .eq. 'Enable_refcut') then
          read(val, *) tmp
          if(tmp .eq. 0) then
              enable_refcut = .false.
          else
              enable_refcut = .true.
          end if
      else if(key .eq. 'Refcut') then
          read(val, *) refcut
      else if(key .eq. 'Position_error') then
          read(val, *) pos_err
      else if(key .eq. 'Enable_uc_scan') then
          read(val, *) tmp
          if(tmp .eq. 0) then
              enable_uc_scan = .false.
          else
              enable_uc_scan = .true.
          end if
      else if(key .eq. 'Enable_gridscan') then
          read(val, *) tmp
          if(tmp .eq. 0) then
              enable_gridscan = .false.
          else
              enable_gridscan = .true.
          end if
      else if(key .eq. 'XCScan') then
          read(val, *) ncentx
      else if(key .eq. 'YCScan') then
          read(val, *) ncenty
      else if(key .eq. 'XYCScanstep') then
          read(val, *) cent_step
      else if(key .eq. 'CellScan') then
          read(val, *) cell_range
      else if(key .eq. 'CellScanstep') then
          read(val, *) cell_step
      else if(key .eq. 'Local_angle_range') then
          read(val, *) local_angle_range
      else if(key .eq. 'Local_angle_step') then
          read(val, *) local_angle_step
      end if

  end do

  close(30)

  !master processor: read in spot within predefined resolution range
  if(node .eq. 0)then

     write(*,*)
     write(*,'(a)') "--------------------------------------------------------------"
     write(*,"(a)") "                 User Input Information"
     write(*,*)
     write(*,"(a,I4)") " Space group code: ", spg_num
     write(*,"(a,a)") " Crystal system: ", crystal_system
     write(*,"(a,3F10.4,3F8.2)") " Unit Cell: ", uc
     write(*,"(a,F8.6)") " Wavelength(A): ", lambda
     write(*,"(a,F6.2)") " Distance(mm): ", distance
     write(*,"(a,F6.4)") " Pixel size(mm): ", pixsize
     write(*,"(a,2F8.2)") " Beam center(pixels): ", centx, centy
     write(*,"(a,F4.2)") " Position Error for indexing(pixels): ", pos_err
     write(*,*) "Resolution cut? {T/F}: ", enable_rescut
     if(enable_rescut)then
         write(*,"(a,F4.1)") " Resolution cutoff(A): ", rescut
     end if
     write(*,*) "Reflection number cut? {T/F}: ", enable_refcut
     if(enable_refcut)then
         write(*,"(a,I4)") " Reflection number cutoff: ", refcut
     end if
     write(*,*) "Unit cell refine? {T/F}: ", enable_uc_scan
     if(enable_uc_scan)then
         write(*,"(a,2F4.1)") " Unit cell refinement range & step: ", cell_range, cell_step
     end if
     write(*,*) "Local refine? {T/F}: ", enable_gridscan
     if(enable_gridscan)then
         write(*,"(a,3I2)") " Beam center refinement range & step: ", ncentx, ncenty, cent_step
         write(*,"(a,2F4.1)") " Local angle range & step: ", local_angle_range,local_angle_step
     end if
     write(*,'(a)') "--------------------------------------------------------------"
     write(*,*)

     !read in spot position
     open(20,file="SPOT-TMP.TXT", status='old')
     nrefl = 0
     do while(.true.)
         read(20,*, iostat=eof) tmp1, tmp2
         if(eof .lt. 0) exit

         if(enable_rescut)then
             call resolution(tmp1, tmp2, dble(centx), dble(centy), dble(distance), res)
             if(res .le. rescut) cycle
         end if

         nrefl = nrefl + 1
         tmp_x(nrefl) = tmp1
         tmp_y(nrefl) = tmp2
     end do

     close(20)

     if(enable_refcut)then

        rescut = 1.0

        do while(nrefl > refcut)

            m = 0
            rescut = rescut + 0.1
            do n = 1, nrefl
               call resolution(tmp_x(n),tmp_y(n), dble(centx), dble(centy), dble(distance), res)
               if(res .le. rescut) cycle
               m = m + 1
               tmp_x(m) = tmp_x(n)
               tmp_y(m) = tmp_y(n)
            end do

            nrefl = m

        end do

     end if


     nrefl_low = 0
     do n = 1, nrefl
        x(n) = tmp_x(n)
        y(n) = tmp_y(n)
        call resolution(x(n), y(n), dble(centx), dble(centy), dble(distance), res)
        if(res .ge. 4.0) nrefl_low = nrefl_low + 1
     end do

     !calculate crystallographic transformation matrix from unitcell
     call transformation_matrix(uc, UB_B)

  end if

  !broadcast to other slave processors from master
  call MPI_Bcast(nrefl, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(nrefl_low, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(x, 10000, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(y, 10000, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(UB_B, 9, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

  if(nrefl .lt. 5 .or. nrefl_low .lt. 3)then
     call MPI_FINALIZE(ierr)
     stop
  end if

  !initial orientation search step and number
  delpsi = 0.03; npsi = nint(pi/2/delpsi)
  deltheta = 0.03; ntheta = nint(pi*2/deltheta)

  !master processor: allocate different jobs to slave processors uniformly
  if(node .eq. 0)then

     num = 0
     do nz = 1, npsi
        nphi = nint(2*pi*sin(nz*delpsi)/delpsi)
        do ny = 1, nphi
           do nx = 1, ntheta
              num = num + 1
           end do
        end do
     end do

     interval = nint(num*1.0/(np-1))

     n = 0; nz_start = 1
     outer:do nz = 1, npsi
        nphi = nint(2*pi*sin(nz*delpsi)/delpsi)
        do ny = 1, nphi
           do nx = 1, ntheta
              n = n + 1
              tmp = mod(n, interval)
              if(tmp .eq. 0 )then
                 id = n/interval
                 split_nz = (/nz_start, nz/)
                 nz_start = nz + 1
                 call MPI_SEND(split_nz, 2, MPI_INTEGER, id, 99, MPI_COMM_WORLD, ierr)
                 if(id .eq. np-2)then
                    split_nz = (/nz_start, npsi/)
                    exit outer
                 end if
              end if
           end do
        end do
     end do outer
 
     call MPI_SEND(split_nz, 2, MPI_INTEGER, np-1, 99, MPI_COMM_WORLD, ierr)

  end if

  !slave processors: perform initial orientation search
  if(node .ne. 0)then

     call MPI_RECV(split_nz, 2, MPI_INTEGER, 0, 99, MPI_COMM_WORLD, ierr)

!     write(*,*) node, split_nz

     score_best = 0.01
     do nz = split_nz(1), split_nz(2)
         nphi = nint(2*pi*sin(nz*delpsi)/delpsi)
         delphi = 2*pi/nphi
         psi = nz*delpsi
         do ny = 1, nphi
            phi = ny*delphi
!            if(phi .lt. pi/3*3 .or. phi .gt. pi/3*4) cycle
            do nx = 1, ntheta
               theta = nx*deltheta
               call rotation_mat_from_rodrigues(psi, phi, theta, rot)
               om = matmul(rot,UB_B)
               call evaluate(om, x, y, dble(centx), dble(centy), dble(distance), nrefl, score1, score2)
!               write(*,*) psi, phi,theta,score2
               if(score2 .gt. score_best)then
                   score_best = score2
                   best = (/score_best, psi, phi, theta, score1/)
               end if
            end do
         end do
     end do

!     write(*,*) node, score_best

  end if

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  call MPI_Gather(best, 5, MPI_REAL8, best_all, 5, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

  if(node .eq. 0)then

      !find the best solution
      score_best = 0.01
      do i = 2, np

          if(best_all(1,i) .gt. score_best)then
              score_best = best_all(1,i)
              psi_init = best_all(2,i)
              phi_init = best_all(3,i)
              theta_init = best_all(4,i)
              score1 = best_all(5,i)
              best = (/score_best, psi_init, phi_init, theta_init, score1/)
          end if

      end do

      write(*,*)
      write(*,'(a)') "--------------------------------------------------------------"
      write(*,'(a)') "                 Initial indexing result"
      write(*,*)
      write(*,'(a,F6.3,a,F8.3)') " Match rate: ", best(1),  "  Positional deviation (pixels): ", best(5)
      write(*,'(a)') "--------------------------------------------------------------"


      psi_best = psi_init; phi_best = phi_init; theta_best = theta_init
      ua_best = uc(1); ub_best = uc(2); uc_best = uc(3)
      centx_best = centx; centy_best = centy; distance_best = distance

      !prior information correction
      if(enable_gridscan) then

         !define neighborhood region about initial orientation
         delpsi = local_angle_step/180*pi; npsi = nint(pi*local_angle_range/180/delpsi)
         delphi = local_angle_step/180*pi; nphi = nint(pi*local_angle_range/180/delphi)
         deltheta = local_angle_step/180*pi; ntheta = nint(pi*local_angle_range/180/deltheta)

         if(.not. enable_uc_scan)then

            !only scan beam center
            score_best = 0.01
            do n2 = -ncentx, ncentx, cent_step
            do n1 = -ncenty, ncenty, cent_step

               do nz = 1, npsi
                   psi = (nz-npsi/2)*delpsi + psi_init
                   do ny = 1, nphi
                      phi = (ny-nphi/2)*delphi + phi_init
!                     if(phi .lt. pi/3*3 .or. phi .gt. pi/3*4) cycle
                      do nx = 1, ntheta
                         theta = (nx-ntheta/2)*deltheta+theta_init
                         call rotation_mat_from_rodrigues(psi, phi, theta, rot)
                         om = matmul(rot,UB_B)
                         call evaluate(om, x, y, dble(centx+n2), dble(centy+n1), dble(distance), nrefl, score1, score2)
                         if(score2 .gt. score_best)then
                             score_best = score2
                             centx_best = centx+n2
                             centy_best = centy+n1
                             best = (/score_best, psi, phi, theta, score1/)
                         end if
                      end do
                   end do
               end do

            end do
            end do

         else

            !scan beam center and unit cell
            score_best = 0.01
            do n2 = -ncentx, ncentx, cent_step
            do n1 = -ncenty, ncenty, cent_step

               do nz = 1, npsi
                   psi = (nz-npsi/2)*delpsi + psi_init
                   do ny = 1, nphi
                      phi = (ny-nphi/2)*delphi + phi_init
!                     if(phi .lt. pi/3*3 .or. phi .gt. pi/3*4) cycle
                      do nx = 1, ntheta

                         theta = (nx-ntheta/2)*deltheta+theta_init
                         call rotation_mat_from_rodrigues(psi, phi, theta, rot)

                         if(trim(adjustl(crystal_system)).eq.'Trigonal' .or. trim(adjustl(crystal_system)).eq.'Cubic')then

                            !only scan one axes
                            do n3 = -nint(cell_range/cell_step), nint(cell_range/cell_step)

                               call transformation_matrix((/uc(1)+n3*cell_step, uc(2)+n3*cell_step, &
                                                            uc(3)+n3*cell_step, uc(4), uc(5), uc(6)/), UB_B)

                               om = matmul(rot,UB_B)
                               call evaluate(om, x, y, dble(centx+n2), dble(centy+n1), dble(distance), nrefl, score1, score2)
                               if(score2 .gt. score_best)then
                                   score_best = score2
                                   centx_best = centx+n2
                                   centy_best = centy+n1
                                   ua_best = uc(1)+n3*cell_step
                                   ub_best = uc(2)+n3*cell_step
                                   uc_best = uc(3)+n3*cell_step
                                   best = (/score_best, psi, phi, theta, score1/)
                               end if

                            end do

                         else if(trim(adjustl(crystal_system)).eq.'Tetragonal' .or. trim(adjustl(crystal_system)).eq.&
                                 'Hexagonal')then

                            !scan two axes
                            do n4 = -nint(cell_range/cell_step), nint(cell_range/cell_step)
                            do n3 = -nint(cell_range/cell_step), nint(cell_range/cell_step)

                               call transformation_matrix((/uc(1)+n3*cell_step, uc(2)+n3*cell_step, &
                                                            uc(3)+n4*cell_step, uc(4), uc(5), uc(6)/), UB_B)

                               om = matmul(rot,UB_B)
                               call evaluate(om, x, y, dble(centx+n2), dble(centy+n1), dble(distance), nrefl, score1, score2)
                               if(score2 .gt. score_best)then
                                   score_best = score2
                                   centx_best = centx+n2
                                   centy_best = centy+n1
                                   ua_best = uc(1)+n3*cell_step
                                   ub_best = uc(2)+n3*cell_step
                                   uc_best = uc(3)+n4*cell_step
                                   best = (/score_best, psi, phi, theta, score1/)
                               end if

                            end do
                            end do

                         else  !for Triclinic, Monoclinic, Orthorhombic

                            !scan all three axes
                            do n5 = -nint(cell_range/cell_step), nint(cell_range/cell_step)
                            do n4 = -nint(cell_range/cell_step), nint(cell_range/cell_step)
                            do n3 = -nint(cell_range/cell_step), nint(cell_range/cell_step)

                               call transformation_matrix((/uc(1)+n3*cell_step, uc(2)+n4*cell_step, &
                                                            uc(3)+n5*cell_step, uc(4), uc(5), uc(6)/), UB_B)

                               om = matmul(rot,UB_B)
                               call evaluate(om, x, y, dble(centx+n2), dble(centy+n1), dble(distance), nrefl, score1, score2)
                               if(score2 .gt. score_best)then
                                   score_best = score2
                                   centx_best = centx+n2
                                   centy_best = centy+n1
                                   ua_best = uc(1)+n3*cell_step
                                   ub_best = uc(2)+n4*cell_step
                                   uc_best = uc(3)+n5*cell_step
                                   best = (/score_best, psi, phi, theta, score1/)
                               end if

                            end do
                            end do
                            end do

                         end if

                      end do
                   end do
               end do

            end do
            end do


         end if

         psi_best = best(2); phi_best = best(3); theta_best = best(4)
         uc = (/ua_best, ub_best, uc_best, uc(4), uc(5), uc(6)/)

         write(*,*)
         write(*,'(a)') "--------------------------------------------------------------"
         write(*,'(a)') "                 After prior information correction"
         write(*,*)
         write(*,'(a,F6.3,a,F8.3)') " Match rate: ", best(1),  "Positional deviation:", best(5)
         write(*,'(a)') "--------------------------------------------------------------"
         write(*,*)

      end if


      call rotation_mat_from_rodrigues(psi_best, phi_best, theta_best, rot)
      call transformation_matrix(uc, UB_B)
      om_best = matmul(rot,UB_B)

      write(*,*)
      write(*,'(a)') "--------------------------------------------------------------"
      write(*,'(a)') "                 Initial Orientation Matrix"
      write(*,*)
      write(*,'(3F10.6)')om_best(1,:)
      write(*,'(3F10.6)')om_best(2,:)
      write(*,'(3F10.6)')om_best(3,:)
      write(*,'(a)') "--------------------------------------------------------------"

      write(*,*)
      write(*,'(a)') "--------------------------------------------------------------"
      write(*,'(a)') "                 Preliminary results before refinement"
      write(*,*)
      write(*,'(18X,a)') "     UC_a    UC_b    UC_c     Psi    Phi    Theta  CenterX CenterY Distance "
      write(*,'(a,9F8.2)') " Inital parameters: ", uc(1), uc(2), uc(3), psi_best/pi*180, phi_best/pi*180,&
                  theta_best/pi*180, centx_best, centy_best, distance_best
      write(*,'(a,F8.3,a)') " Initial Residual Error:  ", best(5), " pixels"
      write(*,'(a)') "--------------------------------------------------------------"
      write(*,*)

      write(*,'(a)') "--------------------------------------------------------------"
      write(*,'(a)') "                 Iterative refinement process"
      write(*,*)

      !iterative refinement
      do niter = 1, 10

         nfound = 0
         do i = 1, nrefl
            call xy2q(x(i), y(i), dble(centx_best), dble(centy_best), dble(distance_best), qx, qy, qz)
            call q2hkl(om_best, qx, qy, qz, h, k, l)
            call hkl2xy(om_best, dble(centx_best), dble(centy_best), dble(distance_best), h, k, l, xo, yo)
            if(sqrt((xo-x(i))**2+(yo-y(i))**2) .le. pos_err)then
                nfound = nfound + 1
                sub_x(nfound) = x(i)
                sub_y(nfound) = y(i)
            end if
         end do

         write(*,'(a,I4)') " Cycle of refinement: ", niter
         write(*,'(I5,a,I5)') nfound, " reflections matched and used for refinement from total of", nrefl

         datas%nrefl = nfound
         datas%x = sub_x; datas%y = sub_y
         ptr = c_loc(datas)

         if(trim(adjustl(crystal_system)).eq.'Trigonal' .or. trim(adjustl(crystal_system)).eq.'Cubic')then

            allocate(xv(7)); allocate(stepv(7))

            min_fslv = fgsl_multimin_fminimizer_alloc(fgsl_multimin_fminimizer_nmsimplex2, 7_8)
            func = fgsl_multimin_function_init(my_minimum_function, 7_8, ptr)

            xv(1:7) = (/uc(1), psi_best, phi_best, theta_best, dble(centx_best), dble(centy_best), dble(distance_best)/)
            x_in = fgsl_vector_init(1.0_fgsl_double)
            status = fgsl_vector_align(xv, 7_8, x_in, 7_8, 0_fgsl_size_t, 1_fgsl_size_t)

            stepv(1:7) = (/0.01_fgsl_double, 0.005_fgsl_double, 0.005_fgsl_double, 0.005_fgsl_double, &
                           0.01_fgsl_double, 0.01_fgsl_double, 0.01_fgsl_double/)
            step_size = fgsl_vector_init(1.0_fgsl_double)
            status = fgsl_vector_align(stepv, 7_8, step_size, 7_8, 0_fgsl_size_t, 1_fgsl_size_t)

            status = fgsl_multimin_fminimizer_set(min_fslv, func, x_in, step_size)

            call fgsl_vector_free(x_in)
            x_in = fgsl_multimin_fminimizer_x(min_fslv)
            status = fgsl_vector_align(xptr, x_in)


            i = 0
            do
                i = i + 1
                status = fgsl_multimin_fminimizer_iterate(min_fslv)
                if(status /= fgsl_success .or. i > itmax)exit

            end do

            write(*,'(a,9F8.2)') " Refined parameters: ", xptr(1), xptr(1), xptr(1), xptr(2:4)/pi*180, xptr(5:7)
            write(*,'(a,F8.3,a)') " Final Residual Error:", fgsl_multimin_fminimizer_minimum(min_fslv), " pixels"
            write(*,*)

            call rotation_mat_from_rodrigues(xptr(2), xptr(3), xptr(4), rot)
            call transformation_matrix((/xptr(1), xptr(1), xptr(1), uc(4), uc(5), uc(6)/), UB_B)
            om_best = matmul(rot,UB_B)
            uc(1) = xptr(1); uc(2)= xptr(1); uc(3) = xptr(1)
            psi_best = xptr(2); phi_best = xptr(3); theta_best = xptr(4)
            centx_best = xptr(5); centy_best = xptr(6)
            distance_best = xptr(7)

            deallocate(xv, stepv)

         else if(trim(adjustl(crystal_system)).eq.'Tetragonal' .or. trim(adjustl(crystal_system)).eq.&
                 'Hexagonal')then

            allocate(xv(8)); allocate(stepv(8))

            min_fslv = fgsl_multimin_fminimizer_alloc(fgsl_multimin_fminimizer_nmsimplex2, 8_8)
            func = fgsl_multimin_function_init(my_minimum_function, 8_8, ptr)

            xv(1:8) = (/uc(1), uc(3), psi_best, phi_best, theta_best, dble(centx_best), dble(centy_best), dble(distance_best)/)
            x_in = fgsl_vector_init(1.0_fgsl_double)
            status = fgsl_vector_align(xv, 8_8, x_in, 8_8, 0_fgsl_size_t, 1_fgsl_size_t)

            stepv(1:8) = (/0.01_fgsl_double, 0.01_fgsl_double, 0.005_fgsl_double, 0.005_fgsl_double, 0.005_fgsl_double, &
                           0.01_fgsl_double, 0.01_fgsl_double, 0.01_fgsl_double/)
            step_size = fgsl_vector_init(1.0_fgsl_double)
            status = fgsl_vector_align(stepv, 8_8, step_size, 8_8, 0_fgsl_size_t, 1_fgsl_size_t)

            status = fgsl_multimin_fminimizer_set(min_fslv, func, x_in, step_size)

            call fgsl_vector_free(x_in)
            x_in = fgsl_multimin_fminimizer_x(min_fslv)
            status = fgsl_vector_align(xptr, x_in)


            i = 0
            do
                i = i + 1
                status = fgsl_multimin_fminimizer_iterate(min_fslv)
                if(status /= fgsl_success .or. i > itmax)exit

            end do

            write(*,'(a,9F8.2)') " Refined parameters: ", xptr(1), xptr(1),xptr(2), xptr(3:5)/pi*180, xptr(6:8)
            write(*,'(a,F8.3,a)') " Final Residual Error:", fgsl_multimin_fminimizer_minimum(min_fslv), " pixels"
            write(*,*)

            call rotation_mat_from_rodrigues(xptr(3), xptr(4), xptr(5), rot)
            call transformation_matrix((/xptr(1), xptr(1), xptr(2), uc(4), uc(5), uc(6)/), UB_B)
            om_best = matmul(rot,UB_B)
            uc(1) = xptr(1); uc(2) = xptr(1); uc(3) = xptr(2)
            psi_best = xptr(3); phi_best = xptr(4); theta_best = xptr(5)
            centx_best = xptr(6); centy_best = xptr(7)
            distance_best = xptr(8)

            deallocate(xv, stepv)

         else  !for Triclinic, Monoclinic, Orthorhombic

            allocate(xv(9)); allocate(stepv(9))

            min_fslv = fgsl_multimin_fminimizer_alloc(fgsl_multimin_fminimizer_nmsimplex2, 9_8)
            func = fgsl_multimin_function_init(my_minimum_function, 9_8, ptr)

            xv(1:9) = (/uc(1), uc(2), uc(3), psi_best, phi_best, theta_best, dble(centx_best), &
                        dble(centy_best), dble(distance_best)/)
            x_in = fgsl_vector_init(1.0_fgsl_double)
            status = fgsl_vector_align(xv, 9_8, x_in, 9_8, 0_fgsl_size_t, 1_fgsl_size_t)

            stepv(1:9) = (/0.01_fgsl_double, 0.01_fgsl_double, 0.01_fgsl_double, 0.005_fgsl_double, 0.005_fgsl_double, &
                           0.005_fgsl_double, 0.01_fgsl_double, 0.01_fgsl_double, 0.01_fgsl_double/)
            step_size = fgsl_vector_init(1.0_fgsl_double)
            status = fgsl_vector_align(stepv, 9_8, step_size, 9_8, 0_fgsl_size_t, 1_fgsl_size_t)

            status = fgsl_multimin_fminimizer_set(min_fslv, func, x_in, step_size)

            call fgsl_vector_free(x_in)
            x_in = fgsl_multimin_fminimizer_x(min_fslv)
            status = fgsl_vector_align(xptr, x_in)


            i = 0
            do
                i = i + 1
                status = fgsl_multimin_fminimizer_iterate(min_fslv)
                if(status /= fgsl_success .or. i > itmax)exit

            end do

            write(*,'(a,9F8.2)') " Refined parameters: ", xptr(1:3), xptr(4:6)/pi*180, xptr(7:9)
            write(*,'(a,F8.3,a)') " Final Residual Error:", fgsl_multimin_fminimizer_minimum(min_fslv), " pixels"
            write(*,*)

            call rotation_mat_from_rodrigues(xptr(4), xptr(5), xptr(6), rot)
            call transformation_matrix((/xptr(1), xptr(2), xptr(3), uc(4), uc(5), uc(6)/), UB_B)
            om_best = matmul(rot,UB_B)
            uc(1) = xptr(1); uc(2) = xptr(2); uc(3) = xptr(3)
            psi_best = xptr(4); phi_best = xptr(5); theta_best = xptr(6)
            centx_best = xptr(7); centy_best = xptr(8)
            distance_best = xptr(9)

            deallocate(xv, stepv)

         end if

      end do

      write(*,'(a)') "--------------------------------------------------------------"
      write(*,*)
      write(*,'(a)') "--------------------------------------------------------------"
      write(*,'(a)') "                 Refined Orientation Matrix"
      write(*,*)
      write(*,'(3F10.6)')om_best(1,:)
      write(*,'(3F10.6)')om_best(2,:)
      write(*,'(3F10.6)')om_best(3,:)
      write(*,'(a)') "--------------------------------------------------------------"

      nfound = 0
      open(20, file='best.txt', status='replace')
      open(21, file='SPOT-TMP.TXT', status='replace')
      do i = 1, nrefl
         call xy2q(x(i), y(i), centx_best, centy_best, distance_best, qx, qy, qz)
         call q2hkl(om_best, qx, qy, qz, h, k, l)
         call hkl2xy(om_best, centx_best, centy_best, distance_best, h, k, l, xo, yo)
         if(sqrt((xo-x(i))**2+(yo-y(i))**2) .le. pos_err)then
             write(20,'(4F10.2,3I5)') x(i), y(i), xo, yo, nint(h), nint(k), nint(l)
             nfound = nfound + 1
         else
             write(21,'(2F10.2)') x(i), y(i)
         end if
      end do
      close(20)
      close(21)

      write(*,*)
      write(*,'(a)') "--------------------------------------------------------------"
      write(*,'(a)') "                 Summary of reflection matching"
      write(*,*)
      write(*,'(I5,a,I5)') nfound, " reflections have been identified after refinement and written to file from total of", nrefl
      write(*,'(a)') "--------------------------------------------------------------"
      write(*,*)

      call fgsl_multimin_fminimizer_free(min_fslv)
      call fgsl_multimin_function_free(func)

  end if



  call MPI_FINALIZE(ierr)

  time2 = MPI_Wtime()

  if(node .eq. 0)then
      write(*,'(a,F10.2)') " Total Time used (seconds): ", time2-time1
      write(*,'(a)') "================================================================"
      write(*,*)
  end if
   
  stop
end program main
