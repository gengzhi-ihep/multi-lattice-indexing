!======================================================================================
! (c) IHEP 2021. All rights reserved.
! Author: Geng Zhi
! Institute of High Energy Physics, Chinese Academy of Science (IHEP, CAS). 
! If you have any problem with the program, please contact author with the
! following 
! email: gengz@ihep.ac.cn
!======================================================================================

module sfx
    use, intrinsic :: iso_c_binding
    use fgsl
    implicit none

    !used for integration
    integer, parameter :: inr_r = 4, mid_r = 5, out_r = 6
    integer, parameter :: Ndx = 4150, Ndy = 4371
    integer, parameter :: BM_IG = 1
    integer, parameter :: BM_BH = 2
    integer, parameter :: BM_BG = 3
    integer, parameter :: BM_PK = 4

    !struct body for panel
    type panel
       integer min_fs, min_ss, max_fs, max_ss !vertex of this panel
       real(8) fsx, fsy, ssx, ssy             !transformation matrix from panel system to cartesian system  
       real(8) corner_x, corner_y             !locations of this panel in cartesian (pixel)
       real(8) phi                            !rotation angle of this panel
    end type panel

    !locations of peaks based on panel system and centerd on the whole detector
    !with origin as left-bottom in unit of pixel
    type peak
       real(8) fs, ss                       !observed locations of peaks on fs-ss(in unit of pixel)
       real(8) xs, ys                       !observed locations in cartesian(in unit of pixel)
    end type peak

    !structure body of one image containing peaks and to-be-refined parameters
    !NOTE : only peaks with existed predicted miller indexes will be used for
    !geometry refinement
    type crystal
       character(len=80) image_name         !name of this image
       integer num_peaks                    !total numbers of observed peaks
       type(peak),allocatable :: obs(:)     !observed peaks on one image
       real(8) lambda                       !wavelength of this image (A)
       real(8) dis                          !distance between crystal and detector
    end type crystal


    contains

    !read in 64 panels to geom
    subroutine read_panels(filename,geom)
       implicit none
       character(len=*) filename
       type(panel) geom(:)
       integer eof, i, j, n
       character(len=50) line
       open(10,file=filename,status='old')
       read(10,'(/////)')  !first omit the first 4 rows

       do n = 1, size(geom) !total panels
       do i = 1, 9          !each panel has 8 parameters

          read(10,'(a)') line
          do j = 1, len_trim(adjustl(line))  !find parameter name
             if(line(j:j).eq.'/')exit
          end do

          if(line(j+1:j+6) .eq. 'min_fs')then
             read(line(j+9:),*) geom(n)%min_fs
          else if(line(j+1:j+6) .eq. 'min_ss')then
             read(line(j+9:),*) geom(n)%min_ss
          else if(line(j+1:j+6) .eq. 'max_fs')then
             read(line(j+9:),*) geom(n)%max_fs
          else if(line(j+1:j+6) .eq. 'max_ss')then
             read(line(j+9:),*) geom(n)%max_ss
          else if(line(j+1:j+8) .eq. 'corner_x')then
             read(line(j+11:),*) geom(n)%corner_x
          else if(line(j+1:j+8) .eq. 'corner_y')then
             read(line(j+11:),*) geom(n)%corner_y
          else if(line(j+1:j+2) .eq. 'fs')then
             read(line(j+6:j+14),*) geom(n)%fsx
             read(line(j+17:j+25),*) geom(n)%fsy
          else if(line(j+1:j+2) .eq. 'ss')then
             read(line(j+6:j+14),*) geom(n)%ssx
             read(line(j+17:j+25),*) geom(n)%ssy
          end if

       end do
       geom(n)%phi = atan2(geom(n)%fsy,geom(n)%fsx)
       end do

       close(10)
       return
    end subroutine read_panels


    !read peaks and braggs from one image
    subroutine read_one_image(fileunit,cryst)
       implicit none
       integer fileunit, n
       type(crystal) :: cryst

       !temporary variables
       character(len=100) line
       real(8) tmp1, tmp2, dis

       !markers for begin and end of crystal and chunk for reading
       character(len=17) :: chunk_begin  = '----- Begin chunk'
       character(len=15) :: chunk_end    = '----- End chunk'
       character(len=22) :: peak_begin   = 'Peaks from peak search'

       outer : do while(.true.)

          read(fileunit,'(a)') line

          !start of one chunk
          if(line(1:17) .eq.chunk_begin)then

             inner:do while(.true.)

                read(fileunit,'(a)') line

                !jump out of the total loops or cycle to find new crystal
                if(line(1:15).eq.chunk_end)then
                    exit outer
                end if

                !image name
                if(line(1:15) .eq. 'Image filename:')then
                   cryst%image_name = line(17:)
                end if

                !read wavelength
                if(line(1:16) .eq. 'photon_energy_eV')then
                   read(line(19:),*) cryst%lambda
                end if

                !read distance if existed
                if(line(1:16) .eq. 'hdf5/distance_mm')then
                   read(line(19:),*) dis
                end if

                !read number of peaks
                if(line(1:9) .eq. 'num_peaks')then
                   read(line(13:),*) cryst%num_peaks
                   allocate(cryst%obs(cryst%num_peaks))
                end if

                !read in strong peaks
                if(line(1:22) .eq. peak_begin)then
                   read(fileunit,'(a)') line
                   do n = 1, cryst%num_peaks
                      read(fileunit,*)cryst%obs(n)%fs, cryst%obs(n)%ss, tmp1, tmp2
                   end do
                end if
             end do inner

          end if

       end do outer

       cryst%lambda = 6.62606896*2.99792458D3/(cryst%lambda)/1.6021773
       cryst%dis = dis/1000

       return
    end subroutine read_one_image

    !turn locations of peaks in fs-ss(in pixel) to cartesian coordinates xy(in pixel)
    subroutine panel2cart(fs,ss,cryst,geom,xs,ys)
       implicit none
       real(8) fs, ss, xs, ys, fs1, ss1
       type(panel) geom(:)
       type(crystal) cryst
       integer i, ind
       logical found

       !first find which panel this peak is in
       found = .false.
       do i = 1, size(geom)
          if(fs .gt. geom(i)%max_fs) cycle
          if(fs .lt. geom(i)%min_fs) cycle
          if(ss .gt. geom(i)%max_ss) cycle
          if(ss .lt. geom(i)%min_ss) cycle
          found = .true.
          ind = i
          exit
       end do
       if(.not. found) then
          write(6,'(2a)')'Error: image has invalid peak not on any panel!',cryst%image_name; stop
       end if
       !then turn location origined in panel corner in unit of pixel
       fs1 = fs - geom(ind)%min_fs
       ss1 = ss - geom(ind)%min_ss
       !now turn fs ss into catesian system xs, ys
       xs = fs1*geom(ind)%fsx + ss1*geom(ind)%ssx
       ys = fs1*geom(ind)%fsy + ss1*geom(ind)%ssy
       !and origined in detector center
       xs = xs + geom(ind)%corner_x
       ys = ys + geom(ind)%corner_y
       return
    end subroutine panel2cart


    subroutine integrate_rings(xi, yi, nref, array, intensities, sigmas)
        implicit none
        real(8) :: xi(10000), yi(10000), array(Ndx, Ndy)
        real(8) :: intensities(10000), sigmas(10000), rs
        integer boxmask(0:2*out_r, 0:2*out_r)
        real(8),target :: bgm_org(3,3), v_org(3), s_tmp(3,3), s_org(3), shifts_org(3)
        real(8) pks_p, pks_q, vmax, intensity, bgmean, sig2_bg, sig2_poisson, sigma
        type(fgsl_matrix) bgm, s_vec
        type(fgsl_vector) v, s_val, shifts
        integer nref, m, n, nn, i, p, q, x, y, cfs, css, status
        real(8),parameter :: max_condition = 1e6

        bgm = fgsl_matrix_init(1.0_fgsl_double)
        status = fgsl_matrix_align(bgm_org,3_8,3_8,3_8,bgm)
        v = fgsl_vector_init(1.0_fgsl_double)
        status = fgsl_vector_align(v_org,3_8,v,3_8,0_fgsl_size_t,1_fgsl_size_t)
        s_val = fgsl_vector_init(1.0_fgsl_double)
        status = fgsl_vector_align(s_org,3_8,s_val,3_8,0_fgsl_size_t,1_fgsl_size_t)
        s_vec = fgsl_matrix_init(1.0_fgsl_double)
        status = fgsl_matrix_align(s_tmp,3_8,3_8,3_8,s_vec)
        shifts = fgsl_vector_init(1.0_fgsl_double)
        status = fgsl_vector_align(shifts_org,3_8,shifts,3_8,0_fgsl_size_t,1_fgsl_size_t)

        !set up box mask
        do p = 0, 2*out_r
        do q = 0, 2*out_r
           rs = (p-out_r)**2 + (q-out_r)**2
           if(rs .gt. out_r**2)then
              boxmask(p,q) = BM_IG      !unused
           else
              if(rs .ge. mid_r**2)then
                 boxmask(p,q) = BM_BG   !Background
              else if(rs .le. inr_r**2)then
                 boxmask(p,q) = BM_PK   !Peak
              else
                 boxmask(p,q) = BM_IG   !unused
              end if
           end if
        end do
        end do

        intensities = 0; sigmas = 0

        !integrate
        outer : do n = 1, nref

           x = nint(xi(n)); y = nint(yi(n))

           do p = 0, 2*out_r
           do q = 0, 2*out_r
               if((x-out_r+p) .le. 0 .or. (y-out_r+q) .le. 0) cycle outer
               if((x-out_r+p) .gt. Ndx .or. (y-out_r+q) .gt. Ndy) cycle outer
           end do
           end do

           bgm_org = 0; v_org = 0; pks_p = 0; pks_q = 0; m = 0

           do p = 0, 2*out_r
           do q = 0, 2*out_r
              if(boxmask(p,q) .eq. BM_BG)then
                 bgm_org(1,1) = bgm_org(1,1) + p*p
                 bgm_org(2,1) = bgm_org(2,1) + p*q
                 bgm_org(3,1) = bgm_org(3,1) + p
                 bgm_org(1,2) = bgm_org(1,2) + p*q
                 bgm_org(2,2) = bgm_org(2,2) + q*q
                 bgm_org(3,2) = bgm_org(3,2) + q
                 bgm_org(1,3) = bgm_org(1,3) + p
                 bgm_org(2,3) = bgm_org(2,3) + q
                 bgm_org(3,3) = bgm_org(3,3) + 1

                 cfs = x - out_r + p
                 css = y - out_r + q
                 v_org(1) = v_org(1) + p*array(cfs, css)
                 v_org(2) = v_org(2) + q*array(cfs, css)
                 v_org(3) = v_org(3) + array(cfs, css)
              else if(boxmask(p,q) .eq. BM_PK)then
                 pks_p = pks_p + p
                 pks_q = pks_q + q
                 m = m + 1
              end if
           end do
           end do

           !then LUV decomposition
           status = fgsl_linalg_sv_decomp_jacobi(bgm,s_vec,s_val)

           !next eigenvalue filtering
           vmax = 0
           do i = 1, 3
              if(s_org(i) .gt. vmax) vmax = s_org(i)
           end do
           do i = 1, 3
              if(s_org(i) .lt. vmax/max_condition) s_org(i) = 0
           end do

           status = fgsl_linalg_sv_solve(bgm,s_vec,s_val,v,shifts)

           !finally integrate intensity subtracted by background
           intensity = 0; bgmean = 0; sig2_bg = 0; nn = 0
           do p = 0, 2*out_r
           do q = 0, 2*out_r
              cfs = x - out_r + p
              css = y - out_r + q
              if(boxmask(p,q) .eq. BM_PK)then
                 intensity = intensity + array(cfs,css)
              else if(boxmask(p,q) .eq. BM_BG)then
                 bgmean = bgmean + array(cfs,css)
                 nn = nn + 1
              end if
           end do
           end do

           intensity=intensity-shifts_org(1)*pks_p-shifts_org(2)*pks_q-shifts_org(3)*m
           bgmean = bgmean / nn

           do p = 0, 2*out_r
           do q = 0, 2*out_r
              cfs = x - out_r + p
              css = y - out_r + q
              if(boxmask(p,q) .eq. BM_BG)then
                 sig2_bg = sig2_bg + (array(cfs,css)-bgmean)**2
              end if
           end do
           end do

           sigma = sqrt(m*sig2_bg / nn)
           intensities(n) = intensity
           sigmas(n) = sigma

        end do outer

        call fgsl_matrix_free(bgm)
        call fgsl_matrix_free(s_vec)
        call fgsl_vector_free(v)
        call fgsl_vector_free(s_val)
        call fgsl_vector_free(shifts)

        return
    end subroutine

end module
