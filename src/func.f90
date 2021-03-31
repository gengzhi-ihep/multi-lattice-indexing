!======================================================================================
! (c) IHEP 2021. All rights reserved.
! Author: Geng Zhi
! Institute of High Energy Physics, Chinese Academy of Science (IHEP, CAS). 
! If you have any problem with the program, please contact author with the
! following 
! email: gengz@ihep.ac.cn
!======================================================================================

module indexing
    use, intrinsic :: iso_c_binding
    use fgsl
    implicit none

    real, parameter :: pi = 3.1415926
    integer spg_num          !space group number
    character(20) crystal_system !crystal system string
    real lambda              !wavelength
    real pixsize             !pixel size
    real centx, centy        !beam center
    real distance            !sample-to-detector distance
    real(8) uc(6)            !unit cell parameter
    real pos_err             !positional tolerance to accept reflections as matched


    !used for fgsl multi-minimization
    type data
        integer(fgsl_size_t) :: nrefl
        real(fgsl_double) :: x(10000), y(10000)
    end type data
    integer,parameter :: itmax = 1000

    !used for integration
    integer, parameter :: inr_r = 4, mid_r = 5, out_r = 6
    integer, parameter :: Ndx = 4150, Ndy = 4371
    integer, parameter :: BM_IG = 1
    integer, parameter :: BM_BH = 2
    integer, parameter :: BM_BG = 3
    integer, parameter :: BM_PK = 4

    contains

    subroutine transformation_matrix(uc, UB_B)
        implicit none
        real(8) UB_B(3,3), uc(6), ax, ay, az, bx, by, bz, cx, cy, cz, vol
        integer is(3), js(3), flag

        vol = uc(1)*uc(2)*uc(3)*sqrt(1-cos(uc(4)/180.*pi)**2-cos(uc(5)/180.*pi)**2-cos(uc(6)/180.*pi)**2+ &
                                     2*cos(uc(4)/180.*pi)*cos(uc(5)/180.*pi)*cos(uc(6)/180.*pi))
        ax = uc(1); ay = 0; az = 0
        bx = uc(2)*cos(uc(6)/180.*pi); by = uc(2)*sin(uc(6)/180.*pi); bz = 0
        cx = uc(3)*cos(uc(5)/180.*pi)
        cy = uc(3)*(cos(uc(4)/180.*pi)-cos(uc(5)/180.*pi)*cos(uc(6)/180.*pi))/sin(uc(6)/180.*pi)
        cz = vol/uc(1)/uc(2)/sin(uc(6)/180.*pi)
        UB_B = reshape((/ax, bx, cx, ay, by, cy, az, bz, cz/), (/3,3/))

        call brinv(UB_B, 3, flag, is , js)
       
        return
    end subroutine

    subroutine rotation_mat_from_rodrigues(psi, phi, theta, mat)
        implicit none
        real(8) psi, phi, theta
        real(8) mat(3,3), nx, ny, nz, sw, cw
        nx = sin(psi)*cos(phi)
        ny = sin(psi)*sin(phi)
        nz = cos(psi)
        sw = sin(theta); cw = cos(theta)
        mat(1,1) = nx**2*(1-cw) + cw
        mat(1,2) = nx*ny*(1-cw) - nz*sw
        mat(1,3) = nx*nz*(1-cw) + ny*sw
        mat(2,1) = nx*ny*(1-cw) + nz*sw
        mat(2,2) = ny**2*(1-cw) + cw
        mat(2,3) = ny*nz*(1-cw) - nx*sw
        mat(3,1) = nx*nz*(1-cw) - ny*sw
        mat(3,2) = ny*nz*(1-cw) + nx*sw
        mat(3,3) = nz**2*(1-cw) + cw
        return
    end subroutine

    subroutine xy2q(xi, yi, centx, centy, distance, qx, qy, qz)
        implicit none
        real(8) xi, yi, qx, qy, qz, x, y, D
        real(8) centx, centy, distance
        x = (xi-centx)*pixsize
        y = (yi-centy)*pixsize
        D = sqrt(x**2+y**2+distance**2)
        qx = x/lambda/D
        qy = y/lambda/D
        qz = distance/lambda/D - 1/lambda
        return
    end subroutine

    subroutine resolution(xi, yi, centx, centy, distance, res)
        implicit none
        real(8) xi, yi, centx, centy, distance, res, theta
        theta = atan2(sqrt((xi-centx)**2+(yi-centy)**2)*pixsize, distance)/2
        res = lambda/(2*sin(theta))
        return
    end subroutine

    subroutine q2hkl(om_in, qx, qy, qz, h, k, l)
        implicit none
        real(8) qx, qy, qz, h, k, l, om_in(3,3), om(3,3)
        integer is(3), js(3), flag

        om = om_in
        call brinv(om, 3, flag, is , js)
        h = om(1,1)*qx + om(1,2)*qy + om(1,3)*qz
        k = om(2,1)*qx + om(2,2)*qy + om(2,3)*qz
        l = om(3,1)*qx + om(3,2)*qy + om(3,3)*qz

        return
    end subroutine

    subroutine hkl2xy(om, centx, centy, distance, hi, ki, li, x, y)
        implicit none
        real(8) om(3,3), hi, ki, li, x, y, qx, qy, qz
        real(8) centx, centy, distance
        integer h, k, l
        h = nint(hi); k = nint(ki); l = nint(li)
        qx = om(1,1)*h + om(1,2)*k + om(1,3)*l
        qy = om(2,1)*h + om(2,2)*k + om(2,3)*l
        qz = om(3,1)*h + om(3,2)*k + om(3,3)*l
        x = (qx*distance/(1/lambda+qz))/pixsize + centx
        y = (qy*distance/(1/lambda+qz))/pixsize + centy
        return
    end subroutine

    subroutine evaluate(om, x, y, centx, centy, distance, nref, score1, score2)
        implicit none
        real(8) om(3,3),x(10000), y(10000), score1, score2, qx, qy, qz, h, k, l, xo, yo
        real(8) centx, centy, distance
        integer nref, n, sum1

        real(8) sum0, err

        sum0 = 0
        sum1 = 0
        do n = 1, nref
            call xy2q(x(n),y(n), centx, centy, distance, qx, qy, qz)
            call q2hkl(om, qx, qy, qz, h, k, l)
            call hkl2xy(om, centx, centy, distance, h, k, l, xo, yo)
            err = sqrt((x(n)-xo)**2+(y(n)-yo)**2)
            if(err .le. pos_err)then
                sum0 = sum0 + err
                sum1 = sum1 + 1
            end if
        end do

        score1 = sum0/sum1
        score2 = sum1*1.0/nref

        return
    end subroutine

!   p : uc_a, uc_c, psi, phi, theta
    function my_minimum_function(x, p) bind(c)
        type(c_ptr), value :: x, p
        type(fgsl_vector) :: f_x
        type(data), pointer :: p_p
        real(8), pointer :: p_x(:)
        real(8) :: my_minimum_function
        integer :: status
        integer :: n, total
        real(8) om(3,3), rot(3,3), UB_B(3,3)
        real(8) qx, qy, qz, h, k, l, xo, yo, err

        call fgsl_obj_c_ptr(f_x, x)
        status = fgsl_vector_align(p_x, f_x)

        call c_f_pointer(p, p_p)

        if(trim(adjustl(crystal_system)).eq.'Trigonal' .or. trim(adjustl(crystal_system)).eq.'Cubic')then

           call rotation_mat_from_rodrigues(p_x(2), p_x(3), p_x(4), rot)
           call transformation_matrix((/p_x(1), p_x(1), p_x(1), uc(4), uc(5), uc(6)/), UB_B)
           om = matmul(rot, UB_B)

           my_minimum_function = 0

           total = 0
           do n = 1, p_p%nrefl
              call xy2q(p_p%x(n), p_p%y(n), p_x(5), p_x(6), p_x(7), qx, qy, qz)
              call q2hkl(om, qx, qy, qz, h, k, l)
              call hkl2xy(om, p_x(5), p_x(6), p_x(7), h, k, l, xo, yo)
              err = sqrt((p_p%x(n)-xo)**2+(p_p%y(n)-yo)**2)
              total = total + 1
              my_minimum_function = my_minimum_function + err
           end do

           my_minimum_function = my_minimum_function / total

        else if(trim(adjustl(crystal_system)).eq.'Tetragonal' .or. trim(adjustl(crystal_system)).eq.&
                'Hexagonal')then

           call rotation_mat_from_rodrigues(p_x(3), p_x(4), p_x(5), rot)
           call transformation_matrix((/p_x(1), p_x(1), p_x(2), uc(4), uc(5), uc(6)/), UB_B)
           om = matmul(rot, UB_B)

           my_minimum_function = 0

           total = 0
           do n = 1, p_p%nrefl
              call xy2q(p_p%x(n), p_p%y(n), p_x(6), p_x(7), p_x(8), qx, qy, qz)
              call q2hkl(om, qx, qy, qz, h, k, l)
              call hkl2xy(om, p_x(6), p_x(7), p_x(8), h, k, l, xo, yo)
              err = sqrt((p_p%x(n)-xo)**2+(p_p%y(n)-yo)**2)
              total = total + 1
              my_minimum_function = my_minimum_function + err
           end do

           my_minimum_function = my_minimum_function / total

        else  !for Triclinic, Monoclinic, Orthorhombic

           call rotation_mat_from_rodrigues(p_x(4), p_x(5), p_x(6), rot)
           call transformation_matrix((/p_x(1), p_x(2), p_x(3), uc(4), uc(5), uc(6)/), UB_B)
           om = matmul(rot, UB_B)

           my_minimum_function = 0

           total = 0
           do n = 1, p_p%nrefl
              call xy2q(p_p%x(n), p_p%y(n), p_x(7), p_x(8), p_x(9), qx, qy, qz)
              call q2hkl(om, qx, qy, qz, h, k, l)
              call hkl2xy(om, p_x(7), p_x(8), p_x(9), h, k, l, xo, yo)
              err = sqrt((p_p%x(n)-xo)**2+(p_p%y(n)-yo)**2)
              total = total + 1
              my_minimum_function = my_minimum_function + err
           end do

           my_minimum_function = my_minimum_function / total

        end if

        return
    end function

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
     
    

   !for matrix inverse
   subroutine brinv(a,n,l,is,js)
      implicit none
      integer n
      real(8) a(n,n)
      integer is(n), js(n)
      real t, d
      integer i, j, k, l
      l = 1

      do k = 1, n
         d = 0

         do i = k, n
         do j = k, n
            if(abs(a(i,j)) .gt. d)then
               d = abs(a(i,j))
               is(k) = i
               js(k) = j
            end if
         end do
         end do

         if(d+1 .eq. 1.0)then
            l = 0
            write(*,*) 'ERR**Not INV'
            return
         end if

         do j = 1, n
            t = a(k,j)
            a(k,j) = a(is(k),j)
            a(is(k),j) = t
         end do

         do i = 1, n
            t = a(i,k)
            a(i,k) = a(i,js(k))
            a(i,js(k)) = t
         end do
         a(k,k) = 1/a(k,k)

         do j = 1, n
            if(j .ne. k)then
               a(k,j) = a(k,j)*a(k,k)
            end if
         end do

         do i = 1, n
            if(i .ne. k)then
               do j = 1, n
                  if(j .ne. k)then
                     a(i,j)=a(i,j)-a(i,k)*a(k,j)
                  end if
               end do
            end if
         end do

         do i = 1, n
            if(i .ne. k)then
               a(i,k) = -a(i,k)*a(k,k)
            end if
         end do
         
      end do

      do k = n,1,-1
         do j = 1, n
            t = a(k,j)
            a(k,j) = a(js(k),j)
            a(js(k),j) = t
         end do
         do i = 1, n
            t = a(i,k)
            a(i,k) = a(i,is(k))
            a(i,is(k)) = t
         end do
      end do

      return
   end subroutine brinv

    !split line into key-value from input configure text
    subroutine extract_values(line, key, val)
        implicit none
        character(200) line
        character(20)  val, key
        integer start, endl, n

        endl = len_trim(line)

        do n = 1, len_trim(line)
           if(line(n:n) .eq. '=')then
               start = n
           else if(line(n:n) .eq. '#')then
               endl = n - 1
               exit
           end if
        end do

        key = trim(adjustl(line(1:start-1)))

        val = trim(adjustl(line(start+1:endl)))

        return
    end subroutine


end module indexing
