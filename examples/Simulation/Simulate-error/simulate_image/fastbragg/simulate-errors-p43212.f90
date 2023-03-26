program main
    implicit none
    integer,parameter :: npats = 5
    real,parameter :: wave = 1.5
    integer :: pat_tot(2048,2048)
    real(8) :: om_tot(3,3,npats)
    character(50) :: command, tmp1
    character(100) :: line
    character(20) :: num, ctmp2, ctmp3, cnpats
    integer n, m, i, j, tmp, err
    
    real(8),external :: ZBQLUAB
    real(8) tmp2, tmp3

    pat_tot = 0

    call ZBQLINI(0)

    open(30,file='parameters.txt',status='replace')

    do n = 1, npats

       write(*,*) "Pattern simulation : ", n
       
       tmp2 = 1024 + ZBQLUAB(dble(-4), dble(4))
       tmp3 = 1024 + ZBQLUAB(dble(-4), dble(4))
       write(ctmp2,'(F10.2)') tmp2
       write(ctmp3,'(F10.2)') tmp3

       open(20,file='2sim-run-cp.sh',status='old')
       open(21,file='2sim-run.sh',status='replace')
       do m = 1, 5
           read(20,'(a)',iostat=err) line
           if (err .ne. 0) exit
           if (m .eq. 3) then
              line = "-ORGX "//trim(adjustl(ctmp2))//" -ORGY "//trim(adjustl(ctmp3))//" \"
           end if
           write(21,'(a)') trim(line)
       end do
       close(20)
       close(21)

       write(30,*) "Crystal", n
       write(30,*) "Beamcenter ", tmp2, tmp3


       tmp2 = 77.51 + ZBQLUAB(dble(-0.8), dble(0.8))
       tmp3 = 37.42 + ZBQLUAB(dble(-0.8), dble(0.8))
       write(ctmp2,'(F10.2)') tmp2
       write(ctmp3,'(F10.2)') tmp3
       open(20,file='1generate-OM-cp.sh',status='old')
       open(21,file='1generate-OM.sh',status='replace')
       do m = 1, 5
           read(20,'(a)',iostat=err) line
           if (err .ne. 0) exit
           if (m .eq. 2) then
              line = "CELL "//trim(adjustl(ctmp2))//" "//trim(adjustl(ctmp2))//" "//trim(adjustl(ctmp3))//" 90.0 90.0 90.0"
           end if
           write(21,'(a)') trim(line)
       end do
       close(20)
       close(21)


       write(30,*) "CELL ", tmp2, tmp2, tmp3, 90.0, 90.0, 90.0


       command = "./1generate-OM.sh >> om.log"
       call system(command)

       open(20,file='A.mat',status='old')
       do m = 1, 3
           read(20,*) om_tot(m,:,n)
       end do
       close(20)

       command = "./2sim-run.sh >> sim.log"
       call system(command)

       command = "./3noise.sh >> noise.log"
       call system(command)

       write(num,'(I4)') n
       write(cnpats,'(I4)') npats
       tmp1 = "sim"//trim(adjustl(cnpats))//"-"//trim(adjustl(num))//".img"
       command = "cp noiseimage.img  "// tmp1 
       call system(command)

       command = "./read_img"
       call system(command)

       open(20,file='new.txt',status='old')
       do j = 1, 2048
       do i = 1, 2048
           read(20,*) tmp
           pat_tot(i,j) = pat_tot(i,j) + tmp - 30
       end do
       end do
       close(20)

    end do
    close(30)

    open(20,file='om.txt',status='replace')
    do n = 1, npats
        write(20,*) "Pattern :", n
    do m = 1, 3
        write(20,*) om_tot(m,:,n)/wave
    end do
    end do
    close(20)

    open(20,file='new.txt',status='replace')
    do j = 1, 2048
    do i = 1, 2048
          write(20,*) pat_tot(i,j)
    end do
    end do
    close(20)

    command = "./write_hdf5"
    call system(command)

    stop
end program
