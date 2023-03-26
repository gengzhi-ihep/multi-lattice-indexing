program main
    implicit none
    integer,parameter :: npats = 10
    real,parameter :: wave = 1.5
    integer :: pat_tot(2048,2048)
    real(8) :: om_tot(3,3,npats)
    character(50) :: command
    integer n, m, i, j, tmp

    pat_tot = 0

    do n = 1, npats

       write(*,*) "Pattern simulation : ", n

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
