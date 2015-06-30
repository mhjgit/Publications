             implicit real*8(a-h,o-z)
             open(unit=6,file='chemf_nn.dat')
             open(unit=7,file='chemf_pn.dat')
             open(unit=8,file='chemf_smn.dat')
             open(unit=9,file='chemf_ln.dat')
             open(unit=10,file='chemf_s0n.dat')
             open(unit=11,file='chemf_spn.dat')


             do i=1,100
                read(5,*) dens, x1, x2, x3, x4, x5, x6
                write(6,'(E12.6,2x,e12.6)') dens, x1
                write(7,'(E12.6,2x,e12.6)') dens, x2
                write(8,'(E12.6,2x,e12.6)') dens, x3
                write(9,'(E12.6,2x,e12.6)') dens, x4
                write(10,'(E12.6,2x,e12.6)') dens, x5
                write(11,'(E12.6,2x,e12.6)') dens, x6
             enddo
             end

