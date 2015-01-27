        real a, array(10)

C        data array/1.,2.,3.,4.,5.,6.,7.,8.,9.,0./

        open(12, file='unform.dat', recl=40, access='direct',
     *           form='unformatted')

C        open(6, file='unform_read.dat', recl=40, access='direct',
C     *           form='formatted')

C        a=10.0

C        write(6,*) 'Files opened'
C        write(12,rec=2) a
C        write(12,rec=1) (array(k), k=1,10)

        read(12,rec=2) a
C        write(6,100,rec=2) a
        write(6,100) a
        read(12,rec=1) (array(k), k=1,10)
C        write(6,110,rec=1) (array(k), k=1,10)
        write(6,110) (array(k), k=1,10)

100     format(1x,f8.3)
110     format(1x,10(f3.1,1x))
        END
