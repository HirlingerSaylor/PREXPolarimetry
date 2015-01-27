      program main
      real result
      integer success
      WRITE(6,*) ' Start QFSRAD'
      call qfsrad(6.,1.,30./180.*3.1415926,0,result,success)
      WRITE(6,*) ' Result=',result
      print*,success,result
      end
