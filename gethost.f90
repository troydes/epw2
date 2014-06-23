SUBROUTINE gethost

#include "f_defs.h"
        
        integer st!,sleepi
        character*80 host

        sleepi = 0
        st = hostnm (host)

        print *, " The current host is " , host
        print *, " The current process ID is ", getpid()
        !do while (sleepi.eq.0)
        !        call sleep(5)
        !end do

END SUBROUTINE gethost
