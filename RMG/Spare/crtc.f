c************************** SVN Revision Information **************************
c **    $Id$    **
c******************************************************************************
 
c ##  Real time clock functions
c ##

      REAL function CRTC()
      REAL rtc

      CRTC = rtc()

      return
      end

      integer function CIRTC()
      integer irtc

      CIRTC = irtc()

      return
      end

