      subroutine chknam(str,len)

c-----------------------------------------------------------------------
c
c                   Check for a Valid File Name
c                   ***************************
c
c This subroutine takes the character string "str" of length "len" and
c removes all leading blanks and blanks out all characters after the
c first blank found in the string (leading blanks are removed first).
c
c
c
c-----------------------------------------------------------------------
      character(len=*), intent(inout) :: str
      integer itrim
c
c Remove leading blanks:
      str=adjustl(str)
c
c find first two blanks and blank out remaining characters:
      itrim=index(str,'   ')
      if (itrim > 0) str(itrim:)=' '
c
c Look for "-fi"
      itrim=index(str,'-fi')
      if (itrim > 0) str(itrim:)=' '
c
c Look for "\fi"
      itrim=index(str,'\fi')
      if (itrim > 0) str(itrim:)=' '
c
c Return with modified file name:
      return
      end
