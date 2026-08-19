module read_input_module
  use globals
  use iso_c_binding, only: c_char, c_null_char
  implicit none

  contains

  subroutine read_input(filename, listname, keyword, keyval, errorio, errorstr)
    character(*), intent(in) :: filename, listname, keyword
    logical(1), intent(inout):: errorio
    character(c_char), dimension(*), intent(out) :: errorstr
    character(len=20) :: f_errorstr
    character(len=200) :: line, trimline, val, thislist, key
    integer:: ios, n, i
    real(dp):: keyval
    logical(1) :: strfind

    ! just exit if a previous parameter was missing
    IF (.not. errorio) then

    f_errorstr = "Parameters loaded"
    strfind = .false.
    thislist = ''

    !------------------------------------------------------------------
    !  Open file and parse YAML key-value pairs by section
    !
    !  Format:
    !    section_name:          <- section header: no leading spaces, ends with ':'
    !      key: value  # comment  <- indented key-value pairs
    !    # comment lines         <- lines whose first non-space char is '#'
    !------------------------------------------------------------------
    OPEN(1, file = filename)

      DO WHILE(.TRUE.)
        read(1, '(A)', iostat=ios) line
        IF (ios /= 0) exit

        ! Trimmed line for comment/empty checks
        trimline = TRIM(ADJUSTL(line))

        ! Skip empty lines
        IF (LEN_TRIM(trimline) == 0) CYCLE

        ! Skip comment lines (first non-space character is '#')
        IF (trimline(1:1) == '#') CYCLE

        ! Detect section header: no leading space, ends with ':'
        IF (line(1:1) /= ' ') THEN
          n = LEN_TRIM(trimline)
          IF (trimline(n:n) == ':') THEN
            thislist = TRIM(trimline(1:n-1))
          END IF
          CYCLE
        END IF

        ! Key-value pair (indented line): only process in the matching section
        IF (TRIM(thislist) == TRIM(listname)) THEN
          line = trimline

          ! Remove inline comments
          n = INDEX(line, '#')
          IF (n /= 0) line = TRIM(ADJUSTL(line(1:n-1)))

          ! Find the first colon (key:value separator)
          n = INDEX(line, ':')
          IF (n == 0) CYCLE

          ! Isolate key and value
          key = TRIM(ADJUSTL(line(1:n-1)))
          val = TRIM(ADJUSTL(line(n+1:LEN(line))))

          ! Skip lines with no value (e.g. nested section headers)
          IF (LEN_TRIM(val) == 0) CYCLE

          ! If the key matches the keyword, read in the value
          IF (TRIM(key) == TRIM(keyword)) THEN
            read(val, *) keyval
            strfind = .true.
          END IF
        END IF
      END DO
      CLOSE(1)

    !
    ! If the keyword was not found, register the error
    !
    IF (.not. strfind) THEN
      print*, 'parameter ', TRIM(keyword), ' not defined for ', listname, ' in ', filename, '.'
      open(2, file='errorfile.txt', status='unknown')
      write(2,*) 'parameter ', TRIM(keyword), ' not defined for ', listname, ' in ', filename, '.'
      f_errorstr = TRIM(keyword) // "                        "
      errorio = .true.
    END IF

    ! Transfer the characters from Fortran to C
    do i = 1, len(f_errorstr)
      if (i <= len_trim(f_errorstr)) then
        errorstr(i) = f_errorstr(i:i)
      else
        errorstr(i) = c_null_char
      end if
    end do

    end if

  end subroutine read_input


end module read_input_module
