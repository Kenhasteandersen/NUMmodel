module read_input_module
  use globals
  use iso_c_binding, only: c_char, c_null_char
  implicit none
  
  contains
        
  subroutine read_input(filename, listname, keyword,keyval,errorio,errorstr)
    character(*), intent(in) :: filename, listname,keyword
    logical(1), intent(inout):: errorio 
    character(c_char), dimension(*), intent(out) :: errorstr
    character(len=20) :: f_errorstr
    character(len=100) :: line, val, thislist, key
    character(len=100) :: str_General
    character(len=100) :: str_Generalists_simple
    character(len=100) :: str_Generalists
    character(len=100) :: str_Diatoms_simple
    character(len=100) :: str_Diatoms
    character(len=100) :: str_Copepods_passive
    character(len=100) :: str_Copepods_active
    character(len=100) :: str_POM
    integer:: ios, n, i
    integer, parameter :: DRK = selected_real_kind (20)
    real(dp):: keyval
    real(dp):: io
    logical(1) :: strfind 
    
    ! just exit if a previos parameter was missing
    IF (.not. errorio) then

    
    f_errorstr ="Parameters loaded             "
    strfind=.false. 
    
 
    !------------------------------------------------------------------
    !  Open file and find keywords and parameter values
    !------------------------------------------------------------------
    str_General='! GENERAL PARAMETERS'
    str_Generalists_simple='! GENERALISTS SIMPLE INPUT PARAMETERS'
    str_Generalists='! GENERALISTS INPUT PARAMETERS'
    str_Diatoms_simple='! DIATOMS SIMPLE INPUT PARAMETERS'
    str_Diatoms='! DIATOMS INPUT PARAMETERS'
    str_Copepods_passive='! COPEPODS PASSIVE INPUT PARAMETERS'
    str_Copepods_active='! COPEPODS ACTIVE INPUT PARAMETERS'
    str_POM='! PARTICULATE ORGANIC MATTER (POM) INPUT PARAMETERS'
    thislist='no list defined yet'
    
    OPEN(1, file = filename)
    
      DO WHILE(.TRUE.)
        read(1, '(A)', iostat=ios) line
        IF (ios /= 0) exit
        
     !
    ! Check if we are on a group define statement
    !
          if (line(1:len(str_General)).eq.str_General) then
            thislist='general'
          else if (line(1:len(str_Generalists_simple)).eq.str_Generalists_simple) then
            thislist='generalists_simple'
          else if (line(1:len(str_Generalists)).eq.str_Generalists) then
            thislist='generalists'
          else if (line(1:len(str_Diatoms_simple)).eq.str_Diatoms_simple) then
            thislist='diatoms_simple'
          else if (line(1:len(str_Diatoms)).eq.str_Diatoms) then
            thislist='diatoms'
          else if (line(1:len(str_Copepods_passive)).eq.str_Copepods_passive) then
            thislist='copepods_passive'
          else if (line(1:len(str_Copepods_active)).eq.str_Copepods_active) then
            thislist='copepods_active'
          else if (line(1:len(str_POM)).eq.str_POM) then
            thislist='POM'
          end if   
    !
    ! Remove leading blank spaces
    !
          line=TRIM(ADJUSTL(line(1:LEN(line))))
    !
    ! If the line is not a comment and it is on the right list then:
    !
           IF ((line(1:1).ne.'!').and.(LEN_TRIM(line).gt.0).and.(TRIM(thislist).eq.TRIM(listname))) THEN
    !
    ! remove comments, if any
    !
            n=index(line,'!')
            IF (n /=0) THEN
              line=TRIM(ADJUSTL(line(1:n-1)))
            END IF
    !
    ! Find equal sign
    !
            n=index(line,'=')
    !
    ! Isolate keywords and values
    !
            key=TRIM(ADJUSTL(line(1:n-1)))
            val=TRIM(ADJUSTL(line(n+1:LEN(line))))
            
    !
    ! If the key matches the keyword, read in the value
    !
         
            IF ( trim(key) == trim(keyword) ) then
              read(val,*) keyval
              !print*, '    ', trim(keyword), '= ', keyval
              strfind=.true.
            END IF
           END IF
      END DO
      CLOSE(1)
      
    !
    ! If the keyword was not found, register the error
    !      
     IF (.not. strfind) THEN
       print*, 'parameter ', TRIM(keyword), ' not defined for ', listname, ' in ',filename, '.'
       open(2, file='errorfile.txt',status='unknown')
       write(2,*) 'parameter ', TRIM(keyword), ' not defined for ', listname, ' in ',filename, '.'
       f_errorstr =TRIM(keyword) // "                        "
       errorio=.true.
       
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


