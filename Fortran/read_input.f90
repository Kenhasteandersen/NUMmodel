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
    character(len=200) :: line, val, thislist, key
    character(len=100) :: str_General
    character(len=100) :: str_Generalists_simple
    character(len=100) :: str_Generalists
    character(len=100) :: str_Diatoms_simple
    character(len=100) :: str_Diatoms
    character(len=100) :: str_Copepods_passive
    character(len=100) :: str_Copepods_active
    character(len=100) :: str_POM
    character(len=100) :: str_Prokaryote
    integer:: ios, n, i
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
    str_Prokaryote='! PROKARYOTE INPUT PARAMETERS'
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
          else if (line(1:len(str_Prokaryote)).eq.str_Prokaryote) then
            thislist='prokaryote'
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
  
  subroutine read_inputX(filename, listname, keyword,keyval,k,errorio,errorstr)
    character(*), intent(in) :: filename, listname,keyword
    integer, intent(in):: k
    logical(1), intent(inout):: errorio 
    character(c_char), dimension(*), intent(out) :: errorstr
    character(len=20) :: f_errorstr
    character(len=100) :: line, val, thislist, key, val2
    character(len=100) :: str_General
    character(len=100) :: str_Generalists_simple
    character(len=100) :: str_Generalists
    character(len=100) :: str_Diatoms_simple
    character(len=100) :: str_Diatoms
    character(len=100) :: str_Copepods_passive
    character(len=100) :: str_Copepods_active
    character(len=100) :: str_POM
    character(len=100) :: str_Prokaryote
    integer:: ios, n, i
    real(dp):: keyval
    real(dp):: io
    logical(1) :: strfind 
    integer :: nargs, start, theend,kk
    character(len=100) :: word,repeatedValue
    integer :: spacePos, tabPos, repeatCount
    
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
    str_Prokaryote='! PROKARYOTE INPUT PARAMETERS'
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
          else if (line(1:len(str_Prokaryote)).eq.str_Prokaryote) then
            thislist='prokaryote'
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
            val=line(n+1:)
            
            
    !
    ! If the key matches the keyword, read in the value number k
    !
            
            IF ( trim(key) == trim(keyword) ) then
            !print*, 'valStart=', val    
            
            ! Check if the word is a repetition notation (e.g., 3*0.8d)
            

              if (LEN_TRIM(val) > 0 .and. INDEX(val, '*') > 0) then
               n = INDEX(val, '*')
                READ(val(1:n-1), *) repeatCount   ! Read the number of repetitions into a temporary variable
                
                repeatedValue = trim(adjustl(val(n+1:)))
                ! Construct the final value by repeating the repeatedValue
               
                val=repeatedValue
                
      
                do i = 2, repeatCount
                 val = char(32) // trim(val) // char(32) // trim(repeatedValue) // char(32)
                !  
               end do
               !print*, 'valNy=', val    
              end if
              
                     
              ! Loop through the different words
              start = 1
              do kk = 1, k
                if (kk.gt.1) then
                  start=theend
                end if
                !Find the first occurence of a character that is not space or tab
                do while ((val(start+1:start+1) == CHAR(9) .or. val(start+1:start+1) == CHAR(32)) .and. (start <= LEN(val)))
                    start = start + 1
                end do
                    start = start+ 1
                !print*, 'startnr=',start
                ! Find the first space or tab after that
                spacePos = INDEX(val(start:len(val)), CHAR(32))
                tabPos = INDEX(val(start:len(val)), CHAR(9))
                !print*, 'spacePos=',spacePos
                !print*, 'tabPos=',tabPos
             
                if (spacePos == 0) then
                   theend = tabPos
                elseif (tabPos == 0) then
                   theend = spacePos
                else
                   theend = MIN(spacePos, tabPos)
                endif
                theend=theend+start-1
                
                !print*, 'theend=',theend
              end do
              
    
              ! Extract the desired word
              word = trim(ADJUSTL(val(start:theend)))
              
              ! If a space or tab is stuck at the end of the word, remove it
              if (LEN_TRIM(word) > 0 .and. word(LEN_TRIM(word):LEN_TRIM(word)) == CHAR(9)) then
                word=word(1:LEN_TRIM(word)-1)
              end if
              
              
              if (LEN_TRIM(word) > 0 .and. word(LEN_TRIM(word):LEN_TRIM(word)) == 'd') then
                
                print*, 'missing "0 after d in parameter ', trim(keyword)
                
              else
               !print*, 'last word=', word
                read(word,*) keyval
                !print*, '    ', trim(keyword), '= ', keyval
                strfind=.true.
              end if
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
       
  end subroutine read_inputX
  
	
 
end module read_input_module


