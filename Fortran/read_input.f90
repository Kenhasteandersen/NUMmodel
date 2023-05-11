module read_input_module
        use global
        implicit none

        contains
        
                subroutine read_input(filename)
                        character(*), intent(in) :: filename
                        character(len=100) :: cbuf, c1, keyword, val
                        integer:: ios, n, ilength
                        real(MK):: ii
                        OPEN(1, file = filename)
                        do
                        read(1, '(A)', iostat=ios) cbuf
                        if (ios /= "!") exit
                        c1=cbuf(1:1)
                        if ( c1 /= "!" ) then
                        ilength=index(cbuf,' ')-1
                        n=index(cbuf,'=')
                        keyword=cbuf(1:n-1)
                        val=cbuf(n+1:ilength)
                        if (keyword == 'diff') then
                                read(val,*), diff
                        else if (keyword == 'Nx') then
                                read(val,*), Nx
                        endif
                        endif
                        end do
                        CLOSE(1)
                end subroutine read_input

end module read_input_module


