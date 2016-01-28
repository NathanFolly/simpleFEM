subroutine readBC(bcfile, bc, nbc)
	use mytypes
	implicit none




	character*20, intent(in) :: bcfile
	type(bcobject), allocatable, intent(inout):: bc(:)
	integer, intent(inout) :: nbc
	integer :: filestatus, i
	character*100 :: line

	open(unit=33, file=bcfile, status='old', iostat=filestatus)

	!look for the $start tag
	rewind(33)

100	read(33,'(A)') line
	if (filestatus.ne.0) then
		print *,'Error reading boundary condition file with the specified name', bcfile
	else if(index(line,'$start').ne.0) then
101		read(33, '(A)') line
		if (index(line,'!') .ne. 0) then
			goto 101
		else
			read(line, *) nbc
		end if


		 allocate(bc(nbc))

		do i=1,nbc
102			read(33,'(A)') line
			if (index(line,'!').ne.0) then
				goto 102
			else if(index(line,'$end') .ne. 0 ) then
				write(*,*) 'Error: The boudary condition file ended before the specified number of boundary conditions could be read.'
			else
				read(line, *) bc(i)%boundaryname, bc(i)%bcnature, bc(i)%conditions(:)
			end if
		end do
	else
		goto 100
	end if

close(33)






end subroutine readBC