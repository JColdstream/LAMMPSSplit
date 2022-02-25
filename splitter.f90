program splitter

implicit none
include 'splitter.inc'

integer(sp) :: i

call start

call skipframe

write(*,*) 'Reading remaining', nframe, 'frames.'

do i = 1, nframe
  call readheader(i)
  call readcoordinates(i)
enddo

do i=1,nframe
  call writecoordinates(i)
enddo

write(*,*) 'Dump file written to ', outputfile

call avgboxlength

contains

! calls readinput, allocates relevant arrays
subroutine start
integer :: i

call get_command_argument(1,inputfile)
if (len_trim(inputfile) == 0) then
  write(*,*) 'Provide an input file.'
  call exit
endif

call readinput

nframe = nframe-nskip

! allocate arrays
allocate(nindex(ntotal))
allocate(atom_type(nanalyse))
allocate(x(nframe,nanalyse))
allocate(y(nframe,nanalyse))
allocate(z(nframe,nanalyse))
allocate(vx(nframe,nanalyse))
allocate(vy(nframe,nanalyse))
allocate(vz(nframe,nanalyse))
allocate(xmin(nframe))
allocate(xmax(nframe))
allocate(ymin(nframe))
allocate(ymax(nframe))
allocate(zmin(nframe))
allocate(zmax(nframe))
allocate(Lx(nframe))
allocate(Ly(nframe))
allocate(Lz(nframe))
allocate(timestep(nframe))
allocate(headertext(nframe,4))

! open lammps trajectory file
open(11, file=trajfile, status='old')
write(*,*) 'Traectory file opened.'
open(21, file=outputfile, status='new')

end subroutine start

! reads analysis.input input file
subroutine readinput
open(10, file=inputfile,status='old')
read(10,*)
read(10,*)
read(10,*) trajfile
read(10,*) outputfile
read(10,*) nframe
read(10,*) nskip
read(10,*) nmol
read(10,*) molsize
read(10,*) ntotal
read(10,*) reset_timestep
write(*,*) trajfile

nanalyse = nmol * molsize

nframecalc = nframe - nskip

end subroutine readinput

! skips reading frames that don't need to be analysed
subroutine skipframe
integer :: i
do i=1, nskip*(ntotal+9)
  read(11,*)
enddo
write(*,*) "Skipping ", nskip, 'frames.'
end subroutine skipframe

! reads in box length, timestep and natoms from the LAMMPS trajectory headers
subroutine readheader(i)

integer :: i

read(11,*) headertext(i,1)  
read(11,*) timestep(i)
read(11,*) headertext(i,2)
read(11,*) natom
read(11,*) headertext(i,3)
read(11,*) xmin(i), xmax(i)
read(11,*) ymin(i), ymax(i)
read(11,*) zmin(i), zmax(i)
read(11,*) headertext(i,4)

!if(natom  .ne. ntotal) stop 'mismatched number of atoms'

Lx(i) = xmax(i) - xmin(i)
Ly(i) = ymax(i) - ymin(i)
Lz(i) = zmax(i) - zmin(i)

end subroutine readheader

! reads in the atom types and coordinates from the LAMMPS trajectory file
subroutine readcoordinates(iframe)
  integer(sp) :: i, iframe
  do i = 1, nanalyse 
    read(11,*) nindex(i), atom_type(i), x(iframe,i), y(iframe,i), z(iframe,i), vx(iframe,i), vy(iframe,i), vz(iframe,i)
  enddo
  do i = 1, ntotal-nanalyse
    read(11,*)
  enddo
end subroutine readcoordinates

! writes frame iframe, to a file.
subroutine writecoordinates(iframe)
  integer :: iframe, i
  if (reset_timestep .eq. 1) timestep = timestep - timestep(1)
  if (reset_timestep .gt. 1) write(*,*) 'Use 0 or 1 as a value for resetting the timestep.'
  open(21, file=outputfile, status='old')
  write(21,*) 'ITEM: TIMESTEP'
  write(21,*) timestep(iframe) 
  write(21,*) 'ITEM: NUMBER OF ATOMS'
  write(21,*) nanalyse
  write(21,*) 'ITEM: BOX BOUNDS pp pp pp' 
  write(21,*) xmin(iframe), xmax(iframe)
  write(21,*) ymin(iframe), ymax(iframe)
  write(21,*) zmin(iframe), zmax(iframe)
  write(21,*) 'ITEM: ATOMS id type x y z vx vy vz'
  do i=1, nanalyse
    write(21,*) nindex(i), atom_type(i), x(iframe,i), y(iframe,i), z(iframe,i), vx(iframe,i), vy(iframe,i), vz(iframe,i)
  enddo
end subroutine writecoordinates

subroutine avgboxlength

  real(sp) :: avglx, avgly, avglz
  avglx = 0
  avgly = 0
  avglz = 0
  do i=1, nframe
    avglx = avglx + Lx(i)
    avgly = avgly + Ly(i)
    avglz = avglz + Lz(i)
  enddo
  avglx = avglx / real(nframe)
  avgly = avgly / real(nframe)
  avglz = avglz / real(nframe)
  write(*,*) 'Average length x :', avglx
  write(*,*) 'Average length y :', avgly
  write(*,*) 'Average length z :', avglz

end subroutine avgboxlength

end program splitter
