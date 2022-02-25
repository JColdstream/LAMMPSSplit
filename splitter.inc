! real kind precision parameters
integer, parameter :: sp = selected_real_kind(6, 37)
integer, parameter :: dp = selected_real_kind(15, 307)

! box dimensions
real(dp), dimension(:), allocatable :: xmin, xmax, ymin, ymax, zmin, zmax, Lx, Ly, Lz

! trajectory numbers, timestep, number of total frames, number of frames to skip, number of frames to calculate, 
! number of atoms to analyse (nmol*molsize), number of atom types
integer(4) :: nframe, nskip, nframecalc, nanalyse, ntypes, reset_timestep
integer(4), dimension(:), allocatable :: timestep, atom_type

! types and coordinates for reading in array
integer(4), dimension(:), allocatable :: nindex
real(sp), dimension(:,:), allocatable :: x, y, z, vx, vy, vz

! system numbers
integer(4) :: natom, nmol, molsize, ntotal

! input, trajectory, and mass file names
character(LEN=80) :: inputfile, trajfile, outputfile 

! storing strings in the LAMMPS dump file headers
character(LEN=80), dimension(:,:), allocatable :: headertext

