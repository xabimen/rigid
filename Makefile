CC=gfortran

main : angles.f90 check_min.f90 io.f90 pair_dist.f90 read_write.f90 main_full_traj.f90
	$(CC) angles.f90 check_min.f90 io.f90 pair_dist.f90 read_write.f90 main_full_traj.f90 -o main 

main_lammps : angles.f90 check_min.f90 io.f90 pair_dist.f90 read_write.f90 main_lammps.f90
	$(CC) angles.f90 check_min.f90 io.f90 pair_dist.f90 read_write.f90 main_lammps.f90 -o main 
	
clean:
	rm -f *.dat
	rm -f *.mod
	rm -f main
	rm -f main_lammps