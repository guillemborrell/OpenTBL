sub-omp:
	make -f makefile-OMP-bg ITPC	
	qsub -n 128 -t 10 --mode smp -q prod-devel -O out-128 ITPC

sub-omp-32:
	make -f makefile-OMP-bg ITPC	
	qsub -n 32 -t 10 --mode smp -q prod-devel -O out-32 ITPC

sub-omp-64:
	make -f makefile-OMP-bg ITPC	
	qsub -n 64 -t 12 --mode smp -q prod-devel -O out-64 ITPC

sub-omp-128:
	make -f makefile-OMP-bg ITPC
	rename ITPC ITPC-128 ITPC	
	qsub -n 128 -t 15 --mode smp -q prod-devel -O out-128 ITPC-128

sub-omp-256:
	make -f makefile-OMP-bg ITPC
	rename ITPC ITPC-256 ITPC	
	qsub -n 256 -t 40 --mode smp -q prod-devel -O out-256 ITPC-256

sub-omp-512:
	make -f makefile-OMP-bg ITPC
	rename ITPC ITPC-512 ITPC	
	qsub -n 512 -t 150 --mode smp -q prod -O out-512 ITPC-512

sub-omp-1024:
	make -f makefile-OMP-bg ITPC
	rename ITPC ITPC-1024 ITPC	
	qsub -n 1024 -t 20 --mode smp -q prod -O out-1024 ITPC-1024	

sub-omp-2048:
	make -f makefile-OMP-bg ITPC
	rename ITPC ITPC-2048-2bls ITPC	
	qsub -n 2048 -t 720 --mode smp -q prod -O out-2048-2bls ITPC-2048-2bls	

sub-omp-4096:
	make -f makefile-OMP-bg ITPC
	rename ITPC ITPC-4096 ITPC	
	qsub -n 4096 -t 720 --mode smp -q prod -O out-4096 ITPC-4096

sub-omp-8192:
	make -f makefile-OMP-bg ITPC
	rename ITPC ITPC-8192 ITPC	
	qsub -n 8192 -t 720 --mode smp -q prod -O out-8192 ITPC-8192

compilar-omp:
	make -f makefile-OMP-bg ITPC
compilar-mpi:
	make -f makefile-mpi ITPC
subv:
	make -f makefile-OMP-vul
	qsub qfile-OMP
# sub-omp:
# 	make -f makefile-OMP-bg ITPC
# 	export np=128
#         export tp=10
# 	qsub -n $np -t $tp --mode smp -q prod-devep -O resultado-$np ITPC

clean: 
	find . \( -name '*.o' \) -exec rm -rf {} \;
	find . \( -name '*.mod' \) -exec rm -rf {} \;
	find . \( -name '*~' \) -exec rm -rf {} \;
