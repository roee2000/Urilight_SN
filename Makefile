# ---- Compiler settings ---------
COMPILER=gfortran
OPTIMIZATION=-g -fcheck=mem -cpp -ffree-line-length-none -fno-range-check -fopenmp 
# -pg: add time profiling
LDFLAGS=-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -fopenmp

# Build directory
BUILD_DIR=build
EXE=$(BUILD_DIR)/UriLight.exe

# SLURM configuration
slurm_partition = cluster
slurm_account = rcl
slurm_ncores = 4

# PBS configuration (any free node, like opacity tables)
pbs_node = 1
pbs_ppn = 4
pbs_mem = 300gb
pbs_walltime = 24:00:00
pbs_jobname = urilight_run

# Create build directory if it doesn't exist
$(shell mkdir -p $(BUILD_DIR))

all: $(EXE)

# Object files with build directory prefix
OBJ= $(addprefix $(BUILD_DIR)/, globals.o logger.o physical_constants.o atomic_data_types.o read_nist_db.o arrays.o general_functions.o transport_general_functions.o radioactive_decay.o randomnumbers.o mesh.o gammatransfer.o uvoirtransfer.o diagnostics.o gamma_physics.o atomic_physics.o uvoir_physics.o gray_opacity.o)

$(EXE): $(OBJ) $(BUILD_DIR)/UriLight.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -o $(EXE) $(BUILD_DIR)/UriLight.o $(OBJ)

$(BUILD_DIR)/globals.o: src/globals.f90 
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/logger.o: src/logger.f90 $(BUILD_DIR)/globals.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/physical_constants.o: src/physical_constants.f90 
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/atomic_data_types.o: src/atomic_data_types.f90
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/read_nist_db.o: src/read_nist_db.f90 $(BUILD_DIR)/atomic_data_types.o $(BUILD_DIR)/physical_constants.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/arrays.o: src/arrays.f90 $(BUILD_DIR)/globals.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/general_functions.o: src/general_functions.f90 
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/randomnumbers.o: src/randomnumbers.f90 $(BUILD_DIR)/physical_constants.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/radioactive_decay.o: src/radioactive_decay.f90 $(BUILD_DIR)/randomnumbers.o $(BUILD_DIR)/physical_constants.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/atomic_physics.o: src/atomic_physics.f90 $(BUILD_DIR)/physical_constants.o $(BUILD_DIR)/general_functions.o $(BUILD_DIR)/globals.o $(BUILD_DIR)/arrays.o $(BUILD_DIR)/atomic_data_types.o $(BUILD_DIR)/read_nist_db.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/mesh.o: src/mesh.f90 $(BUILD_DIR)/general_functions.o $(BUILD_DIR)/randomnumbers.o $(BUILD_DIR)/physical_constants.o $(BUILD_DIR)/globals.o $(BUILD_DIR)/atomic_physics.o $(BUILD_DIR)/logger.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/transport_general_functions.o: src/transport_general_functions.f90 $(BUILD_DIR)/mesh.o $(BUILD_DIR)/general_functions.o $(BUILD_DIR)/randomnumbers.o $(BUILD_DIR)/physical_constants.o $(BUILD_DIR)/globals.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/diagnostics.o: src/diagnostics.f90 $(BUILD_DIR)/mesh.o $(BUILD_DIR)/general_functions.o $(BUILD_DIR)/arrays.o $(BUILD_DIR)/globals.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/gamma_physics.o: src/gamma_physics.f90 $(BUILD_DIR)/physical_constants.o $(BUILD_DIR)/randomnumbers.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/uvoir_physics.o: src/uvoir_physics.f90 $(BUILD_DIR)/physical_constants.o $(BUILD_DIR)/atomic_physics.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/gray_opacity.o: src/gray_opacity.f90 $(BUILD_DIR)/globals.o $(BUILD_DIR)/physical_constants.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/gammatransfer.o: src/gammatransfer.f90 $(BUILD_DIR)/globals.o $(BUILD_DIR)/arrays.o $(BUILD_DIR)/general_functions.o $(BUILD_DIR)/transport_general_functions.o $(BUILD_DIR)/radioactive_decay.o $(BUILD_DIR)/randomnumbers.o $(BUILD_DIR)/mesh.o $(BUILD_DIR)/diagnostics.o $(BUILD_DIR)/gamma_physics.o $(BUILD_DIR)/logger.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/uvoirtransfer.o: src/uvoirtransfer.f90 $(BUILD_DIR)/globals.o $(BUILD_DIR)/arrays.o $(BUILD_DIR)/general_functions.o $(BUILD_DIR)/transport_general_functions.o $(BUILD_DIR)/radioactive_decay.o $(BUILD_DIR)/randomnumbers.o $(BUILD_DIR)/mesh.o $(BUILD_DIR)/diagnostics.o $(BUILD_DIR)/atomic_physics.o $(BUILD_DIR)/physical_constants.o $(BUILD_DIR)/uvoir_physics.o $(BUILD_DIR)/gray_opacity.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/UriLight.o: src/UriLight.f90 $(OBJ)
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)

rclean:
	rm -fr output
	rm -f pbs-urilight.out pbs-urilight.err qsub_out

run:
	$(MAKE) rclean
	./$(EXE)

run-only:
	./$(EXE)

srun:
	rm -fr slurm-*.out slurm-*.err;
	sbatch -n 1 --cpus-per-task=$(slurm_ncores) \
		--account=$(slurm_account) \
		--partition=$(slurm_partition) \
		--job-name="urilight_$(PWD)" \
		--error=slurm-%j.err \
		--output=slurm-%j.out \
		--wrap="export OMP_NUM_THREADS=$(slurm_ncores) && make -C $(PWD) run" > sbatch_out ;\
	cat sbatch_out

prun:
	$(MAKE) rclean
	@echo '#!/bin/bash' > .urilight_run.pbs
	@echo '#PBS -N $(pbs_jobname)' >> .urilight_run.pbs
	@echo '#PBS -l nodes=$(pbs_node):ppn=$(pbs_ppn),mem=$(pbs_mem)' >> .urilight_run.pbs
	@echo '#PBS -l walltime=$(pbs_walltime)' >> .urilight_run.pbs
	@echo '#PBS -S /bin/bash' >> .urilight_run.pbs
	@echo '#PBS -j oe' >> .urilight_run.pbs
	@echo '#PBS -o /dev/null' >> .urilight_run.pbs
	@echo '#PBS -e /dev/null' >> .urilight_run.pbs
	@echo 'cd $(CURDIR)' >> .urilight_run.pbs
	@echo 'OUT="$(CURDIR)/pbs-urilight.out"' >> .urilight_run.pbs
	@echo 'ERR="$(CURDIR)/pbs-urilight.err"' >> .urilight_run.pbs
	@echo ': > "$$OUT"; : > "$$ERR"' >> .urilight_run.pbs
	@echo 'exec >> "$$OUT" 2>> "$$ERR"' >> .urilight_run.pbs
	@echo 'echo "=== job start $$(date) ==="' >> .urilight_run.pbs
	@echo 'module load gcc/gcc-fortran-12.2.0' >> .urilight_run.pbs
	@echo 'export OMP_NUM_THREADS=$(pbs_ppn)' >> .urilight_run.pbs
	@echo 'export OMP_PROC_BIND=spread' >> .urilight_run.pbs
	@echo 'export OMP_PLACES=cores' >> .urilight_run.pbs
	@echo 'export OMP_WAIT_POLICY=passive' >> .urilight_run.pbs
	@echo 'export GFORTRAN_UNBUFFERED_ALL=YES' >> .urilight_run.pbs
	@echo 'stdbuf -oL -eL make run-only' >> .urilight_run.pbs
	@echo 'echo "=== job end $$(date) ==="' >> .urilight_run.pbs
	qsub .urilight_run.pbs > qsub_out
	cat qsub_out

stop:
	@jobname="urilight_$(PWD)"; \
	jobids=$$(squeue -u "$$USER" -h -o "%A %j" | awk -v name="$$jobname" '$$2==name {print $$1}'); \
	if [ -n "$$jobids" ]; then \
		echo "Stopping job(s): $$jobids"; \
		echo "$$jobids" | xargs -r scancel; \
	else \
		echo "No matching running jobs found for $$jobname"; \
	fi

pstop:
	@jobids=$$(qstat -u "$$USER" 2>/dev/null | awk '$$4=="$(pbs_jobname)" {print $$1}'); \
	if [ -n "$$jobids" ]; then \
		echo "Stopping job(s): $$jobids"; \
		echo "$$jobids" | xargs qdel; \
	else \
		echo "No matching jobs for $(pbs_jobname)"; \
	fi

.PHONY: clean all rclean run run-only srun prun stop pstop

