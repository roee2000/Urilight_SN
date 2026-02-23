# ---- Compiler settings ---------
COMPILER=gfortran
OPTIMIZATION=-g -pg -fcheck=bounds -fcheck=mem -cpp -ffree-line-length-none -fno-range-check -fopenmp
LDFLAGS=-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -fopenmp

# Build directory
BUILD_DIR=build
EXE=$(BUILD_DIR)/UriLight.exe

# SLURM configuration
slurm_partition = cluster
slurm_account = rcl
slurm_ncores = 32

# Create build directory if it doesn't exist
$(shell mkdir -p $(BUILD_DIR))

all: $(EXE)

# Object files with build directory prefix
OBJ= $(addprefix $(BUILD_DIR)/, globals.o logger.o physical_constants.o arrays.o general_functions.o transport_general_functions.o radioactive_decay.o randomnumbers.o mesh.o gammatransfer.o uvoirtransfer.o diagnostics.o gamma_physics.o atomic_physics.o uvoir_physics.o gray_opacity.o)

$(EXE): $(OBJ) $(BUILD_DIR)/UriLight.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -o $(EXE) $(BUILD_DIR)/UriLight.o $(OBJ)

$(BUILD_DIR)/globals.o: src/globals.f90 
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/logger.o: src/logger.f90 $(BUILD_DIR)/globals.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/physical_constants.o: src/physical_constants.f90 
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/arrays.o: src/arrays.f90 $(BUILD_DIR)/globals.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/general_functions.o: src/general_functions.f90 
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/randomnumbers.o: src/randomnumbers.f90 $(BUILD_DIR)/physical_constants.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/radioactive_decay.o: src/radioactive_decay.f90 $(BUILD_DIR)/randomnumbers.o $(BUILD_DIR)/physical_constants.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/atomic_physics.o: src/atomic_physics.f90 $(BUILD_DIR)/physical_constants.o $(BUILD_DIR)/general_functions.o $(BUILD_DIR)/globals.o $(BUILD_DIR)/arrays.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/mesh.o: src/mesh.f90 $(BUILD_DIR)/general_functions.o $(BUILD_DIR)/randomnumbers.o $(BUILD_DIR)/physical_constants.o $(BUILD_DIR)/globals.o $(BUILD_DIR)/atomic_physics.o $(BUILD_DIR)/logger.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/transport_general_functions.o: src/transport_general_functions.f90 $(BUILD_DIR)/mesh.o $(BUILD_DIR)/general_functions.o $(BUILD_DIR)/randomnumbers.o $(BUILD_DIR)/physical_constants.o $(BUILD_DIR)/globals.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/diagnostics.o: src/diagnostics.f90 $(BUILD_DIR)/mesh.o $(BUILD_DIR)/general_functions.o $(BUILD_DIR)/arrays.o $(BUILD_DIR)/globals.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/gamma_physics.o: src/gamma_physics.f90 $(BUILD_DIR)/physical_constants.o
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
	@outdir=$$(awk -F"'" '/output_dir/{print $$2}' data_file | tail -1); \
	if [ -z "$$outdir" ]; then outdir=output; fi; \
	rm -fr fort.* run.e* run.o* sbatch_out "$$outdir"

run:
	$(MAKE) rclean;
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

stop:
	@jobname="urilight_$(PWD)"; \
	jobids=$$(squeue -u "$$USER" -h -o "%A %j" | awk -v name="$$jobname" '$$2==name {print $$1}'); \
	if [ -n "$$jobids" ]; then \
		echo "Stopping job(s): $$jobids"; \
		echo "$$jobids" | xargs -r scancel; \
	else \
		echo "No matching running jobs found for $$jobname"; \
	fi

.PHONY: clean all rclean run srun stop

