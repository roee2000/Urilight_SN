# ---- Compiler settings ---------
COMPILER=gfortran
OPTIMIZATION=-g -fcheck=bounds -fcheck=mem -cpp -ffree-line-length-none -fno-range-check
LDFLAGS=-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk

# Build directory
BUILD_DIR=build
EXE=$(BUILD_DIR)/UriLight.exe

# Create build directory if it doesn't exist
$(shell mkdir -p $(BUILD_DIR))

all: $(EXE)

# Object files with build directory prefix
OBJ= $(addprefix $(BUILD_DIR)/, globals.o physical_constants.o arrays.o general_functions.o transport_general_functions.o radioactive_decay.o randomnumbers.o mesh.o gammatransfer.o uvoirtransfer.o diagnostics.o gamma_physics.o atomic_physics.o uvoir_physics.o)

$(EXE): $(OBJ) $(BUILD_DIR)/UriLight.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -o $(EXE) $(BUILD_DIR)/UriLight.o $(OBJ)

$(BUILD_DIR)/globals.o: src/globals.f90 
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

$(BUILD_DIR)/mesh.o: src/mesh.f90 $(BUILD_DIR)/general_functions.o $(BUILD_DIR)/randomnumbers.o $(BUILD_DIR)/physical_constants.o $(BUILD_DIR)/globals.o $(BUILD_DIR)/atomic_physics.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/transport_general_functions.o: src/transport_general_functions.f90 $(BUILD_DIR)/mesh.o $(BUILD_DIR)/general_functions.o $(BUILD_DIR)/randomnumbers.o $(BUILD_DIR)/physical_constants.o $(BUILD_DIR)/globals.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/diagnostics.o: src/diagnostics.f90 $(BUILD_DIR)/mesh.o $(BUILD_DIR)/general_functions.o $(BUILD_DIR)/arrays.o $(BUILD_DIR)/globals.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/gamma_physics.o: src/gamma_physics.f90 $(BUILD_DIR)/physical_constants.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/uvoir_physics.o: src/uvoir_physics.f90 $(BUILD_DIR)/physical_constants.o $(BUILD_DIR)/atomic_physics.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/gammatransfer.o: src/gammatransfer.f90 $(BUILD_DIR)/globals.o $(BUILD_DIR)/arrays.o $(BUILD_DIR)/general_functions.o $(BUILD_DIR)/transport_general_functions.o $(BUILD_DIR)/radioactive_decay.o $(BUILD_DIR)/randomnumbers.o $(BUILD_DIR)/mesh.o $(BUILD_DIR)/diagnostics.o $(BUILD_DIR)/gamma_physics.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/uvoirtransfer.o: src/uvoirtransfer.f90 $(BUILD_DIR)/globals.o $(BUILD_DIR)/arrays.o $(BUILD_DIR)/general_functions.o $(BUILD_DIR)/transport_general_functions.o $(BUILD_DIR)/radioactive_decay.o $(BUILD_DIR)/randomnumbers.o $(BUILD_DIR)/mesh.o $(BUILD_DIR)/diagnostics.o $(BUILD_DIR)/atomic_physics.o $(BUILD_DIR)/physical_constants.o $(BUILD_DIR)/uvoir_physics.o
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/UriLight.o: src/UriLight.f90 $(OBJ)
	$(COMPILER) $(OPTIMIZATION) $(LDFLAGS) -J$(BUILD_DIR) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)

rclean:
	rm -fr fort.* out profile* totals spectrum_uvoir luminocity ejecta_profile run.e* run.o* magnitude opacity_cell_*

.PHONY: clean all rclean

