METAL_FLAGS=
all: shaders
debug: METAL_FLAGS += -fdebug-info-for-profiling
debug: shaders
shaders: 
	xcrun -sdk macosx metal $(METAL_FLAGS) -c src/metal/shaders.metal -o src/metal/shaders.air
	xcrun -sdk macosx metallib src/metal/shaders.air -o src/metal/shaders.metallib