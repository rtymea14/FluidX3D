# FluidX3D - Conjugate Heat Transfer Implementation

This document summarizes the recent implementation of **conjugate heat transfer** capabilities in FluidX3D, enabling simulation of heat transfer between different materials with varying thermal properties.

## Overview

FluidX3D now supports conjugate heat transfer for non-homogeneous media, allowing accurate thermal simulations across fluid-solid interfaces with different material properties. This implementation extends the existing TEMPERATURE extension with material-specific thermal properties.

## New Features

### Enhanced Thermal Properties
- **Material-specific heat capacity**: `rhocP[n]` - Heat capacity (ρ×cp) for each cell
- **Material-specific thermal conductivity**: `k[n]` - Thermal conductivity for each cell  
- **Material type classification**: `matType[n]` - Distinguishes between fluid (0) and solid (1) materials
- **Local thermal diffusivity**: Automatically calculated as `α = k/(ρ×cp)`

### Implementation Details

#### Core Changes
1. **Extended Memory Allocation**: Additional 5 bytes per cell for conjugate heat transfer properties
2. **Material-Aware Thermal Solver**: Modified temperature kernels to handle varying material properties
3. **Interface Heat Transfer**: Proper handling of heat conduction across fluid-solid boundaries

#### Key Properties
```cpp
lbm.k[n] = 1.0f;      // Thermal conductivity (fluid: 1.0, solid: 5.0)
lbm.rhocP[n] = 1.0f;  // Heat capacity (fluid: 1.0, solid: 2.0)
lbm.matType[n] = 0;   // Material type (0=fluid, 1=solid)
```

## Usage Example

### Basic Setup
```cpp
void main_setup() { // Requires TEMPERATURE extension in defines.hpp
    LBM lbm(64u, 64u, 64u, 1u, 1u, 1u, 0.02f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    
    parallel_for(lbm.get_N(), [&](ulong n) {
        uint x, y, z;
        lbm.coordinates(n, x, y, z);
        
        // Initialize default fluid properties
        lbm.k[n] = 1.0f;      // Thermal conductivity of fluid
        lbm.rhocP[n] = 1.0f;  // Heat capacity of fluid
        lbm.matType[n] = 0;   // 0 = fluid
        lbm.T[n] = 1.0f;      // Initial temperature
        
        // Define solid regions with different thermal properties
        if(/* solid region condition */) {
            lbm.flags[n] = TYPE_S;  // Solid boundary
            lbm.k[n] = 5.0f;        // Higher thermal conductivity
            lbm.rhocP[n] = 2.0f;    // Different heat capacity
            lbm.matType[n] = 1;     // 1 = solid
            lbm.T[n] = 2.0f;        // Higher initial temperature
        }
        
        // Set temperature boundaries
        if(/* boundary condition */) {
            lbm.flags[n] = TYPE_T;  // Temperature boundary
            lbm.T[n] = /* boundary temperature */;
        }
    });
}
```

### Example Simulation
The implementation includes a complete example demonstrating:
- **Heterogeneous geometry**: Central solid block surrounded by fluid
- **Material contrast**: 5× higher thermal conductivity and 2× heat capacity in solid
- **Thermal boundaries**: Hot top (T=1.5) and cold bottom (T=0.5) surfaces
- **Heat transfer analysis**: Periodic temperature field export for analysis

## Technical Implementation

### Thermal Diffusion Solver
The conjugate heat transfer solver uses:
- **Local thermal diffusivity**: `α_local = k[n] / rhocP[n]`
- **Heat capacity ratio**: `σ = rhocP[n] / rho_cp_ref`
- **Adaptive relaxation**: `w_T_local = 1/(2×α_local + 0.5)`

### Kernel Modifications
- Enhanced `stream_collide()` kernel with material-aware thermal calculations
- Automatic interface treatment for seamless heat transfer across material boundaries
- Gradient calculations accounting for varying thermal properties

## Requirements

### Required Extensions
- `#define TEMPERATURE` in `defines.hpp`

### Memory Usage
- Additional **5 bytes per cell** for conjugate heat transfer properties:
  - `rhocP`: 4 bytes (float)
  - `matType`: 1 byte (uchar)
  - `k`: Already part of TEMPERATURE extension

## Applications

This conjugate heat transfer implementation enables simulation of:
- **Heat exchangers** with fluid-solid thermal interactions
- **Electronics cooling** with chip-fluid interfaces  
- **Building thermal analysis** with wall-air heat transfer
- **Geothermal systems** with rock-fluid thermal coupling
- **Manufacturing processes** involving thermal treatment

## Performance

The conjugate heat transfer implementation maintains FluidX3D's high performance characteristics:
- Minimal memory overhead (5 bytes/cell)
- GPU-optimized thermal solver
- Compatible with all existing FluidX3D features (multi-GPU, raytracing, etc.)

## Getting Started

1. Enable the TEMPERATURE extension in `defines.hpp`
2. Use the provided conjugate heat transfer example in `setup.cpp` as a starting point
3. Modify material properties and geometry according to your specific application
4. Run simulation and analyze temperature field evolution

---

**Note**: This readme.md file specifically documents the conjugate heat transfer implementation. For complete FluidX3D documentation, performance benchmarks, and general usage instructions, refer to:
- [README.md](README.md) - Main project documentation with comprehensive feature overview
- [DOCUMENTATION.md](DOCUMENTATION.md) - Detailed usage guide and setup instructions