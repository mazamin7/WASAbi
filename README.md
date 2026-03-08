<p align="center">
  <img src="logo.png">
</p>

# WASAbi - a Wave-based Acoustic Simulator using ARD

C++ implementation of Adaptive Rectangular Decomposition (ARD) with frequency dependent atmospheric absorption in 2.5D (3D with constant height).

Theory:
> Gerardo Cicalese, Gabriele Ciaramella, Ilario Mazzieri; Addressing atmospheric absorption in adaptive rectangular decomposition. J. Acoust. Soc. Am. 1 October 2024; 156 (4): 2328–2339. [https://doi.org/10.1121/10.0030468](https://doi.org/10.1121/10.0030468)

Extended from ARD-simulator by [@jinnsjj](https://github.com/jinnsjj).
> [https://github.com/jinnsjj/ARD-simulator](https://github.com/jinnsjj/ARD-simulator)

## Experiment Workflow
WASAbi is structured around **Experiments**. Each experiment is a self-contained directory within `source/experiments/` containing all necessary configurations and definitions.

### 📁 Directory Structure
```text
experiments/[experiment_name]/
├── config.json       # Simulation & Visualization parameters
├── asset.json        # Geometry (Partitions), Sources, and Recorders
└── output/           # Binary simulation results (.bin)
```

### ⚙️ Configuration (`config.json`)
Defines the physics and technical parameters of the simulation.
```json
{
    "simulation": {
        "duration": 0.2,            // Total simulation time (s)
        "dh": 0.5,                 // Spatial step (m)
        "dt": 0.000625,            // Temporal step (s)
        "c0": 343.5,               // Speed of sound (m/s)
        "boundary_absorption": 1.0, // Global wall absorption [0, 1]
        "n_pml_layers": 5          // Absorbing boundary layers
    },
    "visualization": {
        "viz_skip": 10,            // Steps to skip between renders
        "max_viz_gain": 100        // Visual intensity scaling
    }
}
```

### 🧱 Assets (`asset.json`)
Defines the spatial layout and transducer positions.
```json
{
    "partitions": [
        { "x": 0, "y": 0, "z": 0, "w": 30, "h": 10, "d": 10 }
    ],
    "sources": [
        { "type": "gaussian", "x": 15, "y": 25, "z": 3 }
    ],
    "recorders": [
        { "x": 15, "y": 15, "z": 3 }
    ]
}
```

---

## Running the Simulator

Launch the simulator using the provided batch scripts or directly from the command line.

### 🕹️ CLI Modes (`--mode`)
- **`sim-viz`**: Standard real-time visualization mode.
- **`sim-record-field`**: Records the 3D pressure field over time to `record_data_*.bin`.
- **`sim-record-response`**: Records the pressure at recorder positions to `response_data_*.bin`.
- **`viz-record`**: Benchmarking mode for visualization performance.

### 🚀 Commands
**Using the batch scripts (easiest):**
```cmd
.\run_sim_viz.bat hall            # Run 'hall' experiment in viz mode
.\run_sim_record_field.bat room   # Record 3D field for 'room'
```

**Direct Execution:**
```cmd
.\WASAbiApp.exe --experiment hall --mode sim-record-response
```

---

## Post-Processing (MATLAB)
Tools are provided in `postprocessing/` to analyze the binary outputs:
- **`visualize_field.m`**: Standard 3D field visualizer. Automatically orchestrates `config.json` and `asset.json` to decode binary field data.
- **`rir.m`**: Computes Energy Decay Curves (EDC) and RIR statistics.

---

## Building WASAbi

### 🟢 Prerequisites
- **NVIDIA GPU** (Compute Capability 6.1+).
- **CUDA Toolkit** (Tested with v12.4).
- **Visual Studio 2022** with MSVC.

### 🪟 Windows Build
Execute the automated build script from the project root:
```cmd
.\build_cuda.bat
```
This will initialize the environment, configure CMake, and build the `WASAbiApp.exe` executable into `source/build/`.
