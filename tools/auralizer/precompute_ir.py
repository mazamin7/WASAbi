import numpy as np
import json
import os
import sys

def deconvolve(recorded, source, eps=1e-5):
    """
    Deconvolve recorded signal with source signal in frequency domain.
    """
    N = len(recorded)
    # Ensure power of 2 for better performance
    N_fft = 2**int(np.ceil(np.log2(N)))
    
    REC = np.fft.fft(recorded, n=N_fft)
    SRC = np.fft.fft(source, n=N_fft)
    
    # Regularized division
    IR_freq = REC / (SRC + eps)
    ir = np.fft.ifft(IR_freq).real
    
    return ir[:N]

def main():
    if len(sys.argv) < 2:
        print("Usage: python precompute_ir.py <experiment_dir> [z_height_meters]")
        return

    exp_dir = sys.argv[1]
    z_target = float(sys.argv[2]) if len(sys.argv) > 2 else 1.5
    
    config_path = os.path.join(exp_dir, "config.json")
    asset_path = os.path.join(exp_dir, "asset.json")
    record_path = os.path.join(exp_dir, "output", "record_data_0.bin")
    source_path = os.path.join(exp_dir, "output", "source_data_0.bin")

    if not all(os.path.exists(p) for p in [config_path, asset_path, record_path, source_path]):
        print("Error: Missing required files in experiment directory.")
        return

    with open(config_path, 'r') as f:
        config = json.load(f)
    
    with open(asset_path, 'r') as f:
        asset = json.load(f)

    dh = config['simulation'].get('dh', 0.2)
    dt = config['simulation'].get('dt', 0.0002)
    
    # Determine grid dimensions from partitions
    rxs, rys, rzs = 1e9, 1e9, 1e9
    rxe, rye, rze = -1e9, -1e9, -1e9
    
    for p in asset['partitions']:
        rxs = min(rxs, p['x'])
        rys = min(rys, p['y'])
        rzs = min(rzs, p['z'])
        rxe = max(rxe, p['x'] + p['w'])
        rye = max(rye, p['y'] + p['h'])
        rze = max(rze, p['z'] + p['d'])
        
    gs_x = int(round((rxe - rxs) / dh))
    gs_y = int(round((rye - rys) / dh))
    gs_z = int(round((rze - rzs) / dh))
    
    z_idx = int(round((z_target - rzs) / dh))
    z_idx = max(0, min(gs_z - 1, z_idx))
    actual_z = rzs + z_idx * dh
    
    print(f"Grid: {gs_x}x{gs_y}x{gs_z}")
    print(f"Target Z: {z_target}m -> Actual Z: {actual_z}m (idx: {z_idx})")

    # Load source
    source_signal = np.fromfile(source_path, dtype=np.float64)
    n_steps = len(source_signal)
    
    # record_data_0.bin is time-major: [time_step][z][y][x]
    # We want to extract a slice [time_step][z_idx][y][x]
    
    points_per_frame = gs_x * gs_y * gs_z
    slice_size = gs_x * gs_y
    
    # We'll read frame by frame to extract the specific Z slice
    # To be memory efficient we use memmap
    field_data = np.memmap(record_path, dtype=np.float64, mode='r', shape=(n_steps, gs_z, gs_y, gs_x))
    
    slice_data = field_data[:, z_idx, :, :] # shape: (n_steps, gs_y, gs_x)
    
    ir_len = int(0.5 / dt) # 0.5 seconds
    ir_len = min(ir_len, n_steps)
    
    print(f"Generating IRs for {gs_x}x{gs_y} points...")
    
    ir_grid = np.zeros((gs_y, gs_x, ir_len), dtype=np.float32)
    
    for y in range(gs_y):
        for x in range(gs_x):
            recorded = slice_data[:, y, x]
            # Handle potential silence or boundaries
            if np.max(np.abs(recorded)) < 1e-15:
                continue
            
            ir = deconvolve(recorded, source_signal)
            ir_grid[y, x, :] = ir[:ir_len].astype(np.float32)
            
        if y % 10 == 0:
            print(f"Progress: {y}/{gs_y}")

    output_grid_path = os.path.join(exp_dir, "output", "ir_grid_2d.bin")
    output_meta_path = os.path.join(exp_dir, "output", "ir_metadata.json")
    
    ir_grid.tofile(output_grid_path)
    
    meta = {
        "nx": gs_x,
        "ny": gs_y,
        "dh": dh,
        "dt": dt,
        "z": actual_z,
        "ir_len": ir_len,
        "rxs": rxs,
        "rys": rys
    }
    
    with open(output_meta_path, 'w') as f:
        json.dump(meta, f, indent=4)
        
    print(f"Done! Saved to {output_grid_path}")

if __name__ == "__main__":
    main()
