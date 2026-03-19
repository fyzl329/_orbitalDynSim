"""
Example: Run all integrators with unified driver and compare conserved quantities.

This demonstrates:
1. Unified interface: run_simulation(method="...", dt=..., steps=...)
2. Automatic energy + angular momentum tracking
3. Error analysis for conserved quantities
4. Data export using pandas DataFrames
"""
import sys
import numpy as np
from pathlib import Path

# Add simulators to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from simulators import run_simulation, SimulationResults


def main():
    """Run comparison of all integrators under identical conditions."""
    
    print("=" * 70)
    print("ORBITAL DYNAMICS SIMULATOR - UNIFIED DRIVER DEMO")
    print("=" * 70)
    print()
    
    # Identical initial conditions for all methods
    initial_position = np.array([1.0, 0.0])      # Start at r=1
    initial_velocity = np.array([0.0, 1.0])      # Circular orbit velocity
    dt = 0.001                                   # Time step
    steps = 100000                               # Run for 100k steps
    
    methods = ["euler", "verlet", "rk4"]
    
    print(f"Simulation Parameters:")
    print(f"  Initial Position: {initial_position}")
    print(f"  Initial Velocity: {initial_velocity}")
    print(f"  Time Step (dt):   {dt}")
    print(f"  Total Steps:      {steps}")
    print(f"  Total Time:       {steps * dt:.1f}")
    print()
    print("=" * 70)
    
    # Initialize results container with pandas support
    results = SimulationResults()
    
    for method in methods:
        print(f"\nRunning {method.upper()} integrator...")
        
        times, positions, velocities, tracking = run_simulation(
            method=method,
            initial_position=initial_position,
            initial_velocity=initial_velocity,
            dt=dt,
            steps=steps,
        )
        
        # Add to pandas-based results container
        results.add_simulation(method, times, positions, velocities, tracking)
        
        # Extract conserved quantities
        energy = tracking['energy']
        energy_error = tracking['energy_error']
        l_momentum = tracking['angular_momentum']
        l_error = tracking['angular_momentum_error']
        
        # Compute orbit statistics
        distances = np.linalg.norm(positions, axis=1)
        r_min, r_max = distances.min(), distances.max()
        
        # Print results
        print(f"\n  Orbit Statistics:")
        print(f"    r_min:        {r_min:.6f}")
        print(f"    r_max:        {r_max:.6f}")
        print(f"    r_mean:       {distances.mean():.6f}")
        print(f"    Eccentricity: {(r_max - r_min) / (r_max + r_min):.6f}")
        
        print(f"\n  Energy Conservation:")
        print(f"    Initial E:    {energy[0]:.10f}")
        print(f"    Final E:      {energy[-1]:.10f}")
        print(f"    Max Error:    {energy_error.max():.2e}")
        print(f"    Mean Error:   {energy_error.mean():.2e}")
        
        print(f"\n  Angular Momentum Conservation:")
        print(f"    Initial L:    {l_momentum[0]:.10f}")
        print(f"    Initial L:    {l_momentum[0]:.10f}")
        print(f"    Final L:      {l_momentum[-1]:.10f}")
        print(f"    Max Error:    {l_error.max():.2e}")
        print(f"    Mean Error:   {l_error.mean():.2e}")
    
    print("\n" + "=" * 70)
    print("COMPARISON SUMMARY (pandas DataFrame)")
    print("=" * 70)
    print(results.comparison_table())
    
    print("\n" + "=" * 70)
    print("Exporting results to CSV...")
    print("=" * 70)
    
    # Export summary statistics
    analysis_dir = Path(__file__).parent / "analysis" / "data"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    
    summary_csv = analysis_dir / "comparison_summary.csv"
    results.export_summary_csv(summary_csv)
    
    # Export full time-series data
    results.export_all_csv(analysis_dir)
    
    print("\n" + "=" * 70)
    print("Simulation complete. All methods ran under identical conditions.")
    print("=" * 70)


if __name__ == "__main__":
    main()
