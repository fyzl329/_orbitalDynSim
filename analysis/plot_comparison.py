"""
Visualization: Plot trajectories and conserved quantities comparison.

Generates:
- Orbital trajectories for all methods
- Energy error comparison
- Angular momentum error comparison
- CSV exports of all results using pandas
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from simulators import run_simulation, SimulationResults


def plot_comparison():
    """Generate separate comparison plots for each integrator."""
    
    # Run simulations
    initial_position = np.array([1.0, 0.0])
    initial_velocity = np.array([0.0, 1.0])
    dt = 0.001
    steps = 50000
    
    # Use pandas-based results container
    results_pd = SimulationResults()
    results = {}
    
    for method in ["euler", "verlet", "rk4"]:
        print(f"Simulating {method}...")
        times, positions, velocities, tracking = run_simulation(
            method=method,
            initial_position=initial_position,
            initial_velocity=initial_velocity,
            dt=dt,
            steps=steps,
        )
        results[method] = {
            'times': times,
            'positions': positions,
            'velocities': velocities,
            'tracking': tracking,
        }
        # Add to pandas results container
        results_pd.add_simulation(method, times, positions, velocities, tracking)
    
    colors = {'euler': 'red', 'verlet': 'green', 'rk4': 'blue'}
    output_dir = Path(__file__).parent / 'data'
    output_dir.mkdir(exist_ok=True)
    
    # Generate separate plot for each method
    for method in ["euler", "verlet", "rk4"]:
        fig, axes = plt.subplots(1, 3, figsize=(16, 5))
        color = colors[method]
        
        pos = results[method]['positions']
        times = results[method]['times']
        tracking = results[method]['tracking']
        
        # Calculate distances from origin and average deviance
        distances = np.linalg.norm(pos, axis=1)
        avg_deviance = np.mean(np.abs(distances - 1.0))
        max_deviance = np.max(np.abs(distances - 1.0))
        
        # Plot 1: Trajectory with reference circle
        ax = axes[0]
        
        # Draw perfect reference circle (dashed)
        theta = np.linspace(0, 2*np.pi, 1000)
        ref_circle_x = np.cos(theta)
        ref_circle_y = np.sin(theta)
        ax.plot(ref_circle_x, ref_circle_y, 'k--', linewidth=2, label='Reference Circle (r=1)', alpha=0.5)
        
        # Plot trajectory
        ax.plot(pos[:, 0], pos[:, 1], color=color, linewidth=2, label=f'{method.upper()}: avg dev {avg_deviance:.2e}')
        ax.plot(pos[0, 0], pos[0, 1], 'o', color=color, markersize=8, label='Start')
        
        ax.set_xlabel('x', fontsize=11)
        ax.set_ylabel('y', fontsize=11)
        ax.set_title(f'{method.upper()} - Orbital Trajectory', fontsize=12, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.axis('equal')
        ax.plot(0, 0, 'k*', markersize=15, label='Central Body')
        
        # Plot 2: Energy Error
        ax = axes[1]
        energy_error = tracking['energy_error']
        ax.semilogy(times, energy_error, color=color, linewidth=2)
        ax.fill_between(times, energy_error, alpha=0.2, color=color)
        
        ax.set_xlabel('Time', fontsize=11)
        ax.set_ylabel('Relative Energy Error', fontsize=11)
        ax.set_title(f'{method.upper()} - Energy Conservation', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, which='both')
        
        # Add statistics text
        final_energy_error = energy_error[-1]
        max_energy_error = energy_error.max()
        stats_text = f'Final: {final_energy_error:.2e}\nMax: {max_energy_error:.2e}'
        ax.text(0.98, 0.05, stats_text, transform=ax.transAxes, 
                fontsize=9, verticalalignment='bottom', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Plot 3: Angular Momentum Error
        ax = axes[2]
        l_error = tracking['angular_momentum_error']
        ax.semilogy(times, l_error, color=color, linewidth=2)
        ax.fill_between(times, l_error, alpha=0.2, color=color)
        
        ax.set_xlabel('Time', fontsize=11)
        ax.set_ylabel('Relative Angular Momentum Error', fontsize=11)
        ax.set_title(f'{method.upper()} - Angular Momentum Conservation', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, which='both')
        
        # Add statistics text
        final_l_error = l_error[-1]
        max_l_error = l_error.max()
        stats_text = f'Final: {final_l_error:.2e}\nMax: {max_l_error:.2e}'
        ax.text(0.98, 0.05, stats_text, transform=ax.transAxes, 
                fontsize=9, verticalalignment='bottom', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        
        # Save individual figure
        output_path = output_dir / f'comparison_{method}.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_path}")
        plt.close()
    
    # Also create combined comparison figure
    print("\nGenerating combined comparison figure...")
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    # Plot 1: All trajectories on one plot
    ax = axes[0]
    theta = np.linspace(0, 2*np.pi, 1000)
    ref_circle_x = np.cos(theta)
    ref_circle_y = np.sin(theta)
    ax.plot(ref_circle_x, ref_circle_y, 'k--', linewidth=2, label='Reference Circle (r=1)', alpha=0.5)
    
    for method in ["euler", "verlet", "rk4"]:
        pos = results[method]['positions']
        distances = np.linalg.norm(pos, axis=1)
        avg_deviance = np.mean(np.abs(distances - 1.0))
        ax.plot(pos[:, 0], pos[:, 1], label=f'{method.upper()} (avg dev: {avg_deviance:.2e})', 
                color=colors[method], linewidth=1.5, alpha=0.7)
        ax.plot(pos[0, 0], pos[0, 1], 'o', color=colors[method], markersize=6)
    
    ax.set_xlabel('x', fontsize=11)
    ax.set_ylabel('y', fontsize=11)
    ax.set_title('Orbital Trajectories Comparison', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.axis('equal')
    ax.plot(0, 0, 'k*', markersize=12)
    
    # Plot 2: Energy Error comparison
    ax = axes[1]
    for method in ["euler", "verlet", "rk4"]:
        energy_error = results[method]['tracking']['energy_error']
        ax.semilogy(energy_error, label=method.upper(), color=colors[method], linewidth=1.5)
    
    ax.set_xlabel('Step', fontsize=11)
    ax.set_ylabel('Relative Energy Error', fontsize=11)
    ax.set_title('Energy Conservation Comparison', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    
    # Plot 3: Angular Momentum Error comparison
    ax = axes[2]
    for method in ["euler", "verlet", "rk4"]:
        l_error = results[method]['tracking']['angular_momentum_error']
        ax.semilogy(l_error, label=method.upper(), color=colors[method], linewidth=1.5)
    
    ax.set_xlabel('Step', fontsize=11)
    ax.set_ylabel('Relative Angular Momentum Error', fontsize=11)
    ax.set_title('Angular Momentum Conservation Comparison', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    
    # Save combined figure
    output_path = output_dir / 'comparison.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()
    
    # Export all results to CSV using pandas
    print("\nExporting results to CSV...")
    print("-" * 50)
    
    # Export summary statistics
    summary_csv = output_dir / 'comparison_summary.csv'
    results_pd.export_summary_csv(summary_csv)
    
    # Export full time-series data
    results_pd.export_all_csv(output_dir)
    
    print("-" * 50)
    print("Comparison complete!")


if __name__ == "__main__":
    plot_comparison()
