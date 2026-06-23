"""
CBSE Class 12 Informatics Practices / Computer Science Project
Title: Orbital Dynamics Simulator & Conservation Analysis
Author: M. Fayazul Haque
School: Shining Star International School
Academic Session: 2026-27
Description: A simple program to simulate and analyze satellite orbits using
             different mathematical integration methods (Euler, Verlet, and RK4).
             Analyzes data using Pandas and plots graphs using Matplotlib.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =====================================================================
# 1. CORE PHYSICS FUNCTIONS
# =====================================================================

def get_acceleration(x, y, mu=1.0):
    """
    Calculates gravitational acceleration based on Newton's Law of Gravitation.
    Formula: a = -mu * r_vector / r^3
    """
    # Calculate distance from origin (r = sqrt(x^2 + y^2))
    r = np.sqrt(x**2 + y**2)
    
    # Avoid division by zero if satellite crashes into the center
    if r < 1e-10:
        return 0.0, 0.0
        
    ax = -mu * x / (r**3)
    ay = -mu * y / (r**3)
    return ax, ay

def calculate_energy(x, y, vx, vy, mu=1.0):
    """
    Calculates total mechanical energy (Kinetic Energy + Potential Energy).
    Total Energy = 0.5 * v^2 - mu / r
    """
    r = np.sqrt(x**2 + y**2)
    v_square = vx**2 + vy**2
    
    ke = 0.5 * v_square
    pe = -mu / r if r > 1e-10 else 0.0
    return ke + pe

def calculate_angular_momentum(x, y, vx, vy):
    """
    Calculates 2D angular momentum (L = r x v = x*vy - y*vx).
    """
    return x * vy - y * vx

# =====================================================================
# 2. NUMERICAL INTEGRATION METHODS (SIMULATORS)
# =====================================================================

def euler_step(x, y, vx, vy, dt, mu=1.0):
    """
    Euler Method (1st Order): Simplest method.
    New Position = Old Position + Velocity * dt
    New Velocity = Old Velocity + Acceleration * dt
    """
    ax, ay = get_acceleration(x, y, mu)
    
    new_x = x + vx * dt
    new_y = y + vy * dt
    
    new_vx = vx + ax * dt
    new_vy = vy + ay * dt
    
    return new_x, new_y, new_vx, new_vy

def verlet_step(x, y, vx, vy, dt, mu=1.0):
    """
    Velocity Verlet Method (2nd Order): Much better at preserving energy.
    It calculates position first, then finds new acceleration, and updates velocity.
    """
    ax, ay = get_acceleration(x, y, mu)
    
    # Update position: r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2
    new_x = x + vx * dt + 0.5 * ax * dt**2
    new_y = y + vy * dt + 0.5 * ay * dt**2
    
    # Get acceleration at the new position
    new_ax, new_ay = get_acceleration(new_x, new_y, mu)
    
    # Update velocity: v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
    new_vx = vx + 0.5 * (ax + new_ax) * dt
    new_vy = vy + 0.5 * (ay + new_ay) * dt
    
    return new_x, new_y, new_vx, new_vy

def rk4_step(x, y, vx, vy, dt, mu=1.0):
    """
    Runge-Kutta 4th Order (RK4) Method: Highly accurate.
    Calculates 4 intermediate steps (k1, k2, k3, k4) to estimate the next state.
    Written simply without complex vector arrays so it's easy to read.
    """
    # Helper to return derivatives of position and velocity: [vx, vy, ax, ay]
    def get_derivatives(px, py, vx_val, vy_val):
        ax_val, ay_val = get_acceleration(px, py, mu)
        return vx_val, vy_val, ax_val, ay_val

    # k1 (at start of step)
    dx1, dy1, dvx1, dvy1 = get_derivatives(x, y, vx, vy)
    
    # k2 (at midpoint using k1)
    x2 = x + 0.5 * dt * dx1
    y2 = y + 0.5 * dt * dy1
    vx2 = vx + 0.5 * dt * dvx1
    vy2 = vy + 0.5 * dt * dvy1
    dx2, dy2, dvx2, dvy2 = get_derivatives(x2, y2, vx2, vy2)
    
    # k3 (at midpoint using k2)
    x3 = x + 0.5 * dt * dx2
    y3 = y + 0.5 * dt * dy2
    vx3 = vx + 0.5 * dt * dvx2
    vy3 = vy + 0.5 * dt * dvy2
    dx3, dy3, dvx3, dvy3 = get_derivatives(x3, y3, vx3, vy3)
    
    # k4 (at end of step using k3)
    x4 = x + dt * dx3
    y4 = y + dt * dy3
    vx4 = vx + dt * dvx3
    vy4 = vy + dt * dvy3
    dx4, dy4, dvx4, dvy4 = get_derivatives(x4, y4, vx4, vy4)
    
    # Combine steps: weighted average of slopes
    new_x = x + (dt / 6.0) * (dx1 + 2*dx2 + 2*dx3 + dx4)
    new_y = y + (dt / 6.0) * (dy1 + 2*dy2 + 2*dy3 + dy4)
    new_vx = vx + (dt / 6.0) * (dvx1 + 2*dvx2 + 2*dvx3 + dvx4)
    new_vy = vy + (dt / 6.0) * (dvy1 + 2*dvy2 + 2*dvy3 + dvy4)
    
    return new_x, new_y, new_vx, new_vy

# Registry mapping integrator names to solver functions.
# To add a custom integrator, simply define its step function and register it here!
INTEGRATORS = {
    "euler": euler_step,
    "verlet": verlet_step,
    "rk4": rk4_step
}

# =====================================================================
# 3. SIMULATION RUNNER
# =====================================================================

def run_simulation(method, x0=1.0, y0=0.0, vx0=0.0, vy0=1.0, dt=0.005, steps=20000):
    """
    Runs the orbital simulation loop and stores everything in lists.
    Returns a Pandas DataFrame containing the full orbit data.
    """
    # Lists to store the time-series data
    times = []
    x_positions = []
    y_positions = []
    vx_velocities = []
    vy_velocities = []
    energies = []
    energy_errors = []
    angular_momenta = []
    momentum_errors = []
    
    # Set initial state
    x, y = x0, y0
    vx, vy = vx0, vy0
    t = 0.0
    
    # Get initial conserved quantities
    initial_energy = calculate_energy(x, y, vx, vy)
    initial_L = calculate_angular_momentum(x, y, vx, vy)
    
    # Simulation loop
    for i in range(steps + 1):
        # Calculate physics values
        energy = calculate_energy(x, y, vx, vy)
        L = calculate_angular_momentum(x, y, vx, vy)
        
        # Calculate relative errors (fractional deviance from starting values)
        e_err = abs(energy - initial_energy) / abs(initial_energy) if abs(initial_energy) > 1e-10 else 0.0
        l_err = abs(L - initial_L) / abs(initial_L) if abs(initial_L) > 1e-10 else 0.0
        
        # Append data to lists
        times.append(t)
        x_positions.append(x)
        y_positions.append(y)
        vx_velocities.append(vx)
        vy_velocities.append(vy)
        energies.append(energy)
        energy_errors.append(e_err)
        angular_momenta.append(L)
        momentum_errors.append(l_err)
        
        # Step forward using the registered solver
        if method in INTEGRATORS:
            x, y, vx, vy = INTEGRATORS[method](x, y, vx, vy, dt)
        else:
            raise ValueError(f"Unknown integration method: '{method}'. Registered: {list(INTEGRATORS.keys())}")
            
        t += dt
        
    # Convert lists to Pandas DataFrame
    df = pd.DataFrame({
        'time': times,
        'x': x_positions,
        'y': y_positions,
        'vx': vx_velocities,
        'vy': vy_velocities,
        'energy': energies,
        'energy_error': energy_errors,
        'angular_momentum': angular_momenta,
        'momentum_error': momentum_errors
    })
    
    # Calculate additional columns using Pandas
    df['distance'] = np.sqrt(df['x']**2 + df['y']**2)
    df['speed'] = np.sqrt(df['vx']**2 + df['vy']**2)
    
    return df

# =====================================================================
# 4. MAIN DRIVER AND VISUALIZATION
# =====================================================================

def main():
    print("==================================================================")
    # Aim and project introduction printed on startup
    print("PROJECT: Orbital Dynamics Simulator & Conservation Law Analysis")
    print("Student Name: M. Fayazul Haque")
    print("School: Shining Star International School")
    print("Academic Session: 2026-27")
    print("Class XII Physics & Computer Science/IP Project")
    print("==================================================================")
    
    # Simulation Parameters (Standard circular orbit starting values)
    dt = 0.005
    steps = 15000
    
    print(f"Running simulation with:")
    print(f"  Time Step (dt):   {dt}")
    print(f"  Total Steps:      {steps}")
    print(f"  Initial State:    Pos=(1.0, 0.0), Vel=(0.0, 1.0)\n")
    
    # Dict to hold results
    results = {}
    
    # Run simulation for all three methods
    for method in ["euler", "verlet", "rk4"]:
        print(f"-> Generating data for {method.upper()} Method...")
        df = run_simulation(method=method, dt=dt, steps=steps)
        results[method] = df
        
        # Save to CSV using Pandas
        filename = f"{method}_data.csv"
        df.to_csv(filename, index=False)
        print(f"   Saved simulation table to: {filename}")
        
    # Create directory for plots if not existing
    os.makedirs("plots", exist_ok=True)
    
    print("\n--------------------------------------------------------------")
    print("DATA ANALYSIS REPORT (Using Pandas Summary Statistics)")
    print("--------------------------------------------------------------")
    
    summary_list = []
    
    for method, df in results.items():
        # Compute mean distance deviance from perfect circle (r=1)
        mean_distance = df['distance'].mean()
        max_energy_error = df['energy_error'].max()
        max_momentum_error = df['momentum_error'].max()
        
        # Classify the orbit stability based on distance deviance
        std_dev_distance = df['distance'].std()
        if std_dev_distance > 0.1:
            orbit_type = "Unstable (Spiral)"
        else:
            orbit_type = "Stable (Circular)"
            
        summary_list.append({
            'Method': method.upper(),
            'Mean Radius': round(mean_distance, 4),
            'Radius Std Dev': round(std_dev_distance, 6),
            'Max Energy Error': f"{max_energy_error:.2e}",
            'Max L Error': f"{max_momentum_error:.2e}",
            'Orbit Type': orbit_type
        })
        
    # Display Summary Table using Pandas DataFrame display
    summary_df = pd.DataFrame(summary_list)
    print(summary_df.to_string(index=False))
    print("--------------------------------------------------------------\n")
    
    # =====================================================================
    # 5. MATPLOTLIB PLOTTING
    # =====================================================================
    print("Generating Matplotlib plots...")
    
    # Color scheme for the graphs
    colors = {'euler': 'tab:red', 'verlet': 'tab:green', 'rk4': 'tab:blue'}
    
    # Graph 1: Trajectory Comparison Plot
    plt.figure(figsize=(8, 8))
    
    # Draw reference circle of radius = 1 (Sun at the center)
    theta = np.linspace(0, 2 * np.pi, 1000)
    plt.plot(np.cos(theta), np.sin(theta), 'k--', label='Reference Orbit (r=1.0)', alpha=0.5)
    plt.plot(0, 0, 'y*', markersize=15, label='Star / Central Mass')
    
    # Plot satellite trajectories
    for method, df in results.items():
        plt.plot(df['x'], df['y'], label=f"{method.upper()}", color=colors[method], linewidth=1.5)
        # Mark the starting position
        plt.plot(df['x'].iloc[0], df['y'].iloc[0], 'ko', markersize=5)
        
    plt.title("Orbital Trajectories: Comparison of Integration Methods", fontsize=12, fontweight='bold')
    plt.xlabel("X Coordinate", fontsize=10)
    plt.ylabel("Y Coordinate", fontsize=10)
    plt.axis('equal')
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(loc='upper right')
    plt.savefig("plots/trajectory_comparison.png", dpi=150)
    print("   Saved plot: plots/trajectory_comparison.png")
    
    # Graph 2: Energy Conservation Plot
    plt.figure(figsize=(10, 5))
    for method, df in results.items():
        plt.plot(df['time'], df['energy_error'], label=f"{method.upper()}", color=colors[method], linewidth=1.5)
        
    plt.yscale('log') # Log scale is best to display tiny errors vs massive Euler errors
    plt.title("Relative Energy Error Over Time (Conservation of Energy)", fontsize=12, fontweight='bold')
    plt.xlabel("Simulation Time (seconds)", fontsize=10)
    plt.ylabel("Relative Energy Error (Log Scale)", fontsize=10)
    plt.grid(True, which="both", linestyle=':', alpha=0.6)
    plt.legend()
    plt.savefig("plots/energy_conservation.png", dpi=150)
    print("   Saved plot: plots/energy_conservation.png")
    
    # Graph 3: Angular Momentum Conservation Plot
    plt.figure(figsize=(10, 5))
    for method, df in results.items():
        plt.plot(df['time'], df['momentum_error'], label=f"{method.upper()}", color=colors[method], linewidth=1.5)
        
    plt.yscale('log')
    plt.title("Angular Momentum Error Over Time (Conservation of L)", fontsize=12, fontweight='bold')
    plt.xlabel("Simulation Time (seconds)", fontsize=10)
    plt.ylabel("Relative Momentum Error (Log Scale)", fontsize=10)
    plt.grid(True, which="both", linestyle=':', alpha=0.6)
    plt.legend()
    plt.savefig("plots/momentum_conservation.png", dpi=150)
    print("   Saved plot: plots/momentum_conservation.png")
    
    print("\nSimulation and plotting successfully completed! Go to the 'plots/' folder to view output graphs.")
    print("==================================================================")

if __name__ == "__main__":
    main()
