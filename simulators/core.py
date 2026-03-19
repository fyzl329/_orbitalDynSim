"""
Core orbital dynamics physics and unified simulation driver.
"""
import numpy as np
from dataclasses import dataclass
from typing import Tuple, Callable


@dataclass
class SimulationState:
    """State of a 2D orbital system."""
    time: float
    position: np.ndarray  # (2,) - [x, y]
    velocity: np.ndarray  # (2,) - [vx, vy]
    energy: float = 0.0
    angular_momentum: float = 0.0


class OrbitalSystem:
    """Core physics for 2-body orbital system (GM = 1 by convention)."""
    
    def __init__(self, mu: float = 1.0):
        """
        Initialize orbital system.
        
        Args:
            mu: Standard gravitational parameter (GM). Default 1.0 for normalized units.
        """
        self.mu = mu
    
    def acceleration(self, position: np.ndarray) -> np.ndarray:
        """
        Compute gravitational acceleration at position.
        
        Acceleration: a = -μ r / r³
        
        Args:
            position: (2,) array [x, y]
            
        Returns:
            (2,) array of acceleration [ax, ay]
        """
        r = np.linalg.norm(position)
        if r < 1e-10:
            return np.array([0.0, 0.0])
        return -self.mu * position / (r ** 3)
    
    def kinetic_energy(self, velocity: np.ndarray) -> float:
        """Compute kinetic energy (2D, m=1): KE = ½||v||²"""
        return 0.5 * np.dot(velocity, velocity)
    
    def potential_energy(self, position: np.ndarray) -> float:
        """Compute potential energy (2D, m=1): PE = -μ/r"""
        r = np.linalg.norm(position)
        if r < 1e-10:
            return 0.0
        return -self.mu / r
    
    def total_energy(self, position: np.ndarray, velocity: np.ndarray) -> float:
        """Total mechanical energy: E = ½||v||² - μ/r"""
        return self.kinetic_energy(velocity) + self.potential_energy(position)
    
    def angular_momentum(self, position: np.ndarray, velocity: np.ndarray) -> float:
        """
        Compute angular momentum magnitude (2D system).
        L = r × v = x·vy - y·vx (scalar in 2D)
        """
        return position[0] * velocity[1] - position[1] * velocity[0]


class UnifiedSimulationDriver:
    """
    Unified driver that runs all integrators under identical conditions.
    """
    
    def __init__(self, system: OrbitalSystem, integrator: Callable):
        """
        Initialize driver.
        
        Args:
            system: OrbitalSystem instance
            integrator: Integration method function (see integrators.py)
        """
        self.system = system
        self.integrator = integrator
    
    def run(
        self,
        initial_position: np.ndarray,
        initial_velocity: np.ndarray,
        dt: float,
        steps: int,
        track_energy: bool = True,
        track_angular_momentum: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
        """
        Run simulation with given parameters.
        
        Args:
            initial_position: Starting [x, y]
            initial_velocity: Starting [vx, vy]
            dt: Time step
            steps: Number of steps to integrate
            track_energy: Track total energy over time
            track_angular_momentum: Track angular momentum over time
            
        Returns:
            (times, positions, velocities, tracking_data)
        """
        times = np.zeros(steps + 1)
        positions = np.zeros((steps + 1, 2))
        velocities = np.zeros((steps + 1, 2))
        
        # Storage for tracked quantities
        energies = np.zeros(steps + 1) if track_energy else None
        angular_momenta = np.zeros(steps + 1) if track_angular_momentum else None
        
        # Initial conditions
        t = 0.0
        pos = initial_position.copy()
        vel = initial_velocity.copy()
        
        times[0] = t
        positions[0] = pos
        velocities[0] = vel
        
        if track_energy:
            energies[0] = self.system.total_energy(pos, vel)
        if track_angular_momentum:
            angular_momenta[0] = self.system.angular_momentum(pos, vel)
        
        # Integration loop
        for i in range(steps):
            pos, vel = self.integrator(
                self.system,
                pos,
                vel,
                t,
                dt
            )
            t += dt
            
            times[i + 1] = t
            positions[i + 1] = pos
            velocities[i + 1] = vel
            
            if track_energy:
                energies[i + 1] = self.system.total_energy(pos, vel)
            if track_angular_momentum:
                angular_momenta[i + 1] = self.system.angular_momentum(pos, vel)
        
        # Prepare tracking data
        tracking_data = {}
        if track_energy and energies is not None:
            tracking_data['energy'] = energies
            tracking_data['energy_error'] = np.abs(energies - energies[0]) / np.abs(energies[0])
        
        if track_angular_momentum and angular_momenta is not None:
            tracking_data['angular_momentum'] = angular_momenta
            tracking_data['angular_momentum_error'] = np.abs(
                angular_momenta - angular_momenta[0]
            ) / np.abs(angular_momenta[0])
        
        return times, positions, velocities, tracking_data


def run_simulation(
    method: str,
    initial_position: np.ndarray = None,
    initial_velocity: np.ndarray = None,
    dt: float = 0.01,
    steps: int = 10000,
    mu: float = 1.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    """
    High-level unified interface to run simulations.
    
    Args:
        method: Integration method ('euler', 'verlet', 'rk4')
        initial_position: [x, y] - defaults to circular orbit
        initial_velocity: [vx, vy] - defaults to circular orbit
        dt: Time step
        steps: Number of integration steps
        mu: Gravitational parameter
        
    Returns:
        (times, positions, velocities, tracking_data)
    """
    from . import integrators
    
    # Default to circular orbit at r=1
    if initial_position is None:
        initial_position = np.array([1.0, 0.0])
    if initial_velocity is None:
        initial_velocity = np.array([0.0, 1.0])
    
    # Get integrator function
    integrator_map = {
        'euler': integrators.euler_step,
        'verlet': integrators.verlet_step,
        'rk4': integrators.rk4_step,
    }
    
    if method not in integrator_map:
        raise ValueError(f"Unknown method: {method}. Choose from {list(integrator_map.keys())}")
    
    integrator = integrator_map[method]
    
    # Create system and driver
    system = OrbitalSystem(mu=mu)
    driver = UnifiedSimulationDriver(system, integrator)
    
    # Run simulation
    return driver.run(
        initial_position,
        initial_velocity,
        dt=dt,
        steps=steps,
        track_energy=True,
        track_angular_momentum=True,
    )
