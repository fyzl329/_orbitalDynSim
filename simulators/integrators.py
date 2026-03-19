"""
Integration methods for orbital dynamics.
All integrators follow the same interface: (system, pos, vel, t, dt) -> (new_pos, new_vel)
"""
import numpy as np
from typing import Tuple


def euler_step(
    system,
    position: np.ndarray,
    velocity: np.ndarray,
    t: float,
    dt: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Euler forward integration step (1st order).
    
    r(t+dt) = r(t) + v(t)·dt
    v(t+dt) = v(t) + a(t)·dt
    
    Simple but least accurate. Poor for long-term orbital integrations.
    """
    accel = system.acceleration(position)
    new_pos = position + velocity * dt
    new_vel = velocity + accel * dt
    return new_pos, new_vel


def verlet_step(
    system,
    position: np.ndarray,
    velocity: np.ndarray,
    t: float,
    dt: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Velocity Verlet integration step (2nd order).
    
    r(t+dt) = r(t) + v(t)·dt + ½a(t)·dt²
    v(t+dt) = v(t) + ½[a(t) + a(t+dt)]·dt
    
    Better energy conservation than Euler, commonly used in MD simulations.
    """
    accel = system.acceleration(position)
    half_dt = 0.5 * dt
    
    # Update position
    new_pos = position + velocity * dt + 0.5 * accel * dt ** 2
    
    # Compute new acceleration
    new_accel = system.acceleration(new_pos)
    
    # Update velocity
    new_vel = velocity + (accel + new_accel) * half_dt
    
    return new_pos, new_vel


def rk4_step(
    system,
    position: np.ndarray,
    velocity: np.ndarray,
    t: float,
    dt: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Runge-Kutta 4th order integration step (4th order).
    
    y(t+dt) = y(t) + (dt/6)·(k₁ + 2k₂ + 2k₃ + k₄)
    where k_i are stage coefficients with intermediate states.
    
    Classic spectral accuracy method. Excellent for orbital mechanics.
    """
    # State vector: y = [x, y, vx, vy]
    def derivatives(state: np.ndarray) -> np.ndarray:
        """Returns [vx, vy, ax, ay]"""
        pos = state[:2]
        vel = state[2:4]
        accel = system.acceleration(pos)
        return np.concatenate([vel, accel])
    
    y = np.concatenate([position, velocity])
    
    # RK4 coefficients
    k1 = derivatives(y)
    k2 = derivatives(y + 0.5 * dt * k1)
    k3 = derivatives(y + 0.5 * dt * k2)
    k4 = derivatives(y + dt * k3)
    
    y_new = y + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
    
    new_pos = y_new[:2]
    new_vel = y_new[2:4]
    
    return new_pos, new_vel
