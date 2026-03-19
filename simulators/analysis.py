"""
Data analysis module for orbital simulations using pandas.

Provides utilities for organizing, summarizing, and exporting simulation results.
"""
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple


class SimulationResults:
    """Organize and analyze simulation results using pandas DataFrames."""
    
    def __init__(self):
        """Initialize results container."""
        self.data: Dict[str, pd.DataFrame] = {}
        self.summaries: Dict[str, Dict] = {}
    
    def add_simulation(
        self,
        method: str,
        times: np.ndarray,
        positions: np.ndarray,
        velocities: np.ndarray,
        tracking: Dict,
    ) -> None:
        """
        Add simulation results to the collection.
        
        Args:
            method: Integrator name ("euler", "verlet", "rk4")
            times: Time array
            positions: N x 2 array of positions
            velocities: N x 2 array of velocities
            tracking: Dict with 'energy', 'energy_error', 'angular_momentum', 'angular_momentum_error'
        """
        # Calculate orbital metrics
        distances = np.linalg.norm(positions, axis=1)
        speed = np.linalg.norm(velocities, axis=1)
        
        # Create DataFrame
        df = pd.DataFrame({
            'time': times,
            'x': positions[:, 0],
            'y': positions[:, 1],
            'distance_from_origin': distances,
            'vx': velocities[:, 0],
            'vy': velocities[:, 1],
            'speed': speed,
            'energy': tracking['energy'],
            'energy_error': tracking['energy_error'],
            'angular_momentum': tracking['angular_momentum'],
            'angular_momentum_error': tracking['angular_momentum_error'],
        })
        
        self.data[method] = df
        
        # Calculate summary statistics
        self._compute_summary(method, df)
    
    def _compute_summary(self, method: str, df: pd.DataFrame) -> None:
        """Compute and store summary statistics for a method."""
        distances = df['distance_from_origin']
        energy_error = df['energy_error']
        ang_momentum_error = df['angular_momentum_error']
        
        self.summaries[method] = {
            'method': method,
            'avg_distance': distances.mean(),
            'avg_deviance': np.mean(np.abs(distances - 1.0)),
            'max_deviance': np.max(np.abs(distances - 1.0)),
            'min_distance': distances.min(),
            'max_distance': distances.max(),
            'avg_energy_error': energy_error.mean(),
            'max_energy_error': energy_error.max(),
            'min_energy_error': energy_error.min(),
            'final_energy_error': energy_error.iloc[-1],
            'avg_ang_momentum_error': ang_momentum_error.mean(),
            'max_ang_momentum_error': ang_momentum_error.max(),
            'final_ang_momentum_error': ang_momentum_error.iloc[-1],
        }
    
    def get_summary_dataframe(self) -> pd.DataFrame:
        """
        Get all summaries as a DataFrame for easy comparison.
        
        Returns:
            DataFrame with one row per method and summary statistics as columns
        """
        summaries_list = list(self.summaries.values())
        return pd.DataFrame(summaries_list).set_index('method')
    
    def get_method_data(self, method: str) -> pd.DataFrame:
        """Get full time-series data for a single method."""
        if method not in self.data:
            raise ValueError(f"No data for method: {method}")
        return self.data[method]
    
    def export_summary_csv(self, filepath: str) -> None:
        """Export summary statistics to CSV file."""
        summary_df = self.get_summary_dataframe()
        summary_df.to_csv(filepath)
        print(f"Summary exported to: {filepath}")
    
    def export_all_csv(self, directory: str) -> None:
        """Export full time-series data for all methods to separate CSV files."""
        from pathlib import Path
        
        output_dir = Path(directory)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        for method, df in self.data.items():
            filepath = output_dir / f"{method}_timeseries.csv"
            df.to_csv(filepath, index=False)
            print(f"Exported: {filepath}")
    
    def comparison_table(self) -> str:
        """
        Get formatted comparison table as string.
        
        Returns:
            Formatted string suitable for printing or logging
        """
        summary_df = self.get_summary_dataframe()
        # Select key columns for display
        display_cols = [
            'avg_deviance',
            'max_deviance',
            'max_energy_error',
            'final_energy_error',
            'max_ang_momentum_error',
        ]
        
        display_df = summary_df[display_cols]
        return display_df.to_string()
    
    def filter_by_threshold(self, method: str, column: str, threshold: float) -> pd.DataFrame:
        """
        Filter time-series data by a threshold on a column.
        
        Useful for finding timesteps where error exceeds a limit.
        
        Args:
            method: Integrator name
            column: Column name to filter on
            threshold: Threshold value
            
        Returns:
            Filtered DataFrame
        """
        df = self.get_method_data(method)
        return df[df[column] > threshold]
    
    def get_statistics(self, method: str) -> Dict:
        """Get dictionary of statistics for a method."""
        if method not in self.summaries:
            raise ValueError(f"No summary for method: {method}")
        return self.summaries[method].copy()


def compare_methods(results: SimulationResults, methods: List[str] = None) -> pd.DataFrame:
    """
    Compare multiple methods side-by-side.
    
    Args:
        results: SimulationResults object
        methods: List of method names. If None, compare all.
        
    Returns:
        Comparison DataFrame
    """
    if methods is None:
        methods = list(results.data.keys())
    
    comparison = results.get_summary_dataframe().loc[methods]
    return comparison
