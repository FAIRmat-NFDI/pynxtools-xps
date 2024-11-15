# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
"""
Backgrounds for peak fitting.
"""

from typing import Optional, final
import numpy as np


class LinearBackground:
    """Linear background model for XPS spectra that connects the first and last points in the y array."""

    def __init__(self):
        """
        Initialize the linear background model without predefined parameters.
        The slope and intercept will be calculated based on the input data.
        """
        self.slope = None
        self.intercept = None

    @final
    def calc_background(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """
        Calculate the linear background based on the first and last points of the y array.

        Parameters:
            x (np.ndarray): The energy values at which to evaluate the background.
            y (np.ndarray): The observed intensity values at the corresponding energy points in x.

        Returns:
            np.ndarray: The background intensity at each energy point.
        """
        if len(x) != len(y):
            raise ValueError("x and y arrays must have the same length.")

        # Calculate slope and intercept using the first and last points
        self.slope = (y[-1] - y[0]) / (x[-1] - x[0])
        self.intercept = y[0] - self.slope * x[0]

        # Return the background values using the linear equation
        return self.slope * x + self.intercept

    @final
    def formula(self) -> str:
        """
        Returns the formula used for the linear background model.

        Returns:
            str: The formula for the linear background.
        """
        return (
            "Linear Background Formula:\n"
            "B(x) = slope * x + intercept\n"
            "Where:\n"
            "  slope: (y[-1] - y[0]) / (x[-1] - x[0])  # Calculated from the first and last points of y\n"
            "  intercept: y[0] - slope * x[0]           # Calculated from the first point of y\n"
            "  x: Energy values\n"
            "  y: Intensity values\n"
        )


class Shirley:
    """
    Shirley background subtraction using the Sherwood method.

    The Shirley method is used to compute the background of an x-ray or
    spectroscopy data set by iteratively subtracting the baseline until
    convergence. This class implements this method with options for
    tolerance and maximum iterations.
    """

    def __init__(self) -> None:
        """Initialize the Shirley background subtraction object."""
        pass

    @final
    def calc_background(
        self, x: np.ndarray, y: np.ndarray, tol: float = 1e-5, maxit: int = 15
    ) -> np.ndarray:
        """
        Calculate the Shirley background using the Sherwood method.

        Parameters:
        - x (np.ndarray): The x-axis data (energy or time, etc.).
        - y (np.ndarray): The y-axis data (intensity or counts, etc.).
        - tol (float, optional): The tolerance for the convergence criterion. Defaults to 1e-5.
        - maxit (int, optional): The maximum number of iterations to prevent infinite loops. Defaults to 15.

        Raises:
        - ValueError: If input arrays have mismatched dimensions, incorrect types, or the fit does not converge.

        Returns:
        - np.ndarray: The calculated Shirley background.
        """

        # Validate input
        if not isinstance(x, np.ndarray) or not isinstance(y, np.ndarray):
            raise ValueError(
                f"Parameters x and y must be numpy arrays, not {type(x)} and {type(y)}"
            )

        if len(x) != len(y):
            raise ValueError("x and y arrays must have the same length.")

        if x.size == 0 or y.size == 0:
            raise ValueError("x or y array is empty.")

        if x.ndim > 1 or y.ndim > 1:
            raise ValueError(
                f"Data arrays must be one-dimensional. Found shapes x: {x.shape}, y: {y.shape}."
            )

        # Reverse the data if it is in ascending order
        is_reversed = False
        if x[0] < x[-1]:
            is_reversed = True
            x = np.flip(x)
            y = np.flip(y)

        # Initialize background arrays
        background = np.zeros_like(x)
        background_next = np.zeros_like(x)

        # Iterative loop to compute Shirley background
        for iters in range(maxit):
            k = (y[0] - y[-1]) / np.trapz(y - background, x=x)

            # Update the background using trapezoidal integration
            for energy in range(len(x)):
                background_next[energy] = k * np.trapz(
                    y[energy:] - background[energy:], x=x[energy:]
                )

            # Check for convergence
            diff = np.linalg.norm(background_next - background)
            background = np.copy(background_next)
            if diff < tol:
                break
        else:
            # Raise an error if the maximum iterations are reached
            raise ValueError(
                "Maximum number of iterations exceeded before convergence."
            )

        # Reverse back the result if the data was originally reversed
        if is_reversed:
            return np.flip(y[-1] + background)
        return y[-1] + background

    def formula(self) -> str:
        """
        Returns the iterative formula used in the Shirley background subtraction.

        Returns:
        - str: The formula used for the Shirley background subtraction.
        """
        return (
            "Shirley Background Subtraction Formula:\n"
            "1. Initial background B_i = 0 for all i.\n"
            "2. In each iteration, calculate:\n"
            "   k = (y_0 - y_n) / integral(x_0, x_n, (y_j - B_j) dx_j)\n"
            "3. Update the background using the formula:\n"
            "   B_i = k * integral(x_i, x_n, (y_j - B_j) dx_j)\n"
            "4. Repeat steps 2 and 3 until convergence (|B_next - B| < tol).\n"
            "5. If convergence fails after maxit iterations, raise an error.\n"
            "Where:\n"
            "  B_i: Background at point i\n"
            "  k: Scaling factor based on the integral\n"
            "  y_i: Original data at point i\n"
            "  x_i: X-axis values\n"
            "  B: Current background estimate\n"
            "  tol: Convergence tolerance\n"
            "  maxit: Maximum number of iterations"
        )


class TougaardU3:
    """U3 Tougaard background model for XPS spectra.

    Parameters:
        E0 (float): Energy onset for the background, typically near the core-level peak energy.
        A (float): Scaling factor for the background intensity.
        C (float): Shape parameter that determines the background curvature.
    """

    def __init__(self, E0: float, A: float, C: float):
        """
        Initialize the U3 Tougaard model with the required parameters.

        Parameters:
            E0 (float): Energy onset for the background, typically near the core-level peak energy.
            A (float): Scaling factor for the background intensity.
            C (float): Shape parameter that determines the background curvature.
        """
        self.E0 = E0
        self.A = A
        self.C = C

    @final
    def calc_background(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Calculate the Tougaard background at each energy point.

        Args:
            x (np.ndarray): Array of energy values.

        Returns:
            np.ndarray: Tougaard background values corresponding to each energy value in x.
        """
        return self.A / ((x - self.E0) ** 2 + self.C)

    @final
    def formula(self) -> str:
        """
        Returns the formula used for the Tougaard background model.

        Returns:
            str: The formula for the Tougaard background.
        """
        return (
            "Tougaard U3 Background Formula:\n"
            "B(x) = A / ((x - E0)^2 + C)\n"
            "Where:\n"
            "  B(x): Background at energy x\n"
            "  A: Scaling factor\n"
            "  E0: Energy onset (near core-level peak energy)\n"
            "  C: Shape parameter (determines background curvature)\n"
        )

    @final
    def __repr__(self) -> str:
        """
        Returns a string representation of the GaussianLorentzianSum object, including details
        for position, width, intensity, and fraction_gauss.

        Returns:
        - str: The string representation of the object.
        """
        return f"TougaardU3(E0={self.E0}, A={self.A}, " f"C={self.C})"


class TougaardU4:
    """U4 Tougaard background model for XPS spectra."""

    def __init__(self, B: float, C: float, D: float, Eg: float, temp: float = 300.0):
        """
        Initialize the U4 Tougaard model with the required parameters.

        Parameters:
            B (float): A scaling factor that adjusts the overall amplitude of the background.
            C (float): A parameter influencing the background's shape, specifically the width of the energy region where the background decays.
            D (float): Another shape parameter that controls the fall-off of the background, influencing the tail of the background.
            Eg (float): The energy onset for the background, typically near the core-level peak energy, controlling where the background starts.
            temp (float): Temperature in Kelvin (default is 300 K).
        """
        self.B = B
        self.C = C
        self.D = D
        self.Eg = Eg
        self.temp = temp

    @final
    def calc_background(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """
        Calculate the modified Tougaard background at each energy point using a modified expression.

        Parameters:
            x (np.ndarray or float): The energy values at which to evaluate the background.

        Returns:
            np.ndarray or float: The background intensity at each energy point.
        """
        kb = 0.000086  # Boltzmann constant in eV/K

        return (
            (self.B * x)
            / ((self.C - x**2) ** 2 + self.D * x**2)
            * 1
            / (np.exp((self.Eg - x) / (self.temp * kb)) + 1)
        )

    @final
    def formula(self) -> str:
        """
        Returns the formula used for the U4 Tougaard background model.

        Returns:
            str: The formula for the U4 Tougaard background.
        """
        return (
            "U4 Tougaard Background Formula:\n"
            "B(x) = (B * x) / ((C - x^2)^2 + D * x^2) * 1 / (exp((Eg - x) / (kB * T)) + 1)\n"
            "Where:\n"
            "  B: Scaling factor\n"
            "  C: Shape parameter (affects width of energy decay region)\n"
            "  D: Shape parameter (controls the tail of the background)\n"
            "  Eg: Energy onset (near core-level peak energy)\n"
            "  kB: Boltzmann constant (0.000086 eV/K)\n"
            "  T: Temperature (default 300 K, can be modified)\n"
            "  x: Energy values\n"
        )

    @final
    def __repr__(self) -> str:
        """
        Returns a string representation of the GaussianLorentzianSum object, including details
        for position, width, intensity, and fraction_gauss.

        Returns:
        - str: The string representation of the object.
        """
        return (
            f"TougaardU4(B={self.B}, C={self.C}, D={self.D}, "
            f"Eg={self.Eg}, temp={self.temp})"
        )
