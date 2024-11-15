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
Line shapes for peak fitting.
"""

from typing import Optional, final
import numpy as np


class Peak:
    def __init__(self, position: float, width: float, area: float) -> None:
        """
        Base class for a peak with position, width, and intensity.

        Parameters:
        - position (float): Position of the peak on the energy axis.
        - width (float): Full width at half maximum (FWHM) of the peak.
        - area (float): Area of the peak.
        """
        if width <= 0:
            raise ValueError("Width must be positive.")
        if area <= 0:
            raise ValueError("Intensity must be positive.")
        self.position = position
        self.width = width
        self.area = area


class Gaussian(Peak):
    """
    Gaussian peak with specified position, width, and area.
    """

    @final
    def calc_lineshape(self, x: float) -> Optional[float]:
        """
        Calculate the Gaussian lineshape at a given energy value x.

        Parameters:
        - x (float): Energy value on the energy axis.

        Returns:
        - float: Value of the Gaussian lineshape at x, or None if width is zero.
        """
        if self.width != 0:
            # Calculate intensity from area
            intensity = self.area / (np.sqrt(2 * np.pi) * self.width)
            exponent = -4 * np.log(2) * (x - self.position) ** 2 / self.width**2
            return intensity * np.exp(exponent)
        return None

    @final
    def formula(self) -> str:
        """
        Returns a string representation of the Gaussian formula.

        Returns:
        - str: The formula used for the Gaussian lineshape.
        """
        return "G(x) = (area / (sqrt(2 * pi) * width)) * exp[-(4 * ln(2) * (x - position)^2) / width^2]"

    @final
    def __repr__(self) -> str:
        return (
            f"Gaussian(position={self.position}, width={self.width}, area={self.area})"
        )


class Lorentzian(Peak):
    """Lorentzian peak with specified position, width, and area."""

    @final
    def calc_lineshape(self, x: float) -> Optional[float]:
        """
        Calculate the Lorentzian lineshape at a given energy value x.

        Parameters:
        - x (float): Energy value on the energy axis.

        Returns:
        - float: Value of the Lorentzian lineshape at x, or None if width is zero.
        """
        if self.width != 0:
            # Calculate intensity from area
            intensity = self.area / (np.pi * self.width)
            return intensity / (1 + (4 * (x - self.position) ** 2) / self.width**2)
        return None

    @final
    def formula(self) -> str:
        """
        Returns a string representation of the Lorentzian formula.

        Returns:
        - str: The formula used for the Lorentzian lineshape.
        """
        return "L(x) = (area / (pi * width)) / (1 + (4 * (x - position)^2) / width^2)"

    @final
    def __repr__(self) -> str:
        return f"Lorentzian(position={self.position}, width={self.width}, area={self.area})"


class LorentzianAsymmetric:
    def __init__(
        self, position: float, width: float, area: float, alpha: float, beta: float
    ):
        """
        Initialize the LorentzianAsymmetric profile with the required parameters.

        Parameters:
            position (float): The peak position (center).
            width (float): Full width at half maximum (FWHM) of the peak.
            area (float): Area under the peak.
            alpha (float): Scaling factor for width on the right.
            beta (float): Scaling factor for width on the left.
        """
        self.position = position
        self.width = width
        self.area = area
        self.alpha = alpha
        self.beta = beta

    @final
    def calc_lineshape(self, x: np.ndarray) -> np.ndarray:
        """
        Calculate the asymmetric Lorentzian lineshape for an array of energy values x.

        Parameters:
        - x (np.ndarray): Array of energy values on the energy axis.

        Returns:
        - np.ndarray: Array of values of the Lorentzian lineshape for each x.
        """
        # Calculate intensity from area
        intensity = self.area / (np.pi * self.width)

        # Calculate the width based on the position and scaling factors alpha and beta
        width = np.where(
            x < self.position, self.width / self.alpha, self.width / self.beta
        )

        # Calculate the Lorentzian lineshape for each value in x
        return intensity / (1 + (4 * (x - self.position) ** 2) / width**2)

    @final
    def formula(self) -> str:
        """
        Returns a string representation of the asymmetric Lorentzian formula.

        Returns:
        - str: The formula used for the asymmetric Lorentzian lineshape.
        """
        return (
            "L(x) = (area / (pi * width)) / (1 + ((x - position) / (width / alpha))^2) for x <= position, "
            "(area / (pi * width)) / (1 + ((x - position) / (width / beta))^2) for x > position"
        )

    @final
    def __repr__(self) -> str:
        """
        Return a string representation of the LorentzianAsymmetric object.
        """
        return (
            f"LorentzianAsymmetric(position={self.position}, width={self.width}, "
            f"area={self.area}, alpha={self.alpha}, beta={self.beta})"
        )


class LorentzianFinite(LorentzianAsymmetric):
    """Finite Lorentzian peak with specified position, width, area, asymmetry parameters, and damping."""

    def __init__(
        self,
        position: float,
        width: float,
        area: float,
        alpha: float,
        beta: float,
        w: float,
        gauss_contribution: float,
        no_of_convolutions: int,
    ):
        super().__init__(position, width, area, alpha, beta)
        self.w = w
        self.gauss_contribution = gauss_contribution
        self.no_of_convolutions = no_of_convolutions

    @final
    def calc_lineshape(self, x: np.ndarray) -> Optional[np.ndarray]:
        """
        Calculate the finite Lorentzian lineshape with damping for an array of energy values x.

        Parameters:
        - x (np.ndarray): Array of energy values on the energy axis.

        Returns:
        - np.ndarray: Array of finite Lorentzian lineshape values for each x.
        """
        # Get the Lorentzian lineshape for the array of x
        lorentzian = super().calc_lineshape(x)

        # If lorentzian is not None, apply damping
        if lorentzian is not None:
            damping_factor = 1 / (1 + 4 * ((x - self.position) / self.w) ** 2)
            return lorentzian * damping_factor

        # Return None if lorentzian is None
        return None

    @final
    def formula(self) -> str:
        """
        Returns a string representation of the finite Lorentzian formula with damping.

        Returns:
        - str: The formula used for the finite Lorentzian lineshape.
        """
        return (
            "L(x) = (area / (pi * width)) / (1 + (4 * (x - position)^2) / width^2) * "
            "(1 / (1 + 4 * ((x - position) / w)^2))"
        )

    @final
    def __repr__(self) -> str:
        """
        Return a string representation of the LorentzianFinite object.
        """
        return (
            f"LorentzianFinite(position={self.position}, width={self.width}, "
            f"area={self.area}, alpha={self.alpha}, beta={self.beta}, w={self.w})"
        )


class GaussianLorentzianSum(Peak):
    """Combined Gaussian and Lorentzian peak using the existing Gaussian and Lorentzian classes."""

    @final
    def __init__(
        self, position: float, width: float, area: float, fraction_gauss: float = 0.5
    ) -> None:
        """
        Combined Gaussian and Lorentzian sum peak.

        Parameters:
        - position (float): Position of the peak.
        - width (float): Width of the peak.
        - area (float): Area of the peak (instead of intensity).
        - fraction_gauss (float): Fraction of the Gaussian contribution (between 0 and 1).
        """
        super().__init__(position, width, area)
        self.fraction_gauss = fraction_gauss

    @final
    def calc_lineshape(self, x: float) -> Optional[float]:
        """
        Calculate the combined lineshape of the Gaussian and Lorentzian at x.

        Parameters:
        - x (float): Energy value on the energy axis.

        Returns:
        - float: Combined lineshape value at x, or None if width is zero.
        """
        if self.width != 0:
            # Calculate intensity from area for Gaussian and Lorentzian
            intensity = self.area / (np.pi * self.width)

            # Gaussian part (1 - fraction_gauss) contribution
            gauss_part = (1 - self.fraction_gauss) * Gaussian(
                self.position, self.width, intensity
            ).calc_lineshape(x)
            # Lorentzian part (fraction_gauss) contribution
            lorentz_part = self.fraction_gauss * Lorentzian(
                self.position, self.width, intensity
            ).calc_lineshape(x)
            return gauss_part + lorentz_part
        return None

    @final
    def formula(self) -> str:
        """
        Returns a detailed string representation of the combined Gaussian-Lorentzian formula.

        Returns:
        - str: The formula used for the combined lineshape.
        """
        # Using the formula for both Gaussian and Lorentzian
        gauss_formula = "G(x) = (area / (pi * width)) * exp[-(4 * ln(2) * (x - position)^2) / width^2]"
        lorentz_formula = (
            "L(x) = area / (pi * width) / (1 + (4 * (x - position)^2) / width^2)"
        )
        combined_formula = (
            f"SGL(x): G(x) + L(x) = (1 - fraction_gauss) * ({gauss_formula}) + "
            f"fraction_gauss * ({lorentz_formula})"
        )
        return combined_formula

    @final
    def __repr__(self) -> str:
        """
        Returns a string representation of the GaussianLorentzianSum object, including details
        for position, width, area, and fraction_gauss.

        Returns:
        - str: The string representation of the object.
        """
        return (
            f"GaussianLorentzianSum(position={self.position}, width={self.width}, "
            f"area={self.area}, fraction_gauss={self.fraction_gauss})"
        )


class GaussianLorentzianProduct(Peak):
    def __init__(
        self, position: float, width: float, area: float, fraction_gauss: float = 0.5
    ) -> None:
        """
        Combined Gaussian and Lorentzian product peak.

        Parameters:
        - position (float): Position of the peak.
        - width (float): Width of the peak.
        - area (float): Area of the peak (instead of intensity).
        - fraction_gauss (float): Fraction of the Gaussian contribution.
        """
        super().__init__(position, width, area)
        self.fraction_gauss = fraction_gauss

    @final
    def calc_lineshape(self, x: float) -> Optional[float]:
        """
        Calculate the combined lineshape of the Gaussian and Lorentzian product at x.

        Parameters:
        - x (float): Energy value on the energy axis.

        Returns:
        - float: Combined lineshape value at x, or None if width is zero.
        """
        if self.width != 0:
            # Calculate intensity from area for Gaussian and Lorentzian
            intensity = self.area / (np.pi * self.width)

            # Gaussian part (1 - fraction_gauss) contribution
            gauss_part = (1 - self.fraction_gauss) * Gaussian(
                self.position, self.width, intensity
            ).calc_lineshape(x)
            # Lorentzian part (fraction_gauss) contribution
            lorentz_part = self.fraction_gauss * Lorentzian(
                self.position, self.width, intensity
            ).calc_lineshape(x)
            return gauss_part * lorentz_part
        return None

    @final
    def formula(self) -> str:
        """
        Returns a detailed string representation of the combined Gaussian-Lorentzian formula.

        Returns:
        - str: The formula used for the combined lineshape.
        """
        # Using the formula for both Gaussian and Lorentzian
        gauss_formula = "G(x) = (area / (pi * width)) * exp[-(4 * ln(2) * (x - position)^2) / width^2]"
        lorentz_formula = (
            "L(x) = area / (pi * width) / (1 + (4 * (x - position)^2) / width^2)"
        )
        combined_formula = (
            f"GL(x): G(x) * L(x) = (1 - fraction_gauss) * ({gauss_formula}) * "
            f"fraction_gauss * ({lorentz_formula})"
        )
        return combined_formula

    @final
    def __repr__(self) -> str:
        """
        Returns a string representation of the GaussianLorentzianProduct object, including details
        for position, width, area, and fraction_gauss.

        Returns:
        - str: The string representation of the object.
        """
        return (
            f"GaussianLorentzianProduct(position={self.position}, width={self.width}, "
            f"area={self.area}, fraction_gauss={self.fraction_gauss})"
        )


class DoniachSunjic(Peak):
    """Doniach-Sunjic profile for XPS peaks with asymmetry."""

    def __init__(self, position: float, width: float, area: float, beta: float) -> None:
        """
        Initialize the Doniach-Sunjic profile with the required parameters.

        Parameters:
            position (float): The peak position (center).
            width (float): Full width at half maximum (FWHM) of the peak.
            area (float): Area under the peak (instead of intensity).
            beta (float): Asymmetry parameter (1 for symmetric Lorentzian, <1 for left skew, >1 for right skew).
        """
        # Initialize the parent Peak class with area
        super().__init__(position, width, area)
        self.beta = beta

    @final
    def calc_lineshape(self, x: np.ndarray) -> np.ndarray:
        """
        Calculate the Doniach-Sunjic profile at each energy point.

        Parameters:
            x (np.ndarray or float): The energy values at which to evaluate the profile.

        Returns:
            np.ndarray or float: The intensity values at each energy point.
        """
        # Calculate the intensity from area for normalization
        intensity = self.area / (np.pi * self.width)

        # Calculate the Doniach-Sunjic profile
        return intensity / ((1 + ((x - self.position) / self.width) ** 2) ** self.beta)

    @final
    def formula(self) -> str:
        """
        Returns the formula used for the Doniach-Sunjic profile.

        Returns:
            str: The formula for the Doniach-Sunjic profile.
        """
        return (
            "Doniach-Sunjic Profile Formula:\n"
            "f(x) = area / (pi * Gamma) * ((1 + ((x - x0) / Gamma)^2)^beta)\n"
            "Where:\n"
            "  area: Area under the peak\n"
            "  x0: Peak position (center)\n"
            "  Gamma: FWHM (full width at half maximum)\n"
            "  beta: Asymmetry parameter (1 for symmetric Lorentzian)\n"
            "  x: Energy values\n"
        )

    @final
    def __repr__(self) -> str:
        """
        Returns a string representation of the DoniachSunjic object, including details
        for position, width, area, and beta.

        Returns:
        - str: The string representation of the object.
        """
        return (
            f"DoniachSunjic(position={self.position}, width={self.width}, "
            f"area={self.area}, beta={self.beta})"
        )
