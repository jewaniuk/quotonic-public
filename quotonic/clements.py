# Import Python Packages & Modules
from typing import Optional

import numpy as np


class Mesh:
    """Model of a linear Mach-Zehnder interferometer mesh arranged in the Clements configuration.

    Each mesh of Mach-Zehnder interferometers (MZIs) is classified by a number of optical modes $m$.
    Also, fabrication imperfections can optionally be modelled by providing the mean and standard deviation
    of the propagation losses $\\alpha_\mathrm{WG}$, the standard deviation on the splitting ratio of the
    nominally 50:50 directional couplers, and the lengths of components (MZI, phase shifters, flat sections
    in parallel with MZIs). The default lengths correspond to the components considered in J. Ewaniuk *et al.*,
    *Advanced Quantum Technologies* **X**, XX-XX (2023).

    This class features methods to manipulate a mesh once it is constructed, including the generation of its
    matrix representation (only unitary when $\\alpha_\mathrm{WG} = 0\\text{ dB}/\\text{cm}$) from MZI phase
    shifts and the decomposition of a matrix representation to identify the MZI phase shifts required to
    realize it. This decomposition follows the scheme of Clements *et al.*, *Optica* **3**, 1460-1465 (2016).
    """

    def __init__(
        self,
        numModes: int,
        alphaWG: Optional[float] = None,
        std_alphaWG: Optional[float] = None,
        std_SR: Optional[float] = None,
        ellMZI: float = 0.028668,
        ellPS: float = 0.0050,
        ellF: float = 0.028668,
    ):
        """Initialization of a MZI mesh arranged in the Clements configuration.

        The properties of the mesh are first saved, and its phase shifts are initialized. If losses are
        provided, then component-by-component loss probabilities are computed and stored either by
        selecting them randomly from a normal distribution (if a standard deviation is provided), or by
        computing the uniform losses for each component. The loss probabilities are concatenated into a
        single list that is arranged to match the arrangement of phase shifts during mesh encoding.
        If a standard deviation on the directional coupler splitting ratios is provided, then specific
        splitting ratios are selected randomly from a normal distribution for each directional coupler
        in the mesh.

        Args:
            numModes: Number of optical modes $m$
            alphaWG: Mean propagation losses in dB/cm $\\alpha_\mathrm{WG}$.
            std_alphaWG: Standard deviation of the propagation losses in dB/cm $\\alpha_\mathrm{WG}$.
            std_SR: Standard deviation of the directional coupler splitting ratio $t$.
            ellMZI: Characteristic length of a MZI $\ell_\mathrm{MZI}$
            ellPS: Characteristic length of a phase shifter $\ell_\mathrm{PS}$
            ellF: Characteristic length of a flat section in parallel with a MZI $\ell_\mathrm{F}$
        """

        # Store the provided properties of the mesh and initialize the phases
        self.numModes = numModes
        self.alphaWG = alphaWG  # dB/cm
        self.std_alphaWG = std_alphaWG  # dB/cm
        self.std_SR = std_SR
        self.ellMZI = ellMZI  # cm
        self.ellPS = ellPS  # cm
        self.ellF = ellF  # cm
        self.phases = np.zeros(numModes * numModes)

        # If losses are provided, compute the loss probability contributed by each component in the mesh
        if alphaWG is not None:
            # If losses are non-uniform, loss probabilities are selected randomly from a normal distribution
            if std_alphaWG is not None:
                # MZI losses
                mu_alphaMZI = 1.0 - np.power(10, -1.0 * ellMZI * alphaWG / 10)
                sigma_alphaMZI = std_alphaWG * ellMZI * np.log(10) * np.power(10, -1.0 * alphaWG * ellMZI / 10) / 10
                alphaMZI = np.random.normal(mu_alphaMZI, sigma_alphaMZI, numModes * (numModes - 1))

                # Phase shifter losses
                mu_alphaPS = 1.0 - np.power(10, -1.0 * ellPS * alphaWG / 10)
                sigma_alphaPS = std_alphaWG * ellPS * np.log(10) * np.power(10, -1.0 * alphaWG * ellPS / 10) / 10
                alphaPS = np.random.normal(mu_alphaPS, sigma_alphaPS, numModes)

                # Flat section losses
                mu_alphaF = 1.0 - np.power(10, -1.0 * ellF * alphaWG / 10)
                sigma_alphaF = std_alphaWG * ellF * np.log(10) * np.power(10, -1.0 * alphaWG * ellF / 10) / 10
                alphaF = np.random.normal(mu_alphaF, sigma_alphaF, numModes)

            # If losses are uniform, loss probabilities are computed just once for each component
            else:
                alphaMZI = 1.0 - np.power(10, -1.0 * ellMZI * alphaWG / 10) * np.ones(numModes * (numModes - 1))
                alphaPS = 1.0 - np.power(10, -1.0 * ellPS * alphaWG / 10) * np.ones(numModes)
                alphaF = 1.0 - np.power(10, -1.0 * ellF * alphaWG / 10) * np.ones(numModes)

            # Arrange the loss probabilities to match the arrangement of the encoding scheme
            if numModes > 2:
                self.alpha = np.zeros(len(alphaMZI) + len(alphaPS) + len(alphaF))
                ind = 0
                indMZI = 0
                indF = 0
                if numModes % 2 == 0:
                    for i in range(numModes):
                        if i % 2 == 0:
                            self.alpha[ind : ind + numModes] = alphaMZI[indMZI : indMZI + numModes]
                            indMZI += numModes
                        else:
                            self.alpha[ind : ind + 1] = alphaF[indF : indF + 1]
                            self.alpha[ind + 1 : ind + numModes - 1] = alphaMZI[indMZI : indMZI + numModes - 2]
                            self.alpha[ind + numModes - 1 : ind + numModes] = alphaF[indF + 1 : indF + 2]
                            indMZI += numModes - 2
                            indF += 2
                        ind += numModes
                else:
                    for i in range(numModes):
                        if i % 2 == 0:
                            self.alpha[ind : ind + numModes - 1] = alphaMZI[indMZI : indMZI + numModes - 1]
                            self.alpha[ind + numModes - 1 : ind + numModes] = alphaF[indF : indF + 1]
                        else:
                            self.alpha[ind : ind + 1] = alphaF[indF : indF + 1]
                            self.alpha[ind + 1 : ind + numModes] = alphaMZI[indMZI : indMZI + numModes - 1]
                        indMZI += numModes - 1
                        indF += 1
                        ind += numModes
            else:
                self.alpha = np.zeros(len(alphaMZI) + len(alphaPS))
                self.alpha[0:2] = alphaMZI
                ind = 2
            self.alpha[ind::] = alphaPS

        # If directional coupler splitting ratios vary, select them randomly from a normal distribution
        if std_SR is not None:
            self.SR = np.random.normal(0.5, std_SR, numModes * (numModes - 1))
        # If directional coupler splitting ratios do not vary, then they are all 50:50
        else:
            self.SR = 0.5 * np.ones(numModes * (numModes - 1))

    def set_phases(self, phases):
        self.phases = phases

    def set_alpha(self, alpha):
        self.alpha = alpha

    def set_SR(self, SR):
        self.SR = SR

    def dc_column(self, placementSpecifier, SR, pos=True):
        if pos:
            coeff = 1j
        else:
            coeff = -1j

        dc = np.eye(self.numModes, dtype=complex)
        for i, ii in enumerate(np.arange(placementSpecifier, self.numModes - 1, 2)):
            dc[ii : ii + 2, ii : ii + 2] = np.array(
                [[np.sqrt(SR[i]), coeff * np.sqrt(1.0 - SR[i])], [coeff * np.sqrt(1.0 - SR[i]), np.sqrt(SR[i])]],
                dtype=complex,
            )

        return dc

    def ps_column(self, placementSpecifier, phases):
        # Initialize the diagonal of the unitary matrix as an array of ones
        ps_diagonal = np.ones(self.numModes, dtype=complex)

        # Insert phase shifts along the diagonal acting on mode m only
        ps_diagonal[placementSpecifier : self.numModes - 1 : 2] = np.exp(1j * phases)

        # Convert the diagonal to a numModes x numModes matrix and return it
        return np.diag(ps_diagonal)

    def encode(self):
        if self.numModes > 2:
            # Initialize single-photon unitary as numModes x numModes identity matrix
            U = np.eye(self.numModes, dtype=complex)

            indP = 0
            indL = 0
            indSR = 0
            for i in range(self.numModes):
                # Compute placement specifier to act on specific adjacent modes
                placementSpecifier = i % 2

                # Compute number of phase shifts per transformation
                pspt = (self.numModes - placementSpecifier) // 2

                # Compute number of splitting ratios per transformation
                srpt = pspt

                # Compute all layers of the Tmn transformation
                dc_neg = self.dc_column(placementSpecifier, self.SR[indSR : indSR + srpt], pos=False)
                dc_pos = self.dc_column(placementSpecifier, self.SR[indSR + srpt : indSR + (2 * srpt)])
                ps_phi = self.ps_column(placementSpecifier, self.phases[indP : indP + pspt])
                ps_twotheta = self.ps_column(placementSpecifier, self.phases[indP + pspt : indP + (2 * pspt)])

                # Multiply layers with the previous single-photon unitary
                U = dc_neg.dot(ps_twotheta).dot(dc_pos).dot(ps_phi).dot(U)

                # LOSSES
                if self.alphaWG is not None:
                    L = np.sqrt(1.0 - self.alpha[indL : indL + self.numModes])
                    L = np.reshape(L, (self.numModes, 1))
                    U = np.multiply(L, U)
                    indL += self.numModes

                indP += 2 * pspt
                indSR += 2 * srpt

            # Finish encoding by multiplying the output phase shifts to each optical mode
            if self.alphaWG is not None:
                L = np.sqrt(1.0 - self.alpha[indL : indL + self.numModes])
                D = np.reshape(np.exp(1j * self.phases[indP : indP + self.numModes]) * L, (self.numModes, 1))
            else:
                D = np.reshape(np.exp(1j * self.phases[indP : indP + self.numModes]), (self.numModes, 1))
            U = np.multiply(D, U)

            return U
        else:
            SR1 = self.SR[0]
            SR2 = self.SR[1]
            phi = self.phases[0]
            twotheta = self.phases[1]
            U11 = np.sqrt(SR1 * SR2) * np.exp(1j * (phi + twotheta)) + np.sqrt((1.0 - SR1) * (1.0 - SR2)) * np.exp(
                1j * phi
            )
            U12 = 1j * np.sqrt(SR1 * (1.0 - SR2)) * np.exp(1j * twotheta) - 1j * np.sqrt(SR2 * (1.0 - SR1))
            U21 = -1j * np.sqrt(SR2 * (1.0 - SR1)) * np.exp(1j * (phi + twotheta)) + 1j * np.sqrt(
                SR1 * (1.0 - SR2)
            ) * np.exp(1j * phi)
            U22 = np.sqrt((1.0 - SR1) * (1.0 - SR2)) * np.exp(1j * twotheta) + np.sqrt(SR1 * SR2)
            U = np.array([[U11, U12], [U21, U22]])

            if self.alphaWG is not None:
                L = np.sqrt(1.0 - self.alpha[0:2])
                L = np.reshape(L, (2, 1))
                U = np.multiply(L, U)

                L = np.sqrt(1.0 - self.alpha[2:4])
                D = np.reshape(np.exp(1j * self.phases[2::]) * L, (2, 1))
            else:
                D = np.reshape(np.exp(1j * self.phases[2::]), (2, 1))
            U = np.multiply(D, U)

            return U
