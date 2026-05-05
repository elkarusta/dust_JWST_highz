import numpy as np
from scipy.optimize import root_scalar
from numpy.typing import NDArray


def stellar_mass_z_to_metallicity(
        redshift: float, 
        stellar_mass: float | NDArray[np.floating],
    ) -> float | NDArray[np.floating]:
    
    """
    Parameters
    ----------
    redshift: float
        Redshift of the galaxy 
    stellar_mass: float
        Stellar mass of the galaxy [log10(M/M_sun)]
        

    Returns
    -------
    metallicity: float
        Metallicity of the galaxy [log10(Z/Z_sun)]
        

    Notes
    -----
    Using the scaling relation from the FIRE-2 simulations or the Santa Cruz SAM

    References
    ----------
    Marszewsi et al. (2024), ...

    """

    a = 0.37
    b = -4.3

    metallicity = a * stellar_mass + b 

    #redshift-dependent case
    #a = 0.37
    #b = -4.25
    #alpha = 0.21
    #beta = 0.15
    #metallicity = a * ( (1+z)/9 )**alpha *stellar_mass + b * ( (1+z)/9 )**beta

    return metallicity



def mc_evolving_IMF(
        redshift: float, 
        metallicity: float, 
        imf_minimum_mass: float, 
        imf_maximum_mass: float, 
) -> float:
    
    """
    Parameters
    ----------
    redshift: float
        Redshift of the galaxy
    metallicity: float
        Metallicity of the galaxy [log10(Z/Z_sun)]
    imf_minimum_mass: float
        Minimum initial mass of stars (M_sun)
    imf_maximum_mass: float
        Maximum initial mass of stars (M_sun)

        
    Returns
    -------
    Mc: float
        Cutoff mass, defined as the initial stellar mass of the IMF slope change (M_sun)
        

    Notes
    -----
    Calculation of the cutoff mass for a two-component IMF: 
    - first component: salpeter slope (2.35)
    - second component: log-flat (1)
    The cutoff mass depends on redshift and stellar metallicity

    
    References
    ---------
    Evolving IMF recipe from Cueto et al. (2024), based on Chon et al. (2022) simulation results


    """
    
    # fraction of stellar mass in the log-flat IMF component, depends on redshift and metallicity
    x = 1 + Z
    f_massive =  1.07*(1-(2**x)) + 0.04*redshift*(2.67**x)
        
    # imposing 0 < f_massive < 1
    if f_massive <= 0: 
        f_massive = 1e-6 #Mc=100
    elif f_massive >= 1: 
       f_massive = 1 - 1e-6 #Mc=0.1
       
    Mi = imf_minimum_mass
    Mf = imf_maximum_mass

    # cutoff mass estimation from imposing the fraction of stellar mass in the log-flat IMF component
    def equation_for_Mc(Mc):
         # integral of first component
         I1 = (Mi**-0.35 - Mc**-0.35) / 0.35 
         # integral of second component, imposing continuity between the two 
         I2 = Mc**-1.35 * (Mf - Mc) 
         return I2 / (I1 + I2) - f_massive
    
    sol = root_scalar(equation_for_Mc, bracket=[Mi + 1e-6, Mf - 1e-6])
    Mc = sol.root
    
    return Mc



def select_SB99_tables(
        metallicity: float, 
        imf_cutoff_mass: float, 
        imf_maximum_mass: float, 
        SN_maximum_mass: float, 
) -> tuple[str, str, str, str]:

    """
    Parameters
    ----------
    metallicity: float
        Metallicity of the stellar population [log10(Z/Zsun)]
    imf_cutoff_mass: float
        Cutoff mass, defined as the initial mass of the IMF slope change (M_sun)
    imf_maximum_mass: float
        Maximum initial mass of stars (M_sun)
    SN_maximum_mass: float
        Maximum initial mass of stars that explode as supernovae (M_sun)

    Returns
    -------
    filename_spectra: string
        Name of the SB99 table with the spectrum
    filename_Ni: string
        Name of the SB99 table with the number of ionizing photons
    filename_L1500: string
        Name of the SB99 table with the UV luminosity at 1500A
    filename_snr: string
        Name of the SB99 table with the supernova rate

    Notes
    -----
    Selects filenames of pyStarburst99 (SB99) tables for input parameters

    """

    Zsun = 0.02 #Anders&Grevesse89, used in FIRE-2

    abs_Z = 10**(metallicity*Zsun)

    Z_values = np.array([0.0004, 0.001, 0.002, 0.1])
    Z_labels = ["0004", "001", "0002", "1"]

    idx = np.argmin(np.abs(Z_values - abs_Z))  # find closest
    Z = Z_labels[idx]

    Mc = np.round(imf_cutoff_mass, 1)
    Mf = int(imf_maximum_mass)
    SN = int(SN_maximum_mass)

    filename_spectra = f"directory/spectrum_inst_Z{Z}_Mc{Mc}_Mf{Mf}_SN{SN}.csv"
    filename_Ni = f"directory/Ni_inst_Z{Z}_Mc{Mc}_Mf{Mf}_SN{SN}.csv"
    filename_L1500 = f"directory/L1500_inst_Z{Z}_Mc{Mc}_Mf{Mf}_SN{SN}.csv"
    filename_snr = f"directory/snr_inst_Z{Z}_Mc{Mc}_Mf{Mf}_SN{SN}.csv"

    return filename_spectra, filename_Ni, filename_L1500, filename_snr

