import sys
import os
import warnings
import pandas as pd
import numpy as np
import statistics
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import sem
from scipy.signal import butter, filtfilt, savgol_filter
from scipy.optimize import curve_fit, differential_evolution
import matplotlib
from dataclasses import dataclass
import Constants
from Constants import ISO_Percents, END_ROM, pCas, Stiffness_Time, ktr_start_time, ktr_end_time, SampleRate
from Colors import Firebrick, DeepBlue, SkyBlue, Sienna, SeaGreen, Charcoal, LightGray, Black


@dataclass
class FileInfo:

    # Possible information that can be grabbed from filename
    Subject: dict = None
    Muscle: dict = None
    Fibre: dict = None
    pCa: dict = None
    Test: str = None
    rFE_Method: str = None
    Filename: str = None
    Loads: list = None
    LongitudinalTimepoint: dict = None

    OrganizedData: dict = None

    # # Dataframe created from raw text file
    # Data: pd.DataFrame = None

    # # Variables for all single fibre tests
    # FibreLength: dict = None
    # SarcomereLength: dict = None
    # CSA: dict = None
    
    # # Variables from for single fibre pCa test
    # PeakForces: dict = None
    # SpecificForces: dict = None
    # PassiveForces: dict = None

    # # Variables for single fibre ktr test
    # ActiveStiffness: dict = None
    # ktr: dict = None
    # ktrForce: dict = None
    # ktrSpecificForce: dict = None
    # GoodnessFit: dict = None

    # # Variables for single fibre power test
    # Forces: dict = None
    # Velocities: dict = None

    # # Arrays used for ktr modeling
    # Xmodel: np.array = None
    # Ymodel: np.array = None

    # Variables for invivo tests
    Condition: str = None
    Timepoint: str = None
    Frequency: int = None
    PeakTorques: dict = None
    FatigueContractions: int = None
    IsotonicContractionData: pd.DataFrame = None
    ISO_Percent: str = None
    ISO_Torques: dict = None
    ISO_Velocities: dict = None
    ISO_Powers: dict = None
