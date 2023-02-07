import sys
import os
import re
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
from FileInfoClass import FileInfo
from dataclasses import dataclass
import Constants
import Colors
from Constants import ISO_Percents, END_ROM, Stiffness_Time, ktr_start_time, ktr_end_time, SampleRate
from Colors import Firebrick, DeepBlue, SkyBlue, Sienna, SeaGreen, Charcoal, LightGray, Black

class Error(Exception):
    pass

def CreateDataFrame(Model: str = None, Test: str = None, Timepoint: str = None, pCas: list = Constants.pCas) -> pd.DataFrame:
    Results = pd.DataFrame(index = [0]).fillna(0)

    def GenerateColumnNames(Columns: dict = None):
        temp = pd.DataFrame()
        for key in Columns.keys():
            ColumnNames = Columns[key]
            for SpecificColumnName in ColumnNames:
                temp[SpecificColumnName] = SpecificColumnName
        OrganizedColumns = pd.DataFrame(columns=temp.keys())
        return OrganizedColumns

    if Model == 'Single Fibre':
        if Test == 'pCa':
            for ColumnNames in ['Subject', 'Muscle', 'Fibre', 'Fibre Length', 'Sarcomere Length', 'CSA']:
                Results[ColumnNames] = None
            
            ActiveForce = pd.DataFrame()
            SpecificForce = pd.DataFrame()

            for pCa in pCas:
                ActiveForce[pCa] = [f'Active Force (pCa {pCa})']
                SpecificForce[pCa] = [f'Specific Force (pCa {pCa})']

            ActiveForce = GenerateColumnNames(ActiveForce)
            SpecificForce = GenerateColumnNames(SpecificForce)

            Results = pd.concat([Results, ActiveForce, SpecificForce], axis = 1)

        if Test == 'ktr':
            # for ColumnNames in ['Subject', 'Muscle', 'Fibre', 'Fibre Length', 'Sarcomere Length', 'CSA']:
            #     Results[ColumnNames] = None
            
            ActiveStiffness = pd.DataFrame()
            SpecificStiffness = pd.DataFrame()
            ktrColumns = pd.DataFrame()
            
            for pCa in pCas:    
                ActiveStiffness[pCa] = [f'Active Stiffness (pCa {pCa})']
                SpecificStiffness[pCa] = [f'Specific Stiffness (pCa {pCa})']
                ktrColumns[pCa] = ([f'ktr (/s) (pCa {pCa})'], [f'Goodness of fit (r^{2}) (pCa {pCa})'])
            

            ActiveStiffness = GenerateColumnNames(ActiveStiffness)
            SpecificStiffness = GenerateColumnNames(SpecificStiffness)
            ktrColumns = GenerateColumnNames(ktrColumns)

            Results = pd.concat([Results, ActiveStiffness, SpecificStiffness, ktrColumns], axis = 1)

        if Test == 'pCa and ktr':
            for ColumnNames in ['Subject', 'Timepoint', 'Muscle', 'Fibre', 'Fibre Length', 'Sarcomere Length', 'CSA']:
                Results[ColumnNames] = None
            
            ActiveForce = pd.DataFrame()
            SpecificForce = pd.DataFrame()
            ActiveStiffness = pd.DataFrame()
            SpecificStiffness = pd.DataFrame()
            ktrColumns = pd.DataFrame()
            GoodnessFitColumns = pd.DataFrame()
            ktrForceColumns = pd.DataFrame()
            ktrSpecificForceColumns = pd.DataFrame()

            for pCa in pCas:
                ActiveForce[pCa] = [f'Active Force (pCa {pCa})']
                SpecificForce[pCa] = [f'Specific Force (pCa {pCa})']
                ActiveStiffness[pCa] = [f'Active Stiffness (pCa {pCa})']
                SpecificStiffness[pCa] = [f'Specific Stiffness (pCa {pCa})']
                ktrColumns[pCa] = [f'ktr (/s) (pCa {pCa})']
                GoodnessFitColumns[pCa] = [f'Goodness of fit (r^{2}) (pCa {pCa})']
                ktrForceColumns[pCa] = [f'ktr Force (pCa {pCa})']
                ktrSpecificForceColumns[pCa] = [f'ktr Specific Force (pCa {pCa})']
            
            ActiveForce = GenerateColumnNames(ActiveForce)
            SpecificForce = GenerateColumnNames(SpecificForce)
            ActiveStiffness = GenerateColumnNames(ActiveStiffness)
            SpecificStiffness = GenerateColumnNames(SpecificStiffness)
            ktrColumns = GenerateColumnNames(ktrColumns)
            GoodnessFitColumns = GenerateColumnNames(GoodnessFitColumns)
            ktrForceColumns = GenerateColumnNames(ktrForceColumns)
            ktrSpecificForceColumns = GenerateColumnNames(ktrSpecificForceColumns)


            Results = pd.concat([Results, ActiveForce, SpecificForce, ActiveStiffness, SpecificStiffness, ktrColumns, GoodnessFitColumns, ktrForceColumns, ktrSpecificForceColumns], axis = 1)
            pass
    
        if Test == 'Power':
            for ColumnNames in ['Subject', 'Muscle', 'Fibre', 'Fibre Length', 'Sarcomere Length', 'CSA']:
                Results[ColumnNames] = None

            AbsForce_Columns = pd.DataFrame()
            NormForce_Columns = pd.DataFrame()
            AbsVelocity_Columns = pd.DataFrame()
            NormVelocity_Columns = pd.DataFrame()

            for i in Constants.PowerLoads:
                AbsForce_Columns[i] = [f'{i}% Force']
                NormForce_Columns[i] = [f'{i}% Normalized Force']
                AbsVelocity_Columns[i] = [f'{i}% Velocity']
                NormVelocity_Columns[i] = [f'{i}% Normalized Velocity']
            
            AbsForce_Columns = GenerateColumnNames(AbsForce_Columns)
            NormForce_Columns = GenerateColumnNames(NormForce_Columns)
            AbsVelocity_Columns = GenerateColumnNames(AbsVelocity_Columns)
            NormVelocity_Columns = GenerateColumnNames(NormVelocity_Columns)

            Results = pd.concat([Results, AbsForce_Columns, NormForce_Columns, AbsVelocity_Columns, NormVelocity_Columns], axis = 1)

        if Test == 'rFE':
            for ColumnNames in ['Subject', 'Muscle', 'Fibre', 'Fibre Length', 'Sarcomere Length', 'CSA']:
                Results[ColumnNames] = None

            rFEColumns = pd.DataFrame()

            for Method in ['FM', 'STAND']:
                for Protocol in ['ISO', 'rFE']:
                    rFEColumns[Method + Protocol] = [f'{Method} {Protocol} Force', f'{Method} {Protocol} Stiffness', f'{Method} {Protocol} Passive Force']

            rFEColumns = GenerateColumnNames(rFEColumns)

            Results = pd.concat([Results, rFEColumns], axis = 1)

            pass

    if Model == 'in Vivo':
        Hz_Columns = pd.DataFrame()
        Isotonic_Columns = pd.DataFrame()
        Recovery_Columns = pd.DataFrame()
        
        for i in Constants.HzFrequencies:
            Hz_Columns[i] = ([f'{Timepoint} {i} Hz TQ'])
        for i in Constants.ISO_Percents:
            Isotonic_Columns[i] = ([f'{Timepoint} {i}% ISO Velocity'],
                                   [f'{Timepoint} {i}% ISO Torque'])
        for i in Constants.RecoveryTimepoints:
            if i in {1, 2, 5, 10}:
                i = (f'{i} min')
            Recovery_Columns[i] = ([f'{Timepoint} {i} PLFFD Control Twitch TQ'],
                                [f'{Timepoint} {i} PLFFD Control Twitch RTD'],
                                [f'{Timepoint} {i} PLFFD Control Twitch HRT'],
                                [f'{Timepoint} {i} PLFFD 10 Hz TQ'],
                                [f'{Timepoint} {i} PLFFD 10 Hz RTD'],
                                [f'{Timepoint} {i} PLFFD 10 Hz HRT'],
                                [f'{Timepoint} {i} PLFFD 100 Hz TQ'],
                                [f'{Timepoint} {i} PLFFD 100 Hz RTD'],
                                [f'{Timepoint} {i} PLFFD 100 Hz HRT'],
                                [f'{Timepoint} {i} PLFFD Potentiated Twitch TQ'],
                                [f'{Timepoint} {i} PLFFD Potentiated Twitch RTD'],
                                [f'{Timepoint} {i} PLFFD Potentiated Twitch HRT'],
                                [f'{Timepoint} {i} PLFFD 10:100 Hz TQ Ratio'],
                                [f'{Timepoint} {i} ISO Torque'],
                                [f'{Timepoint} {i} ISO Velocity'])
        
        # for Columns in [Hz_Columns, Isotonic_Columns, Recovery_Columns]:
        #     GenerateColumnNames(Columns)

        Hz_Columns = GenerateColumnNames(Hz_Columns)
        Isotonic_Columns = GenerateColumnNames(Isotonic_Columns)
        Recovery_Columns = GenerateColumnNames(Recovery_Columns)
        Fatigue_Columns = pd.DataFrame(columns = [f'{Timepoint} # of Contractions'])

        Results = pd.concat([Results, Hz_Columns, Isotonic_Columns, Fatigue_Columns, Recovery_Columns], axis = 1)

    return Results

def ReadFile(File: str = None, Model: str = None, Test: str = None):
    PowerLoads = []
    Filename: str = os.path.basename(File)
    if Model == 'Single Fibre':
        # Get characteristic info from filename
    
        if Filename.__contains__('.dat'):
            Filename: list = Filename.split("_")
            Filename = Filename[0:4]
        else:
            Filename = Filename.split("_")

        # Filenames are different based on test
        # Grab correct information based off type of test
        if Test in ['pCa', 'ktr']:
            Subject: str = Filename[0]
            Muscle: str = Filename[2]
            Fibre: int = Filename[1]
            pCa = Filename[3]
            pCa: str = pCa[3:]

        if Test == 'pCa and ktr':
            Filename: str = os.path.basename(File)
            FilenameLower = Filename.lower()

            if re.match(r'post_|pre_', FilenameLower):
                Filename = FilenameLower
                pass
            if re.match(r'pre\w', FilenameLower):
                Filename = re.sub(r'pre', r'pre_', FilenameLower)
            if re.match(r'post\w', FilenameLower):
                Filename = re.sub(r'post', r'post_', FilenameLower)

            # Filename = re.sub(r'-', r'_', Filename)
            SplitFilename: list = Filename.split("_")
                

            Timepoint: str = SplitFilename[0]
            Subject: str = SplitFilename[1]
            Muscle: str = SplitFilename[2]
            Fibre: str = SplitFilename[3]
            Protocol: str = SplitFilename[4][0:3]
            pCa: str = SplitFilename[4][3:6]  

        if Test == 'Power':
            Subject: str = Filename[0]
            Muscle: str = Filename[1]
            Fibre: int = Filename[2]
            pCa = []
            if Filename[3].__contains__('LC'):
                Filename[3] = Filename[3][3:14]

            # Split end of filename by comma and append each isotonic load as int to list for future reference
            All_Loads = Filename[3].split(',')
            
            for l in All_Loads: 
                if l.__contains__('.dat'):
                    l = l.split('.')
                    Load = int(l[0])
                else:
                    Load = int(l)
            
                PowerLoads.append(Load)
            
        if Test == 'rFE':
            Subject: str = Filename[0]
            Fibre: str = Filename[1]
            Method: str = Filename[2]
            ContractionType: str = Filename[3]
            pass
        
        # Read first 200 rows of text file to get fibre characteristics from beginning of text file
        # Load columns [0, 1, 2, 3] as they contain info we need and speeds up the process
        Data: pd.DataFrame = pd.read_table(
            File,
            encoding = 'latin1',
            header = None,
            delim_whitespace = True,
            low_memory = False,
            on_bad_lines='skip',
            memory_map = True,
            usecols=[0, 1, 2, 3]
        )
        

        Data.columns = ['Time','Length', 'Stim','Force']

        # Grab desired fibre characteristic info
        # If any of the following fail, an error message will be printed to screen to warn user
        # Example of failure:
        ## if df['Stim'][10] doesn't equal a number (i.e., 1.50, 0.90, etc) 
        ## then Fibre_Length will fail and a warning message will appear
        try:
            FibreLength = float(Data['Stim'][10])
            SarcomereLength = float(Data['Force'][11])
            Diameter = float(Data['Length'][12])
            CSA = np.pi*(Diameter/2)**2
        except:
            print(Error(f"Physical characteristics for {File} not calculated correctly. Check text file for weird lines of data at beginning of file."))
            pass
        
        # Find start of actual data in text file
        DataStart = Data.index[Data['Time'] == '0.00'][0]

        Data = pd.DataFrame(Data[DataStart:], dtype = float).reset_index(drop = True)

        Data['Time'] = Data['Time'].div(1000)
        
        Data['mm Lo Ratio'] = Data['Length'].div(FibreLength)


        if Test in ['pCa', 'ktr']:
            return Data, Subject, Muscle, Fibre, pCa, FibreLength, SarcomereLength, CSA

        if Test == 'pCa and ktr':
            return Data, Timepoint, Subject, Muscle, Fibre, Protocol, pCa, FibreLength, SarcomereLength, CSA

        if Test == 'Power':
            return Data, Subject, Fibre, Muscle, FibreLength, SarcomereLength, CSA, PowerLoads

        if Test == 'rFE':
            return Data, Subject, Fibre, Method, ContractionType, FibreLength, SarcomereLength, CSA

    if Model == 'In Vivo':

        def ButterworthFilter(Data: pd.DataFrame = None, cutoff: int = None, order: int = None) -> pd.Series:

            b, a = butter(order, cutoff, btype='low', analog=False, fs = 10000)
            FilteredTorque = filtfilt(b, a, Data)

            return FilteredTorque

        # Get characteristic info from filename
        Filename = Filename.split("-")

        Condition: str = Filename[0]
        Subject: str = f'{Condition}_{Filename[1]}'
        Timepoint: str = Filename[2]

        if Test == 'Torque-Frequency':
            Hz: int = int(Filename[4])
        
        if Test == 'Torque-Velocity-Power':
            ISO_Percent: int = Filename[4]

        if Test == 'Recovery':
            RecoveryTimepoint: str = Filename[4]
            RecoveryTimepoint = f'{RecoveryTimepoint} min' if RecoveryTimepoint in {'1', '2', '5', '10'} else RecoveryTimepoint

        Data: pd.DataFrame = pd.read_table(
            File,
            encoding = 'latin1',
            header = None,
            delim_whitespace = True,
            low_memory = False,
            on_bad_lines='skip',
            memory_map = True,
            skiprows = 8,
            usecols=[0, 1, 2, 3, 11]
        )
        
        Data.columns = ['Time', 'Length', 'Raw Torque', 'Other', 'Stim']

        # Aurora Scientific includes length and torque scale factors in their data files
        LengthScale: float = float(Data['Raw Torque'][0])
        TorqueScale: float = float(Data['Stim'][0])

        # Find start of actual data in text file
        DataStart = Data.index[Data['Time'] == '0'][0]

        Data = pd.DataFrame(Data[['Time', 'Length', 'Raw Torque', 'Stim']][DataStart:], dtype = float).reset_index(drop = True)

        Data['Time'] = Data['Time'].div(10000)

        # Baseline torque values 
        if File.__contains__('TF' or 'Isotonic'):
            BaselineTorque = Data['Raw Torque'].iloc[15000:16000].mean()
        if File.__contains__('PLFFD'):
            BaselineTorque = Data['Raw Torque'].iloc[21000:22000].mean()
        else:
            BaselineTorque = Data['Raw Torque'].iloc[0:100].mean()

        Data['Raw Torque'] -= BaselineTorque

        # Scale length and torque channels based off Aurora Scientific values
        Data['Length'] = Data['Length'] * LengthScale
        Data['Raw Torque'] = Data['Raw Torque'] * TorqueScale

        # Filter torque signal with highpass (> 100 Hz) filter to eliminate stim artifacts
        Data['Filtered Torque'] = ButterworthFilter(Data['Raw Torque'], 100, 2)

        if Test == 'Torque-Frequency':
            return Data, Filename, Subject, Condition, Timepoint, Hz
        
        if Test == 'Torque-Velocity-Power':
            return Data, Filename, Subject, Condition, Timepoint, ISO_Percent
        
        if Test == 'Fatigue':
            return Data, Filename, Subject, Condition, Timepoint

        if Test == 'Recovery':
            return Data, Filename, Subject, Condition, Timepoint, RecoveryTimepoint



def pCaAnalysis(Data: pd.DataFrame = None, CSA: float = None, Graph: bool = False) -> tuple[pd.DataFrame, float, float]:
    # Subtract baseline force (first 50 ms of test) from force signal
    BaselineForce: float = Data['Force'][0:500].mean()
    Data['Force'] = Data['Force'] - BaselineForce

    # Use a subset of the test (times corresponding to 15-30s following test start) to find peak. 
    # First 15s are ignored so that bath changes aren't captured
    SubsetData = Data[150000:300000]
    
    # Find highest rolling 500 ms window in force
    PeakForce = SubsetData['Force'].rolling(window = 5000, center = True).mean().max()
    SpecificForce = PeakForce / CSA

    return Data, PeakForce, SpecificForce

def ktrAnalysis(Data: pd.DataFrame = None, Graph: bool = False) -> tuple[pd.DataFrame, float, float, float, pd.Series, pd.Series, float]:
    
    def ktr_model(x, a, kt, c):
        return a * (1-np.exp(-kt*x)) + c

    def generate_Initial_Parameters(x_data: pd.Series = None, y_data: pd.Series = None):
        # min and max used for bounds
        def sumOfSquaredError(parameterTuple):
            # do not print warnings by genetic algorithm
            warnings.filterwarnings("ignore")
            val = ktr_model(x_data, *parameterTuple)
            return(np.sum((y_data - val) ** 2.0))

        maxX = max(x_data)
        minX = min(x_data)
        maxY = max(y_data)
        minY = min(y_data)
        
        # find max force during ktr window. To account for files where force drops substantially
        # we compare whether the last 500 ms force is lower than 90% of peak force
        # if it is, we use 90% of peak force as our new maximal force to fix curve to
        Max_Force_param: float = maxY - minY
        if y_data[:-500].mean() < Max_Force_param * 0.9:
            Max_Force_param = Max_Force_param * 0.9
            # print(Error(f"The last 500 ms of force was < 90% of max force in file:{Filename}. Curve fit to 90% max force"))

        # Force at ktr start
        Force_at_T0: float = y_data[0] 
        
        parameterBounds = []
        # search bounds for a (force when at plateau)
        parameterBounds.append([Max_Force_param, Max_Force_param])
        # search bounds for kt (range of values software uses to find ktr)
        parameterBounds.append([0, 30])
        # searh bounds for c (force at t=0)
        parameterBounds.append([Force_at_T0, Force_at_T0])

        # "seed" the numpy random number generator for repeatable results
        result = differential_evolution(sumOfSquaredError, parameterBounds, seed=3)
        
        return result.x

    # Stiffness measurements
    StiffnessWindow = range(int(Stiffness_Time) - 100, int(Stiffness_Time) + 200)
    ForceWindow = range(int(Stiffness_Time) - 5001, int(Stiffness_Time) - 1)
    dF = (Data['Force'][StiffnessWindow]).max() - (Data['Force'][ForceWindow]).mean()
    dLo = (Data['mm Lo Ratio'][StiffnessWindow]).max() - (Data['mm Lo Ratio'][ForceWindow]).mean()
    Stiffness = dF/dLo

    # Fit ktr only to subset of data
    # ktr_start_time = following restretch (i.e., 10160)
    # ktr_end_time = 60000
    ModelData = pd.DataFrame(Data[['Time', 'Force']][Constants.ktr_start_time:Constants.ktr_end_time]).reset_index(drop = True)

    # Find min force value after restretch occurs
    # Becomes real start to ModelData
    min_force_index = ModelData[['Force']].idxmin()
    ktr_start: int = min_force_index[0]+80

    # Cutoff end of raw data to avoid any potential bath moves/movements at end of test
    # Would negatively affect curve fitting
    ktr_end: int = ktr_start + 50000

    # Put time and force data into numpy arrays for curve fitting
    x_data = np.array(ModelData['Time'][ktr_start:ktr_end])
    y_data = np.array(ModelData['Force'][ktr_start:ktr_end])
    
    # Find initial parameters for curve fitting
    ktr_Parameters = generate_Initial_Parameters(x_data, y_data)

    # maxfev = number of iterations code will attempt to find optimal curve fit
    maxfev:int = 1000
    try:
        fittedParameters, pcov = curve_fit(ktr_model, x_data, y_data, ktr_Parameters, maxfev=maxfev)
    except:
        try:
            maxfev = 5000
            fittedParameters, pcov = curve_fit(ktr_model, x_data, y_data, ktr_Parameters, maxfev=maxfev)
        except:
            # print(Error(f"ktr parameters were not fit after {maxfev} iterations for file: {os.path.basename(File)}. Added to 'Files to Check'"))
            pass
            
    
    # Generate model predictions
    modelPredictions = ktr_model(x_data, *fittedParameters)

    # Calculate error
    Max_Force_error: float = np.sqrt(pcov[0, 0])
    ktr_error: float = np.sqrt(pcov[1, 1])
    Force_at_T0_error: float = np.sqrt(pcov[2, 2])
    ktr: float = fittedParameters[1]

    absError = modelPredictions - y_data

    SE: float = np.square(absError)  # squared errors
    MSE: float = np.mean(SE)  # mean squared errors
    RMSE: float = np.sqrt(MSE)  # Root Mean Squared Error, RMSE
    GoodnessFit: float = 1.0 - (np.var(absError) / np.var(y_data))

    x_model: np.array = np.linspace(min(x_data), max(x_data), 100)
    y_model: np.array = ktr_model(x_model, *fittedParameters)
    
    # 
    Xmodel = np.linspace(min(x_data), max(x_data), 100)
    Ymodel = ktr_model(x_model, *fittedParameters)

    ktrForce = Data['Force'][-500:].mean()

    return Data, Stiffness, ktr, GoodnessFit, Xmodel, Ymodel, ktrForce

def PowerAnalysis(Data: pd.DataFrame = None, PowerLoads: list = None, CSA: float = None, FibreLength: float = None, Graph: bool = None) -> tuple[float, dict]:
    ForceVelocityData = {}
    

    MaxForce = Data['Force'][(195000 - 5000):195000].max()
    for idx, Load in enumerate(PowerLoads):

        Rows = range(int(Constants.ForceClampTimes[idx][0] + 500), int(Constants.ForceClampTimes[idx][1]))
        Rows2 = range(int(Constants.ForceClampTimes[idx][0] + 510), int(Constants.ForceClampTimes[idx][1]))
        if Load in [10, 20]:
            Length1 = Data['Length'][Rows]
            Length2 = Data['Length'][Rows2]
            temp = pd.DataFrame()

            PlateauWindow = 50
            
            temp['Length1'] = pd.DataFrame(Length1.rolling(PlateauWindow, center = True).mean()).dropna().reset_index(drop = 'true')
            temp['Length2'] = pd.DataFrame(Length2.rolling(PlateauWindow, center = True).mean()).dropna().reset_index(drop = 'true')

            
            temp['Compare'] =  temp.apply(lambda x: 'True' if x['Length2'] >= x['Length1'] else "", axis=1 )
            Comparison = temp.loc[(temp['Compare'] == 'True')]

        try:
            ClampEnd = Comparison.index[Comparison['Compare'] == 'True'][0] + Rows[0] 
        except:
            ClampEnd = int(Constants.ForceClampTimes[idx][1])

        Rows = range(int(Constants.ForceClampTimes[idx][0] + 500), ClampEnd)

        Start = Rows[0]
        End = Rows[-1]

        Force: float = Data['Force'][Start:End].mean()
    
        try:
            Velocity:float = ((Data['Length'][End] - Data['Length'][Start]) / (Data['Time'][End] - Data['Time'][Start]) * -1)
        except:
            Velocity = 0


        NormForce: float = Force / CSA
        NormVelocity: float = Velocity / FibreLength
        
        ForceVelocityData[Load] = {
            f'{Load}% Force': Force,
            f'{Load}% Normalized Force': NormForce,
            f'{Load}% Velocity': Velocity,
            f'{Load}% Normalized Velocity': NormVelocity
            }
        
    if Graph == True:
        fig, power = plt.subplots(nrows = 1, ncols = 2, figsize=(width := 10, height := 3), layout='constrained', dpi=100)
        power[0].plot(Data['Time'][180000:230000], Data['Force'][180000:230000], color = 'black', label = 'Force')
        power[0].set_xlabel('Time (s)')
        power[1].plot(Data['Time'][180000:230000], Data['Length'][180000:230000], color = 'black', label = 'Length')
        power[1].set_xlabel('Time (s)')

        Force_Inset = fig.add_axes([.05, .7, .20, .30])
        Force_Inset.plot(Data['Time'][195000:215000], Data['Force'][195000:215000], color = 'black', label = 'Force')
        Force_Inset.spines.right.set_visible(True)
        Force_Inset.spines.top.set_visible(True)
        Force_Inset.set_title('Force')

        Length_Inset = fig.add_axes([.6, .5, .20, .30])
        Length_Inset.plot(Data['Time'][195000:215000], Data['Length'][195000:215000], color = 'black', label = 'Length')
        Length_Inset.spines.right.set_visible(True)
        Length_Inset.spines.top.set_visible(True)
        Length_Inset.set_title('Velocity')
        
        Colors = [DeepBlue.Hex_Value, Firebrick.Hex_Value, Sienna.Hex_Value, SeaGreen.Hex_Value]
        for idx, Clamp in enumerate(Constants.ForceClampTimes):

            Calculation_Start = int(Clamp[0] + 500)
            Calculation_End = int(Clamp[1])

            if Calculation_Start == 210500:
                if ClampEnd is not None:
                    Calculation_End = int(ClampEnd)
            Mid_Point = round(Calculation_Start + 1000)

            Force:float = (Data['Force'][Calculation_Start:Calculation_End].mean())
            Velocity:float = ((Data['Length'][Calculation_End] - Data['Length'][Calculation_Start]) / (Data['Time'][Calculation_End] - Data['Time'][Calculation_Start]))

            color = Colors[idx]

            Force_Inset.plot(Data['Time'][Calculation_Start:Calculation_End], Data['Force'][Calculation_Start:Calculation_End], color = color, label = 'Calculations')
            Length_Inset.plot(Data['Time'][Calculation_Start:Calculation_End], Data['Length'][Calculation_Start:Calculation_End], color = color, label = 'Calculations')

            Force_Inset.text(x = Data['Time'][Mid_Point], y = (Force + (Force * 0.15)), s = f'Abs Force = {Force:.4f}')
            Length_Inset.text(x = Data['Time'][Mid_Point], y = (Data['Length'][Calculation_Start] + (Data['Length'][Calculation_Start] * 0.05)), s = f'{Velocity * -1:.4f}')
        plt.show()

    return MaxForce, ForceVelocityData

def rFEAnalysis(Data: pd.DataFrame = None, StiffnessTime: int = 65, SampleRate: int = Constants.SampleRate) -> tuple[float, float, float]:
    StiffnessTime = StiffnessTime * SampleRate
    Baseline_Window = range(int(40000), int(45000))
    Stiffness_Window = range(int(StiffnessTime) - 100, int(StiffnessTime) + 200)
    Force_Window = range(int(StiffnessTime) - 5001, int(StiffnessTime) - 1)
    Passive_Window = range(int(945000), int(950000))

    # Calculadora
    BaselineForce = Data['Force'][Baseline_Window].mean()

    dF = (Data['Force'][Stiffness_Window]).max() - (Data['Force'][Force_Window]).mean()
    dLo = (Data['mm Lo Ratio'][Stiffness_Window]).max() - (Data['mm Lo Ratio'][Force_Window]).mean()
    Stiffness = dF/dLo

    PeakForce = Data['Force'][Force_Window].mean() - BaselineForce

    PassiveForce = Data['Force'][Passive_Window].mean() - BaselineForce

    return PeakForce, PassiveForce, Stiffness

def Find_Peak_Force(Force: pd.DataFrame):
    temp = pd.DataFrame()
    counter = 0
    Peak_Force: float = ()
    Plateau_Window = 5000
    
    temp['roll'] = pd.DataFrame(Force.rolling(Plateau_Window, center = True).mean())
    PeakForce = Force.rolling(Plateau_Window, center = True).mean().max()
    for counter, i in enumerate(temp['roll']):
        if i == PeakForce:
            return counter 
            break
        else:
            counter += 1
    return PeakForce


def TorqueFrequencyCurveFitting(Files: list[FileInfo] = None, AnalyzeTimepoint: str = None, Model: str = None, Test: str = None):
    if Model == 'In Vivo' and Test == 'Torque-Frequency':
        AllTorqueFrequencyData = {}
        VCDMeans = {}
        VCDSEM = {}
        CONMeans = {}
        CONSEM = {}
        SubjectData = {}
        Torques = {}
        Frequencies = {}
        CurveFitResults = {}

        def TF_CurveFit_Equation(x: float, min_TQ: float, max_TQ: float,EC50: float, n: float):
            return min_TQ + (max_TQ - min_TQ) / (1 + ((x / EC50)**n)) 

        # Grab information from each FileInfo object class
        for File in Files:
            for Variable in File.OrganizedData:
                if Variable == 'Peak Torque Data':
                    for SpecificFrequency in File.OrganizedData[Variable].keys():
                            AllTorqueFrequencyData[File.Subject + File.Timepoint + str(File.Frequency)] = {
                                'Condition': File.Condition,
                                'Subject': File.Subject,
                                'Timepoint': File.Timepoint,
                                'Frequency': File.Frequency,
                                'Torque': File.OrganizedData[Variable][SpecificFrequency]
                            }
        
        # Transform dictionary into dataframe
        AllTorqueFrequencyData = pd.DataFrame.from_dict(AllTorqueFrequencyData, orient = 'index')

        # for Subject in AllTorqueFrequencyData['Subject']:
        #     Indices = AllTorqueFrequencyData.index[(AllTorqueFrequencyData['Subject'] == Subject) & (AllTorqueFrequencyData['Timepoint'] == Timepoint)]

        #     Torques[Subject] = AllTorqueFrequencyData['Torque'][Indices]
        #     Frequencies[Subject] = AllTorqueFrequencyData['Frequency'][Indices]

        #     print(Torques)
        
        #     plt.plot(list(Frequencies.values()), list(Torques.values()))
        #     plt.show()
            # plt.plot(list(Frequencies.values()), list(Torques.values()))

            # SubjectData[Subject] = {
            #     'Torque': AllTorqueFrequencyData.loc[(AllTorqueFrequencyData['Subject'] == Subject) & (AllTorqueFrequencyData['Timepoint'] == Timepoint),'Torque'],
            #     'Frequency': AllTorqueFrequencyData.loc[(AllTorqueFrequencyData['Subject'] == Subject) & (AllTorqueFrequencyData['Timepoint'] == Timepoint),'Frequency']
            #     }

        # for Subject in SubjectData.keys():
        #     for Variable in SubjectData[Subject].keys():
        #         if Variable == 'Torque':
        #             Torques[Variable] = SubjectData[Subject][Variable]
        #         if Variable == 'Frequency':
        #             pass
                    # print(f'{SpecificFrequency} = {SubjectData[Subject][SpecificFrequency]}')
            
        for Condition in ['VCD', 'CON']:
            for Frequency in Constants.HzFrequencies:
                if Frequency >= 20:
                    if Condition == 'VCD':
                        # Find average of torques for each frequency, condition, and timepoint
                        Data = AllTorqueFrequencyData.loc[
                            (AllTorqueFrequencyData['Frequency'] == Frequency) & 
                            (AllTorqueFrequencyData['Condition'] == Condition) & 
                            (AllTorqueFrequencyData['Timepoint'] == AnalyzeTimepoint), 
                            'Torque'
                        ]

                        VCDMeans[Frequency] = Data.mean()
                        VCDSEM[Frequency] = np.std(Data) / np.sqrt(len(Data))
        
                    if Condition == 'CON':
                        Data = AllTorqueFrequencyData.loc[
                            (AllTorqueFrequencyData['Frequency'] == Frequency) & 
                            (AllTorqueFrequencyData['Condition'] == Condition) & 
                            (AllTorqueFrequencyData['Timepoint'] == AnalyzeTimepoint), 
                            'Torque'
                        ]

                        CONMeans[Frequency] = Data.mean()
                        CONSEM[Frequency] = np.std(Data) / np.sqrt(len(Data))
        
        fig, RawTorqueCurve = plt.subplots(figsize = (width:= 5, height:= 3), dpi = 300, layout = 'constrained')
        RawTorqueCurve.errorbar(list(VCDMeans.keys()), list(VCDMeans.values()), yerr = list(VCDSEM.values()), marker = 'o', markersize = 4, linestyle = ' ', color = Colors.Black.Hex_Value, label = 'VCD')
        RawTorqueCurve.errorbar(list(CONMeans.keys()), list(CONMeans.values()), yerr = list(CONSEM.values()), marker = 's', markersize = 4, linestyle = ' ', color = Colors.LightGray.Hex_Value, label = 'CON')
        RawTorqueCurve.set_xlabel(Constants.Torque)
        RawTorqueCurve.set_ylabel('Frequency (Hz)')
        RawTorqueCurve.legend(loc = 4)

        for Condition in [VCDMeans, CONMeans]:
            popt, pcov = curve_fit(TF_CurveFit_Equation, list(Condition.keys()), list(Condition.values()))
            maxT, minT, Frequency_at_50percent, SlopeCoeff = popt
            xline = np.arange(min(Condition.keys()), max(Condition.keys()), 1)
            yline = TF_CurveFit_Equation(xline, maxT, minT, Frequency_at_50percent, SlopeCoeff)

            if Condition == VCDMeans:
                RawTorqueCurve.plot(xline, yline, color = Colors.Black.Hex_Value)

                CurveFitResults['VCD'] = {
                    'T50': Frequency_at_50percent,
                    'Slope Coeff': SlopeCoeff
                }

            if Condition == CONMeans:
                RawTorqueCurve.plot(xline, yline, Colors.LightGray.Hex_Value)

                CurveFitResults['CON'] = {
                    'T50': Frequency_at_50percent,
                    'Slope Coeff': SlopeCoeff
                }

    print(CurveFitResults)
    return CurveFitResults

def Torque_Frequency(Data: pd.DataFrame = None, Timepoint: str = None, Subject: str = None, Hz: int = None) -> float:
    # Find index of first stimulation (i.e., contraction start)
    StimIndex = Data.index[Data['Stim'] == 1][0]
    
    # Find the greatest 500 ms rolling average of torque signal
    # The footplate is moved to starting position at beginning of test so first portion of test (i.e., 20000 samples; 2 sec) is ignored
    PeakTorque = Data['Filtered Torque'][StimIndex:].rolling(window = 500, center = True).mean().max()

    if Timepoint == 'D120' and Subject in ['VCD_04', 'VCD_05', 'VCD_15', 'VCD_16']:
        PeakTorque = PeakTorque / 10

    return PeakTorque

def Isotonic_Contractions(Data: pd.DataFrame = None, Subject: str = None, Timepoint: str = None, UserInput: bool = False) -> pd.DataFrame:
    BaselineLength = Data['Length'].iloc[0:100].mean()

    # Find the first sample when footplate starts to move
    # Define as start of contraction
    # When UserInput == False, code attempts to find start of contraction based off criteria (i.e., i <= Baseline_Length * 0.99)
    # When UserInput == True, user will choose start of isotonic contraction so ISO_Start no longer needs to be found
    if UserInput == False:
        try:
            ISO_Start = Data.index[Data['Length'] <= BaselineLength * 0.99][0]
        except:
            ISO_Start: int = 0
        
        Data = Data.iloc[ISO_Start:]

    # Find first sample when footplate crosses end ROM (i.e., -18.99 degrees)
    # Define as end of contraction
    END_ROM = -20 if Timepoint == 'D80' or 'D120' else -18.99

    try:
        if Data['Length'].min() > END_ROM:
            ISO_End = Data.index[Data['Length'] == Data['Length'].min()][0]

        # If end ROM isn't achieved (possible during higher isotonic loads), then end ROM is final length sample in window
        if Data['Length'].min() < END_ROM:

            ISO_End = Data.index[Data['Length'] <= END_ROM][0] - ISO_Start
           
    except:
        print(Error(f'End ROM not found for {Subject}, {Timepoint}'))

    # Return dataframe containing only data during contraction
    return Data

def Torque_Velocity(Data: pd.DataFrame = None, ISO_Percent: str = None, Subject: str = None, Timepoint: str = None) -> tuple[float, float, float]:

    # Contractions occured at different times depending on test
    ContractionData = Data[20000:26000].reset_index(drop = True)

    # Isolate data relevant to isotonic contractions 
    try:
        ContractionData = Isotonic_Contractions(Data = ContractionData, Subject = Subject, Timepoint = Timepoint)
    except: 
        # Contractions where the footplate doesn't move at all have their velocity set at 0
        # May happen when testing against high isotonic loads, following fatigue protocol, etc.
        ISO_Velocity = 0
    else: 
        # Velocity calculated as difference between last and first samples of length and time channels, respectively
        LengthStart =  ContractionData['Length'].iloc[1]
        LengthEnd = ContractionData['Length'].iloc[-1]

        TimeStart = ContractionData['Time'].iloc[1]
        TimeEnd = ContractionData['Time'].iloc[-1]

        ISO_Velocity = ((LengthEnd - LengthStart) / (TimeEnd- TimeStart) * -1)

    # Mean of torque during contraction is used as torque value for power calculation
    # Torque should remain constant because isotonic load is set by researcher 
    ISO_Torque = ContractionData['Filtered Torque'].mean()

    return ISO_Velocity, ISO_Torque, ContractionData

def Recovery_Isotonics(Data: pd.DataFrame = None, File: str = None, Subject: str = None, Timepoint: str = None, Recovery_Timepoint: str = None, Graph: bool = False) -> tuple[float, float, float]:
    # Grab subset of data where isotonic contractions were performed 
    # Make sure to include at least 1000 samples (i.e., 100 ms) prior to contraction starting so baseline signals can be obtained
    Rec1_Data = Data[20000:26000]
    Rec2_Data = Data[55000:61000]

    try: 
        Rec1_Data = Isotonic_Contractions(Data = Rec1_Data, Subject = Subject, Timepoint = Timepoint)
    except:
        Rec1_ISO_Velocity = 0
    else:
        Rec1_ISO_Velocity = ((Rec1_Data['Length'].iloc[-1] - Rec1_Data['Length'].iloc[1]) / (Rec1_Data['Time'].iloc[-1] - Rec1_Data['Time'].iloc[1]) * -1)
    
    Rec1_ISO_Torque = Rec1_Data['Filtered Torque'].mean()
    Rec1_ISO_Power = Rec1_ISO_Torque * Rec1_ISO_Velocity


    try: 
        Rec2_Data = Isotonic_Contractions(Data = Rec2_Data, Subject = Subject, Timepoint = Timepoint)
    except:
        Rec2_ISO_Velocity = 0
    else:
        Rec2_ISO_Velocity = ((Rec2_Data['Length'].iloc[-1] - Rec2_Data['Length'].iloc[1]) / (Rec2_Data['Time'].iloc[-1] - Rec2_Data['Time'].iloc[1]) * -1)
    
    Rec2_ISO_Torque = Rec2_Data['Filtered Torque'].mean()
    Rec2_ISO_Power = Rec2_ISO_Torque * Rec2_ISO_Velocity

    if Rec1_ISO_Power > Rec2_ISO_Power:
        ISO_Velocity = Rec1_ISO_Velocity 
        ISO_Torque = Rec1_ISO_Torque
    else:
        ISO_Velocity = Rec2_ISO_Velocity 
        ISO_Torque = Rec2_ISO_Torque

    if Graph == True:
        fig, Recovery_ISOs = plt.subplots(nrows = 1, ncols = 2, figsize = (width:= 6, height:= 4), layout = 'constrained')
        Recovery_ISOs[0].plot(Rec1_Data['Time'], Rec1_Data['Length'])
        Recovery_ISOs[1].plot(Rec2_Data['Time'], Rec2_Data['Length'])
        plt.title(File)
        plt.show()
    # 
 
    return ISO_Velocity, ISO_Torque

def Fatigue(Data: pd.DataFrame = None, Subject: str = None, Timepoint: str = None) -> int:

    # The end of the actual data must be found for fatigue tests
    # The same protocol is used for each subject so there may be several/many extra contractions written into Aurora Scientific protocol
    # If the end of data isn't found, then there could be several thousand rows of N/A values in data file
    DataEnd = Data.index[Data['Length'] == 0][0]
    Data = Data[:DataEnd]

    # Find number of contractions performed during fatigue test
    IndFatigueContractions = []
    Contraction_Start = 20000
    Contraction_End = Contraction_Start + 5000 

    while Contraction_End <= DataEnd:
        Contraction = Data['Filtered Torque'][Contraction_Start:Contraction_End].max()
        Contraction_Start += 20000
        Contraction_End = Contraction_Start + 5000 
        IndFatigueContractions.append(Contraction)

    FatigueContractions: int = len(IndFatigueContractions)
    First_Fatigue_Contraction: float = IndFatigueContractions[0]
    Last_Fatigue_Contraction: float = IndFatigueContractions[-1]
    Percent_Drop: float = ((Last_Fatigue_Contraction - First_Fatigue_Contraction)/First_Fatigue_Contraction) * 100

    return FatigueContractions

def RTD_Analysis(Data: pd.DataFrame = None, PeakTorque: float = None):
    # Contraction onset defined as the point torque exceeds 3 standard deviations of baseline
    BaselineTorque_Mean = Data['Filtered Torque'].iloc[15000:16000].mean()
    BaselineTorque_STDEV = statistics.stdev(Data['Filtered Torque'][15000:16000])
    RTDStart_Criteria = BaselineTorque_Mean + (3 * BaselineTorque_STDEV)

    # Find the sample where stimulations began to be delivered
    StimIndex = Data.index[Data['Stim'] == 1][0]

    # Search for the point torque exceeds defined contraction onset only in data following stimulation onset
    # So we are confident contraction onset will be defined properly
    ContractionStart = Data[StimIndex:].index[Data['Filtered Torque'][StimIndex:] >= RTDStart_Criteria][0]

    # Then earch for first sample where torque exceeds defined percentage of peak torque
    # The goal is to only capture portions of the contraction with steep incline in force production
    # Slope tends to begin to plateau around ~90% peak torque
    ContractionEnd = Data[ContractionStart:].index[Data['Filtered Torque'][ContractionStart:] >= .90 * PeakTorque][0]

    TimetoPeak = (Data['Time'][ContractionEnd] - Data['Time'][ContractionStart])
    RTD = (Data['Filtered Torque'][ContractionEnd] - Data['Filtered Torque'][ContractionStart]) / (Data['Time'][ContractionEnd] - Data['Time'][ContractionStart])

    plt.plot(Data['Time'], Data['Filtered Torque'], color = Black.Hex_Value)
    plt.plot(Data['Time'][ContractionStart:ContractionEnd], Data['Filtered Torque'][ContractionStart:ContractionEnd], color = Firebrick.Hex_Value)
    plt.title(f'Absolute RTD = {RTD:.4f}, Normalized RTD {RTD / PeakTorque:.4f} Time = {TimetoPeak * 1000:.1f} ms')
    plt.show()

    return RTD

def Twitch_Characteristics(Data: pd.DataFrame = None):
    # Find the sample where stimulations began to be delivered
    StimIndex = Data.index[Data['Stim'] == 1][0]

    # Contraction onset defined as the point torque exceeds 3 standard deviations of baseline
    BaselineTorque_Mean = Data['Filtered Torque'][StimIndex - 500:StimIndex - 1].mean()
    BaselineTorque_STDEV = statistics.stdev(Data['Filtered Torque'][StimIndex - 500:StimIndex - 1])
    RTDStart_Criteria = BaselineTorque_Mean + (3 * BaselineTorque_STDEV)
    
    counter = 0 
    for counter, i in enumerate(Data['Filtered Torque'][StimIndex:]):
        if i > RTDStart_Criteria:
            RTD_Start: int = counter
            break
        else:
            counter += 1
    
    PeakTwitch = Data['Filtered Torque'][StimIndex + 10:].max()

    counter = 0 
    for counter, i in enumerate(Data['Filtered Torque']):
        if i == PeakTwitch:
            Peak: int = counter
            break
        else:
            counter += 1
    
    Post_Twitch_Baseline = Data['Filtered Torque'].iloc[-300:].mean()

    counter = 0
    try: 
        for counter, i in enumerate(Data['Filtered Torque'][Peak:]):
            if i < 0.99 * Post_Twitch_Baseline:
                Twitch_End = Peak + counter
                break
            else:
                counter += 1
    except: 
        Twitch_End = Data['Filtered Torque'][Peak + 500]
    
    Half_Relax_Force: float = Data['Filtered Torque'][Peak] / 2

    counter = 0 
    Diffs = []
    for counter, i in enumerate(Data['Filtered Torque'][Peak:Twitch_End]):
        Diff = (i - Half_Relax_Force)**2
        Diffs.append(Diff)
        Min_Diff = min(Diffs)
    Half_Force = Diffs.index(Min_Diff) + Peak

    Half_Relaxation_Time = (Data['Time'][Half_Force] - Data['Time'][Peak]) * 1000
    RTD_Start = StimIndex + RTD_Start
    RTD = (Data['Filtered Torque'][Peak] - Data['Filtered Torque'][RTD_Start]) / (Data['Time'][Peak] - Data['Time'][RTD_Start])

    return Half_Relaxation_Time, Half_Relax_Force, Half_Force, Peak, StimIndex, RTD_Start, RTD, PeakTwitch, Twitch_End

def PLFFD(Data: pd.DataFrame = None, File: str = None, RecoveryTimepoint: str = None, Subject: str = None, Timepoint: str = None, Graph: bool = False) -> tuple[float, float, float, float, float]:
    if Timepoint == 'D120' and Subject in ['VCD_04', 'VCD_05', 'VCD_15', 'VCD_16']:
        Data['Filtered Torque'] = Data['Filtered Torque'] / 10
    
    CT_Start = 29000
    Control_Twitch_df = Data[CT_Start:32000].reset_index(drop = True)

    LowHz_Data = Data[39000:48000].reset_index(drop = True)

    HighHz_Data = Data[49000:58000].reset_index(drop = True)

    PT_Start = 59000
    Potentiated_Twitch_df = Data[PT_Start:62000].reset_index(drop = True)

    try:
        # CT = Control Twitch
        CT_HRT, CT_Half_Relax_Force, CT_Half_Force, CT_Peak, CT_Twitch_Artifact, CT_RTD_Start, CT_RTD, CT_TQ, CT_Twitch_End = Twitch_Characteristics(Control_Twitch_df)
    except: 
        print(Error(f'Failed control twitch in {File}'))
        pass

    try: 
        # LowHz = 10 Hz
        LowHz_HRT, LowHz_Half_Relax_Force, LowHz_Half_Force, LowHz_Peak, LowHz_Twitch_Artifact, LowHz_RTD_Start, LowHz_RTD, LowHz_TQ, LowHz_Twitch_End = Twitch_Characteristics(LowHz_Data)
    except:
        print(Error(f'Failed 10 Hz in {File}'))
        pass

    try:
        # PT = Potentiated Twitch
        PT_HRT, PT_Half_Relax_Force, PT_Half_Force, PT_Peak, PT_Twitch_Artifact, PT_RTD_Start, PT_RTD, PT_TQ, PT_Twitch_End = Twitch_Characteristics(Potentiated_Twitch_df)
    except: 
        print(Error(f'Failed potentiated twitch in {File}'))
        pass

    LowHz_TQ: float = LowHz_Data['Filtered Torque'].max()

    High_Hz_Baseline_Torque = HighHz_Data['Filtered Torque'][0:100].mean()
    HighHz_Data['Filtered Torque'] = HighHz_Data['Filtered Torque'] - High_Hz_Baseline_Torque
    HighHz_TQ: float = HighHz_Data['Filtered Torque'].max()
    Low_High_Hz_Ratio: float = LowHz_TQ / HighHz_TQ

    if Graph == True:
        
        fig, PLFFD_Fig = plt.subplots(figsize = (width:= 5, height:= 3), layout = 'constrained')
        PLFFD_Fig.plot(Data['Time'], Data['Filtered Torque'], color = Charcoal.Hex_Value)
        PLFFD_Fig.set_xlabel('Time (s)')
        PLFFD_Fig.set_ylabel(Constants.Torque)
        PLFFD_Fig.set_title(f'File: {os.path.basename(File)} \n'
                            f'Recovery Timepoint: {RecoveryTimepoint} \n',
                            fontdict = dict(fontsize = 4),
                            pad = 30)

        CT_Inset = fig.add_axes([0.1, 0.5, 0.25, 0.25])
        CT_Inset.plot(Control_Twitch_df['Time'], Control_Twitch_df['Filtered Torque'])
        CT_Inset.annotate('Control twitch RTD start',
                           xy = (Data['Time'][CT_Start + CT_RTD_Start], Data['Filtered Torque'][CT_Start+ CT_RTD_Start]),
                           xycoords = 'data',
                           xytext = (-10, 10),
                           textcoords = 'offset points',
                           arrowprops = dict(arrowstyle = '->'),
                           fontsize = 4)

        CT_Inset.annotate('Control twitch Peak',
                           xy = (Data['Time'][CT_Start + CT_Peak], Data['Filtered Torque'][CT_Start + CT_Peak]),
                           xycoords = 'data',
                           xytext = (10, 10),
                           textcoords = 'offset points',
                           arrowprops = dict(arrowstyle = '->'),
                           fontsize = 4)
        
        CT_Inset.annotate('CT HRT Force',
                           xy = (Data['Time'][CT_Start + CT_Half_Force], Data['Filtered Torque'][CT_Start + CT_Half_Force]),
                           xycoords = 'data',
                           xytext = (10, 10),
                           textcoords = 'offset points',
                           arrowprops = dict(arrowstyle = '->'),
                           fontsize = 4)
        
        PT_Inset = fig.add_axes([0.8, 0.5, 0.25, 0.25])
        PT_Inset.plot(Potentiated_Twitch_df['Time'], Potentiated_Twitch_df['Filtered Torque'])

        PT_Inset.annotate('Potentiated twitch RTD start',
                           xy = (Data['Time'][PT_Start + PT_RTD_Start], Data['Filtered Torque'][PT_Start + PT_RTD_Start]),
                           xycoords = 'data',
                           xytext = (-10, 10),
                           textcoords = 'offset points',
                           arrowprops = dict(arrowstyle = '->'),
                           fontsize = 4)

        PT_Inset.annotate('Potentiated twitch Peak',
                           xy = (Data['Time'][PT_Start + PT_Peak], Data['Filtered Torque'][PT_Start + PT_Peak]),
                           xycoords = 'data',
                           xytext = (10, 10),
                           textcoords = 'offset points',
                           arrowprops = dict(arrowstyle = '->'),
                           fontsize = 4)

        PT_Inset.annotate('PT HRT Force',
                           xy = (Data['Time'][PT_Start + PT_Half_Force], Data['Filtered Torque'][PT_Start + PT_Half_Force]),
                           xycoords = 'data',
                           xytext = (10, 10),
                           textcoords = 'offset points',
                           arrowprops = dict(arrowstyle = '->'),
                           fontsize = 4)

        PLFFD_Fig.text(x = 2, 
                       y = 1,
                       s = f'Control twitch TQ = {CT_TQ:.3f} (mN x m) \n'
                           f'Control twtich RTD = {CT_RTD:.3f} (mN x m / s) \n'
                           f'Control twitch 1/2 relax time = {CT_HRT:.3f} ms',
                        fontdict = dict(fontsize = 4))

        PLFFD_Fig.text(x = 2, 
                       y = 3,
                       s = f'10 Hz TQ = {LowHz_TQ:.3f} (mN x m) \n'
                           f'10 Hz RTD = {LowHz_RTD:.3f} (mN x m / s) \n'
                           f'10 Hz 1/2 relax time = {LowHz_HRT:.3f} ms',
                        fontdict = dict(fontsize = 4))

        PLFFD_Fig.text(x = 7,
                       y = 1,
                       s = f'Potentiated twitch TQ = {PT_TQ:.3f} (mN x m) \n'
                           f'Potentiated twtich RTD = {PT_RTD:.3f} (mN x m / s) \n'
                           f'Potentiated twitch 1/2 relax time = {PT_HRT:.3f} ms',
                       fontdict = dict(fontsize = 4))
     

        # PLFFD_Fig.annotate('Control twitch RTD start',
        #                    xy = (Data['Time'][29000 + CT_RTD_Start], Data['Filtered Torque'][29000 + CT_RTD_Start]),
        #                    xycoords = 'data',
        #                    xytext = (-100,-10),
        #                    textcoords = 'offset points',
        #                    arrowprops = dict(arrowstyle = '->'),
        #                    fontsize = 4)
        
        # PLFFD_Fig.annotate('Potentiated twitch RTD start',
        #                    xy = (Data['Time'][59000 + CT_RTD_Start], Data['Filtered Torque'][59000 + CT_RTD_Start]),
        #                    xycoords = 'data',
        #                    xytext = (-100,-10),
        #                    textcoords = 'offset points',
        #                    arrowprops = dict(arrowstyle = '->'),
        #                    fontsize = 4)
        
        
        plt.show()

    return CT_TQ, CT_HRT, CT_RTD, LowHz_TQ, HighHz_TQ, Low_High_Hz_Ratio, PT_TQ, PT_HRT, PT_RTD

def Means_and_SEM(Timepoint: str, Test: str, Results: pd.DataFrame) -> tuple[dict, dict]:
    Means: dict = {}
    SEM: dict = {}

    if Test == 'Torque-Frequency':
        print(f'Calculating means and SEM for {Test}')
        for Hz in Constants.HzFrequencies:
            if Hz >= 10:
                for column in Results.columns:
                    if f'{Timepoint} {Hz} Hz TQ' == column:
                        Means[Hz] = Results[column].mean()
                        SEM[Hz] = np.std(Results[column]) / np.sqrt(len(Results[column]))
    
    elif Test == 'Torque-Velocity-Power':
        print(f'Calculating means and SEM for {Test}')

        for Percent in ISO_Percents:
            for column in Results.columns:
                if f'{Timepoint} {Percent}% ISO Velocity' == column:
                    Means[Percent] = Results[column].mean()
                if f'{Timepoint} {Percent}% ISO Torque' == column:
                    SEM[Percent] = Results[column].mean()
                    
    elif Test == 'Torque Velocity (TQ)':
        print(f'Calculating means and SEM for {Test}')

        for Percent in ISO_Percents:
            for column in Results.columns:
                if f'{Timepoint} {Percent}% ISO Torque' == column:
                    Means[Percent] = Results[column].mean()
                    SEM[Percent] = np.std(Results[column]) / np.sqrt(len(Results[column]))

    elif Test == 'PLFFD Control Twitch':
        print(f'Calculating means and SEM for {Test}')

        for Recovery_Timepoint in Constants.RecoveryTimepoints:
            if Recovery_Timepoint in {1, 2, 5, 10}:
                Recovery_Timepoint = f'{Recovery_Timepoint} min'
            for column in Results.columns:
                if f'{Timepoint} {Recovery_Timepoint} PLFFD Control Twitch TQ' == column:
                    Means[Recovery_Timepoint] = Results[column].mean()
                    SEM[Recovery_Timepoint] = np.std(Results[column]) / np.sqrt(len(Results[column]))
  
    elif Test == 'Recovery PLFFD 100 Hz':
        print(f'Calculating means and SEM for {Test}')

        for Recovery_Timepoint in Constants.RecoveryTimepoints:
            if Recovery_Timepoint in {1, 2, 5, 10}:
                Recovery_Timepoint = f'{Recovery_Timepoint} min'
            for column in Results.columns:
                if f'{Timepoint} {Recovery_Timepoint} PLFFD 100 Hz TQ' == column:
                    Means[Recovery_Timepoint] = Results[column].mean()
                    SEM[Recovery_Timepoint] = np.std(Results[column]) / np.sqrt(len(Results[column]))
    
    elif Test == 'Recovery Velocity':
        print(f'Calculating means and SEM for {Test}')
        
        for Recovery_Timepoint in Constants.RecoveryTimepoints:
            if Recovery_Timepoint in {1, 2, 5, 10}:
                Recovery_Timepoint = f'{Recovery_Timepoint} min'
            for column in Results.columns:
                if f'{Timepoint} {Recovery_Timepoint} ISO Velocity' == column:
                    Means[Recovery_Timepoint] = Results[column].mean()
                    SEM[Recovery_Timepoint] = np.std(Results[column]) / np.sqrt(len(Results[column]))
    
    elif Test == 'Recovery Torque':
        print(f'Calculating means and SEM for {Test}')
        
        for Recovery_Timepoint in Constants.RecoveryTimepoints:
            if Recovery_Timepoint in {1, 2, 5, 10}:
                Recovery_Timepoint = f'{Recovery_Timepoint} min'
            for column in Results.columns:
                if f'{Timepoint} {Recovery_Timepoint} ISO Torque' == column:
                    Means[Recovery_Timepoint] = Results[column].mean()
                    SEM[Recovery_Timepoint] = np.std(Results[column]) / np.sqrt(len(Results[column]))

    return Means, SEM

