import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from math import ceil

# Set parameters
## General parameters
t0 = 4  # Start time of analysis (adjust as needed)
t1 = 20  # End time of analysis (adjust as needed)
threshold = 0  # Threshold for peak detection (keep at 0 since i take the mean of the signal)

## Parrameters for FFT (not really needed for wave elevation data)
Tmin = 3  # Minimum period of the signal
fmax = 1 / Tmin  # Hz
timeout = 0.05  # Time out in seconds

# 1. Read and process the data
def read_and_preprocess_data(path, gauge_indices = [0,1]):
    with open(path, 'r') as f:
        first_line = f.readline().strip()
    gauge_locations = first_line.split(';')[2:]
    new_column_names = ['Time [s]'] + [f'{gauge_locations[i]}m' for i in gauge_indices]
    df = pd.read_csv(path, skiprows=4, delimiter=';', names=new_column_names, usecols=[1] + [i + 2 for i in gauge_indices])
    # Print the head of the dataframe
    print(df.head())

    df.iloc[:, 1:] = df.iloc[:, 1:].astype(float)
    df = df[(df['Time [s]'] >= t0) & (df['Time [s]'] <= t1)]
    return df


# 2. Define FFT for filtering out high and low frequencies (not really needed for wave elevation data)
def fft_filter(data):
    for column in data.columns[1:]:
        fft = np.fft.fft(data[column])
        freq = np.fft.fftfreq(len(data[column]), d=timeout)
        fft[(freq > 2 * fmax * timeout) | (freq < -2 * fmax * timeout)] = 0
        data[column] = np.fft.ifft(fft).real
    return data

# 3. Subtract the mean value of each column from the column
def get_eta(df):
    for col in df.columns[1:]:
        df[col] = df[col] - df[col].mean()
    return df

# 4. Get the peaks of the signal using upward crossing (main function here)
def find_significant_peaks(data):
    results = {}
    for column_name in data.columns[1:]:
        upward_crossings = []
        peaks = []
        troughs = []

        for i in range(1, len(data)):
            if data[column_name].iloc[i - 1] < threshold and data[column_name].iloc[i] >= threshold:
                upward_crossings.append(i)

        for i in range(len(upward_crossings) - 1):
            start = upward_crossings[i]
            end = upward_crossings[i + 1]
            peak_index = np.argmax(data[column_name][start:end]) + start
            trough_index = np.argmin(data[column_name][start:end]) + start
            peaks.append(peak_index)
            troughs.append(trough_index)
        peaks, troughs = crop_lists(peaks, troughs)
        results[column_name + '_peaks'] = peaks
        results[column_name + '_troughs'] = troughs
    return results

# Helper function to crop lists (making sure peaks and troughs are of the same length)
def crop_lists(list1, list2):
    if len(list1) > len(list2):
        list1 = list1[:len(list2)]
    elif len(list2) > len(list1):
        list2 = list2[:len(list1)]
    return list1, list2

# 5. Get the heights of the waves (difference between peaks and troughs)
def get_wave_heights(data, peaks_and_troughs):
    results = {}
    for column_name in data.columns[1:]:
        peaks = peaks_and_troughs[column_name + '_peaks']
        troughs = peaks_and_troughs[column_name + '_troughs']
        peak_heights = [data[column_name].iloc[peak] for peak in peaks]
        trough_heights = [data[column_name].iloc[trough] for trough in troughs]
        wave_heights = [peak_height - trough_height for peak_height, trough_height in zip(peak_heights, trough_heights)]
        results[column_name + '_wave_heights'] = sorted(wave_heights, reverse=True)
    return results


# 6. Get significant wave height (average of the top third of the waves; it should not matter for regular waves, but the code is written to handle irregular waves as well). 
def get_significant_wave_height(wave_heights):
    results = {}
    for column_name, heights in wave_heights.items():
        top_third = heights[:ceil(len(heights)/3)]
        significant_wave_height = np.mean(top_third) if top_third else 0
        results[column_name + '_significant_wave_height'] = significant_wave_height
    return np.mean(list(results.values()))

# 7. Plot the data with peaks and troughs (optional, but useful to keep your plots for each iteration)
def plot_data_with_peaks(df, peaks_and_troughs, plot=True, directory='Plot_Directory'):
    if plot:
        if not os.path.exists(directory):
            os.makedirs(directory)
        existing_files = [f for f in os.listdir(directory) if f.endswith('.png')]
        if existing_files:
            max_num = max([int(f.split('_')[-1].split('.')[0]) for f in existing_files])
        else:
            max_num = 0
        file_name = f'iteration_{max_num + 1}.png'
        full_path = os.path.join(directory, file_name)

        fig, axs = plt.subplots(len(df.columns[1:]), 1, figsize=(8, 6 * len(df.columns[1:])))
        for i, column_name in enumerate(df.columns[1:]):
            axs[i].plot(df['Time [s]'], df[column_name])
            axs[i].scatter(df['Time [s]'].iloc[peaks_and_troughs[column_name + '_peaks']], 
                           df[column_name].iloc[peaks_and_troughs[column_name + '_peaks']], color='red', label='Peaks')
            axs[i].scatter(df['Time [s]'].iloc[peaks_and_troughs[column_name + '_troughs']], 
                           df[column_name].iloc[peaks_and_troughs[column_name + '_troughs']], color='green', label='Troughs')
            axs[i].set_xlabel('Time [s]')
            axs[i].set_ylabel('$\eta$ [m]')
            axs[i].set_title('Gauge at ' + column_name)
            axs[i].legend()
        plt.tight_layout()
        plt.savefig(full_path)
        plt.close()
        print(f"Plot saved as {full_path}")


# Main function to process data
def get_wave_height(path, gauge_indices = [0, 1], apply_fft=False, plot=True):
    data = read_and_preprocess_data(path, gauge_indices)
    if apply_fft:
        data = fft_filter(data)
    data = get_eta(data)
    peaks_and_troughs = find_significant_peaks(data)
    wave_heights = get_wave_heights(data, peaks_and_troughs)
    significant_wave_height = get_significant_wave_height(wave_heights)
    plot_data_with_peaks(data, peaks_and_troughs, plot)
    return significant_wave_height

if __name__ == "__main__":
    # Example usage
    file_path = 'MeasureSWL_Elevation.csv'
    gauge_indices = [0, 1]  # Specify which gauges to use
    result = get_wave_height(file_path, gauge_indices = [0, 1])
    print("Significant Wave Height:", result)