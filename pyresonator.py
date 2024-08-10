import numpy as np
from enum import Enum
from scipy.signal import lfilter, butter
import soundfile as sf
import math
import re
import tkinter as tk
from tkinter import filedialog, ttk
import os
import sys

# DEFAULT PARAMS
#---------------------------------------------------

keyNote = "F3"
mode = "MAJOR7"

Q = 25
filterCutoff = 2000
wetDryMix = 0.7

# Data structures representing Notes and Chords
#---------------------------------------------------

# Enum for chords indices
class Chords(Enum):
    SINGLE_NOTE = [0]
    MAJOR = [0, 4, 7]
    MINOR = [0, 3, 7]
    MAJOR7 = [0, 4, 7, 11]
    SUSPENDED4 = [0, 5, 7]

# Enum for relative note index (relative to A)
class Notes(Enum):
    A = 0
    A_SHARP = 1
    B = 2
    C = -9
    C_SHARP = -8
    D = -7
    D_SHARP = -6
    E = -5
    F = -4
    F_SHARP = -3
    G = -2
    G_SHARP = -1

# GUI to select input file
#---------------------------------------------------

# Simple File Explorer for Convenience
def selectFile():
    # Get the directory of the current Python script
    homeDir = os.path.expanduser("~")
    # Create a root window and hide it
    FileExplorerView = tk.Tk()
    FileExplorerView.withdraw()

    # Open file dialog with specified initial directory and file types
    filePath = filedialog.askopenfilename(
        initialdir=homeDir,
        title="Select a .wav file",
        filetypes=(("WAV files", "*.wav"), ("All files", "*.*"))
    )

    # Optional: Close the root window
    FileExplorerView.destroy()

    # Check if the file selection was canceled
    if not filePath:
        sys.exit("File selection was cancelled. Exiting the program.")

    return filePath

inputFile = selectFile()

# GUI to select the input paramaters
#---------------------------------------------------

# Function to update the parameters and print them
def updateParams():
    global keyNote, mode, Q, filterCutoff, wetDryMix

    keyNote = f"{noteVar.get()}{octaveVar.get()}"
    mode = modeVar.get()
    Q = qSlider.get()
    filterCutoff = cutoffSlider.get()
    wetDryMix = dryWetSlider.get()
    
    print(f"KeyNote: {keyNote}")
    print(f"Mode: {mode}")
    print(f"Q: {Q}")
    print(f"Filter Cutoff: {filterCutoff}")
    print(f"Wet/Dry Mix: {wetDryMix}")
    
    ParamSelView.quit()

# Initialize the main window
ParamSelView = tk.Tk()
ParamSelView.title("Audio Parameter Selector")

# Note selector
noteVar = tk.StringVar(value="F")
noteLabel = ttk.Label(ParamSelView, text="Note")
noteLabel.pack(pady=5)
noteSelector = ttk.Combobox(ParamSelView, textvariable=noteVar, values=["A", "B", "C", "D", "E", "F", "G"], state="readonly")
noteSelector.pack()

# Octave selector
octaveVar = tk.StringVar(value="3")
octaveLabel = ttk.Label(ParamSelView, text="Octave")
octaveLabel.pack(pady=5)
octaveSelector = ttk.Combobox(ParamSelView, textvariable=octaveVar, values=[str(i) for i in range(0, 9)], state="readonly")
octaveSelector.pack()

# Mode selector
modeVar = tk.StringVar(value="MAJOR7")
modeLabel = ttk.Label(ParamSelView, text="Mode")
modeLabel.pack(pady=5)
modeSelector = ttk.Combobox(ParamSelView, textvariable=modeVar, values=[
    "SINGLE_NOTE",
    "MAJOR",
    "MINOR",
    "MAJOR7",
    "SUSPENDED4"
], state="readonly")
modeSelector.pack()

# Q slider
qLabel = ttk.Label(ParamSelView, text="Q")
qLabel.pack(pady=5)
qSlider = ttk.Scale(ParamSelView, from_=1.0, to=80.0, orient=tk.HORIZONTAL)
qSlider.set(Q)
qSlider.pack()

# Filter Cutoff slider
cutoffLabel = ttk.Label(ParamSelView, text="Filter Cutoff")
cutoffLabel.pack(pady=5)
cutoffSlider = ttk.Scale(ParamSelView, from_=20.0, to=18000.0, orient=tk.HORIZONTAL)
cutoffSlider.set(filterCutoff)
cutoffSlider.pack()

# Dry/Wet Mix slider
dryWetLabel = ttk.Label(ParamSelView, text="Dry/Wet Mix")
dryWetLabel.pack(pady=5)
dryWetSlider = ttk.Scale(ParamSelView, from_=0.0, to=1.0, orient=tk.HORIZONTAL)
dryWetSlider.set(wetDryMix)
dryWetSlider.pack()

# Submit button
submitButton = ttk.Button(ParamSelView, text="Submit", command=updateParams)
submitButton.pack(pady=20)

# Exit program on close
def onClose():
    sys.exit("Parameter selection was cancelled. Exiting the program.")

# Bind the window close event to the on_close function
ParamSelView.protocol("WM_DELETE_WINDOW", onClose)

# Run the main loop
ParamSelView.mainloop()

# Generate the output path
def generateOutputPath(inputFile, keyNote, mode):
    # Extract directory, filename, and extension
    inputDir = os.path.dirname(inputFile)
    inputName = os.path.splitext(os.path.basename(inputFile))[0]  # Filename without extension
    
    # Generate output filename with the specified format
    outputFilename = f"{inputName}_{keyNote}_{mode}.wav"
    
    # Combine directory and output filename
    outputPath = os.path.join(inputDir, outputFilename)
    
    return outputPath

outputFile = generateOutputPath(inputFile, keyNote, mode)

# Calculate all needed resonant frequencies
#---------------------------------------------------

# Parses a note of format C#5 to it's NoteStep enum and its octave using regular expressions
def parseNoteFromString(NoteString):
    match = re.match(r"([A-Ga-g]#?)(\d+)", NoteString)
    if match:
        note = match.group(1).upper().replace('#', '_SHARP')
        octave = int(match.group(2))
        return note, octave
    else:
        raise ValueError("Invalid note format. Use a format like 'C#5'.")

# Calculates the frequency of a any given Note
def calcFreqFromString(NoteString):
    refFreq = 440.0  # Reference is the frequency of A4
    
    note, octave = parseNoteFromString(NoteString)
    noteIndex = Notes[note].value + (octave - 4) * 12
    
    return calculateFrequency(refFreq, noteIndex)

# Calculate the frequency using the formula: f = f0 * (2 ^ (n / 12))
def calculateFrequency(f0, n):
    return f0 * math.pow(2, n/12)

# Retrieves a list of note indices from a root note and a chord
def getChordOfNote(rootNote:str, chord):
    rootFreq = calcFreqFromString(rootNote)
    relativeIndices = Chords[chord].value
    
    frequencies = [
        calculateFrequency(rootFreq, relIndex) for relIndex in relativeIndices
    ]
    
    return frequencies

# Calculates all of the frequencies
resonantFreqs = getChordOfNote(keyNote, mode)

# Apply effect chain
#---------------------------------------------------

# Resonator filter function
def resonatorFilter(data, sr, freq, Q=5):
    bw = freq / Q
    low = (freq - bw/2) / (sr/2)
    high = (freq + bw/2) / (sr/2)
    b, a = butter(N=2, Wn=[low, high], btype='band')
    return lfilter(b, a, data)

# Read the input file
data, samplerate = sf.read(inputFile)
if data.ndim > 1:
    data = data.mean(axis=1)

# Apply the resonator filter to the data
processedData = np.zeros_like(data)
for freq in resonantFreqs:
    processedData += resonatorFilter(data, samplerate, freq, Q)

# Normalize the output
maxVal = np.max(np.abs(processedData))
if maxVal > 0:
    processedData = processedData / maxVal

# Dry/Wet Mix
processedData = (wetDryMix * processedData) + ((1 - wetDryMix) * data)

# Add a low-pass filter to smooth the sound
def lowPassFilter(data, sr, cutoff):
    norm_cutoff = cutoff / (sr / 2)
    b, a = butter(N=2, Wn=norm_cutoff, btype='low')
    return lfilter(b, a, data)

outputData = lowPassFilter(processedData, samplerate, filterCutoff)

# Generate Output File
#---------------------------------------------------

# Write the processed data to an output file
sf.write(outputFile, outputData, samplerate)