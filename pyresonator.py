import numpy as np
from enum import Enum
from scipy.signal import lfilter, butter
import soundfile as sf
import math
import re

# DEFINE ALL THE PARAMETERS HERE:
#---------------------------------

inputFile = 'input.wav'
outputFile = 'output.wav'

keyNote = "F3"
mode = "MAJOR7"

Q = 25
filterCutoff = 2000
wetDryMix = 0.7

#---------------------------------

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



# Write the processed data to an output file
sf.write(outputFile, outputData, samplerate)