#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  8 15:05:13 2017

@author: wangdayuan
"""

import array
import contextlib
import wave
import numpy as np
import matplotlib.pyplot as plt
from cs591Utilities import *
from scipy.interpolate import interp1d
from scipy import signal
from cmath import rect
from cmath import polar
from scipy import *
from pylab import *
from scipy.io import wavfile


def P2R(radii, angles):
    return radii * np.exp(1j*angles)

'''
P1_A = readWaveFile("Bach.Brandenburg.2.3.wav")

P1_A_polar_spectrum = polarFFT(P1_A)
P1_A_newspectrum = []
for i in range(len(P1_A_polar_spectrum)):
    P1_A_newspectrum += [P2R(P1_A_polar_spectrum[i][0], 0.0)]
P1_A_newSignal = np.fft.irfft(P1_A_newspectrum)
displaySignal(P1_A_newSignal)
writeWaveFile('newBach.Brandenburg.2.3.wav', P1_A_newSignal)
'''

def Problem1(X):
    Spectrum = np.fft.fft(X)
    for i in range(len(Spectrum)):
        Spectrum[i] = np.abs(Spectrum[i])
    result = np.fft.ifft(Spectrum)
    return result

def Problem1b(X):
    Spectrum = np.fft.fft(X)
    for i in range(len(Spectrum)):
        Spectrum[i] = (-1j) * Spectrum[i].imag + Spectrum[i].real
    result = np.fft.ifft(Spectrum)
    return result

P1 = readWaveFile("Bach.Brandenburg.2.3.wav")
P1_f = Problem1(P1)
displaySignal(P1_f)
writeWaveFile('newBach.Brandenburg.2.3.wav', P1_f)

P1_2 = readWaveFile("Bach.Brandenburg.2.3.wav")
P1_f_2 = Problem1b(P1_2)
displaySignal(P1_f_2)
writeWaveFile('newBach_2.Brandenburg.2.3.wav', P1_f)

def timeStretch(X,P):
    skip = 1.0/P
    lo = 0 
    hi = 100
    result = []
    while( hi <= len(X) ):
        a = lo + 5
        b = hi -5
        F = range(lo,hi)
        f = interp1d(F, X[(lo):(hi)], kind="cubic")
        F = np.arange(lo+5,hi-5, skip)
        result += [f(x) for x in F]
        lo = lo + 90
        if hi + 90 > len(X) and hi != len(X):
            hi = len(X)
        hi = hi + 90
    return result


def timeStretch_1(X, P):
    Y = list(X)
    #print(len(X))
    resampleRate = int(len(Y) * P)
    Y = signal.resample(Y, resampleRate)
    print(len(Y))
    return Y

def PhaseVocoder(X, P):
    N = 2048       # window size for FFT
    H = N//4
    sr = 44100
    tscale = 1.0/P
    phi  = zeros(N)
    out = zeros(N, dtype=complex)
    Y = zeros(int(len(X)/tscale+N))
    amp = max(X)
    win = hanning(N)
    p = 0
    pp = 0
    while(p < len(X)-(N+H)):
        spec1 =  fft(win*X[p:p+N])
        spec2 =  fft(win*X[p+H:p+N+H])
        phi += (angle(spec2) - angle(spec1))
        for i in range(len(phi)):
            while(phi[i] < -pi): 
                phi[i] += 2*pi
            while(phi[i] >= pi): 
                phi[i] -= 2*pi
        out.real, out.imag = cos(phi), sin(phi)
        Y[pp:pp+N] += (win*ifft(spec2*out)).real
        pp += H
        p += int(H*tscale)
    return Y

def pitchStretch(X, P):
    Y = list(X)
    p = 1/P
    #result = PhaseVocoder(Y, P)
    result = timeStretch_1(Y, p)
    #p = 1/P
    result = PhaseVocoder(result, P)
    return result
    
    
    
'''   
P3 = readWaveFile("BluesGuitar.wav")
displaySignal(P3)
P3_a = timeStretch(P3, 2.3)
writeWaveFile("test_3.wav", P3_a)
displaySignal(P3_a)
'''

P3 = readWaveFile("BluesGuitar.wav")
displaySignal(P3)
P3_a = timeStretch_1(P3, 2.3)
writeWaveFile("test_3_1.wav", P3_a)
displaySignal(P3_a)

'''
P5 = readWaveFile("Bach.Brandenburg.2.3.wav")
P3_r = pitchStretch(P5, 1/2)
writeWaveFile("test_5.wav", P3_r)
displaySignal(P3_r)
analyzeSpectrum(P3_r)

P5_1 = readWaveFile("Bach.Brandenburg.2.3.wav")
P3_f = pitchStretch(P5_1, 1/(2**(17/12)) )
writeWaveFile("test_5_1.wav", P3_f)
displaySignal(P3_f)
analyzeSpectrum(P3_f)
'''