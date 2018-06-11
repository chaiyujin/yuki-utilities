from python_speech_features import mfcc, fbank
from python_speech_features.sigproc import framesig, magspec
import scipy.io.wavfile as wav
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.signal as sig
from scipy.fftpack import fft, ifft


def lpc_ref(signal, order):
    """Compute the Linear Prediction Coefficients.
    Return the order + 1 LPC coefficients for the signal. c = lpc(x, k) will
    find the k+1 coefficients of a k order linear filter:
      xp[n] = -c[1] * x[n-2] - ... - c[k-1] * x[n-k-1]
    Such as the sum of the squared-error e[i] = xp[i] - x[i] is minimized.
    Parameters
    ----------
    signal: array_like
        input signal
    order : int
        LPC order (the output will have order + 1 items)
    Notes
    ----
    This is just for reference, as it is using the direct inversion of the
    toeplitz matrix, which is really slow"""
    if signal.ndim > 1:
        raise ValueError("Array of rank > 1 not supported yet")
    if order > signal.size:
        raise ValueError("Input signal must have a lenght >= lpc order")

    if order > 0:
        p = order + 1
        r = np.zeros(p, signal.dtype)
        # Number of non zero values in autocorrelation one needs for p LPC
        # coefficients
        nx = np.min([p, signal.size])
        x = np.correlate(signal, signal, 'full')
        r[:nx] = x[signal.size-1:signal.size+order]
        phi = np.dot(sp.linalg.inv(sp.linalg.toeplitz(r[:-1])), -r[1:])
        return np.concatenate(([1.], phi))
    else:
        return np.ones(1, dtype=signal.dtype)


def nextpow2(n):
    """Return the next power of 2 such as 2^p >= n.
    Notes
    -----
    Infinite and nan are left untouched, negative values are not allowed."""
    if np.any(n < 0):
        raise ValueError("n should be > 0")

    if np.isscalar(n):
        f, p = np.frexp(n)
        if f == 0.5:
            return p-1
        elif np.isfinite(f):
            return p
        else:
            return f
    else:
        f, p = np.frexp(n)
        res = f
        bet = np.isfinite(f)
        exa = (f == 0.5)
        res[bet] = p[bet]
        res[exa] = p[exa] - 1
        return res


def _acorr_last_axis(x, nfft, maxlag):
    a = np.real(ifft(np.abs(fft(x, n=nfft) ** 2)))
    print("fft", np.abs(fft(x, n=nfft) ** 2))
    print("ifft", a)
    return a[..., :maxlag+1] / x.shape[-1]


def acorr_lpc(x, axis=-1):
    """Compute autocorrelation of x along the given axis.
    This compute the biased autocorrelation estimator (divided by the size of
    input signal)
    Notes
    -----
        The reason why we do not use acorr directly is for speed issue."""
    if not np.isrealobj(x):
        raise ValueError("Complex input not supported yet")

    maxlag = x.shape[axis]
    nfft = 2 ** nextpow2(2 * maxlag - 1)
    if axis != -1:
        x = np.swapaxes(x, -1, axis)
    a = _acorr_last_axis(x, nfft, maxlag)
    if axis != -1:
        a = np.swapaxes(a, -1, axis)
    return a


def time_ticks(time_len, bins, ticks, include_end=True):
    if not include_end:
        bins -= 1
    locs = np.float32(np.linspace(0, bins, ticks))
    return locs, ["%.02f" % l for l in ((locs * time_len / bins))]


def draw_raw_wave(signal, time_len):
    mid = int(len(signal) / 2)
    signal = signal[mid: mid + 256]
    # signal *= np.hamming(len(signal))
    plt.plot(signal)
    plt.xlabel("time (s)")
    plt.ylabel("sample")
    plt.xlim([0, len(signal) - 1])
    plt.xticks(*time_ticks(time_len, len(signal), 5, False))


rate, audio = wav.read("../asset/test1.wav")
audio = audio[:, 0]
# audio = audio / 32768.0

print("-------------lpc--------------------")
audio_frames = framesig(audio, rate * 0.025, rate * 0.010)
print(audio_frames[0].shape)
# audio_frames[0] = np.zeros_like(audio_frames[0])
# audio_frames[0][0] = 1e-8
print(lpc_ref(audio_frames[0], 12))
print(lpc_ref(audio_frames[10], 12))
the_frame = audio_frames[0: 1]

A = lpc_ref(the_frame[0], 12)
B = 1
w, h = sig.freqz(B, A, 512)
w = rate / 2 * w / (2 * np.pi)
plt.plot(w, 20 * np.log10(np.abs(h)))
print('A', A)
print('w', w[:10])
print('h', 20 * np.log10(np.abs(h))[:10])
spec = 20 * np.log10(magspec(the_frame[:1], 512)[0])
spec -= spec.mean()
xx = np.linspace(0, 2000, len(spec))
plt.plot(xx, spec)

A = [3.2, 2, 1]
rts = np.roots(A)
# formants
rts = np.roots(A)
# rts = [r for r in rts if np.imag(r) >= 0]
angz = np.arctan2(np.imag(rts), np.real(rts))
# Get frequencies.
frqs = angz * (rate / 2 / (2 * np.pi))
bw = -0.5 * (rate / 2 / (2 * (np.pi))) * np.log(np.abs(rts))
formants = [frqs[i] for i in range(len(rts))
            if np.imag(rts[i]) >= 0 and frqs[i] > 90 and bw[i] < 400]
print('----')
print(A)
print(rts)
print(angz)
print(frqs)
print(bw)
print(formants)
print('----')

for f in formants:
    xx = []
    yy = []
    xx.append(f)
    xx.append(f)
    yy.append(-30)
    yy.append(30)
    plt.plot(xx, yy)
plt.show()

print(np.exp(-1j * 1))

quit()

print("-------------mfcc-------------------")
feat, energy = fbank(audio, rate)
feat = feat[:64]
feat = mfcc(audio, rate)[:64]
print(feat.shape)


def draw_mfcc_energy(feat, time_len, colorbar=True):
    plt.imshow(np.transpose(feat), origin="lower",
               aspect="auto", cmap="jet", interpolation="none")
    if colorbar:
        plt.colorbar()
    plt.xlabel("time (s)")
    plt.xlim([0, len(feat) - 1])
    plt.xticks(*time_ticks(time_len, len(feat), 5, False))


# feat = 20 * np.log10(np.abs(feat) + 1e-8)
print(feat.min(), feat.max())

draw_mfcc_energy(feat, 64)
plt.show()
