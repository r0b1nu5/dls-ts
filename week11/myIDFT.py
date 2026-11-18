import numpy as np

def my_IDFT(myF):
    N = 2 * len(myF) - 1            # recover odd N from coefficient length
    half = (N + 1) // 2             # positive frequencies stored: k = 1:(N-1)/2

    F = np.empty(N, dtype=complex)

    # invert the scalings
    F[0] = myF[0] * np.sqrt(N)                                  # F[0] = myF[0]/(1/√N)
    F[1:half] = np.sqrt(N / 2) * (myF[1::2] + 1j * myF[2::2])   # Recreating the complex values

    # fill the conjugate-symmetric half: F[N-k] = conj(F[k])
    F[half:] = np.conj(F[N - half:0:-1])

    return np.fft.ifft(F).real


