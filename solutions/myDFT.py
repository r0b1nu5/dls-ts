import numpy as np
import matplotlib.pyplot as plt

def my_DFT(X):
    N = len(X)

    if N%2 == 0:
        X = X[:-1]
        N -= 1
    
    F = np.fft.fft(X)

    myF = np.array([F[0].real/np.sqrt(N),])

    for k in range(1,int((N+1)/2)):
        myF = np.append(myF, F[k].real*np.sqrt(2/N))
        myF = np.append(myF, F[k].imag*np.sqrt(2/N))

    return myF

def my_DFT_amplitude(X):
    N = len(X)

    if N%2 == 0:
        X = X[:-1]
        N -= 1
    
    F = np.fft.fft(X)

    myF = np.array([F[0].real/np.sqrt(N),])

    for k in range(1,int((N+1)/2)):
        myF = np.append(myF, np.sqrt((F[k].real*np.sqrt(2/N))**2 + (F[k].imag*np.sqrt(2/N))**2))

    return myF

def plot_DFT(myF):
    N = len(myF)

    c = myF[1::2]
    s = myF[2::2]
    x = [myF[0]] + [np.sqrt(c[i]**2 + s[i]**2) for i in range((N-1)/2)]

    fig, ax = plt.subplots(figsize=(7,3))

    ax.plot(range((N-1)/2), x, 'o', color='C2')
    ax.set_xlabel('k')
    ax.set_ylabel('magnitude')

    plt.tight_layout()
    plt.show()

    