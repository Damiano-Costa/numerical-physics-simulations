import matplotlib.pyplot as plt
import numpy as np
#è necessario che il file si chiami staticfluct.txt
x,y=np.loadtxt("staticfluct.txt",unpack=True)
plt.plot(x,y,'o',label="Valori per D")
M=np.mean(y)
sig=np.std(y,ddof=1)
plt.grid()
plt.title('fluttuazioni statistiche')
plt.xlabel('Storie')
plt.ylabel('D(t) al tempo finale')
plt.axhline(M,color='red',linestyle='-',label="mean")
plt.axhline(M-3*sig,color='green',linestyle='-',label="3$\sigma$")
plt.axhline(M+3*sig,color='green',linestyle='-',label="3$\sigma$")
plt.legend()
plt.savefig('StaicFluct.png')
plt.show()
