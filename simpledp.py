import numpy as np
pi=np.pi
ne=2e17; ni=ne; z=1; mu=1; wpe=5.64e4*(ne)**0.5; wpi= 1.32e3* z* mu**-0.5 * ni**0.5
B=215;wce=1.76e7*B
w=2*pi*5e9;k1=1-(wpe/w)**2; k2=1-(wpi/w)**2 + (wpe/wce)**2
r12=abs(k1/k2)**0.5
