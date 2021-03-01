import numpy as np
import matplotlib.pyplot as plt

year=np.array([1978,1980,1982,1984,1986,1988,1990,1992])
toxin=np.array([12,12.7,13,15.2,18.2,19.8,24.1,28.1])
Allyear = np.linspace( 1977.0 , 1997.0 , 201 )

gamma_multiply_yi=np.zeros(len(year))
for xi,x in enumerate(year):
    gamma_multiply_yi[xi]=toxin[xi]/(np.prod(np.ones(len(year)-1)*year[xi]-np.delete(year,[xi])))

Y=np.array([])
for Xi,X in enumerate(Allyear):
    Temp=np.array([])
    for xi,x in enumerate(gamma_multiply_yi):
        Temp=np.append(Temp,np.prod(np.ones(len(year)-1)*X-np.delete(year,[xi]))*gamma_multiply_yi[xi])
    Y=np.append(Y,np.sum(Temp))
    

plt.plot(Allyear,Y)
plt.plot(year,toxin,color='blue',marker='o',linestyle='')
plt.plot(Allyear[np.where(Allyear==1994)],Y[np.where(Allyear==1994)],color='red',marker='o',linestyle='')
plt.xlim(1977,1995)
plt.ylim(-40,29)
plt.xlabel('Year')
plt.ylabel('Toxin concentration')
plt.title('Lagrange polynomial')
plt.grid()
plt.savefig("Lagrange polynomial_1994",dpi=300)
plt.show()

#%%
import numpy as np
import matplotlib.pyplot as plt

year=np.array([1978,1980,1986,1988,1990,1992])
toxin=np.array([12,12.7,18.2,19.8,24.1,28.1])
Allyear = np.linspace( 1977.0 , 1997.0 , 201 )

gamma_multiply_yi=np.zeros(len(year))
for xi,x in enumerate(year):
    gamma_multiply_yi[xi]=toxin[xi]/(np.prod(np.ones(len(year)-1)*year[xi]-np.delete(year,[xi])))

Y=np.array([])
for Xi,X in enumerate(Allyear):
    Temp=np.array([])
    for xi,x in enumerate(gamma_multiply_yi):
        Temp=np.append(Temp,np.prod(np.ones(len(year)-1)*X-np.delete(year,[xi]))*gamma_multiply_yi[xi])
    Y=np.append(Y,np.sum(Temp))
    

plt.plot(Allyear,Y)
plt.plot(year,toxin,color='blue',marker='o',linestyle='')
plt.plot(Allyear[np.where((Allyear==1994) | (Allyear==1984) | (Allyear==1982))],Y[np.where((Allyear==1994) | (Allyear==1984) | (Allyear==1982))],color='red',marker='o',linestyle='')
plt.xlim(1977,1995)
plt.ylim(-40,29)
plt.xlabel('Year')
plt.ylabel('Toxin concentration')
plt.title('Lagrange polynomial')
plt.grid()
plt.savefig("Lagrange polynomial_1984_1982",dpi=300)
plt.show()

#%%
from scipy.interpolate import CubicSpline
import numpy as np
import matplotlib.pyplot as plt
year=np.array([1978,1980,1986,1988,1990,1992])
toxin=np.array([12,12.7,18.2,19.8,24.1,28.1])
c=CubicSpline(year, toxin, axis=0, bc_type='not-a-knot', extrapolate=None)
Allyear = np.linspace( 1977.0 , 1997.0 , 201 )

plt.plot(Allyear,c(Allyear))
plt.plot(year,toxin,color='blue',marker='o',linestyle='')