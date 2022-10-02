import numpy as np
import matplotlib.pyplot as plt 
from math import *
from random import *

## Constantes
L = np.array ([125, 100, 75, 75, 75])
N = 5
p = 0
q = 1-p
NombreMouches = np.zeros(N)
NombreMouches[0] = 100
NombreMouches[1] = 100
NombreMouches[2] = 100

A = np.zeros((N,N))

## Saisie de la matrice de transition
for i in range (N):
    A[(i+1)%N,i] = (q*L[i])/(L[i]+L[(i-1)%N])
    A[(i-1)%N,i] = (q*L[(i-1)%N])/(L[i]+L[(i-1)%N])
    A[i,i] = p 
 
print("Voici la matrice transition : ")  
print (A)                           ## Affichage de la matrice

## Fonction qui permet de calculer une somme
def somme(x):
        s = 0
        for i in range(N):
                s+=x[i]
        return s
    
s = somme(NombreMouches)
x = np.linalg.norm(NombreMouches)

## Application de la méthode des puissances inverses
def methodePuissanceInverse(A):
        pres = 10**(-6)
        z = 0.99
        (n,m) = np.shape(A)
        I = np.zeros((n,n))
        
        for i in range(n):          ## Saisie de la matrice identité
                I[i][i] = 1
        X = np.zeros(n)
        evol = [[]]                 ## Enregistrer l'évolution des proportions de mouches
        for i in range(N-1):
            evol.append([])
        for i in range(n):
                X[i] = NombreMouches[i]
                proportion = X[i]/somme(NombreMouches)
                evol[i].append(proportion)

        R = np.linalg.inv(A-z*I)    ## Calcul de l'inverse de la matrice A-z*I

        convergence = False
        
        print(X)
        print(R)
        
        
        
        while not(convergence):
                Y = np.dot(R,X)
                y = np.linalg.norm(Y)
                X = 1/y*Y
                Ax = np.dot(A,X)
                q = np.dot(X,Ax)
                norme = np.linalg.norm(Ax-q*X)
                convergence = (norme<=pres)
                
                for i in range (len(evol)):
                    proportion = X[i]
                    evol[i].append(proportion)
                
        return(q,X,evol)
    


q,X,evol = methodePuissanceInverse(A)
    
print(q, X)
s2 = somme(X)
X = s/s2*X

print(X)

print("Le vecteur propre associé à la valeur 1 : ")
print(X)                                ## Affichage du vecteur propre 

## Diagramme qui illustre le nombre de mouches dans les salles

abs = [i for i in range (1,np.shape(X)[0]+1)]
plt.bar(abs,X)
plt.ylim(0,100)
plt.xlabel("Pièces")
plt.ylabel("Nombre de mouches")
plt.title("Répartition des mouches dans les salles")

## Affichage de l'evolution de la repartition des mouches

plt.figure()
legend=[]
for i, liste in enumerate(evol):
    plt.plot(liste)
    print(legend)
    legend.append(f'salle{i+1}')
plt.legend(legend)
plt.xlabel("Temps")
plt.ylabel("Proportion des mouches ")
plt.title("Evolution de la répartition des mouches dans chaque salle")
plt.show()