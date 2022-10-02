# Imports
import numpy as np

# Constantes

# Vecteur largeur des portes
L = [125, 100, 75, 75, 75]
# Nombre de pieces
n = len(L)
# Vecteur repartition initiale des mouches
M = [100, 100, 100, 0, 0]
# Probabilite qune mouche ne change pas de piece
reste_p = 0
# Probabilite  qune mouche change de piece
bouge_p = 1 - reste_p

# Intialisation

def initMatriceTransition():
    T = np.zeros((n, n))
    for i in range(n):
        T[i, i] = reste_p
        T[i, (i + 1) % n] = bouge_p * L[i] / (L[i] + L[(i - 1) % n])
        T[i, (i - 1) % n] = bouge_p * L[(i - 1) % n] / (L[i] + L[(i - 1) % n])
        
    return T

# Algo

def puissanceInverse(A):
    precision = 1e-6
    z = 0.99
    x = M
    k = 0
    resolution_aboutit = True
    q = 0
    norme = 0
    I = np.identity(n)
    R = A - z * I
    while resolution_aboutit:
        k += 1
        if k > 5:
            R = A - q * I
        try:
            y = np.linalg.solve(R, x)
            old_norme = norme
            norme = np.linalg.norm(y)
            print(norme, old_norme)
            if abs(norme - old_norme) > precision:
                x = y / norme
                q = np.dot(x, A * x)
            else:
                resolution_aboutit = False
        except np.linalg.LinAlgError:
            resolution_aboutit = False
            
    return q, x
        
    
    
A = initMatriceTransition()
q, x = puissanceInverse(A)

print(q, x)