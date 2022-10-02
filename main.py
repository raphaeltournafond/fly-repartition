# Imports
import decimal
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Constantes

# Vecteur largeur des portes
L = np.array([125, 100, 75, 75, 75])
# Nombre de pieces
n = len(L)
# Vecteur repartition initiale des mouches
M = np.array([100, 100, 100, 0, 0])
# Probabilite qune mouche ne change pas de piece
reste_p = 0
# Probabilite  qune mouche change de piece
bouge_p = 1 - reste_p
# Valeur initiale z
z = 0.78
# Type d'affichage des lignes
courbe = True

# Intialisation

def initMatriceTransition():
    T = np.zeros((n, n))
    for i in range(n):
        T[i, i] = reste_p
        T[(i + 1) % n, i] = bouge_p * L[i] / (L[i] + L[(i - 1) % n])
        T[(i - 1) % n, i] = bouge_p * L[(i - 1) % n] / (L[i] + L[(i - 1) % n])
        
    return T

def puissanceInverse(A, z):
    precision = 1e-6
    r = abs(decimal.Decimal(str(precision)).as_tuple().exponent) 
    I = np.identity(n)
    try:
        R = np.linalg.inv(A - z * I)
        q = z
        x = M
        convergence = False
        evolution = dict()
        evolution[q] = x / np.linalg.norm(x)
        while not convergence:
            y = np.dot(R, x)
            x = y / np.linalg.norm(y)
            prec_q = q
            q = np.dot(x, np.dot(A, x))
            evolution[q] = np.around(x, r)
            cauchy = np.linalg.norm(q - prec_q)
            convergence = cauchy <= precision
        return q, x, evolution
    except np.linalg.LinAlgError:
        print(z, 'est valeur propre')
        raise
    
def numerise(x):
    return np.sum(M) / np.sum(x) * x

# function to add value labels
def addlabels(x,y):
    for i in range(len(x)):
        plt.text(i + 1, y[i], y[i], ha = 'center')
        
def construitStationaire(q, x):
    values = numerise(x)
    rooms = np.arange(1, n + 1)
    plt.bar(rooms, values)
    addlabels(rooms, np.around(values, 4))
    plt.title('Répartition stationaire des mouches par piece\n' +
              'init = ' + str(M) + ' n = ' + str(np.sum(M)) + ' z = ' + str(z))
    plt.xlabel("Piece")
    plt.ylabel("Nombre de mouches")
    
def construitEvolution(evolution):
    plt.figure()
    legende = []
    for i in range(n):
        evolution_salle_i = []
        for q in evolution:
            evolution_salle_i.append(numerise(evolution[q])[i])
        if courbe:
            cubic_interploation_model = interp1d(np.arange(len(evolution)), evolution_salle_i, kind = "cubic")
            # Plotting the Graph
            X = np.linspace(0, len(evolution) - 1, 500)
            Y = cubic_interploation_model(X)
            plt.plot(X, Y)
        else:
            plt.plot(evolution_salle_i)
        legende.append('Piece ' + str(i + 1))
    plt.legend(legende)
    plt.title('Evolution du nombre de mouche par piece\n' +
              'init = ' + str(M) + ' n = ' + str(np.sum(M)) + ' z = ' + str(z))
    plt.xlabel('Itération')
    plt.ylabel('Nombre de mouches')
    
    
def afficheResultats(q, x, evolution):
    print('\nRegime stationaire :\n', 'Le vecteur propre associé à la valeur propre', q, 'est :\n', x)
    print('\nEn évolution :')
    for q in evolution:
        x = evolution[q]
        print('q =', q, '=> x =', numerise(x))
    construitStationaire(q, x)
    construitEvolution(evolution)
    plt.show()

        
A = initMatriceTransition()
print('Matrice de transition A :\n', A)

try:
    q, x, evolution = puissanceInverse(A, z)
    afficheResultats(q, x, evolution)
except np.linalg.LinAlgError:
    print('le z choisis est déjà une valeur propre')