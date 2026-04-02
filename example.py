import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

def simular_injecao_co2():
    r_w = 0.1
    r_e = 2000.0
    h = 50.0
    phi = 0.25
    k = 100e-15
    mu = 0.06e-3
    ct = 5e-10
    q_inj = 0.02
    P_init = 15e6

    N = 150
    dr = (r_e - r_w) / N
    r = np.linspace(r_w + dr/2, r_e - dr/2, N)
    r_faces = np.linspace(r_w, r_e, N + 1)

    dias_simulacao = 365
    dt = 24 * 3600
    num_steps = int((dias_simulacao * 24 * 3600) / dt)

    T = np.zeros(N + 1)
    T[1:-1] = (2 * np.pi * k * h / mu) * r_faces[1:-1] / dr
    T[0] = 0 
    T[-1] = 0 

    V = np.pi * ((r + dr/2)**2 - (r - dr/2)**2) * h
    C = V * phi * ct / dt

    diag_principal = C + T[:-1] + T[1:]
    diag_inferior = -T[1:-1]
    diag_superior = -T[1:-1]
    
    A = diags([diag_inferior, diag_principal, diag_superior], [-1, 0, 1], format='csr')

    P = np.ones(N) * P_init
    resultados_P = {}
    dias_para_salvar = [1, 30, 90, 180, 365]

    for step in range(1, num_steps + 1):
        B = C * P
        B[0] += q_inj
        
        P = spsolve(A, B)
        
        if step in dias_para_salvar:
            resultados_P[step] = P.copy()

    plotar_resultados(r, resultados_P, P_init)

def plotar_resultados(r, resultados_P, P_init):
    plt.figure(figsize=(10, 6))
    plt.axhline(y=P_init / 1e6, color='black', linestyle='--', label='Pressão Inicial')
    
    cores = plt.cm.plasma(np.linspace(0, 0.8, len(resultados_P)))
    
    for (dia, P), cor in zip(resultados_P.items(), cores):
        plt.plot(r, P / 1e6, label=f'Dia {dia}', color=cor)

    plt.title('Perfil de Pressão Radial - Injeção de CO2')
    plt.xlabel('Raio (m)')
    plt.ylabel('Pressão (MPa)')
    plt.xlim(0, 500) 
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()

if __name__ == "__main__":
    simular_injecao_co2()
