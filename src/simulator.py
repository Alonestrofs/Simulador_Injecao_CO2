# Exemplo do simulador feito em python

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
import os
import json

def simular_injecao_co2():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(base_dir, 'data')
    reports_dir = os.path.join(base_dir, 'reports')
    
    # Cria a pasta reports caso não exista
    os.makedirs(reports_dir, exist_ok=True)
    
    params_file = os.path.join(data_dir, 'parametros.json')
    
    # Parâmetros Padrão (Caso a pasra data esteja vazia)
    params = {
        "r_w": 0.1, "r_e": 2000.0, "h": 50.0, "phi": 0.25,
        "k": 100e-15, "mu": 0.06e-3, "ct": 5e-10, "q_inj": 0.02, "P_init": 15e6
    }
    
    # Tenta carregar dados da pasta data
    if os.path.exists(params_file):
        with open(params_file, 'r') as f:
            params.update(json.load(f))
        print(f"Parâmetros carregados de: {params_file}")
    else:
        print("Arquivo de parâmetros não encontrado, usando valores padrão.")

    N = 150
    dr = (params["r_e"] - params["r_w"]) / N
    r = np.linspace(params["r_w"] + dr/2, params["r_e"] - dr/2, N)
    r_faces = np.linspace(params["r_w"], params["r_e"], N + 1)

    dias_simulacao = 365
    dt = 24 * 3600
    num_steps = int((dias_simulacao * 24 * 3600) / dt)

    T = np.zeros(N + 1)
    T[1:-1] = (2 * np.pi * params["k"] * params["h"] / params["mu"]) * r_faces[1:-1] / dr
    T[0] = 0 
    T[-1] = 0 

    V = np.pi * ((r + dr/2)**2 - (r - dr/2)**2) * params["h"]
    C = V * params["phi"] * params["ct"] / dt

    diag_principal = C + T[:-1] + T[1:]
    diag_inferior = -T[1:-1]
    diag_superior = -T[1:-1]
    
    A = diags([diag_inferior, diag_principal, diag_superior], [-1, 0, 1], format='csr')

    P = np.ones(N) * params["P_init"]
    resultados_P = {}
    dias_para_salvar = [1, 30, 90, 180, 365]

    for step in range(1, num_steps + 1):
        B = C * P
        B[0] += params["q_inj"]
        
        P = spsolve(A, B)
        
        if step in dias_para_salvar:
            resultados_P[step] = P.copy()

    plotar_resultados(r, resultados_P, params["P_init"], reports_dir)

def plotar_resultados(r, resultados_P, P_init, reports_dir):
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

  
    plot_path = os.path.join(reports_dir, 'grafico_ccus.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Gráfico exportado com sucesso para: {plot_path}")
    
    plt.show()

if __name__ == "__main__":
    simular_injecao_co2()
