# Simulador Numérico de Injeção de CO₂

Este projeto implementa um modelo numérico unidimensional radial para simular a elevação da pressão em um aquífero salino durante a injeção de dióxido de carbono (CO₂).

Esse tipo de simulação é um componente crítico na avaliação da integridade de projetos de CCUS (Carbon Capture, Utilization, and Storage).

## Descrição do Modelo

O simulador resolve a equação da difusividade para o escoamento de fluido levemente compressível em meios porosos:
$$\frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{\partial P}{\partial r} \right) = \frac{\phi \mu c_t}{k} \frac{\partial P}{\partial t}$$

## Metodologia Numérica
Discretização espacial: Diferenças Finitas em coordenadas radiais
Discretização temporal: Método totalmente implícito (Backward Euler), garantindo estabilidade incondicional
Solução do sistema: Matriz esparsa tridiagonal resolvida via scipy.sparse.linalg

## Parâmetros do projeto:

| Variável | Descrição                    | Valor Padrão            |
| -------- | ---------------------------- | ----------------------- |
| `r_w`    | Raio do poço injetor         | 0.1 m                   |
| `r_e`    | Raio externo do reservatório | 2000 m                  |
| `k`      | Permeabilidade da rocha      | 100 mD ($10^{-13} m^2$) |
| `phi`    | Porosidade                   | 25%                     |
| `mu`     | Viscosidade do CO₂           | 0.06 cP                 |
| `q_inj`  | Taxa de injeção volumétrica  | 0.02 $m^3/s$            |
| `P_init` | Pressão original da formação | 15 MPa                  |


## Estrutura do Código
- Definição da malha:
Criação de blocos radiais (logarítmicos ou lineares) para capturar adequadamente o gradiente de pressão próximo ao poço

- Cálculo de transmissibilidade:
Define a facilidade de fluxo entre blocos adjacentes com base na geometria radial

- Matriz de acúmulo:
Representa a capacidade de armazenamento do meio poroso (porosidade + compressibilidade)

- Loop temporal:
Resolve o sistema linear a cada passo de tempo (Δt = 1 dia): $A \cdot P^{n+1} = B(P^n)$

## Requisitos

### Python
```
* Python 3.x
* NumPy
* SciPy
* Matplotlib
```
### C
```
Compilador GCC
Biblioteca matemática padrão (math.h)
```

## Como Executar

### Executando em Python
- Instale as dependências:
```
pip install numpy scipy matplotlib
```
- Execute o script principal:
```
python simulator.py
```

### Executando em C
- Compile o código:
```
gcc Simulator.c -o Simulator -lm
```
- Execute o programa:
```
./Simulator.c
```
