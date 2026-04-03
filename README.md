Simulador Numérico de Injeção de CO₂ (CCUS)

Este projeto implementa um modelo numérico unidimensional radial para simular a elevação da pressão em um aquífero salino durante a injeção de dióxido de carbono (CO₂).

Esse tipo de simulação é um componente crítico na avaliação da integridade de projetos de CCUS (Carbon Capture, Utilization, and Storage).

## Descrição do Modelo

O simulador resolve a equação da difusividade para o escoamento de fluido levemente compressível em meios porosos:
$$\frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{\partial P}{\partial r} \right) = \frac{\phi \mu c_t}{k} \frac{\partial P}{\partial t}$$

## Metodologia Numérica
Discretização espacial: Diferenças Finitas em coordenadas radiais
Discretização temporal: Método totalmente implícito (Backward Euler), garantindo estabilidade incondicional
Solução do sistema: Matriz esparsa tridiagonal resolvida via scipy.sparse.linalg
