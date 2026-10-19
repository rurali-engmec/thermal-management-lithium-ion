# thermal-management-lithium-ion
# Modelagem Computacional do Gerenciamento Térmico de Células de Íon-Lítio (LFP)

## 📌 Sobre o Projeto
Este projeto foi desenvolvido no âmbito da disciplina **SEM 0588 - Fenômenos de Transporte Computacional** no Departamento de Engenharia Mecânica da **EESC-USP (Escola de Engenharia de São Carlos)**. 

O objetivo principal é a simulação numérica bidimensional do comportamento térmico transiente em células de bateria de íon-lítio do tipo *pouch* submetidas a regimes severos de carga rápida (*Extreme Fast Charging*), avaliando a eficácia de diferentes sistemas de arrefecimento.

---

## 🔬 Metodologia e Fundamentação Matemática
* **Equação Governante:** Condução de calor bidimensional transiente com termo fonte volumétrico.
* **Modelo Eletroquímico:** Baseado no balanço de energia de Bernardi et al., contabilizando o aquecimento irreversível por Efeito Joule ($I^2R_{int}$).
* **Método Numérico:** Implementação do **Método de Diferenças Finitas (MDF)** utilizando o esquema implícito de **Direções Alternadas (ADI)** para garantir estabilidade incondicional.
* **Resolução Matricial:** Utilização otimizada do **Algoritmo de Thomas (TDMA)** para sistemas tridiagonais, reduzindo o custo computacional de $O(N^3)$ para $O(N)$.
* **Condições de Contorno:** Aplicação de condições de terceira espécie (Robin / convecção) com o rigoroso **Método do Nó Fictício** para preservar a precisão de segunda ordem nas fronteiras.
* **Refinamento Avançado:** Modelagem não-homogênea da geração de calor utilizando uma **distribuição espacial Gaussiana** para capturar o fenômeno de concentração de corrente (*current crowding*) nos terminais elétricos (*tabs*).

---

## 📊 Principais Resultados
1. **Independência de Malha:** Validação da convergência numérica para malhas superiores a $60 \times 60$ nós, com erro relativo inferior a 0,05%.
2. **Análise Paramétrica (Convecção):** 
   * *Modelo Homogêneo:* Exigia um coeficiente convectivo rigoroso de $h = 166 \, \text{W/m}^2\cdot\text{K}$ para manter a célula abaixo do limite crítico de $60 \, ^\circ\text{C}$.
   * *Modelo Não-Homogêneo (Gaussiano):* Demonstrou que, ao focar o calor nos terminais próximos às fronteiras de resfriamento, um coeficiente de $h = 64 \, \text{W/m}^2\cdot\text{K}$ já é suficiente para garantir a segurança operacional, justificando estratégias modernas de *tab cooling*.

---

## 🛠️ Tecnologias e Ferramentas Utilizadas
* **Linguagem / Ambiente:** GNU Octave / MATLAB (implementação dos solvers numéricos e pós-processamento gráfico).
* **Engenharia Computacional:** Elementos finitos/diferenças finitas, métodos numéricos implícitos e análise de transferência de calor.

---

## 📂 Estrutura do Repositório
- `/codigo/` — Scripts em Octave/MATLAB com os algoritmos ADI e TDMA.
- `/relatorio/` — Relatório técnico completo em formato PDF.
- `/resultados/` — Mapas de calor (isotermas), curvas transientes e gráficos de convergência.
