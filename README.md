# Introduction
Algoritmo para correção automática de linha de base em espectros, utilizando mínimos quadrados penalizados com ajuste ótimo do parâmetro de suavização.

# baseline-als-correction

Este repositório contém a implementação de um algoritmo para **correção automática de linha de base em espectros experimentais**, com base em mínimos quadrados penalizados assimetricamente ponderados (Asymmetric Least Squares - ALS) e seleção automática do parâmetro de suavização.

## 📌 Descrição

O algoritmo realiza a remoção da linha de base de sinais espectrais por meio de um ajuste iterativo que penaliza resíduos negativos. O parâmetro de suavização (`λ`) é otimizado automaticamente com base na minimização conjunta de critérios de ajuste (resíduos) e suavidade.

## ⚙️ Funcionalidades

- Correção robusta de linha de base em espectros com ruído ou deriva de fundo;
- Seleção automática do parâmetro de suavização;
- Visualização gráfica da função de custo para escolha de `λ`.

## 🧪 Exemplo de uso

```python
from baseline_correction import automatic_baseline

# y = espectro experimental (array 1D)
baseline = automatic_baseline(y)
corrected_spectrum = y - baseline
