import numpy as np
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt
import pandas as pd






def baseline_als(y, lam=1e8, ratio=0.0001):
    y = np.asarray(y)  # Converte a entrada para um array numpy
    N = len(y)  # Obtém o comprimento do array
    D = scipy.sparse.csc_matrix(np.diff(np.eye(N), 2))  # Calcula a segunda derivada discreta
    H = lam * D.dot(D.transpose())  # Calcula a matriz H, que é uma função da segunda derivada discreta
    w = np.ones(N)  # Inicializa um vetor de pesos com valores iguais a 1
    i = 0  # Inicializa o contador de iterações

    while True:
        W = scipy.sparse.spdiags(w, 0, N, N)  # Cria uma matriz diagonal esparsa com os pesos
        Z = W + H  # Calcula a matriz Z adicionando a matriz de pesos à matriz H
        z = scipy.sparse.linalg.spsolve(Z, w*y)  # Resolve o sistema linear Zz = wy para obter z
        d = y-z  # Calcula a diferença entre o espectro original e o espectro suavizado
        dn = d[d<0]  # Seleciona apenas os valores negativos da diferença
        m = np.mean(dn)  # Calcula a média dos valores negativos
        s = np.std(dn)  # Calcula o desvio padrão dos valores negativos
        wt = np.zeros(len(y))  # Inicializa um novo vetor de pesos
        for n in np.arange(len(y)):  # Loop sobre os elementos do vetor
            wt[n] = 1. / (1. + np.exp((2.*(d[n]-((2.*s)-m))/s)))  # Calcula os novos pesos com base na diferença
        if scipy.linalg.norm(w-wt)/scipy.linalg.norm(w) < ratio:  # Verifica a convergência dos pesos
            break  # Se convergir, sai do loop
        elif i==10000:  # Se o número máximo de iterações for atingido
            print(i)  # Imprime o número máximo de iterações
            break  # Sai do loop
        w = wt  # Atualiza os pesos
        i += 1  # Incrementa o contador de iterações

    return z, w  # Retorna o espectro suavizado e os pesos
    


def automatic_baseline(y, lambda_range=[1, 10], grid=21):
    """Algoritmo baseado na Seleção Automática do Parâmetro Ótimo para Correção da Linha de Base
    usando Mínimos Quadrados Penalizados assimetricamente ponderados. dx.doi.org/10.5573/ieie.2016.53.3.124"""
    
    y = np.asarray(y)  # Converte a entrada para um array numpy
    N = len(y)  # Obtém o comprimento do array
    D = np.transpose(np.diff(np.eye(N), 2))  # Calcula a segunda derivada discreta
    z, w = [], []  # Inicializa listas para armazenar os espectros suavizados e os pesos
    F, T = [], []  # Inicializa listas para armazenar os valores de fitness e suavização
    
    lambda_ = np.logspace(lambda_range[0], lambda_range[1], grid)  # Gera valores de lambda
    
    # Loop sobre os valores de lambda
    for lam in lambda_:
        z_, w_ = baseline_als(y, lam)  # Aplica a função baseline_als para obter o espectro suavizado e os pesos
        z.append(z_)  # Adiciona o espectro suavizado à lista
        w.append(w_)  # Adiciona os pesos à lista
    
    # Loop sobre os valores calculados
    for i in range(grid):
        F_, T_ = [], []  # Inicializa listas temporárias para armazenar os valores de fitness e suavização
        
        YZ = y - z[i]  # Calcula a diferença entre o espectro original e o espectro suavizado
        
        # Fitness Score: calcula o produto interno dos resíduos ponderados
        F_ = YZ.T.dot(np.diag(w[i]).dot(YZ))
        F.append(F_)  # Adiciona o valor de fitness à lista
        
        # Smooth Score: calcula a suavização penalizada
        T_ = z[i].T.dot(D.T).dot(D).dot(z[i])
        T.append(T_)  # Adiciona o valor de suavização à lista
    
    # Normaliza os valores de fitness e suavização
    F_norm = (F - np.min(F)) / (np.max(F) - np.min(F))
    T_norm = (T - np.min(T)) / (np.max(T) - np.min(T))
    
    # Calcula a soma normalizada de fitness e suavização
    S = F_norm + T_norm
    
    # Plota os resultados para análise
    plt.plot(np.log10(lambda_), S)
    
    # Encontra o índice do mínimo da soma normalizada
    opt = S.argmin()
    lambda_optim = lambda_[opt]  # Obtém o valor de lambda correspondente ao mínimo
    z_optim = z[opt]  # Obtém o espectro suavizado correspondente ao mínimo
    
    return z_optim  # Retorna o espectro suavizado ótimo





