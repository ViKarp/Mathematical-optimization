import random as rd
import pandas as pd
import numpy as np
from scipy import optimize

def golden_ratio(func, a, b, eps, df):
    '''метод золотого сечения одномерной минимизации'''
    c = a + (3 - (5)**0.5) / 2 * (b - a)
    d = a + ((5)**0.5 - 1) / 2 * (b - a)
    while (b - a) / 2 > eps:
        df.loc[ len(df.index )] = [a, b, (b - a) / 2, c, d, func(c), func(d)]
        if func(c) > func(d):
            a, c = c, d
            d = a + ((5)**0.5 - 1) / 2 * (b - a)
        else:
            b, d = d, c
            c = a + (3 - (5)**0.5) / 2 * (b - a)
    return a+ (b - a) / 2

n=2
global_epsilon = 0.000000001

def e(j):
    '''выбор координаты в покоординатном спуске'''
    ar=np.zeros(n)
    ar[j-1]=1
    return ar

def f1(x):
    '''искомая функция'''
    a, b = x
    return 5 * b**2 + 20 * a * b + 30 * a**2 + 10

def f2(x):
    '''искомая функция'''
    a, b = x
    return abs(b + 2) + abs(a + b - 2) + abs(2 * (a) + b - 1) + abs(3 * (a) + b) + abs(4 * (a) + b +1)
def coord_boost(func, a, b, eps, df):
    '''покоординатный спуск'''
    x0 = np.array((rd.uniform(-10,10), rd.uniform(-10,10))) #выбор начальной точки
    xk = np.array((rd.uniform(-10,10), rd.uniform(-10,10)))
    df.loc[ len(df.index )] = [xk, np.linalg.norm(x0-xk), func(xk), abs(func(xk)-func(x0))] #запись в таблицу
    k=0
    while np.linalg.norm(x0-xk)>eps:
        for _ in range(2): #покоординатно
            pk=e(k-int(k/n)*n+1) #выбор координаты
            k+=1
            alpha=golden_ratio(lambda x: func(xk+x*pk),-30,30, eps, pd.DataFrame(columns=["a", "b", "eps", "c", "d", "f(c)", "f(d)"]))
            x0=xk
            xk=xk+alpha*pk
        df.loc[ len(df.index )] = [xk, np.linalg.norm(x0-xk), func(xk), abs(func(xk)-func(x0))]
    
    print('Final point: ', (xk, func(xk)))
    return(xk)


def derivative_x(epsilon, arg, f):
    '''производная по х'''
    return (f((global_epsilon + epsilon, arg)) -
            f((epsilon, arg))) / global_epsilon


def derivative_y(epsilon, arg, f):
    '''производная по у'''
    return (f((arg, epsilon + global_epsilon)) -
            f((arg, epsilon))) / global_epsilon

def count_alpha(mode, x, gradient, alpha, eps, func):
    '''подсчет альфа в зависимости от метода'''
    if mode == 1:
        return alpha
    if mode == 2:
        while func(x-alpha*gradient) > func(x) - 0.5 * alpha *  (np.linalg.norm(gradient))**2:
             alpha *= 0.5
        return alpha
    if mode == 3:
        alpha = golden_ratio(func=lambda l: func(x - l * gradient), a=-30, b=30, eps=eps, df=pd.DataFrame(columns=["a", "b", "eps", "c", "d", "f(c)", "f(d)"]))
        return alpha

def Gr_m(func, eps, mode, df, cnt=100000):
    alpha = 0.02    # Шаг сходимости, подсчитано вручную "< 2/L"
    X_prev = np.array([rd.uniform(-10,10), rd.uniform(-10,10)])   #Начальная точка
    gradient = np.array((derivative_x(X_prev[0], X_prev[1], func), derivative_y(X_prev[1], X_prev[0], func))) #Градиент
    X = X_prev - alpha * gradient
    k = 0
    df.loc[ len(df.index )] = [X_prev, func(X_prev), np.linalg.norm(gradient)]
    while np.linalg.norm(X - X_prev) > eps and k < cnt:
        k += 1
        X_prev = X.copy()   
        gradient = np.array((derivative_x(X_prev[0], X_prev[1], func), derivative_y(X_prev[1], X_prev[0], func)))
        df.loc[ len(df.index )] = [X_prev, func(X_prev), np.linalg.norm(gradient)]
        X = X_prev - count_alpha(mode, X_prev, gradient, alpha, eps, func) * gradient #подсчет следующей точки
    
    return X
    
 
 
dfcoordf1 = pd.DataFrame(columns=["xk", "||x0-xk||", "f(xk)", "|f(xk)-f(xk-1)|"])
m1f1 = coord_boost(f1, -10, 10, 10**(-6), dfcoordf1)
dfcoordf1.to_csv("dfcoordf1.csv")

dfgrad1f1 = pd.DataFrame(columns=["xk", "f(xk)", "||f'(xk)||"])
min_x, min_y = Gr_m(f1, 10**(-6), 1, dfgrad1f1)
minimum = (min_x, min_y, f1((min_x, min_y)))
print("const alpha: ", minimum)
dfgrad1f1.to_csv("dfgrad1f1.csv")

dfgrad2f1 = pd.DataFrame(columns=["xk", "f(xk)", "||f'(xk)||"])
min_x, min_y = Gr_m(f1, 10**(-6), 2, dfgrad2f1)
minimum = (min_x, min_y, f1((min_x, min_y)))
print("split alpha: ", minimum)
dfgrad2f1.to_csv("dfgrad2f1.csv")
dfgrad3f1 = pd.DataFrame(columns=["xk", "f(xk)", "||f'(xk)||"])
min_x, min_y = Gr_m(f1, 10**(-6), 3, dfgrad3f1)
minimum = (min_x, min_y, f1((min_x, min_y)))
print("fastest descent: ", minimum)
dfgrad3f1.to_csv("dfgrad3f1.csv")



dfcoordf2 = pd.DataFrame(columns=["xk", "||x0-xk||", "f(xk)", "|f(xk)-f(xk-1)|"])
m1f2 =coord_boost(f2, -10, 10, 10**(-6), dfcoordf2) 
dfcoordf2.to_csv("dfcoordf2.csv")