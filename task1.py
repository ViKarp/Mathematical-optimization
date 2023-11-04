import random as rd
import pandas as pd
from scipy import optimize
def f1(gamma):
    global a0, b0, g1, g2
    return 5 * (b0 + gamma * g2)**2 + 20 * (a0 + gamma * g1) * (b0 + gamma * g2) + 30 * (a0 + gamma * g1)**2 + 10

def f2(gamma):
    global a0, b0, g1, g2
    return abs(b0 + gamma * g2 +2) + abs(a0 + gamma * g1 + b0 + gamma * g2-2) + abs(2 * (a0 + gamma * g1) + b0 + gamma * g2 - 1) + abs(3 * (a0 + gamma * g1) + b0 + gamma * g2) + abs(4 * (a0 + gamma * g1) + b0 + gamma * g2 +1)

def dichotomy(func, a, b, eps, df):
    while (b - a) / 2 > eps:
        h = a + (b - a) / 2
        c = h - eps / 2
        d = h + eps / 2
        df.loc[ len(df.index )] = [a, b, (b - a) / 2, c, d, func(c), func(d)]
        if func(c) >= func(d):
            a = c
        else:
            b = d
    return h

def golden_ratio(func, a, b, eps, df):
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


a0 = rd.uniform(-10,10)
b0 = rd.uniform(-10,10)
g1 = rd.uniform(-1,1)
g2 = (1 - g1**2)**0.5
print("a0, b0 = ({}, {})".format(a0, b0))
print("g1, g2 = ({}, {})".format(g1, g2))

dfm1f1 = pd.DataFrame(columns=["a", "b", "eps", "c", "d", "f(c)", "f(d)"])
dfm2f1 = pd.DataFrame(columns=["a", "b", "eps", "c", "d", "f(c)", "f(d)"])
m1f1 = dichotomy(f1, -10, 10, 10**(-6), dfm1f1)
m2f1 = golden_ratio(f1, -10, 10, 10**(-6), dfm2f1)            
print("method: dichotomy; function 1: ", m1f1)
print("method: golden ratio; function 1: ", m2f1)
print("check: ", optimize.fminbound(f1, -10, 10))

dfm1f2 = pd.DataFrame(columns=["a", "b", "eps", "c", "d", "f(c)", "f(d)"])
dfm2f2 = pd.DataFrame(columns=["a", "b", "eps", "c", "d", "f(c)", "f(d)"])
m1f2 = dichotomy(f2, -10, 10, 10**(-6), dfm1f2)
m2f2 = golden_ratio(f2, -10, 10, 10**(-6), dfm2f2)
print("method: dichotomy; function 2: ", m1f2)
print("method: golden ratio; function 2: ",m2f2)
print("check: ", optimize.fminbound(f2, -10, 10))   

dfm1f1.to_csv("dfm1f1.csv")
dfm2f1.to_csv("dfm2f1.csv")
dfm1f2.to_csv("dfm1f2.csv")
dfm2f2.to_csv("dfm2f2.csv")