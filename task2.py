from scipy import optimize
from scipy.optimize import Bounds
import pandas as pd
from random import randint
#Метод ломаных
#Константа Липшица на всем отрезке
L = 12+3/34
def F(x):
    '''Целевая функция'''
    return (abs(-4-x)*(abs(x-4)/2-1)*(3*abs(x/4+4)+1))/17

def broken_line(eps, L, start_x, start_y):
    points = []
    points.append((start_x, start_y))
    dflom = pd.DataFrame(columns=["n","current min x", "current min y", "2L*delta", "new x1", "new x2", "new y"])
    delta = 1 / (2 * L) * (F(start_x) - start_y)
    n = 1
    while delta * 2 * L > eps:
        cur_min = points[0]
        new_x_1 = cur_min[0] - delta
        new_x_2 = cur_min[0] + delta
        new_y = 1 / 2 * (F(cur_min[0]) + cur_min[1])

        if (n - 1) % 5 == 0:
            dflom.loc[ len(dflom.index )] = [n, *cur_min, round(delta * 2 * L,4), new_x_1, new_x_2, new_y]

        points.pop(0)
        points.append((new_x_1, new_y))
        points.append((new_x_2, new_y))
        points.sort(key= lambda x: x[1])

        delta = 1 / (2 * L) * (F(points[0][0]) - points[0][1])
        n += 1
    dflom.loc[ len(dflom.index )] = [n, *points[0], 0, '-', '-', '-']
    dflom.to_csv('dflom.csv',index=False)
x1 = -20
x2 = 10

eps = 10 ** (-6)


# Начальная точка
start_x = 1 / (2 * L) * (F(x1) - F(x2) + L * (x1 + x2))
start_y = 1 / 2 * (F(x1) + F(x2) + L * (x1 - x2))

broken_line(eps,L,start_x, start_y)


bounds = Bounds ([x1], [x2])
print('real min = ', optimize.minimize(F,start_x, bounds=bounds).x)

#Метод Ньютона-Рафсона

def F1(x):
    return [((i+1)**2+3)*(i+2)**2 for i in x]
def F2(x):
    return 4*x**3+18*x**2+32*x+24
def F3(x):
    return 12*x**2+36*x+32
def newton(x0, x1, x2, eps):
    x_last = x0
    n = 1
    alpha = 0.5
    x_next = x_last - alpha*F2(x_last)/F3(x_last)
    if x_next > x2:
        x_next = x2
    elif x_next < x1:
        x_next = x1
    df_tangent = pd.DataFrame(columns=["n","x", "first derivative", "second derivative", "alpha"])
    df_tangent.loc[ len(df_tangent.index )] = [n, x_last, F2(x_last), F3(x_last), alpha]
    while abs(x_next - x_last) > eps:
        n+=1
        x_last = x_next
        alp_pret = [x_last - 0.1*alp*F2(x_last)/F3(x_last) for alp in range(1,11)]
        x_next = alp_pret[F1(alp_pret).index(min(F1(alp_pret)))]
        if x_next > x2:
            x_next = x2
        elif x_next < x1:
            x_next = x1
        alpha = 0.1*(alp_pret.index(min(alp_pret))+1)
        df_tangent.loc[ len(df_tangent.index )] = [n, x_last, round(F2(x_last),3), round(F3(x_last),3), alpha]
    df_tangent.loc[ len(df_tangent.index )] = [n+1, x_next, '-', '-', '-']
    df_tangent.to_csv('df_tangent.csv',index=False)
x1 = -10
x2 = 10
x0 = randint(x1,x2)
newton(x0,x1,x2,eps)
bounds = Bounds ([x1], [x2])
print('real min = ', optimize.minimize(F1,start_x, bounds=bounds).x)
