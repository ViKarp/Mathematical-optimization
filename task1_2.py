import matplotlib.pyplot as plt
import numpy as np
from sortedcontainers import SortedList
# Вариант 1

a = 5
b = 5
c = 4
d = 14
h = -1

def phi(x):
    return 1 / d * abs(a - x) * (1 / 2 * abs(x + b) + h) * (3 * abs(1 / 4 * x + c) + 1)

# График исходной функции
x = np.linspace(-20., 10., 300, endpoint=True)
y = phi(x)
plt.plot(x, y)
plt.tight_layout()
plt.show()

# Задание 1

# Модули меняют знак (ф-я меняет поведение) в точках x = a (x = 5), x = -b (x = -5), x = -4c (x = -16)

# На уч-ке [-20, -16] phi(x) = (5-x)/14 * (x+7)/2 * (3x+44)/4
# phi'(x) = 1/112 (-9x^2-100x+17) - парабола, ветви вниз,
# вершина (x0 = -100/18 ~ -5.5) вне рассматриваемого диапазона,
# max|phi'(x)| на левой границе при x = -20, L1 = 14.14

# На уч-ке [-16, -5] phi(x) = (x-5)/14 * (x+7)/2 * (3x+52)/4
# phi'(x) = 1/112 (9x^2+116x-1), ветви вверх, x0 = -116/18 ~ -6.44,
# максимум на левой границе при x = -16, L2 = 4

# На уч-ке [-5, 5] phi(x) = (5-x)/14 * (x+3)/2 * (3x+52)/4
# phi'(x) = 1/112 (-9x^2-92x+149), ветви вниз, x0 = -92/18 ~ -5.11
# максимум на правой границе при x = 5, L3 = 4.79

# На уч-ке [5, 10] phi(x) = (x-5)/14 * (x+3)/2 * (3x+52)/4
# phi'(x) = 1/112 (9x^2+92x-149),ветви вверх, x0 = -92/18 ~ -5.11
# максимум на правой границе при x = 10, L4 = 14.92

# L = max(14.92, 14.14, 4.79, 4) = 14.92 - константа Липшица

L = 14.92

# Задание 2
x1 = -20
x2 = 10

eps = 10 ** (-6)

# Начальная точка
start_x = 1 / (2 * L) * (phi(x1) - phi(x2) + L * (x1 + x2))
start_y = 1 / 2 * (phi(x1) + phi(x2) + L * (x1 - x2))

# Самосортирующийся список, первый элемент всегда будет тот что нужен
points = SortedList(key=lambda p: p[1])
points.add((start_x, start_y))

delta = 1 / (2 * L) * (phi(start_x) - start_y)
n = 1
while delta * 2 * L > eps:
    current_min = points[0]
    new_x_1 = current_min[0] - delta
    new_x_2 = current_min[0] + delta
    new_y = 1 / 2 * (phi(current_min[0]) + current_min[1])

    if (n - 1) % 5 == 0:
        print(n, *current_min, delta * 2 * L, new_x_1, new_x_2, new_y)

    points.pop(0)
    points.add((new_x_1, new_y))
    points.add((new_x_2, new_y))

    delta = 1 / (2 * L) * (phi(points[0][0]) - points[0][1])
    n += 1

print('min = ', points[0])