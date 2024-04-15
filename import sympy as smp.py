import sympy as smp
from sympy import *
import numpy as np
# Определяем символы
x = smp.symbols('x')
w = smp.symbols('w', cls=smp.Function)
theta = smp.symbols('theta', cls=smp.Function)
q = smp.symbols('q')
E, I = smp.symbols('E I')
# Определяем дифференциальное уравнение
diff_eq = Eq(Derivative(w(x), x, x, x, x), q/(E*I))
solution = dsolve(diff_eq, w(x))
solution1 = Derivative(solution, x)
print("Решение дифференциального уравнения с граничными условиями:")
solution
#Параметры материала
E1 = 22000000000
E2 = 700000000
G12 = 500000000
nu21 = 0.23
ro = 1580
h0 = 0.00027
n = 5
alpha = [30, 30, -30, -30, 0]
alpha_rad = [smp.rad(a) for a in alpha]
sin_alpha = [smp.sin(a) for a in alpha_rad]
nu12 = nu21*E2/E1
print(nu12)
E12 = E2 * nu12 + 2 * G12
print(E12)
E_1 = E1/(1 + nu12 * nu21)
E_2 = E2/(1 + nu12 * nu21)
#коэффициенты жесткости
A11 = np.zeros(n)
A12 = np.zeros(n)
A22 = np.zeros(n)
A44 = np.zeros(n)

# Calculate the elements of the matrices
for i in range(n):
    A11[i] = E_1 * np.cos(alpha[i])**4 + E_2 * np.sin(alpha[i])**4 + 2*E12*(np.cos(alpha[i])**2) * np.sin(alpha[i])**2
    A12[i] = E_1 * nu12 + (E_1 + E_2 + 2 * E12) * (np.cos(alpha[i])**2) * np.sin(alpha[i])**2
    A22[i] = E_1 * np.sin(alpha[i])**4 + E_2 * np.cos(alpha[i])**4 + 2*E12*(np.cos(alpha[i])**2) * np.sin(alpha[i])**2
    A44[i] = (E_1 + E_2 + 2 * E1 * nu12) * (np.cos(alpha[i])**2) * np.sin(alpha[i])**2 + G12*((np.cos(alpha[i])**2) - np.sin(alpha[i])**2)

# Calculate the sums and multiply by h0
B11 = sum(A11) * h0 * 0.25
B12 = sum(A12) * h0 * 0.25 
B22 = sum(A22) * h0 * 0.25
B44 = sum(A44) * h0 * 0.25

print(B11, B12, B22)

E_x = B11 - (B12**2)/B22
print(E_x)

G_xy = B44
#Параметры сечения
b = 0.05
delta = n * h0
E_I = (2/3)*(b**3)*delta*E_x
G_F = (5/3)*b*delta*G_xy
#Записываем условие
L = 1.2
q = 1200
M = 450


boundary_conditions = [
Eq(w(0), 0), # w(0) = 0
Eq(Derivative(w(x), x).subs(x, 0), 0),
Eq(Derivative(w(x), x, x).subs(x, L), M/(E_I)),
Eq(Derivative(w(x), x, x, x).subs(x, L), 0)

]
diff_eq = Eq(Derivative(w(x), x, x, x, x), q/(E_I))
print("Граничные условия:")
for bc in boundary_conditions:
    print(bc)
solution = dsolve(diff_eq, w(x), ics={bc.lhs: bc.rhs for bc in boundary_conditions})


print("Решение дифференциального уравнения с граничными условиями:")
solution