import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from sympy import solve, Eq  # S, Rational, Mul, Pow,
from IPython.display import display  # Math, Latex
matplotlib.rcParams.update({'font.size': 12,
                            'font.family': 'calibri',
                            'mathtext.fontset': 'stix'})
# In[377]:
k, No, mu, B, U, T = sp.symbols('k N mu B U T'.split())
Nu, Nd = sp.var('N_u N_d')
E, invT, M, CB, U = sp.symbols('S (1/T) M C_B U'.split())
equations = [
    Eq(E/k, sp.log(sp.factorial(No)) -
        sp.log(sp.factorial(Nu)) -
        sp.log(sp.factorial(No-Nu))),   # [0]
    Eq(E/k, No*sp.log(No) -
        Nu*sp.log(Nu) -
        (No - Nu)*sp.log(No-Nu)),       # [1]
    Eq(invT, k/(2*mu*B)*sp.log((No - U/(mu*B))/(No + U/(mu*B)))),   # [2]
    Eq(M, No*mu*sp.tanh((mu*B)/(k*T))),                             # [3]
    Eq(M, -U/B),                # [4]
    Eq(CB, No*k*(((mu*B)/(k*T))**2/(sp.cosh((mu*B)/(k*T))**2))),    # [5]
    Eq(U, (mu*B*(No - 2*Nu))),  # [6]
    Eq(U/(mu*B), No - 2*Nu),    # [7]
    Eq(k, 8.62*10**(-5)),       # [8]
    Eq(mu, 5*10**(-8)),         # [9]
    Eq(T, 300),                 # [10]
    Eq(B, 1),                   # [11]
    Eq(No, 100)]                # [12]
equations


# In[378]:

newM = equations[4].rhs.subs(U, equations[6].rhs)  # [13]
equations.append(Eq(M, newM))
equations


# In[379]:

# from sympy.parsing.sympy_parser import parse_expr
newNu = (Eq(Nu, (solve(equations[13], Nu))[0]))
equations.append(newNu)  # [14]
equations


# In[380]:

newSk = Eq(E/k, (equations[1].rhs.subs(Nu, equations[14].rhs)))  # [15]
equations.append(newSk)
equations


# In[381]:

y = equations[15].rhs.subs([(No, equations[12].rhs),
                           (mu, equations[9].rhs),
                            (M, 1)])


# In[382]:

x = equations[7].rhs.subs(No, equations[12].rhs)


# In[384]:

print(type(x), type(y))
display(x, y)


# In[364]:

fig, ax = plt.subplots(figsize=(20, 8))
ax.set_ylabel("S/k", fontsize=18)
ax.set_xlabel('U/mu*B', fontsize=18)
# ax.spines["top"].set_visible(False)
# ax.spines["bottom"].set_visible(False)
# ax.spines["right"].set_visible((False)
# ax.spines["left"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# plt.xlabel('$\Delta t$ $(s)$',fontsize=20)
# plt.ylabel('$\Delta p$ $(hPa)$',fontsize=20)
plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
plt.title('Entropy as a Function of Energy', fontsize=22)
plt.plot((x, (Nu, 0, 100)), y)
plt.show()
# In[ ]:
# get_ipython().magic('pinfo sympy.evaluate')
