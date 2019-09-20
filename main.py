import numpy as np
import math
import matplotlib.pyplot as plt
import random


def mu(t):
    return 1+t;


def calculate_u(l, T, phi_u, gamma, h_, tau_, u0):
    values_u = np.zeros((l+1,T+1))
    values_a = np.zeros((l+1,T+1))

    # for i in range(l):
    #     values_a[i,0]=0
    for j in range(T+1):
        values_u[0,j]=mu(j*tau_)
    for i in range(1,l+1):
        for j in range(0,T+1):
            b1=0
            if j != 0:
                b1=gamma[2 * j - 1] * (phi_u[j] + phi_u[j - 1]) / 2 *\
                   (u0[i, j] + u0[i, j - 1]) / 2 * math.exp(-gamma[2 * j - 1] * tau_) * tau_
                values_a[i, j] = values_a[i, j - 1] * math.exp(-gamma[2 * j - 1] * tau_) + b1
            values_u[i, j] = values_u[i - 1, j] + gamma[2*j]*(values_a[i,j]-phi_u[j]*u0[i,j]
                                                            + values_a[i-1,j]-phi_u[j]*u0[i-1,j])*h_/2

            # print(values_u[i,j], values_a[i, j], b1, flag)

    return values_u


def calculate_phi(l, T, phi0, gamma, h_, tau_, u0, g):
    values_phi = np.zeros(T+1)
    values_a = np.zeros((l+1, T+1))
    # values_phi[0]=1/(l*h_*gamma[0])*(math.log(mu(0)/g[0]))
    # print(mu(0),g[0])
    integral=0
    integral2 = 0
    for j in range(0,T+1):
        integral=0
        integral2 = 0
        for i in range(0,l+1):
            if j == 0:
                values_a[i,j]=0
            else:
                values_a[i, j] = values_a[i, j - 1] * math.exp(-gamma[2*j-1] * tau_) + \
                             gamma[2 * j - 1] * (phi0[j] + phi0[j - 1]) / 2 * (
                                     u0[i, j] + u0[i, j - 1]) / 2 * math.exp(-gamma[2*j-1] * tau_)*tau_
            # if i==0 or i==l:
            #     integral += gamma[2*j+1]*values_a[i,j]*math.exp(-gamma[j]*phi0[j]*(l*h-(i+0.5)*h))*h/2
            # else:
            #     integral += gamma[j]*values_a[i,j]*math.exp(-gamma[j]*phi0[j]*(l*h-i*h))*h
            # if i != 0:
            #     integral += gamma[2 * j] * (
            #                 values_a[i, j] + values_a[i - 1, j]) / 2 * math.exp(-gamma[j] * phi0[j] * (l * h - (i-0) * h)) * h
            if i != 0:
                integral+=(values_a[i, j] + values_a[i - 1, j]) / 2*h_
                integral2+=(u0[i, j] + u0[i-1, j]) / 2 * h_
        # print(g[j]-integral,g[j],integral)
        # values_phi[j]=1/(l*h*gamma[j])*(math.log(mu(j*tau))-math.log(g[j]-integral))
        values_phi[j] = (gamma[2*j]*integral+mu(j*tau)-g[j])/(integral2*gamma[2*j])
    return values_phi

def calculate_phi2(l, T, phi0, gamma, h_, tau_, u0, g):
    values_phi = np.zeros(T+1)
    values_a = np.zeros((l+1, T+1))
    # values_phi[0]=1/(l*h_*gamma[0])*(math.log(mu(0)/g[0]))
    # print(mu(0),g[0])
    integral=0
    integral2 = 0
    for j in range(0,T+1):
        integral=0
        integral2 = 0
        for i in range(0,l+1):
            if j == 0:
                values_a[i,j]=0
            else:
                values_a[i, j] = values_a[i, j - 1] * math.exp(-gamma[2*j-1] * tau_) + \
                             gamma[2 * j - 1] * (phi0[j] + phi0[j - 1]) / 2 * (
                                     u0[i, j] + u0[i, j - 1]) / 2 * math.exp(-gamma[2*j-1] * tau_)*tau_
            # if i==0 or i==l:
            #     integral += gamma[2*j+1]*values_a[i,j]*math.exp(-gamma[j]*phi0[j]*(l*h-(i+0.5)*h))*h/2
            # else:
            #     integral += gamma[j]*values_a[i,j]*math.exp(-gamma[j]*phi0[j]*(l*h-i*h))*h
            # if i != 0:
            #     integral += gamma[2 * j] * (
            #                 values_a[i, j] + values_a[i - 1, j]) / 2 * math.exp(-gamma[j] * phi0[j] * (l * h - (i-0) * h)) * h
            if i != 0:
                integral+=(values_a[i, j]*math.exp(gamma[2 * j]*phi0[j]*i*h) + values_a[i - 1, j]*math.exp(gamma[2 * j]*phi0[j]*(i-1)*h)) / 2*h_
                # integral2+=(u0[i, j] + u0[i-1, j]) / 2 * h_
        # print(g[j]-integral,g[j],integral)
        # values_phi[j]=1/(l*h*gamma[j])*(math.log(mu(j*tau))-math.log(g[j]-integral))
        # values_phi[j] = (gamma[2*j]*integral+mu(j*tau)-g[j])/(integral2*gamma[2*j])
        if j == 0:
            print(gamma[2 * j], math.exp(gamma[2 * j]*phi0[j]*l*h), integral, mu(j * tau), g[j], h*l * gamma[2 * j])
        values_phi[j] = math.log((gamma[2 * j] *integral + mu(j * tau)) / g[j]) / (h*l * gamma[2 * j])
    return values_phi


def calculate_gamma(l, T, phi, gamma0, h_, tau_, u0, g):
    values_gamma = np.zeros(2*T+1)
    values_a = np.zeros((l+1, T+1))

    for j in range(0,T+1):
        integral=0
        integral2 = 0
        for i in range(0,l+1):
            if j == 0:
                values_a[i,j]=0
            else:
                values_a[i, j] = values_a[i, j - 1] * math.exp(-gamma0[2*j-1] * tau_) + \
                             gamma0[2 * j - 1] * (phi[j] + phi[j - 1]) / 2 * (
                                     u0[i, j] + u0[i, j - 1]) / 2 * math.exp(-gamma0[2*j-1] * tau_)*tau_

            if i != 0:
                integral+=(values_a[i, j] + values_a[i - 1, j]) / 2*h_
                integral2+=(u0[i, j] + u0[i-1, j]) / 2 * h_
        values_gamma[2*j] = (-mu(j*tau_)+g[j])/(integral-phi[j]*integral2)
        if j!=0:
            values_gamma[2*j-1]=(values_gamma[2*j]+values_gamma[2*j-2])/2
    return values_gamma

def calculate_gamma1(l, T, phi, gamma0, h_, tau_, u0, g):


def get_u(l, T, phi_gu, gamma, h_, tau_, u0):
    u_get = calculate_u(l, T, phi_gu, gamma, h_, tau_, u0)
    delta_u = u_get - u0
    delta_u = delta_u ** 2
    delta = np.sum(delta_u)
    while delta > 0.0001:
        # print('here')
        u0 = u_get.copy()
        u_get = calculate_u(l, T, phi_gu, gamma, h_, tau_, u0)
        delta_u = u_get - u0
        delta_u = delta_u ** 2
        delta = np.sum(delta_u)
    return u_get


def get_phi(l, T, phi0, gamma, h_, tau_, u0, g):
    phi_0 = phi0.copy()
    u_g = get_u(l, T, phi_0, gamma, h_, tau_, u0)
    phi_v = calculate_phi2(l, T, phi_0, gamma, h_, tau_, u_g, g)
    phi_1 = phi_v.copy()
    delta_phi = phi_v - phi_0
    delta_phi = delta_phi ** 2
    delta = np.sum(delta_phi)
    # while delta > 0.001:
    #     phi_0 = phi_v.copy()
    #     u_g = get_u(l, T, phi_0, gamma, h_, tau_, u0)
    #     phi_v = calculate_phi2(l, T, phi_0, gamma, h_, tau_, u_g, g)
    #     delta_phi = phi_v - phi_0
    #     delta_phi = delta_phi ** 2
    #     delta = np.sum(delta_phi)
    # for i in range(4):
    i=0
    while delta > 0.001:
        phi_0 = phi_v.copy()
        u_g = get_u(l, T, phi_0, gamma, h_, tau_, u0)
        phi_v = calculate_phi2(l, T, phi_0, gamma, h_, tau_, u_g, g)
        delta_phi = phi_v - phi_0
        delta_phi = delta_phi ** 2
        delta = np.sum(delta_phi)
        i=i+1
        print(i)
    return (phi_v,phi_1, i)


def get_gamma(l, T, phi, gamma0, h_, tau_, u0, g):
    gamma_0 = gamma0.copy()
    u_g = get_u(l, T, phi, gamma_0, h_, tau_, u0)
    gamma_v = calculate_gamma(l, T, phi, gamma_0, h_, tau_, u_g, g)
    delta_phi = gamma_v - gamma_0
    delta_phi = delta_phi ** 2
    delta = np.sum(delta_phi)
    while delta > 0.01:
        gamma_0 = gamma_v.copy()
        u_g = get_u(l, T, phi, gamma_0, h_, tau_, u0)
        gamma_v = calculate_gamma(l, T, phi, gamma_0, h_, tau_, u_g, g)
        delta_phi = gamma_v - gamma_0
        delta_phi = delta_phi ** 2
        delta = np.sum(delta_phi)
    # for i in range(10):
    #     gamma_0 = gamma_v.copy()
    #     u_g = get_u(l, T, phi, gamma_0, h_, tau_, u0)
    #     gamma_v = calculate_gamma(l, T, phi, gamma_0, h_, tau_, u_g, g)
    #     delta_phi = gamma_v - gamma_0
    #     delta_phi = delta_phi ** 2
    #     delta = np.sum(delta_phi)
    return gamma_v



l = 100
T = 100
h = 0.01
tau = 0.01
phi_=np.zeros(T+1)
gamma_=np.zeros(2*T+1)
gamma0_=np.zeros(2*T+1)+1
u0_=np.zeros((l+1,T+1))+1
phi0_=np.zeros(T+1)+1
for j in range(T+1):
    # phi_[j]=2-math.exp(-tau*j)
    phi_[j] = math.cos(tau*j)
    # gamma_[2 * j] =2-math.exp(-tau*j)
    # gamma_[2*j]=j*tau+1   #+math.cos((j)/2*tau)
    # gamma_[2 * j] = j * tau+5
    gamma_[2*j]=5#2 - math.sin(tau * j)
    if j !=T:
        # gamma_[2 * j + 1] = (j + 0.5) * tau+5
        # gamma_[2 * j+1] = 2 - math.exp(-tau * (j+0.5))
        # gamma_[2*j+1]=(j+0.5)*tau+1  #*tau+math.cos((j+1)/2*tau)
        gamma_[2 * j+1] = 5#2 - math.sin(tau * (j+0.5))

u1 = get_u(l, T, phi_, gamma_, h, tau, u0_)
# u2 = get_u(l, T, phi0_, gamma_, h, tau, u0_)
# u1 = calculate_u(l, T, phi_, gamma_, h, tau, u0_, 0)
# u2 = calculate_u(l, T, phi0_, gamma_, h, tau, u0_, 1)
# print()
# u3=u1-u2
g_ = u1[l,:]

# for j in range(T):
#     r = random.randint(-1000000, 1000000)/1000000000*5
#     print(r)
#     g_[j]+=r

phi1 = get_phi(l,T,phi0_, gamma_, h, tau, u0_, g_)
# phi2=phi1-phi_
# gamma1=get_gamma(l,T,phi_, gamma0_, h, tau, u0_, g_)
# phi_v = calculate_phi(l,T,phi0_, gamma_,h, tau, u0_, u[l,:])
x=np.zeros(l+1)
for i in range(l+1):
    x[i]=i*h
t=np.zeros(T+1)
for j in range(T+1):
    t[j]= j*tau
t2=np.zeros(2*T+1)

fig, ax = plt.subplots()

ax.plot(t, phi_,
        linestyle = '-',
        linewidth = 1,
        color = 'red')

# Пунктирная линия ('--' или 'dashed'):
ax.plot(t, phi1[0],
        linestyle = '--',
        linewidth = 1,
        color = 'green')

# Точка-тире ('-.' или 'dashdot'):
ax.plot(t, phi1[1],
        linestyle = '-.',
        linewidth = 1,
        color = 'blue')
ax.plot(t, phi0_,
        linestyle = ':',
        linewidth = 2,
        color = 'gray')


plt.legend((r'Исходная функция $\varphi$',r'Вычисленная функция $\varphi$ на '+str(phi1[2])+r'-й итерации ($\varphi_{'+str(phi1[2])+r'}$)',r'Вычисленная функция $\varphi$ на 1-й итерации ($\varphi_1$)',r'$\varphi_0$'))

# graph = plt.plot(t,phi_)
# # graph = plt.plot(t, phi2)
# graph2 = plt.plot(t, phi1[0])
# graph3 = plt.plot(t, phi1[1])
# plt.ylabel('phi(t)')
# # plt.title('phi(t)=cos(t), gamma(t)=2-sin(t), g=g+eps, eps = [-0.01,0.01]')
# plt.title('phi(t)=cos(t), gamma(t)=2-sin(t)')

# print(gamma_)


# graph = plt.plot(t, gamma1[0::2])
# graph2 = plt.plot(t, gamma_[0::2])


# graph = plt.plot(x,u3[:,5])
# graph2 = plt.plot(x,u2[:,5])
# graph2 = plt.plot(x, u1)
plt.show()

