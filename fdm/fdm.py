# Написать программу для численного решения дифференциального уравнения:
# 
# du/dt + a*d²u/dx² = 0, 0<=x<=1, 0<=t<=T, u=u(t,x);
# 
# с начальным условием:
# 
# u(0, x) = phi(x);
# 
# с граничными условиями:
# u(t, 0) = u0(t);
# u(t, 1) = u1(t).
# 
# Здесь
# phi(x) = sin(pi*n*x);
# u0(t)=0;
# u1(t)=sin(t) / (t^2+1);
# T = 20 / n;
# a = 2;
# где n>0 - номер вашего варианта.
# 
# Использовать сеточный метод для решения.
# Решить как при помощи явной схемы, так и при помощи неявной схемы.
# Найти профили решения для t=T и вывести в файл.
# (скоро здесь появится ссылка на гитхаб)

import sys
import numpy as np
import matplotlib.pyplot as plt

def compute_u_explicit_scheme(phi,u0,u1,n,a,t_step,x_step):
    T = 20 / n

    xs = np.arange(0,1,x_step)
    ts = np.arange(0,T,t_step)

    n_x = len(xs)
    n_t = len(ts)

    u = np.ndarray(shape=(n_x,n_t),dtype=float)

    # initial conditions
    u[:,0] = phi(n, xs)
    # boundary conditions
    u[0,:] = u0(ts)
    u[-1,:] = u1(ts)

    for i in range(1,n_x-1):
        for k in range(0,n_t-1):
            u[i,k+1] = u[i,k] + (a*t_step/(x_step**2))*(u[i+1,k] - 2*u[i,k] + u[i-1,k])

    return (xs, u)


def compute_u_implicit_scheme(phi,u0,u1,n,a,t_step,x_step):
    T = 20 / n

    xs = np.arange(0,1,x_step)
    ts = np.arange(0,T,t_step)

    n_x = len(xs)
    n_t = len(ts)

    u = np.ndarray(shape=(n_x,n_t),dtype=float)

    # initial conditions
    u[:,0] = phi(n, xs)
    # boundary conditions
    u[0,:] = u0(ts)
    u[-1,:] = u1(ts)

    # if we denote the main diagonal coefficients as b_i,
    # the upper diagonal coefficients as c_i,
    # and the lower diagonal coefficients as a_i,
    # then the Thomas algorithm is stable iff |b_i| > |c_i| + |a_i| (for every i)
    # ref: https://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)

    uld = -a/(x_step**2) # upper and lower diagonal coefficient (uld = a_i = c_i for every i)
    md = 1/t_step + 2*a/(x_step**2) # main diagonal coefficient (md = b_i for every i)
    # here we have |md| - 2*|uld| = 1/m > 0. as such, the algorithm will converge.

    uld_arr = np.zeros(n_x-1)
    md_arr = np.zeros(n_x-1)
    d = np.zeros(n_x-1)

    for k in range(n_t-1):
        d[1:-1] = u[1:-2,k] / t_step
        d[0] += uld * u[0,k+1]
        d[-1] += uld * u[-1,k+1]
        uld_arr[:] = uld
        md_arr[:] = md

        for i in range(1,n_x-2):
            m = uld_arr[i-1]/md_arr[i-1]
            md_arr[i] -= m * uld_arr[i-1]
            d[i] -= m * d[i-1]
        
        u[i-2,k+1] = d[-2] / md_arr[-2]
        for i in range(n_x-2, 0, -1):
            u[i,k+1] = (d[i] - uld_arr[i]*u[i+1,k+1]) / md_arr[i]

    return (xs, u)


def main(argv):
    n = 20
    a = 2
    def phi(n,x):
        return np.sin(np.pi*n*x)

    def u0(t):
        return 0

    def u1(t):
        return np.sin(t)/(1+t**2)

    def closure(mode, png_filename, txt_filename):
        if mode == "explicit":
            computation = compute_u_explicit_scheme
        elif mode == "implicit":
            computation = compute_u_implicit_scheme
        else:
            assert False, '"mode" arg must be either "explicit" or "implicit"'

        (xs, u_res) = computation(phi,u0,u1,n,a,1/500,1/500)
        
        plt.plot(xs,u_res[-1,:])
        plt.savefig(png_filename)
        plt.close()

        with open(txt_filename, 'w+') as outfile:
            for i in range(len(xs)):
                outfile.write(str(xs[i]) +" "+ str(u_res[i,:])+"\n")
            outfile.close()

    closure("explicit", "u_expl_T.png", "u_expl_T.txt")
    closure("implicit", "u_impl_T.png", "u_impl_T.txt")


if __name__ == '__main__':
    main(sys.argv)
