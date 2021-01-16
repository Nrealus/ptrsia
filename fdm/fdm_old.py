import sys
import numpy as np
import matplotlib.pyplot as plt

def ode_forward_euler(f, u0, dt, T):
    N_t = int(T/dt)
    u = np.zeros(N_t+1, len(u0))
    t = np.arange(0, len(u0), N_t*dt)
    u[0,:] = u0
    for n in range(0,N_t):
        u[n+1,:] = u[n,:] + dt*f(u[n,:]n, t[n])
    
    return u, t

def solve_diffusion_pde_forward_euler(phi,u0,u1,n,a,t_step,x_step):
    
    def ode_rhs(uu, tt):

        gamma = lambda ttt : 

    T = 20 / n

    ts = np.arange(0,T,t_step)
    xs = np.arange(0,1,x_step)

    u = np.ndarray(shape=(len(ts),len(xs)),dtype=float)

    # initial conditions
    u[0,:] = phi(n, xs)

    # boundary conditions
    u[:,0] = u0(ts)
    u[:,-1] = u1(ts)


    return xs, u


def main(argv):
    n = 20
    a = 2

    def phi(n,x):
        return np.sin(np.pi*n*x)
    def u0(t):
        return 0
    def u1(t):
        return np.sin(t)/(1+t**2)

    (xs, u_expl) = solve_diffusion_pde_forward_euler(phi,u0,u1,n,a,1/500,1/500)
    plt.plot(xs,u_expl[-1,:])
    plt.savefig("u_expl_T.png")
    plt.close()

    #(xs, u_impl) = solve_diffusion_pde_forward_euler(phi,u0,u1,n,a,1/500,1/500)
    #plt.plot(xs,u_impl[-1,:])
    #plt.savefig("u_impl_T.png")

    with open("u_expl_T.txt", 'w+') as outfile_expl:
        for i in range(len(xs)):
            outfile_expl.write(str(xs[i]) +" "+ str(u_expl[:,i])+"\n")
        outfile_expl.close()

if __name__ == '__main__':
    main(sys.argv)
