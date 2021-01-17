import math

g = 9.84
T = 2
tau = 0.1
m1 = 1
m2 = 1
l1 = 1
l2 = 1

helper_const1 = (m1 / 2 + m2) * g * l1
helper_const2 = m2 * l1 * l2 / 4
helper_const3 = m2 * g * l2 / 2

helper_const4 = (m1/3 + m2) * l1**2
helper_const5 = (m2/3) * l2**2

def helper_get_coefficients(theta1, theta2, omega1, omega2):
    _temp = helper_const2 * math.sin(theta1 - theta2)
    f1 = helper_const1 * math.sin(theta1) - (omega2**2) * _temp
    f2 = helper_const3 * math.sin(theta2) + (omega1**2) * _temp 

    m11 = helper_const4
    m12 = helper_const2 * math.cos(theta1 - theta2)
    m21 = m12
    m22 = helper_const5

    detM = m11*m22 - m12*m21

    return (f1, f2, m11, m12, m21, m22, detM)

def epsilon1(theta1, theta2, omega1, omega2):
    (f1, f2, m11, m12, m21, m22, detM) = helper_get_coefficients(theta1, theta2, omega1, omega2)
    return (m22 * f1 - m12 * f2) / detM


def epsilon2(theta1, theta2, omega1, omega2):
    (f1, f2, m11, m12, m21, m22, detM) = helper_get_coefficients(theta1, theta2, omega1, omega2)
    return (-m21 * f1 + m11 * f2) / detM


def one_step_runge_kutta(theta1, theta2, omega1, omega2):
    """
    y[i+1] = y[i] + tau/6 * (k1 + 2*k2 + 2*k3 + k4)

    k1 = f(t[i], y[i])
    k2 = f(t[i] + tau/2, y[i] + tau/2 * k1)
    k3 = f(t[i] + tau/2, y[i] + tau/2 * k2)
    k4 = f(t[i] + tau  , y[i] + tau   * k3)
    """
    k1theta1 = omega1
    k1theta2 = omega2
    k1omega1 = epsilon1(theta1, theta2, omega1, omega2)
    k1omega2 = epsilon2(theta1, theta2, omega1, omega2)


    k2theta1 = omega1 + tau / 2 * k1omega1
    k2theta2 = omega2 + tau / 2 * k1omega2
    k2omega1 = epsilon1(
        theta1 + tau/2*k1theta1,
        theta2 + tau/2*k1theta2,
        omega1 + tau/2*k1omega1,
        omega2 + tau/2*k1omega2
        )
    k2omega2 = epsilon2(
        theta1 + tau/2*k1theta1,
        theta2 + tau/2*k1theta2,
        omega1 + tau/2*k1omega1,
        omega2 + tau/2*k1omega2
    )


    k3theta1 = omega1 + tau / 2 * k2omega1
    k3theta2 = omega2 + tau / 2 * k2omega2
    k3omega1 = epsilon1(
        theta1 + tau/2*k2theta1,
        theta2 + tau/2*k2theta2,
        omega1 + tau/2*k2omega1,
        omega2 + tau/2*k2omega2
        )
    k3omega2 = epsilon2(
        theta1 + tau/2*k2theta1,
        theta2 + tau/2*k2theta2,
        omega1 + tau/2*k2omega1,
        omega2 + tau/2*k2omega2
    )


    k4theta1 = omega1 + tau * k3omega1
    k4theta2 = omega2 + tau * k3omega2
    k4omega1 = epsilon1(
        theta1 + tau*k3theta1,
        theta2 + tau*k3theta2,
        omega1 + tau*k3omega1,
        omega2 + tau*k3omega2
        )
    k4omega2 = epsilon2(
        theta1 + tau*k3theta1,
        theta2 + tau*k3theta2,
        omega1 + tau*k3omega1,
        omega2 + tau*k3omega2
    )


    dtheta1 = (k1theta1 + 2 * k2theta1 + 2 * k3theta1 + k4theta1) * tau / 6
    dtheta2 = (k1theta2 + 2 * k2theta2 + 2 * k3theta2 + k4theta2) * tau / 6
    domega1 = (k1omega1 + 2 * k2omega1 + 2 * k3omega1 + k4omega1) * tau / 6
    domega2 = (k1omega2 + 2 * k2omega2 + 2 * k3omega2 + k4omega2) * tau / 6
    return (theta1 + dtheta1,
            theta2 + dtheta2,
            omega1 + domega1,
            omega2 + domega2
           )

t = 0
theta1 = 1
theta2 = 0.5
omega1 = 0
omega2 = 0
print(t, theta1, theta2, omega1, omega2)

with open("output.txt", 'w+') as outfile:
    outfile.write("t\ttheta1\ttheta2\tomega1\tomega2\n")
    for i in range(int(T / tau)):
        t = i * tau

        theta1, theta2, omega1, omega2 = one_step_runge_kutta(theta1, theta2, omega1, omega2)

        print(t, theta1, theta2, omega1, omega2)
        outfile.write("{}\t{}\t{}\t{}\t{}\n".format(t, theta1, theta2, omega1, omega2))
    outfile.close()
