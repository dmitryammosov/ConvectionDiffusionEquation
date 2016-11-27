import matplotlib.pyplot as plt
import math

# problem:
# du/dt - du/dx = 0, 0 <= x <= L, 0 <= t <= T
# u(x, 0) = u0(x), 0 <= x <= L
# u(L, t) = nu2(t), 0 <= t <= T

N = 100
M = 100
A = 1.
L = 1.
T = 1.
h = L / N
tau = T / M


def u0(x):
    return A * math.exp(-x)


def nu2(t):
    return A * math.exp(-L - t)


def exact_u(x, t):
    return A * math.exp(-x - t)


def print_errors(u):
    for j in range(N + 1):
        print math.fabs(u[M][j] - exact_u(j * h, M * tau))


def plot_graph(u):
    t_index = N
    x = [0.0] * (N + 1)

    for i in range(N + 1):
        x[i] = i * h

    plt.plot(x, u[t_index])
    plt.xlabel('x')
    plt.ylabel('u')
    plt.title('t =' + str(M * tau))
    plt.show()


def implicit_scheme_second_order():
    a = [0.0] * (N + 1)
    b = [0.0] * (N + 1)
    c = [0.0] * (N + 1)
    f = [0.0] * (N + 1)
    alpha = [0.0] * (N + 1)
    beta = [0.0] * (N + 1)
    u = [[0] * (N + 1) for i in range(M + 1)]

    for j in range(N + 1):
        u[0][j] = u0(j * h)

    for i in range(M):
        for j in range(1, N):
            a[j] = 1.0 / (2 * h)
            b[j] = -1.0 / (2 * h)
            c[j] = -1.0 / tau
            f[j] = -u[i][j] / tau

        alpha[1] = 1.0 / (h / tau + 1.0)
        beta[1] = u[i][0] / (1.0 + tau / h)

        for j in range(1, N):
            alpha[j + 1] = b[j] / (c[j] - a[j] * alpha[j])
            beta[j + 1] = (a[j] * beta[j] + f[j]) / (c[j] - a[j] * alpha[j])

        u[i + 1][N] = nu2((i + 1) * tau)

        for j in reversed(range(N)):
            u[i + 1][j] = alpha[j + 1] * u[i + 1][j + 1] + beta[j + 1]

    return u


def main():
    u = implicit_scheme_second_order()
    plot_graph(u)
    print_errors(u)


main()