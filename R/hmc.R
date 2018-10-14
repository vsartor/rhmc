#' Numerical Gradient
#'
#' Performs numerical differentiation of a function at a specific point.
#' Uses some numerical tricks to always achieve a reliable, though not
#' necessarily optimal, error.
#'
#' @param f The function for which the gradient is desired.
#' @param x The point at which the gradient should be approximated.
#'
#' @return The gradient of the function `f` at `x`.
#'
#' @examples
#' func = function(x) exp(-0.5 * x ^ 2) / sqrt(2 * pi)
#' grad = function(x) -x * exp(-0.5 * x ^ 2) / sqrt(2 * pi)
#' num_grad(func, -2)
#' abs(num_grad(func, -2) - grad(-2))
#'
#' @export
num_grad = function(f, x) {
    # Faz uma aproximação numérica do gradiente da função f no ponto x
    d = length(x)
    g = numeric(d)

    for (i in 1:d) {
        # Calcula diferença parcial para o i-ésimo parâmetro da função

        # Truques numéricos
        h  = sqrt(.Machine$double.eps) * if (x[i] != 0) abs(x[i]) else 1e-8
        xh = x[i] + h
        dx = xh - x[i]

        # Cálculo da diferença parcial
        if (dx == 0) next # Evita divisão por 0 quando gradiente é 0
        Xh    = x
        Xh[i] = xh
        g[i]  = (f(Xh) - f(x)) / dx
    }

    g
}

#' Hamiltonian Dynamics
#'
#' Approximates Hamiltonian dynamics for some potential function and a L2-norm
#' kinectic funcion, assuming H(q,p) = U(q) + K(p).
#'
#' @param U Potential function of the system.
#' @param q Initial position vector.
#' @param p Initial momentum vector.
#' @param L Number of steps.
#' @param eps Size of each step.
#' @param m Mass vector.
#'
#' @return A list with the position `q` and momentum `p` at the end of the
#'         trajectory.
#'
#' @examples
#' U = function(x) exp(-0.5 * x^2) / sqrt(2 * pi)
#' hamiltonian_dynamics(U, -2, 0.8, 100, 0.1, 1)
#' hamiltonian_dynamics(U, -2, 0.85, 100, 0.1, 1)
#'
#' @export
hamiltonian_dynamics = function(U, q, p, L, eps, m) {
    p = p - .5 * eps * num_grad(U, q)
    for (i in 1:(L - 1)) {
        q = q + eps * p / m
        p = p - eps * num_grad(U, q)
    }
    q = q + eps * p / m
    p = p - .5 * eps * num_grad(U, q)

    list(q = q, p = p)
}

#' Hamiltonian Monte Carlo
#'
#' Performs Hamiltonian Monte Carlo for a desired target function.
#'
#' @param f Minus log-density function of interest.
#' @param init Initial point for the algorithm.
#' @param numit Number of iterations.
#' @param L Leapfrog parameter: number of steps.
#' @param eps Leapfrog parameter: size of each step.
#' @param mass Mass vector.
#'
#' @return A list with the chain with the samples of interest, the values of
#'         the log-density calculated at each step and the acceptance rate.
#'
#' @importFrom stats rnorm runif
#'
#' @examples
#' f = function(x) -dnorm(x, 20, 10, log = TRUE)
#' hmc(f, 19, 1000, 16, 0.3, 0.1)
#'
#' @export
hmc = function(f, init, numit, L, eps, mass) {
    d = length(init)

    # Cadeias
    q = matrix(nrow = d, ncol = numit)
    U = numeric(numit)

    # Inicializar valores
    q[ ,1] = init
    U[1]   = f(init)
    ar     = 0

    # Monte Carlo Hamiltoniano
    for (i in 2:numit) {
        p = rnorm(d, 0, sqrt(mass))
        K = sum(p^2 / (2 * mass))

        prop = hamiltonian_dynamics(f, q[ ,i-1], p, L, eps, mass)

        U_prop = f(prop$q)
        K_prop = sum(prop$p^2 / (2 * mass))

        if (runif(1) < exp(K - K_prop + U[i-1] - U_prop)) {
            q[ ,i] = prop$q
            U[i]   = U_prop
            ar     = ar + 1
        } else {
            q[ ,i] = q[,i - 1]
            U[i]   = U[i-1]
        }
    }

    list(chain = q, U = U, ar = ar / (numit - 1))
}
