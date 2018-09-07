
#' Leapfrog Integrator for Hamiltonian Dynamics given the gradient
#'
#' Performs leapfrog steps for the Hamiltonian dynamics system given the
#' gradient function.
#'
#' @param g Gradient function for the system of interest.
#' @param q Initial position of the system.
#' @param p Initial momentum of the system.
#' @param L Number of leapfrog steps.
#' @param eps Size of each leapfrog step.
#'
#' @return A list containing both the position `q` and momentum `p` at the end
#'         of the integrator.
#'
#' @export
leapfrog.grad = function(g, q, p, L, eps) {
    # Start performing a half-step for the momentum
    p = p - .5 * eps * g(q)

    # Performs the first L-1 steps
    for (i in 1:(L-1)) {
        q = q + eps * p
        p = p - eps * g(q)
    }

    # The last step is a half-step for the momentum
    q = q + eps * p
    p = p - .5 * eps * g(q)

    # Return the final position and momentum
    list(q = q, p = p)
}


#' Numerical Approximation of the Gradient
#'
#' Performs simple numerical differentiation at a specific point.
#'
#' @param f Function to differentiate.
#' @param x Point to calculate differential.
#'
#' @return Approximation of the gradient.
#'
#' @export
num.grad = function(f, x) {
    d = length(x)
    g = numeric(d)
    for (i in 1:d) {
        # Calculate partial difference for i-th parameter
        h = sqrt(.Machine$double.eps) * x[i]
        xh = x[i] + h
        dx = xh - x[i]
        # Evaluate differential
        if (dx == 0) next() # g[i] is already set to 0
        Xh = x
        Xh[i] = xh
        g[i] = (f(Xh) - f(x)) / dx
    }
    g
}


#' Leapfrog Integrator for Hamiltonian Dynamics
#'
#' Performs leapfrog steps for the Hamiltonian dynamics system given the
#' integral of the gradient function for the system of interest, using
#' numerical derivatives for the gradient.
#'
#' @param f Function of the system of interest.
#' @param q Initial position of the system.
#' @param p Initial momentum of the system.
#' @param L Number of leapfrog steps.
#' @param eps Size of each leapfrog step.
#'
#' @return A list containing both the position `q` and momentum `p` at the end
#'         of the integrator.
#'
#' @export
leapfrog = function(f, q, p, L, eps) {
    # Start performing a half-step for the momentum
    p = p - .5 * eps * num.grad(f, q)

    # Performs the first L-1 steps
    for (i in 1:(L-1)) {
        q = q + eps * p
        p = p - eps * num.grad(f, q)
    }

    # The last step is a half-step for the momentum
    q = q + eps * p
    p = p - .5 * eps * num.grad(f, q)

    # Return the final position and momentum
    list(q = q, p = p)
}
