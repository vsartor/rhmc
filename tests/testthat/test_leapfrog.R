context("Leapfrog Implementations")

test_that(
    "leapfrog.grad 1-dimension",
    {
        g = function(x) -x * exp(-1/2 * x^2) / sqrt(2*pi)
        l = leapfrog.grad(g, -2, .8, 500, 0.01)
        expect_equal(l$q, -0.3920266161027324836)
        expect_equal(l$p, -0.095468083724177782434)
    }
)

test_that(
    "num.grad 1-dimension",
    {
        f = function(x) exp(-1/2 * x^2) / sqrt(2*pi)
        g = function(x) -x * exp(-1/2 * x^2) / sqrt(2*pi)
        expect_equal(num.grad(f, -2), g(-2))
        expect_equal(num.grad(f, -1), g(-1))
        expect_equal(num.grad(f, -.1), g(-.1))
        expect_equal(num.grad(f, 0), g(0))
        expect_equal(num.grad(f, .1), g(.1))
        expect_equal(num.grad(f, 1), g(1))
        expect_equal(num.grad(f, 2), g(2))
    }
)

test_that(
    "num.grad 2-dimension, linear case",
    {
        set.seed(123)
        y = rgamma(50, 5, 6)
        f = function(theta) theta[1] * sum(log(y)) - theta[2] * sum(y)
        g = function(theta) c(sum(log(y)), - sum(y))
        expect_equal(num.grad(f, c(4,5)), g(c(4,5)))
        expect_equal(num.grad(f, c(5,6)), g(c(5,6)))
        expect_equal(num.grad(f, c(6,7)), g(c(6,7)))
    }
)

test_that(
    "num.grad 3-dimension, non-linear case",
    {
        f = function(X) {
            x = X[1]; y = X[2]; z = X[3];
            x^2 * y * log(z) + exp(x * y *z)
        }
        g = function(X) {
            x = X[1]; y = X[2]; z = X[3];
            c(2 * x * y * log(z) + y * z * exp(x * y * z),
              x^2 * log(z) + x * z * exp(x * y * z),
              x^2 * y / z + x * y * exp(x * y * z))
        }
        expect_equal(num.grad(f, c(.5, -.4, 0.5)), g(c(.5, -.4, 0.5)))
    }
)

test_that(
    "leapfrog 1-dimension",
    {
        f = function(x) exp(-1/2 * x^2) / sqrt(2*pi)
        g = function(x) -x * exp(-1/2 * x^2) / sqrt(2*pi)
        lg = leapfrog.grad(g, -2, .8, 500, 0.01)
        ln = leapfrog(f, -2, .8, 500, 0.01)
        expect_equal(lg$q, ln$q)
        expect_equal(lg$p, ln$p)
    }
)
