# The Lasserre hierarchy

We refer the reader to the mentioned papers for exact description and results on the SDP hierarchy. We only recall the principle and the form of the problems to detail the implementation more precisely.

!!! warning

    To be completed

## A measure problem

The starting problem is a polynomial optimization problem in complex variables, that is a finite dimensional non linear problem:

```math
\inf_{z\in\mathbb{C}^n} \left\{ f(z) : g_i(z) \ge 0,  i=1, ..., m \right\}
```

This problem can be rewritten as an infinite dimensional linear optimization problem:

```math
\inf_{\mu \ge 0} \left\{ \int_K f(z)d\mu(z) \right\}
```

where $K = \left\{ z\in\mathbb{C}^n : g_i(x) \ge 0,  i=1, ..., \right\}$

At this point two elements make this problem intractable:

- the variable belongs to an infinite dimensional space,
- the variable itself is an infinite dimensional object.

A moment $\mu$ is exactly described by the sequence of its moments $\int z^\alpha d\mu(z), \forall \alpha \in \mathbb N^n$. Then the first point can be dealt with by selecting a finite dimensional subspace of ??, such as the subspace of measures having all moments of degree greater than $d$ null. This gives a creasing sequence of sub-spaces of ??:

```math
\left\{ \mu > 0 : \int z^\alpha d\mu(z) = 0 \forall \alpha \in \mathbb N^n, |\alpha| > d\right\}
```

Given an order $d$, one can therefore build the order $d$ moment relaxation, that is the following problem $\inf_{\mu \ge 0} \left\{ \int_K f(z)d\mu(z) \right\}$

Note:

- difference between real and complex.
- Equivalent, with moment description of the measure to: Moment Relaxation.
- symmetries
- multiorder

## The moment relaxation

## The Sum Of Squares relaxation

## Generic SDP problem form

The general form of primal SDP problem accepted by most solvers can be described as follows:

```math
\begin{aligned}
& \min_{Z_i\succeq 0, x_i\in \mathbb{R}}
& \sum_i A_{0i} \cdot Z_i + \sum_i b_{0i} x_i + c_0 \\
& \text{s.t.}
& \sum_i A_{ji} \cdot Z_i + \sum_i b_{ji} x_i + c_j =0\\
&&Z_i \succeq 0, \; i=1,\ldots, n_{sdp}&\\
&&x_i\in\mathbb{R}, \; i=1,\ldots,n_{scal}
\end{aligned}
```