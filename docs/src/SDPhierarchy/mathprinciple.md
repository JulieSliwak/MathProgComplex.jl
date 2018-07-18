# The Lasserre hierarchy

Implementing Lasserre hierarchy for complex numbers, CJ, DM paper.

## Initial problem

## Deriving the moment relaxation

## Deriving the Sum Of Squares relaxation

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
