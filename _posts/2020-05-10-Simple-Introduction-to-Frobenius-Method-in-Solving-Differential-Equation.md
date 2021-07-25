---
layout: post
title: Simple Introduction to Frobenius Method in Solving Differential Equation
author: Homam BETAR
categories: journal
tags: [Mathematics, ODEs, PDEs, Frobenius, Cauchy, Euler]
image: Taylor_series_1.png
---
## <ins> Introduction</ins>
We are not going to introduce this method with all the detailed proofs. You can find them in many references, such as those listed at the end. Instead our main objective in this blog is to gain an understanding of the mathematics required in physics, and also in numerical analysis. In this article, we will focus on explaining the main principle behind this method through a discussion of an example. Although providing all the proofs for the theories described here is far behind the scope of this blog, I always try to explain these theories and methods relying on intuition, which does not always work either because I fail to do that, or It won't be possible. _Is it sometimes impossible to explain intuitively an idea you understand?_

In many fields in physics and other domains, second-order ordinary differential equations (ODEs) with variable coefficients can result when studying some phenomena. Although there is no general method to solve this kind of equations, one can use many techniques hoping they will help in finding the solution. The difficulty spectre of solving them is very wide with equations at its far end without known analytical solutions. One way to solve the second-order ODE is to guess one of its two independent solutions, and then use it to lower the order of the equation to a first-order ODE which might be easier to treat (remember that the general solution of an ODE is a linear combination of its linearly independent solutions). Another way is the Frobenius method, named after mathematician [Ferdinand Georg Frobenius](http://mathshistory.st-andrews.ac.uk/Biographies/Frobenius.html), which is without any doubt much more general. When applying the Frobenius method, we can obtain the general formulae, as [power series](http://localhost:4000/journal/What-Is-Taylor-Expansion.html), for many special functions, such as Bessel, Legendre, or Hermite functions. However, this does not mean that applying this method will always lead us to a known function. Finally, we stress the fact that when we use this technique, we expand the solution of the ODE around its [regular singular points](http://localhost:4000/journal/What-are-the-Ordinary-and-Singular-Points-of-an-ODE.html).   

## <ins> Frobenius method</ins>
In this section, we focus on the following second-order ODE:

$$
\begin{equation}
x^2\frac{d^2y}{dx^2} + x a(x) \frac{dy}{dx}+b(x)y=0,
\tag{1}\label{eq:gen_form}
\end{equation}
$$

Where the variable coefficients $$a(x)$$ and $$b(x)$$ may be expanded as power series in $$x$$. Thus, 

$$
\begin{equation}
a(x) = \sum_{i=0}^{\infty}a_i x^i, \quad b(x) = \sum_{i=0}^{\infty}b_i x^i
\tag{2}\label{eq:coeffs}
\end{equation}
$$

By looking carefully at Eqs[\eqref{eq:gen_form},\eqref{eq:coeffs}], it is very logic to suggest that the solution may also be expanded as power series, then

$$
\begin{equation}
y(x) = \sum_{j=0}^{\infty}c_jx^{j+r},
\tag{3}\label{eq:sol_power_series}
\end{equation}
$$
Where $$c_0\neq 0$$.
Because we assume that all the terms in the ODE may be expanded, then after substituting the power series in Eq\eqref{eq:gen_form}, one can obtain all $$c_j$$ by requiring that the coefficient of each power of $$x$$ equal to zero if we want the equation to be correct for all values of $$x$$. $$r$$ in the previous relationship represents the lowest power of $$x$$ in the solution. The first and second derivatives of $$y$$ are then

$$
\begin{equation}
y' = \sum_{j=0}^{\infty}c_j(j+r)x^{j+r-1}, \quad y'' = \sum_{j=0}^{\infty}c_j(j+r)(j+r-1)x^{j+r-2},
\tag{4}\label{eq:derivatives}
\end{equation}
$$

Now substituting Eqs[\eqref{eq:derivatives}, \eqref{eq:sol_power_series}, \eqref{eq:coeffs}] in Eq\eqref{eq:gen_form}, we obtain the _indicial equation_ after equating to $$0$$ the cooefficient of $$x^{0}$$, which results when $$j=0$$ and $$i=0$$. This equation gives the values of the lowest power of $$r$$ in the solution. It reads, 

$$
\begin{equation}
r(r-1) + a_0 r + b_0 = 0, \quad c_0\neq 0,
\tag{5}\label{eq:indicial_eq}
\end{equation}
$$

Depending on the roots indicial equation $$r_1$$ and $$r_2$$ ($$r_1\leq r_2$$), we end with one of the following four possible scenarios (_Frobenius Algorithm_):
+ $$r_1\neq r_2$$, and the difference $$d=r_2-r_1$$ is not a positive integer. In this case, one obtains the two independent solutions $$y_1$$, and $$y_2$$ from Eq\eqref{eq:sol_power_series} after substituting $$r_1$$ and $$r_2$$ in it. The general solution is then, $$y = C_1 y_1 + C_2 y_2$$.
+ $$d$$ is a positive integer. In this case, $$y_1$$ and $$y_2$$ are given by,
$$
\begin{equation}
y_1 = \{(r-r_1)y(x,r)\}_{r=r_1}, \ y_2 = \{\frac{d}{dr}[(r-r_1)y(x,r)]\}_{r=r_1}, 
\tag{6}\label{eq:d_integer}
\end{equation}
$$
We can understand that as following. When $$d=r_2-r_1 \in \mathbb{Z^*}$$, then the terms $c_{j\geq d}$ in $$y(x,r_1)$$ have zeros in their denominators. Then, to eleminate this problem we multibly both sides of the ODE by $$(r-r_1)$$, which then leads us to Eq\eqref{eq:d_integer} -Notice that $$d\neq 0$$. _Even though it is no more than some algebra, I advise you to do the derivations because it empowers your understanding._
+ $$d=0$$. In this case, the two independent solutions are 
$$
\begin{equation}
y_1 = y(x,r_1), \ y_2 = \{\frac{d}{dr}y(x,r)\}_{r=r_1}, 
\tag{7}\label{eq:d_0}
\end{equation}
$$
+ In some cases we may have, not only $$c_0\neq 0$$, but also another coefficient $$c_m\neq0$$ in $$y(x,r_1)$$. Here the two independent solutions are obtained from $$y(x,r_1)$$. We can understand that since the general solution of the second-order ODE must contain two linearly independent solutions, each of them multiplied by arbitrary constants which are $$c_0$$ and $$c_m$$ here.

## <ins> Hold your Pen and start an Example </ins>
I decided to choose _Cauchy-Euler_ equation as an example to apply the Frobenius method because it represents a particular case. The second-order version of this equation reads,

$$
\begin{equation}
x^2y''+a x y' + b y = 0, \quad a, b \in \mathbb{C},
\tag{8}\label{eq:C_E}
\end{equation}
$$

In our example $$a=1$$, and $$b=-1$$. 

### Direct substitution of $$x^r$$:
One strategy to solve the previous equation is by substituting $$y=x^{r}$$, and then finding the roots of the resulting second-order $$r$$- polynomial. In our case,
$$
\begin{equation}
y=x^r \implies y'=rx^{r-1} \implies y''=r(r-1)x^{r-2},
\tag{9}\label{eq:C_E_der}
\end{equation}
$$
Then Eq\eqref{eq:C_E} gives, 

$$
\begin{equation}
r(r-1) + r - 1  = 0, \implies r_{1,2} = \pm 1,
\tag{10}\label{eq:Ind_C_E}
\end{equation}
$$

Thus, $$y_{1} = x$$, and $$y_2 = x^{-1}=\frac{1}{x}$$. Since these two functions are linearly independent, the general solution of Eq\eqref{eq:C_E} with $$a=1$$, and $$b=-1$$ is then

$$
\begin{equation}
y(x) = Ax + \frac{B}{x},
\tag{11}\label{eq:gen_sol_C_E_ge}
\end{equation}
$$

In Fig1, I plot Eq\eqref{eq:gen_sol_C_E_ge} for $$A=B=1$$. We observe two things, first when $$x\gg1$$, we can approximate the solution to $$y(x)\approx x$$, while in the opposite situation $$x\ll1$$, one writes $$y\approx 1/x$$. The following _matlab_ code is used to solve Eq\eqref{eq:C_E}.

```matlab
% Matlab code used to produce Fig.1
%		Homam BETAR 
%	<-- Institut Jean Lamour -->
%		10/05/2020
% Remark: Codes in this blog are not the most efficient ones.  
%         They only help to profound the understanding of the
%		       article's subject.  ^_^ 

syms y(x)
ode = x^2 * diff(y,x,2) + x * diff(y,x) - y == 0;
y_s = dsolve(ode);

disp(y_s); % it prints C2*x + C1/(2*x)
y_s = x + 1/x % I took C2 = C1 = 1 

fplot(x,yy, [-10, 10]); grid on;
```
### Frobenius method:
Here, we start by assuming that the solution to the equation can expanded as power series given by Eq\eqref{eq:sol_power_series}. Then substituting Eqs[\eqref{eq:sol_power_series}, \eqref{eq:derivatives}] in Eq\eqref{eq:C_E}, one fidns

$$
\begin{equation}
\begin{split}
& \sum_{j=0}^{\infty} c_j (j+r)(j+r-1)x^{j+r} + \sum_{j=0}^{\infty} c_j (j+r)x^{j+r} -\sum_{j=0}^{\infty}c_j x^{j+r} = 0, \\
& \sum_{j=0}^{\infty} c_j\Big\{\ (j+r)(j+r-1) + (j+r) -1 \Big\} x^{j+r} = 0,
\end{split}
\tag{12}\label{eq:Frob_Cauchy_1}
\end{equation}
$$

The previous equation must be correct for all values $$x$$. Thus, following the Frobenius method with $$c_0=0$$, one obtains the following two equations for the coefficients, 

$$
\begin{equation}
\begin{split}
& r(r-1)+ r - 1 = 0 \implies r_{1,2} = \pm 1, \\
& (j+r)(j+r-1) + (j+r) -1 = 0, \ (a_j\neq0, \ j\geq 1) ??
\end{split}
\tag{13}\label{eq:Frob_Cauchy_2}
\end{equation}
$$

It is obvious from the previous equation that the first relationship gives the same roots obtained in Eq\eqref{eq:Ind_C_E}. Therefore, when $$r=r_{1,2}$$ the second relationship in Eq\eqref{eq:Frob_Cauchy_2} is impossible. Then, the only possible solution to this situation is that $$c_j=0$$ for $$j\geq 1$$, which means either y$$y_1=c_0x$$ or $$y_2=c_0/x$$. Because these two functions are linearly independent, the general solution of the equation, taking into account that $$c_0$$ is arbitrary, reads

$$
\begin{equation}
y(x) = A'x + \frac{B'}{x},
\tag{14}\label{eq:Frob_Cauchy_3}
\end{equation}
$$

This soluation, as expected, is identical to Eq\eqref{eq:gen_sol_C_E_ge}.
### Yep, Another method: Lowering the order of the ODE:
If you can guess in this example that $$y_1(x)=A_1x$$ is a solution to Eq\eqref{eq:C_E}, then we can take advantage of knowing this solution, and lowering the order of the equation to becoming a first-order ODE. The second solution is written as 

$$
\begin{equation}
y_2 = y_1 u(x) = A_1xu(x),
\tag{15}\label{eq:lower_1}
\end{equation}
$$ 

The basic idea here is that $$y_1$$ will allow us to eliminate the last term (in this example $$-y$$) from the ODE, which converts it into a first-order one. After substituting Eq\eqref{eq:lower_1} into the ODE, it is easy to show that the resulting equation is

$$
\begin{equation}
\begin{split}
& xu'' + 3u' = 0, \quad take \ v = u' \implies \\
& xv' + 3v = 0, \implies v = \frac{du}{dx}= \frac{b}{x^3}, \implies \\
& u = \frac{-1}{2}\frac{b}{x^2} =  \frac{b_1}{x^2}, \implies  y_2 = A_1x * \frac{b_1}{x^2} = \frac{B_1}{x},
\end{split}
\tag{16}\label{eq:lower_2}
\end{equation}
$$ 

Then the general solution of the equation will be, 

$$
\begin{equation}
y(x) = A_1x + \frac{B_1}{x},
\tag{14}\label{eq:lower_3}
\end{equation}
$$

Which is again identical to Eqs[\eqref{eq:gen_sol_C_E_ge},\eqref{eq:Frob_Cauchy_3}].

## <ins> Conclusion </ins>
The Frobenius method is a reliable choice to find solutions of ODEs when relatively more straightforward ways fail to do so. On the one hand, applying this method helps us to find answers in closed form to many ODEs, such as Bessel, Legendre, Laguerre, Hermite equations. In fact, following this method enables us to write explicit formulae (in the form of series power) for these functions, which in turns makes solving the equation is a matter of recognizing its class. On the other hand, in most cases, it is not an easy task at all to obtain a closed-form of the power series. Understanding requires a lot of exercising and patience. So, _hold your pen and start exercising!_    
