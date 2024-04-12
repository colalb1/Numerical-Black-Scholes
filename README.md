# Numerical Black-Scholes-Merton
Numerical approximation of Black-Scholes-Merton PDE using Crank-Nicolson with Rannacher smoothing and TR-BDF2 for European and Barrier options.

## Precursor
This problem deals with spurious oscillations for payoff functions with discontinuities for numerical approximations of the Black-Scholes-Merton PDE (below) for European and barrier options. The point of this mini-app is to show various methods of handling these discontinuities via two second-order accurate finite differencing methods. This notebook will highlight how Rannacher time-stepping (smoothing) handles these oscillations for the Crank-Nicolson (CN) discretization by applying the $L$-stable Backward Euler for the first few timesteps (Crank-Nicolson is only $A$-stable, leading to oscillations given a discontinuity). I then extend the first-order accurate Backward Euler method to second-degree accuracy by implementing the Trapezoidal Rule with Second-Order Backward Difference Formula (TR-BDF2) to maintain $L$-stability while simultaneously increasing accuracy.

$$\frac{\partial V}{\partial t} + \frac{\sigma ^ 2 S ^ 2}{2} \frac{\partial ^ 2 V}{\partial S ^ 2} + (r - q)S \frac{\partial V}{\partial S} - rV = 0$$

I will define $A$-stability and $L$-stability as these terms are non-trivial and used often.

For $u_t = \lambda u$, a discretization is $A$-stable if $\lambda$ lies in the left half of the complex plane, $\mathbb{C}$. 

$A$-stable discretization $u^{(\ell)}$ where $\lim_{|z|\to\infty}{\frac{u^{(\ell + 1)}}{u^{(\ell)}}}\to 0$ for $z = h_t\lambda$ is $L$-stable.




## Crank-Nicolson with Rannacher Smoothing
Crank-Nicolson is second-order accurate but has spurious oscillations due to its $A$-stability. Thus, the Backward Euler method is applied to the first few (depends on how many you choose) timesteps to reduce oscillations before applying Crank-Nicolson. Applying the Backward Euler method for the first few timesteps does not affect the second-order accuracy of Crank-Nicolson. [This article](https://people.maths.ox.ac.uk/gilesm/files/NA-05-16.pdf) explains this phenomenon.

The last paragraph of page 107 in this [article](https://www.researchgate.net/publication/228524629_Convergence_analysis_of_Crank-Nicolson_and_Rannacher_time-marching) explains why I used quarter steps for ONE iteration of Rannacher smoothing.
</br>
</br>
The following system of equations solves the Crank-Nicolson discretization:
$$Pu_{\ell - 1} = Qu_\ell + \Delta tb_\ell$$ 

where $P = \Delta t * tridiag\left(-u, \frac{1}{\Delta t} + v, -w\right)$, $Q = \Delta t * tridiag\left(u, \frac{1}{\Delta t} - v, w\right)$

and $u = \frac{S}{4\Delta S}\left(\frac{\sigma ^ 2 S}{\Delta S} - (r - q)\right)$, $v = \frac{1}{2}\left(\frac{\sigma ^ 2 S ^ 2}{(\Delta S) ^ 2} + r\right)$, and $w = \frac{S}{4\Delta S}\left(\frac{\sigma ^ 2 S}{\Delta S} + (r - q)\right)$. $S$ is a vector of underlying prices, and the other terms are constant.

It appears as two parts in the implementation, but $b_{\ell}$ simplifies to the following:

$$b_\ell = \Delta t\begin{bmatrix} u[1] * (V[0, \ell] + V[0, \ell - 1])  \\\ 0 \\\ \vdots \\\ 0 \\\ w[-1] * (V[-1, \ell] + V[-1, \ell - 1]) \end{bmatrix}$$

where $V[:, \ell - 1]$ is the option price grid at the $\ell - 1$ timestep. I solved by iterating backward, so $V[:, \ell]$ is solved before $V[:, \ell - 1]$.
</br>
</br>
</br>
Following the same form ($Pu_{\ell - 1} = Qu_{\ell} + \Delta tb_\ell$), the system of equations for Backward Euler follows as such:

$P = \Delta t * tridiag\left(\alpha, \frac{1}{\Delta t} + \beta, \gamma\right)$ and $Q = I$

where $\alpha = \frac{S}{2\Delta S}\left(-\frac{\sigma ^ 2 S}{\Delta S} + (r - q)\right)$, $\beta = \frac{1}{\Delta t} + r + \frac{\sigma ^ 2 S ^ 2}{(\Delta S) ^ 2}$, and $\gamma = \frac{S}{2\Delta S}\left(-\frac{\sigma ^ 2 S}{\Delta S} - (r - q)\right)$.

Also, $$b_\ell = -\Delta t\begin{bmatrix} \alpha[1] * V[0, \ell] \\\ 0 \\\ \vdots \\\ 0 \\\ \gamma[-1] * V[-1, \ell] \end{bmatrix}$$



## Trapezoidal Rule with Second-Order Backward Differentiation (TR-BDF2) 
This method was originally implemented for Black-Scholes in [this](https://chasethedevil.github.io/lefloch_trbdf2_draft.pdf) paper. This was part of Fabien Le Floc’h's Ph.D. dissertation.

The trapezoidal steps are helper values for the BDF2 method and are ghosts on the final grid. I have satisfied all of the main method requirements, which include preprocessing half the grid with Crank-Nicolson, assembling the second-order backward method with the preprocessed grid, and implementing the boundary conditions for the TR and BDF2 steps. The output (without any manipulation) is very similar to the empirical surface. 

There is an issue at the underlying price cross sections at $0$ and $2K$ as they each display an oscillation; I imputed these edges with the second-to-last cross-section. Admittedly, this is a hacky fix. I tweaked the boundary conditions and checked my implementation against the formulas written in the paper but was not able to converge on the issue. My assumption is that it stems from the offset values.

I do not believe in "close enough," so I will instead say that the TR-BDF2 implementation is accurate up to (exclusive) the first and last space iterations. This imputation showed a negligible difference during error analysis, but there *are* some trivial workarounds to these oscillations.

One may simply increase the upper bound of the grid and add more iterations if they desire a more robust estimate (such as if the underlying is predicted to increase many times its original price). One could also make the grid finer to increase the method's accuracy. Hence, the imputation yields a practically infinitesimal difference (if the underlying likely will not increase more than twice its current price before $T$). I mention these ideas to show that an oscillation at **two** gridpoints is nearly insignificant (since $2 \ll M$).

To the readers who may have found the issue causing these oscillations, please let me know, and I will fix it.
<br>
<br>
I will defer a detailed explanation of TR-BDF2 to the article at the top of this section and will provide a basic explanation instead. TR-BDF2 uses a second-order accurate, fully implicit Runge Kutta timestepping technique for which ghost/intermediate time steps are used to quell accuracy loss using variable timesteps in the standard BDF2 method. Given $N$ time steps, the time grid is structured as follows:

$$\{0, \alpha, 1, 1 + \alpha, 2, \dots, N - 1, N - 1 + \alpha, N\}$$

The mentioned "intermediate grid" is the sequence of terms with $\alpha$: $\{\alpha, 1 + \alpha, \dots, N - 1 + \alpha\}$. The popular choices for $\alpha$ are $\frac{1}{2}$ and $2 - \sqrt{2}$. According to Le'Floch, $\alpha = \frac{1}{2}$ simplifies the BDF2 discretization, but $\alpha = 2 - \sqrt{2}$ decreases truncation error as this yields proportional Jacobians for the TR and BDF2 steps.

Crank-Nicolson is used to solve for the price at these timesteps before using them for BDF2. One may find the explanation of the Crank-Nicolson discretization in the CN section.

The discretization of BDF2 yields a system of equations of the form $Pu_{\ell - 1} = Qu_\ell + \frac{1}{\alpha(2 - \alpha)}u_{\ell - \alpha} + \alpha\Delta tb_\ell$ where

$P = \alpha\Delta t * tridiag\left(-u, \frac{1}{\alpha\Delta t} + v, -w\right)$ and $Q = -\frac{(1 - \alpha) ^ 2}{\alpha(2 - \alpha)}I$. 

$u, v, w$ and $b_\ell$ are defined the CN section.

Modifications to $u, v,$ and $w$ are made to agree with the boundary conditions:
* $u[0] = 0$
* $u[-1] = -\frac{(r - q)S[-1]}{2\Delta S}$
* $v[0] = \frac{1}{2}\left(\frac{(r - q)S[0]}{\Delta S} + r\right)$
* $v[-1] = \frac{1}{2}\left(r - \frac{(r - q)S[-1]}{\Delta S}\right)$
* $w[0] = \frac{(r - q)S[0]}{2\Delta S}$
* $w[-1] = 0$

This yields the solution on the standard $0, 1, \dots, N$ time grid.


## Conclusion
CN with Rannacher smoothing is more accurate than CN because it removes the oscillations from CN. This is because Backward Euler is $L$-stable instead of $A$-stable like CN.

TR-BDF2 is also $L$-stable, but shows larger errors for most implied volatility values. A method for which one changes the approximation method based on the implied volatility could be implemented to maintain the highest degree of accuracy.

The results seem to highlight the goal of my mini-app as the oscillations in the CN implementation became mild/negligible after Rannacher smoothing. The TR-BDF2 method seems to be without oscillations as it is more stable than the plain Crank-Nicolson scheme.

In creating this mini-app, I learned how to apply discretization schemes to second-order linear partial differential equations and how various dynamics of the Black-Scholes-Merton equation affect the payoff functions. I have also become much more familiar with novel finite differencing techniques for pricing said equation through literature analysis. I also learned the lesson of patience when implementing techniques described in academic papers.

Among the issues that could be improved, the oscillations at the edges of the space grid could be smoothed out in the TR-BDF2 method. TR-BDF2 is a rather uncommon method for pricing options so it was difficult to determine whether the oscillations were caused by an implementation issue or whether they were an artifact of the algorithm. At the time of this writing, there is no other open-source TR-BDF2 implementation that I am aware of. I would also like to extend these methods to different, more exotic, option types to test their accuracy in a variety of contexts.

Future directions for this project include but are not limited to: fixing TR-BDF2 to remove the barrier oscillations, implementing other exotic options to determine which scenarios each method works best in, adding varying risk-free rate and volatility dynamics, and implementing higher-order schemes to achieve greater accuracy.
