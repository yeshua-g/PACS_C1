# Pacs Challenge 1 (Academic Year 2023-2024)
This README.md is written in Markdown so I suggest to read it [online](https://github.com/yeshua-g/PACS_C1/blob/main/README.md).

The aim of this code is to solve the following minimization problem:\
Let $f\colon R^n\to R$ be a given function that has a minimum. Find $x\in R^n$ such that\
$$x = \min_{y \in R^n} f(y).$$
To address the problem four methods have been implemented:
1. **Gradient descent method**
2. **Hevy ball method**
3. **Nesterov method**
4. **Adam method**

For all methods, the same stopping criteria have been adopted:
* *control on the step length:* $\lVert x_{k+1} - x_k \rVert < \varepsilon_s$
* *control on the residual:*  $|f(x_{k+1}) - f(x_k) | < \varepsilon_r$
* *limit number of iterations reached:* $k > k_{\max}$

To let the program work the User can just type `make` in the terminal and this will generate an executable called **minimization**, that can be executed by tiping `./minimization`.

The Method choice and parameters setting is done trough a GetPot file called **dataGetPot**, whose options are described in the following, together with specific parameters for each method.

Before starting it is **EXTREMELY IMPORTANT** to say that the code is written using the *muparser* library, thath is not included in this folder. So the User is required to change the first two lines of the makefile in order to set the proper path on his local machine

## Choosing the method
To choose the Method the User have to open the **dataGetPot** file with an editor and modify the [Method] section. 
```
[Method]
choice='your_preferred_choice'

// Change 'your_prefferred_choice with:
//'Gradient' or 'HeavyBall' or 'Nesterov' or 'Adam'
```
In order to define the function the User have to modify the [Function] section of the **dataGetPot** file. Note that the insertion of the Exact Gradient is not compulsory since for every method it has been given the possibility to use a numerical approximation of it by finite differences.
```
[Function]

#Function to minimize
Expression = 'x1*x2 + 4*x1^4 + x2^2 + 3*x1'

#Components of the gradient 
#.... if the problem is n-dimensional just add the
# needed components with the same notation up to
# ExactGrad_n = '...'
ExactGrad_0 = 'x2 + 16*x1^3 +3'
ExactGrad_1 = 'x1 + 2*x2'

# Dimension of the problem
Dimension = 2

#Initial condition
#.... if the problem is n-dimensional just add the
# needed components with the same notation up to
# x_0_n = '...'
x_0_0 = 0.
x_0_1 = 0.
```
Let's take a look at specific instructions for each method

### Gradient descent
The minimization problem is solved here by taking the following update rule $$x_{k+1} = x_k - \alpha_k \nabla f(x_k)$$
Here the User can choose among three strategies for computing $\alpha_k$ given $\alpha_0$
* *Exponential decay:* $\alpha_k=\alpha_0 e^{-\mu k}$
* *Inverse decay:* $\alpha_k=\frac{\alpha_0 }{1+\mu k}$
* Approximate line search with *Armijo rule* that checks the following condition $f(x_k)-f(x_k-\alpha_0 \nabla f(x_k)) \ge \sigma \alpha_0 || \nabla f(x_k) ||^2$ and if not satisfied sets $\alpha_0=\frac{\alpha_0}{2}$ and repeats. The value that finally satisfies the condition is the desired $\alpha_k$

To deal with Gradient descent parameters the User have to go in [GradientParameters] section of the **dataGetPot** file that looks like this:
```
[GradientParameters]
Alpha0 = 0.05
Mu = 0.2
EpsilonR = 1e-6
EpsilonS = 1e-6
MaxIterations = 1000
Sigma= 0.2
# Possible choices are 'ExponentialDecay',
# 'InverseDecay' and 'ApproximateLineSearch'
StepSizeMethod = 'ApproximateLineSearch'
# Set this to 0 if you want to use the Exact Gradient
# set it to 1 if you want to use the fd one
UseNumericalGradient = 1
```
### Heavy-ball 
The minimization problem is solved here by taking the following update rule:\
given $d_0 = -\alpha_0 \nabla f(x_0)$ $$x_{k+1} = x_k + d_k$$ $$d_{k+1} = \eta d_k - \alpha_{k+1} \nabla f(x_{k+1})$$
Here the User cannot longer choose *Armijo rule* since $d_k$ cannot be guaranteed to be a descent direction. So alongside *Exponential Decay* and *Inverse Decay* it has been added a possibility to use a *Fixed Alpha* ( $\alpha_0$ )

To deal with Heavy-ball parameters the User have to go in [HeavyballParameters] section of the **dataGetPot** file that looks like this:
```
[HeavyballParameters]
Alpha0 = 0.05
Eta = 0.9
Mu = 0.2
EpsilonR = 1e-6
EpsilonS = 1e-6
MaxIterations = 1000
# Possible choices are 'ExponentialDecay',
# 'InverseDecay' and 'FixedAlpha'
StepSizeMethod = 'FixedAlpha'
# Set this to 0 if you want to use the Exact Gradient
# set it to 1 if you want to use the fd one
UseNumericalGradient = 1
```
### Nesterov
The minimization problem is solved here by taking the following update rule:\
 $$y = x_k + \eta (x_k - x_{k-1})$$ $$x_{k+1} = y - \alpha_{k} \nabla f(y)$$
Stepsize strategies are the same as for the Heavy-ball method.
To deal with Nesterov parameters the User have to go in [Nesterov] section of the **dataGetPot** file that looks like this:
```
[Nesterov]
Alpha0 = 0.05
Eta = 0.9
Mu = 0.2
EpsilonR = 1e-6
EpsilonS = 1e-6
MaxIterations = 1000
# Possible choices are 'ExponentialDecay',
# 'InverseDecay' and 'FixedAlpha'
StepSizeMethod = 'FixedAlpha'
# Set this to 0 if you want to use the Exact Gradient
# set it to 1 if you want to use the fd one
UseNumericalGradient = 1
```
### Adam
This method works in the following way ( abuse of notation for vectorial operations)
$$m_k = \beta_1 m_{k-1} + (1-\beta_1) \nabla f(x_k)$$
$$v_k = \beta_2 v_{k-1} + (1-\beta_2) (\nabla f(x_k))^2$$
$$\hat{m_k} = \frac{m_k}{1-\beta_1^k}$$
$$\hat{v_k} = \frac{v_k}{1-\beta_2^k}$$
$$x_k = x_{k-1} -\frac{\eta}{\sqrt{\hat{v_k} + \epsilon}} \hat{m_k}$$
To deal with Adam parameters the User have to go in [Adam] section of the **dataGetPot** file that looks like this:
```
[Adam]
Beta_1 = 0.9
Beta_2 = 0.999
Eta = 0.01
Epsilon = 1e-8
EpsilonR = 1e-6
EpsilonS = 1e-6
MaxIterations = 1000
# Set this to 0 if you want to use the Exact Gradient
# set it to 1 if you want to use the fd one
UseNumericalGradient = 0
```



   
