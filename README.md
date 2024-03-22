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







   
