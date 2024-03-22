# Pacs Challenge 1 (Academic Year 2023-2024)
This README.md is written in Markdown so it is bettere to read it [online]()
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

   
