## Exercise 4

**Solution:**

##### a)

A point $x$ is a fixed point of the function $\phi$ should satisfy
$$
\phi(x) = x \ \ \ \text{if and only if} \ \ \ f(x) = 0
$$
 For $f(x) = x^3 - 2$, we can easily get when $x = \sqrt[3]{2}$, $f(x) = 0$.

Next, we should let $\phi(\sqrt[3]{2}) = \sqrt[3]{2}$.
$$
\phi(x) &=& x(1-\frac{\omega}{3}) + x^3(1-\omega) + \frac{2\omega}{3x^2} + 2(\omega-1)\\
\phi(\sqrt[3]{2}) &=& \sqrt[3]{2}(1-\frac{\omega}{3}) + 2(1-\omega) + \frac{2\omega}{3(\sqrt[3]{2})^2} + 2(\omega-1)\\
$$
Since $\phi(\sqrt[3]{2}) = \sqrt[3]{2}$, we can get
$$
\sqrt[3]{2}(1-\frac{\omega}{3}) + 2(1-\omega) + \frac{2\omega}{3(\sqrt[3]{2})^2} + 2(\omega-1) &=& \sqrt[3]{2} \\
\sqrt[3]{2}(1-\frac{\omega}{3}) + 2(1-\omega) + \frac{2\omega}{3(\sqrt[3]{2})^2} - 2(1-\omega) &=& \sqrt[3]{2} \\
\sqrt[3]{2}(1-\frac{\omega}{3}) + \frac{2\omega}{3(\sqrt[3]{2})^2} &=& \sqrt[3]{2} \\
2(1-\frac{\omega}{3}) + \frac{2\omega}{3} &=& 2 \\
2-\frac{2\omega}{3} + \frac{2\omega}{3} &=& 2 \\
2 &=& 2 \\
$$
Therefore, $\omega$ can be any value, and the root of $f(x) = 0$ is a fixed point.





##### b)

The method is locally convergent when $\phi : [a, b] \to \mathbb{R}$ is continuously differentiable and $\alpha \in (a, b)$ is a fixed point of $\phi$ such that $|\phi'(\alpha)| < 1$.

Therefore, for root $x$, we should find $\omega$ that makes $|\phi'(x)| < 1$.

Since $\phi(x) = x(1-\frac{\omega}{3}) + x^3(1-\omega) + \frac{2\omega}{3x^2} + 2(\omega-1)\\$, then we can get
$$
\phi'(x) = (1-\frac{\omega}{3}) + 3x^2(1-\omega) - \frac{4\omega}{3x^3}
$$
Because the root $x = \sqrt[3]{2}$, and $|\phi'(x)| < 1$, so
$$
-1 < \phi'(x) < 1 \\
-1 < (1-\frac{\omega}{3}) + 3x^2(1-\omega) - \frac{4\omega}{3x^3} < 1 \\
-1 < 1-\frac{\omega}{3} + 3(\sqrt[3]{2})^2(1-\omega) - \frac{4\omega}{6} < 1 \\
-1 < 1-\omega + 3(\sqrt[3]{2})^2(1-\omega) < 1 \\
- 1 < (1-\omega) [1 + 3(\sqrt[3]{2})^2] < 1 \\
- \frac{1}{1 + 3(\sqrt[3]{2})^2} < 1-\omega < \frac{1}{1 + 3(\sqrt[3]{2})^2} \\
- \frac{1}{1 + 3(\sqrt[3]{2})^2} - 1 < -\omega < \frac{1}{1 + 3(\sqrt[3]{2})^2} - 1 \\
\frac{2 + 3(\sqrt[3]{2})^2}{1 + 3(\sqrt[3]{2})^2} > \omega > \frac{3(\sqrt[3]{2})^2}{1 + 3(\sqrt[3]{2})^2} \\
$$
Therefore, for the value $\frac{3(\sqrt[3]{2})^2}{1 + 3(\sqrt[3]{2})^2} < \omega < \frac{2 + 3(\sqrt[3]{2})^2}{1 + 3(\sqrt[3]{2})^2}$, the method is locally convergent.





##### c)

The method of second order when satisfied $\phi : [a, b] \to \mathbb{R}$ is twice continuously differentiable and $\alpha \in (a, b)$ is a fixed point of $\phi$ satisfying $\phi'(\alpha) = 0$.

Therefore, for root $x$, we should find $\omega$ that makes $\phi'(x) = 0$.

Since $\phi'(x) = (1-\frac{\omega}{3}) + 3x^2(1-\omega) - \frac{4\omega}{3x^3}$, and the root $x = \sqrt[3]{2}$, we can get
$$
\phi'(x) = (1-\frac{\omega}{3}) + 3x^2(1-\omega) - \frac{4\omega}{3x^3} = 0 \\ 
\phi'(\sqrt[3]{2}) = (1-\frac{\omega}{3}) + 3(\sqrt[3]{2})^2(1-\omega) - \frac{4\omega}{3(\sqrt[3]{2})^3} = 0 \\ 
1-\frac{\omega}{3} + 3(\sqrt[3]{2})^2(1-\omega) - \frac{2\omega}{3} = 0 \\
1-\omega + 3(\sqrt[3]{2})^2(1-\omega) = 0 \\
\omega - 3(\sqrt[3]{2})^2(1-\omega) = 1 \\
\omega\sqrt[3]{2} - 6(1-\omega) = \sqrt[3]{2} \\
\omega\sqrt[3]{2} + 6\omega = \sqrt[3]{2} + 6 \\
\omega = 1 \\
$$
Therefore, for the value $\omega = 1$, the method of second order.





##### d)

Firstly, we want to find in which condition the method is of order higher than 2.

A Taylor expansion of $\phi(x_k)$ around $x = \alpha$ gives
$$
x_{k+1} - \alpha &=& \phi(x_k) - \phi(\alpha) \\
&=& \phi'(\alpha)(x_k - \alpha) + \frac{\phi''(\alpha)}{2}(x_k - \alpha)^2 + ... + \frac{\phi^{(p)}(\xi_k)}{p!}(x_k - \alpha)^p
$$
for some point $\xi_k$ between $x_k$ and $\alpha$. And we can know that if we want to let the method of order higher than 2, it should satisfying $\phi'(\alpha) = 0$, $\phi''(\alpha) = 0$, and $\phi^{(j)}$ is continuous.

Since $\phi'(x) = (1-\frac{\omega}{3}) + 3x^2(1-\omega) - \frac{4\omega}{3x^3}$, then we can get
$$
\phi''(x) = 6x(1-\omega) + \frac{4\omega}{x^4}
$$

We already get when $\omega = 1$, $\phi'(x) = 0$ before. However, when $\omega = 1 \\$, $\phi''(x) \ne 0$. 

Therefore, there is no value of $\omega$ make the method is of order higher than 2.
