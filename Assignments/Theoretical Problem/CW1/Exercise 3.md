## Exercise 3

**Solution:**

##### a)

It is **NOT** possible to use the bisection method to find both roots of the function $$f(x)$$, and we can only find the root from the right side.

Because for the right side root $$r_{right}$$, we can certainly get $$f(a)f(b) \leq 0$$, and the function $$f$$ changes sign on any interval $$[a, b]$$ that contain the root $$r_{right}$$. Therefore, we can find the root $$r_{right}$$ by applying the bisection method.

However, for the left side root $$r_{left}$$, we cannot find it using the bisection method. Because there is no sign changes on the interval $$[a, b]$$ that contain the root $$r_{left}$$, and $$f(a)f(b) > 0$$, so it is impossible to find the root $$r_{left}$$.

(Personally, I think there is still a slim possibility that we can find the root $$r_{left}$$. If we can get a perfect interval $$(r_{left}-\delta, r_{left}+\delta)$$ for $$\delta \in \mathbb{R}$$ initially or after some iterations, then $$x_k$$ will be $$\frac{(r_{left}-\delta) + (r_{left}+\delta)}{2} = r_{left}$$, and $$f(x_k) = f(r_{left}) = 0$$, so that we can find the root $$r_{left}$$. But this probability is too small, so we only consider to find the root $$r_{right}$$ in the following.)





Assume that the error between $x_k$ and root $\alpha$ is $e_k$, and 
$$
|e_k| = |x_k - \alpha| \leq tol
$$
According to the priori criteria, we can know that the bisection method produces iterates satisfy
$$
|e_k| = |x_k - \alpha| \leq \frac{b - a}{2^{k+1}} = tol \\
k \ge log_2\Big(\frac{b - a}{tol}\Big) - 1
$$
where $a$ and $b$ are the bound of interval.

To find the root $$\alpha$$ by bisection with a relative accuracy $$tol = 10^{-10}$$, suppose we have a starting interval $$[-\pi, \pi]$$, then we can get
$$
k \ge log_2\Big(\frac{b - a}{tol}\Big) - 1 \\
k \ge log_2\Big(\frac{\pi - (-\pi)}{10^{-10}}\Big) - 1 \\
k \ge 34.87
$$
Therefore, after 35 iterations, the bisection method can found the root $\alpha$ to the accuracy $tol = 10^{-10}$.





##### b)

The formula for the Newton iteration is:

$$
x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)}, \ k = 0,1,2,...
$$
and the function $$f$$ is:
$$
f(x) = \frac{x}{2} - sinx + \frac{\pi}{6} - \frac{\sqrt{3}}{2}
$$
So, we can get $$f'(x) = \frac{1}{2} - cosx$$.

Therefore, the Newton's method for the problem is:
$$
x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)} = x_k - \frac{\frac{x}{2} - sinx + \frac{\pi}{6} - \frac{\sqrt{3}}{2}}{\frac{1}{2} - cosx}
$$




To determine the order of convergence for the Newton method, we should calculate the value $$p$$ from:
$$
||x_{n+1} - x|| \leq C ||x_n - x|| ^ p, \ p \geq 1
$$
Suppose that this function $f$ will finally converge to a root $x$, so $f(x) = 0$.

Since $x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)}$, we can get $x_{k+1} - x = x_k - x + \frac{f(x_k)}{f'(x_k)}$.

Next, according to the Taylor's theorem, we can get
$$
f(x) = f(x_k)+f'(x_k)(x-x_k)+\frac{f''(x_k)}{2}(x-x_k)^2+...+o(|x-x_k|^p)
$$

Because Taylor approximation is accurate enough such that we can ignore higher order terms, so we neglect third and higher powers of $(x-x_k)$, and now get
$$
f(x) \approx f(x_k)+f'(x_k)(x-x_k)+\frac{f''(x_k)}{2}(x-x_k)^2
$$
Since $f(x)=0$, then we can get
$$
0 &=& f(x_k)+f'(x_k)(x-x_k)+\frac{f''(x_k)}{2}(x-x_k)^2 \\
f(x_k)+f'(x_k)(x-x_k) &=& - \frac{f''(x_k)}{2}(x-x_k)^2 \\
\frac{f(x_k)}{f'(x_k)}+(x-x_k) &=& - \frac{f''(x_k)}{2f'(x_k)}(x-x_k)^2 \\
x - x_k + \frac{f(x_k)}{f'(x_k)} &=& - \frac{f''(x_k)}{2f'(x_k)}(x-x_k)^2 \\
x - \Big[x_k - \frac{f(x_k)}{f'(x_k)}\Big] &=& - \frac{f''(x_k)}{2f'(x_k)}(x-x_k)^2 \\
$$

Since $x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)}$, the equation will become
$$
x - x_{k+1} = - \frac{f''(x_k)}{2f'(x_k)}(x-x_k)^2 \\
||x_{n+1} - x|| \leq - \frac{f''(x_k)}{2f'(x_k)} ||x_n - x|| ^ 2
$$
Therefore, we can get $p=2$ from
$$
||x_{n+1} - x|| \leq C ||x_n - x|| ^ 2
$$

So, for both of the two zeros, the order of convergence of the Newton method is 2.





##### c)

According to the Contraction mapping theorem, we can know that if there exists a constant $0 < \Lambda < 1$ such that $|\phi'(x^0)|\leq \Lambda$ for all $x^0 \in [\frac{\pi}{2}, \pi]$, and $\phi(x) \in [\frac{\pi}{2}, \pi]$ for all $x \in [\frac{\pi}{2}, \pi]$, then we can prove the point iteration converges linearly to $\alpha$.

Since $\phi(x) = sinx + \frac{x}{2} - \Big(\frac{\pi}{6} - \frac{\sqrt{3}}{2} \Big)$, we can get
$$
\phi'(x^0) = cosx^0 + \frac{1}{2}
$$
According to the property of $cosx$, we can know that $\phi'(x^0)$ monotonic decreasing in interval $[\frac{\pi}{2}, \pi]$, and value range is $[\frac{1}{2}, -\frac{1}{2}]$. Thus, there exists a constant $0 < \Lambda < 1$ such that $|\phi'(x^0)|\leq \Lambda$ for all $x^0 \in [\frac{\pi}{2}, \pi]$.

Next, we should prove $\phi(x) \in [\frac{\pi}{2}, \pi]$ for all $x \in [\frac{\pi}{2}, \pi]$.

From $\phi'(x)$, we can calculate that when $x = \frac{2\pi}{3}$, $\phi'(x) = 0$. Because $\phi'(\frac{\pi}{2}) = \frac{1}{2}$, and $\phi'(\pi) = -\frac{1}{2}$, we can know that $\phi(x)$ monotonic increasing in interval $[\frac{\pi}{2}, \frac{2\pi}{3}]$, and monotonic decreasing in interval $[\frac{2\pi}{3}, \pi]$.

Thus, we can calculate that value range of $\phi(x)$.
$$
\phi\Big( \frac{2\pi}{3} \Big) = sinx + \frac{x}{2} - \Big(\frac{\pi}{6} - \frac{\sqrt{3}}{2} \Big) = \sqrt{3} + \frac{\pi}{6} \\ 
\phi\Big( \frac{\pi}{2} \Big) = sinx + \frac{x}{2} - \Big(\frac{\pi}{6} - \frac{\sqrt{3}}{2} \Big) = \frac{2 +\sqrt{3}}{2} + \frac{\pi}{12} \\
\phi(\pi) = sinx + \frac{x}{2} - \Big(\frac{\pi}{6} - \frac{\sqrt{3}}{2} \Big) = \frac{\pi}{3} + \frac{\sqrt{3}}{2}
$$
Therefore, the value range of $\phi(x)$ is $[\frac{\pi}{3} + \frac{\sqrt{3}}{2}, \sqrt{3} + \frac{\pi}{6}]$, and $\frac{\pi}{3} + \frac{\sqrt{3}}{2} > \frac{\pi}{2}$, and $\sqrt{3} + \frac{\pi}{6} < \pi$.

So, the fixed point iteration converges linearly to $\alpha$ for every initial guess $x^0 \in [\frac{\pi}{2}, \pi]$, with $|x^{k+1}-\alpha| \leq \Lambda |x^k - \alpha|$.
