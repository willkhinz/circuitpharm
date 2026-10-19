# Theorem — Positive Fast Deactivation Amplitude in Two-Open Schemes

**Status:** Mathematical Theorem (Closed-form weight identity, spectral sign partition theorem, sufficient condition on rates, and necessary condition for multi-component positive deactivation).  
**Branch:** `claude/agent-brief-20`  
**Date:** 2026-10-18  

---

## 1. Statement of the Theorems

### Theorem 1 (Spectral Sign Partition for Two-Open Conducting Manifolds)

Let $\mathcal{M}$ be a continuous-time Markov receptor model of GABA-A receptor kinetics possessing the 6-state sequential two-open topology (`kinetic_two_open`):
$$R \underset{k_{\text{off}}}{\overset{2 k_{\text{on}} x}{\rightleftharpoons}} AR \underset{2 k_{\text{off}}}{\overset{k_{\text{on}} x}{\rightleftharpoons}} A_2R \underset{\alpha}{\overset{\beta}{\rightleftharpoons}} A_2O \underset{\gamma}{\overset{\delta}{\rightleftharpoons}} A_2O'$$
with desensitisation branching from $A_2R$:
$$A_2R \underset{r}{\overset{d}{\rightleftharpoons}} A_2D$$
operating under detailed balance, where all microscopic rate constants $(k_{\text{on}}, k_{\text{off}}, \alpha, \beta, \gamma, \delta, d, r)$ are strictly positive. Both $A_2O$ and $A_2O'$ are conducting states with identical unit conductance, so the observable open probability is:
$$P_O(t) = p_{A_2O}(t) + p_{A_2O'}(t)$$

Consider the deactivation protocol following an agonist application pulse of arbitrary duration $\Delta t > 0$ and concentration profile $g(t) \ge 0$, followed by instantaneous washout to zero agonist at $t = 0$.

Let $B$ be the $4 \times 4$ subgenerator governing the autonomous doubly-liganded subnetwork $\mathcal{D} = \{A_2D, A_2R, A_2O, A_2O'\}$ at zero agonist:
$$B = \begin{pmatrix} -r & r & 0 & 0 \\ d & -(2 k_{\text{off}} + \beta + d) & \beta & 0 \\ 0 & \alpha & -(\alpha + \delta) & \delta \\ 0 & 0 & \gamma & -\gamma \end{pmatrix}$$
with real symmetric Jacobi symmetrization $H = D B D^{-1}$, where $D = \text{diag}(\sqrt{d/r}, 1, \sqrt{\beta/\alpha}, \sqrt{(\beta/\alpha)(\delta/\gamma)})$.

Let $s_k = -\lambda_k > 0$ be the decay rate of relaxation mode $k \in \{0, 1, 2, 3\}$.

Then:
1. **Exact Gateway Weight Identity:** The effective projection weight of the gateway shut state $A_2R$ onto the total conducting observable for mode $k$ is given in closed form by:
$$\mathcal{W}^{(k)}_{A_2R} = - (u^{(k)}_{A_2R})^2 \beta \frac{s_k - (\gamma + \delta)}{(s_k - \mu_1)(s_k - \mu_2)}$$
where $\mu_1 > \mu_2 > 0$ are the eigenvalues of the open-manifold submatrix:
$$\begin{pmatrix} \alpha + \delta & -\sqrt{\delta\gamma} \\ -\sqrt{\delta\gamma} & \gamma \end{pmatrix}, \quad \mu_{1,2} = \frac{(\alpha + \delta + \gamma) \pm \sqrt{(\alpha + \delta + \gamma)^2 - 4\alpha\gamma}}{2}$$
and satisfy the strict interlacing order:
$$0 < \mu_2 < \gamma + \delta < \mu_1$$

2. **Spectral Sign Partition:** The gateway weight carries sign:
$$\operatorname{sgn}\left(\mathcal{W}^{(k)}_{A_2R}\right) = \begin{cases} +1 & \text{if } s_k \in (\gamma + \delta, \mu_1) \cup (0, \mu_2) \\ -1 & \text{if } s_k \in (\mu_1, \infty) \cup (\mu_2, \gamma + \delta) \end{cases}$$

3. **Reachable Initial Occupancy Invariance:** Because probability enters $\mathcal{D}$ exclusively through the gateway state $A_2R$, Duhamel's variation-of-parameters formula yields:
$$\text{amp}_k = \mathcal{W}^{(k)}_{A_2R} \int_0^{\Delta t} \sigma(s) e^{-s_k (\Delta t - s)} ds$$
Consequently, for every reachable initial state $p_{\mathcal{D}}(\Delta t)$:
$$\operatorname{sgn}(\text{amp}_k) \equiv \operatorname{sgn}\left(\mathcal{W}^{(k)}_{A_2R}\right)$$
The sign of each mode amplitude is an intrinsic property of the spectrum and does not depend on pulse duration or concentration.

---

### Theorem 2 (Sufficient Condition for Positive Fast Amplitude)

A two-open scheme admits a positive fast relaxation mode $s^* \in (\gamma + \delta, \mu_1)$ with $\mathcal{W}_{A_2R} > 0$ whenever the microscopic rates satisfy the explicit algebraic inequality:
$$\gamma + \delta > r \quad \text{and} \quad \frac{2 k_{\text{off}}}{\gamma + \delta} + \frac{d}{\gamma + \delta - r} > 1$$
Under this condition, the secular characteristic function $f(s) = 1 - \frac{2 k_{\text{off}}}{s} - \frac{d}{s - r} - \beta \frac{s - (\gamma + \delta)}{D(s)}$ satisfies $f(\gamma + \delta) < 0$ and $f(\mu_1^-) = +\infty$. By the Intermediate Value Theorem, $f(s)$ is guaranteed to have a root $s^* \in (\gamma + \delta, \mu_1)$. That root constitutes a physical relaxation mode whose amplitude satisfies $a(s^*) > 0$ unconditionally.

---

### Theorem 3 (Necessary Condition: Single-Open Infeasibility)

For any single-open 5-state Markov scheme ($\mathcal{O} = \{A_2O\}$), the open manifold has dimension 1 ($\mu_1 = \alpha, \delta = \gamma = 0$), and the gateway weight reduces to:
$$\mathcal{W}^{(k)}_{A_2R} = - (u^{(k)}_{A_2R})^2 \frac{\beta}{s_k - \alpha}$$
By Cauchy interlacing, the fastest relaxation rate $s_0$ strictly exceeds $\alpha$ ($s_0 > \alpha + \beta > \alpha$). Therefore, $\mathcal{W}^{(0)}_{A_2R} < 0$ unconditionally for all positive rate constants.

Furthermore, because $\dim(\mathcal{D}) = 3$, the system possesses at most 3 non-zero relaxation modes. Since mode 0 is strictly non-positive, **at most two relaxation modes can carry positive amplitude**. 

Consequently, observing three distinct positive decay components ($a_1 > 0, a_2 > 0, a_3 > 0$, as reported by Barberis et al. 2007: $0.57 \pm 0.04, 0.22 \pm 0.02, 0.21 \pm 0.03$) **strictly necessitates at least two conducting open states in the open manifold**.

---

## 2. Mathematical Proof

### Step 1: Symmetrized Jacobi Form and Eigenvector Projections
Order the states of $\mathcal{D}$ as $1: A_2D, 2: A_2R, 3: A_2O, 4: A_2O'$.
The subgenerator $B$ is tridiagonal:
$$B = \begin{pmatrix} -r & r & 0 & 0 \\ d & -(2 k_{\text{off}} + \beta + d) & \beta & 0 \\ 0 & \alpha & -(\alpha + \delta) & \delta \\ 0 & 0 & \gamma & -\gamma \end{pmatrix}$$
By detailed balance, $\pi = (d/r, 1, \beta/\alpha, (\beta/\alpha)(\delta/\gamma))$.
The similarity transformation $H = D B D^{-1}$ yields the real symmetric Jacobi matrix:
$$H = \begin{pmatrix} -r & \sqrt{dr} & 0 & 0 \\ \sqrt{dr} & -(2 k_{\text{off}} + \beta + d) & \sqrt{\alpha\beta} & 0 \\ 0 & \sqrt{\alpha\beta} & -(\alpha + \delta) & \sqrt{\delta\gamma} \\ 0 & 0 & \sqrt{\delta\gamma} & -\gamma \end{pmatrix}$$
Let $H u^{(k)} = -s_k u^{(k)}$, with $s_k > 0$ and $\|u^{(k)}\| = 1$.
The right eigenvectors of $B$ are $v^{(k)} = D^{-1} u^{(k)}$, and the left eigenvectors are $w^{(k)} = (u^{(k)})^T D$.
For the observable $P_O = p_3 + p_4$, the weight vector is:
$$\mathcal{W}^{(k)} = \left( u^{(k)}_3 \sqrt{\pi_3} + u^{(k)}_4 \sqrt{\pi_4} \right) D^{-1} u^{(k)}$$
Evaluating at the gateway state $j = 2$ (where $\pi_2 = 1$):
$$\mathcal{W}^{(k)}_{A_2R} = u^{(k)}_2 \left( u^{(k)}_3 \sqrt{\pi_3} + u^{(k)}_4 \sqrt{\pi_4} \right)$$

### Step 2: Exact Closed-Form Weight Ratio
From row 4 of $H u^{(k)} = -s_k u^{(k)}$:
$$\sqrt{\delta\gamma} u^{(k)}_3 - \gamma u^{(k)}_4 = -s_k u^{(k)}_4 \implies u^{(k)}_4 = -\frac{\sqrt{\delta\gamma}}{s_k - \gamma} u^{(k)}_3$$
Multiplying by $\sqrt{\pi_4} / \sqrt{\pi_3} = \sqrt{\delta/\gamma}$:
$$u^{(k)}_4 \sqrt{\pi_4} = -\frac{\delta}{s_k - \gamma} u^{(k)}_3 \sqrt{\pi_3}$$
Summing the conducting projections:
$$u^{(k)}_3 \sqrt{\pi_3} + u^{(k)}_4 \sqrt{\pi_4} = u^{(k)}_3 \sqrt{\pi_3} \left[ 1 - \frac{\delta}{s_k - \gamma} \right] = u^{(k)}_3 \sqrt{\pi_3} \left[ \frac{s_k - (\gamma + \delta)}{s_k - \gamma} \right]$$

From row 3 of $H u^{(k)} = -s_k u^{(k)}$:
$$\sqrt{\alpha\beta} u^{(k)}_2 - (\alpha + \delta) u^{(k)}_3 + \sqrt{\delta\gamma} u^{(k)}_4 = -s_k u^{(k)}_3$$
Substituting $u^{(k)}_4$:
$$\sqrt{\alpha\beta} u^{(k)}_2 = - \left[ s_k - \alpha - \delta - \frac{\delta\gamma}{s_k - \gamma} \right] u^{(k)}_3 = - \left[ s_k - \alpha - \frac{\delta s_k}{s_k - \gamma} \right] u^{(k)}_3$$
$$= - \frac{s_k(s_k - \gamma - \delta) - \alpha(s_k - \gamma)}{s_k - \gamma} u^{(k)}_3 = - \frac{D(s_k)}{s_k - \gamma} u^{(k)}_3$$
where $D(s) = s(s - \gamma - \delta) - \alpha(s - \gamma) = \det(s I - (-H)_{\{3,4\}})$.
Solving for $u^{(k)}_3$:
$$u^{(k)}_3 = - \frac{\sqrt{\alpha\beta}(s_k - \gamma)}{D(s_k)} u^{(k)}_2$$
Multiplying by $u^{(k)}_2$:
$$u^{(k)}_2 u^{(k)}_3 = - (u^{(k)}_2)^2 \frac{\sqrt{\alpha\beta}(s_k - \gamma)}{D(s_k)}$$
Substituting into $\mathcal{W}^{(k)}_{A_2R}$, and noting $\sqrt{\pi_3} = \sqrt{\beta/\alpha}$:
$$\mathcal{W}^{(k)}_{A_2R} = \left[ - (u^{(k)}_2)^2 \frac{\sqrt{\alpha\beta}(s_k - \gamma)}{D(s_k)} \right] \sqrt{\frac{\beta}{\alpha}} \left[ \frac{s_k - (\gamma + \delta)}{s_k - \gamma} \right]$$
The factor $(s_k - \gamma)$ cancels out identically, and $\sqrt{\alpha\beta} \sqrt{\beta/\alpha} = \beta$:
$$\mathcal{W}^{(k)}_{A_2R} = - (u^{(k)}_2)^2 \beta \frac{s_k - (\gamma + \delta)}{D(s_k)} \equiv - (u^{(k)}_2)^2 \beta \frac{s_k - (\gamma + \delta)}{(s_k - \mu_1)(s_k - \mu_2)}$$
This establishes Theorem 1(1).

### Step 3: Exact Secular Bracket Connection
The characteristic polynomial $\det(s I - (-H)) = 0$ expands by Schur complement on node 2:
$$\det(s I - (-H)) = (s - r) D(s) \left[ s - (2 k_{\text{off}} + \beta + d) - \frac{dr}{s - r} - \alpha\beta \frac{s - \gamma}{D(s)} \right] = 0$$
Using $d + \frac{dr}{s-r} = \frac{ds}{s-r}$ and $\beta + \frac{\alpha\beta(s-\gamma)}{D(s)} = \beta \frac{s(s-\gamma-\delta)}{D(s)}$:
$$s - 2 k_{\text{off}} - \frac{ds}{s - r} - \beta \frac{s(s - \gamma - \delta)}{D(s)} = 0$$
Dividing by $s > 0$:
$$1 - \frac{2 k_{\text{off}}}{s} - \frac{d}{s - r} \equiv \beta \frac{s - (\gamma + \delta)}{D(s)}$$
Substituting this directly into $\mathcal{W}^{(k)}_{A_2R}$:
$$\mathcal{W}^{(k)}_{A_2R} = - (u^{(k)}_2)^2 \left[ 1 - \frac{2 k_{\text{off}}}{s_k} - \frac{d}{s_k - r} \right]$$
Therefore:
$$\mathcal{W}^{(k)}_{A_2R} > 0 \iff 1 - \frac{2 k_{\text{off}}}{s_k} - \frac{d}{s_k - r} < 0 \iff \beta \frac{s_k - (\gamma + \delta)}{D(s_k)} < 0$$
Evaluating $D(\gamma + \delta)$:
$$D(\gamma + \delta) = (\gamma + \delta)(0) - \alpha(\gamma + \delta - \gamma) = -\alpha \delta < 0$$
Because $D(s) = (s - \mu_1)(s - \mu_2)$ and $D(\gamma + \delta) < 0$, $\gamma + \delta$ lies strictly between $\mu_2$ and $\mu_1$:
$$\mu_2 < \gamma + \delta < \mu_1$$
Therefore:
- In $(\gamma + \delta, \mu_1)$: $s_k - (\gamma + \delta) > 0$ and $D(s_k) < 0 \implies \frac{s_k - (\gamma + \delta)}{D(s_k)} < 0 \implies \mathcal{W}^{(k)}_{A_2R} > 0$.
- In $(0, \mu_2)$: $s_k - (\gamma + \delta) < 0$ and $D(s_k) > 0 \implies \frac{s_k - (\gamma + \delta)}{D(s_k)} < 0 \implies \mathcal{W}^{(k)}_{A_2R} > 0$.
- In $(\mu_2, \gamma + \delta)$: $s_k - (\gamma + \delta) < 0$ and $D(s_k) < 0 \implies \frac{s_k - (\gamma + \delta)}{D(s_k)} > 0 \implies \mathcal{W}^{(k)}_{A_2R} < 0$.
- In $(\mu_1, \infty)$: $s_k - (\gamma + \delta) > 0$ and $D(s_k) > 0 \implies \frac{s_k - (\gamma + \delta)}{D(s_k)} > 0 \implies \mathcal{W}^{(k)}_{A_2R} < 0$.
This establishes Theorem 1(2).

### Step 4: Proof of the Sufficient Condition
Consider the continuous function $f(s) = 1 - \frac{2 k_{\text{off}}}{s} - \frac{d}{s - r} - \beta \frac{s - (\gamma + \delta)}{D(s)}$ on the interval $(\gamma + \delta, \mu_1)$.
At the left endpoint $s = \gamma + \delta$:
$$f(\gamma + \delta) = 1 - \frac{2 k_{\text{off}}}{\gamma + \delta} - \frac{d}{\gamma + \delta - r}$$
If $\frac{2 k_{\text{off}}}{\gamma + \delta} + \frac{d}{\gamma + \delta - r} > 1$, then $f(\gamma + \delta) < 0$.
At the right endpoint as $s \to \mu_1^-$:
$(s - \mu_1) \to 0^-$ from below, while $s - (\gamma + \delta) > 0$ and $s - \mu_2 > 0$.
Thus $D(s) = (s - \mu_1)(s - \mu_2) \to 0^-$, so $-\beta \frac{s - (\gamma + \delta)}{D(s)} \to +\infty$.
Therefore $\lim_{s \to \mu_1^-} f(s) = +\infty$.
Since $f(s)$ is continuous on $(\gamma + \delta, \mu_1)$, $f(\gamma + \delta) < 0$, and $\lim_{s \to \mu_1^-} f(s) = +\infty$, by the Intermediate Value Theorem there exists at least one root $s^* \in (\gamma + \delta, \mu_1)$ such that $f(s^*) = 0$.
Because $s^* \in (\gamma + \delta, \mu_1)$, Theorem 1 guarantees $\mathcal{W}_{A_2R}(s^*) > 0$.
By Duhamel's formula, the mode amplitude satisfies $a(s^*) > 0$ unconditionally.
Q.E.D.
