# On computing the probability of Shor's original factoring algorithm splitting a given integer of known factorization
This repository contains a [script](sfa-success-probability.sage) for computing the probability of Shor's factoring algorithm — as it was originally described by Shor in [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) — succeeding in splitting a given odd composite integer $N$ that is not a perfect prime power into two non-trivial factors, given the order $r$ of an element $g$ selected uniformly at random from $\mathbb Z_N^\ast$.
To perform the computation, the script requires the complete factorization of $N$ to be known.

This repository furthermore contains a [test script](sfa-success-probability-test.sage) that may be used to verify the correctness for small $N$.

The aforementioned scripts were originally written back in 2020 with students taking introductory lectures to Shor's algorithms in mind.
They are made public in the hope that they may be useful also to others.
They are distributed "as is" without warranty of any kind, either expressed or implied.
For further details, see the [license](LICENSE.md).

## Understanding why these scripts are of limited practical relevance
In practice, one would not factor integers using Shor's algorithm as originally described in [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700).
This is because there are better ways to factor integers quantumly:

- In [[E21b]](https://doi.org/10.1007/s11128-021-03069-1), it is shown that <i>all</i> prime factors of any integer $N$ can be found efficiently classically with very high probability given the order $r$ of an element $g$ selected uniformly at random from $\mathbb Z_N^\ast$.
This by  using a more involved classical post-processing algorithm that is essentially due to Miller [[Miller76]](https://doi.org/10.1016/S0022-0000(76)80043-8).

   More specifically, a lower bound on the probability of completely factoring $N$ is given in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1): Asymptotically, in the limit as $N$ tends to infinity, the lower bound on the success probability tends to one if properly parameterized.
   Already for moderate $N$, a very high success probability can be guaranteed.

   For a reference implementation of the classical post-processing algorithm, see [this repository](https://github.com/ekera/factoritall).

- In [[E24]](https://doi.org/10.1145/3655026), the lower bound in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1) is extended to also encompass the quantum order-finding part of the algorithm, so that only the availability of a fault-tolerant quantum computer has to be assumed.

   More specifically, a lower bound on the success probability of first finding the order $r$ of $g$, and of then completely factoring $N$ given $r$, is given in [[E24]](https://doi.org/10.1145/3655026): Again, asymptotically in the limit as  $N$ tends to infinity, the lower bound on the success probability tends to one if properly parameterized.
   Already for moderate $N$, a high success probability can be guaranteed.
   This conclusively shows that a single run of the quantum part is usually sufficient.

   This fact may also be seen in simulations, such as those in App. A of [[E21]](https://doi.org/10.1515/jmc-2020-0006).

- The above results are for factoring <i>any</i> integer $N$ into all of its prime factors.

   An important special class of integers in cryptography are RSA integers $N = pq$ with two large distinct random prime factors $p$, $q$ of similar bit lengths.
   To factor such integers, it is better to use Ekerå-Håstad's derivative of Shor's algorithm [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20) with the post-processing in [[E20]](https://doi.org/10.1007/s10623-020-00783-2) or [[E23p]](https://doi.org/10.48550/arXiv.2309.01754).
   A single run of the quantum part is then also usually sufficient, but this single run imposes less requirements on the quantum computer.

- For earlier works, see also e.g.
   [Knill95],
   [[McAnally01]](https://doi.org/10.48550/arXiv.quant-ph/0112055),
   [[Seifert01]](https://doi.org/10.1007/3-540-45353-9_24),
   [[Leander02]](https://doi.org/10.48550/arXiv.quant-ph/0208183),
   [[Gerjuoy05]](https://doi.org/10.1119/1.1891170),
   [[BW07]](https://doi.org/10.26421/QIC7.5-6-7),
   [[M-L+12]](https://doi.org/10.1038/nphoton.2012.259),
   [[Lawson15]](https://doi.org/10.1007/s11128-014-0910-z),
   [[GLMS15]](https://doi.org/10.48550/arXiv.1511.04385) and
   [[MRS18]](https://doi.org/10.48550/arXiv.1802.08444),
   along with e.g.
   [[Pollard74]](https://doi.org/10.1017/S0305004100049252),
   [[Miller76]](https://doi.org/10.1016/S0022-0000(76)80043-8),
   [[Rabin80]](https://doi.org/10.1016/0022-314X(80)90084-0) and
   [Long].

   For more comprehensive literature surveys, see e.g.
   [[E21b]](https://doi.org/10.1007/s11128-021-03069-1), Sect. 2,
   and
   [[E24]](https://doi.org/10.1145/3655026), Sect. 1.4–1.5.

## Prerequisites
To install [Sage](https://www.sagemath.org) under [Ubuntu 20.04 LTS](https://releases.ubuntu.com/20.04), simply execute:

```console
$ sudo apt install sagemath
```
For other Linux and Unix distributions, or operating systems, you may need to [download Sage](https://www.sagemath.org/download) and install it manually.
These scripts were originally developed for an older version of Sage, but then have since been updated and tested to work with Sage 9.4.

### Attaching the scripts
Launch Sage and attach the main script [<code>sfa-success-probability.sage</code>](sfa-success-probability.sage) by executing:

```console
$ sage
(..)
sage: %attach sfa-success-probability.sage
```

You may also attach the test script [<code>sfa-success-probability-test.sage</code>](sfa-success-probability-test.sage) by executing:

```console
sage: %attach sfa-success-probability-test.sage
```

## Computing the success probability
In what follows, let
$N = \prod_{i = 1}^n p_i^{e_i}$,
where
$n \ge 2$,
the
$p_i$
are pairwise distinct odd prime factors, and the
$e_i$
are positive integer exponents.
All odd integers that are not perfect prime powers can be written on this form.

To compute the probability of Shor's original algorithm succeeding in factoring $N$, when implemented exactly as described in the next subsection, call <code>sfa_success_probability(N, factors = None, printouts = False)</code>.

The arguments to this function are as follows:

- <code>N</code> — The integer $N$ to be factored.

  It is required — as in [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) — that $N$ be an odd positive integer and not a perfect prime power.

- <code>factors</code> — Either the prime factors of $N = \prod_{i = 1}^n p_i^{e_i}$, on the form <code>[[p1, e1], .., [pn, en]]</code>, or <code>None</code> in which case $N$ will be factored using functions native to Sage.
Note that using the latter option is only computationally feasible for small to moderate size $N$, or otherwise easily factorizable N.

- <code>printouts</code> — A boolean flag that may be set to <code>True</code> to print intermediary success probabilities of passing the above steps.
Defaults to <code>False</code> for no printouts.

The function returns the exact success probability as a rational.

### Assumptions
More specifically, the <code>sfa_success_probability()</code> function assumes that the below steps are performed:

1. Select an integer $g$ uniformly at random from $[1, N)$.

   If $d = \text{gcd}(g, N) \neq 1$, return $d$ as a non-trivial factor of $N$. Otherwise, continue.

   Note that if the algorithm continues, $g$ will have been selected uniformly at random from $\mathbb Z_N^\ast$, the multiplicative group of the ring of integers modulo $N$. (To have this property, and to keep the analysis as simple as possible, we do not exclude $1$ from the interval. It is easy to adapt the analysis to instead sample from $(1, N)$  if desired.)

2. Perceive $g$ as an element of $\mathbb Z_N^\ast$ and compute the order $r$ of $g$.

3. If $r$ is odd, the algorithm fails.

4. If $g^{r / 2} \equiv -1 \quad (\text{mod } N)$, the algorithm fails.

5. Otherwise, the algorithm succeeds, and returns $d = \gcd(g^{r / 2} \pm 1, N)$ as the two non-trivial factors of $N$.

The <code>sfa_success_probability()</code> function computes the probability of the above algorithm succeeding for a given $N$; either immediately in step 1, or in step 5, having successfully passed steps 3 and 4.

The function furthermore assumes that the quantum order finding function called in step 2 correctly computes the order $r$.
Failures due to an incorrect order being returned in step 2 are hence not accounted for when the probability of reaching step 5 is computed.

Note that it is reasonable to assume that $r$ will in general be correctly computed:

See for example App. A in [[E21]](https://doi.org/10.1515/jmc-2020-0006) for estimates based on simulations that show that enumerating a limited set of vectors in a lattice will yield the order, or the recent analysis in [[E24]](https://doi.org/10.1145/3655026) that yields a lower bound.

### Test suites
To test the [main script](sfa-success-probability.sage), attach the [test script](sfa-success-probability-test.sage) and call <code>test_sfa_success_probability()</code>.

The <code>test_sfa_success_probability()</code> function tests that the success probability is correctly computed for
- all positive odd integers $N < 10^4$ that are not perfect prime powers, and
- ten integers $N$ selected uniformly at random from $[10^6, 2 \cdot 10^6)$, again with the constraints that $N$ must be odd and not a perfect prime power.

To check that the values yielded by the <code>sfa_success_probability()</code> function are correct, the test function <code>test_sfa_success_probability()</code> performs an exhaustive search over all integers $g$ in $[1, N)$ to count the number of $g$ for which $N$ is successfully factored by <a href="#Assumptions">the above algorithm</a>.

### Understanding how the script works
To understand how the script works, let us consider each step in <a href="#Assumptions">the above algorithm</a>:

1. There are
   $N-1$
   elements on the interval
   $[1, N)$ in step 1 of the algorithm.

   Of these,
   $N - 1 - \phi(N)$
   elements are not in
   $\mathbb Z_N^\ast$,
   and hence immediately yield a non-trivial divisor of
   $N$.

   The probability of a non-trivial divisor being produced in this manner is
   $1 - \phi(N) / (N - 1)$.
   Otherwise — i.e. with probability
   $\phi(N) / (N - 1)$
   — the algorithm will have succeeded in selecting
   $g$
   uniformly at random from
   $\mathbb Z_N^\ast$.

   Above, and in what follows below,
   $\phi$
   is Euler's totient function.

   - <i>Before continuing, it is useful to introduce some more notation:</i>

      Selecting
      $g$
      uniformly at random from
      $\mathbb Z_N^\ast$
      is equivalent to selecting
      $g_i$
      uniformly at random from
      $\mathbb Z_{p_i^{e_i}}^\ast$
      for all
      $i \in [1, n]$,
      and then computing
      $g$
      from
      $(g_1, \ldots, g_n)$
      by requiring that
      $g_i \equiv g \quad (\text{mod } p_i^{e_i})$
      for all
      $i \in [1, n]$
      and applying the Chinese remainder theorem.

      For
      $r_i$
      the order of
      $g_i$,
      it then furthermore holds that
      $r = \mathrm{lcm}(r_1, \ldots, r_n)$
      is the order of
      $g \simeq (g_1, \ldots, g_n)$.

2. By assumption, the order $r$ of $g \in \mathbb Z_N^\ast$ is successfully computed in step 2 of the algorithm.

   - <i>Before continuing, it is useful to introduce some more notation:</i>:

      Let
      $2^{\kappa_i}$
      be the greatest power of two to divide
      $\phi(p_i^{e_i})$.

      Let
      $t_i \in [0, \kappa_i]$
      be such that
      $2^{t_i}$
      is the greatest power of two to divide
      $r_i$.

   - For each
     $i \in [1, n]$,
     the [main script](sfa-success-probability.sage) then tabulates the probability
     $\mathrm{Pr}((p_i, e_i), t_i)$
     of selecting an element
     $g_i \equiv g \quad (\text{mod } p_i^{e_i})$
     of order
     $r_i$
     such that
     $t_i = 0, 1, \ldots, \kappa_i$.

3. If $r$ is odd, the algorithm fails in step 3.

   - This occurs if and only if
     $t_1 = \ldots = t_n = 0$,
     since
     $r_1, \ldots, r_n$
     are then all odd, so
     $r = \mathrm{lcm}(r_1, \ldots, r_n)$
     is odd.

   - The probability of this occuring is hence $\prod_{i = 1}^n \mathrm{Pr}((p_i, e_i), t_i = 0)$.

4. If $g^{r / 2} \equiv -1 \quad (\text{mod } N)$, the algorithm fails in step 4.

   - This occurs if and only if
     $t_1 = \ldots = t_n \neq 0$.

     When combined with 3, this implies that Shor's algorithm fails to split
     $N$
     if and only if
     $t_1 = \ldots = t_n$.

   - The probability of this, or the situation in 3, occuring, is hence
     $\sum_{(t_1, \ldots, t_n)} \prod_{i = 1}^n \mathrm{Pr}((p_i, e_i), t_i)$,
     where the sum runs over all tuples
     $(t_1, \ldots, t_n)$
     such that
     $t_1 = \ldots = t_n$,
     and where we recall that each
     $t_i \in [0, \kappa_i]$.

5. Otherwise, the algorithm succeeds.

### Notes on the relation to published works
Note that the above procedure for determining whether a given fixed element $g$ will split $N = \prod_{i=1}^n p_i^{e_i}$ is not new: It is described by Shor in the journal version of his original paper as a part of the analysis, see in particular the third paragraph on p. 1498 of [[Shor97]](https://doi.org/10.1137/S0097539795293172). The idea is as follows:

For each
$i \in [1, n]$,
let
$2^{t_i}$
be the greatest power of two to divide the order
$r_i$
of
$g_i = g \quad (\text{mod } p_i^{e_i})$.

Then Shor's algorithm fails to split
$N = \prod_{i=1}^n p_i^{e_i}$
if and only if
$t_1 = \ldots = t_n$.

The [main script](sfa-success-probability.sage) in this repository simply implements this procedure from [[Shor97]](https://doi.org/10.1137/S0097539795293172).

A similar extended analysis was also recently given in [[BK21]](https://doi.org/10.1016/j.procs.2021.11.020), but the difference compared to [[Shor97]](https://doi.org/10.1137/S0097539795293172) is marginal.

Compared to [[Shor97]](https://doi.org/10.1137/S0097539795293172) and [[BK21]](https://doi.org/10.1016/j.procs.2021.11.020), the [main script](sfa-success-probability.sage) in this repository computes the probability of success when $g$ is selected uniformly at random from $\mathbb Z_N^\ast$, as in Shor's algorithm — not only for a given fixed $g$ in $\mathbb Z_N^\ast$.
It does so for any odd $N$ that is not a perfect prime power, and for which the factorization is known.

## About and acknowledgments
These scripts were developed by [Martin Ekerå](mailto:ekera@kth.se), in part at [KTH, the Royal Institute of Technology](https://www.kth.se/en), in Stockholm, [Sweden](https://www.sweden.se).
Valuable comments and advice were provided by Johan Håstad throughout the development process.

Funding and support was provided by the Swedish NCSA that is a part of the [Swedish Armed Forces](https://www.mil.se).