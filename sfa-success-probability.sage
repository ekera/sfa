# ------------------------------------------------------------------------------
# Computes the probability of Shor's factoring algorithm — as it was originally
# described by Shor in [Shor94] [Shor97] — succeeding in finding non-trivial
# factor of a given integer N.
#
# More specifically, it is assumed that the below steps are performed:
#
#  1. Select an integer g uniformly at random from [1, N).
#
#     If d = gcd(g, N) != 1, return d as a non-trivial factor of N.
#
#     Otherwise, continue.
#
#     Note that if the algorithm continues, g will have been selected uniformly
#     at random from Z_N^*, the multiplicative group of the ring of integers
#     modulo N.
#
#  2. Perceive g as an element of Z_N^* and compute the order r of g.
#
#  3. If r is odd, the algorithm fails.
#
#  4. If g^(r / 2) == -1 (mod N), the algorithm fails.
#
#  5. Otherwise, the algorithm succeeds, and returns d = gcd(g^(r / 2) ± 1, N)
#     as the two non-trivial factors of N.
#
# This function computes the probability of the above algorithm succeeding for a
# given N, either immediately in step 1, or in step 5, having successfully
# passed steps 3 and 4.
#
# The script assumes that the quantum order finding function called in step 2
# correctly computes the order r. Failures due to an incorrect order being
# returned in step 2 are hence not accounted for when the probability of
# reaching step 5 is computed.
#
# Note that it is reasonable to assume that the order will in general be
# correctly computed. (See e.g. [E24] and Appendix A to [E20].)
#
# Note furthermore that there are more efficient methods of post-processing the
# order r returned by the order finding algorithm (See e.g. [E21b] and [E24].)
#
# This script is provided only to answer commonly posed questions regarding the
# success probability of Shor's original factoring algorithm; not to propose
# that the algorithm be used in practice as it was originally described.
#
# The arguments to this function are as follows:
#
#  N:          The integer to be factored. It is required that N be an odd
#              positive integer and not a perfect prime power.
#
#  factors:    Either the prime factors of N = p1^e1 * .. * pn^en, on the
#              form [[p1, e1], .., [pn, en]], or None, as is the default.
#
#              If set to None, the function will factor N. Note however that
#              this is only computationally feasible for small to moderate size
#              N, or easily factorizable N.
#
#  printouts:  A boolean flag that may be set to True to print intermediary
#              success probabilities of passing the above steps. Defaults to
#              False for no printouts.
#
# This function returns the success probability.
#
# For more detailed information, please see instead the README.md file.
#
# Bibliography:
#
#  [E20]     Ekerå, M.: Quantum algorithms for computing general discrete
#            logarithms and orders with tradeoffs. J. Math. Cryptol. 15(1),
#            pp. 359–407 (2020).
#            https://doi.org/10.1515/jmc-2020-0006
#
#  [E21b]    Ekerå, M.: Ekerå, M. On completely factoring any integer
#            efficiently in a single run of an order-finding algorithm.
#            Quantum Inf. Process 20:205 (2021).
#            https://doi.org/10.1007/s11128-021-03069-1
#
#  [E24]     Ekerå, M.: On the success probability of quantum order finding.
#            ACM Trans. Quantum Comput. 5(2):11 (2024).
#            https://doi.org/10.1145/3655026
#
#  [Shor94]  Shor, P.W.: Algorithms for Quantum Computation: Discrete Logarithms
#            and Factoring. In: SFCS, Proceedings of the 35th Annual Symposium
#            on Foundations of Computer Science, IEEE Computer Society,
#            Washington, DC, pp. 124–134 (1994).
#            https://doi.org/10.1109/SFCS.1994.365700
#
#  [Shor97]  Shor, P.W.: Polynomial-time algorithms for prime factorization and
#            discrete logarithms on a quantum computer. SIAM J. Comput. 26(5),
#            pp. 1484–1509 (1997).
#            https://doi.org/10.1137/S0097539795293172
def sfa_success_probability(N, factors = None, printouts = False):
  # A real field with 32 bit precision used exclusively for printouts.
  F = RealField(32);

  # Sanity checks to make sure the function is called correctly.
  if (N % 2) == 0:
    raise Exception("Error: The integer N must be odd.");

  if factors == None:
    factors = factor(N); # Compute the factorization if it is not given.
  else:
    # Sanity check the given factorization.
    if len(factors) != len(set([pi for [pi, ei] in factors])):
      raise Exception("Error: A factor was twice specified.");

    for [pi, ei] in factors:
      if not pi.is_prime(proof = False):
        raise Exception("Error: The factor " + str(pi) + " is not prime.");

      if ei < 1:
        raise Exception("Error: All exponents must be positive integers.");

    if prod([pi^ei for [pi, ei] in factors]) != N:
      raise Exception("Error: The given factorization of N is incorrect.");

  # Extract n,.
  n = len(factors);

  # Sanity check n.
  if n < 2:
    raise Exception("Error: The integer N is prime or a pure prime power.");

  # The total number of elements in the group Z_N^* is phi(N).
  phi = prod([pi^(ei - 1) * (pi - 1) for [pi, ei] in factors]);

  # The algorithm selects g uniformly at random from [1, N). There are N - 1
  # integers on this interval. The algorithm may by chance select g such that
  # d = gcd(g, N) != 1, in which case d is a non-trivial factor of N.
  if printouts:
    print("1. The probability of factoring via a non-trivial GCD is:", \
      1 - phi / (N - 1), "≈ 10^" + str(F(log(1 - phi / (N - 1))/log(10))));

  # The main procedure follows below:
  #
  # It computes the number of elements g in Z_N^* with even order r, and the
  # number of good elements with even order r such that g^(r/2) != -1.

  even_elements = 0; # The number of elements g in Z_N^* with even order r.
  good_elements = 0; # The number of elements with even r and g^(r/2) != -1.

  # Selecting g uniformly at random from Z_N^* is equivalent to selecting g_i
  # uniformly at random from Z_{p_i^e_i}^* and then computing g using the
  # Chinese remainder theorem. As Shor points out on p. 1498 of [Shor97], if
  # r_i is the multiplicative order of g_i, then r = lcm(r_1, .., r_n) is the
  # order of g, and Shor's algorithm only fails if 2^t for some t is the
  # largest power of two to divide all of r_1, .., r_n.

  # The below function computes the probability of an element g_i selected
  # uniformly at random from Z_{p_i^e_i}^* having an order r_i that is divisible
  # by 2^{t_i}. It returns a dictionary of Pr((p_i, e_i), t_i) indexed by ti,
  # where t_i runs from 0 up to and including kappa_i, for 2^{kappa_i} the
  # largest power of two to divide \phi(p_i^e_i).
  def tabulate_probabilities_for_all_ti(pi, ei):
    order = pi^(ei - 1) * (pi - 1); # = \phi(p_i^e_i).

    # Let kappa_i be such that order = 2^kappa_i o for o odd.
    #
    # Note: This implies that kappa_i >= 1, since pi - 1 divides order, where
    #       p_i is an odd prime, so the order must be even.

    kappai = 0;
    while (order % 2^(kappai + 1)) == 0:
      kappai += 1;

    results = dict();
    results[kappai] = 1/2;

    for ti in range(1, kappai):
      results[kappai - ti] = (1/2) * (1/2^ti);

    results[0] = 1/2^kappai;

    return results;

  # Tabulate probabilities in t_i for each subgroup modulo p_i^e_i.
  tabulated_probabilities_for_all_ti = \
    [tabulate_probabilities_for_all_ti(pi, ei) for [pi, ei] in factors];

  # Compute the number of possible (t_1, .., t_n) tuples.
  combinations_of_tis = \
    prod([len(table) for table in tabulated_probabilities_for_all_ti]);

  # Compute the number of elements g in Z_N^* with odd and even order r.
  odd_elements = \
    phi * prod([table[0] for table in tabulated_probabilities_for_all_ti]);
  even_elements = phi - odd_elements;

  # Compute the number of good and bad elements g in Z_N^*, where the element is
  # said to be good iff Shor's algorithm splits N when g is selected.
  max_ti = \
    min([max(table.keys()) for table in tabulated_probabilities_for_all_ti]);

  bad_elements = 0;
  for i in range(max_ti + 1):
    bad_elements += \
      phi * prod([table[i] for table in tabulated_probabilities_for_all_ti]);
  good_elements = phi - bad_elements;

  # Compute the success probability.
  probability = (1 - phi / (N - 1)) + (phi / (N - 1)) * (good_elements / phi);

  if printouts:
    print("\nAssume this *did not* happen:");

    print("\n 2. Assume the order r of g is correctly computed:");

    print("\n  3. The probability of selecting g of even order r is:", \
      even_elements / phi);

    print("\n  Assume this *did* happen:");

    print("\n   4. The probability of g being such that g^(r/2) != -1 is:", \
      good_elements / even_elements);

    print("\n   Assume this *did* happen:");

    print("\n    5. Non-trivial factors are produced by gcd(g^(r/2) ± 1, N).");

    print("\nThe total success probability is hence:", \
      probability, "≈", F(probability));

    print("");

  # Return the probability.
  return probability;