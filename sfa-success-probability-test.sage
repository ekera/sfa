# ------------------------------------------------------------------------------
# Performs basic unit tests for the sfa_success_probability() function.
#
# This function first tests that the success probability is correctly computed
# for all odd positive N < 10^4 with at least two distinct prime factors.
#
# It then tests that the success probability is correctly computed for ten
# integers N selected uniformly at random from [10^6, 2 * 10^6), again with the
# constraints that N must be odd and have at least two distinct prime factors.
#
# To check that the values yielded by sfa_success_probability() are correct,
# this function performs an exhaustive search over all integers g in [1, N) to
# count the number of g for which N is successfully factored.
def test_sfa_success_probability():
  # Exhaustively computes the success probability for small N for comparison.
  def exhaust_sfa_success_probability(N):
    R = IntegerModRing(N);

    success = 0;
    fail = 0;

    for g in range(1, N):
      if gcd(g, N) != 1:
        success += 1;
      else:
        g = R(g);
        r = g.multiplicative_order();

        if r % 2 == 0:
          if g^(r / 2) != -1:
            success += 1;
          else:
            fail += 1;
        else:
          fail += 1;

    return success / (success + fail);

  # Test all positive odd N < 10^4 with at least two distinct prime factors.
  print("Testing all odd N < 10^4 with at least two distinct prime factors:");

  for N in range(10^4):
    if N % 2 == 0:
      continue;

    if len(factor(N)) < 2:
      continue;

    print(" Testing N:", N);

    if exhaust_sfa_success_probability(N) != \
      sfa_success_probability(N, printouts = False):
      raise Exception("Error: Test failed for N:", N);

  # Tests ten integers N selected uniformly at random from all odd integers on
  # the interval [10^6, 2 * 10^6) with at least two distinct prime factors.
  print("\nTesting ten random N selected from the interval [10^6, 2 * 10^6):");

  for _ in range(10):
    while True:
      N = 10^6 + IntegerModRing(10^6).random_element().lift();

      if N % 2 == 0:
        continue;

      if len(factor(N)) >= 2:
        break;

    print(" Testing N:", N);

    if exhaust_sfa_success_probability(N) != \
      sfa_success_probability(N, printouts = False):
      raise Exception("Error: Test failed for N:", N);

  # Report success.
  print("\nAll tests passed.");