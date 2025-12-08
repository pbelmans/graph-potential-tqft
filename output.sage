load("tqft.sage")

for g in range(2, 11):
    print(f"Genus {g} and odd parity")
    print(regularized_period(g, 1))
    print(f"Genus {g} and even parity")
    print(regularized_period(g, 0))
    print()
# note that no caching of powers of `A` is done between computations for different genus
