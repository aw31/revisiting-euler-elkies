# Revisiting $A^4 + B^4 + C^4 = D^4$ in 2023

In 1769, Euler conjectured that the equation $A^4 + B^4 + C^4 = D^4$ has no non-trivial integer solutions, in a similar spirit to Fermat's last theorem. Two hundred years later, [Elkies (1988)](https://www.ams.org/journals/mcom/1988-51-184/S0025-5718-1988-0930224-9/S0025-5718-1988-0930224-9.pdf) refuted this conjecture using elliptic curve theory, and [Frye (1988)](https://ieeexplore.ieee.org/document/74138) found the smallest counterexample through a supercomputer search that took over 100 hours. With today's hardware and an improved algorithm, we show that Euler's conjecture can be refuted in seconds using a consumer laptop and a smidgen of elementary number theory.

## Algorithm

To solve the Diophantine equation $A^4 + B^4 + C^4 = D^4$, we implement an efficient version of the approach sketched in [Elkies (1988)](https://www.ams.org/journals/mcom/1988-51-184/S0025-5718-1988-0930224-9/S0025-5718-1988-0930224-9.pdf). (A naive search would need $\sim 4\cdot 400000^4\approx 10^{23}$ integer operations before finding a solution, far exceeding even the limits of modern hardware.)

Our algorithm performs a "meet-in-the-middle" search over pairs $(A, B)$ and $(C, D)$. It runs in two phases, each taking quadratic time in the upper bound on $D$. This improves over both the naive quartic time algorithm and the cubic time algorithm deployed by [Frye (1988)](https://ieeexplore.ieee.org/document/74138).

The two phases of our algorithm are as follows:

1. Generate a list of "candidate differences" $D^4 - C^4$ by iterating over pairs $(C, D)$ and pruning pairs that fail to satisfy certain modular constraints. Specifically, we consider the equation modulo $5$, $2^8$, $3^6$, $13^2$, and $29^2$.
2. Check pairs $(A, B)$ until we find a pair such that $A^4 + B^4$ equals a candidate difference $D^4 - C^4$. For fast, cache-friendly lookups, we store candidate differences using a Bloom filter in front of a hash map.

Running this algorithm uncovers the solution $414560^4 + 95800^4 + 217519^4 = 422481^4$.

## Usage

```
git clone https://github.com/aw31/revisiting-euler-elkies.git
cd revisiting-euler-elkies && make && ./a.out
```

## Output

On a 2020 MacBook Pro, the program (with parallelized search) produces the following output:
```
clang++ -std=c++20 -fopenmp -Ofast *.cc -lomp
Searching up to D = 500000 with 16 threads

Found 1764000 good pairs (0.154864%)
Found 26415362 candidate differences (0.0105661%)

=== Compute differences ===
Time: 4.126s
Total: 4.126s

=== Populate Bloom filter and hash map ===
Time: 2.442s
Total: 6.568s

=== Check pairwise sums ===
Time: 10.174s
Total: 16.742s

Solution found: 414560^4 + 95800^4 + 217519^4 = 422481^4
```
