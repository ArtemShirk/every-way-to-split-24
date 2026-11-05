# Every Way to Split 24

Part of **Collapse Theory**. Watch the series: [Collapse Theory playlist](https://www.youtube.com/watch?v=E33IrIBLFfA&list=PLbrl13LRQB2U) · [@shirkvids on YouTube](https://www.youtube.com/@shirkvids)

Split 24 into two whole numbers. Use them to stretch two copies of Euler's pentagonal product, multiply, and read the numbers that come out. At every prime, that number tells you whether the prime fits a simple formula.

## The rule

Pick any `a` from 1 to 12 and let `b = 24 − a`. Multiply out

```
q · (1 − qᵃ)(1 − q²ᵃ)(1 − q³ᵃ)··· · (1 − qᵇ)(1 − q²ᵇ)(1 − q³ᵇ)···
```

For a prime `p` that doesn't divide `a·b`, look at the number in front of `qᵖ`:

| Number | Meaning |
|---|---|
| **2** | p = x² + axy + 6ay² for some whole numbers x, y |
| **−1 or −2** | p = 4x² + (a−4)xy + (a+1)y² instead |
| **0** | p fits neither |

The same two formulas work for all twelve splits. You only change `a`.

**Example.** Split 24 as 6 + 18. The number in front of q³¹ is 2, and 31 = 2² + 27·1². (With a = 6 the first formula is x² + 6xy + 36y², which is x² + 27y² rearranged.)

For mathematicians:

$$\eta(az)\,\eta\big((24-a)z\big)=\tfrac12\left(\theta_{[1,\,a,\,6a]}-\theta_{[4,\,a-4,\,a+1]}\right)$$

## The twelve splits

| Split | Primes with a 2 |
|---|---|
| 1 + 23 | 59, 101, 167, 173, 211 |
| 2 + 22 | 47, 53, 103, 163, 199 |
| 3 + 21 | 67, 79, 127, 163, 277 |
| 4 + 20 | 29, 89, 101, 181, 229 |
| 5 + 19 | 131, 239, 389, 419, 461 |
| 6 + 18 | 31, 43, 109, 127, 157 |
| 7 + 17 | 263, 443, 557, 701, 1019 |
| 8 + 16 | 41, 113, 137, 257, 313 |
| 9 + 15 | 139, 151, 199, 331, 541 |
| 10 + 14 | 71, 149, 179, 331, 359 |
| 11 + 13 | 179, 467, 653, 719, 797 |
| 12 + 12 | 37, 61, 157, 193, 313 |

## Where it came from

This comes out of **Collapse Theory**, a framework I'm developing about how systems of rules break and rebuild. I explain it in my [Collapse Theory video series](https://www.youtube.com/watch?v=E33IrIBLFfA&list=PLbrl13LRQB2U) on [my YouTube channel](https://www.youtube.com/@shirkvids). Its test bed is the **Oasis**: a system takes in one newcomer per day. When too many newcomers can't be explained, the system collapses, and it has to find a new rule that keeps everything it already knew.

The **Prime Oasis** runs that protocol on multiplication. It finds the primes as its own collapses and hits the failure of unique factorization in Z[√−5] on day 33. Later it runs out of remainder-based rules and climbs to Gauss's cube test for x² + 27y². Asked to search every way of splitting 24, it found the rule above.

## Status

- **Proved.** Both sides are modular forms of the same type, and for these forms matching the first 24 coefficients forces a match forever. We checked 100,000.
- **Tested** at every prime below one million for all twelve splits: 941,953 checks, 0 failures.
- **Possibly known.** Two cases are classical: 1 + 23 (van der Blij, 1952) and 6 + 18 (Gauss's cube test). I haven't found all twelve stated together, or the single formula for the second form. If you know a reference, please open an issue.

## Run it

```
pip install numpy scipy mpmath
python sim/run_all.py
```

Every stage runs in about 10 seconds. To run one stage: `python sim/run_all.py 24` (or `1`, `2`, `3`, `predict`, `zeros`).

## What's here

| Folder | Contents |
|---|---|
| `paper/` | The proof, with a plain-language half for non-mathematicians |
| `report/` | The Prime Oasis visual report |
| `sim/` | The simulator, every stage, and the result data |

## Credit

**Artem Shirk**, more at [youtube.com/@shirkvids](https://www.youtube.com/@shirkvids).

**Claude (Anthropic)** wrote the simulator code, ran the computations and helped draft the proof, working from Artem's framework.

## Cite

Shirk, A. C. (2026). *Every Way to Split 24: Twelve Eta Products and the Principal Form.* https://github.com/artemshirk/every-way-to-split-24

## License

Code (`sim/`): MIT, see `LICENSE`.
Paper, report and text: CC BY 4.0, see `LICENSE-DOCS`. Reuse is welcome with credit to Artem C. Shirk.
