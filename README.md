# Montgomery

## Overview

The "Montgomery" project offers a comprehensive library for performing arithmetic operations within the Montgomery domain, particularly over a prime field. This project aims to support efficient hardware (HW) implementations of these operations, which include addition, subtraction, multiplication, division, and inversion. The focus is on educational purposes and illustration, and it is **not intended for use in production environments under any circumstances**.


## Montgomery Domain

The Montgomery domain refers to a transformation of integers modulo `N`, usually a prime, into a different representation which can make certain modular arithmetic operations, particularly modular multiplication, more efficient to compute, especially on hardware platforms.

### Montgomery Representation

A number `x` in the conventional representation is transformed into the Montgomery domain by computing `xR mod N`, where `R` is typically a power of 2 greater than `N`. This transformation is beneficial because it allows modular reductions to be performed more efficiently.

### Operations in the Montgomery Domain

- **Addition/Subtraction**: These operations are performed as usual followed by a modular reduction.
- **Multiplication**: Instead of performing a standard modular multiplication, the numbers are first transformed 
  into their Montgomery forms, multiplied, and then transformed back, utilizing the Montgomery reduction process. 
  Note: single operations inside the Montgomery domain might not be beneficial taken consideration the cost of field 
  changing.
- **Division/Inversion**: These operations involve complex steps using extended Euclidean algorithms, adapted to the 
  Montgomery form for efficiency. We adopt a method described in reference \[1\] here to boost its performance.

## Installation

```bash
git clone https://github.com/r5by/montgomery.git
cd montgomery
pip install -r requirements.txt
```

## Example Usage

```python

# Example Montgomery domain with realization 5 as the multiplication (refer [2]) operator
M = Montgomery.factory(mod=31, mul_opt='real5').build(m=5)

x_, y_ = 7, 21  # x', y' as elements in prime field F_{31}
# enter the Montgomery domain over F_{31}
x, y = M(7), M(21)
a = x + y
b = x - y
c = x * y
d = 1 / x
e = x / y

# exit the Montomery domain
_x, _y = int(x), int(y)
assert _x == x_
assert _y == y_
```

## Known Issues

Realization 8 in reference [2] only works for sequential compressor not ternary one, finding out why...

## References
1. [Montgomery inversion](https://link.springer.com/article/10.1007/s13389-017-0161-x)
2. [High-Radix Design of a Scalable Montgomery Modular Multiplier With Low Latency](https://ieeexplore.ieee.org/abstract/document/9328560)
3. [Topics in Computational Number Theory Inspired by Peter L. Montgomery](https://www.cambridge.org/core/books/topics-in-computational-number-theory-inspired-by-peter-l-montgomery/4F7A9AE2CE219D490B7D253558CF6F00)

