# ! /usr/bin/env python
# -*- coding: utf-8 -*-
# ===================================================================
# eccwrapper.py
#
# 26.06.2019, Benjamin Kromer <benjamin.kromer1@web.de>
#
# @desc: Wrapper class for fastecdsa class. For Elliptic Curve and
#        Point representation and manipulation.
# ===================================================================
from fastecdsa.point import Point as FastecdsaPoint
from ECCGT import random_generator


class ShortPoint:
    """Elliptic Curve Point representation.

    Attributes:
        x (int): the x coordinate of the point
        y (int): the y coordinate of the point
    """
    x: int = None
    y: int = None

    def __init__(self):
        pass

    def __eq__(self, other):
        """Compare two elliptic curve points for equality.

        Args:
            other (ShortPoint): second elliptic curve point

        Returns:
            bool: True if points are the same, False else
        """
        if not isinstance(other, ShortPoint):
            # don't attempt to compare against unrelated types
            return NotImplemented

        return self.x == other.x and self.y == other.y


class Curve:
    name = None
    order = None
    generator = ShortPoint()

    def __init__(self):
        raise NotImplementedError('Abstract method __init__')


class Fastecdsa(Curve):
    """Wrapper class for fastecdsa library.

    Attributes:
        _curve: curve object from fastecdsa
        name: curve name
        generator: the base point of the curve
        order: the order of the base point of the curve
        rand_gen: random generator object for random numbers
    """
    def __init__(self, curve):
        """
        Args:
            curve: curve object from fastecdsa
        """
        self._curve = curve
        self.name = curve.name
        self.generator.x = self._curve.gx
        self.generator.y = self._curve.gy
        self.order = self._curve.q
        self.rand_gen = random_generator.RandomGenerator(self.order)

    def multiplication(self, k, P):
        """Multiply a elliptic curve point P by a integer k

        Args:
            k (int): integer in range of 1,...,order-1
            P (ShortPoint): elliptic curve point

        Returns:
            k*P (ShortPoint)
        """
        if P.x == 0 and P.y == 1:
            return P
        product = k * self.shortpoint_to_point(P)
        val = ShortPoint()
        val.x = product.x
        val.y = product.y
        return val

    def addition(self, P, Q):
        """Add two elliptic curve points P and Q

        Args:
            P (ShortPoint): elliptic curve point
            Q (ShortPoint): elliptic curve point

        Returns:
            P+Q (ShortPoint)
        """
        sum1 = self.shortpoint_to_point(P)
        sum2 = self.shortpoint_to_point(Q)
        sum1 = sum1 + sum2
        val = ShortPoint()
        val.x = sum1.x
        val.y = sum1.y
        return val

    def subtraction(self, P, Q):
        """Subtract two elliptic curve points P and Q

        Args:
            P (ShortPoint): elliptic curve point
            Q (ShortPoint): elliptic curve point

        Returns:
            P-Q (ShortPoint)
        """
        point1 = self.shortpoint_to_point(P)
        point2 = self.shortpoint_to_point(Q)
        diff1 = point1 - point2
        val = ShortPoint()
        val.x = diff1.x
        val.y = diff1.y
        return val

    def negation(self, P):
        """Negate elliptic curve point P

        Args:
            P (ShortPoint): elliptic curve point

        Returns:
            -P (ShortPoint)
        """
        point = self.shortpoint_to_point(P)
        neg = - point
        val = ShortPoint()
        val.x = neg.x
        val.y = neg.y
        return val

    def isoncurve(self, P):
        """Check if point P is on curve _curve

        Args:
            P (ShortPoint): elliptic curve point

        Returns:
            True if point is on curve, False else
        """
        point = FastecdsaPoint(P.x, P.y, self._curve)
        return self._curve.is_point_on_curve(point)

    def shortpoint_to_point(self, P):
        """Transform ShortPoint to fastecdsa point

        Args:
            P (ShortPoint): elliptic curve point

        Returns:
            fastecdsa point
        """
        if P.x == 0 and P.y == 1:
            return FastecdsaPoint(P.x, P.y, None)
        else:
            return FastecdsaPoint(P.x, P.y, self._curve)
