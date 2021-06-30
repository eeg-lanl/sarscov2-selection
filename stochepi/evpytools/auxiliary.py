#!/usr/bin/env python3
# coding: utf-8


def unzip(xs):
    return [x[0] for x in xs], [x[1] for x in xs]



def unique(seq):
    """preserves the order (unlike list(set(seq)))"""
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def flatten(xss):
    if len(xss) > 0 and type(xss[0]) is list:
        return flatten([x for xs in xss for x in xs])
    else:
        return xss