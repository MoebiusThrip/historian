# racquetballs.py to test particle in a box hypotheses

# import system tools
import os, sys
import json

# import time
from time import time, sleep

# import trig functions
from math import sin, cos, pi, exp, log

# import numpy
import math
import numpy
from numpy.random import random, rand, randint
from numpy import sin, cos, pi, exp, log, sqrt

# import matplotlib for plots
from matplotlib import pyplot
from matplotlib import gridspec
from matplotlib import colorbar
from matplotlib import colors as Colors
from matplotlib import style as Style
Style.use('fast')


# play racquetball
def racquet(steps=10000, frequency=2.0, left=0.0, right=20.0, teleportation=10.0):
    """Simulate a session of quantum racquetball.

    Arguments:
        steps: int, number of steps
        frequency: float, maximum random step distance
        left: float, left wll position
        right: float, right wall position
        teleportation: float, teleport distance

    Returns:
        None
    """

    # begin segments
    segments = []

    # set beginning at far left, with rightward movement
    first = left
    rightward = True

    # create random lengths
    lengths = numpy.random.rand(steps,) * frequency

    # create directions at random
    directions = numpy.random.choice(a=[False, True], size=(steps, ))

    # for each step
    for length, direction in zip(lengths, directions):

        # first non teleport, then teleport
        for teleport in (False, True):

            # change length
            if teleport:

                # change length
                length = teleportation
                rightward = direction

            # if rightward
            if rightward:

                # add
                second = first + length

                # check for crossing boundary
                if second > right:

                    # add first segment
                    segment = (first, right, teleport)
                    segments.append(segment)

                    # create other portion
                    portion = length - (right - first)
                    segment = (right, right - portion, teleport)
                    segments.append(segment)
                    rightward = False

                # otherwise
                else:

                    # just add the segment
                    segment = (first, second, teleport)
                    segments.append(segment)

            # otherwise
            else:

                # substart
                second = first - length

                # check for crossing boundary
                if second < left:

                    # add first segment
                    segment = (first, left, teleport)
                    segments.append(segment)

                    # create other portion
                    portion = length - (first - left)
                    segment = (left, left + portion, teleport)
                    segments.append(segment)
                    rightward = True

                # otherwise
                else:

                    # just add the segment
                    segment = (first, second, teleport)
                    segments.append(segment)

            # reassign first
            first = segments[-1][1]

    # use matplotlib
    pyplot.clf()
    for index, segment in enumerate(segments):

        # set color
        color = 'b-'
        if segment[2]:

            # set to red for teleport
            color = 'r--'

        # plot segment
        pyplot.plot([segment[0], segment[1]], [index, index], color)

    # saveplot
    pyplot.savefig('racquetball.png')

    # make histogram
    chunks = 1000
    size = (right - left) / chunks
    bins = [(right + (size * index), right + size + (size * index)) for index in range(chunks)]

    # got through each segment
    histogram = []
    for segment in segments:

        # only reals
        if not segment[2]:

            # determine number of bins
            number = (segment[0] - segment[1]) / size
            sign = numpy.sign((segment[0] - segment[1]))
            histogram += [segment[0] + index * size * sign for index in range(int(abs(number)))]

    # make histogram
    pyplot.clf()
    pyplot.hist(histogram, bins=100)
    pyplot.savefig('histogram.png')

    return None