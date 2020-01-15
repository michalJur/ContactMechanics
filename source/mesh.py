"""
Created at 29.10.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

from source.point import Point
from source.edge import Edge
from source.element import Element
from source.smart_array import SmartArray
import numpy as np


class Mesh:

    # TODO: rename point -> points etc.
    def __init__(self, point: Point, edge: Edge, element: Element,  subarea: dict):
        self.point = point
        self.edge = edge
        self.element = element
        self.subarea = subarea
        # self.independent = SmartArray(self.point, idx)
