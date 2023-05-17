import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
import glider_sizing

import numpy as np

def test_balloon_mass():
	assert 1 > glider_sizing.balloon_mass(1.361) > 1.5