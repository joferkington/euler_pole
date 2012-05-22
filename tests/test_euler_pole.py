import itertools
import numpy as np
import euler_pole

class TestCartesianSphericalConversions:
    def setup_method(self, method):
        self.cart_sph = [[(0, 0, 1), (90, 0, 1)],
                         [(0, 1, 0), (0, 90, 1)],
                         [(1, 0, 0), (0, 0, 1)],
                         [(0, 0, -1), (-90, 0, 1)],
                         [(0, -1, 0), (0, -90, 1)],
                         [(0, -1, 0), (0, 270, 1)],
                         [(-1, 0, 0), (0, 180, 1)],
                         [(-1, 0, 0), (0, -180, 1)],
                         [(1, 0, 0), (0, 360, 1)],
                         [(3, 4, 0), (0, 53.13, 5)],
                        ]

    def positive_longitude(self, lon):
        if lon < 0:
            lon += 360
        while lon >= 360:
            lon -= 360
        return lon

    def equal_sph(self, coord1, coord2):
        coord1, coord2 = list(coord1), list(coord2)
        for coord in [coord1, coord2]:
            coord[1] = self.positive_longitude(coord[1])
        return np.allclose(coord1, coord)

    def test_cart2sph(self):
        for cart, sph in self.cart_sph:
            assert self.equal_sph(sph, euler_pole.cart2sph(*cart))

    def test_sph2cart(self):
        for cart, sph in self.cart_sph:
            assert np.allclose(cart, euler_pole.sph2cart(*sph))

    def test_chain_cart2sph(self):
        for cart, sph in self.cart_sph:
            trans_cart = euler_pole.sph2cart(*euler_pole.cart2sph(*cart))
            trans_sph = euler_pole.cart2sph(*euler_pole.sph2cart(*sph))
            assert np.allclose(cart, trans_cart)
            assert self.equal_sph(sph, trans_sph)

    def test_magnitude(self):
        x, y, z = euler_pole.sph2cart(0, 0, 1)
        assert np.allclose(np.sqrt(x**2 + y**2 + z**2), 1)

class TestLocalCoords:
    def test_simple(self):
        lat, lon = 0, 0
        x, y, z = 0, 1, 0
        assert np.allclose([0, 1, 0], euler_pole.local_coords(lat, lon, x,y,z))

    def test_zero_down(self):
        positions = itertools.product([45, -45], [45, -45, 135, 215])
        for lat, lon in positions:
            north_pole = [0, 0, 1]
            x, y, z = np.cross(north_pole, euler_pole.sph2cart(lat, lon))
            assert np.allclose(z, 0)
            assert np.allclose(np.hypot(x, y), np.sqrt(2) / 2)


class TestEulerPole:
    def test_velocity(self):
        pass
        

