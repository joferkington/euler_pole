import numpy as np

class EulerPole(object):
    earth_radius = 6371.009 # Km (Mean radius from IUGG)

    def __init__(self, lat, lon, rate):
        """
        Parameters:
        -----------
            lat: The latitude of the pole in degrees
            lon: The longitude of the pole in degrees
            rate: The rotation velocity of the pole in deg/myr
                          (counterclockwise-positive)
        """
        self.lat = lat
        self.lon = lon
        self.rot_velocity = -rate

    def __sub__(self, other):
        lat, lon, vel = cart2sph(*(self.omega - other.omega))
        return EulerPole(lat, lon, np.degrees(vel))

    def __add__(self, other):
        lat, lon, vel = cart2sph(*(self.omega + other.omega))
        return EulerPole(lat, lon, np.degrees(vel))

    def __neg__(self):
        lat, lon, vel = cart2sph(*-self.omega)
        return EulerPole(lat, lon, np.degrees(vel))

    def __repr__(self):
        template = 'EulerPole(lat={0}, lon={1}, rot_velocity={2})'
        return template.format(self.lat, self.lon, -self.rot_velocity)

    @property
    def coord_basis(self):
        """The coordinate basis for the pole's reference frame. 
        (a 3x3 numpy array)"""
        x, y, z = sph2cart(self.lat, self.lon)
        pole = np.array([x,y,z])
        vec1 = np.cross(pole, [0, 0, 1])
        vec2 = np.cross(pole, vec1)

        basis = np.vstack([pole, vec1, vec2]).T
        return basis

    @property
    def pole_transform(self):
        """Transformation matrix to go between earth reference coordinates
        and the euler pole reference coordinates (a 3x3 numpy array)"""
        return np.linalg.inv(self.coord_basis)

    @property
    def inv_pole_transform(self):
        """Transformation matrix to go between the euler pole's reference
        coordinates and the earth's reference coordinates (a 3x3 numpy
        array)"""
        return self.coord_basis

    def rotate(self, lat, lon, angle):
        """Rotates a feature (given by arrays of "lat" and "lon") around 
        the Euler pole by the given angle. Returns 2 arrays of lat and lon
        with the same length as the input arrays."""
        # Convert the input into Euler pole basis coordinates (with the 
        # "north pole" at the euler pole)
        x, y, z = sph2cart(lat, lon)
        inp_xyz = np.vstack([x,y,z]).T
        xx, yy, zz = np.dot(inp_xyz, self.pole_transform).T
        lat_prime, lon_prime, r_prime = cart2sph(xx, yy, zz)

        # Add the rotation angle to the longitude...
        lon_prime += angle
        xx, yy, zz = sph2cart(lat_prime, lon_prime, r_prime)

        # ...And convert back into lat, long.
        xyz = np.vstack([xx,yy,zz]).T
        xx, yy, zz = np.dot(xyz, self.inv_pole_transform).T
        new_lat, new_lon, _ = cart2sph(xx, yy, zz)

        return new_lat, new_lon
        
    def move(self, lat, lon, time):
        """Moves a feature _back_ in time by "time" million years.
        Use a negative time to move the feature "into the future"."""
        angle = time * self.rot_velocity
        new_lat, new_lon = self.rotate(lat, lon, angle)
        return new_lat, new_lon

    @property
    def omega(self):
        """The Euler vector for the pole in geocentric Cartesian Coordinates."""
        vec = sph2cart(self.lat, self.lon, np.radians(self.rot_velocity))
        return np.array(vec)

    def velocity(self, lat, lon):
        """
        Calculates the azimuth (in degrees) and rate of plate motion 
        (in millimeters per year) at a given point.

        Parameters:
        -----------
            lat : The latitude of the point in degrees
            lon : The longitude of the point in degrees

        Returns:
            azimuth : The azimuth in degrees clockwise from north
            rate : The rate in mm/yr

        """
        east, north, down = self.velocity_components(lat, lon)
        azi = azimuth(east, north)
        rate = np.sqrt(north**2 + east**2 + down**2)
        return azi, rate

    def velocity_components(self, lat, lon):
        """
        Calculates the eastward, northward, and downward componenents (in mm/yr) 
        of plate velocity at a given point.

        Parameters:
        -----------
            lat : The latitude of the point in degrees
            lon : The longitude of the point in degrees

        Returns:
        --------
            east, north : The eastward and northward components of the plate
                velocity in millimeters per year at the given point.
        """
        # Convert position (from lat, lon) into geocentric cartesian coords
        r = sph2cart(lat, lon, self.earth_radius)

        # Velocity at the earth's surface is then just omega x r
        v = np.cross(self.omega, r)

        # We can then convert this back to local (north, east, down) coordinates
        east, north, down = local_coords(lat, lon, v[0], v[1], v[2])
        return east, north, down


#-- Utility Functions --------------------------------------------------------
def sph2cart(lat, lon, r=1):
    """Convert spherical coordinates to cartesian.  Default raduis is 1 (unit
    length).  Input is in degrees, output is in km."""
    lat, lon = np.radians(lat), np.radians(lon)
    x = r * np.cos(lat) * np.cos(lon)
    y = r * np.cos(lat) * np.sin(lon)
    z = r * np.sin(lat)
    return x,y,z

def cart2sph(x,y,z):
    """Convert cartesian geocentric coordinates to spherical coordinates.
    Output is in degrees (for lat & lon) and whatever input units are for the
    radius.""" 
    r = np.sqrt(x**2 + y**2 + z**2)
    lat = np.arcsin(z / r)
    lon = np.arctan2(y, x)
    return np.degrees(lat), np.degrees(lon), r

def local_coords(lat, lon, x,y,z):
    """Calculate local east,north,down components of x,y,z at lat,lon"""
    lat, lon = np.radians(lat), np.radians(lon)

    north = - np.sin(lat) * np.cos(lon) * x \
            - np.sin(lat) * np.sin(lon) * y \
            + np.cos(lat) * z

    east = - np.sin(lon) * x + np.cos(lon) * y

    down = - np.cos(lat) * np.cos(lon) * x \
           - np.cos(lat) * np.sin(lon) \
           - np.sin(lat) * z
    return east, north, down

def azimuth(east, north):
    """Returns azimuth in degrees counterclockwise from North given north and
    east components"""
    azi = np.degrees(np.arctan2(north, east))
    azi = 90 - azi
    if azi <= 0: 
        azi +=360
    return azi


