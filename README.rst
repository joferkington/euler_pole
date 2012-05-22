euler_pole
============
`euler_pole` provides an `EulerPole` class for plate motion calculations using
Euler poles.  This allows the calculation of velocity at specified points on
the Earth's surface as well as calculating the relative motion between two
defined Euler poles.

Basic Usage
-----------
The `EulerPole` class allows addition and subtraction of euler poles and the calculation of velocity


As an example::

        from euler_pole import EulerPole

        # Here we have two euler poles from Loveless & Meade, 2010 given
        # relative to the stable Eurasian reference frame of Apel, 2006.
        # We want to calculate the relative motion between them to determine
        # the motion of the Philippine Sea Plate relative to the Amur Plate.

        # Movement of the Amur plate relative to the Eurasian reference frame
        #                degrees    degrees     deg / myr
        amur = EulerPole(lat=66.80, lon=143.93, rate=0.26)

        # Movement of the Philippine Sea Plate relative to the Eurasian frame
        philippine = EulerPole(-44.85, 336.45, 1.07)

        # Movement of the Philippine Sea Plate relative to the Amur Plate
        phil_amur = philippine - amur

        # Velocity at a given point:
        print phil_amur.velocity(33.273, 136.781)


