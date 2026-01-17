from typing import Optional, Tuple
import warnings

from scipy.spatial.transform import Rotation

OBLIQUITY = 23.4392911

class ElementsConverter:
    """Converts orbital elements to Celestia convention."""

    def __init__(self, ra: float, dec: float) -> None:
        """Creates a converter for given RA and Dec (both in degrees)."""
        self._convert = Rotation.from_euler('yzx', [-dec-90, ra, -OBLIQUITY], degrees=True)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')  # ignore gimbal lock warning
            self._arg_peri, self._inclination, self._node = self._convert.inv().as_euler('zxz', degrees=True)
        self._inclination = -self._inclination
        self._node += 180
        if self._node >= 360:
            self._node -= 360

    def convert(
        self,
        arg_peri: Optional[float] = None,
        inclination: Optional[float] = None,
        node: Optional[float] = None,
    ) -> Tuple[float, float, float]:
        """Converts orbital elements to Celestia convention. Angles in degrees."""
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')  # ignore gimbal lock warning
            orbit = Rotation.from_euler(
                'zxz',
                [
                    arg_peri if arg_peri is not None else self._arg_peri,
                    inclination if inclination is not None else self._inclination,
                    node if node is not None else self._node,
                ],
                degrees=True,
            )

            arg_peri, inclination, node = (self._convert * orbit).as_euler('zxz', degrees=True)

        if arg_peri < 0:
            arg_peri += 360
        if node < 0:
            node += 360

        return arg_peri, inclination, node


if __name__ == '__main__':
    print('70 Vir b')
    ec = ElementsConverter(
        (13 + 28/60 + 25.80819/3600) * 15,
        (13 + 46/60 + 43.6430/3600),
    )

    arg_peri, inclination, node = ec.convert(arg_peri=358.71-180.0)  # exoplanets use omega_1 so apply 180 degree correction
    print(f'arg_peri    = {arg_peri}')
    print(f'inclination = {inclination}')
    print(f'node        = {node}')

    print('\nHD 80606 b')
    ec = ElementsConverter(
        (9 + 22/60 + 37.5769/3600) * 15,
        (50 + 36.60/60 + 13.430/3600),
    )

    arg_peri, inclination, node = ec.convert(
        arg_peri=300.4977-180.0,  # exoplanets use omega_1 so apply 180 degree correction
        inclination=89.285,
    )
    print(f'arg_peri    = {arg_peri}')
    print(f'inclination = {inclination}')
    print(f'node        = {node}')

    print('\nalf Cen B')
    ec = ElementsConverter(
        (14 + 39/60 + 39.49400/3600) * 15,
        (-60 - 50/60 - 2.3737/3600),
    )
    arg_peri, inclination, node = ec.convert(
        arg_peri=231.65,  # no correction needed for visual binary
        inclination=79.205,
        node=204.85
    )
    print(f'arg_peri    = {arg_peri}')
    print(f'inclination = {inclination}')
    print(f'node        = {node}')
