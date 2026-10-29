import numpy as np
def latlon_to_xy(lat, lon, proj):
    """GOES-R PUG 5.1.2.8.1: geodetic lat/lon -> fixed-grid scan angles (rad)."""
    req = proj['semi_major_axis']; rpol = proj['semi_minor_axis']
    H = proj['perspective_point_height'] + req
    lon0 = np.deg2rad(proj['longitude_of_projection_origin'])
    e = np.sqrt((req**2 - rpol**2) / req**2)
    phi = np.deg2rad(lat); lam = np.deg2rad(lon)
    phic = np.arctan((rpol**2 / req**2) * np.tan(phi))
    rc = rpol / np.sqrt(1 - e**2 * np.cos(phic)**2)
    sx = H - rc * np.cos(phic) * np.cos(lam - lon0)
    sy = -rc * np.cos(phic) * np.sin(lam - lon0)
    sz = rc * np.sin(phic)
    y = np.arctan(sz / sx)
    x = np.arcsin(-sy / np.sqrt(sx**2 + sy**2 + sz**2))
    return x, y
def xy_to_latlon(x, y, proj):
    req = proj['semi_major_axis']; rpol = proj['semi_minor_axis']
    H = proj['perspective_point_height'] + req
    lon0 = np.deg2rad(proj['longitude_of_projection_origin'])
    x, y = np.meshgrid(x, y)
    a = np.sin(x)**2 + np.cos(x)**2 * (np.cos(y)**2 + (req**2 / rpol**2) * np.sin(y)**2)
    b = -2 * H * np.cos(x) * np.cos(y)
    c = H**2 - req**2
    with np.errstate(invalid='ignore'):
        rs = (-b - np.sqrt(b**2 - 4 * a * c)) / (2 * a)
    sx = rs * np.cos(x) * np.cos(y); sy = -rs * np.sin(x); sz = rs * np.cos(x) * np.sin(y)
    lat = np.rad2deg(np.arctan((req**2 / rpol**2) * sz / np.sqrt((H - sx)**2 + sy**2)))
    lon = np.rad2deg(lon0 - np.arctan(sy / (H - sx)))
    return lat, lon
