Epsilon = 1e-12
Ellipsoid_WGS_84 = {"a": 6378137.0, "b": 6356752.314245179}
Ellipsoid_PZ90_11 = {"a": 6378136, "reciprocalAlpha": 298.25784}
Ellipsoid_Krasovsky = {"a": 6378245, "b": 6356863.0188, "reciprocalAlpha": 298.3}
Ellipsoid = Ellipsoid_WGS_84

if "b" not in Ellipsoid:
    Ellipsoid["b"] = Ellipsoid["a"] * (1 -  1 / Ellipsoid["reciprocalAlpha"])


