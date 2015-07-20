(function(geoUtil) {
  var a = 6378137,
    b = 6356752.31424518,
    f = 1 / 298.257223563, // WGS84 ellipsoid params
    toRad = Math.PI / 180;

  //private functions
  function _toEqualAreaPoint(pt) {
    var eSq = 0.00669437999019741354678198566736,
      e = 0.08181919084296430236105472696748,
      sinY = Math.sin(pt[1] * toRad),
      q = (1 - eSq) * ((sinY / (1 - eSq * (sinY * sinY)) - (1 / (2 * e)) * Math.log((1 - e * sinY) / (1 + e * sinY)))),
      x = a * pt[0] * toRad,
      y = a * q * 0.5,
      equalAreaCynlindricalProjectedPt = [x, y];
    return equalAreaCynlindricalProjectedPt;
  };

  //direct and inverse geodesic solver algorithm is based on http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
  //http://en.wikipedia.org/wiki/Vincenty's_formulae
  function _directGeodeticSolver( /*radians*/ lat, /*radians*/ lon, /*radians*/ alpha, /*meters*/ s) {
    var sinA = Math.sin(alpha),
      cosA = Math.cos(alpha),
      tanU1 = (1 - f) * Math.tan(lat),
      cosU1 = 1 / Math.sqrt((1 + tanU1 * tanU1)),
      sinU1 = tanU1 * cosU1,
      sigma1 = Math.atan2(tanU1, cosA),
      sinsqAlpha = (cosU1 * sinA) * (cosU1 * sinA),
      cosSqAlpha = 1 - sinsqAlpha,
      uSq = cosSqAlpha * (a * a - b * b) / (b * b),
      coef1 = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq))),
      coef2 = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq))),
      sigma = s / (b * coef1),
      sigmaP = 2 * Math.PI,
      sinSigma, cosSigma, cos2SigmaM, deltaSigma;
    while (Math.abs(sigma - sigmaP) > 1e-12) {
      cos2SigmaM = Math.cos(2 * sigma1 + sigma);
      sinSigma = Math.sin(sigma);
      cosSigma = Math.cos(sigma);
      deltaSigma = coef2 * sinSigma * (cos2SigmaM + coef2 / 4 * (cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) - coef2 / 6 * cos2SigmaM * (-3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM)));
      sigmaP = sigma;
      sigma = s / (b * coef1) + deltaSigma;
    }
    var temp = sinU1 * sinSigma - cosU1 * cosSigma * cosA,
      lat2 = Math.atan2(sinU1 * cosSigma + cosU1 * sinSigma * cosA, (1 - f) * Math.sqrt(sinsqAlpha + temp * temp)),
      lambda = Math.atan2(sinSigma * sinA, cosU1 * cosSigma - sinU1 * sinSigma * cosA),
      coef3 = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha)),
      adjustedLon = lambda - (1 - coef3) * f * Math.sqrt(sinsqAlpha) * (sigma + coef3 * sinSigma * (cos2SigmaM + coef3 * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM))),
      lat2Deg = lat2 / (Math.PI / 180),
      lon2Deg = (lon + adjustedLon) / (Math.PI / 180),
      pt = [lon2Deg, lat2Deg];
    return pt;
  }

  function _inverseGeodeticSolver( /*radians*/ lat1, /*radians*/ lon1, /*radians*/ lat2, /*radians*/ lon2) {
    var d = (lon2 - lon1),
      u1 = Math.atan((1 - f) * Math.tan(lat1)),
      u2 = Math.atan((1 - f) * Math.tan(lat2)),
      sinU1 = Math.sin(u1),
      cosU1 = Math.cos(u1),
      sinU2 = Math.sin(u2),
      cosU2 = Math.cos(u2),
      lambda = d,
      lambdaP, iterLimit = 1000,
      cosSqAlpha, sinSigma, cos2SigmaM, cosSigma, sigma, sinLambda, cosLambda, sinAlpha, temp;
    do {
      sinLambda = Math.sin(lambda);
      cosLambda = Math.cos(lambda);
      sinSigma = Math.sqrt((cosU2 * sinLambda) * (cosU2 * sinLambda) + (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) * (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda));
      if (sinSigma === 0) {
        return 0;
      }
      cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;
      sigma = Math.atan2(sinSigma, cosSigma);
      sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
      cosSqAlpha = 1 - sinAlpha * sinAlpha;
      cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;
      if (isNaN(cos2SigmaM)) {
        cos2SigmaM = 0;
      }
      temp = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));
      lambdaP = lambda;
      lambda = d + (1 - temp) * f * sinAlpha * (sigma + temp * sinSigma * (cos2SigmaM + temp * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM)));
    }
    while (Math.abs(lambda - lambdaP) > 1e-12 && --iterLimit > 0);
    if (iterLimit === 0) {
      //As Vincenty pointed out, when two points are nearly antipodal, the formula may not converge
      //It's time to switch to other formula, which may not as highly accurate as Vincenty's. Just for the special case.
      //Here implements Haversine formula
      var haversine_R = 6371009, // km
        haversine_d = Math.acos(Math.sin(lat1) * Math.sin(lat2) + Math.cos(lat1) * Math.cos(lat2) * Math.cos(lon2 - lon1)) * haversine_R,
        dLon = lon2 - lon1,
        y = Math.sin(dLon) * Math.cos(lat2),
        x = Math.cos(lat1) * Math.sin(lat2) - Math.sin(lat1) * Math.cos(lat2) * Math.cos(dLon),
        brng = Math.atan2(y, x);
      return {
        "azimuth": brng,
        "geodesicDistance": haversine_d
      };
    }
    var uSq = cosSqAlpha * (a * a - b * b) / (b * b),
      coef1 = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq))),
      coef2 = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq))),
      deltaSigma = coef2 * sinSigma * (cos2SigmaM + coef2 / 4 * (cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) - coef2 / 6 * cos2SigmaM * (-3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM))),
      geodesicDistance = b * coef1 * (sigma - deltaSigma),
      azimuth = Math.atan2(cosU2 * Math.sin(lambda), cosU1 * sinU2 - sinU1 * cosU2 * Math.cos(lambda)),
      reverseAzimuth = Math.atan2(cosU1 * Math.sin(lambda), cosU1 * sinU2 * Math.cos(lambda) - sinU1 * cosU2),
      inverseResult = {
        azimuth: azimuth,
        geodesicDistance: geodesicDistance,
        reverseAzimuth: reverseAzimuth
      };
    return inverseResult;
  }

  function geodesicDensify(pts, maxSegmentLength) {
    var radius = 6371008.771515059,
      result = [];
    if (!maxSegmentLength) {
      maxSegmentLength = 10000;
    }
    if (maxSegmentLength < radius / 10000) {
      maxSegmentLength = radius / 10000;
    }

    var i, j, result = [],
      lon1, lat1, lon2, lat2, i, j;
    result.push([pts[0][0], pts[0][1]]);
    lon1 = pts[0][0] * toRad;
    lat1 = pts[0][1] * toRad;
    for (i = 0; i < pts.length - 1; i++) {
      lon2 = pts[i + 1][0] * toRad;
      lat2 = pts[i + 1][1] * toRad;
      if (lon1 === lon2 && lat1 === lat2) {
        continue;
      }
      var inverseGeodeticResult = _inverseGeodeticSolver(lat1, lon1, lat2, lon2);
      var azimuth = inverseGeodeticResult.azimuth; //radians
      var geodesicDist = inverseGeodeticResult.geodesicDistance; //meters
      var numberOfSegment = geodesicDist / maxSegmentLength;
      if (numberOfSegment > 1) {
        for (j = 1; j <= numberOfSegment - 1; j++) {
          var length = j * maxSegmentLength;
          var pt = _directGeodeticSolver(lat1, lon1, azimuth, length);
          result.push([pt[0], pt[1]]);
        }
        var lastDensifiedLength = (geodesicDist + Math.floor(numberOfSegment - 1) * maxSegmentLength) / 2;
        var lastSecondPt = _directGeodeticSolver(lat1, lon1, azimuth, lastDensifiedLength);
        result.push([lastSecondPt[0], lastSecondPt[1]]);
      }
      var endPt = _directGeodeticSolver(lat1, lon1, azimuth, geodesicDist);
      result.push([endPt[0], endPt[1]]);
      lon1 = endPt[0] * toRad;
      lat1 = endPt[1] * toRad;
    }

    return result;
  };

  var unitConversion = {
    //length unit conversion from miles
    "miles": 1,
    "kilometers": 1.609344,
    "feet": 5280,
    "meters": 1609.34,
    "yards": 1760,
    "nauticalMiles": 0.869,
    "centimeters": 160934,
    "decimeters": 16093.4,
    "inches": 63360,
    "millimeters": 1609340,
    //area unit conversion from acres
    "acres": 1,
    "ares": 40.4685642,
    "squareKilometers": 0.00404685642,
    "squareMiles": 0.0015625,
    "squareFeet": 43560,
    "squareMeters": 4046.85642,
    "hectares": 0.404685642,
    "squareYards": 4840,
    "squareInches": 6272640,
    "squareMillimeters": 4046856420,
    "squareCentimeters": 40468564.2,
    "squareDecimeters": 404685.642
  };

  geoUtil.geodesicDensify = geodesicDensify;

  //geodesic distance
  geoUtil.computeDistance = function(pt1, pt2, unit) {
    //pt1 and pt2 are array with lng and lat value
    var distance, lng1 = pt1[0] * toRad,
      lng2 = pt2[0] * toRad,
      lat1 = pt1[1] * toRad,
      lat2 = pt2[1] * toRad;

    if (!(lat1 === lat2 && lng1 === lng2)) {
      inverseGeodeticResult = _inverseGeodeticSolver(lat1, lng1, lat2, lng2);
      distance = inverseGeodeticResult.geodesicDistance / 1609.344; //miles
    }

    if (unit) {
      distance *= unitConversion[unit];
    }

    return distance;
  };

  geoUtil.computeArea = function(pts, unit) {
    //pts is array of pt, which is an array with 2 values
    if (!pts.length || pts.length < 3) {
      return null;
    }
    var ring = geodesicDensify(pts, 10000),
      point1 = _toEqualAreaPoint([ring[0][0], ring[0][1]]),
      point2 = _toEqualAreaPoint([ring[ring.length - 1][0], ring[ring.length - 1][1]]),
      area = point2[0] * point1[1] - point1[0] * point2[1],
      i;

    for (i = 0; i < ring.length - 1; i++) {
      point1 = _toEqualAreaPoint([ring[i + 1][0], ring[i + 1][1]]);
      point2 = _toEqualAreaPoint([ring[i][0], ring[i][1]]);
      area += point2[0] * point1[1] - point1[0] * point2[1];
    }
    area /= 4046.87; //acres
    if (unit) {
      area *= unitConversion[unit];
    }
    area /= (-2);

    return Math.abs(area);
  };

  geoUtil.createBuffer = function(pts, radius) {
    //create buffer around point, polyline or polygon

  };
})(window.geoUtil = window.geoUtil || {});