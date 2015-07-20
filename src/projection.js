//conversion between geographic and web mercator
//theory: https://en.wikipedia.org/wiki/Web_Mercator

//conversion between UTM and geographic 
//theory: https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system#From_UTM_coordinates_.28E.2C_N.2C_Zone.2C_Hemi.29_to_latitude.2C_longitude_.28.CF.86.2C_.CE.BB.29
//tool: http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html

(function(geoUtil) {
  var pi = 3.14159265358979,
    /* Ellipsoid model constants (actual values here are for WGS84) */
    axisA = 6378137.0,
    axisB = 6356752.314,
    eccSquared = 6.69437999013e-03,
    utmScaleFactor = 0.9996;

  function footpointLatitude(y) {
    var y_, alpha_, beta_, gamma_, delta_, epsilon_, n;
    var result;
    /* Precalculate n (Eq. 10.18) */
    n = (axisA - axisB) / (axisA + axisB);
    /* Precalculate alpha_ (Eq. 10.22) */
    /* (Same as alpha in Eq. 10.17) */
    alpha_ = ((axisA + axisB) / 2.0) * (1 + (Math.pow(n, 2.0) / 4) + (Math.pow(n, 4.0) / 64));
    /* Precalculate y_ (Eq. 10.23) */
    y_ = y / alpha_;
    /* Precalculate beta_ (Eq. 10.22) */
    beta_ = (3.0 * n / 2.0) + (-27.0 * Math.pow(n, 3.0) / 32.0) + (269.0 * Math.pow(n, 5.0) / 512.0);
    /* Precalculate gamma_ (Eq. 10.22) */
    gamma_ = (21.0 * Math.pow(n, 2.0) / 16.0) + (-55.0 * Math.pow(n, 4.0) / 32.0);
    /* Precalculate delta_ (Eq. 10.22) */
    delta_ = (151.0 * Math.pow(n, 3.0) / 96.0) + (-417.0 * Math.pow(n, 5.0) / 128.0);
    /* Precalculate epsilon_ (Eq. 10.22) */
    epsilon_ = (1097.0 * Math.pow(n, 4.0) / 512.0);
    /* Now calculate the sum of the series (Eq. 10.21) */
    result = y_ + (beta_ * Math.sin(2.0 * y_)) + (gamma_ * Math.sin(4.0 * y_)) + (delta_ * Math.sin(6.0 * y_)) + (epsilon_ * Math.sin(8.0 * y_));

    return result;
  }

  function mapXYToLatLon(x, y, lambda0, philambda) {
    var phif, Nf, Nfpow, nuf2, ep2, tf, tf2, tf4, cf;
    var x1frac, x2frac, x3frac, x4frac, x5frac, x6frac, x7frac, x8frac;
    var x2poly, x3poly, x4poly, x5poly, x6poly, x7poly, x8poly;

    /* Get the value of phif, the footpoint latitude. */
    phif = footpointLatitude(y);
    /* Precalculate ep2 */
    ep2 = (Math.pow(axisA, 2.0) - Math.pow(axisB, 2.0)) / Math.pow(axisB, 2.0);
    /* Precalculate cos (phif) */
    cf = Math.cos(phif);
    /* Precalculate nuf2 */
    nuf2 = ep2 * Math.pow(cf, 2.0);
    /* Precalculate Nf and initialize Nfpow */
    Nf = Math.pow(axisA, 2.0) / (axisB * Math.sqrt(1 + nuf2));
    Nfpow = Nf;
    /* Precalculate tf */
    tf = Math.tan(phif);
    tf2 = tf * tf;
    tf4 = tf2 * tf2;
    /* Precalculate fractional coefficients for x**n in the equations
       below to simplify the expressions for latitude and longitude. */
    x1frac = 1.0 / (Nfpow * cf);
    Nfpow *= Nf; /* now equals Nf**2) */
    x2frac = tf / (2.0 * Nfpow);
    Nfpow *= Nf; /* now equals Nf**3) */
    x3frac = 1.0 / (6.0 * Nfpow * cf);
    Nfpow *= Nf; /* now equals Nf**4) */
    x4frac = tf / (24.0 * Nfpow);
    Nfpow *= Nf; /* now equals Nf**5) */
    x5frac = 1.0 / (120.0 * Nfpow * cf);
    Nfpow *= Nf; /* now equals Nf**6) */
    x6frac = tf / (720.0 * Nfpow);
    Nfpow *= Nf; /* now equals Nf**7) */
    x7frac = 1.0 / (5040.0 * Nfpow * cf);
    Nfpow *= Nf; /* now equals Nf**8) */
    x8frac = tf / (40320.0 * Nfpow);
    /* Precalculate polynomial coefficients for x**n.
       -- x**1 does not have a polynomial coefficient. */
    x2poly = -1.0 - nuf2;
    x3poly = -1.0 - 2 * tf2 - nuf2;
    x4poly = 5.0 + 3.0 * tf2 + 6.0 * nuf2 - 6.0 * tf2 * nuf2 - 3.0 * (nuf2 * nuf2) - 9.0 * tf2 * (nuf2 * nuf2);
    x5poly = 5.0 + 28.0 * tf2 + 24.0 * tf4 + 6.0 * nuf2 + 8.0 * tf2 * nuf2;
    x6poly = -61.0 - 90.0 * tf2 - 45.0 * tf4 - 107.0 * nuf2 + 162.0 * tf2 * nuf2;
    x7poly = -61.0 - 662.0 * tf2 - 1320.0 * tf4 - 720.0 * (tf4 * tf2);
    x8poly = 1385.0 + 3633.0 * tf2 + 4095.0 * tf4 + 1575 * (tf4 * tf2);
    /* Calculate latitude */
    philambda[0] = phif + x2frac * x2poly * (x * x) + x4frac * x4poly * Math.pow(x, 4.0) + x6frac * x6poly * Math.pow(x, 6.0) + x8frac * x8poly * Math.pow(x, 8.0);
    /* Calculate longitude */
    philambda[1] = lambda0 + x1frac * x + x3frac * x3poly * Math.pow(x, 3.0) + x5frac * x5poly * Math.pow(x, 5.0) + x7frac * x7poly * Math.pow(x, 7.0);

    return;
  }

  function radToDeg(rad) {
    return (rad / pi * 180.0)
  }

  function degToRad(deg) {
    return (deg / 180.0 * pi)
  }

  function validateLngLat(lng, lat){
    var valid = true;
    if (lng < -180 || lng > 180 || lat < -90 || lat > 90) {
      valid = false; 
    }
    return valid; 
  }

  function utmCentralMeridian(zone) {
    var cmeridian;
    cmeridian = degToRad(-183.0 + (zone * 6.0));
    return cmeridian;
  }

  function utmXYToLatLon(x, y, zone, southhemi, latlon) {
    var cmeridian;
    x -= 500000.0;
    x /= utmScaleFactor;
    /* If in southern hemisphere, adjust y accordingly. */
    if (southhemi)
      y -= 10000000.0;
    y /= utmScaleFactor;
    cmeridian = utmCentralMeridian(zone);
    mapXYToLatLon(x, y, cmeridian, latlon);

    return;
  }

  function utm2Lnglat(x, y, zone, southhemi) {
    var latlon = new Array(2);
    utmXYToLatLon(x, y, zone, southhemi, latlon);
    return {
      "longitude": radToDeg(latlon[1]),
      "latitude": radToDeg(latlon[0])
    };
  }

  function arcLengthOfMeridian (phi)
    {
        var alpha, beta, gamma, delta, epsilon, n;
        var result;

        /* Precalculate n */
        n = (axisA - axisB) / (axisA + axisB);

        /* Precalculate alpha */
        alpha = ((axisA + axisB) / 2.0)
           * (1.0 + (Math.pow (n, 2.0) / 4.0) + (Math.pow (n, 4.0) / 64.0));

        /* Precalculate beta */
        beta = (-3.0 * n / 2.0) + (9.0 * Math.pow (n, 3.0) / 16.0)
           + (-3.0 * Math.pow (n, 5.0) / 32.0);

        /* Precalculate gamma */
        gamma = (15.0 * Math.pow (n, 2.0) / 16.0)
            + (-15.0 * Math.pow (n, 4.0) / 32.0);
    
        /* Precalculate delta */
        delta = (-35.0 * Math.pow (n, 3.0) / 48.0)
            + (105.0 * Math.pow (n, 5.0) / 256.0);
    
        /* Precalculate epsilon */
        epsilon = (315.0 * Math.pow (n, 4.0) / 512.0);
    
    /* Now calculate the sum of the series and return */
    result = alpha
        * (phi + (beta * Math.sin (2.0 * phi))
            + (gamma * Math.sin (4.0 * phi))
            + (delta * Math.sin (6.0 * phi))
            + (epsilon * Math.sin (8.0 * phi)));

    return result;
    }

  function mapLatLonToXY (phi, lambda, lambda0, xy)
    {
        var N, nu2, ep2, t, t2, l;
        var l3coef, l4coef, l5coef, l6coef, l7coef, l8coef;
        var tmp;

        /* Precalculate ep2 */
        ep2 = (Math.pow (axisA, 2.0) - Math.pow (axisB, 2.0)) / Math.pow (axisB, 2.0);
    
        /* Precalculate nu2 */
        nu2 = ep2 * Math.pow (Math.cos (phi), 2.0);
    
        /* Precalculate N */
        N = Math.pow (axisA, 2.0) / (axisB * Math.sqrt (1 + nu2));
    
        /* Precalculate t */
        t = Math.tan (phi);
        t2 = t * t;
        tmp = (t2 * t2 * t2) - Math.pow (t, 6.0);

        /* Precalculate l */
        l = lambda - lambda0;
    
        /* Precalculate coefficients for l**n in the equations below
           so a normal human being can read the expressions for easting
           and northing
           -- l**1 and l**2 have coefficients of 1.0 */
        l3coef = 1.0 - t2 + nu2;
    
        l4coef = 5.0 - t2 + 9 * nu2 + 4.0 * (nu2 * nu2);
    
        l5coef = 5.0 - 18.0 * t2 + (t2 * t2) + 14.0 * nu2
            - 58.0 * t2 * nu2;
    
        l6coef = 61.0 - 58.0 * t2 + (t2 * t2) + 270.0 * nu2
            - 330.0 * t2 * nu2;
    
        l7coef = 61.0 - 479.0 * t2 + 179.0 * (t2 * t2) - (t2 * t2 * t2);
    
        l8coef = 1385.0 - 3111.0 * t2 + 543.0 * (t2 * t2) - (t2 * t2 * t2);
    
        /* Calculate easting (x) */
        xy[0] = N * Math.cos (phi) * l
            + (N / 6.0 * Math.pow (Math.cos (phi), 3.0) * l3coef * Math.pow (l, 3.0))
            + (N / 120.0 * Math.pow (Math.cos (phi), 5.0) * l5coef * Math.pow (l, 5.0))
            + (N / 5040.0 * Math.pow (Math.cos (phi), 7.0) * l7coef * Math.pow (l, 7.0));
    
        /* Calculate northing (y) */
        xy[1] = arcLengthOfMeridian (phi)
            + (t / 2.0 * N * Math.pow (Math.cos (phi), 2.0) * Math.pow (l, 2.0))
            + (t / 24.0 * N * Math.pow (Math.cos (phi), 4.0) * l4coef * Math.pow (l, 4.0))
            + (t / 720.0 * N * Math.pow (Math.cos (phi), 6.0) * l6coef * Math.pow (l, 6.0))
            + (t / 40320.0 * N * Math.pow (Math.cos (phi), 8.0) * l8coef * Math.pow (l, 8.0));
    
        return;
    }

  function lnglat2Utm(lng, lat)
    {
      if (!validateLngLat(lng, lat)) {
        return ;
      }
      var xy = new Array(2), zone = Math.floor ((lng + 180.0) / 6) + 1, south = lat < 0;
        mapLatLonToXY (degToRad(lat), degToRad(lng), utmCentralMeridian (zone), xy);

        /* Adjust easting and northing for UTM system. */
        xy[0] = xy[0] * utmScaleFactor + 500000.0;
        xy[1] = xy[1] * utmScaleFactor;
        if (xy[1] < 0.0)
            xy[1] = xy[1] + 10000000.0;

        return {
          "zone": zone,
          "south": south,
          "x": xy[0],
          "y": xy[1]
        };
    }

  //convert epsg to utm
  function getUtm(epsg) {
    var south = true,
      zone = epsg - 32700;
    if (zone < 0) {
      south = false;
      zone += 100;
    }
    return {
      "zone": zone,
      "south": south
    }
  }

  function getEpsg(/*number*/zone, /*bool*/south) {
    var epsg;
    if (south) {
      epsg = 32700 + zone;
    }
    else {
      epsg = 32600 + zone;
    }
    return epsg;
  }

  function lnglat2WebMercator(lng, lat) {
    if (!validateLngLat(lng, lat)) {
      return ;
    }
    var lat_rad = degToRad(lat);
    return [degToRad(lng) * axisA, axisA / 2.0 * Math.log((1.0 + Math.sin(lat_rad)) / (1.0 - Math.sin(lat_rad)))];
  }

  function webMercator2Lnglat(x, y, isLinear) {
    var lng_deg = radToDeg(x / axisA);

    if (isLinear) {
      return [lng_deg, radToDeg((PI / 2) - (2 * Math.atan(Math.exp(-1.0 * y / axisA))))];
    }
    return [lng_deg - (Math.floor((lng_deg + 180) / 360) * 360), radToDeg((pi / 2) - (2 * Math.atan(Math.exp(-1.0 * y / axisA))))];
  };

  /*************************************************/
  /* public mehods *********************************/
  /*************************************************/
  //convert between lat lon and web mercator
  geoUtil.lnglat2WebMercator = lnglat2WebMercator;
  geoUtil.webMercator2Lnglat = webMercator2Lnglat;
  //convert between UTM and lat lon
  geoUtil.lnglat2Utm = function(lng, lat){
    var utm = lnglat2Utm(lng, lat);
    utm.epsg = getEpsg(utm.zone, utm.south);
    return utm;
  };
  geoUtil.utm2Lnglat = function (x, y, epsg) {
    utmZone = getUtm(epsg);
    return utm2Lnglat(x, y, utmZone.zone, utmZone.south);
  };
})(window.geoUtil = window.geoUtil || {});