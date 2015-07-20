# geo-util
JavaScript utitlities for online mapping. It includes several useful tools

## [projection](https://github.com/J-Huang/geo-util/blob/master/src/projection.js)
- Conversion between UTM and Lat/Lng, and conversion between Lat/Lng and Web Mercator
UTM is probably the most used projection in GIS. In order to plot it on web maps, which are under Web Mercator, it needs to reproject it from UTM to Lat/Lng, then convert Lat/Lng to Web Mercator

### samples
[<img src="https://github.com/J-Huang/geo-util/blob/master/examples/utm_conversion.png" style="width:300px;" alt="UTM conversion" />](https://rawgit.com/J-Huang/geo-util/master/examples/projection.html)

### usage
- geoUtil.utm2Lnglat(x, y, epsg). It returns an object with latitude and longitude properties.
- geoUtil.lnglat2Utm(lng, lat). It returns an object, which has several properties, including "zone", "south", "x", "y" and "epsg".
- geoUtil.webMercator2Lnglat(x, y). It returns an array with longitude and latitude value.
- geoUtil.lnglat2WebMercator(lng, lat). It returns an array with x and y values in web mercator.

## [geodesic](https://github.com/J-Huang/geo-util/blob/master/src/geodesic.js)
- geodesic util tools can densify between two lng/lat value pairs. It provides methods to measure distance and area.

### samples
- geodesic densify and distance measurement
[<img src="https://github.com/J-Huang/geo-util/blob/master/examples/geodesic_densify.png" style="width:300px;" alt="geodesic distance" />](https://rawgit.com/J-Huang/geo-util/master/examples/geodesic_distance_densify.html)

- area measurement
[<img src="https://github.com/J-Huang/geo-util/blob/master/examples/geodesic_area.png" style="width:300px;" alt="geodesic area" />](https://rawgit.com/J-Huang/geo-util/master/examples/geodesic_area.html)
