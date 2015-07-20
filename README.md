# geo-util
JavaScript utitlities for online mapping. It includes several useful tools

## [projection](https://github.com/J-Huang/geo-util/blob/master/src/projection.js)
- Conversion between UTM and Lat/Lng, and conversion between Lat/Lng and Web Mercator
UTM is probably the most used projection in GIS. In order to plot it on web maps, which are under Web Mercator, it needs to reproject it from UTM to Lat/Lng, then convert Lat/Lng to Web Mercator<br/>
### samples
[<img src="https://github.com/J-Huang/geo-util/blob/master/examples/utm_conversion.png" alt="UTM conversion" />](https://rawgit.com/J-Huang/geo-util/master/examples/projection.html)
### usage
- add script tag linking to the projection.js
- geoUtil.utm2Lnglat(x, y, epsg). It returns an object with latitude and longitude properties.
- geoUtil.lnglat2Utm(lng, lat). It returns an object, which has several properties
