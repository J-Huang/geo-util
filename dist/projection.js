!function(t){function a(t){var a,o,h,M,n,w,r,p;return r=(g-m)/(g+m),o=(g+m)/2*(1+Math.pow(r,2)/4+Math.pow(r,4)/64),a=t/o,h=3*r/2+-27*Math.pow(r,3)/32+269*Math.pow(r,5)/512,M=21*Math.pow(r,2)/16+-55*Math.pow(r,4)/32,n=151*Math.pow(r,3)/96+-417*Math.pow(r,5)/128,w=1097*Math.pow(r,4)/512,p=a+h*Math.sin(2*a)+M*Math.sin(4*a)+n*Math.sin(6*a)+w*Math.sin(8*a)}function o(t,o,h,M){var n,w,r,p,u,e,i,c,s,f,v,l,d,z,x,y,U,Z,b,q,A,L,I,P;n=a(o),u=(Math.pow(g,2)-Math.pow(m,2))/Math.pow(m,2),s=Math.cos(n),p=u*Math.pow(s,2),w=Math.pow(g,2)/(m*Math.sqrt(1+p)),r=w,e=Math.tan(n),i=e*e,c=i*i,f=1/(r*s),r*=w,v=e/(2*r),r*=w,l=1/(6*r*s),r*=w,d=e/(24*r),r*=w,z=1/(120*r*s),r*=w,x=e/(720*r),r*=w,y=1/(5040*r*s),r*=w,U=e/(40320*r),Z=-1-p,b=-1-2*i-p,q=5+3*i+6*p-6*i*p-3*p*p-9*i*p*p,A=5+28*i+24*c+6*p+8*i*p,L=-61-90*i-45*c-107*p+162*i*p,I=-61-662*i-1320*c-720*c*i,P=1385+3633*i+4095*c+1575*c*i,M[0]=n+v*Z*t*t+d*q*Math.pow(t,4)+x*L*Math.pow(t,6)+U*P*Math.pow(t,8),M[1]=h+f*t+l*b*Math.pow(t,3)+z*A*Math.pow(t,5)+y*I*Math.pow(t,7)}function h(t){return t/l*180}function M(t){return t/180*l}function n(t,a){var o=!0;return(-180>t||t>180||-90>a||a>90)&&(o=!1),o}function w(t){var a;return a=M(-183+6*t)}function r(t,a,h,M,n){var r;t-=5e5,t/=d,M&&(a-=1e7),a/=d,r=w(h),o(t,a,r,n)}function p(t,a,o,M){var n=new Array(2);return r(t,a,o,M,n),{longitude:h(n[1]),latitude:h(n[0])}}function u(t){var a,o,h,M,n,w,r;return w=(g-m)/(g+m),a=(g+m)/2*(1+Math.pow(w,2)/4+Math.pow(w,4)/64),o=-3*w/2+9*Math.pow(w,3)/16+-3*Math.pow(w,5)/32,h=15*Math.pow(w,2)/16+-15*Math.pow(w,4)/32,M=-35*Math.pow(w,3)/48+105*Math.pow(w,5)/256,n=315*Math.pow(w,4)/512,r=a*(t+o*Math.sin(2*t)+h*Math.sin(4*t)+M*Math.sin(6*t)+n*Math.sin(8*t))}function e(t,a,o,h){var M,n,w,r,p,e,i,c,s,f,v,l,d;w=(Math.pow(g,2)-Math.pow(m,2))/Math.pow(m,2),n=w*Math.pow(Math.cos(t),2),M=Math.pow(g,2)/(m*Math.sqrt(1+n)),r=Math.tan(t),p=r*r,d=p*p*p-Math.pow(r,6),e=a-o,i=1-p+n,c=5-p+9*n+4*n*n,s=5-18*p+p*p+14*n-58*p*n,f=61-58*p+p*p+270*n-330*p*n,v=61-479*p+179*p*p-p*p*p,l=1385-3111*p+543*p*p-p*p*p,h[0]=M*Math.cos(t)*e+M/6*Math.pow(Math.cos(t),3)*i*Math.pow(e,3)+M/120*Math.pow(Math.cos(t),5)*s*Math.pow(e,5)+M/5040*Math.pow(Math.cos(t),7)*v*Math.pow(e,7),h[1]=u(t)+r/2*M*Math.pow(Math.cos(t),2)*Math.pow(e,2)+r/24*M*Math.pow(Math.cos(t),4)*c*Math.pow(e,4)+r/720*M*Math.pow(Math.cos(t),6)*f*Math.pow(e,6)+r/40320*M*Math.pow(Math.cos(t),8)*l*Math.pow(e,8)}function i(t,a){if(n(t,a)){var o=new Array(2),h=Math.floor((t+180)/6)+1,r=0>a;return e(M(a),M(t),w(h),o),o[0]=o[0]*d+5e5,o[1]=o[1]*d,o[1]<0&&(o[1]=o[1]+1e7),{zone:h,south:r,x:o[0],y:o[1]}}}function c(t){var a=!0,o=t-32700;return 0>o&&(a=!1,o+=100),{zone:o,south:a}}function s(t,a){var o;return o=a?32700+t:32600+t}function f(t,a){if(n(t,a)){var o=M(a);return[M(t)*g,g/2*Math.log((1+Math.sin(o))/(1-Math.sin(o)))]}}function v(t,a,o){var M=h(t/g);return o?[M,h(PI/2-2*Math.atan(Math.exp(-1*a/g)))]:[M-360*Math.floor((M+180)/360),h(l/2-2*Math.atan(Math.exp(-1*a/g)))]}var l=3.14159265358979,g=6378137,m=6356752.314,d=.9996;t.lnglat2WebMercator=f,t.webMercator2Lnglat=v,t.lnglat2Utm=function(t,a){var o=i(t,a);return o.epsg=s(o.zone,o.south),o},t.utm2Lnglat=function(t,a,o){return utmZone=c(o),p(t,a,utmZone.zone,utmZone.south)}}(window.geoUtil=window.geoUtil||{});