var gulp = require('gulp'),
  uglify = require('gulp-uglify');

gulp.task('default', function() {
  return gulp.src('src/*.js')
    .pipe(uglify())
    .pipe(gulp.dest('dist'));
});