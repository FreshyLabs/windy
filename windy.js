

/* Global class for simulating wind 
// params
    url to wind data .json
    canvas container el Id 

*/ 

var Windy = function( params ){
  var VELOCITY_SCALE = 1/60000;             // scale for wind velocity (completely arbitrary--this value looks nice)
  var OVERLAY_ALPHA = Math.floor(0.4*255);  // overlay transparency (on scale [0, 255])
  var INTENSITY_SCALE_STEP = 10;            // step size of particle intensity color scale
  var MAX_WIND_INTENSITY = 17;              // wind velocity at which particle intensity is maximum (m/s)
  var MAX_PARTICLE_AGE = 100;               // max number of frames a particle is drawn before regeneration
  var PARTICLE_LINE_WIDTH = .5;            // line width of a drawn particle
  var PARTICLE_MULTIPLIER = 7;              // particle count scalar (completely arbitrary--this values looks nice)
  var PARTICLE_REDUCTION = 0.75;            // reduce particle count to this much of normal for mobile devices
  var FRAME_RATE = 30;                      // desired milliseconds per frame

  var NULL_WIND_VECTOR = [NaN, NaN, null];  // singleton for no wind in the form: [u, v, magnitude]
  var TRANSPARENT_BLACK = [255, 0, 0, 0];



    var bilinearInterpolateScalar = function(x, y, g00, g10, g01, g11) {
        var rx = (1 - x);
        var ry = (1 - y);
        return g00 * rx * ry + g10 * x * ry + g01 * rx * y + g11 * x * y;
    };

    var bilinearInterpolateVector = function(x, y, g00, g10, g01, g11) {
        var rx = (1 - x);
        var ry = (1 - y);
        var a = rx * ry,  b = x * ry,  c = rx * y,  d = x * y;
        var u = g00[0] * a + g10[0] * b + g01[0] * c + g11[0] * d;
        var v = g00[1] * a + g10[1] * b + g01[1] * c + g11[1] * d;
        return [u, v, Math.sqrt(u * u + v * v)];
    };

    var createScalarBuilder = function(record) {
        var data = record.data;
        return {
            header: record.header,
            //recipe: recipeFor(""),
            data: function(i) {
                return data[i];
            },
            interpolate: bilinearInterpolateScalar
        }
    };

    var createWindBuilder = function(uComp, vComp) {
        var uData = uComp.data, vData = vComp.data;
        return {
            header: uComp.header,
            //recipe: recipeFor("wind-" + uComp.header.surface1Value),
            data: function(i) {
                return [uData[i], vData[i]];
            },
            interpolate: bilinearInterpolateVector
        }
    };

    var createBuilder = function(data) {
        var uComp = null, vComp = null, scalar = null;

        data.forEach(function(record) {
            switch (record.header.parameterCategory + "," + record.header.parameterNumber) {
                case "2,2": uComp = record; break;
                case "2,3": vComp = record; break;
                default:
                    scalar = record;
            }
        });

        return uComp ? createWindBuilder(uComp, vComp) : createScalarBuilder(scalar);
    };

    var buildGrid = function(data, callback) {
        var builder = createBuilder(data);

        var header = builder.header;
        var λ0 = header.lo1, φ0 = header.la1;  // the grid's origin (e.g., 0.0E, 90.0N)
        var Δλ = header.dx, Δφ = header.dy;    // distance between grid points (e.g., 2.5 deg lon, 2.5 deg lat)
        var ni = header.nx, nj = header.ny;    // number of grid points W-E and N-S (e.g., 144 x 73)
        var date = new Date(header.refTime);
        date.setHours(date.getHours() + header.forecastTime);

        // Scan mode 0 assumed. Longitude increases from λ0, and latitude decreases from φ0.
        // http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table3-4.shtml
        var grid = [], p = 0;
        var isContinuous = Math.floor(ni * Δλ) >= 360;
        for (var j = 0; j < nj; j++) {
            var row = [];
            for (var i = 0; i < ni; i++, p++) {
                row[i] = builder.data(p);
            }
            if (isContinuous) {
                // For wrapped grids, duplicate first column as last column to simplify interpolation logic
                row.push(row[0]);
            }
            grid[j] = row;
        }

        function interpolate(λ, φ) {
            var i = floorMod(λ - λ0, 360) / Δλ;  // calculate longitude index in wrapped range [0, 360)
            var j = (φ0 - φ) / Δφ;                 // calculate latitude index in direction +90 to -90

            var fi = Math.floor(i), ci = fi + 1;
            var fj = Math.floor(j), cj = fj + 1;


            var row;
            if ((row = grid[fj])) {
                var g00 = row[fi];
                var g10 = row[ci];
                if (isValue(g00) && isValue(g10) && (row = grid[cj])) {
                    var g01 = row[fi];
                    var g11 = row[ci];
                    if (isValue(g01) && isValue(g11)) {
                        // All four points found, so interpolate the value.
                        return builder.interpolate(i - fi, j - fj, g00, g10, g01, g11);
                    }
                }
            }
            return null;
        }

        callback( {
            date: date,
            //recipe: builder.recipe,
            interpolate: interpolate
        });
    };


  
    /**
     * @returns {Boolean} true if the specified value is not null and not undefined.
     */
    function isValue(x) {
        return x !== null && x !== undefined;
    }


  /**
     * @returns {Object} the first argument if not null and not undefined, otherwise the second argument.
     */
    function coalesce(a, b) {
        return isValue(a) ? a : b;
    }

    /**
     * @returns {Number} returns remainder of floored division, i.e., floor(a / n). Useful for consistent modulo
     *          of negative numbers. See http://en.wikipedia.org/wiki/Modulo_operation.
     */
    function floorMod(a, n) {
        return a - n * Math.floor(a / n);
    }

    /**
     * @returns {Number} distance between two points having the form [x, y].
     */
    function distance(a, b) {
        var Δx = b[0] - a[0];
        var Δy = b[1] - a[1];
        return Math.sqrt(Δx * Δx + Δy * Δy);
    }

    /**
     * @returns {Number} the value x clamped to the range [low, high].
     */
    function clamp(x, range) {
        return Math.max(range[0], Math.min(x, range[1]));
    }

    /**
     * Pad number with leading zeros. Does not support fractional or negative numbers.
     */
    function zeroPad(n, width) {
        var s = n.toString();
        var i = Math.max(width - s.length, 0);
        return new Array(i + 1).join("0") + s;
    }

    /**
     * Calculate distortion of the wind vector caused by the shape of the projection at point (x, y). The wind
     * vector is modified in place and returned by this function.
     */
    function distort(projection, λ, φ, x, y, scale, wind) {
        var u = wind[0] * scale;
        var v = wind[1] * scale;
        var d = distortion(projection, λ, φ, x, y);

        // Scale distortion vectors by u and v, then add.
        wind[0] = d[0] * u + d[2] * v;
        wind[1] = d[1] * u + d[3] * v;
        return wind;
    }

    function proportion(i, bounds) {
        return (clamp(i, bounds) - bounds[0]) / (bounds[1] - bounds[0]);
    }


   function distortion(projection, λ, φ, x, y) {
        var τ = 2 * Math.PI;
        var H = Math.pow(10, -5.2);
        var hλ = λ < 0 ? H : -H;
        var hφ = φ < 0 ? H : -H;
        var pλ = project([λ + hλ, φ]);
        var pφ = project([λ, φ + hφ]);
        //var pλ = projection([λ + hλ, φ]);
        //var pφ = projection([λ, φ + hφ]);
        
        // Meridian scale factor (see Snyder, equation 4-3), where R = 1. This handles issue where length of 1º λ
        // changes depending on φ. Without this, there is a pinching effect at the poles.
        var k = Math.cos(φ / 360 * τ);
        return [
            (pλ[0] - x) / hλ / k,
            (pλ[1] - y) / hλ / k,
            (pφ[0] - x) / hφ,
            (pφ[1] - y) / hφ
        ];
    }


  
  function createField(columns, bounds, callback) {

        /**
         * @returns {Array} wind vector [u, v, magnitude] at the point (x, y), or [NaN, NaN, null] if wind
         *          is undefined at that point.
         */
        function field(x, y) {
            var column = columns[Math.round(x)];
            return column && column[Math.round(y)] || NULL_WIND_VECTOR;
        }

        // Frees the massive "columns" array for GC. Without this, the array is leaked (in Chrome) each time a new
        // field is interpolated because the field closure's context is leaked, for reasons that defy explanation.
        field.release = function() {
            columns = [];
        };

        field.randomize = function(o) {  // UNDONE: this method is terrible
            var x, y;
            var safetyNet = 0;
            do { 
                x = Math.round(Math.floor(Math.random() * bounds.width) + bounds.x);
                y = Math.round(Math.floor(Math.random() * bounds.height) + bounds.y)
            } while (field(x, y)[2] === null && safetyNet++ < 30);
            o.x = x;
            o.y = y;
            return o;
        };

        //field.overlay = mask.imageData;
        //return field;
        callback( bounds, field );
  }



  function buildBounds( bounds, width, height ) {
      var upperLeft = bounds[0];
      var lowerRight = bounds[1];
      var x = Math.round(upperLeft[0]); //Math.max(Math.floor(upperLeft[0], 0), 0);
      var y = Math.max(Math.floor(upperLeft[1], 0), 0);
      var xMax = Math.min(Math.ceil(lowerRight[0], width), width - 1);
      var yMax = Math.min(Math.ceil(lowerRight[1], height), height - 1);

      return {x: x, y: y, xMax: width, yMax: yMax, width: width, height: height};
  }

  function currentPosition() {
        var λ = floorMod(new Date().getTimezoneOffset() / 4, 360);  // 24 hours * 60 min / 4 === 360 degrees
        return [λ, 0];
  }

  var invert = function(x) {
    var point = map.pointLocation(new MM.Point(x[1], x[0]));
    return [point.lon, point.lat];
  };

  var project = function(x) {
    var point = map.locationPoint(new MM.Location(x[1], x[0]));
    return [point.x, point.y];
  };

  
  function ensureNumber(num, fallback) {
    return _.isFinite(num) || num === Infinity || num === -Infinity ? num : fallback;
  }

  function interpolateField( grid, bounds, callback ) {

    var projection = d3.geo.mercator().precision(.1); 
    var velocityScale = bounds.height * VELOCITY_SCALE;

    var columns = [];
    var point = [];
    var x = bounds.x;
    var scale = { bounds: bounds }; //grid.recipe.scale; //, gradient = scale.gradient;

    function interpolateColumn(x) {
        var column = [];
        for (var y = bounds.y; y <= bounds.yMax; y += 2) {
                point[1] = x; point[0] = y;
                var coord = invert(point);
                var color = TRANSPARENT_BLACK;
                if (coord) {
                    var λ = coord[0], φ = coord[1];
                    if (isFinite(λ)) {
                        var wind = grid.interpolate(λ, φ);
                        if (wind) {
                            wind = distort(projection, λ, φ, x, y, velocityScale, wind);
                            column[y+1] = column[y] = wind;
                            //color = gradient(proportion(wind[2], scale.bounds), OVERLAY_ALPHA);
                        }
                    }
                }
                //mask.set(x, y, color).set(x+1, y, color).set(x, y+1, color).set(x+1, y+1, color);
        }
        columns[x+1] = columns[x] = column;
    }

    (function batchInterpolate() {
        try {
                var start = Date.now();
                //console.log('start while', x, bounds.width);
                while (x < bounds.width) {
                    interpolateColumn(x);
                    x += 2;
                    if ((Date.now() - start) > 400) { //MAX_TASK_TIME) {
                        // Interpolation is taking too long. Schedule the next batch for later and yield.
                        //report.progress((x - bounds.x) / (bounds.xMax - bounds.x));
                        setTimeout(batchInterpolate, 25);
                        return;
                    }
                }
          createField(columns, bounds, callback);
        }
        catch (e) {
            console.log('error in batch interp', e);
        }
    })();
  }


  var animate = function(bounds, field) {

    function asColorStyle(r, g, b, a) {
        return "rgba(" + r + ", " + g + ", " + b + ", " + a + ")";
    }
  
    function windIntensityColorScale(step, maxWind) {
        var result = [];
        for (var j = 200; j >= 0; j = j - step) {
            result.push(asColorStyle(j, j, j, .5));
        }
        result.indexFor = function(m) {  // map wind speed to a style
            return Math.floor(Math.min(m, maxWind) / maxWind * (result.length - 1));
        };
        return result;
    }

    var colorStyles = windIntensityColorScale(INTENSITY_SCALE_STEP, MAX_WIND_INTENSITY);
    var buckets = colorStyles.map(function() { return []; });
    var particleCount = Math.round(bounds.width * PARTICLE_MULTIPLIER);
    //if (µ.isMobile()) {
      //particleCount *= .75; //PARTICLE_REDUCTION;
    //}
    
    var fadeFillStyle = "rgba(0, 0, 0, 0.97)";

    var particles = [];
    for (var i = 0; i < particleCount; i++) {
        particles.push(field.randomize({age: Math.floor(Math.random() * MAX_PARTICLE_AGE) + 0}));
    }

    function evolve() {
        buckets.forEach(function(bucket) { bucket.length = 0; });
        particles.forEach(function(particle) {
            if (particle.age > MAX_PARTICLE_AGE) {
                field.randomize(particle).age = 0;
            }
            var x = particle.x;
            var y = particle.y;
            var v = field(x, y);  // vector at current position
            var m = v[2];
            if (m === null) {
                particle.age = MAX_PARTICLE_AGE;  // particle has escaped the grid, never to return...
            }
            else {
                var xt = x + v[0];
                var yt = y + v[1];
                if (field(xt, yt)[2] !== null) {
                    // Path from (x,y) to (xt,yt) is visible, so add this particle to the appropriate draw bucket.
                    particle.xt = xt;
                    particle.yt = yt;
                    buckets[colorStyles.indexFor(m)].push(particle);
                }
                 else {
                    // Particle isn't visible, but it still moves through the field.
                    particle.x = xt;
                    particle.y = yt;
                }
            }
            particle.age += 1;
        });
    }

    var g = params.canvas.getContext("2d");
    g.lineWidth = PARTICLE_LINE_WIDTH;
    g.fillStyle = fadeFillStyle;

    function draw() {
        // Fade existing particle trails.
        var prev = g.globalCompositeOperation;
        g.globalCompositeOperation = "destination-in";
        g.fillRect(bounds.x, bounds.y, bounds.width, bounds.height);
        g.globalCompositeOperation = prev;

        // Draw new particle trails.
        buckets.forEach(function(bucket, i) {
            if (bucket.length > 0) {
                g.beginPath();
                g.strokeStyle = colorStyles[i];
                bucket.forEach(function(particle) {
                    g.moveTo(particle.x, particle.y);
                    g.lineTo(particle.xt, particle.yt);
                    particle.x = particle.xt;
                    particle.y = particle.yt;
                });
                g.stroke();
            }
        });
    }

    (function frame() {
        try {
            evolve();
            draw();
            setTimeout(frame, FRAME_RATE);
        }
        catch (e) {
            console.error(e);
        }
    })();
  }

  function start( bounds, width, height ){
    stop();
    // build grid 
    buildGrid( params.data, function(grid){
      // interpolateField
      // create Field
      interpolateField( grid, buildBounds( bounds, width, height), function( bounds, field ){
        // animate the canvas with random points 
        windy.field = field;
        animate( bounds, field );
      });
    });
  }

  function stop(){
    if (windy.field) windy.field.release();
  }

  var windy = {
    params: params,
    start: start,
    stop: stop
  };

  return windy;
}

//
if ( typeof(module) != "undefined" ){
  module.exports = Windy;
}
