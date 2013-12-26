

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
  var PARTICLE_LINE_WIDTH = 1.0;            // line width of a drawn particle
  var PARTICLE_MULTIPLIER = 7;              // particle count scalar (completely arbitrary--this values looks nice)
  var PARTICLE_REDUCTION = 0.75;            // reduce particle count to this much of normal for mobile devices
  var FRAME_RATE = 40;                      // desired milliseconds per frame

  var NULL_WIND_VECTOR = [NaN, NaN, null];  // singleton for no wind in the form: [u, v, magnitude]
  var TRANSPARENT_BLACK = [0, 0, 0, 0];



  var grid = {

    bilinearInterpolateScalar: function(x, y, g00, g10, g01, g11) {
        var rx = (1 - x);
        var ry = (1 - y);
        return g00 * rx * ry + g10 * x * ry + g01 * rx * y + g11 * x * y;
    },

    bilinearInterpolateVector: function(x, y, g00, g10, g01, g11) {
        var rx = (1 - x);
        var ry = (1 - y);
        var a = rx * ry,  b = x * ry,  c = rx * y,  d = x * y;
        var u = g00[0] * a + g10[0] * b + g01[0] * c + g11[0] * d;
        var v = g00[1] * a + g10[1] * b + g01[1] * c + g11[1] * d;
        return [u, v, Math.sqrt(u * u + v * v)];
    },

    createScalarBuilder: function(record) {
        var data = record.data;
        return {
            header: record.header,
            recipe: recipeFor(""),
            data: function(i) {
                return data[i];
            },
            interpolate: bilinearInterpolateScalar
        }
    },

    createWindBuilder: function(uComp, vComp) {
        var uData = uComp.data, vData = vComp.data;
        return {
            header: uComp.header,
            recipe: recipeFor("wind-" + uComp.header.surface1Value),
            data: function(i) {
                return [uData[i], vData[i]];
            },
            interpolate: bilinearInterpolateVector
        }
    },

    createBuilder: function(data) {
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
    },

    buildGrid: function(data) {
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
            var i = µ.floorMod(λ - λ0, 360) / Δλ;  // calculate longitude index in wrapped range [0, 360)
            var j = (φ0 - φ) / Δφ;                 // calculate latitude index in direction +90 to -90

            var fi = Math.floor(i), ci = fi + 1;
            var fj = Math.floor(j), cj = fj + 1;

            var row;
            if ((row = grid[fj])) {
                var g00 = row[fi];
                var g10 = row[ci];
                if (µ.isValue(g00) && µ.isValue(g10) && (row = grid[cj])) {
                    var g01 = row[fi];
                    var g11 = row[ci];
                    if (µ.isValue(g01) && µ.isValue(g11)) {
                        // All four points found, so interpolate the value.
                        return builder.interpolate(i - fi, j - fj, g00, g10, g01, g11);
                    }
                }
            }
            return null;
        }

        return {
            date: date,
            recipe: builder.recipe,
            interpolate: interpolate
        };
    }

  }


  var animate = function(globe, field) {
    if (!globe || !field) return;

//    var cancel = this.cancel;
//    var bounds = globe.bounds(view);
    var colorStyles = µ.windIntensityColorScale(INTENSITY_SCALE_STEP, MAX_WIND_INTENSITY);
    var buckets = colorStyles.map(function() { return []; });
    var particleCount = Math.round(bounds.width * PARTICLE_MULTIPLIER);
    if (µ.isMobile()) {
        particleCount *= PARTICLE_REDUCTION;
    }
    var fadeFillStyle = µ.isFF() ? "rgba(0, 0, 0, 0.95)" : "rgba(0, 0, 0, 0.97)";  // FF Mac alpha behaves oddly

    log.debug("particle count: " + particleCount);
    var particles = [];
    for (var i = 0; i < particleCount; i++) {
        particles.push(field.randomize({age: _.random(0, MAX_PARTICLE_AGE)}));
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

    var g = d3.select("#animation").node().getContext("2d");
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
            if (cancel.requested) {
                field.release();
                return;
            }
            evolve();
            draw();
            setTimeout(frame, FRAME_RATE);
        }
        catch (e) {
            report.error(e);
        }
    })();
  }

  var windy = {
    params: params,
    grid: grid,
    animate: animate
  };
  return windy;
}

//
if ( typeof(module) != "undefined" ){
  module.exports = Windy;
}
