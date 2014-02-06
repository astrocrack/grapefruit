// Based on:
//  http://devblog.phillipspiess.com/2010/02/23/better-know-an-algorithm-1-marching-squares/
//  https://github.com/ralphbean/js-experiments/blob/master/marching-squares/index.html

/**
 * Marching squares movement chart: http://devblog.phillipspiess.com/wp-content/uploads/2010/02/MarchingStates.png
 */

var inherit = require('../utils/inherit'),
    math = require('./math'),
    Vector = require('./Vector'),
    Polygon = require('../geom/Polygon'),
    C = require('../constants'),
    cp = require('../vendor/cp');

function TiledGeom(layer) {
    this.data = layer.tileIds;
    this.size = layer.map.size;
    //this.testFn = testFn;
    this.map = layer.map;

    this._startPoint = new Vector();
    this._previousStep = null;
    this._nextStep = null;
    this._state = null;

    this.alpha = 0.5;

    this._verts = [
        new Vector(),
        new Vector(),
        new Vector(),
        new Vector()
    ];
    this._cell = {
        x: 0,
        y: 0,
        v: [
            new Vector(),
            new Vector(),
            new Vector(),
            new Vector()
        ]
    };
}

inherit(TiledGeom, Object, {
    create: function() {
        var poly,
            verts = {},
            points = [],
            isolevel = 0.0008,
            tid, x, y,
            szx = this.size.x,
            tsx = this.map.tileSize.x,
            tsy = this.map.tileSize.y,
            pts,
            nPts,
            pt;

        //iterate through each tile and accumulate their verticies
        //eliminating duplicates along the way.
        for(var i = 0, len = this.data.length; i < len; ++i) {
            //get the tile
            tid = this.data[i];
            poly = this.getTilePoly(tid);

            //if it has a poly, add it to the segment pool
            if(poly) {
                x = (i % szx);// * tsx;
                y = math.floor(i / szx);// * tsy;

                if(poly._shapetype === C.SHAPE.RECTANGLE) {
                    poly = poly.toPolygon();
                }

                pts = poly.points;
                nPts = pts.length;

                //add the vertices, filtering duplicates
                for(var p = 0; p < nPts; ++p) {
                    pt = pts[p];

                    if(!pt._worldValues) {
                        pt.x += x;
                        pt.y += y;
                        pt._worldValues = true;
                    }

                    var id = pt.toHashkey();

                    if(!verts[id]) {
                        verts[id] = 1;
                        points.push(pt);
                    }
                }
            }
        }

        //TODO: Convert points to an array of triangle polygons via https://github.com/r3mi/poly2tri.js

        console.log(points.map(function(v) {
            return v.x + ',' + v.y;
        }).join(','));

        //reduce points to a convex hull
        var newPoints = [];

        cp.convexHull(points, newPoints);

        console.log(newPoints);

        return new Polygon(0, 0, newPoints);
    },
    getTilePoly: function(tid) {
        if(!tid) return false;

        var set = this.map.getTileset(tid),
            props;

        if(set) {
            props = set.getTileProperties(tid);

            if(props.mass) {
                return props.hitArea || new Polygon(0, 0, [
                    new Vector(0, 0),
                    new Vector(set.tileSize.x, 0),
                    new Vector(set.tileSize.x, set.tileSize.y),
                    new Vector(0, set.tileSize.y)
                ]);
            }
        }

        return false;
    },
    simplify: function(segments, tolerance) {
        if(!segments.length) return [];

        tolerance = tolerance || 0;

        var deduped = [segments[0]],
            prev = deduped[0],
            seg;

        if(segments.length === 1) return deduped;

        //any segments that are the same coors

        //combine adjacent segments, if previous b === current a then extend
        //the previous segment to include both.
        for(var i = 1, len = segments.length; i < len; ++i) {
            seg = segments[i];

            //if duplicated, combine
            if(math.abs(prev.b.x - seg.a.x) < tolerance && math.abs(prev.b.y - seg.a.y) < tolerance) {
                prev.b = seg.b;
            }
            //otherwise add this segment
            else {
                deduped.push(seg);
                prev = seg;
            }
        }

        //finally check the wrap point
        var f = deduped[0],
            l = deduped[deduped.length - 1];

        if(math.abs(l.b.x - f.a.x) < tolerance && math.abs(l.b.y - f.a.y) < tolerance) {
            l.b = f.b;
            deduped.shift();
        }

        return deduped;
    },
    /*
    //public-facing function that returns an array of segments from the march
    march: function() {
        this._findStartPoint();

        return this._walkPerimeter(this._startPoint.x, this._startPoint.y);
    },
    //remove duplicated segments from a previous march
    simplify: function(segments, tolerance) {
        if(!segments.length) return [];

        tolerance = tolerance || 0;

        var deduped = [segments[0]],
            prev = deduped[0],
            seg;

        if(segments.length === 1) return deduped;

        //combine adjacent segments, if previous b === current a then extend
        //the previous segment to include both.
        for(var i = 1, len = segments.length; i < len; ++i) {
            seg = segments[i];

            //if duplicated, combine
            if(math.abs(prev.b.x - seg.a.x) < tolerance && math.abs(prev.b.y - seg.a.y) < tolerance) {
                prev.b = seg.b;
            }
            //otherwise add this segment
            else {
                deduped.push(seg);
                prev = seg;
            }
        }

        //finally check the wrap point
        var f = deduped[0],
            l = deduped[deduped.length - 1];

        if(math.abs(l.b.x - f.a.x) < tolerance && math.abs(l.b.y - f.a.y) < tolerance) {
            l.b = f.b;
            deduped.shift();
        }

        return deduped;
    },
    //find a non-empty point in the data
    _findStartPoint: function() {
        //iterate to find a non-empty point
        for(var i = 0; i < this.data.length; ++i) {
            if(this.data[i] !== 0) {
                return this._startPoint.set((i % this.size.x) - 1, math.floor(i / this.size.x) - 1);
            }
        }

        return Vector.ZERO;
    },
    _walkPerimeter: function(startX, startY) {
        //sanity checking
        if(startX < 0) startX = 0;
        else if(startX > this.size.x) startX = this.size.x;

        if(startY < 0) startY = 0;
        else if(startY > this.size.y) startY = this.size.y;

        //setup return list
        var segments = [];

        //current x/y position during walk
        var x = startX,
            y = startY;

        //main while loop, steps until reaching initial points
        do {
            //evaluate state and setup next direction
            this._step(x, y, segments);

            switch(this._nextStep) {
                case C.DIRECTION.UP:    y--; break;
                case C.DIRECTION.LEFT:  x--; break;
                case C.DIRECTION.DOWN:  y++; break;
                case C.DIRECTION.RIGHT: x++; break;
            }
        } while(x !== startX || y !== startY);

        return segments;
    },
    //Determines and sets the state of the 4 tiles that represent our current state,
    //and sets our current and previous directions
    _step: function(x, y, segments) {
        //store previous step
        this._previousStep = this._nextStep;

        //determine state by scanning the 4 tile area
        this._state = 0;

        if(this.testFn(x, y)) //up, left
            this._state |= 1;
        if(this.testFn(x+1, y)) //up, right
            this._state |= 2;
        if(this.testFn(x+1, y+1)) //down, right
            this._state |= 4;
        if(this.testFn(x, y+1)) //down, left
            this._state |= 8;

        //state is now in the range [0, 15]. Use the switch statement to determine
        //what our next step is based on the Marching Squares Movement Chart (below).
        switch(this._state) {
            case 1: this._nextStep = C.DIRECTION.UP; break;
            case 2: this._nextStep = C.DIRECTION.RIGHT; break;
            case 3: this._nextStep = C.DIRECTION.RIGHT; break;
            case 4: this._nextStep = C.DIRECTION.DOWN; break;
            case 5: this._nextStep = C.DIRECTION.UP; break;
            case 6:
                if(this._previousStep === C.DIRECTION.UP) {
                    this._nextStep = C.DIRECTION.LEFT;
                }
                else {
                    this._nextStep = C.DIRECTION.RIGHT;
                }
                break;
            case 7: this._nextStep = C.DIRECTION.DOWN; break;
            case 8: this._nextStep = C.DIRECTION.LEFT; break;
            case 9:
                if(this._previousStep === C.DIRECTION.RIGHT) {
                    this._nextStep = C.DIRECTION.UP;
                }
                else {
                    this._nextStep = C.DIRECTION.DOWN;
                }
                break;
            case 10: this._nextStep = C.DIRECTION.LEFT; break;
            case 11: this._nextStep = C.DIRECTION.RIGHT; break;
            case 12: this._nextStep = C.DIRECTION.LEFT; break;
            case 13: this._nextStep = C.DIRECTION.UP; break;
            case 14: this._nextStep = C.DIRECTION.LEFT; break;
            default: this._nextStep = C.DIRECTION.NONE; break;
        }

        if(this._state === 0)
            return segments || [];

        var edge = TiledGeom.edgeTable[this._state],
            seg = TiledGeom.segmentTable[this._state],
            alpha = this.alpha;

        this._verts[0].set(0, 0);
        this._verts[1].set(0, 0);
        this._verts[2].set(0, 0);
        this._verts[3].set(0, 0);

        this._cell[0].set(x, y);
        this._cell[1].set(x+1, y);
        this._cell[2].set(x+1, y+1);
        this._cell[3].set(x, y+1);

        if(edge & 1) {
            this._verts[0].copy(this._cell[0]).lerp(this._cell[1], alpha);
        }

        if(edge & 2) {
            this._verts[1].copy(this._cell[1]).lerp(this._cell[2], alpha);
        }

        if(edge & 4) {
            this._verts[2].copy(this._cell[2]).lerp(this._cell[3], alpha);
        }

        if(edge & 8) {
            this._verts[3].copy(this._cell[3]).lerp(this._cell[0], alpha);
        }

        segments = segments || [];
        for(var i = 0; seg[i] !== -1; i += 2) {
            segments.push({
                a: this._verts[seg[i]],
                b: this._verts[seg[i+1]]
            });
        }

        return segments;
    }
    */
    /*
    //Simplify is based on PathFitter.fit() from paper.js
    // https://github.com/paperjs/paper.js/blob/master/src/path/PathFitter.js
    simplify: function(segments, tolerance) {
        if(segments.length < 2) return segments;

        var points = [],
            a, b,
            prev;

        //copy points and filter out adjacent duplicates
        for(var i = 0, len = segments.length; i < len; ++i) {
            a = segments[i].a;
            b = segments[i].b;

            if(!prev || !prev.equals(a)) {
                points.push(a);
                prev = a;
            }

            if(!prev.equals(b)) {
                points.push(b);
                prev = b;
            }
        }

        var plen = points.length;

        if(plen < 2) return points;

        tolerance = tolerance || 0;

        var simplified = [];
        this.simplifyCubic(
            points,
            simplified,
            0,
            plen - 1,
            points[1].clone().sub(points[0]).normalize(),
            points[plen - 2].clone().sub(points[plen - 1]).normalize()
        );

        return simplified;
    },
    simplifyCubic: function(points, out, first, last, tan1, tan2, tolerance) {
        out = out || [];

        //Use heuristic if region only has two points in it
        if(last - first == 1) {
            var pt1 = points[first],
                pt2 = points[last],
                dist = pt1.distanceTo(pt2) / 3;

            this.addCurve([
                pt1,
                pt1.add(tan1.normalize(dist)),
                pt2.add(tan2.normalize(dist)),
                pt2
            ], out);
            return out;
        }

        // Parameterize points, and attempt to fit curve
        var uPrime = this.chordLengthParameterize(points, first, last),
            maxError = Math.max(tolerance, tolerance * tolerance),
            split,
            curve = [0, 0, 0, 0],
            max = {
                error: maxDist,
                index: index
            };

        // Try 4 iterations
        for(var i = 0; i <= 4; i++) {
            this.generateBezier(points, curve, first, last, uPrime, tan1, tan2);

            //  Find max deviation of points to fitted curve
            this.findMaxError(points, max, first, last, curve, uPrime);

            if(max.error < tolerance) {
                this.addCurve(curve, out);
                return out;
            }

            split = max.index;

            // If error not too large, try reparameterization and iteration
            if(max.error >= maxError)
                break;

            this.reparameterize(points, first, last, uPrime, curve);
            maxError = max.error;
        }

        // Fitting failed -- split at max error point and fit recursively
        var V1 = points[split - 1].clone().sub(points[split]),
            V2 = points[split].clone().sub(points[split + 1]),
            tanCenter = V1.add(V2).divideScalar(2).normalize();

        this.fitCubic(first, split, tan1, tanCenter);
        this.fitCubic(split, last, tanCenter.negate(), tan2);
    },
    addCurve: function(curve, out) {
        out.push({
            a: curve[3],
            b: curve[2].sub(curve[3])
        });
    },
    generateBezier: function(points, out, first, last, uPrime, tan1, tan2) {
        out = out || [0, 0, 0, 0];

        var epsilon = 10e-12,
            p1 = points[first],
            p2 = points[last],
            C = [[0, 0], [0,0]],
            X = [0,0],
            tmp = new Vector(),
            ttn1 = new Vector(),
            ttn2 = new Vector(),
            tpt1 = new Vector(),
            tpt2 = new Vector();

        for(var i = 0, len = last - first + 1; i < len; ++i) {
            var u = uPrime[i],
                t = 1 - u,
                b = 3 * u * t,
                b0 = t * t * t,
                b1 = b * t,
                b2 = b * u,
                b3 = u * u * u,
                a1 = ttn1.copy(tan1).normalize(b1),
                a2 = ttn2.copy(tan2).normalize(b2),
                tmp = tmp.copy(points[first + i])
                    .sub(tpt1.copy(pt1).multiplyScalar(b0 + b1))
                    .sub(tpt2.copy(pt2).multiplyScalar(b2 + b3));

            C[0][0] += a1.dot(a1);
            C[0][1] += a1.dot(a2);
            // C[1][0] += a1.dot(a2);
            C[1][0] = C[0][1];
            C[1][1] += a2.dot(a2);
            X[0] += a1.dot(tmp);
            X[1] += a2.dot(tmp);
        }

        // Compute the determinants of C and X
        var detC0C1 = C[0][0] * C[1][1] - C[1][0] * C[0][1],
            alpha1, alpha2;

        if(Math.abs(detC0C1) > epsilon) {
            // Kramer's rule
            var detC0X  = C[0][0] * X[1]    - C[1][0] * X[0],
                detXC1  = X[0]    * C[1][1] - X[1]    * C[0][1];

            // Derive alpha values
            alpha1 = detXC1 / detC0C1;
            alpha2 = detC0X / detC0C1;
        } else {
            // Matrix is under-determined, try assuming alpha1 == alpha2
            var c0 = C[0][0] + C[0][1],
                c1 = C[1][0] + C[1][1];

            if(Math.abs(c0) > epsilon) {
                alpha1 = alpha2 = X[0] / c0;
            } else if(Math.abs(c1) > epsilon) {
                alpha1 = alpha2 = X[1] / c1;
            } else {
                // Handle below
                alpha1 = alpha2 = 0;
            }
        }

        // If alpha negative, use the Wu/Barsky heuristic (see text)
        // (if alpha is 0, you get coincident control points that lead to
        // divide by zero in any subsequent NewtonRaphsonRootFind() call.
        var segLength = pt2.getDistance(pt1);
        epsilon *= segLength;

        if(alpha1 < epsilon || alpha2 < epsilon) {
            // fall back on standard (probably inaccurate) formula,
            // and subdivide further if needed.
            alpha1 = alpha2 = segLength / 3;
        }

        // First and last control points of the Bezier curve are
        // positioned exactly at the first and last data points
        // Control points 1 and 2 are positioned an alpha distance out
        // on the tangent vectors, left and right, respectively
        out[0] = pt1;
        out[1] = tpt1.copy(pt1).add(ttn1.copy(tan1).normalize(alpha1));
        out[2] = pt2;
        out[3] = tpt2.copy(pt2).add(ttn2.copy(tan2).normalize(alpha2));
        return out;
    },
    reparameterize: function(points, first, last, u, curve) {
        for(var i = first; i <= last; ++i) {
            u[i - first] = this.findRoot(curve, points[i], u[i - first]);
        }
    },
    // Use Newton-Raphson iteration to find better root.
    findRoot: function(curve, point, u) {
        var curve1 = [],
            curve2 = [],
            i = 0;

        // Generate control vertices for Q'
        for(i = 0; i <= 2; ++i) {
            curve1[i] = curve[i + 1].clone().sub(curve[i]).multiply(3);
        }

        // Generate control vertices for Q''
        for(i = 0; i <= 1; ++i) {
            curve2[i] = curve1[i + 1].clone().sub(curve1[i]).multiply(2);
        }

        // Compute Q(u), Q'(u) and Q''(u)
        var pt = this.evaluate(3, curve, u),
            pt1 = this.evaluate(2, curve1, u),
            pt2 = this.evaluate(1, curve2, u),
            diff = pt.sub(point),
            df = pt1.dot(pt1) + diff.dot(pt2);

        // Compute f(u) / f'(u)
        if (Math.abs(df) < 10e-6) //tolerance
            return u;

        // u = u - f(u) / f'(u)
        return u - diff.dot(pt1) / df;
    },
    // Evaluate a Bezier curve at a particular parameter value
    evaluate: function(degree, curve, t) {
        // Copy array
        var tmp = curve.slice(),
            tv = new Vector();

        // Triangle computation
        for(var i = 1; i <= degree; i++) {
            for(var j = 0; j <= degree - i; j++) {
                tmp[j] = tmp[j].clone().multiplyScalar(1 - t).add(tv.copy(tmp[j + 1]).multiplyScalar(t));
            }
        }

        return tmp[0];
    },
    // Assign parameter values to digitized points
    // using relative distances between points.
    chordLengthParameterize: function(points, first, last) {
        var u = [0],
            i = 0;

        for(i = first + 1; i <= last; i++) {
            u[i - first] = u[i - first - 1] + points[i].distanceTo(points[i - 1]);
        }

        for(i = 1, m = last - first; i <= m; i++) {
            u[i] /= u[m];
        }

        return u;
    },
    // Find the maximum squared distance of digitized points to fitted curve.
    findMaxError: function(points, out, first, last, curve, u) {
        out = out || {
            error: maxDist,
            index: index
        };

        var index = Math.floor((last - first + 1) / 2),
            maxDist = 0;

        for(var i = first + 1; i < last; i++) {
            var P = this.evaluate(3, curve, u[i - first]),
                v = P.sub(points[i]),
                dist = v.lengthSq();

            if(dist >= maxDist) {
                maxDist = dist;
                index = i;
            }
        }

        out.error = maxDist;
        out.index = index;
        return out;
    }
    */
});

TiledGeom.edgeTable = [
    0x0,     //0000
    0x9,     //1001
    0x3,     //0011
    0xa,     //1010
    0x6,     //0110
    0xf,     //1111
    0x5,     //0101
    0xc,     //1100
    0xc,     //1100
    0x5,     //0101
    0xf,     //1111
    0x6,     //0110
    0xa,     //1010
    0x3,     //0011
    0x9,     //1001
    0x0,     //0000
];

TiledGeom.segmentTable = [
    [-1,-1,-1,-1,-1],
    [0,3,-1,-1,-1],
    [1,0,-1,-1,-1],
    [1,3,-1,-1,-1],
    [2,1,-1,-1,-1],
    [2,1,0,3,-1],
    [2,0,-1,-1,-1],
    [2,3,-1,-1,-1],
    [3,2,-1,-1,-1],
    [0,2,-1,-1,-1],
    [1,0,3,2,-1],
    [1,2,-1,-1,-1],
    [3,1,-1,-1,-1],
    [0,1,-1,-1,-1],
    [3,0,-1,-1,-1],
    [-1,-1,-1,-1,-1]
];

module.exports = TiledGeom;
