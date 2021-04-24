const area = 300; // width and height of an SVG
const scale = 100; // radius of a beachball
function round(x) { return Math.round(x * 100) / 100; }

const N = 12; // number of points per 180 degrees
const EPSIL = 1e-4;

// color settings
var penW = { stroke: "red", fill: "none" };
var penL = penW;
var penT = penW;
var fillE = "lightyellow";
var fillG = "gray";
var penFp = { stroke: "none", fill: fillG }; // -Fp & -Fg
var penFt = { stroke: "none", fill: fillE }; // -Ft & -Fe

var draw = SVG().addTo('#svg-area').size(area, area);
var group = draw.group().transform({ tx: area / 2, ty: area / 2 }).flip("y");

function deg2rad(theta) { return theta * (Math.PI / 180); }
function rad2deg(theta) { return theta / (Math.PI / 180); }

function polar2xy(radius, theta) {
	theta = deg2rad(theta);
	return [radius * Math.cos(theta), radius * Math.sin(theta)];
}

function polar2zne(radius, azimuth, plunge) {
	if (radius == 0) return [0, 0, 0];
	let [rh, rz] = polar2xy(radius, plunge);
	let [rn, re] = polar2xy(rh, azimuth);
	return [rz, rn, re]; // z (depth), northward, eastward
}

class Axis {
	// z * z + n * n + e * e == 1
	constructor(z, n, e, value = undefined) {
		this.z = z; // z (depth)
		this.n = n; // northward
		this.e = e; // eastward
		this.value = value;
	}

	static fromPolar(azimuth, plunge, value = undefined) {
		let [z, n, e] = polar2zne(1, azimuth, plunge);
		return new Axis(z, n, e, value);
	}

	getPolar() {
		let { z, n, e } = this;

		let h2 = n * n + e * e;
		if (h2 == 0) return [0, z >= 0 ? 90 : -90];

		return [
			rad2deg(Math.atan2(e, n)),
			rad2deg(Math.atan2(z, Math.sqrt(h2))),
		];
	}

	getOpposite() {
		return new Axis(-this.z, -this.n, -this.e, this.value);
	}

	static sumsub(a1, a2) {
		// assume that a1 and a2 are orthonormal
		return [
			new Axis(
				(a1.z + a2.z) * Math.SQRT1_2,
				(a1.n + a2.n) * Math.SQRT1_2,
				(a1.e + a2.e) * Math.SQRT1_2
			),
			new Axis(
				(a1.z - a2.z) * Math.SQRT1_2,
				(a1.n - a2.n) * Math.SQRT1_2,
				(a1.e - a2.e) * Math.SQRT1_2
			)
		];
	}
}

class NodalPlane {
	constructor(normalVector, slipVector = undefined) {
		this.nv = normalVector;
		this.sv = slipVector;
	}

	static fromParams(str, dip, rake = undefined) {
		let nv = Axis.fromPolar(str + 90, dip - 90); // nv.z <= 0
		let sv;

		if (rake != null) {
			let [cs, ss] = polar2xy(1, str);
			let [cd, sd] = polar2xy(1, dip);
			let [cr, sr] = polar2xy(1, rake);
			sv = new Axis(-sd * sr, cs * cr + ss * cd * sr, ss * cr - cs * cd * sr);
		}

		return new NodalPlane(nv, sv);
	}

	static fromPartial(str1, dip1, str2, fault) {
		let nv = NodalPlane.fromParams(str1, dip1).nv; // nv.z <= 0
		let sv;

		let sz = 0, [sn, se] = polar2xy(1, str2);
		let [z, n, e] = [nv.n * se - nv.e * sn, nv.e * sz - nv.z * se, nv.z * sn - nv.n * sz];
		let norm = Math.hypot(z, n, e);

		if (norm > 0) {
			sv = new Axis(z / norm, n / norm, e / norm);
			if (sv.z * fault < 0) sv = sv.getOpposite();
		} else {
			console.warn("cannot determine NP2 (any dip is possible)");
		}

		return new NodalPlane(nv, sv);
	}

	static fromAxes(tAxis, pAxis) {
		let [nv, sv] = Axis.sumsub(tAxis, pAxis);

		return (nv.z <= 0)
			? new NodalPlane(nv, sv)
			: new NodalPlane(nv.getOpposite(), sv.getOpposite());
	}

	getParams() {
		let { nv, sv } = this;
		let str, dip, rake;
		let sz, sn, se; // vector pointing to strike

		let { z, n, e } = nv;
		let h2 = n * n + e * e;
		if (h2 == 0) {
			str = 0; // arbitrary
			dip = 0;
			[sz, sn, se] = [0, 1, 0];
		} else {
			str = rad2deg(Math.atan2(-n, e));
			dip = rad2deg(Math.atan2(Math.sqrt(h2), -z));
			[sz, sn, se] = [0, e, -n];
		}

		if (sv != null) {
			// assert that nv.norm == 1
			let [tz, tn, te] = [n * se - e * sn, e * sz - z * se, z * sn - n * sz];
			({ z, n, e } = sv);
			rake = rad2deg(Math.atan2(tz * z + tn * n + te * e, sz * z + sn * n + se * e));
		}

		return [str, dip, rake];
	}

	getOtherPlane() {
		if (this.sv == null) throw "rake is undefined";

		// swap two vectors, new nv.z <= 0
		return (this.sv.z <= 0)
			? new NodalPlane(this.sv, this.nv)
			: new NodalPlane(this.sv.getOpposite(), this.nv.getOpposite());
	}

	getPrincipalAxes() {
		if (this.sv == null) throw "rake is undefined";

		// assert that nv.norm == sv.norm == 1 and dot(nv,sv) == 0
		let { nv, sv } = this;

		let [tAxis, pAxis] = Axis.sumsub(nv, sv);
		tAxis.value = +1;
		pAxis.value = -1;
		let nAxis = new Axis(
			nv.n * sv.e - nv.e * sv.n,
			nv.e * sv.z - nv.z * sv.e,
			nv.z * sv.n - nv.n * sv.z,
			0
		);

		return [
			(tAxis.z >= 0) ? tAxis : tAxis.getOpposite(),
			(nAxis.z >= 0) ? nAxis : nAxis.getOpposite(),
			(pAxis.z >= 0) ? pAxis : pAxis.getOpposite(),
		];
	}
}

class FocalMechanism {
	constructor(np, axes) {
		this.np = np;
		this.axes = axes;
	}

	static fromSa(str, dip, rake, ...rest) {
		let np = NodalPlane.fromParams(str, dip, rake);
		let axes = np.getPrincipalAxes();
		return new FocalMechanism(np, axes);
	}

	static fromSc(str1, dip1, rake1, str2, dip2, rake2, ...rest) {
		let np1 = NodalPlane.fromParams(str1, dip1, rake1);
		let np2 = NodalPlane.fromParams(str2, dip2, rake2);
		let np = np1;
		let axes = np.getPrincipalAxes();
		return new FocalMechanism(np, axes);
	}

	static fromSp(str1, dip1, str2, fault, ...rest) {
		let np = NodalPlane.fromPartial(str1, dip1, str2, fault);
		let axes = np.getPrincipalAxes();
		return new FocalMechanism(np, axes);
	}

	static fromSm(mrr, mtt, mff, mrt, mrf, mtf, ...rest) {
		let a = [
			mrr, mrt, mrf,
			mrt, mtt, mtf,
			mrf, mtf, mff,
		];
		let d = new Array(3);
		let b = new Array(3);
		let z = new Array(3);
		let v = new Array(9);

		jacobi3(a, d, v, b, z);

		let axes = new Array(3);
		for (let j = 0; j < 3; j++) {
			// from (r, t, f) to (z, n, e)
			let axis = new Axis(-v[j*3], -v[j*3+1], v[j*3+2], d[j]);
			axes[j] = (axis.z >= 0) ? axis : axis.getOpposite();
		}
		let np = NodalPlane.fromAxes(axes[0], axes[2]);

		return new FocalMechanism(np, axes);
	}

	static fromSx(...params) {
		let axes = new Array(3);
		for (let j = 0; j < 3; j++) {
			let axis = Axis.fromPolar(params[j*3+1], params[j*3+2], params[j*3]);
			axes[j] = (axis.z >= 0) ? axis : axis.getOpposite();
		}
		let np = NodalPlane.fromAxes(axes[0], axes[2]);

		return new FocalMechanism(np, axes);
	}
}

class Beachball {
	static zne2proj(rz, rn, re) { // assume rz >= 0 (lower-hemisphere)
		let rh2 = rn * rn + re * re;
		if (rh2 == 0) return [0, 0];

		let r = Math.hypot(rz, rn, re);
		let f = 1 / Math.sqrt(r * (r + rz));
		return [re * f, rn * f];
	}

	static xy(n, rz, rn, re) {
		let arr = new Array(n * 2);
		for (let i = 0; i < n; i++) {
			let [x, y] = this.zne2proj(rz[i], rn[i], re[i]);
			arr[i*2+0] = round(x * scale);
			arr[i*2+1] = round(y * scale);
		}
		return arr;
	}

	drawCircle(canvas, x0, y0, radius) {
		let left = (x0 - radius) * scale;
		let top  = (y0 - radius) * scale;
		return canvas.circle(radius * scale * 2).transform({ tx: left, ty: top });
	}

	drawLine(canvas, n, rz, rn, re) {
		return canvas.polyline(Beachball.xy(n, rz, rn, re));
	}

	drawPolygon(canvas, n, rz, rn, re) {
		return canvas.polygon(Beachball.xy(n, rz, rn, re));
	}

	drawTensor(canvas, axes, fillComp, fillExt, plotZeroTrace = false) {
		axes = Array.from(axes); // shallow copy
		let v = Array.from(axes, x => x.value);

		if (v[0] < v[1] || v[1] < v[2]) {
			throw "principal axes must satisfy T.value >= N.value >= P.value";
		}

		if (plotZeroTrace) {
			let vi = (v[0] + v[1] + v[2]) / 3; // isotropic
			for (let j = 0; j < 3; j++) v[j] -= vi;
		}

		if (v[2] >= 0) {
			// explosive
			this.drawCircle(canvas, 0, 0, 1).fill(fillComp);
			return;
		} else if (v[0] <= 0) {
			// implosive
			this.drawCircle(canvas, 0, 0, 1).fill(fillExt);
			return;
		}

		let v0, v1, v2;
		let fillBase, fillClosed;
		if (v[1] <= 0) {
			[v0, v1, v2] = [v[0], -v[1], -v[2]];
			[fillBase, fillClosed] = [fillExt, fillComp];
		} else {
			axes.reverse();
			[v0, v1, v2] = [-v[2], v[1], v[0]];
			[fillBase, fillClosed] = [fillComp, fillExt];
		}
		/* push-pull boundaries (r^{T} M r == 0) surround axes[0] */

		let q = new Array(9); // orthogonal matrix
		for (let j = 0; j < 3; j++) {
			let axis = axes[j];
			let { z, n, e } = (axis.z >= 0) ? axis : axis.getOpposite();
			[q[j], q[j+3], q[j+6]] = [z, n, e];
		}

		let pz = new Array(N * 2 + 1);
		let pn = new Array(N * 2 + 1);
		let pe = new Array(N * 2 + 1);
		for (let i = 0; i < N * 2; i++) {
			let [c, s] = polar2xy(1, i * 180 / N);

			let h2 = v1 * c * c + v2 * s * s;
			let r2 = v0 + h2;
			let [t, n, p] = [
				Math.sqrt(h2 / r2),
				Math.sqrt(v0 / r2) * c,
				Math.sqrt(v0 / r2) * s,
			];

			pz[i] = q[0] * t + q[1] * n + q[2] * p;
			pn[i] = q[3] * t + q[4] * n + q[5] * p;
			pe[i] = q[6] * t + q[7] * n + q[8] * p;
		}
		[pz[N*2], pn[N*2], pe[N*2]] = [pz[0], pn[0], pe[0]];

		// detect the start of upper-/lower-hemisphere
		let i1 = null, i2 = null;
		let prev = pz[0]; // assert that pz[0] >= 0
		for (let i = 1; i <= N * 2; i++) {
			if (pz[i] < 0) {
				if (prev >= 0) i1 = i;
			} else {
				if (prev < 0) i2 = i;
			}
			prev = pz[i];
		}

		if (i1 == null) {
			this.drawCircle(canvas, 0, 0, 1).fill(fillBase);
			this.drawPolygon(canvas, N * 2, pz, pn, pe).fill(fillClosed);
			return;
		}

		let intersection = (z1, n1, e1, z2, n2, e2) => {
			let t = (0 - z1) / (z2 - z1);
			let [n, e] = [n1 * (1 - t) + n2 * t, e1 * (1 - t) + e2 * t];
			let h = Math.hypot(n, e);
			return [n / h, e / h];
		};
		let [pn1, pe1] = intersection(pz[i1-1], pn[i1-1], pe[i1-1], pz[i1], pn[i1], pe[i1]);
		let [pn2, pe2] = intersection(pz[i2-1], pn[i2-1], pe[i2-1], pz[i2], pn[i2], pe[i2]);
		let theta0 = rad2deg(Math.atan2(pn1 * pe2 - pe1 * pn2, pn1 * pn2 + pe1 * pe2));
		let m = Math.ceil(Math.abs(theta0) / 180 * N);

		let rn = new Array(m + 1);
		let re = new Array(m + 1);
		for (let i = 0; i <= m; i++) {
			let [c, s] = polar2xy(1, theta0 * i / m);
			[rn[i], re[i]] = [pn1 * c - pe1 * s, pn1 * s + pe1 * c];
		}

		this.drawCircle(canvas, 0, 0, 1).fill(fillBase);

		let bz = new Array(N * 3 + 1);
		let bn = new Array(N * 3 + 1);
		let be = new Array(N * 3 + 1);
		let j;

		j = 0;
		for (let i = i2; i < N * 2; i++, j++) [bz[j], bn[j], be[j]] = [pz[i], pn[i], pe[i]];
		for (let i = 0 ; i < i1   ; i++, j++) [bz[j], bn[j], be[j]] = [pz[i], pn[i], pe[i]];
		for (let i = 0 ; i <= m   ; i++, j++) [bz[j], bn[j], be[j]] = [0,     rn[i], re[i]];
		this.drawPolygon(canvas, j, bz, bn, be).fill(fillClosed);

		j = 0;
		for (let i = i1; i < i2; i++, j++) [bz[j], bn[j], be[j]] = [-pz[i], -pn[i], -pe[i]];
		for (let i = m ; i >= 0; i--, j++) [bz[j], bn[j], be[j]] = [-0,     -rn[i], -re[i]];
		this.drawPolygon(canvas, j, bz, bn, be).fill(fillClosed);

		return;
	}

	drawPlane(canvas, normalVector, pen) {
		let { z, n, e } = normalVector;
		let h2 = n * n + e * e;
		if (h2 == 0) return;

		let [sz, sn, se] = [0, e, -n];
		let [tz, tn, te] = [n * se - e * sn, e * sz - z * se, z * sn - n * sz];

		let pz = new Array(N + 1);
		let pn = new Array(N + 1);
		let pe = new Array(N + 1);
		for (let i = 0; i < N + 1; i++) {
			let [c, s] = polar2xy(1, -i * 180 / N); // rake: from 0 to -180
			[pz[i], pn[i], pe[i]] = [sz * c + tz * s, sn * c + tn * s, se * c + te * s];
		}

		this.drawLine(canvas, N + 1, pz, pn, pe).attr(pen);
		return;
	}

	drawOutline(canvas, pen) {
		this.drawCircle(canvas, 0, 0, 1).attr(pen);
		return;
	}
}

function jacobi3(a, d, v, b, z) {
	/* For the original implementation, see gmt_jacobi() in src/gmt_vector.c */

	const MAX_SWEEPS = 50;

	/* Begin by initializing v, b, d, and z.  v = identity matrix,
	        b = d = diag(a), and z = 0:  */

	v.splice(0, 9, ...[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
	b.splice(0, 3, ...[a[0], a[4], a[8]]);
	d.splice(0, 3, ...[a[0], a[4], a[8]]);
	z.fill(0.0, 0, 3);

	/* End of initializations.  Set counters and begin:  */

	nrots = 0;
	nsweeps = 0;

	const iter = [ /* p, q, pq (== p+q*3), ri, rj */
		[0, 1, 3, 6, 7],
		[0, 2, 6, 3, 7],
		[1, 2, 7, 3, 6],
	];

	while (nsweeps < MAX_SWEEPS) {
		/* Sum off-diagonal elements of upper triangle.  */
		let sum = Math.abs(a[3]) + Math.abs(a[6]) + Math.abs(a[7]);

		/* Exit this loop (converged) when sum == 0.0  */
		if (sum == 0.0) break;

		/* If (nsweeps < 3) do only bigger elements;  else all  */
		let threshold = (nsweeps < 3) ? 0.2 * sum / 9 : 0.0;

		/* Now sweep whole upper triangle doing Givens rotations:  */
		iter.forEach(([p, q, pq, ri, rj]) => {
			if (a[pq] == 0.0) return; /* continue */

			let g = 100.0 * Math.abs(a[pq]);
			let gIsSmallTo = (x) => { let z = Math.abs(x); return z + g == z };

			/* After four sweeps, if g is small relative
			        to a(p,p) and a(q,q), skip the
			        rotation and set a(p,q) to zero.  */

			if ((nsweeps > 3) && gIsSmallTo(d[p]) && gIsSmallTo(d[q])) {
				a[pq] = 0.0;
				return; /* continue */
			}

			if (Math.abs(a[pq]) <= threshold) return; /* continue */

			let h = d[q] - d[p];
			let t;
			if (h == 0.0) t = 1.0;
			else if (gIsSmallTo(h)) t = a[pq] / h;
			else {
				let theta = 0.5 * h / a[pq];
				t = 1.0 / (Math.abs(theta) + Math.sqrt(1.0 + theta * theta));
				if (theta < 0.0) t = -t;
			}

			let c = 1.0 / Math.sqrt(1.0 + t * t);
			let s = t * c;
			let tau = s / (1.0 + c);

			h = t * a[pq];
			z[p] -= h; z[q] += h;
			d[p] -= h; d[q] += h;
			a[pq] = 0.0;

			let rotate = (x, i, j) => {
				let g = x[i];
				let h = x[j];
				x[i] = g - s * (h + g * tau);
				x[j] = h + s * (g - h * tau);
			};

			rotate(a, ri, rj);
			for (j = 0; j < 3; j++) rotate(v, j + p * 3, j + q * 3);

			nrots++;
		});

		/* End of one sweep of the upper triangle.  */

		nsweeps++;

		for (let p = 0; p < 3; p++) {
			b[p] += z[p];   /* Update the b copy of diagonal  */
			d[p] = b[p];    /* Replace d with b to reduce round-off error  */
			z[p] = 0.0;     /* Clear z.  */
		}
	}

	/* Get here via break when converged, or when nsweeps == MAX_SWEEPS.
	        Sort eigenvalues by insertion:  */

	for (let i = 0; i < 2; i++) {
		let k = i;
		let g = d[i];
		for (let j = i + 1; j < 3; j++) {  /* Find max location  */
			if (d[j] > g) [k, g] = [j, d[j]];
		}
		if (k == i) continue;

		/*  Need to swap value and vector  */
		[d[k], d[i]] = [d[i], g];
		let p = i * 3;
		let q = k * 3;
		for (let j = 0; j < 3; j++) {
			[v[j + p], v[j + q]] = [v[j + q], v[j + p]];
		}
	}

	/* Return 0 if converged; else print warning and return -1:  */

	if (nsweeps == MAX_SWEEPS) {
		console.log(`jacobi3 failed to converge in ${nsweeps} sweeps`);
		return -1;
	}
	return 0;
}


const DATA_FORMATS = {
	a: { code: "a", argc: 3, legend: "strike dip rake [magnitude ...]" },
	c: { code: "c", argc: 6, legend: "strike dip rake (of plane 1, 2) [mantissa exponent ...]" },
	m: { code: "m", argc: 6, legend: "mrr mtt mff mrt mrf mtf [exponent ...]" },
	p: { code: "p", argc: 4, legend: "strike1 dip1 strike2 \u00b11 [magnitude ...]" },
	x: { code: "x", argc: 9, legend: "value azimuth plunge (of T, N, P axis) [exponent ...]" },
};
let selectedFormat = DATA_FORMATS["a"];

/* Define events */
const elementDataFormat = document.querySelector("#data-format");
// const elementPlotMode = document.querySelector("#plot-mode");
const elementParams = document.querySelector("#params");

elementDataFormat.addEventListener("change", event => {
	let fmt = DATA_FORMATS[event.target.value];
	if (fmt == null) return;

	selectedFormat = fmt;
	document.querySelector("#legend").textContent = fmt.legend;
});

elementParams.addEventListener("change", event => {
	let args = event.target.value.trim().split(/\s+/);
	if (args.length < selectedFormat.argc) return;

	let params = args.slice(0, selectedFormat.argc).map(s => Number(s));
	let fm;
	switch (selectedFormat.code) {
		case "a":
			fm = FocalMechanism.fromSa(...params);
			break;
		case "c":
			fm = FocalMechanism.fromSc(...params);
			break;
		case "p":
			fm = FocalMechanism.fromSp(...params);
			break;
		case "m":
			fm = FocalMechanism.fromSm(...params);
			break;
		case "x":
			fm = FocalMechanism.fromSx(...params);
			break;
		default:
			return;
	}

	console.log(fm);

	group.clear();
	let b = new Beachball();
	b.drawTensor(group, fm.axes, fillG, fillE, false);
	b.drawPlane(group, fm.np.nv, penT);
	b.drawPlane(group, fm.np.sv, penT);
	b.drawOutline(group, penL);
});
