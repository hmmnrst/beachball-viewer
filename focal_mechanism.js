/**
 * Manage a focal mechanism.
 * Interpret various data formats suppoted by GMT.
 *
 * Copyright: 2021 Masahiro Nomoto
 * This code is licensed under GNU LGPL.
 */

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
