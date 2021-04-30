/**
 * Draw a focal mechanism.
 *
 * Copyright (C) 2021 hmmnrst
 * This code is licensed under GNU LGPL.
 */

const scale = 100; // radius of a beachball
function round(x) { return Math.round(x * 100) / 100; }

const N = 12; // number of points per 180 degrees

class Beachball {
	static xy(n, rz, rn, re) {
		let arr = new Array(n * 2);
		for (let i = 0; i < n; i++) {
			let [x, y] = zne2proj(rz[i], rn[i], re[i]);
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
