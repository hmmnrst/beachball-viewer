/**
 * Manage a principal axis (or any vector).
 * Optionally, `value` can represent a principal value.
 *
 * This code is licensed under CC0.
 */

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

	static sumsub(a1, a2, v1 = undefined, v2 = undefined) {
		// assume that a1 and a2 are orthonormal
		return [
			new Axis(
				(a1.z + a2.z) * Math.SQRT1_2,
				(a1.n + a2.n) * Math.SQRT1_2,
				(a1.e + a2.e) * Math.SQRT1_2,
				v1
			),
			new Axis(
				(a1.z - a2.z) * Math.SQRT1_2,
				(a1.n - a2.n) * Math.SQRT1_2,
				(a1.e - a2.e) * Math.SQRT1_2,
				v2
			)
		];
	}
}
