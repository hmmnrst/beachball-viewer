/**
 * Manage a nodal plane of a focal mechanism, i.e. a fault (and its slip direction).
 *
 * This code is licensed under CC0.
 */

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

		let [tAxis, pAxis] = Axis.sumsub(nv, sv, +1, -1);
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
