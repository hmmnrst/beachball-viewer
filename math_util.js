/**
 * Math utility functions.
 *
 * This code is licensed under CC0.
 */

function deg2rad(theta) { return theta * (Math.PI / 180); }
function rad2deg(theta) { return theta / (Math.PI / 180); }

function polar2xy(radius, theta) {
	if (radius == 0) return [0, 0];
	theta = deg2rad(theta);
	return [radius * Math.cos(theta), radius * Math.sin(theta)];
}

function polar2zne(radius, azimuth, plunge) {
	if (radius == 0) return [0, 0, 0];
	let [rh, rz] = polar2xy(radius, plunge);
	let [rn, re] = polar2xy(rh, azimuth);
	return [rz, rn, re]; // z (depth), northward, eastward
}

function zne2proj(rz, rn, re) { // assume rz >= 0 (lower-hemisphere)
	let rh2 = rn * rn + re * re;
	if (rh2 == 0) return [0, 0];

	let r = Math.hypot(rz, rn, re);
	// let f = 1 / Math.sqrt(r * (r + rz));
	let f = 1 / r;
	return [re * f, rn * f];
}
