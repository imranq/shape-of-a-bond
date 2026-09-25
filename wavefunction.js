// CPU side of the wavefunction: Slater-type AOs, a numerical overlap matrix
// (so any MO or blend of MOs can be normalised) and |ψ|² point sampling.
//
// Every AO is a Slater-type orbital  φ(r) = |d|^k · (A + B·d) · e^(−ζ|d|),  d = r − center.
// s: k = n−1, B = 0.  p: k = n−2, A = 0, B along the axis.
// The GPU raymarcher evaluates the exact same expression (k is packed as zeta + 100k).

const AXES = { px: [1, 0, 0], py: [0, 1, 0], pz: [0, 0, 1] };
const CUTOFF = 18; // e^-18 ≈ 1.5e-8
const FACT = [1, 1, 2, 6, 24, 120, 720, 5040];
const radialNorm = (n, z) => Math.pow(2 * z, n + 0.5) / Math.sqrt(FACT[2 * n]);

// basis entries: { atom, t: 's' | 'px' | 'py' | 'pz', n, z (Å⁻¹) }
export function buildAOs(mol, pos) {
  return mol.basis.map((b) => {
    const [cx, cy, cz] = pos[b.atom];
    const N = radialNorm(b.n, b.z);
    if (b.t === 's') {
      return { cx, cy, cz, z: b.z, k: b.n - 1, A: N / Math.sqrt(4 * Math.PI), bx: 0, by: 0, bz: 0 };
    }
    const f = N * Math.sqrt(3 / (4 * Math.PI)), d = AXES[b.t];
    return { cx, cy, cz, z: b.z, k: b.n - 2, A: 0, bx: d[0] * f, by: d[1] * f, bz: d[2] * f };
  });
}

const rpow = (r, k) => (k === 0 ? 1 : k === 1 ? r : r * r);

export function evalAOs(aos, x, y, z, out) {
  for (let i = 0; i < aos.length; i++) {
    const a = aos[i];
    const dx = x - a.cx, dy = y - a.cy, dz = z - a.cz;
    const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
    const zr = a.z * r;
    out[i] = zr > CUTOFF ? 0 : rpow(r, a.k) * (a.A + a.bx * dx + a.by * dy + a.bz * dz) * Math.exp(-zr);
  }
}

export function psiAt(aos, w, x, y, z) {
  let s = 0;
  for (let i = 0; i < aos.length; i++) {
    if (w[i] === 0) continue;
    const a = aos[i];
    const dx = x - a.cx, dy = y - a.cy, dz = z - a.cz;
    const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
    const zr = a.z * r;
    if (zr > CUTOFF) continue;
    s += w[i] * rpow(r, a.k) * (a.A + a.bx * dx + a.by * dy + a.bz * dz) * Math.exp(-zr);
  }
  return s;
}

export function bounds(pos, margin) {
  const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
  for (const p of pos) for (let k = 0; k < 3; k++) {
    min[k] = Math.min(min[k], p[k] - margin);
    max[k] = Math.max(max[k], p[k] + margin);
  }
  return { min, max };
}

// S_ij = ∫ φ_i φ_j dV on a uniform grid (midpoint rule; converges fast for these smooth decaying functions)
export function overlapMatrix(aos, pos, h = 0.1, margin = 2.6) {
  const n = aos.length, S = new Float64Array(n * n), v = new Float64Array(n);
  const { min, max } = bounds(pos, margin);
  for (let x = min[0] + h / 2; x < max[0]; x += h)
    for (let y = min[1] + h / 2; y < max[1]; y += h)
      for (let z = min[2] + h / 2; z < max[2]; z += h) {
        evalAOs(aos, x, y, z, v);
        for (let i = 0; i < n; i++) {
          const vi = v[i];
          if (vi === 0) continue;
          for (let j = i; j < n; j++) S[i * n + j] += vi * v[j];
        }
      }
  const h3 = h * h * h;
  for (let i = 0; i < n; i++) for (let j = i; j < n; j++) {
    S[i * n + j] *= h3;
    S[j * n + i] = S[i * n + j];
  }
  return S;
}

export function normOf(c, S) {
  const n = c.length;
  let s = 0;
  for (let i = 0; i < n; i++) {
    if (c[i] === 0) continue;
    for (let j = 0; j < n; j++) s += c[i] * c[j] * S[i * n + j];
  }
  return Math.sqrt(Math.max(s, 1e-12));
}

// Rejection-sample points from |ψ|². Each point is stored relative to its nearest nucleus,
// so the cloud can follow the atoms while the bond stretches without resampling.
export function sampleCloud(aos, w, pos, count, margin = 2.4) {
  const { min, max } = bounds(pos, margin);
  const span = [max[0] - min[0], max[1] - min[1], max[2] - min[2]];
  const rnd = (k) => min[k] + Math.random() * span[k];

  let peak = 0;
  for (let i = 0; i < 8000; i++) {
    const v = psiAt(aos, w, rnd(0), rnd(1), rnd(2));
    peak = Math.max(peak, v * v);
  }
  const probes = [0, 0.12, 0.25, 0.4, 0.6];
  for (const p of pos) for (const d of probes) for (let k = 0; k < 3; k++) for (const s of [-1, 1]) {
    const q = [p[0], p[1], p[2]];
    q[k] += s * d;
    const v = psiAt(aos, w, q[0], q[1], q[2]);
    peak = Math.max(peak, v * v);
  }
  peak *= 1.05;

  const offsets = new Float32Array(count * 3), near = new Uint8Array(count), sign = new Float32Array(count);
  let k = 0, tries = 0, best = -Infinity, worst = Infinity, iMax = -1, iMin = -1;
  while (k < count && tries < 2e6) {
    tries++;
    const x = rnd(0), y = rnd(1), z = rnd(2);
    const v = psiAt(aos, w, x, y, z);
    if (Math.random() * peak >= v * v) continue;
    let ai = 0, bd = Infinity;
    for (let a = 0; a < pos.length; a++) {
      const d = (x - pos[a][0]) ** 2 + (y - pos[a][1]) ** 2 + (z - pos[a][2]) ** 2;
      if (d < bd) { bd = d; ai = a; }
    }
    near[k] = ai;
    offsets[k * 3] = x - pos[ai][0];
    offsets[k * 3 + 1] = y - pos[ai][1];
    offsets[k * 3 + 2] = z - pos[ai][2];
    sign[k] = v >= 0 ? 1 : -1;
    // label anchor: strong |ψ| but away from any nucleus, i.e. the body of a lobe
    const score = v * Math.min(Math.sqrt(bd), 0.9);
    if (score > best) { best = score; iMax = k; }
    if (score < worst) { worst = score; iMin = k; }
    k++;
  }
  return { count: k, offsets, near, sign, iMax: best > 0 ? iMax : -1, iMin: worst < 0 ? iMin : -1 };
}
