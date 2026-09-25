// Extended Hückel theory (Hoffmann 1963).
//   Basis: valence Slater orbitals; overlaps from 3-Gaussian fits of each Slater function
//   (fitted offline, overlap with the exact STO > 0.9998).
//   H_ii = valence-state ionisation energy,  H_ij = K·S_ij·(H_ii + H_jj)/2,  K = 1.75
//   Solve HC = SCE by Löwdin orthogonalisation.
import { EL, BOHR } from './elements.js';

// ζ = 1 fits: [exponents], [coefficients of normalised primitives]
const FIT = {
  '1s': [[2.2276608, 0.4057712, 0.1098175], [0.154329, 0.5353281, 0.4446346]],
  '2s': [[2.5815781, 0.1567622, 0.0601833], [-0.0599448, 0.5960384, 0.4581788]],
  '2p': [[0.9192377, 0.2359194, 0.0800981], [0.1623949, 0.5661709, 0.4223071]],
  '3s': [[0.5641488, 0.0692442, 0.0326953], [-0.1782578, 0.8612763, 0.2261841]],
  '3p': [[2.6928778, 0.148936, 0.0573958], [-0.0106195, 0.5218564, 0.5450015]],
};
const K_WH = 1.75;
const AX = { px: 0, py: 1, pz: 2 };

// basis entry: { atom, t: 's'|'px'|'py'|'pz', n, z (Å⁻¹), zAU, H }
export function valenceBasis(atoms) {
  const basis = [];
  atoms.forEach((a, i) => {
    const e = EL[a.el];
    const zA = e.zeta / BOHR;
    basis.push({ atom: i, t: 's', n: e.n, z: zA, zAU: e.zeta, H: e.hs });
    if (a.el !== 'H') for (const t of ['px', 'py', 'pz']) basis.push({ atom: i, t, n: e.n, z: zA, zAU: e.zeta, H: e.hp });
  });
  return basis;
}

function primitives(b) {
  const key = b.n + (b.t === 's' ? 's' : 'p');
  const [al, co] = FIT[key];
  const z2 = b.zAU * b.zAU;
  return al.map((a0, k) => {
    const a = a0 * z2;
    const ns = Math.pow((2 * a) / Math.PI, 0.75);
    return { a, c: co[k] * (b.t === 's' ? ns : ns * 2 * Math.sqrt(a)) };
  });
}

// analytic overlap of contracted s/p Cartesian Gaussians (positions in bohr)
export function overlapMatrix(basis, posA) {
  const n = basis.length;
  const pos = posA.map((p) => p.map((x) => x / BOHR));
  const prims = basis.map(primitives);
  const S = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    const bi = basis[i], A = pos[bi.atom], li = AX[bi.t];
    for (let j = i; j < n; j++) {
      const bj = basis[j], B = pos[bj.atom], lj = AX[bj.t];
      const AB2 = (A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2;
      if (AB2 > 400) continue; // > ~10 Å
      let s = 0;
      for (const p of prims[i]) for (const q of prims[j]) {
        const g = p.a + q.a, mu = (p.a * q.a) / g;
        const ss = Math.pow(Math.PI / g, 1.5) * Math.exp(-mu * AB2);
        let v;
        if (li === undefined && lj === undefined) v = ss;
        else {
          const P = [0, 1, 2].map((k) => (p.a * A[k] + q.a * B[k]) / g);
          const PA = li === undefined ? 0 : P[li] - A[li];
          const PB = lj === undefined ? 0 : P[lj] - B[lj];
          if (li === undefined) v = PB * ss;
          else if (lj === undefined) v = PA * ss;
          else v = (PA * PB + (li === lj ? 1 / (2 * g) : 0)) * ss;
        }
        s += p.c * q.c * v;
      }
      S[i * n + j] = S[j * n + i] = s;
    }
  }
  // exact normalisation on the diagonal
  const d = new Float64Array(n);
  for (let i = 0; i < n; i++) d[i] = 1 / Math.sqrt(S[i * n + i]);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) S[i * n + j] *= d[i] * d[j];
  return S;
}

// symmetric eigen-decomposition (cyclic Jacobi). Returns ascending values, vectors as columns in V (row-major n×n)
export function eigh(Ain, n) {
  const A = Float64Array.from(Ain), V = new Float64Array(n * n);
  for (let i = 0; i < n; i++) V[i * n + i] = 1;
  for (let sweep = 0; sweep < 60; sweep++) {
    let off = 0;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += A[p * n + q] ** 2;
    if (off < 1e-20) break;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
      const apq = A[p * n + q];
      if (Math.abs(apq) < 1e-14) continue;
      const th = (A[q * n + q] - A[p * n + p]) / (2 * apq);
      const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < n; k++) {
        const akp = A[k * n + p], akq = A[k * n + q];
        A[k * n + p] = c * akp - s * akq;
        A[k * n + q] = s * akp + c * akq;
      }
      for (let k = 0; k < n; k++) {
        const apk = A[p * n + k], aqk = A[q * n + k];
        A[p * n + k] = c * apk - s * aqk;
        A[q * n + k] = s * apk + c * aqk;
      }
      for (let k = 0; k < n; k++) {
        const vkp = V[k * n + p], vkq = V[k * n + q];
        V[k * n + p] = c * vkp - s * vkq;
        V[k * n + q] = s * vkp + c * vkq;
      }
    }
  }
  const order = [...Array(n).keys()].sort((a, b) => A[a * n + a] - A[b * n + b]);
  const values = order.map((k) => A[k * n + k]);
  const vecs = order.map((k) => { const v = new Float64Array(n); for (let r = 0; r < n; r++) v[r] = V[r * n + k]; return v; });
  return { values, vecs };
}

export function solveEHT(atoms, pos, charge = 0, basis = valenceBasis(atoms)) {
  const n = basis.length;
  const S = overlapMatrix(basis, pos);
  const H = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    H[i * n + j] = i === j ? basis[i].H : (K_WH * S[i * n + j] * (basis[i].H + basis[j].H)) / 2;
  }
  // X = S^-1/2
  const { values: sv, vecs: U } = eigh(S, n);
  const X = new Float64Array(n * n);
  for (let k = 0; k < n; k++) {
    const f = 1 / Math.sqrt(Math.max(sv[k], 1e-8));
    const u = U[k];
    for (let i = 0; i < n; i++) { const ui = u[i] * f; if (ui) for (let j = 0; j < n; j++) X[i * n + j] += ui * u[j]; }
  }
  const T = mul(X, H, n), Hp = mul(T, X, n);
  const { values: E, vecs: Vp } = eigh(Hp, n);
  const C = Vp.map((v) => {
    const c = new Float64Array(n);
    for (let i = 0; i < n; i++) { let s = 0; for (let j = 0; j < n; j++) s += X[i * n + j] * v[j]; c[i] = s; }
    // deterministic phase: largest coefficient positive
    let m = 0;
    for (let i = 1; i < n; i++) if (Math.abs(c[i]) > Math.abs(c[m])) m = i;
    if (c[m] < 0) for (let i = 0; i < n; i++) c[i] = -c[i];
    return c;
  });
  const nElec = atoms.reduce((s, a) => s + EL[a.el].ve, 0) - charge;
  const occ = E.map((_, k) => Math.max(0, Math.min(2, nElec - 2 * k)));
  const homo = Math.ceil(nElec / 2) - 1;
  return { basis, S, E, C, occ, nElec, homo, lumo: homo + 1 };
}

function mul(A, B, n) {
  const C = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let k = 0; k < n; k++) {
    const a = A[i * n + k];
    if (!a) continue;
    for (let j = 0; j < n; j++) C[i * n + j] += a * B[k * n + j];
  }
  return C;
}

// ---------- orbital analysis for the plain-language descriptions ----------
export function analyseMO(res, k, atoms, bonds, opts = {}) {
  const { basis, S, C } = res;
  const n = basis.length, c = C[k];
  const share = new Float64Array(atoms.length);
  const pShare = { x: 0, y: 0, z: 0, s: 0 };
  for (let i = 0; i < n; i++) {
    let g = 0;
    for (let j = 0; j < n; j++) g += c[j] * S[i * n + j];
    const q = c[i] * g;
    share[basis[i].atom] += q;
    pShare[basis[i].t === 's' ? 's' : basis[i].t[1]] += q;
  }
  const bondPop = bonds.map((b) => {
    let p = 0;
    for (let i = 0; i < n; i++) if (basis[i].atom === b.a) for (let j = 0; j < n; j++) if (basis[j].atom === b.b) p += 2 * c[i] * c[j] * S[i * n + j];
    return p;
  });
  return { share, bondPop, pShare };
}
