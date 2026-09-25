// 3D structure from a molecular graph.
// VSEPR-style model: every atom's bonds *and lone pairs* are electron domains that
// spread to ideal angles (2 → linear, 3 → trigonal, 4 → tetrahedral). Lone pairs are
// ghost particles. A small force field (bond lengths, 1-3 distances, steric repulsion,
// π-planarity) is relaxed with FIRE, first in 4D to escape tangles and then squeezed to 3D.
import { EL } from './elements.js';
import { neighbors, lonePairs } from './smiles.js';

const IDEAL = { 2: Math.PI, 3: (2 * Math.PI) / 3, 4: Math.acos(-1 / 3) };

export function bondLength(ea, eb, order, domA, domB) {
  const shrink = (e, d) => (e === 'H' ? 0 : d === 2 ? 0.07 : d === 3 ? 0.03 : 0);
  if (ea === 'H' && eb === 'H') return 0.74;
  let r = EL[ea].r + EL[eb].r - shrink(ea, domA) - shrink(eb, domB);
  if (order === 1.5) r -= 0.06;
  else if (order === 2) r -= 0.12;
  else if (order === 3) r -= 0.18;
  return r;
}

// opts: { long: Set(bondIndex) stretched ×1.35, linear: [[i, c, j]], init: [[x,y,z]], iters }
export function buildModel(g, opts = {}) {
  const nb = neighbors(g);
  const N = g.atoms.length;
  const lp = g.atoms.map((_, i) => (opts.lpOverride ? opts.lpOverride[i] : lonePairs(g, i, nb)));
  const dom = g.atoms.map((_, i) => nb[i].length + lp[i]);

  // particles: atoms then lone-pair ghosts
  const owner = [];
  const pairs = []; // [i, j, d0, k, kind] kind: 0 spring, 1 min-distance
  g.bonds.forEach((b, k) => {
    let d = bondLength(g.atoms[b.a].el, g.atoms[b.b].el, b.order, dom[b.a], dom[b.b]);
    if (opts.long?.has(k)) d *= 1.35;
    pairs.push([b.a, b.b, d, 1.0, 0]);
  });
  const bondLen = new Map();
  for (const p of pairs) { bondLen.set(`${p[0]},${p[1]}`, p[2]); bondLen.set(`${p[1]},${p[0]}`, p[2]); }

  let M = N;
  const domains = g.atoms.map(() => []);
  for (let i = 0; i < N; i++) {
    for (const { j } of nb[i]) domains[i].push({ p: j, d: bondLen.get(`${i},${j}`) });
    for (let k = 0; k < lp[i]; k++) {
      owner[M] = i;
      pairs.push([i, M, 0.55, 1.0, 0]);
      domains[i].push({ p: M, d: 0.55 });
      M++;
    }
  }

  const bonded = new Set();
  const key = (a, b) => (a < b ? `${a},${b}` : `${b},${a}`);
  for (const p of pairs) bonded.add(key(p[0], p[1]));
  const linear = new Set((opts.linear || []).map(([i, , j]) => key(i, j)));

  // 1-3 terms around each centre
  for (let c = 0; c < N; c++) {
    const D = domains[c];
    const theta = IDEAL[D.length];
    for (let u = 0; u < D.length; u++) for (let v = u + 1; v < D.length; v++) {
      const a = D[u], b = D[v];
      const kk = key(a.p, b.p);
      bonded.add(kk);
      if (linear.has(kk)) { pairs.push([a.p, b.p, a.d + b.d, 0.6, 0]); continue; }
      if (theta) {
        const d = Math.sqrt(a.d * a.d + b.d * b.d - 2 * a.d * b.d * Math.cos(theta));
        pairs.push([a.p, b.p, d, 0.5, 0]);
      } else {
        // 5+ domains: just keep them at least ~90° apart
        pairs.push([a.p, b.p, Math.SQRT2 * Math.min(a.d, b.d) * 0.98, 0.5, 1]);
      }
    }
  }
  // steric repulsion between everything else (atoms only)
  for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
    if (bonded.has(key(i, j))) continue;
    const hi = g.atoms[i].el === 'H', hj = g.atoms[j].el === 'H';
    const dmin = hi && hj ? 2.3 : hi || hj ? 2.5 : 2.85;
    pairs.push([i, j, dmin, 0.25, 1]);
  }

  // π planarity: across bonds joining two trigonal centres, all substituents share a plane
  const quads = [];
  const trig = (i) => dom[i] === 3 && g.atoms[i].el !== 'H';
  g.bonds.forEach((b) => {
    if (!(trig(b.a) && trig(b.b))) return;
    for (const { j: i } of nb[b.a]) if (i !== b.b) for (const { j: l } of nb[b.b]) if (l !== b.a) quads.push([i, b.a, b.b, l]);
  });

  const energyGrad = (X, G, dim, w4) => {
    G.fill(0);
    let E = 0;
    for (const [i, j, d0, k, kind] of pairs) {
      let d2 = 0;
      for (let c = 0; c < dim; c++) { const t = X[i * 4 + c] - X[j * 4 + c]; d2 += t * t; }
      const d = Math.sqrt(d2) + 1e-9;
      const diff = d - d0;
      if (kind === 1 && diff > 0) continue;
      E += k * diff * diff;
      const f = (2 * k * diff) / d;
      for (let c = 0; c < dim; c++) {
        const t = (X[i * 4 + c] - X[j * 4 + c]) * f;
        G[i * 4 + c] += t;
        G[j * 4 + c] -= t;
      }
    }
    if (dim === 4 && w4 > 0) for (let i = 0; i < M; i++) { const w = X[i * 4 + 3]; E += w4 * w * w; G[i * 4 + 3] += 2 * w4 * w; }
    if (dim === 3) {
      const kp = 0.6;
      for (const [a, b, c, d] of quads) {
        const A = [0, 1, 2].map((t) => X[a * 4 + t]), B = [0, 1, 2].map((t) => X[b * 4 + t]);
        const C = [0, 1, 2].map((t) => X[c * 4 + t]), Dd = [0, 1, 2].map((t) => X[d * 4 + t]);
        const u = sub(B, A), v = sub(C, A), w = sub(Dd, A);
        const vol = dot(u, cross(v, w));
        E += kp * vol * vol;
        const s = 2 * kp * vol;
        const gb = cross(v, w), gc = cross(w, u), gd = cross(u, v);
        for (let t = 0; t < 3; t++) {
          G[b * 4 + t] += s * gb[t];
          G[c * 4 + t] += s * gc[t];
          G[d * 4 + t] += s * gd[t];
          G[a * 4 + t] -= s * (gb[t] + gc[t] + gd[t]);
        }
      }
    }
    return E;
  };

  const run = (X, dim, iters, w4sched) => {
    const V = new Float64Array(M * 4), G = new Float64Array(M * 4);
    let dt = 0.03, alpha = 0.1, npos = 0, E = 0;
    for (let it = 0; it < iters; it++) {
      const w4 = w4sched ? w4sched(it / iters) : 0;
      E = energyGrad(X, G, dim, w4);
      let P = 0, fn = 0, vn = 0;
      for (let q = 0; q < M * 4; q++) { P -= G[q] * V[q]; fn += G[q] * G[q]; vn += V[q] * V[q]; }
      fn = Math.sqrt(fn) + 1e-12; vn = Math.sqrt(vn);
      if (P > 0) {
        for (let q = 0; q < M * 4; q++) V[q] = (1 - alpha) * V[q] - (alpha * vn * G[q]) / fn;
        if (++npos > 5) { dt = Math.min(dt * 1.1, 0.12); alpha *= 0.99; }
      } else { V.fill(0); dt *= 0.5; alpha = 0.1; npos = 0; }
      for (let q = 0; q < M * 4; q++) {
        V[q] -= G[q] * dt;
        let step = V[q] * dt;
        if (step > 0.2) step = 0.2; else if (step < -0.2) step = -0.2;
        X[q] += step;
      }
      if (dim === 3) for (let i = 0; i < M; i++) X[i * 4 + 3] = 0;
    }
    return E;
  };

  let best = null, bestE = Infinity;
  const tries = opts.init ? 1 : N < 30 ? 4 : 2;
  const iters = opts.iters || (N < 40 ? 1600 : 1000);
  for (let t = 0; t < tries; t++) {
    const X = new Float64Array(M * 4);
    const box = 1.2 * Math.cbrt(M) + 1;
    for (let i = 0; i < M; i++) {
      const base = opts.init ? (i < N ? opts.init[i] : opts.init[owner[i]]) : null;
      for (let c = 0; c < 4; c++) {
        X[i * 4 + c] = base
          ? (c < 3 ? base[c] + (i >= N ? (Math.random() - 0.5) * 0.4 : 0) : 0)
          : (Math.random() - 0.5) * box * (c === 3 ? 0.5 : 1);
      }
    }
    if (!opts.init) run(X, 4, iters, (f) => (f < 0.4 ? 0 : 3 * (f - 0.4) / 0.6));
    const E = run(X, 3, Math.round(iters * (opts.init ? 0.6 : 0.7)));
    if (E < bestE) { bestE = E; best = X; }
  }
  const pos = [];
  for (let i = 0; i < N; i++) pos.push([best[i * 4], best[i * 4 + 1], best[i * 4 + 2]]);
  return { pos, lp, dom, energy: bestE };
}

const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];

// Centre and rotate so the widest spread lies along x, the flattest along y
// (planar molecules lie flat on the stage with their π lobes pointing up).
export function orient(pos, weights) {
  const n = pos.length;
  const w = weights || pos.map(() => 1);
  const W = w.reduce((a, b) => a + b, 0);
  const c = [0, 0, 0];
  pos.forEach((p, i) => { for (let k = 0; k < 3; k++) c[k] += (p[k] * w[i]) / W; });
  const q = pos.map((p) => sub(p, c));
  const cov = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  q.forEach((p, i) => { for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) cov[a][b] += p[a] * p[b] * w[i]; });
  const { values, vectors } = jacobi3(cov);
  const idx = [0, 1, 2].sort((a, b) => values[b] - values[a]); // big, mid, small
  const ax = vectors[idx[0]], az = vectors[idx[1]];
  let ay = cross(az, ax);
  const out = q.map((p) => [dot(p, ax), dot(p, ay), dot(p, az)]);
  return n ? out : [];
}

function jacobi3(A) {
  const a = A.map((r) => r.slice());
  const v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (let sweep = 0; sweep < 30; sweep++) {
    for (let p = 0; p < 2; p++) for (let r = p + 1; r < 3; r++) {
      if (Math.abs(a[p][r]) < 1e-12) continue;
      const th = 0.5 * Math.atan2(2 * a[p][r], a[r][r] - a[p][p]);
      const c = Math.cos(th), s = Math.sin(th);
      for (let k = 0; k < 3; k++) {
        const x = a[k][p], y = a[k][r];
        a[k][p] = c * x - s * y; a[k][r] = s * x + c * y;
      }
      for (let k = 0; k < 3; k++) {
        const x = a[p][k], y = a[r][k];
        a[p][k] = c * x - s * y; a[r][k] = s * x + c * y;
      }
      for (let k = 0; k < 3; k++) {
        const x = v[k][p], y = v[k][r];
        v[k][p] = c * x - s * y; v[k][r] = s * x + c * y;
      }
    }
  }
  return { values: [a[0][0], a[1][1], a[2][2]], vectors: [0, 1, 2].map((k) => [v[0][k], v[1][k], v[2][k]]) };
}
