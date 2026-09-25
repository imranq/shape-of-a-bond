// Resonance structures: the ways to place π bonds (and any charge) over a conjugated
// system without moving atoms. Covers Kekulé structures (benzene, pyridine, naphthalene)
// and charge-delocalised ions (carboxylate, allyl, nitrate-like). Each structure is a
// list of per-bond orders plus per-atom formal charges.
import { EL } from './elements.js';
import { neighbors } from './smiles.js';

const NEUTRAL_VALENCE = { B: 3, C: 4, N: 3, O: 2, P: 3, S: 2, F: 1, Cl: 1, H: 1 };

export function resonanceStructures(g, limit = 8) {
  const nb = neighbors(g);
  const N = g.atoms.length;
  const sigma = (i) => nb[i].length;
  const hasPi = (i) => nb[i].some(({ k }) => g.bonds[k].order >= 1.5 && g.bonds[k].order < 3);
  const inTriple = (i) => nb[i].some(({ k }) => g.bonds[k].order === 3);

  // atoms that could carry one π bond: neutral valence leaves room for exactly one more bond
  const capable = (i) => {
    const a = g.atoms[i];
    if (a.el === 'H' || inTriple(i)) return false;
    return (NEUTRAL_VALENCE[a.el] ?? 0) - sigma(i) >= 1;
  };
  const pi = new Set();
  for (let i = 0; i < N; i++) if (hasPi(i) && capable(i)) pi.add(i);
  // charged atoms next to the π system join it (their charge can move)
  for (let i = 0; i < N; i++) {
    if (!g.atoms[i].charge || pi.has(i) || !capable(i)) continue;
    if (nb[i].some(({ j }) => pi.has(j))) pi.add(i);
  }
  if (pi.size < 3) return null;

  const piAtoms = [...pi];
  const edges = g.bonds.map((b, k) => ({ ...b, k })).filter((b) => pi.has(b.a) && pi.has(b.b));
  const charged = piAtoms.filter((i) => g.atoms[i].charge);
  const q = charged.reduce((s, i) => s + g.atoms[i].charge, 0);
  const kBonds = Math.floor((piAtoms.length - charged.length) / 2);
  if (kBonds < 1) return null;

  // enumerate matchings of size kBonds (backtracking over edges)
  const results = [];
  const used = new Set(), chosen = [];
  const rec = (start) => {
    if (results.length >= 64) return;
    if (chosen.length === kBonds) { results.push([...chosen]); return; }
    for (let e = start; e < edges.length; e++) {
      const { a, b } = edges[e];
      if (used.has(a) || used.has(b)) continue;
      used.add(a); used.add(b); chosen.push(e);
      rec(e + 1);
      used.delete(a); used.delete(b); chosen.pop();
    }
  };
  rec(0);

  const structures = [];
  const seen = new Set();
  for (const m of results) {
    const matched = new Set(m.flatMap((e) => [edges[e].a, edges[e].b]));
    const unmatched = piAtoms.filter((i) => !matched.has(i));
    // unmatched atoms must be able to hold the charge (or the radical)
    if (charged.length && unmatched.length !== charged.length) continue;
    if (q > 0 && unmatched.some((i) => g.atoms[i].el === 'O' || g.atoms[i].el === 'F')) continue;
    const orders = g.bonds.map((b) => b.order);
    for (const e of edges) orders[e.k] = 1;
    for (const e of m) orders[edges[e].k] = 2;
    const charges = g.atoms.map((a, i) => (pi.has(i) ? 0 : a.charge));
    const per = charged.length ? q / charged.length : 0;
    for (const i of unmatched) charges[i] = per;
    const key = orders.join(',') + '|' + charges.join(',');
    if (seen.has(key)) continue;
    seen.add(key);
    structures.push({ orders, charges });
    if (structures.length >= limit) break;
  }
  if (structures.length < 2) return null;
  // hybrid: average π order over structures
  const hybrid = g.bonds.map((b, k) => structures.reduce((s, st) => s + st.orders[k], 0) / structures.length);
  const hCharges = g.atoms.map((_, i) => structures.reduce((s, st) => s + st.charges[i], 0) / structures.length);
  return { structures, hybrid, hybridCharges: hCharges, piAtoms };
}
