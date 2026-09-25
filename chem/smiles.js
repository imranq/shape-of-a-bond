// A compact SMILES reader: organic subset, bracket atoms (charge, H count, atom maps),
// branches, ring closures (incl. %nn), aromatic atoms, dot-separated fragments.
// Stereo marks (@, /, \) are accepted and ignored.
import { EL } from './elements.js';

const ORGANIC = ['Cl', 'Br', 'B', 'C', 'N', 'O', 'P', 'S', 'F', 'I'];
const AROMATIC = ['b', 'c', 'n', 'o', 'p', 's'];
const BRACKET = /^(\d+)?([A-Z][a-z]?|[a-z]{1,2})(@*)(?:H(\d*))?([+-]+\d*)?(?::(\d+))?$/;

function parseBracket(body) {
  const m = BRACKET.exec(body);
  if (!m) throw new Error(`Can't read bracket atom [${body}]`);
  const sym = m[2];
  const aromatic = sym[0] === sym[0].toLowerCase();
  const el = aromatic ? sym[0].toUpperCase() + sym.slice(1) : sym;
  let charge = 0;
  if (m[5]) {
    const sign = m[5][0] === '+' ? 1 : -1;
    const digits = m[5].replace(/[+-]/g, '');
    charge = sign * (digits ? parseInt(digits, 10) : m[5].length);
  }
  return {
    el, aromatic, charge,
    hcount: m[4] === undefined ? 0 : m[4] === '' ? 1 : parseInt(m[4], 10),
    map: m[6] ? parseInt(m[6], 10) : 0,
  };
}

export function parseSmiles(str) {
  str = str.trim();
  if (!str) throw new Error('Type a SMILES string, e.g. CCO');
  const atoms = [], bonds = [], stack = [], rings = new Map();
  let prev = -1, pending = null, i = 0;

  const addBond = (a, b, sym) => {
    if (a === b) throw new Error('Ring closure onto the same atom');
    if (bonds.some((x) => (x.a === a && x.b === b) || (x.a === b && x.b === a))) throw new Error('Duplicate bond');
    let order = 1, arom = false;
    if (sym === '=') order = 2;
    else if (sym === '#') order = 3;
    else if (sym === ':') { order = 1.5; arom = true; }
    else if (!sym || sym === '/' || sym === '\\') {
      if (atoms[a].aromatic && atoms[b].aromatic) { order = 1.5; arom = true; }
    }
    bonds.push({ a, b, order, arom });
  };
  const addAtom = (a) => {
    atoms.push(a);
    const idx = atoms.length - 1;
    if (prev >= 0) addBond(prev, idx, pending);
    pending = null;
    prev = idx;
  };

  while (i < str.length) {
    const ch = str[i];
    if (ch === ' ' || ch === '\t') break;
    if (ch === '(') { if (prev < 0) throw new Error('Branch with no atom before it'); stack.push(prev); i++; continue; }
    if (ch === ')') { if (!stack.length) throw new Error('Unmatched )'); prev = stack.pop(); i++; continue; }
    if ('-=#:/\\'.includes(ch)) { pending = ch; i++; continue; }
    if (ch === '.') { prev = -1; pending = null; i++; continue; }
    if (/\d/.test(ch) || ch === '%') {
      let key;
      if (ch === '%') { key = str.substr(i + 1, 2); i += 3; } else { key = ch; i++; }
      if (prev < 0) throw new Error('Ring number with no atom');
      if (rings.has(key)) {
        const r = rings.get(key);
        rings.delete(key);
        addBond(r.atom, prev, pending || r.bond);
      } else rings.set(key, { atom: prev, bond: pending });
      pending = null;
      continue;
    }
    if (ch === '[') {
      const j = str.indexOf(']', i);
      if (j < 0) throw new Error('Unclosed [');
      addAtom(parseBracket(str.slice(i + 1, j)));
      i = j + 1;
      continue;
    }
    const org = ORGANIC.find((e) => str.startsWith(e, i));
    if (org) { addAtom({ el: org, aromatic: false, charge: 0, hcount: null, map: 0 }); i += org.length; continue; }
    if (AROMATIC.includes(ch)) { addAtom({ el: ch.toUpperCase(), aromatic: true, charge: 0, hcount: null, map: 0 }); i++; continue; }
    throw new Error(`Unexpected "${ch}" at position ${i + 1}`);
  }
  if (rings.size) throw new Error(`Ring ${[...rings.keys()][0]} is never closed`);
  if (stack.length) throw new Error('Unclosed (');
  if (!atoms.length) throw new Error('No atoms found');
  for (const a of atoms) if (!EL[a.el]) throw new Error(`Element ${a.el} isn't supported (H B C N O F P S Cl)`);
  return { atoms, bonds };
}

// Adds implicit hydrogens as explicit atoms. Each H remembers its parent,
// and gets a mapping key when its parent carries an atom-map number.
export function addHydrogens(g) {
  const atoms = g.atoms.map((a) => ({ ...a })), bonds = g.bonds.map((b) => ({ ...b }));
  const heavy = atoms.length;
  for (let i = 0; i < heavy; i++) {
    const a = atoms[i];
    let h = a.hcount;
    if (h === null) {
      const mine = bonds.filter((b) => b.a === i || b.b === i);
      let sum = mine.reduce((s, b) => s + (b.arom ? 1 : b.order), 0);
      if (a.aromatic && !mine.some((b) => b.order === 2)) sum += 1;
      const base = EL[a.el].val;
      // aromatic atoms keep their lowest valence (pyrrole-type N needs an explicit [nH])
      let v = a.aromatic ? base[0] : base.find((x) => x >= sum) ?? base[base.length - 1];
      if (a.el === 'C') v -= Math.abs(a.charge);
      else if (a.el === 'B') v -= a.charge;
      else v += a.charge;
      h = Math.max(0, v - sum);
    }
    a.hcount = h;
    for (let k = 0; k < h; k++) {
      atoms.push({ el: 'H', aromatic: false, charge: 0, hcount: 0, map: 0, parent: i, hkey: a.map ? `${a.map}:${k}` : null });
      bonds.push({ a: i, b: atoms.length - 1, order: 1, arom: false });
    }
  }
  return { atoms, bonds };
}

export function neighbors(g) {
  const nb = g.atoms.map(() => []);
  g.bonds.forEach((b, k) => { nb[b.a].push({ j: b.b, k }); nb[b.b].push({ j: b.a, k }); });
  return nb;
}

// Lone pairs from the valence count; aromatic bonds count 1.5.
export function lonePairs(g, i, nb = neighbors(g)) {
  const a = g.atoms[i];
  const used = nb[i].reduce((s, x) => s + g.bonds[x.k].order, 0);
  return Math.max(0, Math.floor((EL[a.el].ve - a.charge - used) / 2 + 1e-6));
}

export function components(g) {
  const nb = neighbors(g), comp = new Array(g.atoms.length).fill(-1);
  let c = 0;
  for (let s = 0; s < g.atoms.length; s++) {
    if (comp[s] >= 0) continue;
    const q = [s];
    comp[s] = c;
    while (q.length) {
      const u = q.pop();
      for (const { j } of nb[u]) if (comp[j] < 0) { comp[j] = c; q.push(j); }
    }
    c++;
  }
  return { comp, count: c };
}

// Aromatic rings: for every aromatic bond, the shortest aromatic path back gives a ring.
export function aromaticRings(g) {
  const nb = neighbors(g), seen = new Set(), rings = [];
  g.bonds.forEach((b, k) => {
    if (!b.arom) return;
    const prev = new Map([[b.a, -1]]), q = [b.a];
    while (q.length) {
      const u = q.shift();
      if (u === b.b) break;
      for (const { j, k: bk } of nb[u]) {
        if (bk === k || !g.bonds[bk].arom || prev.has(j)) continue;
        prev.set(j, u);
        q.push(j);
      }
    }
    if (!prev.has(b.b)) return;
    const ring = [];
    for (let u = b.b; u !== -1; u = prev.get(u)) ring.push(u);
    if (ring.length > 8) return;
    const key = [...ring].sort((x, y) => x - y).join(',');
    if (!seen.has(key)) { seen.add(key); rings.push(ring); }
  });
  return rings;
}

export function formulaOf(g) {
  const count = {};
  for (const a of g.atoms) count[a.el] = (count[a.el] || 0) + 1;
  const order = ['C', 'H', ...Object.keys(count).filter((e) => e !== 'C' && e !== 'H').sort()];
  const sub = (n) => (n === 1 ? '' : String(n).replace(/\d/g, (d) => '₀₁₂₃₄₅₆₇₈₉'[d]));
  const heavy = Object.keys(count).filter((e) => e !== 'H').sort();
  // Hill order with carbon; otherwise heavy atoms then H (NH₃), but H first for acids/water (H₂O, HCl)
  const hFirst = heavy.length === 1 && ['O', 'S', 'F', 'Cl'].includes(heavy[0]);
  const keys = count.C ? order : hFirst ? ['H', ...heavy] : [...heavy, 'H'];
  let s = keys.filter((e) => count[e]).map((e) => e + sub(count[e])).join('');
  const q = g.atoms.reduce((t, a) => t + a.charge, 0);
  if (q) s += (Math.abs(q) > 1 ? Math.abs(q) : '') + (q > 0 ? '⁺' : '⁻');
  return s;
}
