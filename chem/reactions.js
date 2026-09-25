// Reaction lab. Each reaction is an atom-mapped "reactants>>products" SMILES.
// The path is built from three geometries sharing one atom ordering:
//   TS  – embedded from the union of reactant and product bonds, with changing bonds stretched
//   R/P – the reactant / product graphs relaxed from the TS and pulled apart into fragments
// ξ ∈ [0,1] interpolates R → TS → P, and orbitals are recomputed with EHT along the way.
import { parseSmiles, addHydrogens, components, neighbors, lonePairs } from './smiles.js';
import { buildModel, orient } from './embed.js';

export const REACTIONS = [
  {
    id: 'h2',
    title: 'H· + H· → H₂',
    kind: 'Radical combination',
    smiles: '[H:1].[H:2]>>[H:1][H:2]',
    Ea: 0, dH: -436,
    iso: 0.2,
    story: [
      'Two hydrogen atoms, each with one unpaired electron in a 1s orbital (a SOMO).',
      'The 1s waves overlap in phase. The two singly occupied orbitals mix into a bonding σ and an antibonding σ*.',
      'Both electrons pair up in σ. The energy released (436 kJ/mol) is the H–H bond energy.',
    ],
  },
  {
    id: 'lewis',
    title: 'NH₃ + BF₃ → H₃N–BF₃',
    kind: 'Lewis acid–base',
    smiles: '[NH3:1].[B:2]([F:3])([F:4])[F:5]>>[NH3+:1][B-:2]([F:3])([F:4])[F:5]',
    Ea: 0, dH: -130,
    story: [
      "Ammonia's HOMO is its lone pair. Boron in BF₃ is flat, with an empty 2p orbital as the LUMO.",
      'The lone pair points into the empty p orbital. Their overlap starts a new σ bond, and BF₃ begins to pyramidalise.',
      'A dative N→B bond has formed: both electrons came from nitrogen. Boron is now tetrahedral.',
    ],
  },
  {
    id: 'sn2',
    title: 'HO⁻ + CH₃Cl → CH₃OH + Cl⁻',
    kind: 'SN2 substitution',
    smiles: '[OH-:1].[CH3:2][Cl:3]>>[OH:1][CH3:2].[Cl-:3]',
    Ea: 103, dH: -75,
    story: [
      "Hydroxide's lone-pair HOMO approaches carbon from the back, aiming at the lobe of the C–Cl σ* LUMO.",
      'Transition state: carbon is five-coordinate and the three H atoms are flat. O–C forms as C–Cl breaks.',
      "The umbrella has flipped (Walden inversion). Chloride leaves with the pair that used to be C–Cl's.",
    ],
    energyNote: 'aqueous values',
  },
  {
    id: 'da',
    title: 'Butadiene + Ethylene → Cyclohexene',
    kind: 'Diels–Alder [4+2]',
    smiles: '[CH2:1]=[CH:2][CH:3]=[CH2:4].[CH2:5]=[CH2:6]>>[CH2:1]1[CH:2]=[CH:3][CH2:4][CH2:6][CH2:5]1',
    Ea: 115, dH: -166,
    story: [
      "The diene's HOMO (ψ₂) has matching phases at both ends, the same as the dienophile's π* LUMO. The overlap is symmetry-allowed.",
      'Transition state: both new σ bonds form at the same time (concerted), and the π bonds shift to the middle.',
      'Two π bonds have become two stronger σ bonds, leaving one new C=C in the ring.',
    ],
  },
  {
    id: 'methyl',
    title: '2 CH₃· → C₂H₆',
    kind: 'Radical combination',
    smiles: '[CH3:1].[CH3:2]>>[CH3:1][CH3:2]',
    Ea: 0, dH: -377,
    story: [
      'Each methyl radical is flat, with its unpaired electron in a 2p orbital sticking out of the plane.',
      'The p orbitals point at each other and overlap end-on. The radicals start to pyramidalise.',
      'A C–C σ bond forms and each carbon becomes tetrahedral (sp³).',
    ],
  },
  {
    id: 'proton',
    title: 'H₂O + H⁺ → H₃O⁺',
    kind: 'Proton transfer',
    smiles: '[OH2:1].[H+:2]>>[OH2+:1][H:2]',
    Ea: 0, dH: -691,
    story: [
      "A bare proton has an empty 1s orbital. Water's HOMO is a lone pair on oxygen.",
      'The lone pair reaches toward the proton, and O–H overlap grows as the bond forms.',
      'The hydronium ion is pyramidal. Its three O–H bonds are equivalent and the + charge is shared.',
    ],
    energyNote: 'gas-phase proton affinity',
  },
];

function mappedGraph(smi) {
  const g = addHydrogens(parseSmiles(smi));
  const keys = g.atoms.map((a) => {
    if (a.el === 'H' && a.parent !== undefined) {
      const p = g.atoms[a.parent];
      if (!p.map) throw new Error('Every heavy atom needs an atom map number');
      return `h${a.hkey}`;
    }
    if (!a.map) throw new Error('Every heavy atom needs an atom map number');
    return `m${a.map}`;
  });
  return { g, keys };
}

export function buildReaction(def) {
  const [rs, ps] = def.smiles.split('>>');
  const R = mappedGraph(rs), P = mappedGraph(ps);
  // unify atom order on the reactant side; implicit hydrogens pair up by parent map
  const pIndex = new Map(P.keys.map((k, i) => [k, i]));
  const toP = R.keys.map((k) => {
    if (!pIndex.has(k)) throw new Error(`Atom ${k} has no partner in the products`);
    return pIndex.get(k);
  });
  if (P.keys.length !== R.keys.length) throw new Error('Reactants and products have different atom counts');
  const fromP = new Array(toP.length);
  toP.forEach((p, r) => (fromP[p] = r));

  const atomsR = R.g.atoms;
  const atomsP = toP.map((p) => P.g.atoms[p]);
  const pk = (a, b) => (a < b ? `${a},${b}` : `${b},${a}`);
  const bondsR = R.g.bonds.map((b) => ({ ...b }));
  const bondsP = P.g.bonds.map((b) => ({ ...b, a: fromP[b.a], b: fromP[b.b] }));
  const rMap = new Map(bondsR.map((b) => [pk(b.a, b.b), b]));
  const pMap = new Map(bondsP.map((b) => [pk(b.a, b.b), b]));

  // union graph for the transition state
  const union = [];
  const changing = new Set();
  const forming = [], breaking = [];
  for (const key of new Set([...rMap.keys(), ...pMap.keys()])) {
    const br = rMap.get(key), bp = pMap.get(key);
    const [a, b] = key.split(',').map(Number);
    const order = br && bp ? (br.order + bp.order) / 2 : (br || bp).order;
    union.push({ a, b, order, rOrder: br ? br.order : 0, pOrder: bp ? bp.order : 0, arom: false });
    if (!br || !bp) {
      changing.add(union.length - 1);
      (br ? breaking : forming).push([a, b]);
    }
  }
  const gR = { atoms: atomsR, bonds: bondsR };
  const gP = { atoms: atomsP, bonds: bondsP };
  const gTS = {
    atoms: atomsR.map((a, i) => ({ ...a, charge: 0 })),
    bonds: union,
  };
  // lone pairs at the TS: the smaller of the two ends (a pair being donated is half-used)
  const nbR = neighbors(gR), nbP = neighbors(gP);
  const lpTS = atomsR.map((_, i) => Math.min(lonePairs(gR, i, nbR), lonePairs(gP, i, nbP)));
  // backside attack: a centre that gains one bond and loses another keeps them collinear
  const linear = [];
  for (const [fa, fb] of forming) for (const [ba, bb] of breaking) {
    for (const c of [fa, fb]) for (const d of [ba, bb]) {
      if (c === d) linear.push([fa === c ? fb : fa, c, ba === d ? bb : ba]);
    }
  }
  const ts = buildModel(gTS, { long: changing, linear, lpOverride: lpTS });
  const tsPos = orient(ts.pos, gTS.atoms.map((a) => (a.el === 'H' ? 0.3 : 1)));

  const separate = (g, pos) => {
    const { comp, count } = components(g);
    if (count < 2) return pos;
    const out = pos.map((p) => p.slice());
    const cents = [], sizes = [];
    for (let c = 0; c < count; c++) {
      const idx = comp.map((x, i) => (x === c ? i : -1)).filter((i) => i >= 0);
      sizes.push(idx.reduce((s, i) => s + (g.atoms[i].el === 'H' ? 1 : 12), 0));
      cents.push([0, 1, 2].map((k) => idx.reduce((s, i) => s + pos[i][k], 0) / idx.length));
    }
    const total = sizes.reduce((a, b) => a + b, 0);
    const G = [0, 1, 2].map((k) => cents.reduce((s, c, i) => s + (c[k] * sizes[i]) / total, 0));
    comp.forEach((c, i) => {
      let d = [0, 1, 2].map((k) => cents[c][k] - G[k]);
      let L = Math.hypot(...d);
      if (L < 1e-3) return;
      // fragments end up ~2.3 Å further apart than at the transition state
      const push = 2.3 * (1 - sizes[c] / total);
      for (let k = 0; k < 3; k++) out[i][k] += (d[k] / L) * push;
    });
    return out;
  };

  const relaxed = (g, init) => buildModel(g, { init, iters: 900 }).pos;
  const rPos = relaxed(gR, separate(gR, tsPos));
  const pPos = relaxed(gP, separate(gP, tsPos));

  return {
    def,
    atoms: atomsR,
    atomsP,
    bondsUnion: union,
    forming, breaking,
    gR, gP,
    charge: atomsR.reduce((s, a) => s + a.charge, 0),
    fragR: components(gR),
    frames: { R: rPos, TS: tsPos, P: pPos },
  };
}

const smooth = (t) => t * t * (3 - 2 * t);
export function pathAt(rx, xi) {
  const { R, TS, P } = rx.frames;
  const [A, B, t] = xi < 0.5 ? [R, TS, smooth(xi * 2)] : [TS, P, smooth(xi * 2 - 1)];
  return A.map((p, i) => [0, 1, 2].map((k) => p[k] + (B[i][k] - p[k]) * t));
}

// schematic energy profile (kJ/mol): Hermite through R(0), TS(Ea), P(ΔH)
export function profileAt(def, xi) {
  const { Ea, dH } = def;
  if (Ea <= 0) {
    // barrierless: smooth monotone descent
    const t = smooth(Math.min(1, Math.max(0, (xi - 0.15) / 0.75)));
    return dH * t;
  }
  if (xi < 0.5) return Ea * smooth(xi * 2);
  return Ea + (dH - Ea) * smooth(xi * 2 - 1);
}
