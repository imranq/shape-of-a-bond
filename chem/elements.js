// Element data for geometry building and extended Hückel theory.
//   n      principal quantum number of the valence shell
//   zeta   Slater exponent (a.u.) used by EHT (Hoffmann parameters)
//   hs/hp  valence-state ionisation energies H_ii (eV)
//   ve     valence electrons
//   val    allowed valences for implicit hydrogens (SMILES organic subset)
//   r      single-bond covalent radius (Å)
export const BOHR = 0.529177;

export const EL = {
  H: { n: 1, zeta: 1.3, hs: -13.6, ve: 1, val: [1], r: 0.31, color: 0xf1f3f6, vis: 0.19 },
  B: { n: 2, zeta: 1.3, hs: -15.2, hp: -8.5, ve: 3, val: [3], r: 0.84, color: 0xf4a3a8, vis: 0.27 },
  C: { n: 2, zeta: 1.625, hs: -21.4, hp: -11.4, ve: 4, val: [4], r: 0.76, color: 0x2c3038, vis: 0.28 },
  N: { n: 2, zeta: 1.95, hs: -26.0, hp: -13.4, ve: 5, val: [3, 5], r: 0.71, color: 0x3f6ff0, vis: 0.28 },
  O: { n: 2, zeta: 2.275, hs: -32.3, hp: -14.8, ve: 6, val: [2], r: 0.66, color: 0xe4473e, vis: 0.29 },
  F: { n: 2, zeta: 2.425, hs: -40.0, hp: -18.1, ve: 7, val: [1], r: 0.57, color: 0x7fe07a, vis: 0.26 },
  P: { n: 3, zeta: 1.6, hs: -18.6, hp: -14.0, ve: 5, val: [3, 5], r: 1.07, color: 0xff9a3c, vis: 0.33 },
  S: { n: 3, zeta: 1.817, hs: -20.0, hp: -13.3, ve: 6, val: [2, 4, 6], r: 1.05, color: 0xf0cf3a, vis: 0.33 },
  Cl: { n: 3, zeta: 2.033, hs: -26.3, hp: -14.2, ve: 7, val: [1], r: 1.02, color: 0x46d160, vis: 0.32 },
};

// Average bond enthalpies (kJ/mol) for the potential-energy card
const BOND_KJ = {
  'H-H': 436, 'C-H': 413, 'N-H': 391, 'O-H': 463, 'S-H': 339, 'B-H': 389, 'F-H': 567, 'Cl-H': 431,
  'C-C': 348, 'C=C': 614, 'C#C': 839, 'C:C': 518,
  'C-N': 293, 'C=N': 615, 'C#N': 891, 'C:N': 500,
  'C-O': 358, 'C=O': 745, 'C#O': 1072, 'C:O': 470,
  'C-F': 485, 'C-Cl': 328, 'C-S': 272, 'C=S': 573, 'C:S': 400, 'C-P': 264, 'C-B': 356,
  'N-N': 163, 'N=N': 418, 'N#N': 945, 'N-O': 201, 'N=O': 607, 'N:N': 400, 'N:O': 400,
  'O-O': 146, 'O=O': 495, 'B-F': 613, 'B-N': 389, 'B-O': 536, 'P-O': 335, 'P=O': 544, 'S-O': 265, 'S=O': 522,
  'P-H': 322, 'P-Cl': 326, 'S-S': 266, 'Cl-Cl': 242, 'F-F': 155,
};

export function bondEnthalpy(a, b, order) {
  const sym = order === 1.5 ? ':' : order === 2 ? '=' : order === 3 ? '#' : '-';
  return BOND_KJ[`${a}${sym}${b}`] ?? BOND_KJ[`${b}${sym}${a}`] ?? 350;
}

export function bondSymbol(order) {
  return order === 1.5 ? '⋯' : order === 2 ? '=' : order === 3 ? '≡' : '–';
}
