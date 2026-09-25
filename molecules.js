// Molecule library. Coordinates in ångström; AO exponents (z) in Å⁻¹ (Slater exponents / a₀).
// MO coefficients are qualitative LCAO combinations (Hückel-exact for benzene π);
// they are renormalised against the numerical overlap matrix at runtime.
//
// MO fields:
//   name  – label with _{sub} / ^{sup} markup
//   e     – orbital energy at equilibrium (eV): photoelectron data for occupied levels,
//           estimates for virtual ones
//   ref   – energy the level relaxes to when the bond is pulled apart (atomic parent / nonbonding)
//   g     – degeneracy group id (levels sharing g are drawn side by side)
//   from  – atomic levels it correlates with (diatomics only)

const H = 2.34;  // H 1s in a molecule (ζ = 1.24 a.u.)
const C = 3.07;  // C 2s/2p (ζ = 1.625 a.u.)
const N = 3.69;  // N 2s/2p (ζ = 1.95 a.u.)
const O = 4.30;  // O 2s/2p (ζ = 2.275 a.u.)

export const ELEMENTS = {
  H: { color: 0xf1f3f6, radius: 0.19, rough: 0.18 },
  C: { color: 0x2c3038, radius: 0.28, rough: 0.3 },
  N: { color: 0x3f6ff0, radius: 0.28, rough: 0.25 },
  O: { color: 0xe4473e, radius: 0.3, rough: 0.25 },
};

const valence = (atom, zeta) => [
  { atom, t: 's', n: 2, z: zeta },
  { atom, t: 'px', n: 2, z: zeta },
  { atom, t: 'py', n: 2, z: zeta },
  { atom, t: 'pz', n: 2, z: zeta },
];
const h1s = (atom) => ({ atom, t: 's', n: 1, z: H });

// ---------- benzene (Hückel π system, ring in the xz-plane, π lobes along y) ----------
const RING = 1.397, RING_H = 2.484;
const ang = [0, 1, 2, 3, 4, 5].map((j) => (j * Math.PI) / 3);
const benzeneAtoms = [
  ...ang.map((a, j) => ({ el: 'C', p: [RING * Math.cos(a), 0, RING * Math.sin(a)], label: j === 0 ? 'C · sp²' : null })),
  ...ang.map((a, j) => ({ el: 'H', p: [RING_H * Math.cos(a), 0, RING_H * Math.sin(a)], label: j === 0 ? 'H' : null })),
];
const huckel = (k, fn) => ang.map((a) => +fn(k * a).toFixed(6));
const ALPHA = -6.8, BETA = -2.45;

export const MOLECULES = [
  {
    id: 'h2',
    formula: 'H₂',
    name: 'Hydrogen',
    blurb: 'The simplest bond: two protons sharing two electrons.',
    atoms: [
      { el: 'H', p: [-0.3705, 0, 0], label: 'H · 1s' },
      { el: 'H', p: [0.3705, 0, 0], label: 'H · 1s' },
    ],
    bonds: [[0, 1, 1]],
    rep: [0, 1],
    bondName: 'H–H',
    order: '1',
    morse: { De: 4.52, Re: 0.741, a: 1.94 },
    basis: [h1s(0), h1s(1)],
    atomic: [{ name: '1s', e: -13.6, n: 1, el: 1 }],
    mos: [
      {
        name: 'σ_{g}1s', e: -15.4, ref: -13.6, occ: 2, g: 0, from: ['1s'], c: [1, 1],
        desc: 'Both 1s waves add in phase. Electron density piles up between the nuclei and pulls them together. That shared density is the covalent bond.',
      },
      {
        name: 'σ*_{u}1s', e: -1.0, ref: -13.6, occ: 0, g: 1, from: ['1s'], c: [1, -1],
        desc: 'The same waves, out of phase. A nodal plane cuts between the nuclei, so electrons placed here push the atoms apart. It is destabilised more than σ is stabilised.',
      },
    ],
    defaultMO: 0, iso: 0.22, cam: 7.2,
  },
  {
    id: 'n2',
    formula: 'N₂',
    name: 'Nitrogen',
    blurb: 'A triple bond: one σ and two π, among the strongest bonds in chemistry.',
    atoms: [
      { el: 'N', p: [-0.549, 0, 0], label: 'N · sp' },
      { el: 'N', p: [0.549, 0, 0], label: 'N · sp' },
    ],
    bonds: [[0, 1, 3]],
    rep: [0, 1],
    bondName: 'N≡N',
    order: '3',
    morse: { De: 9.79, Re: 1.098, a: 2.69 },
    basis: [...valence(0, N), ...valence(1, N)],
    atomic: [
      { name: '2s', e: -25.6, n: 1, el: 2 },
      { name: '2p', e: -13.2, n: 3, el: 3 },
    ],
    mos: [
      { name: '2σ_{g}', e: -37.3, ref: -25.6, occ: 2, g: 0, from: ['2s'], c: [0.6, 0.15, 0, 0, 0.6, -0.15, 0, 0],
        desc: 'The 2s orbitals overlap in phase. The orbital is deep and compact and strongly bonding.' },
      { name: '2σ*_{u}', e: -18.8, ref: -25.6, occ: 2, g: 1, from: ['2s'], c: [0.6, -0.2, 0, 0, -0.6, -0.2, 0, 0],
        desc: 'The 2s orbitals combine out of phase. Its node cancels most of what 2σ_g built, so the bond is left to the 2p orbitals.' },
      { name: '1π_{u} (y)', e: -16.7, ref: -13.2, occ: 2, g: 2, from: ['2p'], c: [0, 0, 0.7, 0, 0, 0, 0.7, 0],
        desc: 'The 2p orbitals overlap side-on, above and below the axis. This orbital and its twin make two of the three bonds.' },
      { name: '1π_{u} (z)', e: -16.7, ref: -13.2, occ: 2, g: 2, from: ['2p'], c: [0, 0, 0, 0.7, 0, 0, 0, 0.7],
        desc: 'This is the second π bond, rotated 90°. Together the two π orbitals form a cylinder of density around the axis.' },
      { name: '3σ_{g}', e: -15.6, ref: -13.2, occ: 2, g: 3, from: ['2p'], c: [-0.3, 0.6, 0, 0, -0.3, -0.6, 0, 0], tag: 'HOMO',
        desc: 'The HOMO. The 2p orbitals overlap end-on. s–p mixing moves density to the outside, which gives N₂ its lone-pair-like lobes.' },
      { name: '1π*_{g} (y)', e: 2.3, ref: -13.2, occ: 0, g: 4, from: ['2p'], c: [0, 0, 0.7, 0, 0, 0, -0.7, 0], tag: 'LUMO',
        desc: 'The LUMO. The 2p orbitals combine side-on and out of phase. It is empty in N₂. In O₂ it holds two electrons, which is why O₂ has a double bond.' },
      { name: '1π*_{g} (z)', e: 2.3, ref: -13.2, occ: 0, g: 4, from: ['2p'], c: [0, 0, 0, 0.7, 0, 0, 0, -0.7],
        desc: 'This is the degenerate partner of the π* LUMO. It has a node between the atoms and a node along the axis.' },
      { name: '3σ*_{u}', e: 6.0, ref: -13.2, occ: 0, g: 5, from: ['2p'], c: [0.25, 0.6, 0, 0, -0.25, 0.6, 0, 0],
        desc: 'The 2p orbitals combine end-on and out of phase. It is the most antibonding valence orbital.' },
    ],
    defaultMO: 4, iso: 0.12, cam: 7.8,
  },
  {
    id: 'h2o',
    formula: 'H₂O',
    name: 'Water',
    blurb: 'A bent molecule: two O–H bonds plus two lone pairs on the oxygen.',
    atoms: [
      { el: 'O', p: [0, 0.195, 0], label: 'O · sp³' },
      { el: 'H', p: [0.757, -0.391, 0], label: 'H · 1s' },
      { el: 'H', p: [-0.757, -0.391, 0], label: null },
    ],
    bonds: [[0, 1, 1], [0, 2, 1]],
    rep: [0, 1],
    bondName: 'O–H',
    order: '1',
    morse: { De: 5.1, Re: 0.958, a: 2.27 },
    basis: [...valence(0, O), h1s(1), h1s(2)],
    mos: [
      { name: '2a_{1}', e: -32.2, ref: -30.5, occ: 2, g: 0, c: [0.8, 0, 0, 0, 0.18, 0.18],
        desc: 'Mostly O 2s with a little of both hydrogens. It is low in energy, compact and bonding everywhere.' },
      { name: '1b_{2}', e: -18.7, ref: -15.0, occ: 2, g: 1, c: [0, 0.55, 0, 0, 0.4, -0.4],
        desc: 'O 2pₓ bonds to the two hydrogens with opposite signs. This orbital holds most of the O–H bonding along the H···H direction.' },
      { name: '3a_{1}', e: -14.7, ref: -13.3, occ: 2, g: 2, c: [-0.25, 0, -0.55, 0, 0.3, 0.3],
        desc: 'O 2p along the bisector, mixed with 2s. It bonds to the hydrogens below and bulges into a lone-pair lobe above.' },
      { name: '1b_{1}', e: -12.6, ref: -12.6, occ: 2, g: 3, c: [0, 0, 0, 1, 0, 0], tag: 'HOMO',
        desc: 'The HOMO is a pure O 2p lone pair that sticks out of the molecular plane. It has no hydrogen character, so its energy stays flat as you stretch the bonds.' },
      { name: '4a_{1}', e: 2.0, ref: -10.0, occ: 0, g: 4, c: [0.45, 0, -0.45, 0, -0.6, -0.6], tag: 'LUMO',
        desc: 'The LUMO is antibonding along the bisector. Every O–H bond has a node in it.' },
      { name: '2b_{2}', e: 4.0, ref: -10.0, occ: 0, g: 5, c: [0, 0.5, 0, 0, -0.6, 0.6],
        desc: 'The antibonding partner of 1b₂. The phase flips between oxygen and each hydrogen.' },
    ],
    defaultMO: 1, iso: 0.12, cam: 8.0,
  },
  {
    id: 'ch4',
    formula: 'CH₄',
    name: 'Methane',
    blurb: 'Four equivalent C–H bonds come from delocalised orbitals, not four separate ones.',
    atoms: [
      { el: 'C', p: [0, 0, 0], label: 'C · sp³' },
      { el: 'H', p: [0.6276, 0.6276, 0.6276], label: 'H · 1s' },
      { el: 'H', p: [0.6276, -0.6276, -0.6276], label: null },
      { el: 'H', p: [-0.6276, 0.6276, -0.6276], label: null },
      { el: 'H', p: [-0.6276, -0.6276, 0.6276], label: null },
    ],
    bonds: [[0, 1, 1], [0, 2, 1], [0, 3, 1], [0, 4, 1]],
    rep: [0, 1],
    bondName: 'C–H',
    order: '1',
    morse: { De: 4.55, Re: 1.087, a: 1.8 },
    basis: [...valence(0, C), ...[1, 2, 3, 4].map(h1s)],
    mos: [
      { name: '2a_{1}', e: -23.0, ref: -19.0, occ: 2, g: 0, c: [0.6, 0, 0, 0, 0.25, 0.25, 0.25, 0.25],
        desc: 'Carbon 2s surrounded by all four hydrogens in phase. The result is one large bonding sphere.' },
      { name: '1t_{2} (x)', e: -14.4, ref: -12.0, occ: 2, g: 1, c: [0, 0.55, 0, 0, 0.3, 0.3, -0.3, -0.3],
        desc: 'This is one of three degenerate t₂ orbitals, one per axis. Together with 2a₁ they hold the four C–H bonds, and no single orbital is one bond.' },
      { name: '1t_{2} (y)', e: -14.4, ref: -12.0, occ: 2, g: 1, c: [0, 0, 0.55, 0, 0.3, -0.3, 0.3, -0.3], tag: 'HOMO',
        desc: 'The same shape as the x orbital, rotated. Degenerate orbitals are symmetry copies of each other.' },
      { name: '1t_{2} (z)', e: -14.4, ref: -12.0, occ: 2, g: 1, c: [0, 0, 0, 0.55, 0.3, -0.3, -0.3, 0.3],
        desc: 'The third t₂ orbital. Adding the densities of all three gives a perfectly tetrahedral cloud.' },
      { name: '3a_{1}*', e: 3.5, ref: -15.0, occ: 0, g: 2, c: [0.6, 0, 0, 0, -0.35, -0.35, -0.35, -0.35], tag: 'LUMO',
        desc: 'The LUMO. Carbon 2s is out of phase with all four hydrogens, so every C–H bond has a node.' },
      { name: '2t_{2}* (x)', e: 6.0, ref: -9.0, occ: 0, g: 3, c: [0, 0.6, 0, 0, -0.35, -0.35, 0.35, 0.35],
        desc: 'The antibonding partners of the t₂ set.' },
      { name: '2t_{2}* (y)', e: 6.0, ref: -9.0, occ: 0, g: 3, c: [0, 0, 0.6, 0, -0.35, 0.35, -0.35, 0.35],
        desc: 'The antibonding partners of the t₂ set.' },
      { name: '2t_{2}* (z)', e: 6.0, ref: -9.0, occ: 0, g: 3, c: [0, 0, 0, 0.6, -0.35, 0.35, 0.35, -0.35],
        desc: 'The antibonding partners of the t₂ set.' },
    ],
    defaultMO: 1, iso: 0.1, cam: 8.4,
  },
  {
    id: 'c2h4',
    formula: 'C₂H₄',
    name: 'Ethylene',
    blurb: 'A double bond: a σ framework in the plane plus one π bond above and below it.',
    atoms: [
      { el: 'C', p: [-0.6695, 0, 0], label: 'C · sp²' },
      { el: 'C', p: [0.6695, 0, 0], label: null },
      { el: 'H', p: [-1.2345, 0, 0.929], label: 'H · 1s' },
      { el: 'H', p: [-1.2345, 0, -0.929], label: null },
      { el: 'H', p: [1.2345, 0, 0.929], label: null },
      { el: 'H', p: [1.2345, 0, -0.929], label: null },
    ],
    bonds: [[0, 1, 2, [0, 0, 1]], [0, 2, 1], [0, 3, 1], [1, 4, 1], [1, 5, 1]],
    rep: [0, 1],
    bondName: 'C=C',
    order: '2',
    morse: { De: 7.55, Re: 1.339, a: 2.0 },
    basis: [...valence(0, C), ...valence(1, C), ...[2, 3, 4, 5].map(h1s)],
    mos: [
      { name: 'σ 2a_{g}', e: -23.7, ref: -19.5, occ: 2, g: 0, c: [0.5, 0, 0, 0, 0.5, 0, 0, 0, 0.15, 0.15, 0.15, 0.15],
        desc: 'Both carbon 2s orbitals and all four hydrogens combine in phase. It is the lowest valence orbital and bonds everywhere.' },
      { name: 'σ 2b_{1u}', e: -19.1, ref: -16.0, occ: 2, g: 1, c: [0.45, -0.15, 0, 0, -0.45, -0.15, 0, 0, 0.22, 0.22, -0.22, -0.22],
        desc: 'The 2s orbitals combine out of phase across C–C but in phase with each CH₂ group. It has a node in the middle.' },
      { name: 'σ_{CH} 1b_{2u}', e: -15.9, ref: -13.0, occ: 2, g: 2, c: [0, 0, 0, 0.4, 0, 0, 0, 0.4, 0.3, -0.3, 0.3, -0.3],
        desc: 'In-plane p orbitals bond to the hydrogens. This is C–H bonding spread over the whole molecule.' },
      { name: 'σ_{CC} 3a_{g}', e: -14.7, ref: -11.5, occ: 2, g: 3, c: [-0.12, 0.5, 0, 0, -0.12, -0.5, 0, 0, 0, 0, 0, 0],
        desc: 'The sp² lobes overlap head-on along the axis. This is the σ half of the double bond.' },
      { name: 'σ_{CH} 1b_{3g}', e: -12.9, ref: -11.0, occ: 2, g: 4, c: [0, 0, 0, 0.4, 0, 0, 0, -0.4, 0.3, -0.3, -0.3, 0.3],
        desc: 'C–H bonding that is antisymmetric across the C–C bond. It has a node through the middle.' },
      { name: 'π', e: -10.5, ref: -4.35, occ: 2, g: 5, c: [0, 0, 0.7, 0, 0, 0, 0.7, 0, 0, 0, 0, 0], tag: 'HOMO',
        desc: 'The HOMO. The unhybridised 2p orbitals overlap side-on above and below the plane. This π bond stops the molecule from twisting.' },
      { name: 'π*', e: 1.8, ref: -4.35, occ: 0, g: 6, c: [0, 0, 0.7, 0, 0, 0, -0.7, 0, 0, 0, 0, 0], tag: 'LUMO',
        desc: 'The LUMO. The same p orbitals combine out of phase. Light at about 7 eV moves an electron into it (π→π*), which breaks the π bond.' },
      { name: 'σ*_{CC}', e: 4.5, ref: -11.5, occ: 0, g: 7, c: [0.2, 0.55, 0, 0, -0.2, 0.55, 0, 0, 0, 0, 0, 0],
        desc: 'The sp² lobes combine head-on and out of phase. This is the antibonding partner of the C–C σ bond.' },
    ],
    defaultMO: 5, iso: 0.1, cam: 9.2,
  },
  {
    id: 'c6h6',
    formula: 'C₆H₆',
    name: 'Benzene',
    blurb: 'Six π electrons spread evenly around the ring. This delocalisation is aromaticity.',
    atoms: benzeneAtoms,
    bonds: [
      ...ang.map((_, j) => [j, (j + 1) % 6, 1.5]),
      ...ang.map((_, j) => [j, j + 6, 1]),
    ],
    rep: [0, 1],
    bondName: 'C–C (arom.)',
    order: '1½',
    morse: { De: 5.37, Re: 1.397, a: 2.0 },
    basis: ang.map((_, j) => ({ atom: j, t: 'py', n: 2, z: C })),
    mos: [
      { name: 'π_{1} a_{2u}', e: ALPHA + 2 * BETA, ref: ALPHA, occ: 2, g: 0, c: huckel(0, () => 1),
        desc: 'All six p orbitals are in phase. They form two continuous rings of density above and below the plane. This orbital is why benzene is so stable.' },
      { name: 'π_{2} e_{1g}', e: ALPHA + BETA, ref: ALPHA, occ: 2, g: 1, c: huckel(1, Math.cos), tag: 'HOMO',
        desc: 'One half of the degenerate HOMO pair. A single nodal plane cuts through the ring.' },
      { name: 'π_{3} e_{1g}', e: ALPHA + BETA, ref: ALPHA, occ: 2, g: 1, c: huckel(1, Math.sin),
        desc: 'The other half of the HOMO pair, with its nodal plane rotated 90°. π₁, π₂ and π₃ together hold the 6 π electrons.' },
      { name: 'π_{4}* e_{2u}', e: ALPHA - BETA, ref: ALPHA, occ: 0, g: 2, c: huckel(2, Math.cos), tag: 'LUMO',
        desc: 'The degenerate LUMO pair. Two nodal planes cross the ring, and the orbital is empty in the ground state.' },
      { name: 'π_{5}* e_{2u}', e: ALPHA - BETA, ref: ALPHA, occ: 0, g: 2, c: huckel(2, Math.sin),
        desc: 'This is the second LUMO, a symmetry copy of π₄* rotated by 45°.' },
      { name: 'π_{6}* b_{2g}', e: ALPHA - 2 * BETA, ref: ALPHA, occ: 0, g: 3, c: huckel(3, Math.cos),
        desc: 'Every neighbouring pair is out of phase. Three nodal planes make it the most antibonding π orbital.' },
    ],
    defaultMO: 0, iso: 0.08, cam: 11.5,
  },
];
