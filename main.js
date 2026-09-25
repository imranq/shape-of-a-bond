import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { RoomEnvironment } from 'three/addons/environments/RoomEnvironment.js';
import { EffectComposer } from 'three/addons/postprocessing/EffectComposer.js';
import { RenderPass } from 'three/addons/postprocessing/RenderPass.js';
import { UnrealBloomPass } from 'three/addons/postprocessing/UnrealBloomPass.js';
import { OutputPass } from 'three/addons/postprocessing/OutputPass.js';
import { ShaderPass } from 'three/addons/postprocessing/ShaderPass.js';
import { Pass, FullScreenQuad } from 'three/addons/postprocessing/Pass.js';

import { MOLECULES } from './molecules.js';
import { buildAOs, overlapMatrix as gridOverlap, normOf, sampleCloud, bounds } from './wavefunction.js';
import {
  MAX_AO, MAX_ATOMS, MAX_BONDS,
  orbitalVert, orbitalFrag, cloudVert, cloudFrag, backdropFrag, compositeFrag, finishFrag,
} from './shaders.js';
import { MODiagram, MorseChart, ProfileChart, richHTML } from './ui.js';
import { EL, bondEnthalpy } from './chem/elements.js';
import { parseSmiles, addHydrogens, aromaticRings, formulaOf, neighbors } from './chem/smiles.js';
import { buildModel, orient, bondLength } from './chem/embed.js';
import { solveEHT, valenceBasis, analyseMO } from './chem/eht.js';
import { resonanceStructures } from './chem/resonance.js';
import { REACTIONS, buildReaction, pathAt, profileAt } from './chem/reactions.js';

const S_MIN = 0.6, S_MAX = 2.6;
const CLOUD_N = 4500;
const COL_POS = new THREE.Color('#4fd8ff');
const COL_NEG = new THREE.Color('#ffab3d');
const STAGE_Y = -2.35;
const $ = (id) => document.getElementById(id);
const fmt = (v, d = 2) => v.toFixed(d).replace('-', '−');

// ============================================================ renderer / scene
const canvas = $('gl');
const renderer = new THREE.WebGLRenderer({ canvas, antialias: false, powerPreference: 'high-performance' });
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
renderer.setSize(window.innerWidth, window.innerHeight);
renderer.toneMapping = THREE.ACESFilmicToneMapping;
renderer.toneMappingExposure = 1.0;
renderer.shadowMap.enabled = true;
renderer.shadowMap.type = THREE.PCFSoftShadowMap;
renderer.setClearColor(0x0b0d11, 1);

const scene = new THREE.Scene();
scene.fog = new THREE.Fog(0x0b0d11, 18, 48);
const pmrem = new THREE.PMREMGenerator(renderer);
scene.environment = pmrem.fromScene(new RoomEnvironment(), 0.04).texture;
scene.environmentIntensity = 0.55;

const camera = new THREE.PerspectiveCamera(36, window.innerWidth / window.innerHeight, 0.05, 200);
camera.position.set(9, 6, 16);

const controls = new OrbitControls(camera, canvas);
controls.enableDamping = true;
controls.dampingFactor = 0.06;
controls.minDistance = 2.2;
controls.maxDistance = 30;
controls.maxPolarAngle = Math.PI * 0.62;
controls.autoRotate = true;
controls.autoRotateSpeed = 0.45;
controls.target.set(0, -0.15, 0);

// ============================================================ environment & stage
{
  scene.add(new THREE.Mesh(
    new THREE.SphereGeometry(80, 48, 24),
    new THREE.ShaderMaterial({
      side: THREE.BackSide, depthWrite: false, fog: false,
      vertexShader: 'varying vec3 vWorld; void main(){ vWorld = position; gl_Position = projectionMatrix * modelViewMatrix * vec4(position,1.); }',
      fragmentShader: backdropFrag,
    }),
  ));
  scene.add(new THREE.HemisphereLight(0xbfd6ff, 0x1a120a, 0.35));
  const key = new THREE.DirectionalLight(0xfff1dc, 2.4);
  key.position.set(-5, 9, 6);
  key.castShadow = true;
  key.shadow.mapSize.set(2048, 2048);
  key.shadow.camera.left = key.shadow.camera.bottom = -8;
  key.shadow.camera.right = key.shadow.camera.top = 8;
  key.shadow.camera.near = 1;
  key.shadow.camera.far = 30;
  key.shadow.radius = 6;
  key.shadow.bias = -0.0004;
  scene.add(key);
  const rim = new THREE.DirectionalLight(0x8fb8ff, 1.3);
  rim.position.set(6, 3, -8);
  scene.add(rim);
  const warm = new THREE.PointLight(0xffb45c, 18, 20, 2);
  warm.position.set(0, STAGE_Y + 0.6, -6);
  scene.add(warm);
}

const gold = new THREE.MeshStandardMaterial({ color: 0xd4a24c, metalness: 1, roughness: 0.3 });
const goldBright = new THREE.MeshStandardMaterial({ color: 0xe8b85a, metalness: 1, roughness: 0.22, emissive: 0x3a2508, emissiveIntensity: 0.6 });
const plateMat = new THREE.MeshStandardMaterial({ color: 0x15171c, metalness: 0.7, roughness: 0.38 });

{
  const c = document.createElement('canvas');
  c.width = c.height = 512;
  const g = c.getContext('2d');
  g.fillStyle = '#0e1014';
  g.fillRect(0, 0, 512, 512);
  g.strokeStyle = 'rgba(255,255,255,0.05)';
  g.lineWidth = 2;
  for (let i = 0; i <= 512; i += 64) {
    g.beginPath(); g.moveTo(i, 0); g.lineTo(i, 512); g.stroke();
    g.beginPath(); g.moveTo(0, i); g.lineTo(512, i); g.stroke();
  }
  const tex = new THREE.CanvasTexture(c);
  tex.wrapS = tex.wrapT = THREE.RepeatWrapping;
  tex.repeat.set(20, 20);
  tex.colorSpace = THREE.SRGBColorSpace;
  tex.anisotropy = 8;
  const floor = new THREE.Mesh(new THREE.PlaneGeometry(90, 90), new THREE.MeshStandardMaterial({ map: tex, roughness: 0.82, metalness: 0.1 }));
  floor.rotation.x = -Math.PI / 2;
  floor.position.y = STAGE_Y - 0.42;
  floor.receiveShadow = true;
  scene.add(floor);

  const base = new THREE.Mesh(new THREE.CylinderGeometry(4.75, 4.95, 0.4, 128), plateMat);
  base.position.y = STAGE_Y - 0.2;
  base.receiveShadow = base.castShadow = true;
  scene.add(base);
  const top = new THREE.Mesh(new THREE.CylinderGeometry(4.3, 4.3, 0.02, 128), new THREE.MeshStandardMaterial({ color: 0x0c0d10, roughness: 0.55, metalness: 0.4 }));
  top.position.y = STAGE_Y + 0.005;
  top.receiveShadow = true;
  scene.add(top);
  const ring = new THREE.Mesh(new THREE.TorusGeometry(4.32, 0.025, 12, 200), goldBright);
  ring.rotation.x = Math.PI / 2;
  ring.position.y = STAGE_Y + 0.02;
  scene.add(ring);

  // engraved degree ticks around the rim, like a focus ring
  const m = new THREE.Matrix4(), q = new THREE.Quaternion(), s = new THREE.Vector3(), p = new THREE.Vector3();
  const Y = new THREE.Vector3(0, 1, 0);
  const ticks = new THREE.InstancedMesh(new THREE.BoxGeometry(1, 1, 1), gold, 180);
  for (let i = 0; i < 180; i++) {
    const a = (i / 180) * Math.PI * 2;
    const len = i % 15 === 0 ? 0.3 : i % 5 === 0 ? 0.2 : 0.11;
    p.set(Math.cos(a) * (4.42 + len / 2), STAGE_Y + 0.012, Math.sin(a) * (4.42 + len / 2));
    q.setFromAxisAngle(Y, -a);
    s.set(len, 0.02, i % 15 === 0 ? 0.035 : 0.022);
    ticks.setMatrixAt(i, m.compose(p, q, s));
  }
  scene.add(ticks);
  const knurl = new THREE.InstancedMesh(new THREE.BoxGeometry(1, 1, 1), gold, 240);
  for (let i = 0; i < 240; i++) {
    const a = (i / 240) * Math.PI * 2;
    p.set(Math.cos(a) * 4.93, STAGE_Y - 0.2, Math.sin(a) * 4.93);
    q.setFromAxisAngle(Y, -a);
    s.set(0.06, 0.3, 0.05);
    knurl.setMatrixAt(i, m.compose(p, q, s));
  }
  knurl.castShadow = true;
  scene.add(knurl);
}

// brass placard on the front edge
const plate = (() => {
  const c = document.createElement('canvas');
  c.width = 512; c.height = 128;
  const tex = new THREE.CanvasTexture(c);
  tex.colorSpace = THREE.SRGBColorSpace;
  tex.anisotropy = 8;
  const mat = new THREE.MeshStandardMaterial({ map: tex, metalness: 0.9, roughness: 0.35 });
  const mesh = new THREE.Mesh(new THREE.BoxGeometry(1.9, 0.46, 0.06), [gold, gold, gold, gold, mat, gold]);
  mesh.position.set(0, STAGE_Y - 0.16, 5.02);
  mesh.rotation.x = -0.12;
  mesh.castShadow = true;
  scene.add(mesh);
  let last = null;
  const draw = (title, sub) => {
    last = [title, sub];
    const g = c.getContext('2d');
    const grd = g.createLinearGradient(0, 0, 512, 128);
    grd.addColorStop(0, '#e2b560'); grd.addColorStop(0.5, '#f3cf85'); grd.addColorStop(1, '#c8973f');
    g.fillStyle = grd;
    g.fillRect(0, 0, 512, 128);
    g.strokeStyle = 'rgba(60,35,5,0.55)';
    g.lineWidth = 4;
    g.strokeRect(10, 10, 492, 108);
    g.fillStyle = '#3a2508';
    g.textAlign = 'center';
    let size = 48;
    g.font = `700 ${size}px Lexend, sans-serif`;
    while (g.measureText(title).width > 470 && size > 20) g.font = `700 ${(size -= 2)}px Lexend, sans-serif`;
    g.fillText(title, 256, 66);
    g.font = '500 20px "JetBrains Mono", monospace';
    g.fillText(sub, 256, 100);
    tex.needsUpdate = true;
  };
  return { set: draw, redraw: () => last && draw(...last) };
})();

// ============================================================ molecule model (atoms + sticks)
const molGroup = new THREE.Group();
scene.add(molGroup);
const sphereGeo = new THREE.SphereGeometry(1, 64, 40);
const stickGeo = new THREE.CylinderGeometry(1, 1, 1, 28, 1, true);
const atomMats = Object.fromEntries(Object.entries(EL).map(([k, v]) => [k, new THREE.MeshPhysicalMaterial({
  color: v.color, roughness: k === 'H' ? 0.18 : 0.28, metalness: 0, clearcoat: 1, clearcoatRoughness: 0.08,
  sheen: 0.4, sheenColor: new THREE.Color(0x9fdcff),
})]));
const ringMat = new THREE.MeshStandardMaterial({ color: 0xe8b85a, metalness: 1, roughness: 0.25, emissive: 0x4a3008, emissiveIntensity: 1.2 });
const partialMat = new THREE.MeshStandardMaterial({ color: 0xe8b85a, metalness: 1, roughness: 0.25, emissive: 0x2a1a04, transparent: true, opacity: 0.75 });

let atomMeshes = [], sticks = [], ringMeshes = [], stickSig = '';

function buildAtoms(sys) {
  for (const m of atomMeshes) molGroup.remove(m);
  atomMeshes = sys.atoms.map((a) => {
    const m = new THREE.Mesh(sphereGeo, atomMats[a.el]);
    m.castShadow = true;
    molGroup.add(m);
    return m;
  });
  for (const r of ringMeshes) molGroup.remove(r);
  ringMeshes = sys.rings.map(() => {
    const r = new THREE.Mesh(new THREE.TorusGeometry(1, 0.022, 12, 120), ringMat);
    molGroup.add(r);
    return r;
  });
  stickSig = '';
}

const orderClass = (o) => (o < 1.15 ? 1 : o < 1.85 ? 1.5 : o < 2.5 ? 2 : 3);

function buildSticks(bonds) {
  for (const s of sticks) molGroup.remove(s.m);
  sticks = [];
  bonds.forEach((b, bi) => {
    const oc = orderClass(b.order);
    const slots = oc === 1 ? [0] : oc === 1.5 ? [0, 'p'] : oc === 2 ? [-1, 1] : [0, 1, 2];
    for (const k of slots) {
      const m = new THREE.Mesh(stickGeo, k === 'p' ? partialMat : gold);
      m.castShadow = true;
      molGroup.add(m);
      sticks.push({ m, bi, k, oc });
    }
  });
}

// ============================================================ caliper (bond-length gauge)
const caliper = new THREE.Group();
scene.add(caliper);
const calLine = new THREE.Mesh(new THREE.CylinderGeometry(0.011, 0.011, 1, 10), goldBright);
const calA = new THREE.Mesh(new THREE.CylinderGeometry(0.011, 0.011, 0.22, 10), goldBright);
const calB = calA.clone();
caliper.add(calLine, calA, calB);

// ============================================================ orbital raymarch pass
const orbitalScene = new THREE.Scene();
const vec4s = (n) => Array.from({ length: n }, () => new THREE.Vector4());
const orbitalUniforms = {
  uC: { value: vec4s(MAX_AO) },
  uB: { value: vec4s(MAX_AO) },
  uN: { value: 0 },
  uAtoms: { value: vec4s(MAX_ATOMS) },
  uNA: { value: 0 },
  uBondA: { value: vec4s(MAX_BONDS) },
  uBondB: { value: Array.from({ length: MAX_BONDS }, () => new THREE.Vector3()) },
  uNB: { value: 0 },
  uBoxMin: { value: new THREE.Vector3() },
  uBoxMax: { value: new THREE.Vector3() },
  uIso: { value: 0.1 },
  uOpacity: { value: 0 },
  uGlow: { value: 1.1 },
  uTime: { value: 0 },
  uScan: { value: 0.35 },
  uPos: { value: COL_POS.clone().multiplyScalar(1.15) },
  uNeg: { value: COL_NEG.clone().multiplyScalar(1.15) },
  uLight: { value: new THREE.Vector3(-0.5, 0.9, 0.6) },
};
const orbitalBox = new THREE.Mesh(
  new THREE.BoxGeometry(1, 1, 1),
  new THREE.ShaderMaterial({ uniforms: orbitalUniforms, vertexShader: orbitalVert, fragmentShader: orbitalFrag, side: THREE.BackSide, depthTest: false, depthWrite: false }),
);
orbitalBox.frustumCulled = false;
orbitalScene.add(orbitalBox);

class OrbitalPass extends Pass {
  constructor() {
    super();
    this.scale = 0.75;
    this.rt = new THREE.WebGLRenderTarget(1, 1, { type: THREE.HalfFloatType });
    this.quad = new FullScreenQuad(new THREE.ShaderMaterial({
      uniforms: { tBase: { value: null }, tOrb: { value: null } },
      vertexShader: 'varying vec2 vUv; void main(){ vUv = uv; gl_Position = vec4(position.xy, 0., 1.); }',
      fragmentShader: compositeFrag, depthTest: false, depthWrite: false,
    }));
    this.w = this.h = 1;
    this._c = new THREE.Color();
  }
  setSize(w, h) { this.w = w; this.h = h; this.applyScale(); }
  setScale(s) { this.scale = s; this.applyScale(); }
  applyScale() { this.rt.setSize(Math.max(1, Math.round(this.w * this.scale)), Math.max(1, Math.round(this.h * this.scale))); }
  render(r, writeBuffer, readBuffer) {
    r.getClearColor(this._c);
    const a = r.getClearAlpha();
    r.setRenderTarget(this.rt);
    r.setClearColor(0x000000, 0);
    r.clear();
    r.render(orbitalScene, camera);
    r.setClearColor(this._c, a);
    this.quad.material.uniforms.tBase.value = readBuffer.texture;
    this.quad.material.uniforms.tOrb.value = this.rt.texture;
    r.setRenderTarget(this.renderToScreen ? null : writeBuffer);
    this.quad.render(r);
  }
}

// ============================================================ electron measurements (point cloud)
const cloudGeo = new THREE.BufferGeometry();
const cloudPos = new Float32Array(CLOUD_N * 3);
const cloudSign = new Float32Array(CLOUD_N);
cloudGeo.setAttribute('position', new THREE.BufferAttribute(cloudPos, 3).setUsage(THREE.DynamicDrawUsage));
cloudGeo.setAttribute('aSign', new THREE.BufferAttribute(cloudSign, 1));
cloudGeo.setAttribute('aSeed', new THREE.BufferAttribute(new Float32Array(CLOUD_N).map(() => Math.random()), 1));
cloudGeo.setAttribute('aIdx', new THREE.BufferAttribute(new Float32Array(CLOUD_N).map((_, i) => i), 1));
cloudGeo.setDrawRange(0, 0);
const cloudUniforms = {
  uTime: { value: 0 }, uSize: { value: 60 }, uFade: { value: 0 },
  uCursor: { value: 0 }, uCount: { value: 1 }, uOcc: { value: 2 },
  uPos: { value: COL_POS.clone() }, uNeg: { value: COL_NEG.clone() },
};
const cloud = new THREE.Points(cloudGeo, new THREE.ShaderMaterial({
  uniforms: cloudUniforms, vertexShader: cloudVert, fragmentShader: cloudFrag,
  transparent: true, depthWrite: false, blending: THREE.AdditiveBlending,
}));
cloud.frustumCulled = false;
scene.add(cloud);

const lightPos = new THREE.PointLight(COL_POS, 0, 9, 2);
const lightNeg = new THREE.PointLight(COL_NEG, 0, 9, 2);
scene.add(lightPos, lightNeg);

// ============================================================ post-processing
const composer = new EffectComposer(renderer, new THREE.WebGLRenderTarget(1, 1, { type: THREE.HalfFloatType, samples: 4 }));
composer.setPixelRatio(renderer.getPixelRatio());
composer.setSize(window.innerWidth, window.innerHeight);
composer.addPass(new RenderPass(scene, camera));
const orbitalPass = new OrbitalPass();
composer.addPass(orbitalPass);
const bloom = new UnrealBloomPass(new THREE.Vector2(window.innerWidth, window.innerHeight), 0.62, 0.55, 0.82);
composer.addPass(bloom);
composer.addPass(new OutputPass());
const finish = new ShaderPass(new THREE.ShaderMaterial({
  uniforms: { tDiffuse: { value: null }, uTime: { value: 0 }, uRes: { value: new THREE.Vector2(window.innerWidth, window.innerHeight) } },
  vertexShader: 'varying vec2 vUv; void main(){ vUv = uv; gl_Position = projectionMatrix * modelViewMatrix * vec4(position,1.); }',
  fragmentShader: finishFrag,
}));
composer.addPass(finish);

// ============================================================ helpers
const sub3 = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dist = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
const centroid = (P) => [0, 1, 2].map((k) => P.reduce((s, p) => s + p[k], 0) / P.length);
const splitFactor = (Re, s) => Math.min(2.2, Math.exp(-1.7 * Re * (s - 1)));
const bondSym = (o) => (o >= 2.5 ? '≡' : o >= 1.85 ? '=' : '–');

function atomNames(atoms) {
  const count = {};
  return atoms.map((a) => {
    if (a.el === 'H') return 'H';
    count[a.el] = (count[a.el] || 0) + 1;
    return `${a.el}${count[a.el]}`;
  });
}

function sdot(a, b, S, n) {
  let s = 0;
  for (let i = 0; i < n; i++) { if (!a[i]) continue; for (let j = 0; j < n; j++) s += a[i] * b[j] * S[i * n + j]; }
  return s;
}

// ============================================================ systems
// A system is anything the bench can show: a gallery molecule, a SMILES build or a reaction.
// It provides geometry(s, ξ), a basis, and solve(pos) → orbitals at that geometry.

function makeGallery(mol, index) {
  const c = centroid(mol.atoms.map((a) => a.p));
  const base = mol.atoms.map((a) => sub3(a.p, c));
  const atoms = mol.atoms.map((a) => ({ el: a.el, charge: 0, label: a.label }));
  const graph = { atoms, bonds: mol.bonds.map(([a, b, order]) => ({ a, b, order, arom: order === 1.5 })) };
  const sCache = new Map();
  let lastS = null, lastT = 0;
  const overlapAt = (s, now) => {
    const key = Math.round(s * 100);
    if (sCache.has(key)) return (lastS = sCache.get(key));
    if (lastS && now - lastT < 110) return lastS;
    const q = key / 100;
    const P = base.map((p) => p.map((v) => v * q));
    lastS = gridOverlap(buildAOs(mol, P), P);
    sCache.set(key, lastS);
    lastT = now;
    return lastS;
  };
  return {
    kind: 'gallery', index, id: mol.id, name: mol.name, formula: mol.formula, blurb: mol.blurb,
    plateSub: `${mol.bondName}  ${mol.morse.Re.toFixed(3)} Å  ·  ORDER ${mol.order}`,
    atoms, graph, basis: mol.basis, rings: aromaticRings(graph), nb: neighbors(graph),
    names: atomNames(atoms),
    geometry: (s) => base.map((p) => p.map((v) => v * s)),
    stretchable: true,
    morse: { ...mol.morse, name: mol.bondName },
    rep: mol.rep, bondName: mol.bondName, order: mol.order,
    cam: mol.cam, iso: mol.iso,
    resonance: resonanceStructures(graph),
    view: { mos: mol.mos, atomic: mol.atomic, capL: mol.atoms[0].el, capR: mol.atoms[1]?.el, capM: mol.formula },
    scaleEnergies: mol.mos.map((m) => m.e),
    defaultMO: mol.defaultMO,
    bonds: () => graph.bonds,
    charges: () => atoms.map(() => 0),
    solve(pos, s, now) {
      const f = splitFactor(mol.morse.Re, s);
      const S = overlapAt(s, now);
      return {
        S, version: S,
        coefs: mol.mos.map((m) => m.c),
        energies: mol.mos.map((m) => m.ref + (m.e - m.ref) * f),
        occ: mol.mos.map((m) => m.occ),
      };
    },
    describe: (i) => mol.mos[i].desc,
  };
}

// shared machinery for extended-Hückel systems (SMILES and reactions)
function ehtCore(atoms, charge) {
  const basis = valenceBasis(atoms);
  const core = { basis, last: null, lastKey: null, version: 0, cost: 0, lastT: 0 };
  core.solveAt = (pos, key, now, force = false) => {
    if (core.lastKey === key && core.last) return core.last;
    if (core.last && !force && now - core.lastT < core.cost * 2.5) return core.last;
    const t0 = performance.now();
    const res = solveEHT(atoms, pos, charge, basis);
    // Hund: electrons spread over a degenerate frontier set
    const E = res.E, occ = res.occ.slice();
    if (res.lumo < E.length && Math.abs(E[res.lumo] - E[res.homo]) < 0.02) {
      let lo = res.homo, hi = res.lumo;
      while (lo > 0 && Math.abs(E[lo - 1] - E[res.homo]) < 0.02) lo--;
      while (hi + 1 < E.length && Math.abs(E[hi + 1] - E[res.homo]) < 0.02) hi++;
      let e = 0;
      for (let k = lo; k <= hi; k++) e += occ[k];
      for (let k = lo; k <= hi; k++) occ[k] = Math.min(2, e / (hi - lo + 1));
    }
    res.occ = occ;
    // keep phases continuous with the previous geometry
    if (core.last) {
      const n = basis.length;
      for (let k = 0; k < res.C.length; k++) {
        if (sdot(core.last.res.C[k], res.C[k], res.S, n) < 0) for (let i = 0; i < n; i++) res.C[k][i] = -res.C[k][i];
      }
    }
    core.version++;
    core.last = { res, S: res.S, version: core.version };
    core.lastKey = key;
    core.cost = performance.now() - t0;
    core.lastT = now;
    return core.last;
  };
  return core;
}

function ehtWindow(res, include = res.homo) {
  const lo = Math.max(0, Math.min(res.homo - 4, include - 1)), hi = Math.min(res.E.length - 1, res.lumo + 3);
  const ks = [];
  for (let k = lo; k <= hi; k++) ks.push(k);
  let g = 0;
  return ks.map((k, i) => {
    if (i > 0 && Math.abs(res.E[k] - res.E[ks[i - 1]]) > 0.03) g++;
    const rel = k - res.homo;
    const somo = res.occ[k] > 0 && res.occ[k] < 2;
    const name = rel === 0 ? (somo ? 'SOMO' : 'HOMO') : rel < 0 ? `HOMO−${-rel}` : rel === 1 ? 'LUMO' : `LUMO+${rel - 1}`;
    return { k, name, occ: res.occ[k], g, tag: rel === 0 ? (somo ? 'SOMO' : 'HOMO') : rel === 1 ? 'LUMO' : null };
  });
}

function describeEHT(sys, res, k) {
  const bl = sys.bonds(state.xi);
  const { share, bondPop, pShare } = analyseMO(res, k, sys.atoms, bl);
  const names = sys.names;
  const idx = [...share.keys()].sort((a, b) => share[b] - share[a]);
  const top = idx.filter((i) => share[i] > 0.07).slice(0, 4);
  const spread = [...share].filter((v) => v > 0.05).length;
  const tot = pShare.s + pShare.x + pShare.y + pShare.z || 1;
  const piFrac = sys.planar ? pShare.y / tot : 0;
  const pops = bondPop.map((p, i) => ({ p, i }));
  const bonding = pops.filter((b) => b.p > 0.04).sort((a, b) => b.p - a.p).slice(0, 3);
  const anti = pops.filter((b) => b.p < -0.06).sort((a, b) => a.p - b.p).slice(0, 3);
  const bname = (b) => `${names[bl[b.i].a]}${bondSym(bl[b.i].order)}${names[bl[b.i].b]}`;
  let type;
  if (piFrac > 0.8) type = anti.length && !bonding.length ? 'π* (antibonding)' : bonding.length && !anti.length ? 'π bonding' : 'π';
  else if (share[idx[0]] > 0.6 && sys.atoms[idx[0]].el !== 'C' && sys.atoms[idx[0]].el !== 'H' && (bonding[0]?.p ?? 0) < 0.15) type = `${names[idx[0]]} lone pair (n)`;
  else if (bonding.length && !anti.length) type = 'σ bonding';
  else if (anti.length && !bonding.length) type = 'σ* antibonding';
  else if (!bonding.length && !anti.length) type = 'Nonbonding';
  else type = 'Mixed σ';
  let text = `<b>${type}</b>` + (spread > 2 ? `, spread over ${spread} atoms` : '') + ': ' +
    top.map((i) => `${names[i]} ${Math.round(share[i] * 100)}%`).join(' · ') + '.';
  if (bonding.length) text += ` In phase across ${bonding.map(bname).join(', ')}.`;
  if (anti.length) text += ` Node across ${anti.map(bname).join(', ')}.`;
  if (sys.frags) {
    const fr = sys.frags.map((f) => f.atoms.reduce((s, i) => s + share[i], 0));
    text += ` <span class="frag">Shared: ${sys.frags.map((f, i) => `${f.name} ${Math.round(Math.max(0, fr[i]) * 100)}%`).join(' · ')}</span>`;
  }
  return text;
}

function makeSmiles(str) {
  const g = addHydrogens(parseSmiles(str));
  if (g.atoms.length > MAX_ATOMS) throw new Error(`${g.atoms.length} atoms: live orbitals are limited to ${MAX_ATOMS}`);
  const model = buildModel(g);
  const raw = orient(model.pos, g.atoms.map((a) => (a.el === 'H' ? 0.3 : 1)));
  const charge = g.atoms.reduce((s, a) => s + a.charge, 0);
  const core = ehtCore(g.atoms, charge);
  const heavy = g.atoms.map((a, i) => i).filter((i) => g.atoms[i].el !== 'H');
  const planar = heavy.length >= 3 && heavy.every((i) => Math.abs(raw[i][1]) < 0.12);
  // representative bond: highest order between heavy atoms, else any bond
  const hb = g.bonds.filter((b) => g.atoms[b.a].el !== 'H' && g.atoms[b.b].el !== 'H');
  const rb = (hb.length ? hb : g.bonds).slice().sort((a, b) => b.order - a.order)[0];
  const names = atomNames(g.atoms);
  const Re = rb ? dist(raw[rb.a], raw[rb.b]) : 1;
  const rbName = rb ? `${names[rb.a]}${bondSym(rb.order)}${names[rb.b]}` : '—';
  const De = rb ? bondEnthalpy(g.atoms[rb.a].el, g.atoms[rb.b].el, rb.order) / 96.485 : 4;
  const res0 = core.solveAt(raw, 's100', 0, true).res;
  const win = ehtWindow(res0);
  const radius = Math.max(...raw.map((p) => Math.hypot(...p)));
  const labelled = g.atoms.map((a, i) => ((a.el !== 'C' && a.el !== 'H') ? i : -1)).filter((i) => i >= 0).slice(0, 6);
  const atoms = g.atoms.map((a, i) => ({ el: a.el, charge: a.charge, label: labelled.includes(i) ? names[i] : null }));
  const formula = formulaOf(g);
  return {
    kind: 'smiles', id: 'smi:' + str, name: formula, formula, smiles: str,
    blurb: `${g.atoms.length} atoms · ${res0.nElec} valence electrons · extended Hückel orbitals`,
    plateSub: str.length > 30 ? str.slice(0, 29) + '…' : str,
    atoms, graph: g, basis: core.basis, rings: aromaticRings(g), nb: neighbors(g), names, planar,
    geometry: (s) => raw.map((p) => p.map((v) => v * s)),
    stretchable: true,
    morse: { De, Re, a: 2.0, name: rbName },
    rep: rb ? [rb.a, rb.b] : [0, 0], bondName: rbName,
    cam: Math.max(7, 5 + radius * 2.6), iso: 0.1,
    resonance: resonanceStructures(g),
    view: { mos: win },
    win,
    scaleEnergies: win.map((m) => res0.E[m.k]),
    defaultMO: win.findIndex((m) => m.k === res0.homo),
    bonds: () => g.bonds,
    charges: () => g.atoms.map((a) => a.charge),
    solve(pos, s, now) {
      const r = core.solveAt(pos, 's' + Math.round(s * 100), now);
      return { S: r.S, res: r.res, coefs: win.map((m) => r.res.C[m.k]), energies: win.map((m) => r.res.E[m.k]), occ: win.map((m) => r.res.occ[m.k]), version: r.version };
    },
  };
}

function makeReaction(def, index) {
  const rx = buildReaction(def);
  const core = ehtCore(rx.atoms, rx.charge);
  const names = atomNames(rx.atoms);
  const { comp, count } = rx.fragR;
  const frags = [];
  for (let c = 0; c < count; c++) {
    const idx = comp.map((x, i) => (x === c ? i : -1)).filter((i) => i >= 0);
    frags.push({ atoms: idx, name: formulaOf({ atoms: idx.map((i) => rx.atoms[i]) }) });
  }
  const res0 = core.solveAt(pathAt(rx, 0), 'x0', 0, true).res;
  const scale = [0, 0.5, 1].flatMap((xi) => {
    const r = solveEHT(rx.atoms, pathAt(rx, xi), rx.charge, core.basis);
    return ehtWindow(r).map((m) => r.E[m.k]);
  });
  // start by following the donor fragment's highest occupied orbital (e.g. hydroxide's lone pair in SN2)
  const donor = frags[def.donor ?? 0];
  let startK = res0.homo;
  for (let k = res0.homo; k >= 0; k--) {
    const { share } = analyseMO(res0, k, rx.atoms, []);
    if (donor.atoms.reduce((t, i) => t + share[i], 0) > 0.6) { startK = k; break; }
  }
  const win = ehtWindow(res0, startK);
  const graphU = { atoms: rx.atoms, bonds: rx.bondsUnion };
  const basisByAtom = rx.atoms.map((_, a) => core.basis.map((b, i) => (b.atom === a ? i : -1)).filter((i) => i >= 0));
  const donorBasis = frags[def.donor ?? 0].atoms.flatMap((a) => basisByAtom[a]);
  const radius = Math.max(...['R', 'TS', 'P'].flatMap((f) => rx.frames[f].map((p) => Math.hypot(...p))));
  const rep = rx.forming[0] || rx.breaking[0];
  const changing = new Set([...rx.forming, ...rx.breaking].flat());
  const atoms = rx.atoms.map((a, i) => ({ el: a.el, charge: a.charge, label: changing.has(i) && a.el !== 'H' ? names[i] : null }));
  // bond orders switch at the TS; presence fades with distance
  const bondsAt = (xi) => rx.bondsUnion.map((b) => ({ a: b.a, b: b.b, order: xi < 0.5 ? b.rOrder || b.pOrder : b.pOrder || b.rOrder, arom: false }));
  return {
    kind: 'reaction', index, id: 'rx:' + def.id, basisByAtom, donorBasis, name: def.title, formula: def.kind, def, rx,
    blurb: def.kind,
    plateSub: def.kind.toUpperCase(),
    atoms, graph: graphU, basis: core.basis, rings: [], nb: neighbors(graphU), names, planar: false, frags,
    geometry: (s, xi) => pathAt(rx, xi),
    stretchable: false,
    rep, bondName: rep ? `${names[rep[0]]}···${names[rep[1]]}` : '',
    cam: Math.max(7.5, 4.6 + radius * 2.1), iso: def.iso ?? 0.1,
    resonance: null,
    view: { mos: win },
    win,
    scaleEnergies: scale,
    defaultMO: Math.max(0, win.findIndex((m) => m.k === startK)),
    bonds: bondsAt,
    presence: true,
    charges: (xi) => (xi < 0.5 ? rx.atoms.map((a) => a.charge) : rx.atomsP.map((a) => a.charge)),
    solve(pos, s, now) {
      const r = core.solveAt(pos, 'x' + Math.round(state.xi * 400), now);
      return { S: r.S, res: r.res, coefs: win.map((m) => r.res.C[m.k]), energies: win.map((m) => r.res.E[m.k]), occ: win.map((m) => r.res.occ[m.k]), version: r.version };
    },
  };
}

// ============================================================ state
const state = {
  mode: 'gallery',
  sys: null,
  sol: null,
  sel: 0,
  fromC: null,
  blend: 1,
  shown: null,          // coefficients currently on screen (target of the last frame)
  track: null,          // reaction mode: the orbital being followed along the path
  stretch: 1, stretchV: 0, dragS: false, lastInput: 0,
  xi: 0, playing: false, playHold: 0,
  iso: 0.1,
  view: 'surface',
  vibrate: false, vibAmp: 0,
  pop: 0, orbAlpha: 0, surfAmt: 1, cloudAmt: 0, cloudFade: 0,
  cloudData: null, cloudKey: null, geoChangedAt: 0,
  measure: { cursor: 0, t: 0 },
  resIndex: -1, resCycle: false, resT: 0,
  switching: false, pending: null,
  galleryIndex: 0, rxIndex: 0,
};

// ============================================================ tweens
const tweens = [];
const easeInOut = (t) => (t < 0.5 ? 4 * t * t * t : 1 - (-2 * t + 2) ** 3 / 2);
const easeOutBack = (t) => 1 + 2.2 * (t - 1) ** 3 + 1.2 * (t - 1) ** 2;
function tween(dur, fn, ease = easeInOut) {
  return new Promise((resolve) => tweens.push({ t: 0, dur, fn, ease, resolve }));
}
function runTweens(dt) {
  for (let i = tweens.length - 1; i >= 0; i--) {
    const tw = tweens[i];
    tw.t = Math.min(1, tw.t + dt / tw.dur);
    tw.fn(tw.ease(tw.t));
    if (tw.t >= 1) { tweens.splice(i, 1); tw.resolve(); }
  }
}

// ============================================================ UI wiring
const diagram = new MODiagram($('mo'), (i) => selectMO(i));
const morse = new MorseChart($('morse'), {
  onStart: () => (state.dragS = true),
  onDrag: (s) => setStretch(s),
  onEnd: () => (state.dragS = false),
});
const profile = new ProfileChart($('profile'), {
  onStart: () => setPlaying(false),
  onDrag: (xi) => setXi(xi),
});

MOLECULES.forEach((m, i) => {
  const b = document.createElement('button');
  b.textContent = m.formula;
  b.title = `${m.name} (${i + 1})`;
  b.addEventListener('click', () => loadGallery(i));
  $('molSeg').appendChild(b);
});
function loadGallery(i) {
  state.galleryIndex = i;
  [...$('molSeg').children].forEach((b, k) => b.classList.toggle('on', k === i));
  loadSystem(() => makeGallery(MOLECULES[i], i));
}

const EXAMPLES = [
  ['Pyridine', 'c1ccncc1'], ['Acetate', 'CC(=O)[O-]'], ['Formaldehyde', 'C=O'], ['Ethanol', 'CCO'],
  ['Acetylene', 'C#C'], ['HCN', 'C#N'], ['CO₂', 'O=C=O'], ['Butadiene', 'C=CC=C'],
  ['Naphthalene', 'c1ccc2ccccc2c1'], ['Allyl cation', 'C=C[CH2+]'], ['Aspirin', 'CC(=O)Oc1ccccc1C(=O)O'], ['Caffeine', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C'],
];
for (const [label, smi] of EXAMPLES) {
  const b = document.createElement('button');
  b.textContent = label;
  b.title = smi;
  b.addEventListener('click', () => { $('smiIn').value = smi; buildSmiles(); });
  $('smiEx').appendChild(b);
}
$('smiForm').addEventListener('submit', (e) => { e.preventDefault(); buildSmiles(); });

function buildSmiles() {
  const str = $('smiIn').value.trim();
  const msg = $('smiMsg');
  msg.className = 'msg';
  msg.textContent = 'Building geometry and solving orbitals…';
  setTimeout(() => {
    let sys;
    try { sys = makeSmiles(str); }
    catch (err) { msg.className = 'msg err'; msg.textContent = err.message; return; }
    msg.textContent = `${sys.formula} · ${sys.atoms.length} atoms · ${sys.basis.length} basis orbitals${sys.resonance ? ` · ${sys.resonance.structures.length} resonance structures` : ''}`;
    loadSystem(() => sys);
  }, 20);
}

REACTIONS.forEach((r, i) => {
  const b = document.createElement('button');
  b.innerHTML = `${r.title}<small>${r.kind}</small>`;
  b.addEventListener('click', () => loadReaction(i));
  $('rxSeg').appendChild(b);
});
function loadReaction(i) {
  state.rxIndex = i;
  [...$('rxSeg').children].forEach((b, k) => b.classList.toggle('on', k === i));
  loadSystem(() => makeReaction(REACTIONS[i], i));
}

$('modeSeg').addEventListener('click', (e) => {
  const b = e.target.closest('button');
  if (b) setMode(b.dataset.m);
});
function setMode(m) {
  if (m === state.mode) return;
  state.mode = m;
  for (const b of $('modeSeg').children) b.classList.toggle('on', b.dataset.m === m);
  $('paneGallery').hidden = m !== 'gallery';
  $('paneSmiles').hidden = m !== 'smiles';
  $('paneReaction').hidden = m !== 'reaction';
  $('stretchBox').hidden = m === 'reaction';
  $('sliders').classList.toggle('single', m === 'reaction');
  $('morse').toggleAttribute('hidden', m === 'reaction');
  $('profile').toggleAttribute('hidden', m !== 'reaction');
  setPlaying(false);
  if (m === 'gallery') loadGallery(state.galleryIndex);
  else if (m === 'smiles') buildSmiles();
  else loadReaction(state.rxIndex);
}

const stretchEl = $('stretch'), isoEl = $('iso'), xiEl = $('xi');
const setRangeFill = (el) => el.style.setProperty('--p', `${((el.value - el.min) / (el.max - el.min)) * 100}%`);
stretchEl.addEventListener('pointerdown', () => (state.dragS = true));
window.addEventListener('pointerup', () => (state.dragS = false));
stretchEl.addEventListener('input', () => { state.lastInput = performance.now(); setStretch(stretchEl.value / 100, true); });
const isoFromSlider = (v) => 0.015 * 15 ** (v / 100);
const sliderFromIso = (iso) => (Math.log(iso / 0.015) / Math.log(15)) * 100;
isoEl.addEventListener('input', () => {
  state.iso = isoFromSlider(+isoEl.value);
  $('isoVal').textContent = state.iso.toFixed(3);
  setRangeFill(isoEl);
});
xiEl.addEventListener('input', () => { setPlaying(false); setXi(xiEl.value / 1000, true); });
$('btnPlay').addEventListener('click', () => setPlaying(!state.playing));

function setPlaying(p) {
  state.playing = p;
  state.playHold = 0;
  $('btnPlay').classList.toggle('on', p);
  if (p && state.xi >= 0.999) setXi(0);
}
function setXi(xi, fromSlider = false) {
  state.xi = Math.max(0, Math.min(1, xi));
  if (!fromSlider) xiEl.value = Math.round(state.xi * 1000);
  setRangeFill(xiEl);
  state.geoChangedAt = performance.now();
}
function setStretch(s, fromSlider = false) {
  state.stretch = Math.max(S_MIN, Math.min(S_MAX, s));
  state.stretchV = 0;
  if (!fromSlider) stretchEl.value = Math.round(state.stretch * 100);
  setRangeFill(stretchEl);
  state.geoChangedAt = performance.now();
}

$('viewSeg').addEventListener('click', (e) => {
  const b = e.target.closest('button');
  if (b) setView(b.dataset.v);
});
function setView(v) {
  state.view = v;
  for (const b of $('viewSeg').children) b.classList.toggle('on', b.dataset.v === v);
  const surf = v === 'surface' ? 1 : v === 'both' ? 0.5 : 0;
  const cl = v === 'surface' ? 0 : 1;
  const s0 = state.surfAmt, c0 = state.cloudAmt;
  if (cl && c0 < 0.01) resetMeasure();
  tween(0.6, (t) => { state.surfAmt = s0 + (surf - s0) * t; state.cloudAmt = c0 + (cl - c0) * t; });
}

const btnRotate = $('btnRotate'), btnVibrate = $('btnVibrate');
btnRotate.addEventListener('click', () => toggleRotate());
btnVibrate.addEventListener('click', () => toggleVibrate());
function toggleRotate() { controls.autoRotate = !controls.autoRotate; btnRotate.classList.toggle('on', controls.autoRotate); }
function toggleVibrate() { state.vibrate = !state.vibrate; btnVibrate.classList.toggle('on', state.vibrate); }

const help = $('help');
$('btnHelp').addEventListener('click', () => (help.hidden = false));
$('helpClose').addEventListener('click', () => (help.hidden = true));
help.addEventListener('click', (e) => { if (e.target === help) help.hidden = true; });

window.addEventListener('keydown', (e) => {
  if (e.target.id === 'smiIn') return;
  if (e.target.tagName === 'INPUT' && !e.key.startsWith('Arrow')) return;
  const n = parseInt(e.key, 10);
  if (n >= 1 && n <= 6) {
    if (state.mode === 'gallery' && MOLECULES[n - 1]) loadGallery(n - 1);
    if (state.mode === 'reaction' && REACTIONS[n - 1]) loadReaction(n - 1);
  } else if (e.key === 'ArrowRight' || e.key === 'ArrowUp') { e.preventDefault(); selectMO(state.sel + 1); }
  else if (e.key === 'ArrowLeft' || e.key === 'ArrowDown') { e.preventDefault(); selectMO(state.sel - 1); }
  else if (e.key === 'v' || e.key === 'V') setView({ surface: 'both', both: 'cloud', cloud: 'surface' }[state.view]);
  else if (e.key === 'r' || e.key === 'R') toggleRotate();
  else if (e.key === ' ') { e.preventDefault(); state.mode === 'reaction' ? setPlaying(!state.playing) : toggleVibrate(); }
  else if (e.key === '?' || e.key === 'h' || e.key === 'H') help.hidden = !help.hidden;
  else if (e.key === 'Escape') help.hidden = true;
});

// ---------- resonance
function setupResonance(sys) {
  const row = $('resRow'), seg = $('resSeg');
  seg.replaceChildren();
  state.resCycle = false;
  $('btnResCycle').classList.remove('on');
  const R = sys.resonance;
  row.hidden = !R;
  if (!R) { state.resIndex = -1; return; }
  const roman = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII'];
  R.structures.forEach((_, i) => {
    const b = document.createElement('button');
    b.textContent = roman[i];
    b.title = `Resonance structure ${roman[i]}`;
    b.addEventListener('click', () => setResonance(i));
    seg.appendChild(b);
  });
  const h = document.createElement('button');
  h.textContent = 'Hybrid';
  h.title = 'The real molecule: an average of the structures';
  h.addEventListener('click', () => setResonance(-1));
  seg.appendChild(h);
  setResonance(-1, true);
}
function setResonance(i, keepCycle = false) {
  state.resIndex = i;
  if (!keepCycle) { state.resCycle = false; $('btnResCycle').classList.remove('on'); }
  const btns = [...$('resSeg').children];
  btns.forEach((b, k) => b.classList.toggle('on', i === -1 ? k === btns.length - 1 : k === i));
}
$('btnResCycle').addEventListener('click', () => {
  state.resCycle = !state.resCycle;
  state.resT = 0;
  $('btnResCycle').classList.toggle('on', state.resCycle);
  if (state.resCycle) setResonance(0, true);
});

function displayBonds(sys) {
  const base = sys.bonds(state.xi);
  const R = sys.resonance;
  if (!R) return base;
  const orders = state.resIndex < 0 ? R.hybrid : R.structures[state.resIndex].orders;
  return base.map((b, k) => ({ ...b, order: orders[k] }));
}
function displayCharges(sys) {
  const R = sys.resonance;
  if (!R) return sys.charges(state.xi);
  return state.resIndex < 0 ? R.hybridCharges : R.structures[state.resIndex].charges;
}

// ============================================================ labels
const labelsEl = $('labels');
let labels = [], chargeLabels = [], labelPos, labelNeg, labelBond;
function makeLabel(text, cls = '') {
  const d = document.createElement('div');
  d.className = 'tag3d ' + cls;
  d.innerHTML = `<i></i><span>${text}</span>`;
  labelsEl.appendChild(d);
  return { el: d, span: d.querySelector('span'), pos: new THREE.Vector3(), dx: 12, dy: -12, visible: true };
}
function buildLabels(sys) {
  labelsEl.replaceChildren();
  labels = [];
  sys.atoms.forEach((a, i) => {
    if (!a.label) return;
    const l = makeLabel(a.label);
    l.atom = i;
    labels.push(l);
  });
  chargeLabels = sys.atoms.map((_, i) => {
    const l = makeLabel('', 'charge');
    l.atom = i;
    l.charge = true;
    l.dx = -14; l.dy = -18;
    l.visible = false;
    labels.push(l);
    return l;
  });
  labelPos = makeLabel('ψ &gt; 0', 'pos');
  labelNeg = makeLabel('ψ &lt; 0', 'neg');
  labelBond = makeLabel('', 'gold');
  labelBond.dx = 0; labelBond.dy = 0;
  labels.push(labelPos, labelNeg, labelBond);
}
const _v = new THREE.Vector3();
function placeLabels() {
  const w = window.innerWidth, h = window.innerHeight;
  for (const l of labels) {
    _v.copy(l.pos).project(camera);
    const on = l.visible && _v.z < 1 && Math.abs(_v.x) < 1.1 && Math.abs(_v.y) < 1.1;
    l.el.style.opacity = on ? 1 : 0;
    if (!on) continue;
    const x = (_v.x * 0.5 + 0.5) * w + l.dx, y = (-_v.y * 0.5 + 0.5) * h + l.dy;
    l.el.style.transform = `translate(${x.toFixed(1)}px, ${y.toFixed(1)}px) translate(${l.dx === 0 ? '-50%' : '0'}, -50%)`;
  }
}
const chargeText = (q) => {
  const a = Math.abs(q);
  const mag = Math.abs(a - 1) < 0.01 ? '' : Math.abs(a - 0.5) < 0.01 ? '½' : Math.abs(a - 1 / 3) < 0.01 ? '⅓' : a.toFixed(1);
  return (/^\d/.test(mag) ? 'δ' : '') + mag + (q > 0 ? '+' : '−');
};

// ============================================================ orbital selection
function selectMO(i, instant = false) {
  const sys = state.sys;
  if (!sys || !state.sol) return;
  const n = sys.view.mos.length;
  i = (i + n) % n;
  if (i === state.sel && !instant) return;
  const cur = state.shown;
  state.sel = i;
  if (instant || !cur) { state.fromC = null; state.blend = 1; }
  else {
    state.fromC = Float64Array.from(cur);
    state.blend = 0;
    tween(0.75, (t) => (state.blend = t));
  }
  if (sys.kind === 'reaction') {
    state.track = Float64Array.from(state.sol.coefs[i]);
    state.trackK = sys.win[i].k;
    if (!instant) state.trackMode = 'manual';
  }
  diagram.setSelected(i);
  refreshMOText(true);
  state.cloudKey = null;
  resetMeasure();
}

let descT = 0;
function refreshMOText(force = false) {
  const sys = state.sys, sol = state.sol;
  if (!sys || !sol) return;
  const now = performance.now();
  if (!force && now - descT < 250) return;
  descT = now;
  const m = sys.view.mos[state.sel];
  const occ = sol.occ[state.sel];
  $('moName').innerHTML = richHTML(m.name) + (sys.kind !== 'gallery' ? ` <small>${fmt(sol.energies[state.sel], 1)} eV</small>` : '');
  const el = $('moOcc');
  el.textContent = occ >= 2 ? '2 e⁻ · filled' : occ > 0 ? `${+occ.toFixed(2)} e⁻ · half-filled` : 'empty';
  el.classList.toggle('empty', !occ);
  $('moDesc').innerHTML = sys.kind === 'gallery' ? richHTML(sys.describe(state.sel)) : describeEHT(sys, sol.res, sys.win[state.sel].k);
}

// Reaction tracking. Auto mode follows the "bond-forming orbital": the occupied orbital with the
// largest overlap population across the forming bond(s). Early on, when nothing overlaps yet,
// or after the user picks a level by hand, it follows by continuity (max overlap with the last one).
function formingPop(res, c, sys) {
  const n = res.C[0].length, S = res.S, B = sys.basisByAtom;
  let p = 0;
  for (const [a, b] of sys.rx.forming) for (const i of B[a]) for (const j of B[b]) p += 2 * c[i] * c[j] * S[i * n + j];
  return p;
}
function followTrack(sol) {
  if (!state.track) return;
  const sys = state.sys, res = sol.res, n = res.C[0].length;
  let bk = -1;
  if (state.trackMode === 'auto') {
    // among orbitals still holding ≥25% of the donor's electrons, the one that bonds most across the forming bond
    let best = -Infinity;
    const D = sys.donorBasis;
    for (let k = 0; k < res.C.length; k++) {
      if (!(res.occ[k] > 0)) continue;
      const c = res.C[k];
      let share = 0;
      for (const i of D) { let g = 0; for (let j = 0; j < n; j++) g += res.S[i * n + j] * c[j]; share += c[i] * g; }
      if (share < 0.25) continue;
      const p = formingPop(res, c, sys) + 0.01 * res.E[k] / 10;
      if (p > best) { best = p; bk = k; }
    }
  }
  if (bk < 0) {
    let best = 0;
    for (let k = 0; k < res.C.length; k++) {
      const o = Math.abs(sdot(state.track, res.C[k], res.S, n));
      if (o > best) { best = o; bk = k; }
    }
  }
  const c = Float64Array.from(res.C[bk]);
  if (sdot(state.track, c, res.S, n) < 0) for (let i = 0; i < n; i++) c[i] = -c[i];
  if (bk !== state.trackK && state.trackK !== undefined && state.shown) {
    // the followed orbital changed character: morph instead of jumping
    state.fromC = Float64Array.from(state.shown);
    state.blend = 0;
    tween(0.45, (t) => (state.blend = t));
  }
  state.trackK = bk;
  state.track = c;
  const w = sys.win.findIndex((m) => m.k === bk);
  if (w >= 0 && w !== state.sel) { state.sel = w; diagram.setSelected(w); state.cloudKey = null; refreshMOText(true); }
}

function resetMeasure() {
  state.measure.cursor = 0;
  state.measure.t = 0;
}

// ============================================================ loading systems
async function loadSystem(factory) {
  if (state.switching) { state.pending = factory; return; }
  state.switching = true;
  const first = !state.sys;
  if (!first) {
    const a0 = state.orbAlpha, p0 = state.pop, c0 = state.cloudFade;
    await tween(0.35, (t) => { state.orbAlpha = a0 * (1 - t); state.pop = p0 * (1 - t); state.cloudFade = c0 * (1 - t); });
  }
  let sys;
  try { sys = factory(); }
  catch (err) {
    console.error(err);
    $('smiMsg').className = 'msg err';
    $('smiMsg').textContent = err.message;
    state.switching = false;
    if (state.sys) { state.orbAlpha = 1; state.pop = 1; }
    return;
  }
  state.sys = sys;
  state.sol = null;
  state.track = null;
  state.trackK = undefined;
  state.trackMode = 'auto';
  state.shown = null;
  state.stretch = 1; state.stretchV = 0;
  stretchEl.value = 100;
  setRangeFill(stretchEl);
  state.xi = 0;
  xiEl.value = 0;
  setRangeFill(xiEl);
  state.iso = sys.iso;
  isoEl.value = sliderFromIso(sys.iso);
  $('isoVal').textContent = sys.iso.toFixed(3);
  setRangeFill(isoEl);
  orbitalUniforms.uIso.value = sys.iso;

  // keep big molecules clear of the stage
  const P0 = sys.kind === 'reaction' ? ['R', 'TS', 'P'].flatMap((f) => sys.rx.frames[f]) : sys.geometry(1, 0);
  const minY = Math.min(...P0.map((p) => p[1]));
  sys.lift = Math.max(0, STAGE_Y + 1.5 - minY);

  buildAtoms(sys);
  buildLabels(sys);
  plate.set(sys.kind === 'gallery' ? `${sys.formula}  ·  ${sys.name.toUpperCase()}` : sys.kind === 'reaction' ? sys.name : sys.formula, sys.plateSub);
  diagram.setMolecule(sys.view, sys.scaleEnergies);
  setupResonance(sys);

  $('molName').textContent = sys.name;
  $('molBlurb').textContent = sys.blurb;
  $('eyebrow').textContent = { gallery: 'Molecular orbital bench', smiles: 'SMILES builder · extended Hückel', reaction: 'Reaction lab' }[sys.kind];
  $('diagCap').textContent = sys.kind === 'reaction' ? 'Orbitals along the path' : 'Orbital diagram';
  if (sys.kind === 'reaction') {
    profile.setProfile(sys.def, (xi) => profileAt(sys.def, xi));
    $('leftCap').innerHTML = `Reaction profile · <em>${sys.def.kind}</em>`;
  } else {
    morse.setMorse(sys.morse, S_MIN, S_MAX);
    $('leftCap').innerHTML = `Potential energy · <em>${sys.morse.name}</em>`;
    $('leftFoot').textContent = 'Drag the curve to pull the bond. Let go and it springs back to Rₑ.';
  }

  // prime the orbitals so the first frame has something to draw
  state.sol = sys.solve(sys.geometry(1, 0), 1, 0);
  state.sel = -1;
  selectMO(Math.max(0, sys.defaultMO), true);

  // frame the molecule
  const dir = camera.position.clone().sub(controls.target).normalize();
  const from = camera.position.clone();
  const target = new THREE.Vector3(0, -0.15 + sys.lift, 0);
  const t0 = controls.target.clone();
  if (first) {
    controls.target.copy(target);
    camera.position.copy(target).add(new THREE.Vector3(0.9, 0.62, 1.6).normalize().multiplyScalar(sys.cam * 2.6));
    const f = camera.position.clone();
    const to = target.clone().add(new THREE.Vector3(0.62, 0.36, 1).normalize().multiplyScalar(sys.cam));
    tween(2.6, (t) => camera.position.lerpVectors(f, to, t));
  } else {
    const to = target.clone().add(dir.multiplyScalar(sys.cam));
    tween(1.1, (t) => { camera.position.lerpVectors(from, to, t); controls.target.lerpVectors(t0, target, t); });
  }
  lastUI = '';
  await Promise.all([
    tween(first ? 1.4 : 0.7, (t) => (state.pop = t), easeOutBack),
    tween(first ? 1.8 : 0.8, (t) => (state.orbAlpha = t)),
  ]);
  state.switching = false;
  if (sys.kind === 'reaction' && !state.pending) setPlaying(true);
  if (state.pending) {
    const p = state.pending;
    state.pending = null;
    loadSystem(p);
  }
}

function resampleCloud(pos, sol, c) {
  const aos = buildAOs(state.sys, pos);
  const nrm = normOf(c, sol.S);
  const w = Array.from(c, (x) => x / nrm);
  const data = sampleCloud(aos, w, pos, CLOUD_N);
  state.cloudData = data;
  cloudSign.set(data.sign.subarray(0, data.count));
  cloudGeo.attributes.aSign.needsUpdate = true;
  cloudGeo.setDrawRange(0, data.count);
  cloudUniforms.uCount.value = data.count;
  if (state.cloudFade < 1) tween(0.6, (t) => (state.cloudFade = Math.max(state.cloudFade, t)));
}

// ============================================================ per-frame update
const up = new THREE.Vector3(0, 1, 0);
const zAxis = new THREE.Vector3(0, 0, 1);
const tA = new THREE.Vector3(), tB = new THREE.Vector3(), tD = new THREE.Vector3(), tU = new THREE.Vector3(), tW = new THREE.Vector3(), tC = new THREE.Vector3();
const lobeP = new THREE.Vector3(), lobeN = new THREE.Vector3();
let lastUI = '';

function update(time, dt, now) {
  const sys = state.sys;
  if (!sys || !state.sol) return;

  // ---- bond physics: the bond springs back to Rₑ unless you hold it
  const holding = state.dragS || now - state.lastInput < 250;
  if (sys.stretchable && !holding && (Math.abs(state.stretch - 1) > 1e-4 || Math.abs(state.stretchV) > 1e-4)) {
    const acc = -55 * (state.stretch - 1) - 3.2 * state.stretchV;
    state.stretchV += acc * dt;
    state.stretch += state.stretchV * dt;
    if (Math.abs(state.stretch - 1) < 1e-4 && Math.abs(state.stretchV) < 1e-3) { state.stretch = 1; state.stretchV = 0; }
    stretchEl.value = Math.round(state.stretch * 100);
    setRangeFill(stretchEl);
    state.geoChangedAt = now;
  }
  // ---- reaction playback
  if (sys.kind === 'reaction' && state.playing && !state.switching) {
    if (state.xi >= 1) {
      state.playHold += dt;
      if (state.playHold > 2.4) { state.playHold = 0; setXi(0); }
    } else setXi(state.xi + dt / 7);
  }
  // ---- resonance auto-cycle
  if (state.resCycle && sys.resonance) {
    state.resT += dt;
    if (state.resT > 1.6) {
      state.resT = 0;
      const n = sys.resonance.structures.length;
      setResonance(state.resIndex === -1 ? 0 : state.resIndex + 1 >= n ? -1 : state.resIndex + 1, true);
    }
  }

  state.vibAmp += ((state.vibrate && sys.stretchable ? 0.075 : 0) - state.vibAmp) * 0.06;
  const s = sys.stretchable ? state.stretch * (1 + state.vibAmp * Math.sin(time * Math.PI * 2 * 1.15)) : 1;
  const bob = Math.sin(time * 0.8) * 0.05;
  const raw = sys.geometry(s, state.xi);
  const lift = sys.lift + bob;
  const pos = raw.map((p) => [p[0], p[1] + lift, p[2]]);

  // ---- orbitals at this geometry
  const sol = sys.solve(raw, s, now);
  const fresh = sol.version !== state.sol.version;
  state.sol = sol;
  if (fresh && sys.kind === 'reaction') followTrack(sol);
  const target = sys.kind === 'reaction' && state.track ? state.track : sol.coefs[state.sel];
  let c = target;
  if (state.blend < 1 && state.fromC && state.fromC.length === target.length) {
    c = new Float64Array(target.length);
    for (let i = 0; i < c.length; i++) c[i] = state.fromC[i] * (1 - state.blend) + target[i] * state.blend;
  }
  state.shown = c;
  const nrm = normOf(c, sol.S);

  // ---- atoms
  const pop = Math.max(0, state.pop);
  atomMeshes.forEach((m, i) => {
    m.position.fromArray(pos[i]);
    m.scale.setScalar(EL[sys.atoms[i].el].vis * pop);
  });

  // ---- sticks (rebuilt when the bond pattern changes: resonance, reactions)
  const bl = displayBonds(sys);
  const sig = bl.map((b) => orderClass(b.order)).join('');
  if (sig !== stickSig) { buildSticks(bl); stickSig = sig; }
  const presence = bl.map((b) => {
    if (!sys.presence) return 1;
    const d0 = bondLength(sys.atoms[b.a].el, sys.atoms[b.b].el, b.order, 4, 4);
    const d = dist(pos[b.a], pos[b.b]);
    return Math.max(0, Math.min(1, (1.75 * d0 - d) / (0.55 * d0)));
  });
  for (const st of sticks) {
    const b = bl[st.bi];
    tA.fromArray(pos[b.a]);
    tB.fromArray(pos[b.b]);
    tD.subVectors(tB, tA);
    const len = tD.length();
    tD.normalize();
    const pr = presence[st.bi] * pop;
    let off = 0, r = 0.058;
    if (st.oc > 1) {
      // side direction: toward a ring centre, else toward a neighbouring substituent (in the σ plane)
      const ring = sys.rings.find((Rg) => Rg.includes(b.a) && Rg.includes(b.b));
      if (ring) tU.fromArray(centroid(ring.map((i) => pos[i]))).sub(tA);
      else {
        const nbr = sys.nb[b.a].find(({ j }) => j !== b.b) || sys.nb[b.b].find(({ j }) => j !== b.a);
        if (nbr) tU.fromArray(pos[nbr.j]).sub(tA); else tU.set(0, 1, 0);
      }
      tU.addScaledVector(tD, -tU.dot(tD));
      if (tU.lengthSq() < 1e-6) tU.set(0, 1, 0).addScaledVector(tD, -tD.y);
      if (tU.lengthSq() < 1e-6) tU.set(0, 0, 1);
      tU.normalize();
      tW.crossVectors(tD, tU);
      if (st.oc === 1.5) {
        if (st.k === 'p') { off = 0.13; r = 0.022 * (0.5 + Math.min(1, Math.abs(b.order - 1) * 2) * 0.5); }
      } else if (st.oc === 2) {
        r = 0.042;
        if (ring) off = st.k === 1 ? 0.14 : 0;
        else { off = 0.075; tU.multiplyScalar(st.k); }
      } else {
        const a = (st.k / 3) * Math.PI * 2 + Math.PI / 2;
        tU.multiplyScalar(Math.cos(a)).addScaledVector(tW, Math.sin(a));
        off = 0.1;
        r = 0.04;
      }
    }
    const shorten = off > 0.12 ? 0.72 : 1;
    st.m.position.addVectors(tA, tB).multiplyScalar(0.5);
    if (off) st.m.position.addScaledVector(tU, off);
    st.m.quaternion.setFromUnitVectors(up, tD);
    st.m.scale.set(r * pr, len * shorten * Math.min(1, pop * 1.3), r * pr);
    st.m.visible = pr > 0.02;
  }
  // aromatic rings: the hybrid's delocalised π, hidden while a single Kekulé structure is shown
  const showRings = !sys.resonance || state.resIndex < 0;
  sys.rings.forEach((Rg, i) => {
    const m = ringMeshes[i];
    const P = Rg.map((k) => pos[k]);
    const c0 = centroid(P);
    tA.fromArray(sub3(P[0], c0));
    tB.fromArray(sub3(P[1], c0));
    tC.crossVectors(tA, tB).normalize();
    m.position.fromArray(c0);
    m.quaternion.setFromUnitVectors(zAxis, tC);
    const rad = P.reduce((q, p) => q + dist(p, c0), 0) / P.length;
    m.scale.setScalar(rad * 0.62 * pop);
    m.visible = showRings;
  });

  // ---- caliper on the representative bond
  const [ri, rj] = sys.rep || [0, 0];
  const R = dist(pos[ri], pos[rj]);
  caliper.visible = pop > 0.5 && ri !== rj;
  if (caliper.visible) {
    tA.fromArray(pos[ri]);
    tB.fromArray(pos[rj]);
    tD.subVectors(tB, tA).normalize();
    tU.set(0, -1, 0).addScaledVector(tD, tD.y);
    if (tU.lengthSq() < 1e-3) tU.set(0, 0, 1);
    tU.normalize().multiplyScalar(0.95);
    tA.add(tU); tB.add(tU);
    calLine.position.addVectors(tA, tB).multiplyScalar(0.5);
    calLine.quaternion.setFromUnitVectors(up, tD);
    calLine.scale.set(1, R, 1);
    tW.copy(tU).normalize();
    calA.position.copy(tA); calA.quaternion.setFromUnitVectors(up, tW);
    calB.position.copy(tB); calB.quaternion.setFromUnitVectors(up, tW);
    labelBond.pos.copy(calLine.position).addScaledVector(tW, 0.28);
    labelBond.span.textContent = `${R.toFixed(2)} Å`;
  }
  labelBond.visible = caliper.visible;

  // ---- raymarch uniforms (only the AOs that matter for this orbital)
  const aos = buildAOs(sys, pos);
  const U = orbitalUniforms;
  const wts = aos.map((_, i) => c[i] / nrm);
  const wmax = Math.max(1e-9, ...wts.map(Math.abs));
  const keep = wts.map((_, i) => i).filter((i) => Math.abs(wts[i]) > wmax * 0.02).sort((a, b) => Math.abs(wts[b]) - Math.abs(wts[a])).slice(0, MAX_AO);
  keep.forEach((i, n) => {
    const a = aos[i], w = wts[i];
    U.uC.value[n].set(a.cx, a.cy, a.cz, a.z + 100 * a.k);
    U.uB.value[n].set(a.bx * w, a.by * w, a.bz * w, a.A * w);
  });
  U.uN.value = keep.length;
  const NA = Math.min(pos.length, MAX_ATOMS);
  for (let i = 0; i < NA; i++) U.uAtoms.value[i].set(pos[i][0], pos[i][1], pos[i][2], EL[sys.atoms[i].el].vis * pop * 0.98);
  U.uNA.value = NA;
  let nb = 0;
  bl.forEach((b, k) => {
    if (nb >= MAX_BONDS) return;
    U.uBondA.value[nb].set(pos[b.a][0], pos[b.a][1], pos[b.a][2], 0.055 * pop * presence[k]);
    U.uBondB.value[nb].set(pos[b.b][0], pos[b.b][1], pos[b.b][2]);
    nb++;
  });
  U.uNB.value = nb;
  const bb = bounds(pos, 2.4);
  U.uBoxMin.value.fromArray(bb.min);
  U.uBoxMax.value.fromArray(bb.max);
  orbitalBox.position.set((bb.min[0] + bb.max[0]) / 2, (bb.min[1] + bb.max[1]) / 2, (bb.min[2] + bb.max[2]) / 2);
  orbitalBox.scale.set(bb.max[0] - bb.min[0], bb.max[1] - bb.min[1], bb.max[2] - bb.min[2]);
  U.uIso.value += (state.iso - U.uIso.value) * 0.2;
  U.uOpacity.value = state.orbAlpha * state.surfAmt;
  U.uTime.value = time;

  // ---- electron measurements: resample |ψ|² once the geometry settles
  const geoKey = `${sys.id}|${state.sel}|${Math.round(s * 100)}|${Math.round(state.xi * 400)}`;
  if (state.cloudKey !== geoKey) {
    const idle = now - state.geoChangedAt > 260 && !state.playing && !holding;
    if (state.cloudKey === null || idle) {
      resampleCloud(pos, sol, target);
      state.cloudKey = geoKey;
    }
  }
  const cd = state.cloudData;
  if (cd) {
    for (let k = 0; k < cd.count; k++) {
      const p = pos[cd.near[k]] || pos[0];
      cloudPos[k * 3] = p[0] + cd.offsets[k * 3];
      cloudPos[k * 3 + 1] = p[1] + cd.offsets[k * 3 + 1];
      cloudPos[k * 3 + 2] = p[2] + cd.offsets[k * 3 + 2];
    }
    cloudGeo.attributes.position.needsUpdate = true;
    const at = (i, v) => v.set(cloudPos[i * 3], cloudPos[i * 3 + 1], cloudPos[i * 3 + 2]);
    if (cd.iMax >= 0) at(cd.iMax, tA); else tA.set(0, 99, 0);
    if (cd.iMin >= 0) at(cd.iMin, tB); else tB.set(0, 99, 0);
    lobeP.lerp(tA, lobeP.lengthSq() ? 0.15 : 1);
    lobeN.lerp(tB, lobeN.lengthSq() ? 0.15 : 1);
  }
  // measurements arrive slowly at first (watch the electrons jump), then faster to build the cloud
  const occ = sol.occ[state.sel];
  if (state.cloudAmt > 0.01) {
    state.measure.t += dt;
    const rate = state.measure.t < 4 ? 2.5 : Math.min(420, 2.5 + (state.measure.t - 4) ** 2 * 16);
    state.measure.cursor += rate * dt;
  }
  cloudUniforms.uCursor.value = state.measure.cursor;
  cloudUniforms.uOcc.value = occ > 1.5 ? 2 : 1;
  cloudUniforms.uFade.value = state.orbAlpha * state.cloudAmt * state.cloudFade;
  cloudUniforms.uTime.value = time;
  cloudUniforms.uSize.value = 30 * renderer.getPixelRatio();
  const meas = $('measure');
  if (state.cloudAmt > 0.5) {
    const count = Math.floor(state.measure.cursor);
    const eTxt = occ > 1.5 ? '<b>2 electrons</b> in this orbital' : occ > 0 ? '<b>1 electron</b> in this orbital' : '<b>empty orbital</b>: where an added electron would go';
    const html = `${eTxt} · ${count.toLocaleString()} measured position${count === 1 ? '' : 's'}`;
    if (meas.innerHTML !== html) meas.innerHTML = html;
  } else if (meas.innerHTML) meas.innerHTML = '';

  const glowAmt = state.orbAlpha * Math.max(state.surfAmt, state.cloudAmt * 0.7);
  const hasPos = cd && cd.iMax >= 0, hasNeg = cd && cd.iMin >= 0;
  lightPos.position.copy(lobeP).multiplyScalar(1.2);
  lightNeg.position.copy(lobeN).multiplyScalar(1.2);
  lightPos.intensity = hasPos ? 7 * glowAmt : 0;
  lightNeg.intensity = hasNeg ? 7 * glowAmt : 0;
  labelPos.pos.copy(lobeP);
  labelNeg.pos.copy(lobeN);
  const lobeLabels = sys.atoms.length < 14 && !state.playing;
  labelPos.visible = lobeLabels && hasPos && glowAmt > 0.3 && state.blend > 0.9;
  labelNeg.visible = lobeLabels && hasNeg && glowAmt > 0.3 && state.blend > 0.9;
  const charges = displayCharges(sys);
  for (const l of labels) {
    if (l.atom === undefined) continue;
    l.pos.fromArray(pos[l.atom]);
    if (l.charge) {
      const q = charges[l.atom];
      l.visible = pop > 0.6 && Math.abs(q) > 0.01;
      if (l.visible) {
        const t = chargeText(q);
        if (l.span.textContent !== t) { l.span.textContent = t; l.el.classList.toggle('minus', q < 0); l.el.classList.toggle('plus', q > 0); }
      }
    } else l.visible = pop > 0.6;
  }

  // ---- HUD numbers & charts
  const key = `${sys.id}|${s.toFixed(3)}|${state.xi.toFixed(3)}|${state.sel}|${sol.version}`;
  if (key !== lastUI) {
    lastUI = key;
    diagram.update(sol.energies, sol.occ);
    refreshMOText();
    const eps = fmt(sol.energies[state.sel], 1) + ' eV';
    if (sys.kind === 'reaction') {
      profile.update(state.xi);
      $('xiVal').textContent = `ξ ${state.xi.toFixed(2)}`;
      const f = sys.rx.forming[0], b = sys.rx.breaking[0];
      setChip(1, f ? `${sys.names[f[0]]}···${sys.names[f[1]]}` : 'Bond', f ? `${dist(pos[f[0]], pos[f[1]]).toFixed(2)} Å` : '—');
      setChip(2, b ? `${sys.names[b[0]]}···${sys.names[b[1]]}` : 'Progress', b ? `${dist(pos[b[0]], pos[b[1]]).toFixed(2)} Å` : `${Math.round(state.xi * 100)}%`);
      setChip(3, 'Energy', `${fmt(profileAt(sys.def, state.xi), 0)} kJ/mol`);
      setChip(4, 'Orbital ε', eps);
      const phase = state.xi < 0.33 ? 0 : state.xi < 0.67 ? 1 : 2;
      const st = $('bondState');
      st.textContent = ['reactants', sys.def.Ea > 0 ? 'transition state' : 'bond forming', 'products'][phase];
      st.className = 'state' + (phase === 1 ? ' warn' : '');
      const story = sys.def.story[phase] + (sys.def.energyNote ? ` <span class="muted">(${sys.def.energyNote})</span>` : '');
      if ($('leftFoot').innerHTML !== story) $('leftFoot').innerHTML = story;
    } else {
      morse.update(s);
      const Rrep = sys.morse.Re * s;
      $('stretchVal').textContent = `${Rrep.toFixed(2)} Å`;
      const V = sys.morse.De * (1 - Math.exp(-sys.morse.a * (Rrep - sys.morse.Re))) ** 2 - sys.morse.De;
      setChip(1, sys.bondName, `${Rrep.toFixed(2)} Å`);
      if (sys.kind === 'gallery') setChip(2, 'Bond order', sys.order);
      else {
        const hm = sol.res.homo;
        setChip(2, 'HOMO–LUMO gap', sol.res.E[hm + 1] !== undefined ? `${(sol.res.E[hm + 1] - sol.res.E[hm]).toFixed(2)} eV` : '—');
      }
      setChip(3, 'Bond energy', `${fmt(V)} eV`);
      setChip(4, 'Orbital ε', eps);
      const st = $('bondState');
      const frac = -V / sys.morse.De;
      if (s < 0.92 && V > -0.9 * sys.morse.De) { st.textContent = 'compressed'; st.className = 'state warn'; }
      else if (frac > 0.85) { st.textContent = Math.abs(s - 1) < 0.01 ? 'at Rₑ' : 'bound'; st.className = 'state'; }
      else if (frac > 0.12) { st.textContent = 'stretching'; st.className = 'state warn'; }
      else { st.textContent = 'dissociated'; st.className = 'state off'; }
    }
  }
}

function setChip(n, label, value) {
  const l = $(`c${n}l`), v = $(`c${n}v`);
  if (l.textContent !== label) l.textContent = label;
  if (v.textContent !== value) v.textContent = value;
}

// ============================================================ adaptive resolution for the raymarch
let frameAvg = 16, lastAdjust = 0;
function adaptQuality(dtMs, now) {
  frameAvg += (dtMs - frameAvg) * 0.05;
  if (now - lastAdjust < 900) return;
  lastAdjust = now;
  let sc = orbitalPass.scale;
  if (frameAvg > 21 && sc > 0.42) sc -= 0.08;
  else if (frameAvg < 14.5 && sc < 1) sc += 0.05;
  if (sc !== orbitalPass.scale) orbitalPass.setScale(Math.max(0.42, Math.min(1, sc)));
}

// ============================================================ loop
let last = performance.now(), time = 0;
function frame(now) {
  const dtMs = Math.min(100, now - last);
  last = now;
  const dt = Math.min(0.05, dtMs / 1000);
  time += dt;
  runTweens(dt);
  controls.update();
  try { update(time, dt, now); } catch (err) { console.error(err); }
  finish.uniforms.uTime.value = time;
  composer.render(dt);
  placeLabels();
  adaptQuality(dtMs, now);
  requestAnimationFrame(frame);
}

window.addEventListener('resize', () => {
  const w = window.innerWidth, h = window.innerHeight;
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
  renderer.setSize(w, h);
  composer.setSize(w, h);
  bloom.resolution.set(w, h);
  finish.uniforms.uRes.value.set(w, h);
});

document.fonts?.ready.then(() => plate.redraw());
loadGallery(0);
requestAnimationFrame(frame);
requestAnimationFrame(() => $('fade').classList.add('gone'));
