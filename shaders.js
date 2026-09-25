export const MAX_AO = 64;
export const MAX_ATOMS = 64;
export const MAX_BONDS = 80;

export const orbitalVert = /* glsl */ `
varying vec3 vWorld;
void main() {
  vec4 w = modelMatrix * vec4(position, 1.0);
  vWorld = w.xyz;
  gl_Position = projectionMatrix * viewMatrix * w;
}`;

// Raymarches the molecular orbital ψ(r) = Σ c_i φ_i(r) analytically.
// Draws up to four glassy |ψ| = iso shells, colours them by phase and adds
// volumetric glow inside the lobes. Atoms and bonds occlude the ray analytically.
export const orbitalFrag = /* glsl */ `
precision highp float;
#define MAX_AO ${MAX_AO}
#define MAX_ATOMS ${MAX_ATOMS}
#define MAX_BONDS ${MAX_BONDS}

uniform vec4 uC[MAX_AO];      // center.xyz, zeta + 100·k  (radial factor r^k)
uniform vec4 uB[MAX_AO];      // p-direction * weight, s weight
uniform int uN;
uniform vec4 uAtoms[MAX_ATOMS];
uniform int uNA;
uniform vec4 uBondA[MAX_BONDS]; // a.xyz, radius
uniform vec3 uBondB[MAX_BONDS];
uniform int uNB;
uniform vec3 uBoxMin;
uniform vec3 uBoxMax;
uniform float uIso;
uniform float uOpacity;
uniform float uGlow;
uniform float uTime;
uniform float uScan;
uniform vec3 uPos;
uniform vec3 uNeg;
uniform vec3 uLight;
varying vec3 vWorld;

float psi(vec3 p) {
  float s = 0.0;
  for (int i = 0; i < MAX_AO; i++) {
    if (i >= uN) break;
    vec3 d = p - uC[i].xyz;
    float r = length(d);
    float k = floor(uC[i].w * 0.01);
    float zr = (uC[i].w - 100.0 * k) * r;
    if (zr > 18.0) continue;
    float rad = k < 0.5 ? 1.0 : (k < 1.5 ? r : r * r);
    s += rad * (uB[i].w + dot(uB[i].xyz, d)) * exp(-zr);
  }
  return s;
}

vec3 gradDir(vec3 p) {
  const vec2 k = vec2(1.0, -1.0);
  const float h = 0.0035;
  return k.xyy * psi(p + k.xyy * h) + k.yyx * psi(p + k.yyx * h) +
         k.yxy * psi(p + k.yxy * h) + k.xxx * psi(p + k.xxx * h);
}

vec2 boxHit(vec3 ro, vec3 rd) {
  vec3 inv = 1.0 / rd;
  vec3 t0 = (uBoxMin - ro) * inv, t1 = (uBoxMax - ro) * inv;
  vec3 mn = min(t0, t1), mx = max(t0, t1);
  return vec2(max(max(mn.x, mn.y), mn.z), min(min(mx.x, mx.y), mx.z));
}

float sphHit(vec3 ro, vec3 rd, vec4 s) {
  vec3 oc = ro - s.xyz;
  float b = dot(oc, rd);
  float c = dot(oc, oc) - s.w * s.w;
  float h = b * b - c;
  if (h < 0.0) return 1e9;
  float t = -b - sqrt(h);
  return t > 0.0 ? t : 1e9;
}

float capHit(vec3 ro, vec3 rd, vec3 pa, vec3 pb, float ra) {
  vec3 ba = pb - pa, oa = ro - pa;
  float baba = dot(ba, ba), bard = dot(ba, rd), baoa = dot(ba, oa);
  float rdoa = dot(rd, oa), oaoa = dot(oa, oa);
  float a = baba - bard * bard;
  float b = baba * rdoa - baoa * bard;
  float c = baba * oaoa - baoa * baoa - ra * ra * baba;
  float h = b * b - a * c;
  if (h >= 0.0) {
    float t = (-b - sqrt(h)) / a;
    float y = baoa + t * bard;
    if (y > 0.0 && y < baba && t > 0.0) return t;
    vec3 oc = (y <= 0.0) ? oa : ro - pb;
    b = dot(rd, oc);
    c = dot(oc, oc) - ra * ra;
    h = b * b - c;
    if (h > 0.0) { float t2 = -b - sqrt(h); if (t2 > 0.0) return t2; }
  }
  return 1e9;
}

float opaqueHit(vec3 ro, vec3 rd) {
  float t = 1e9;
  for (int i = 0; i < MAX_ATOMS; i++) { if (i >= uNA) break; t = min(t, sphHit(ro, rd, uAtoms[i])); }
  for (int i = 0; i < MAX_BONDS; i++) { if (i >= uNB) break; if (uBondA[i].w > 0.004) t = min(t, capHit(ro, rd, uBondA[i].xyz, uBondB[i], uBondA[i].w)); }
  return t;
}

vec4 shade(vec3 p, vec3 rd, float v) {
  vec3 n = normalize(-sign(v) * gradDir(p));  // outward: towards decreasing |ψ|
  vec3 V = -rd;
  float nv = dot(n, V);
  bool front = nv > 0.0;
  vec3 nn = front ? n : -n;
  vec3 L = normalize(uLight);
  float fres = pow(1.0 - abs(nv), 3.0);
  float diff = max(dot(nn, L), 0.0) * 0.65 + 0.35;
  float spec = pow(max(dot(nn, normalize(L + V)), 0.0), 90.0);
  float spec2 = pow(max(dot(nn, normalize(vec3(-0.6, 0.3, -0.7) + V)), 0.0), 40.0) * 0.35;
  vec3 R = reflect(rd, nn);
  float env = smoothstep(-0.3, 1.0, R.y);
  vec3 base = v > 0.0 ? uPos : uNeg;

  // travelling contour lines, like light caught in the edge of a lens element
  float band = sin(dot(p, vec3(0.0, 1.0, 0.25)) * 26.0 - uTime * 1.6);
  float scan = smoothstep(0.93, 1.0, band) * uScan;

  vec3 col = base * (0.16 + 0.32 * diff)
           + base * fres * 1.5
           + base * scan * 1.4
           + mix(vec3(0.015), vec3(0.28, 0.32, 0.4), env) * 0.3;
  float a = front ? 0.14 + 0.7 * fres : 0.07 + 0.3 * fres;
  a = clamp(a + scan * 0.25, 0.0, 1.0);
  return vec4(col * a + vec3(1.0, 0.97, 0.92) * (spec * 1.6 + spec2) * (front ? 1.0 : 0.35), a);
}

void main() {
  vec3 ro = cameraPosition;
  vec3 rd = normalize(vWorld - ro);
  vec2 tb = boxHit(ro, rd);
  float t0 = max(tb.x, 0.0);
  float t1 = min(tb.y, opaqueHit(ro, rd));
  if (t1 <= t0 || uOpacity <= 0.001) { gl_FragColor = vec4(0.0); return; }

  float jitter = fract(sin(dot(gl_FragCoord.xy, vec2(12.9898, 78.233))) * 43758.5453);
  float t = t0 + jitter * 0.02;
  float prevF = abs(psi(ro + rd * t)) - uIso;
  vec4 acc = vec4(0.0);
  int hits = 0;

  for (int i = 0; i < 240; i++) {
    if (t >= t1 || acc.a > 0.97) break;
    // big steps far from the surface, small ones near it
    float dt = mix(0.014, 0.085, clamp(abs(prevF) / uIso, 0.0, 1.0));
    dt = min(dt, t1 - t + 1e-4);
    t += dt;
    vec3 p = ro + rd * t;
    float v = psi(p);
    float f = abs(v) - uIso;
    vec3 gc = v > 0.0 ? uPos : uNeg;

    if (f > 0.0) {
      float dens = uGlow * min(f / uIso, 2.5) * dt;
      acc.rgb += (1.0 - acc.a) * gc * dens * 0.16;
      acc.a += (1.0 - acc.a) * dens * 0.05;
    } else {
      float hq = abs(v) / uIso;
      acc.rgb += (1.0 - acc.a) * gc * hq * hq * hq * dt * 0.12 * uGlow;
    }

    if (sign(f) != sign(prevF) && hits < 4) {
      float ta = t - dt, tc = t, fa = prevF;
      for (int k = 0; k < 6; k++) {
        float tm = 0.5 * (ta + tc);
        float fm = abs(psi(ro + rd * tm)) - uIso;
        if (sign(fm) == sign(fa)) { ta = tm; fa = fm; } else { tc = tm; }
      }
      float th = 0.5 * (ta + tc);
      vec3 ph = ro + rd * th;
      vec4 c = shade(ph, rd, psi(ph));
      acc.rgb += (1.0 - acc.a) * c.rgb;
      acc.a += (1.0 - acc.a) * c.a;
      hits++;
    }
    prevF = f;
  }
  gl_FragColor = acc * uOpacity;
}`;

// Electron measurements. Each point is one sampled position from |ψ|² (a "measurement").
// Points appear one at a time in random order: the newest uOcc points are the electrons
// "right now"; older ones fade into the long-exposure cloud that |ψ|² describes.
export const cloudVert = /* glsl */ `
attribute float aSign;
attribute float aSeed;
attribute float aIdx;
uniform float uTime;
uniform float uSize;
uniform float uCursor;
uniform float uCount;
uniform float uOcc;
varying float vSign;
varying float vBright;
varying float vLive;
void main() {
  float age = uCursor - 1.0 - aIdx;
  if (uCursor > uCount) age = mod(age, uCount);
  if (age < 0.0) { gl_Position = vec4(2.0, 2.0, 2.0, 1.0); gl_PointSize = 0.0; return; }
  float live = age < uOcc ? 1.0 : 0.0;
  vec3 p = position + (1.0 - live) * 0.012 * vec3(sin(uTime * 1.3 + aSeed * 21.0), sin(uTime * 1.7 + aSeed * 37.0), sin(uTime * 1.1 + aSeed * 53.0));
  vec4 mv = modelViewMatrix * vec4(p, 1.0);
  gl_Position = projectionMatrix * mv;
  float fresh = exp(-age / 900.0);
  gl_PointSize = uSize * (live > 0.5 ? 3.4 : (0.5 + 0.35 * fract(aSeed * 7.13) + 0.5 * exp(-age / 6.0))) / -mv.z;
  vBright = live > 0.5 ? 3.2 : (0.35 + 0.65 * fresh) * (0.75 + 0.25 * sin(uTime * 2.6 + aSeed * 61.0));
  vLive = live;
  vSign = aSign;
}`;

export const cloudFrag = /* glsl */ `
uniform vec3 uPos;
uniform vec3 uNeg;
uniform float uFade;
varying float vSign;
varying float vBright;
varying float vLive;
void main() {
  float d = length(gl_PointCoord - 0.5);
  float a = smoothstep(0.5, 0.0, d);
  a *= a;
  vec3 c = vSign > 0.0 ? uPos : uNeg;
  c = mix(c, vec3(1.0), vLive * smoothstep(0.25, 0.0, d) * 0.8);
  gl_FragColor = vec4(c * a * vBright * uFade * 1.4, 1.0);
}`;

export const backdropFrag = /* glsl */ `
varying vec3 vWorld;
void main() {
  vec3 d = normalize(vWorld);
  float h = d.y;
  vec3 top = vec3(0.028, 0.032, 0.042);
  vec3 mid = vec3(0.055, 0.058, 0.07);
  vec3 bot = vec3(0.02, 0.02, 0.024);
  vec3 c = h > 0.0 ? mix(mid, top, smoothstep(0.0, 0.7, h)) : mix(mid, bot, smoothstep(0.0, 0.25, -h));
  // warm glow behind the stage
  c += vec3(0.06, 0.035, 0.012) * pow(max(0.0, 1.0 - abs(h - 0.05) * 3.0), 3.0) * smoothstep(-0.2, 0.8, -d.z);
  gl_FragColor = vec4(c, 1.0);
}`;

export const compositeFrag = /* glsl */ `
uniform sampler2D tBase;
uniform sampler2D tOrb;
varying vec2 vUv;
void main() {
  vec4 b = texture2D(tBase, vUv);
  vec4 o = texture2D(tOrb, vUv);
  gl_FragColor = vec4(b.rgb * (1.0 - o.a) + o.rgb, 1.0);
}`;

export const finishFrag = /* glsl */ `
uniform sampler2D tDiffuse;
uniform float uTime;
uniform vec2 uRes;
varying vec2 vUv;
void main() {
  vec4 c = texture2D(tDiffuse, vUv);
  vec2 q = vUv - 0.5;
  q.x *= uRes.x / uRes.y;
  float v = smoothstep(1.05, 0.25, length(q));
  c.rgb *= mix(0.5, 1.0, v);
  float n = fract(sin(dot(vUv * uRes + fract(uTime) * 91.0, vec2(12.9898, 78.233))) * 43758.5453);
  c.rgb += (n - 0.5) * 0.022;
  gl_FragColor = c;
}`;
