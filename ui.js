const NS = 'http://www.w3.org/2000/svg';

function svgEl(tag, attrs = {}, parent) {
  const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(e);
  return e;
}

// "σ*_{u}1s" → σ*<sub>u</sub>1s
export function richHTML(s) {
  return s.replace(/_\{([^}]*)\}/g, '<sub>$1</sub>').replace(/\^\{([^}]*)\}/g, '<sup>$1</sup>');
}

function richSVG(text, s) {
  const re = /([_^])\{([^}]*)\}/g;
  let last = 0, m;
  const push = (str, shift) => {
    if (!str) return;
    const t = svgEl('tspan', shift ? { 'baseline-shift': shift, 'font-size': '72%' } : {}, text);
    t.textContent = str;
  };
  while ((m = re.exec(s))) {
    push(s.slice(last, m.index));
    push(m[2], m[1] === '_' ? 'sub' : 'super');
    last = re.lastIndex;
  }
  push(s.slice(last));
}

// Continuous, monotonic energy → y map. Knots are levels spaced by √ΔE,
// so deep levels don't squash the valence region.
function makeEnergyScale(energies, top, bottom) {
  const es = [...new Set(energies.map((e) => Math.round(e * 10) / 10))].sort((a, b) => a - b);
  const pos = [0];
  for (let i = 1; i < es.length; i++) pos.push(pos[i - 1] + Math.sqrt(es[i] - es[i - 1]) + 0.55);
  const span = pos[pos.length - 1] || 1;
  const toY = (p) => bottom - (p / span) * (bottom - top);
  return (e) => {
    if (es.length === 1) return (top + bottom) / 2;
    let i = 0;
    if (e <= es[0]) i = 0;
    else if (e >= es[es.length - 1]) i = es.length - 2;
    else while (e > es[i + 1]) i++;
    const f = (e - es[i]) / (es[i + 1] - es[i]);
    const y = toY(pos[i] + f * (pos[i + 1] - pos[i]));
    return Math.max(top - 14, Math.min(bottom + 12, y));
  };
}

function arrows(x, y, n) {
  if (n <= 0) return '';
  const up = (cx) => `M${cx} ${y + 6}L${cx} ${y - 6}M${cx - 2.6} ${y - 3}L${cx} ${y - 6.2}L${cx + 2.6} ${y - 3}`;
  const dn = (cx) => `M${cx} ${y - 6}L${cx} ${y + 6}M${cx - 2.6} ${y + 3}L${cx} ${y + 6.2}L${cx + 2.6} ${y + 3}`;
  return n < 1.5 ? up(x) : up(x - 4) + dn(x + 4);
}

// view: { mos: [{ name, occ, g, tag, from }], atomic?, capL?, capR?, capM? }
export class MODiagram {
  constructor(svg, onSelect) {
    this.svg = svg;
    this.onSelect = onSelect;
    this.W = 380;
    this.H = 250;
    svg.setAttribute('viewBox', `0 0 ${this.W} ${this.H}`);
  }

  setMolecule(view, scaleEnergies) {
    const { W, H } = this;
    const svg = this.svg;
    svg.replaceChildren();
    this.view = view;
    const diatomic = !!view.atomic;
    this.y = makeEnergyScale([...scaleEnergies, ...(view.atomic || []).map((a) => a.e)], 22, H - 18);

    const ax = svgEl('g', { class: 'axis' }, svg);
    svgEl('path', { d: `M10 ${H - 10}L10 12M6.5 17L10 11L13.5 17` }, ax);
    svgEl('text', { x: 16, y: 16, class: 'ax-label' }, ax).textContent = 'E';

    this.items = [];
    this.connectors = [];
    this.atomLevels = [];
    const cx = diatomic ? W / 2 : 150;
    const groups = new Map();
    view.mos.forEach((m, i) => {
      if (!groups.has(m.g)) groups.set(m.g, []);
      groups.get(m.g).push(i);
    });
    const connLayer = svgEl('g', { class: 'conn' }, svg);

    if (diatomic) {
      for (const side of [-1, 1]) {
        const colX = side < 0 ? 46 : W - 46;
        for (const lv of view.atomic) {
          const w = lv.n > 1 ? 13 : 30, gap = 4;
          const total = lv.n * w + (lv.n - 1) * gap;
          const x0 = colX - total / 2;
          const y = this.y(lv.e);
          const g = svgEl('g', { class: 'alevel' }, svg);
          let left = lv.el;
          for (let k = 0; k < lv.n; k++) {
            const xs = x0 + k * (w + gap);
            svgEl('line', { x1: xs, x2: xs + w, y1: y, y2: y }, g);
            // Hund's rule: singly fill degenerate levels before pairing
            const nE = lv.n === 1 ? left : (lv.el <= lv.n ? (k < lv.el ? 1 : 0) : (k < lv.el - lv.n ? 2 : 1));
            if (lv.n === 1) left = 0;
            svgEl('path', { d: arrows(xs + w / 2, y, nE), class: 'e' }, g);
          }
          svgEl('text', { x: side < 0 ? x0 - 6 : x0 + total + 6, y: y + 3.5, 'text-anchor': side < 0 ? 'end' : 'start', class: 'alabel' }, g).textContent = lv.name;
          this.atomLevels.push({ lv, side, xIn: side < 0 ? x0 + total : x0, y });
        }
        svgEl('text', { x: colX, y: H - 4, 'text-anchor': 'middle', class: 'col-cap' }, svg).textContent = side < 0 ? view.capL : view.capR;
      }
      svgEl('text', { x: cx, y: H - 4, 'text-anchor': 'middle', class: 'col-cap' }, svg).textContent = view.capM;
    }

    for (const [, members] of groups) {
      const n = members.length;
      const w = n > 1 ? (n > 2 ? 22 : 28) : 44, gap = 9;
      const total = n * w + (n - 1) * gap;
      const x0 = cx - total / 2;
      const first = view.mos[members[0]];
      const gEl = svgEl('g', { class: 'group' }, svg);
      let label = null, tag = null;
      if (!diatomic) {
        label = svgEl('text', { x: x0 + total + 10, class: 'mlabel' }, gEl);
        richSVG(label, first.name.replace(/ \([xyz]\)/, ''));
      }
      const tagged = members.map((i) => view.mos[i].tag).find(Boolean);
      if (tagged) {
        tag = svgEl('text', { x: diatomic ? x0 + total + 8 : x0 - 8, 'text-anchor': diatomic ? 'start' : 'end', class: 'tag ' + tagged.toLowerCase() }, gEl);
        tag.textContent = tagged;
      }
      members.forEach((i, k) => {
        const xs = x0 + k * (w + gap);
        const it = svgEl('g', { class: 'lvl' }, gEl);
        const hit = svgEl('rect', { x: xs - 3, width: w + 6, height: 22, class: 'hit' }, it);
        const line = svgEl('line', { x1: xs, x2: xs + w }, it);
        const e = svgEl('path', { class: 'e' }, it);
        svgEl('title', {}, it).textContent = view.mos[i].name.replace(/[_^]\{([^}]*)\}/g, '$1');
        it.addEventListener('click', () => this.onSelect(i));
        this.items[i] = { it, hit, line, e, xs, w, label: k === 0 ? label : null, tag: k === 0 ? tag : null };
      });
      if (diatomic) {
        for (const al of this.atomLevels) {
          if (!first.from?.includes(al.lv.name)) continue;
          const path = svgEl('path', { class: 'connector' }, connLayer);
          this.connectors.push({ path, al, g: members[0], xOut: al.side < 0 ? x0 : x0 + total });
        }
      }
    }
  }

  update(energies, occs) {
    this.view.mos.forEach((m, i) => {
      const it = this.items[i];
      const y = this.y(energies[i]);
      const occ = occs ? occs[i] : m.occ;
      it.line.setAttribute('y1', y);
      it.line.setAttribute('y2', y);
      it.hit.setAttribute('y', y - 11);
      it.e.setAttribute('d', arrows(it.xs + it.w / 2, y, occ));
      it.it.classList.toggle('occ', occ > 0);
      if (it.label) it.label.setAttribute('y', y + 3.5);
      if (it.tag) it.tag.setAttribute('y', y + 3.5);
    });
    // hide group labels that would collide with the one below
    const labs = this.items.filter((it) => it && it.label).map((it) => ({ el: it.label, y: +it.label.getAttribute('y') })).sort((a, b) => b.y - a.y);
    let lastY = Infinity;
    for (const l of labs) {
      const show = lastY - l.y >= 12.5;
      l.el.style.opacity = show ? 1 : 0;
      if (show) lastY = l.y;
    }
    for (const c of this.connectors) {
      const y = this.y(energies[c.g]);
      c.path.setAttribute('d', `M${c.al.xIn + c.al.side * -3} ${c.al.y}L${c.xOut + c.al.side * 3} ${y}`);
    }
  }

  setSelected(i) {
    this.items.forEach((it, k) => it && it.it.classList.toggle('sel', k === i));
  }
}

// Draggable SVG chart base: maps pointer x to a value and reports drag start/end.
class DragChart {
  constructor(svg, W, H, handlers) {
    this.svg = svg;
    this.W = W;
    this.H = H;
    svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
    let dragging = false;
    const handle = (ev) => {
      const r = svg.getBoundingClientRect();
      handlers.onDrag(this.fromPx(((ev.clientX - r.left) / r.width) * W));
    };
    svg.addEventListener('pointerdown', (ev) => { dragging = true; svg.setPointerCapture(ev.pointerId); handlers.onStart?.(); handle(ev); });
    svg.addEventListener('pointermove', (ev) => dragging && handle(ev));
    const end = () => { if (dragging) { dragging = false; handlers.onEnd?.(); } };
    svg.addEventListener('pointerup', end);
    svg.addEventListener('pointercancel', end);
  }

  frame() {
    this.svg.replaceChildren();
    const defs = svgEl('defs', {}, this.svg);
    const grad = svgEl('linearGradient', { id: 'mfill' + this.constructor.name, x1: 0, y1: 0, x2: 0, y2: 1 }, defs);
    svgEl('stop', { offset: '0', 'stop-color': '#4fd8ff', 'stop-opacity': '0' }, grad);
    svgEl('stop', { offset: '1', 'stop-color': '#4fd8ff', 'stop-opacity': '0.22' }, grad);
    return `url(#mfill${this.constructor.name})`;
  }

  marker() {
    this.guide = svgEl('line', { class: 'guide' }, this.svg);
    this.halo = svgEl('circle', { r: 9, class: 'halo' }, this.svg);
    this.dot = svgEl('circle', { r: 4.2, class: 'dot' }, this.svg);
  }

  place(x, y, yBase) {
    for (const [el, a, b] of [[this.dot, 'cx', 'cy'], [this.halo, 'cx', 'cy']]) { el.setAttribute(a, x); el.setAttribute(b, y); }
    this.guide.setAttribute('x1', x);
    this.guide.setAttribute('x2', x);
    this.guide.setAttribute('y1', y);
    this.guide.setAttribute('y2', yBase);
  }
}

export class MorseChart extends DragChart {
  constructor(svg, handlers) {
    super(svg, 320, 150, handlers);
    this.m = { l: 34, r: 12, t: 12, b: 24 };
  }

  fromPx(px) {
    const R = this.R0 + ((px - this.m.l) / (this.W - this.m.l - this.m.r)) * (this.R1 - this.R0);
    return R / this.morse.Re;
  }

  V(R) {
    const { De, Re, a } = this.morse;
    return De * (1 - Math.exp(-a * (R - Re))) ** 2 - De;
  }

  setMorse(morse, sMin, sMax) {
    this.morse = morse;
    const { W, H, m } = this;
    const { De, Re } = morse;
    this.R0 = Re * sMin;
    this.R1 = Re * sMax;
    this.E0 = -De * 1.12;
    this.E1 = De * 0.62;
    this.x = (R) => m.l + ((R - this.R0) / (this.R1 - this.R0)) * (W - m.l - m.r);
    this.yE = (E) => m.t + (1 - (Math.min(E, this.E1) - this.E0) / (this.E1 - this.E0)) * (H - m.t - m.b);
    const fill = this.frame();
    const svg = this.svg;

    const g = svgEl('g', { class: 'grid' }, svg);
    const span = this.R1 - this.R0;
    const step = span < 1.5 ? 0.25 : span < 3 ? 0.5 : 1;
    for (let R = Math.ceil(this.R0 / step) * step; R <= this.R1 + 1e-6; R += step) {
      svgEl('line', { x1: this.x(R), x2: this.x(R), y1: m.t, y2: H - m.b }, g);
      svgEl('text', { x: this.x(R), y: H - m.b + 13, 'text-anchor': 'middle', class: 'tick' }, svg).textContent = +R.toFixed(2);
    }
    const eStep = De > 6 ? 4 : De > 3 ? 2 : 1;
    for (let E = -Math.floor(De / eStep) * eStep; E <= this.E1; E += eStep) {
      svgEl('line', { x1: m.l, x2: W - m.r, y1: this.yE(E), y2: this.yE(E) }, g);
      svgEl('text', { x: m.l - 5, y: this.yE(E) + 3, 'text-anchor': 'end', class: 'tick' }, svg).textContent = E;
    }
    svgEl('text', { x: W - m.r, y: m.t - 3, 'text-anchor': 'end', class: 'ax-label' }, svg).textContent = 'R (Å) →';
    svgEl('text', { x: 4, y: m.t - 3, class: 'ax-label' }, svg).textContent = 'eV';
    svgEl('line', { x1: m.l, x2: W - m.r, y1: this.yE(0), y2: this.yE(0), class: 'asym' }, svg);
    svgEl('text', { x: W - m.r - 2, y: this.yE(0) - 4, 'text-anchor': 'end', class: 'note' }, svg).textContent = 'separated atoms';

    let d = '';
    for (let i = 0; i <= 160; i++) {
      const R = this.R0 + (i / 160) * (this.R1 - this.R0);
      d += `${i ? 'L' : 'M'}${this.x(R).toFixed(2)} ${this.yE(this.V(R)).toFixed(2)}`;
    }
    svgEl('path', { d: `${d}L${this.x(this.R1)} ${H - m.b}L${this.x(this.R0)} ${H - m.b}Z`, fill }, svg);
    svgEl('path', { d, class: 'curve' }, svg);
    svgEl('line', { x1: this.x(Re), x2: this.x(Re), y1: this.yE(-De), y2: H - m.b, class: 're' }, svg);
    richSVG(svgEl('text', { x: this.x(Re) + 5, y: this.yE(-De) - 6, class: 'note' }, svg), 'R_{e}');
    const bx = this.x(this.R1) - 26;
    svgEl('path', { d: `M${bx} ${this.yE(0)}L${bx} ${this.yE(-De)}`, class: 'de' }, svg);
    richSVG(svgEl('text', { x: bx - 4, y: this.yE(-De / 2) + 3, 'text-anchor': 'end', class: 'note' }, svg), `D_{e} ${De.toFixed(2)} eV`);
    this.marker();
  }

  update(s) {
    const R = this.morse.Re * s;
    this.place(this.x(R), this.yE(this.V(R)), this.H - this.m.b);
  }
}

// Reaction coordinate diagram: schematic enthalpy profile through R, TS, P
export class ProfileChart extends DragChart {
  constructor(svg, handlers) {
    super(svg, 320, 150, handlers);
    this.m = { l: 40, r: 12, t: 16, b: 24 };
  }

  fromPx(px) {
    return Math.max(0, Math.min(1, (px - this.m.l) / (this.W - this.m.l - this.m.r)));
  }

  setProfile(def, fn) {
    this.fn = fn;
    const { W, H, m } = this;
    const lo = Math.min(0, def.dH) * 1.12, hi = Math.max(def.Ea, 0, 40) * 1.25 + 10;
    this.x = (xi) => m.l + xi * (W - m.l - m.r);
    this.yE = (E) => m.t + (1 - (E - lo) / (hi - lo)) * (H - m.t - m.b);
    const fill = this.frame();
    const svg = this.svg;
    const g = svgEl('g', { class: 'grid' }, svg);
    const step = hi - lo > 500 ? 200 : hi - lo > 200 ? 100 : 50;
    for (let E = Math.ceil(lo / step) * step; E <= hi; E += step) {
      svgEl('line', { x1: m.l, x2: W - m.r, y1: this.yE(E), y2: this.yE(E) }, g);
      svgEl('text', { x: m.l - 5, y: this.yE(E) + 3, 'text-anchor': 'end', class: 'tick' }, svg).textContent = E;
    }
    svgEl('text', { x: 4, y: m.t - 5, class: 'ax-label' }, svg).textContent = 'kJ/mol';
    [['reactants', 0, 'start'], ['TS', 0.5, 'middle'], ['products', 1, 'end']].forEach(([t, xi, anchor]) => {
      svgEl('text', { x: this.x(xi), y: H - 6, 'text-anchor': anchor, class: 'tick' }, svg).textContent = def.Ea > 0 || t !== 'TS' ? t : '';
    });
    let d = '';
    for (let i = 0; i <= 160; i++) d += `${i ? 'L' : 'M'}${this.x(i / 160).toFixed(2)} ${this.yE(fn(i / 160)).toFixed(2)}`;
    svgEl('path', { d: `${d}L${this.x(1)} ${H - m.b}L${this.x(0)} ${H - m.b}Z`, fill }, svg);
    svgEl('path', { d, class: 'curve' }, svg);
    svgEl('line', { x1: m.l, x2: W - m.r, y1: this.yE(0), y2: this.yE(0), class: 'asym' }, svg);
    if (def.Ea > 0) {
      const x = this.x(0.5);
      svgEl('path', { d: `M${x - 30} ${this.yE(0)}L${x - 30} ${this.yE(def.Ea)}`, class: 'de' }, svg);
      richSVG(svgEl('text', { x: x - 34, y: this.yE(def.Ea / 2) + 3, 'text-anchor': 'end', class: 'note' }, svg), `E_{a} ${def.Ea}`);
    }
    const xp = this.x(0.93);
    svgEl('path', { d: `M${xp} ${this.yE(0)}L${xp} ${this.yE(def.dH)}`, class: 'de' }, svg);
    richSVG(svgEl('text', { x: xp - 4, y: this.yE(def.dH / 2) + 3, 'text-anchor': 'end', class: 'note' }, svg), `ΔH ${def.dH}`);
    this.marker();
  }

  update(xi) {
    this.place(this.x(xi), this.yE(this.fn(xi)), this.H - this.m.b);
  }
}
