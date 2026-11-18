// Globe renderer (from rtl-sdr-workbench display/globe.js, unchanged but for asset paths): three.js 0.160 (vendor/three, MIT), a NASA Black Marble night
// texture, Natural Earth coastlines (world-atlas, vendor/earth) and a 30° graticule. Markers keep a
// constant size on screen; tracks follow great circles; aircraft-to-ground-station links are arcs.
// The camera orbits the Earth with north kept up: drag rotates, the wheel zooms, and every action
// also has an on-screen button in the page.
import * as THREE from './vendor/three/three.module.min.js';

const R = 1;
const MIN_DIST = 1.08, MAX_DIST = 6;
const rad = d => d * Math.PI / 180;

// The same mapping as SphereGeometry's UVs with an equirectangular texture (u = 0 at 180° W).
export function toVec(lat, lon, r = R) {
 const phi = rad(lon + 180), theta = rad(90 - lat);
 return new THREE.Vector3(-r * Math.sin(theta) * Math.cos(phi), r * Math.cos(theta), r * Math.sin(theta) * Math.sin(phi));
}
function toLatLon(v) {
 const n = v.clone().normalize();
 const lat = 90 - Math.acos(Math.max(-1, Math.min(1, n.y))) * 180 / Math.PI;
 let lon = Math.atan2(n.z, -n.x) * 180 / Math.PI - 180;
 if (lon < -180) lon += 360;
 return {lat, lon};
}
// Points along the great circle from a to b (unit vectors), at most `step` radians apart.
function greatCircle(a, b, step = rad(1)) {
 const ang = a.angleTo(b), n = Math.max(1, Math.ceil(ang / step)), out = [];
 if (ang < 1e-9) return [a.clone()];
 const s = Math.sin(ang);
 for (let i = 0; i <= n; i++) {
  const f = i / n, w1 = Math.sin((1 - f) * ang) / s, w2 = Math.sin(f * ang) / s;
  out.push(a.clone().multiplyScalar(w1).add(b.clone().multiplyScalar(w2)));
 }
 return out;
}

// Round, square or diamond points of a fixed pixel size; `shape` 0 circle, 1 square, 2 diamond, 3 ring.
const POINT_VS = `attribute vec3 color; attribute float size; attribute float shape; varying vec3 vColor; varying float vShape;
void main() { vColor = color; vShape = shape; vec4 mv = modelViewMatrix * vec4(position, 1.0); gl_Position = projectionMatrix * mv; gl_PointSize = size * PIXEL_RATIO; }`;
const POINT_FS = `varying vec3 vColor; varying float vShape;
void main() {
 vec2 p = gl_PointCoord * 2.0 - 1.0; float a = 1.0;
 if (vShape < 0.5) { float d = length(p); if (d > 1.0) discard; a = smoothstep(1.0, 0.75, d); if (d > 0.72) { gl_FragColor = vec4(0.03, 0.06, 0.1, a); return; } }
 else if (vShape < 1.5) { if (max(abs(p.x), abs(p.y)) > 0.78) { gl_FragColor = vec4(0.03, 0.06, 0.1, 1.0); return; } }
 else if (vShape < 2.5) { float d = abs(p.x) + abs(p.y); if (d > 1.0) discard; if (d > 0.74) { gl_FragColor = vec4(0.03, 0.06, 0.1, 1.0); return; } }
 else { float d = length(p); if (d > 1.0 || d < 0.72) discard; a = smoothstep(1.0, 0.9, d); }
 gl_FragColor = vec4(vColor, a);
}`;

export class Globe {
 // container: element to fill. onPick(item | null) when the person clicks a marker or empty space.
 constructor(container, {onPick = () => {}, onView = () => {}} = {}) {
  this.container = container; this.onPick = onPick; this.onView = onView;
  this.items = [];            // [{id, lat, lon, color, size, shape, data}]
  this.view = {lat: 30, lon: -40, dist: 3.2};
  this.spin = false; this.follow = null; this.anim = null;
  this.renderer = new THREE.WebGLRenderer({antialias: true, alpha: false});
  this.renderer.setPixelRatio(Math.min(2, devicePixelRatio || 1));
  this.renderer.setClearColor(0x04080e);
  container.append(this.renderer.domElement);
  this.scene = new THREE.Scene();
  this.camera = new THREE.PerspectiveCamera(40, 1, 0.01, 50);
  this.earth = new THREE.Group(); this.scene.add(this.earth);
  this.layers = {};
  this.buildEarth();
  this.points = this.makePoints(); this.scene.add(this.points);
  this.highlight = this.makePoints(); this.scene.add(this.highlight);
  this.tracks = new THREE.LineSegments(new THREE.BufferGeometry(), new THREE.LineBasicMaterial({vertexColors: true, transparent: true, opacity: 0.85}));
  this.links = new THREE.LineSegments(new THREE.BufferGeometry(), new THREE.LineBasicMaterial({color: 0xffa94d, transparent: true, opacity: 0.45}));
  this.selLink = new THREE.LineSegments(new THREE.BufferGeometry(), new THREE.LineBasicMaterial({color: 0xffd08a, transparent: true, opacity: 0.95}));
  this.scene.add(this.tracks, this.links, this.selLink);
  this.bindInput();
  new ResizeObserver(() => this.resize()).observe(container);
  this.resize();
  const loop = () => { this.frame(); this.raf = requestAnimationFrame(loop); };
  this.raf = requestAnimationFrame(loop);
 }

 buildEarth() {
  const base = new URL('./vendor/earth/', import.meta.url).href;
  const tex = new THREE.TextureLoader().load(base + 'earth-night.jpg');
  tex.colorSpace = THREE.SRGBColorSpace; tex.anisotropy = 8;
  this.surface = new THREE.Mesh(new THREE.SphereGeometry(R, 96, 64), new THREE.MeshBasicMaterial({map: tex, color: 0xb8c8dc}));
  this.earth.add(this.surface);
  // Graticule every 30°.
  const g = [];
  for (let lat = -60; lat <= 60; lat += 30) for (let lon = -180; lon < 180; lon += 2) g.push(toVec(lat, lon, 1.0008), toVec(lat, lon + 2, 1.0008));
  for (let lon = -180; lon < 180; lon += 30) for (let lat = -88; lat < 88; lat += 2) g.push(toVec(lat, lon, 1.0008), toVec(lat + 2, lon, 1.0008));
  this.layers.grid = new THREE.LineSegments(new THREE.BufferGeometry().setFromPoints(g), new THREE.LineBasicMaterial({color: 0x5fb4ff, transparent: true, opacity: 0.14}));
  this.earth.add(this.layers.grid);
  // Atmosphere rim.
  const rim = new THREE.Mesh(new THREE.SphereGeometry(R * 1.035, 64, 48), new THREE.ShaderMaterial({
   transparent: true, side: THREE.BackSide, depthWrite: false,
   vertexShader: 'varying vec3 n; void main(){ n = normalize(normalMatrix * normal); gl_Position = projectionMatrix * modelViewMatrix * vec4(position,1.0); }',
   fragmentShader: 'varying vec3 n; void main(){ float i = pow(0.62 - dot(n, vec3(0.0,0.0,1.0)), 3.0); gl_FragColor = vec4(0.37,0.7,1.0,1.0) * clamp(i, 0.0, 1.0); }'}));
  this.scene.add(rim);
  this.loadCoastlines(base + 'land-50m.json');
 }
 // world-atlas TopoJSON (quantized arcs) → line segments; only the arcs are needed for outlines.
 async loadCoastlines(url) {
  try {
   const topo = await (await fetch(url)).json(), [sx, sy] = topo.transform.scale, [tx, ty] = topo.transform.translate, pts = [];
   for (const arc of topo.arcs) {
    let x = 0, y = 0, prev = null;
    for (const [dx, dy] of arc) {
     x += dx; y += dy;
     const lon = x * sx + tx, lat = y * sy + ty, v = toVec(lat, lon, 1.0012);
     if (prev) pts.push(prev, v);
     prev = v;
    }
   }
   this.layers.coast = new THREE.LineSegments(new THREE.BufferGeometry().setFromPoints(pts), new THREE.LineBasicMaterial({color: 0x9fd4ff, transparent: true, opacity: 0.55}));
   this.layers.coast.visible = this.coastVisible ?? true;
   this.earth.add(this.layers.coast);
  } catch (e) { console.warn('Coastlines failed to load', e); }
 }
 setLayer(name, on) {
  if (name === 'coast') { this.coastVisible = on; if (this.layers.coast) this.layers.coast.visible = on; }
  else if (name === 'grid') this.layers.grid.visible = on;
  else if (name === 'lights') this.surface.material.color.set(on ? 0xb8c8dc : 0x1c2a3a);
 }

 makePoints() {
  const geo = new THREE.BufferGeometry();
  const mat = new THREE.ShaderMaterial({vertexShader: POINT_VS, fragmentShader: POINT_FS, transparent: true,
   defines: {PIXEL_RATIO: this.renderer.getPixelRatio().toFixed(2)}});
  const p = new THREE.Points(geo, mat); p.frustumCulled = false;
  return p;
 }
 fillPoints(points, list) {
  const n = list.length, pos = new Float32Array(n * 3), col = new Float32Array(n * 3), size = new Float32Array(n), shape = new Float32Array(n), c = new THREE.Color();
  list.forEach((it, i) => {
   const v = toVec(it.lat, it.lon, 1.004);
   pos.set([v.x, v.y, v.z], i * 3); c.set(it.color); col.set([c.r, c.g, c.b], i * 3);
   size[i] = it.size; shape[i] = it.shape;
  });
  const g = points.geometry;
  g.setAttribute('position', new THREE.BufferAttribute(pos, 3)); g.setAttribute('color', new THREE.BufferAttribute(col, 3));
  g.setAttribute('size', new THREE.BufferAttribute(size, 1)); g.setAttribute('shape', new THREE.BufferAttribute(shape, 1));
  g.computeBoundingSphere();
 }
 // items: [{id, lat, lon, color, size, shape, data}]
 setItems(items) { this.items = items; this.fillPoints(this.points, items); this.refreshHighlight(); }
 select(id) { this.selected = id; this.refreshHighlight(); }
 refreshHighlight() {
  const it = this.items.find(x => x.id === this.selected);
  this.fillPoints(this.highlight, it ? [{...it, color: '#ffffff', size: it.size + 12, shape: 3}] : []);
 }
 // tracks: [{color, points: [{lat, lon}]}]; consecutive points joined along great circles.
 setTracks(tracks) {
  const pos = [], col = [], c = new THREE.Color();
  for (const tr of tracks) {
   c.set(tr.color);
   for (let i = 1; i < tr.points.length; i++) {
    const seg = greatCircle(toVec(tr.points[i - 1].lat, tr.points[i - 1].lon), toVec(tr.points[i].lat, tr.points[i].lon));
    for (let j = 1; j < seg.length; j++) {
     const a = seg[j - 1].multiplyScalar(1.003), b = seg[j].clone().multiplyScalar(1.003);
     pos.push(a.x, a.y, a.z, b.x, b.y, b.z); col.push(c.r, c.g, c.b, c.r, c.g, c.b);
    }
   }
  }
  this.tracks.geometry.dispose();
  this.tracks.geometry = new THREE.BufferGeometry();
  this.tracks.geometry.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
  this.tracks.geometry.setAttribute('color', new THREE.Float32BufferAttribute(col, 3));
 }
 // links: [{from: {lat, lon}, to: {lat, lon}}] drawn as arcs that rise with distance.
 arcs(links) {
  const pos = [];
  for (const l of links) {
   const a = toVec(l.from.lat, l.from.lon), b = toVec(l.to.lat, l.to.lon), seg = greatCircle(a, b, rad(1.5)), h = 0.06 + a.angleTo(b) * 0.08;
   const lifted = seg.map((v, i) => v.multiplyScalar(1.003 + h * Math.sin(Math.PI * i / (seg.length - 1 || 1))));
   for (let i = 1; i < lifted.length; i++) pos.push(lifted[i - 1].x, lifted[i - 1].y, lifted[i - 1].z, lifted[i].x, lifted[i].y, lifted[i].z);
  }
  const g = new THREE.BufferGeometry(); g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
  return g;
 }
 setLinks(links) { this.links.geometry.dispose(); this.links.geometry = this.arcs(links); }
 setSelectedLink(link) { this.selLink.geometry.dispose(); this.selLink.geometry = this.arcs(link ? [link] : []); }

 // ---- View ----
 zoom(f) { this.view.dist = Math.max(MIN_DIST, Math.min(MAX_DIST, 1 + (this.view.dist - 1) * f)); this.follow = null; this.anim = null; }
 reset() { this.flyTo(30, -40, 3.2); }
 setSpin(on) { this.spin = on; }
 setFollow(id) { this.follow = id; }
 // Animates the view to lat/lon (and distance) over about a second.
 flyTo(lat, lon, dist = this.view.dist) {
  const from = {...this.view};
  let dl = ((lon - from.lon + 540) % 360) - 180;
  this.anim = {from, to: {lat: Math.max(-85, Math.min(85, lat)), lon: from.lon + dl, dist}, t0: performance.now()};
 }
 // Centres on a set of positions with a distance that keeps them in view.
 fit(positions) {
  if (!positions.length) return;
  const sum = new THREE.Vector3();
  for (const p of positions) sum.add(toVec(p.lat, p.lon));
  if (sum.length() < 1e-6) { this.reset(); return; }
  const c = toLatLon(sum), cv = toVec(c.lat, c.lon);
  const spread = Math.max(...positions.map(p => cv.angleTo(toVec(p.lat, p.lon))));
  this.flyTo(c.lat, c.lon, Math.max(1.25, Math.min(3.6, 1.15 + spread * 2.4)));
 }
 resize() {
  const w = this.container.clientWidth || 1, h = this.container.clientHeight || 1;
  this.renderer.setSize(w, h, false); this.camera.aspect = w / h; this.camera.updateProjectionMatrix();
 }
 frame() {
  const now = performance.now();
  if (this.anim) {
   const f = Math.min(1, (now - this.anim.t0) / 900), e = f < 0.5 ? 2 * f * f : 1 - (-2 * f + 2) ** 2 / 2, {from, to} = this.anim;
   this.view = {lat: from.lat + (to.lat - from.lat) * e, lon: from.lon + (to.lon - from.lon) * e, dist: from.dist + (to.dist - from.dist) * e};
   if (f >= 1) { this.anim = null; this.view.lon = ((this.view.lon + 540) % 360) - 180; }
  } else if (this.follow) {
   const it = this.items.find(x => x.id === this.follow);
   if (it) { this.view.lat += (it.lat - this.view.lat) * 0.08; this.view.lon += ((((it.lon - this.view.lon) + 540) % 360) - 180) * 0.08; }
  } else if (this.spin && !this.dragging) this.view.lon = ((this.view.lon + 0.04 + 540) % 360) - 180;
  const v = toVec(this.view.lat, this.view.lon, this.view.dist);
  this.camera.position.copy(v); this.camera.up.set(0, 1, 0); this.camera.lookAt(0, 0, 0);
  this.renderer.render(this.scene, this.camera);
  this.onView(this.view);
 }

 // ---- Input ----
 bindInput() {
  const el = this.renderer.domElement;
  let down = null, moved = false;
  el.addEventListener('pointerdown', e => { down = {x: e.clientX, y: e.clientY, lat: this.view.lat, lon: this.view.lon}; moved = false; el.setPointerCapture(e.pointerId); });
  el.addEventListener('pointermove', e => {
   if (!down) { el.style.cursor = this.pickAt(e.clientX, e.clientY) ? 'pointer' : 'grab'; return; }
   const dx = e.clientX - down.x, dy = e.clientY - down.y;
   if (!moved && Math.hypot(dx, dy) < 4) return;
   moved = true; this.dragging = true; this.anim = null; this.follow = null;
   const k = (this.view.dist - 1) * 0.11;   // degrees per pixel: slower when zoomed in
   this.view.lon = ((down.lon - dx * k + 540) % 360) - 180;
   this.view.lat = Math.max(-85, Math.min(85, down.lat + dy * k));
   el.style.cursor = 'grabbing';
  });
  el.addEventListener('pointerup', e => {
   if (down && !moved) this.onPick(this.pickAt(e.clientX, e.clientY));
   down = null; this.dragging = false; el.style.cursor = 'grab';
  });
  el.addEventListener('wheel', e => { e.preventDefault(); this.zoom(e.deltaY > 0 ? 1.15 : 1 / 1.15); }, {passive: false});
 }
 // The visible marker nearest the pointer, within 14 px.
 pickAt(x, y) {
  const r = this.renderer.domElement.getBoundingClientRect(), cam = this.camera.position, v = new THREE.Vector3();
  let best = null, bestD = 14 * 14;
  for (const it of this.items) {
   const p = toVec(it.lat, it.lon, 1.004);
   if (p.dot(cam) <= p.lengthSq()) continue;    // behind the Earth
   v.copy(p).project(this.camera);
   const sx = (v.x + 1) / 2 * r.width + r.left, sy = (1 - v.y) / 2 * r.height + r.top, d = (sx - x) ** 2 + (sy - y) ** 2;
   if (d < bestD) { bestD = d; best = it; }
  }
  return best;
 }
 // Screen position of a lat/lon, or null when it is behind the Earth.
 screenOf(lat, lon) {
  const p = toVec(lat, lon, 1.004), cam = this.camera.position;
  if (p.dot(cam) <= p.lengthSq()) return null;
  const v = p.project(this.camera), r = this.renderer.domElement.getBoundingClientRect();
  return {x: (v.x + 1) / 2 * r.width, y: (1 - v.y) / 2 * r.height};
 }
}
