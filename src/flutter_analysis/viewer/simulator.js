/* Same strip equations as HStab.rhs; usable in a Web Worker and Node checks. */
"use strict";
(() => {
  function random(seed) {
    let spare = null;
    const uniform = () => {seed = (Math.imul(1664525, seed) + 1013904223) >>> 0; return (seed + .5) / 4294967296;};
    return () => {
      if (spare !== null) {const value = spare; spare = null; return value;}
      const r = Math.sqrt(-2 * Math.log(uniform())), a = 2 * Math.PI * uniform();
      spare = r * Math.sin(a); return r * Math.cos(a);
    };
  }
  function makeRhs(kernel, V) {
    const p = kernel.params, N = p.n_strips;
    const {b, a, phi: ph, psi: ps, dy} = kernel;
    const m = kernel.inverseMass, K = kernel.K, C = kernel.C;
    const lam0 = b.map(v => .0455 * V / v), lam1 = b.map(v => .3 * V / v);
    let nc0 = 0, nc1 = 0;
    for (let i = 0; i < N; i++) {
      nc0 += Math.PI * p.rho * b[i] ** 2 * V * ph[i] * ps[i] * dy[i];
      nc1 -= Math.PI * p.rho * b[i] ** 3 * V * (.5 - a[i]) * ps[i] ** 2 * dy[i];
    }
    return (s, gust, out) => {
      let fh = nc0 * s[3] - K[0] * s[0] - C[0] * s[2];
      let ft = nc1 * s[3] - K[1] * s[1] - C[1] * s[3];
      for (let i = 0; i < N; i++) {
        const Q = V * ps[i] * s[1] - ph[i] * s[2] + b[i] * (.5 - a[i]) * ps[i] * s[3] + gust;
        const cq = .5 * Q + .165 * lam0[i] * s[4+i] + .335 * lam1[i] * s[4+N+i];
        const alpha = p.alpha0_deg * Math.PI / 180 + cq / V;
        const lift = p.rho * V * V * b[i] * p.cl_max * Math.tanh(p.cl_alpha * alpha / p.cl_max);
        fh += lift * ph[i] * dy[i];
        ft += lift * b[i] * (a[i] + .5) * ps[i] * dy[i];
        out[4+i] = -lam0[i] * s[4+i] + Q;
        out[4+N+i] = -lam1[i] * s[4+N+i] + Q;
      }
      out[0] = s[2]; out[1] = s[3];
      out[2] = m[0][0] * fh + m[0][1] * ft;
      out[3] = m[1][0] * fh + m[1][1] * ft;
      return out;
    };
  }
  function simulate(kernel, options) {
    const {speed: V, duration, excitation: e} = options;
    const fs = options.fs || 600, substeps = options.substeps || 4, dt = 1 / (fs * substeps);
    if (!(V > 0 && duration > 0 && Number.isFinite(V) && Number.isFinite(duration))) throw Error("Invalid speed or duration");
    const steps = Math.max(2, Math.ceil(duration / dt) + 1), n = 4 + 2 * kernel.params.n_strips;
    const state = [new Float64Array(n), new Float64Array(n)];
    if (e.kind === "release") {
      for (const s of state) {s[0] = .004; s[1] = .5 * Math.PI / 180;}
      if (e.pattern === "antisymmetric") for (let j = 0; j < n; j++) state[0][j] *= -1;
    }
    const normal = random(e.seed), process = [normal(), normal(), normal()];
    const decay = Math.exp(-V * dt / e.length_scale), innovation = Math.sqrt(1 - decay * decay);
    const common = Math.sqrt(e.correlation), local = Math.sqrt(1 - e.correlation);
    function gustAt(i) {
      if (i) for (let j = 0; j < 3; j++) process[j] = decay * process[j] + innovation * normal();
      let left = 0, right = 0;
      if (e.kind === "gusts" || e.kind === "combined") {
        if (e.pattern === "symmetric") left = right = e.gust_rms * process[0];
        else if (e.pattern === "antisymmetric") {right = e.gust_rms * process[0]; left = -right;}
        else {left = e.gust_rms * (common * process[0] + local * process[1]); right = e.gust_rms * (common * process[0] + local * process[2]);}
      }
      if (e.kind === "periodic" || e.kind === "combined") {
        const tone = e.tone_amplitude * Math.sin(2 * Math.PI * e.tone_frequency * i * dt);
        right += tone; left += e.pattern === "antisymmetric" ? -tone : tone;
      }
      return [left, right];
    }
    const result = {speed:V, time:[],leftHeave:[],leftTwist:[],rightHeave:[],rightTwist:[],gust:[],stopped:false,excitation:e};
    const rhs = makeRhs(kernel, V), work = Array.from({length:5}, () => new Float64Array(n));
    const [k1,k2,k3,k4,tmp] = work;
    let gust = gustAt(0);
    for (let i = 0; i < steps; i++) {
      if (i % substeps === 0) {
        result.time.push(i * dt); result.leftHeave.push(state[0][0]); result.leftTwist.push(state[0][1]);
        result.rightHeave.push(state[1][0]); result.rightTwist.push(state[1][1]); result.gust.push(gust);
        if (state.some(s => Math.abs(s[0]) > .1 * kernel.params.semispan || Math.abs(s[1]) > Math.PI/18)) {result.stopped = true; break;}
      }
      if (i + 1 === steps) break;
      const nextGust = gustAt(i+1);
      for (let side = 0; side < 2; side++) {
        const s = state[side], g0 = gust[side], g1 = nextGust[side], gm = .5 * (g0+g1);
        rhs(s,g0,k1);
        for (let j=0;j<n;j++) tmp[j]=s[j]+dt/2*k1[j]; rhs(tmp,gm,k2);
        for (let j=0;j<n;j++) tmp[j]=s[j]+dt/2*k2[j]; rhs(tmp,gm,k3);
        for (let j=0;j<n;j++) tmp[j]=s[j]+dt*k3[j]; rhs(tmp,g1,k4);
        for (let j=0;j<n;j++) {s[j]+=dt/6*(k1[j]+2*k2[j]+2*k3[j]+k4[j]); if(!Number.isFinite(s[j])) throw Error("Integration became nonfinite");}
      }
      gust = nextGust;
    }
    return result;
  }
  globalThis.TailSimulation = {makeRhs,simulate};
  if (typeof WorkerGlobalScope !== "undefined" && self instanceof WorkerGlobalScope) {
    self.onmessage = ({data}) => {
      try {self.postMessage({id:data.id,result:simulate(data.kernel,data.options)});}
      catch (error) {self.postMessage({id:data.id,error:error.message});}
    };
  }
})();
