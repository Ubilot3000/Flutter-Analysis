/* Optional browser integration check. Uses installed Chrome and Node's native CDP client. */
import {spawn} from "node:child_process";
import {mkdir, mkdtemp, readFile, writeFile, rm} from "node:fs/promises";
import {tmpdir} from "node:os";
import path from "node:path";
import {pathToFileURL} from "node:url";
import assert from "node:assert/strict";

const root = path.resolve(import.meta.dirname, "..");
const output = path.join(root, "output");
await mkdir(output, {recursive:true});
// A unique local temp profile avoids OneDrive locking Chrome's port file.
const profile = await mkdtemp(path.join(tmpdir(), "flutter-viewer-"));
const chrome = process.env.CHROME_PATH || "C:/Program Files/Google/Chrome/Application/chrome.exe";
const processHandle = spawn(chrome, ["--headless=new", "--remote-debugging-port=0", `--user-data-dir=${profile}`,
  "--no-first-run", "--no-default-browser-check", "--disable-background-networking", "--enable-webgl",
  "--no-sandbox", "--use-angle=swiftshader", "--enable-unsafe-swiftshader", "about:blank"], {windowsHide: true, stdio: ["ignore", "ignore", "pipe"]});
let browserLog = "";
processHandle.stderr.on("data", chunk => {browserLog += chunk.toString();});
const delay = ms => new Promise(resolve => setTimeout(resolve, ms));
const errors = [];
let socket, id = 0, pending = new Map();
async function poll(fn, message) {
  for (let i = 0; i < 150; i++) {try {const value = await fn(); if (value) return value;} catch {} await delay(100);}
  throw Error(message);
}
function call(method, params = {}) {
  const request = ++id;
  return new Promise((resolve, reject) => {
    const timeout = setTimeout(() => {pending.delete(request); reject(Error(`CDP timeout: ${method}`));}, 15000);
    pending.set(request, {resolve: value => {clearTimeout(timeout); resolve(value);}, reject});
    socket.send(JSON.stringify({id: request, method, params}));
  });
}
async function evaluate(expression) {
  const result = await call("Runtime.evaluate", {expression, returnByValue: true, awaitPromise: true});
  if (result.exceptionDetails) throw Error(JSON.stringify(result.exceptionDetails));
  return result.result.value;
}
try {
  const pages = await poll(async () => {
    const port = (await readFile(path.join(profile, "DevToolsActivePort"), "utf8")).split("\n")[0];
    const list = await (await fetch(`http://127.0.0.1:${port}/json/list`)).json();
    return list.some(p => p.type === "page") ? list : null;
  }, "Chrome debugging endpoint did not become ready");
  socket = new WebSocket(pages.find(p => p.type === "page").webSocketDebuggerUrl);
  await new Promise((resolve, reject) => {socket.onopen = resolve; socket.onerror = reject;});
  socket.onmessage = event => {
    const message = JSON.parse(event.data);
    if (message.id) {const p = pending.get(message.id); if (p) {pending.delete(message.id); if (message.error) p.reject(Error(JSON.stringify(message.error))); else p.resolve(message.result);}}
    if (message.method === "Runtime.exceptionThrown") errors.push(message.params.exceptionDetails);
    if (message.method === "Runtime.consoleAPICalled" && message.params.type === "error") errors.push(message.params);
  };
  await call("Runtime.enable"); await call("Page.enable");
  await call("Emulation.setDeviceMetricsOverride", {width: 1400, height: 1060, deviceScaleFactor: 1, mobile: false});
  await call("Page.navigate", {url: pathToFileURL(path.join(output, "hstab_3d.html")).href});
  await poll(() => evaluate("!!window.__viewerState"), "3D view failed to initialize");
  await evaluate(`document.getElementById('play').click(); document.getElementById('time').value=.05; document.getElementById('time').dispatchEvent(new Event('input'));`);
  await delay(150);
  let state = await evaluate("window.__viewerState");
  assert.equal(state.playing, false); assert.equal(state.glError, 0); assert(state.triangles > 5000);
  const expected = await evaluate(`(() => {const d=JSON.parse(document.getElementById('model-data').textContent);const c=d.response;return {h:c.rightHeave[30],t:c.rightTwist[30]};})()`);
  assert(Math.abs(state.qh - expected.h) < 1e-10); assert(Math.abs(state.qt - expected.t) < 1e-10);
  await evaluate(`document.getElementById('play').click()`);
  await delay(200);
  assert((await evaluate("window.__viewerState")).time > .05, "Play does not advance the simulation");
  await evaluate(`document.getElementById('play').click();document.getElementById('time').value=.05;document.getElementById('time').dispatchEvent(new Event('input'));`);
  await evaluate(`document.getElementById('motion').value='heave'; document.getElementById('gain').value=9; document.getElementById('gain').dispatchEvent(new Event('input'));`);
  await delay(100); state = await evaluate("window.__viewerState");
  assert.equal(state.qt, 0); assert.equal(state.gain, 9); assert(Math.abs(state.qh - expected.h) < 1e-10);
  await evaluate(`document.getElementById('motion').value='twist'`);
  await delay(100); state = await evaluate("window.__viewerState");
  assert.equal(state.qh, 0); assert(Math.abs(state.qt - expected.t) < 1e-10);
  await evaluate(`document.getElementById('motion').value='coupled'; document.getElementById('gain').value=5; document.getElementById('gain').dispatchEvent(new Event('input')); document.getElementById('restart').click();`);
  await delay(100); assert.equal((await evaluate("window.__viewerState")).time, 0);
  async function input(id,value,event="input") {
    await evaluate(`document.getElementById('${id}').value=${JSON.stringify(value)};document.getElementById('${id}').dispatchEvent(new Event('${event}'));`);
    await poll(()=>evaluate("!window.__simulationBusy"),`Simulation did not finish after ${id}`);
    await delay(100);
  }
  async function scrub(t) {
    await evaluate(`document.getElementById('time').value=${t};document.getElementById('time').dispatchEvent(new Event('input'));`);await delay(80);
    return evaluate("window.__viewerState");
  }
  await input("velocity",17.3);
  state=await scrub(.12);assert.equal(state.speed,17.3);
  assert(Math.abs(state.left[0]-state.right[0])>1e-6,"Mixed gusts must not impose mirror motion");
  await input("pattern","symmetric","change");
  state=await scrub(.12);assert(Math.abs(state.left[0]-state.right[0])<1e-12);
  await input("pattern","antisymmetric","change");
  state=await scrub(.12);assert(Math.abs(state.left[0]+state.right[0])<1e-12);
  await input("gust-rms",0);
  state=await scrub(.12);assert.equal(state.left[0],0);assert.equal(state.right[0],0);
  await input("gust-rms",.25);
  await input("excitation","periodic","change");
  assert.equal(await evaluate("document.getElementById('tone-controls').hidden"),false);
  await input("tone-frequency",12,"change");
  state=await scrub(.12);assert(Math.abs(state.qh)>1e-6);
  await input("excitation","release","change");
  state=await scrub(0);assert(Math.abs(state.qh-.004)<1e-12);
  await input("velocity",60);
  assert.match(await evaluate("document.getElementById('status').textContent"),/Divergence/);
  assert.equal(await evaluate("window.__response.stopped"),true);
  await input("velocity",15);
  await input("excitation","gusts","change");
  await input("pattern","mixed","change");
  state=await scrub(.2);
  await evaluate(`document.querySelector('[data-view="top"]').click();document.getElementById('reset-camera').click();`);
  await delay(100);
  const layout = await call("Page.getLayoutMetrics");
  let screenshot = await call("Page.captureScreenshot", {format:"png", captureBeyondViewport: true,
    clip:{x:0,y:0,width:1400,height:Math.ceil(layout.cssContentSize.height),scale:1}});
  await writeFile(path.join(output, "viewer-preview.png"), Buffer.from(screenshot.data, "base64"));
  await call("Emulation.setDeviceMetricsOverride", {width:390, height:844, deviceScaleFactor:1, mobile:false});
  await delay(150);
  assert(await evaluate("document.documentElement.scrollWidth <= innerWidth+1"), "Mobile layout overflows horizontally");
  screenshot = await call("Page.captureScreenshot", {format:"png", captureBeyondViewport:false});
  await writeFile(path.join(output, "viewer-mobile.png"), Buffer.from(screenshot.data, "base64"));
  assert.deepEqual(errors, []);
  console.log("PASS: WebGL rendering, simulated samples, scrub/pause/restart, modal filters, gain, arbitrary airspeed slider, continuous gusts, wake excitation, independent tips, camera presets, responsive layout, and no browser errors.");
} finally {
  await writeFile(path.join(output, "browser-check.log"), browserLog);
  if (socket?.readyState === 1) {try {await call("Browser.close");} catch {} socket.close();}
  processHandle.kill();
  const resolved=path.resolve(profile),base=path.resolve(tmpdir());
  if(resolved.startsWith(base+path.sep)&&path.basename(resolved).startsWith("flutter-viewer-")) {
    try {await rm(resolved,{recursive:true,force:true,maxRetries:5,retryDelay:200});} catch {}
  }
}
