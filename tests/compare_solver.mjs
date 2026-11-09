/* Called by unittest: independently verify the browser equations against Python. */
import "../src/flutter_analysis/viewer/simulator.js";
let input="";
for await (const chunk of process.stdin) input+=chunk;
const data=JSON.parse(input);
const rhs=new Float64Array(data.state.length);
TailSimulation.makeRhs(data.kernel,17.3)(data.state,.27,rhs);
process.stdout.write(JSON.stringify({rhs:Array.from(rhs),results:data.options.map(o=>TailSimulation.simulate(data.kernel,o))}));
