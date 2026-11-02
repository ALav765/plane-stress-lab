// Runs the verification suite against exact solutions. Exit code 1 on failure.
const { TESTS } = require('../src/solver.js');
let failed = 0;
for (const [name, run] of TESTS) {
  const [fem, exact, tol] = run();
  const err = fem / exact - 1, ok = Math.abs(err) < tol;
  if (!ok) failed++;
  console.log((ok ? 'PASS' : 'FAIL') + '  ' + name + '  FEM ' + fem.toPrecision(4) + '  exact ' + exact.toPrecision(4) + '  error ' + (err * 100).toFixed(2) + '%  (tol ' + tol * 100 + '%)');
}
process.exit(failed ? 1 : 0);
