# Plane Stress Lab

A browser-based finite element analyser for 2D solid mechanics, written from scratch in plain JavaScript with no dependencies. Built for a Mechanics of Solids course project.

Pick a problem, a material and a load. It shows the stress field, the deformed shape, the peak stress and a safety factor against yield, then lets you check the solver against exact solutions.

## Features

- Three parametric problems: plate with a window or hole, cantilever beam, pressurised thick cylinder
- 4-node quadrilateral elements, plane stress or plane strain, linear elastic
- Four materials (steel, aluminium, titanium, ABS plastic)
- Contours of von Mises, maximum principal, sx, sy and txy stress
- Deformed-shape view, click-to-probe, safety factor against yield
- Mesh-convergence study on any problem
- Built-in verification suite against exact solutions
- Two case studies: the de Havilland Comet window corners and the Sleipner A mesh lesson

## Run it

Open `index.html` in a browser (it needs the `src/` folder beside it). For a single offline file:

```
npm run build      # writes dist/plane-stress-lab.html
```

Run the verification suite from the command line (Node 14 or newer):

```
npm test
```

## Deploy on GitHub Pages

Push the repo, then Settings, Pages, deploy from branch `main`, folder `/ (root)`. The app is served from `index.html`.

## Verification

| Test | Exact | Solver | Error |
|---|---|---|---|
| Kirsch plate with a hole, stress concentration | 3.000 | 3.032 | 1.1% |
| Cantilever tip deflection (beam theory) | 0.0200 mm | 0.0200 mm | 0.2% |
| Cantilever mid-span bending stress | 3.000 MPa | 3.014 MPa | 0.5% |
| Thick cylinder bore hoop stress (Lame) | 33.33 MPa | 33.94 MPa | 1.8% |

The beam match is partly lucky: element stiffening in bending and shear deformation missing from beam theory are both small and partly cancel.

## Structure

```
index.html        page layout and styles
src/solver.js     meshes, element stiffness, solver, stress recovery, tests
src/ui.js         rendering, controls, probe, mesh study
test/verify.js    command-line verification runner
build.js          inlines everything into one offline file
docs/             project report
```

## How it works

1. A parametric mesh generator builds a structured quadrilateral mesh for the chosen problem, using symmetry to model a quarter where possible.
2. Each element stiffness matrix comes from 2x2 Gauss integration of B^T D B.
3. The system K u = F is solved matrix-free with Jacobi-preconditioned conjugate gradients.
4. Stresses are evaluated at element corners and averaged at the nodes. Von Mises stress compares against the material yield strength to give the safety factor.

## Limits

Linear elastic only, 2D only, fixed parametric shapes, no free-form drawing. The fixed wall of the cantilever creates an artificial stress singularity, so peak stress is read one depth away. Peak stresses at very sharp corners keep rising as the mesh is refined, so they are lower bounds. The case studies are simplified 2D lessons, not reconstructions of the real structures.

## License

MIT
