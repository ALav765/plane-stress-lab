// Inlines src/*.js into one offline-ready file: dist/plane-stress-lab.html
const fs = require('fs');
let html = fs.readFileSync('index.html', 'utf8');
html = html.replace(/<script src="(.+?)"><\/script>/g, (m, f) => '<script>\n' + fs.readFileSync(f, 'utf8') + '\n</script>');
fs.mkdirSync('dist', { recursive: true });
fs.writeFileSync('dist/plane-stress-lab.html', html);
console.log('Wrote dist/plane-stress-lab.html (' + html.length + ' bytes)');
