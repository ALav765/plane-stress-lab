const $=id=>document.getElementById(id);
const STOPS=[[0,[38,70,150]],[.35,[40,160,170]],[.6,[240,215,90]],[.8,[240,130,50]],[1,[190,30,30]]];
function cm(t){t=Math.max(0,Math.min(1,t));for(let i=1;i<STOPS.length;i++){const a=STOPS[i-1],b=STOPS[i];if(t<=b[0]){const f=(t-a[0])/(b[0]-a[0]);return 'rgb('+[0,1,2].map(k=>Math.round(a[1][k]+f*(b[1][k]-a[1][k]))).join(',')+')';}}return 'rgb(190,30,30)';}
$('leg').style.background='linear-gradient(90deg,'+STOPS.map(s=>'rgb('+s[1].join(',')+') '+s[0]*100+'%').join(',')+')';
function plot(id,xs,ys,o){
  const W=400,H=220,L=46,R=14,T=14,B=40,fx=Math.log;
  const x0=Math.min(...xs.map(fx)),x1=Math.max(...xs.map(fx)),ym=Math.max(...ys,o.ref||0)*1.12;
  const X=v=>L+(fx(v)-x0)/(x1-x0||1)*(W-L-R),Y=v=>H-B-v/ym*(H-B-T);
  let s='<svg viewBox="0 0 '+W+' '+H+'" width="100%" style="max-width:560px" role="img" aria-label="Peak stress against number of elements">';
  for(let k=0;k<=4;k++){const v=ym*k/4;s+='<line x1="'+L+'" x2="'+(W-R)+'" y1="'+Y(v)+'" y2="'+Y(v)+'" stroke="var(--ln)"/><text x="'+(L-6)+'" y="'+(Y(v)+4)+'" text-anchor="end" font-size="11" fill="var(--mu)">'+v.toFixed(0)+'</text>';}
  xs.forEach(v=>{s+='<text x="'+X(v)+'" y="'+(H-B+16)+'" text-anchor="middle" font-size="11" fill="var(--mu)">'+v+'</text>';});
  s+='<text x="'+((L+W-R)/2)+'" y="'+(H-6)+'" text-anchor="middle" font-size="11" fill="var(--mu)">number of elements</text><text x="12" y="'+(T+4)+'" font-size="11" fill="var(--mu)">peak von Mises, MPa</text>';
  s+='<line x1="'+L+'" x2="'+(W-R)+'" y1="'+Y(o.ref)+'" y2="'+Y(o.ref)+'" stroke="var(--dg)" stroke-dasharray="5 4"/><text x="'+(W-R)+'" y="'+(Y(o.ref)-5)+'" text-anchor="end" font-size="11" fill="var(--dg)">finest mesh</text>';
  s+='<polyline fill="none" stroke="var(--ac)" stroke-width="2" points="'+xs.map((v,i)=>X(v)+','+Y(ys[i])).join(' ')+'"/>';
  xs.forEach((v,i)=>{s+='<circle cx="'+X(v)+'" cy="'+Y(ys[i])+'" r="3" fill="var(--ac)"/>';});
  $(id).innerHTML=s+'</svg>';
}
const TPL={
 win:{name:'Plate with a window',ps:0,ns:[4,8,16,32,48],n:32,par:[['rc','Corner radius (window half-width = 1)',.06,1,.01,.06],['q','Hoop stress, MPa',10,150,1,53]],build:(v,n)=>win({a:1,b:1.2,rc:v.rc,W:24,n:n,m:n,lam:.5,q:v.q}),note:'Quarter model, biaxial tension with axial stress at half the hoop stress, as in a pressurised fuselage.'},
 beam:{name:'Cantilever beam',ps:0,ns:[2,4,8,16],n:8,par:[['L','Length, mm',50,200,5,100],['h','Depth, mm',5,30,1,10],['q','Tip load, N per mm of width',.5,20,.5,5]],build:(v,n)=>beam({L:v.L,h:v.h,P:v.q,ny:n}),note:'Fixed at the left, loaded at the right tip. Peak stress is read on the top surface at least one depth from the wall.'},
 cyl:{name:'Pressurised thick cylinder',ps:1,ns:[4,8,16,32],n:32,par:[['r','Outer ÷ inner radius',1.2,4,.1,2],['q','Internal pressure, MPa',1,100,1,20]],build:(v,n)=>cyl({a:50,b:50*v.r,n:n,q:v.q}),note:'Quarter of a ring, 50 mm bore, plane strain.'}
};
const S={t:'win',m:'Aluminium',ps:0,f:'vm',d:.5,v:{}};
function setT(t,v){
  S.t=t;S.ps=TPL[t].ps;S.v={};TPL[t].par.forEach(p=>S.v[p[0]]=p[5]);Object.assign(S.v,v||{});
  $('tpl').value=t;$('ps').value=S.ps;$('mat').value=S.m;
  $('ctl').innerHTML=TPL[t].par.map(p=>'<div><label><span>'+p[1]+'</span><span id="o_'+p[0]+'"></span></label><input type="range" id="s_'+p[0]+'" min="'+p[2]+'" max="'+p[3]+'" step="'+p[4]+'" value="'+S.v[p[0]]+'"></div>').join('');
  TPL[t].par.forEach(p=>$('s_'+p[0]).addEventListener('input',e=>{S.v[p[0]]=+e.target.value;run();}));
  $('note').textContent=TPL[t].note;run();
}
function run(){
  const T=TPL[S.t],mt=MATS[S.m];
  S.M=T.build(S.v,T.n);S.r=fe(S.M,mt[0],mt[1],S.ps);
  T.par.forEach(p=>$('o_'+p[0]).textContent=S.v[p[0]]);
  const pk=peak(S.M,S.r);
  $('pv').textContent=pk.v.toFixed(1);
  const sf=mt[2]/pk.v;$('sf').textContent=sf.toFixed(2);$('sf').style.color=sf<1?'var(--dg)':'';
  S.k=pk.k;draw();
}
function draw(){
  const M=S.M,r=S.r,V=r[S.f],g=$('cv').getContext('2d'),nn=M.X.length;
  let um=0;for(let k=0;k<nn;k++)um=Math.max(um,Math.hypot(r.x[2*k],r.x[2*k+1]));
  $('um').textContent=um<.01?um.toExponential(1):um.toFixed(3);
  const MX=M.mir?[[1,1],[-1,1],[1,-1],[-1,-1]]:[[1,1]],xr=Math.max(...M.X)-Math.min(...M.X);
  const k=um?S.d*.12*(S.t=='win'?4:xr)/um:0;
  $('dv').textContent=S.d?'×'+Math.round(k)+' scale':'off';
  const dx=new Float64Array(nn),dy=new Float64Array(nn);
  for(let i=0;i<nn;i++){dx[i]=M.X[i]+k*r.x[2*i];dy[i]=M.Y[i]+k*r.x[2*i+1];}
  let x0=-3.6,x1=3.6,y0=-2.4,y1=2.4;
  if(S.t!='win'){x0=y0=1e9;x1=y1=-1e9;for(let i=0;i<nn;i++)for(const m of MX){x0=Math.min(x0,m[0]*dx[i]);x1=Math.max(x1,m[0]*dx[i]);y0=Math.min(y0,m[1]*dy[i]);y1=Math.max(y1,m[1]*dy[i]);}}
  const sc=Math.min(760/(x1-x0),480/(y1-y0)),cx=(x0+x1)/2,cy=(y0+y1)/2;
  S.vw={sc:sc,cx:cx,cy:cy,dx:dx,dy:dy,MX:MX};
  const so=Array.from(V).sort((a,b)=>a-b),lo=so[Math.floor(.002*nn)],hi=so[Math.floor(.998*nn)],tx=(x,m)=>400+(m*x-cx)*sc,ty=(y,m)=>260-(m*y-cy)*sc;
  $('lo').textContent=lo.toFixed(0);$('hi').textContent=hi.toFixed(0)+' MPa';$('fl').textContent=$('fld').selectedOptions[0].text;
  g.clearRect(0,0,800,520);g.lineWidth=.6;
  for(const e of M.els){
    let v=0,out=true;for(const i of e){v+=V[i];if(tx(dx[i],1)>-20&&tx(dx[i],1)<820&&ty(dy[i],1)>-20&&ty(dy[i],1)<540)out=false;}
    if(out&&S.t=='win')continue;
    g.fillStyle=g.strokeStyle=cm((v/4-lo)/(hi-lo||1));
    for(const m of MX){g.beginPath();e.forEach((i,j)=>{const x=tx(dx[i],m[0]),y=ty(dy[i],m[1]);j?g.lineTo(x,y):g.moveTo(x,y);});g.closePath();g.fill();g.stroke();}
  }
  g.strokeStyle='#fff';g.lineWidth=2;
  for(const m of MX){g.beginPath();g.arc(tx(dx[S.k],m[0]),ty(dy[S.k],m[1]),8,0,6.2832);g.stroke();}
}
$('cv').addEventListener('click',ev=>{
  const c=$('cv'),b=c.getBoundingClientRect(),px=(ev.clientX-b.left)*800/b.width,py=(ev.clientY-b.top)*520/b.height,w=S.vw;let bd=1e9,bk=0;
  for(let i=0;i<S.M.X.length;i++)for(const m of w.MX){const d=Math.hypot(400+(m[0]*w.dx[i]-w.cx)*w.sc-px,260-(m[1]*w.dy[i]-w.cy)*w.sc-py);if(d<bd){bd=d;bk=i;}}
  const r=S.r;$('pr').textContent='Node at ('+S.M.X[bk].toFixed(2)+', '+S.M.Y[bk].toFixed(2)+'): von Mises '+r.vm[bk].toFixed(1)+', max principal '+r.s1[bk].toFixed(1)+', σx '+r.sx[bk].toFixed(1)+', σy '+r.sy[bk].toFixed(1)+', τxy '+r.sxy[bk].toFixed(1)+' MPa.';
});
$('tpl').innerHTML=Object.keys(TPL).map(k=>'<option value="'+k+'">'+TPL[k].name+'</option>').join('');
$('mat').innerHTML=Object.keys(MATS).map(k=>'<option>'+k+'</option>').join('');
$('tpl').onchange=e=>setT(e.target.value);
$('mat').onchange=e=>{S.m=e.target.value;run();};
$('ps').onchange=e=>{S.ps=+e.target.value;run();};
$('fld').onchange=e=>{S.f=e.target.value;draw();};
$('df').oninput=e=>{S.d=+e.target.value;draw();};
$('c1').onclick=()=>{S.m='Aluminium';setT('win',{rc:.06,q:53});};
$('c2').onclick=()=>{S.m='Aluminium';setT('win',{rc:.6,q:53});};
$('go').onclick=()=>{
  const T=TPL[S.t],mt=MATS[S.m],xs=[],ys=[];$('go').textContent='Running…';
  const step=i=>{
    if(i<T.ns.length){const M=T.build(S.v,T.ns[i]);xs.push(M.els.length);ys.push(peak(M,fe(M,mt[0],mt[1],S.ps)).v);setTimeout(()=>step(i+1),10);return;}
    const fin=ys[ys.length-1];plot('p2',xs,ys,{ref:fin});
    $('msg').textContent='The coarsest mesh ('+xs[0]+' elements) reports '+ys[0].toFixed(1)+' MPa, '+Math.round(Math.abs(1-ys[0]/fin)*100)+'% '+(ys[0]<fin?'below':'above')+' the finest-mesh value of '+fin.toFixed(1)+' MPa. Only trust a result once the curve has flattened.';
    $('go').textContent='Run again';
  };
  step(0);
};
setT('win');
setTimeout(()=>TESTS.forEach((t,i)=>setTimeout(()=>{
  const q=t[1](),e=q[0]/q[1]-1,ok=Math.abs(e)<q[2];
  $('tb').insertAdjacentHTML('beforeend','<tr><td>'+t[0]+'</td><td>'+q[0].toPrecision(4)+'</td><td>'+q[1].toPrecision(4)+'</td><td>'+(e*100).toFixed(1)+'%</td><td style="color:'+(ok?'var(--ac)':'var(--dg)')+'">'+(ok?'Pass':'Fail')+' (±'+q[2]*100+'%)</td></tr>');
},40*i)),50);
