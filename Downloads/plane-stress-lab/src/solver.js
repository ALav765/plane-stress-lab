/* Plane Stress Lab: 2D linear-elastic finite element solver.
   4-node quadrilaterals, plane stress or plane strain, matrix-free
   Jacobi-preconditioned conjugate gradients. No dependencies. */
function rrHit(c,s,a,b,rc){
  if(s*a<=c*(b-rc))return a/c;
  if(c*b<=s*(a-rc))return b/s;
  const cx=a-rc,cy=b-rc,d=c*cx+s*cy;
  return d+Math.sqrt(d*d-cx*cx-cy*cy+rc*rc);
}
function thetas(n,tc,K){
  const T=new Array(n+1),h=n/2,S=300;
  for(let seg=0;seg<2;seg++){
    const lo=seg*Math.PI/4,hi=lo+Math.PI/4,cum=[0];
    for(let k=1;k<=S;k++){const t=lo+(hi-lo)*(k-0.5)/S;cum.push(cum[k-1]+1+K*Math.exp(-(((t-tc)/0.15)**2)));}
    for(let k=0;k<=h;k++){
      const target=cum[S]*k/h;let q=1;
      while(q<S&&cum[q]<target)q++;
      const fr=(target-cum[q-1])/(cum[q]-cum[q-1]||1);
      T[seg*h+k]=lo+(hi-lo)*(q-1+fr)/S;
    }
  }
  T[0]=0;T[n]=Math.PI/2;T[n/2]=Math.PI/4;
  return T;
}
function mesh(p){
  const a=p.a,b=p.b,rc=p.rc,n=p.n,m=p.m,W=p.W,N=n+1,X=[],Y=[];
  const K=Math.max(0,8*(1-rc/Math.min(a,b))),TH=thetas(n,Math.atan2(b,a),K);
  for(let j=0;j<=m;j++){
    const t=Math.pow(j/m,1.8);
    for(let i=0;i<=n;i++){
      const th=TH[i],c=Math.cos(th),s=Math.sin(th),r=rrHit(c,s,a,b,rc);
      let xo,yo;
      if(i<=n/2){xo=W;yo=W*s/c;}else{yo=W;xo=W*c/s;}
      X.push((1-t)*r*c+t*xo);Y.push((1-t)*r*s+t*yo);
    }
  }
  return {X:X,Y:Y};
}
function bm(x,y,r,s){
  const dr=[-(1-s),(1-s),(1+s),-(1+s)],ds=[-(1-r),-(1+r),(1+r),(1-r)];
  let j11=0,j12=0,j21=0,j22=0;
  for(let k=0;k<4;k++){j11+=dr[k]*x[k]/4;j12+=dr[k]*y[k]/4;j21+=ds[k]*x[k]/4;j22+=ds[k]*y[k]/4;}
  const det=j11*j22-j12*j21,B=new Float64Array(24);
  for(let k=0;k<4;k++){
    const nx=(j22*dr[k]-j12*ds[k])/(4*det),ny=(-j21*dr[k]+j11*ds[k])/(4*det);
    B[2*k]=nx;B[8+2*k+1]=ny;B[16+2*k]=ny;B[16+2*k+1]=nx;
  }
  return {B:B,det:det};
}
function dB(D,B){
  const R=new Float64Array(24);
  for(let r=0;r<3;r++)for(let c=0;c<8;c++){let s=0;for(let k=0;k<3;k++)s+=D[r*3+k]*B[k*8+c];R[r*8+c]=s;}
  return R;
}
function keF(x,y,D){
  const K=new Float64Array(64),g=0.5773502691896258,P=[[-g,-g],[g,-g],[g,g],[-g,g]];
  for(const q of P){
    const o=bm(x,y,q[0],q[1]),DB=dB(D,o.B);
    for(let a=0;a<8;a++)for(let b=0;b<8;b++){
      let s=0;for(let c=0;c<3;c++)s+=o.B[c*8+a]*DB[c*8+b];K[a*8+b]+=s*o.det;
    }
  }
  return K;
}
function dot(a,b){let s=0;for(let i=0;i<a.length;i++)s+=a[i]*b[i];return s;}
const MATS={Steel:[200000,.3,250],Aluminium:[70000,.33,270],Titanium:[114000,.34,880],'ABS plastic':[2300,.35,40]};
function win(p){
  const n=p.n,m=p.m,N=n+1,ms=mesh(p),X=ms.X,Y=ms.Y,nd=2*X.length,els=[],F=new Float64Array(nd),fix=new Uint8Array(nd);
  for(let j=0;j<m;j++)for(let i=0;i<n;i++)els.push([j*N+i,(j+1)*N+i,(j+1)*N+i+1,j*N+i+1]);
  for(let j=0;j<=m;j++){fix[2*j*N+1]=1;fix[2*(j*N+n)]=1;}
  for(let i=0;i<n;i++){
    const a=m*N+i,b=a+1,L=Math.hypot(X[b]-X[a],Y[b]-Y[a]);
    if(i<n/2){F[2*a]+=p.lam*p.q*L/2;F[2*b]+=p.lam*p.q*L/2;}else{F[2*a+1]+=p.q*L/2;F[2*b+1]+=p.q*L/2;}
  }
  return {X:X,Y:Y,els:els,F:F,fix:fix,pk:Array.from({length:N},(_,i)=>i),mir:1};
}
function beam(p){
  const ny=p.ny,nx=2*Math.min(100,Math.round(ny*p.L/p.h/2)),N=nx+1,X=[],Y=[],els=[];
  for(let j=0;j<=ny;j++)for(let i=0;i<=nx;i++){X.push(p.L*i/nx);Y.push(p.h*(j/ny-.5));}
  for(let j=0;j<ny;j++)for(let i=0;i<nx;i++)els.push([j*N+i,j*N+i+1,(j+1)*N+i+1,(j+1)*N+i]);
  const nd=2*X.length,F=new Float64Array(nd),fix=new Uint8Array(nd);
  for(let j=0;j<=ny;j++){fix[2*j*N]=1;fix[2*j*N+1]=1;F[2*(j*N+nx)+1]=-p.P*(j%ny?1:.5)/ny;}
  const pk=[];for(let i=0;i<=nx;i++)if(X[i]>=p.h)pk.push(ny*N+i);
  return {X:X,Y:Y,els:els,F:F,fix:fix,pk:pk,N:N,nx:nx,ny:ny};
}
function cyl(p){
  const n=p.n,nr=n,N=n+1,X=[],Y=[],els=[];
  for(let j=0;j<=nr;j++){const r=p.a+(p.b-p.a)*j/nr;for(let i=0;i<=n;i++){const t=Math.PI/2*i/n;X.push(r*Math.cos(t));Y.push(r*Math.sin(t));}}
  for(let j=0;j<nr;j++)for(let i=0;i<n;i++)els.push([j*N+i,(j+1)*N+i,(j+1)*N+i+1,j*N+i+1]);
  const nd=2*X.length,F=new Float64Array(nd),fix=new Uint8Array(nd);
  for(let j=0;j<=nr;j++){fix[2*j*N+1]=1;fix[2*(j*N+n)]=1;}
  for(let i=0;i<n;i++){const L=Math.hypot(X[i+1]-X[i],Y[i+1]-Y[i]);for(const k of [i,i+1]){F[2*k]+=p.q*L/2*X[k]/p.a;F[2*k+1]+=p.q*L/2*Y[k]/p.a;}}
  return {X:X,Y:Y,els:els,F:F,fix:fix,pk:Array.from({length:N},(_,i)=>i),mir:1};
}
function fe(M,E,nu,ps){
  const X=M.X,Y=M.Y,els=M.els,nn=X.length,nd=2*nn,fix=M.fix,F=Float64Array.from(M.F);
  let D;
  if(ps){const f=E/((1+nu)*(1-2*nu));D=[f*(1-nu),f*nu,0,f*nu,f*(1-nu),0,0,0,f*(1-2*nu)/2];}
  else{const f=E/(1-nu*nu);D=[f,f*nu,0,f*nu,f,0,0,0,f*(1-nu)/2];}
  const Ks=els.map(e=>keF(e.map(k=>X[k]),e.map(k=>Y[k]),D)),dg=new Float64Array(nd);
  els.forEach((e,q)=>e.forEach((k,a)=>{dg[2*k]+=Ks[q][2*a*8+2*a];dg[2*k+1]+=Ks[q][(2*a+1)*8+2*a+1];}));
  for(let d=0;d<nd;d++)if(fix[d])F[d]=0;
  const u8=new Float64Array(8);
  function mv(x,y){
    y.fill(0);
    for(let q=0;q<els.length;q++){
      const e=els[q],K=Ks[q];
      for(let a=0;a<4;a++){u8[2*a]=x[2*e[a]];u8[2*a+1]=x[2*e[a]+1];}
      for(let r=0;r<8;r++){let s=0;for(let c=0;c<8;c++)s+=K[r*8+c]*u8[c];y[2*e[r>>1]+(r&1)]+=s;}
    }
    for(let d=0;d<nd;d++)if(fix[d])y[d]=0;
  }
  const x=new Float64Array(nd),r=Float64Array.from(F),z=new Float64Array(nd),pp=new Float64Array(nd),qv=new Float64Array(nd);
  for(let d=0;d<nd;d++)z[d]=fix[d]?0:r[d]/dg[d];
  pp.set(z);
  let rz=dot(r,z),it=0;const r0=Math.sqrt(dot(r,r));
  while(it<20000&&Math.sqrt(dot(r,r))>1e-9*r0){
    mv(pp,qv);const al=rz/dot(pp,qv);
    for(let d=0;d<nd;d++){x[d]+=al*pp[d];r[d]-=al*qv[d];z[d]=fix[d]?0:r[d]/dg[d];}
    const rz2=dot(r,z),be=rz2/rz;rz=rz2;
    for(let d=0;d<nd;d++)pp[d]=z[d]+be*pp[d];
    it++;
  }
  const A=new Float64Array(nn),B=new Float64Array(nn),C=new Float64Array(nn),ct=new Float64Array(nn),CR=[[-1,-1],[1,-1],[1,1],[-1,1]];
  els.forEach(e=>{
    const ex=e.map(k=>X[k]),ey=e.map(k=>Y[k]),ue=new Float64Array(8);
    for(let a=0;a<4;a++){ue[2*a]=x[2*e[a]];ue[2*a+1]=x[2*e[a]+1];}
    for(let a=0;a<4;a++){
      const Bm=bm(ex,ey,CR[a][0],CR[a][1]).B,ep=[0,0,0];
      for(let k=0;k<3;k++)for(let c=0;c<8;c++)ep[k]+=Bm[k*8+c]*ue[c];
      const k=e[a];
      A[k]+=D[0]*ep[0]+D[1]*ep[1];B[k]+=D[3]*ep[0]+D[4]*ep[1];C[k]+=D[8]*ep[2];ct[k]++;
    }
  });
  const vm=new Float64Array(nn),s1=new Float64Array(nn);
  for(let k=0;k<nn;k++){
    const a=A[k]/ct[k],b=B[k]/ct[k],c=C[k]/ct[k],sz=ps?nu*(a+b):0;
    A[k]=a;B[k]=b;C[k]=c;
    vm[k]=Math.sqrt(a*a+b*b+sz*sz-a*b-b*sz-sz*a+3*c*c);
    s1[k]=(a+b)/2+Math.sqrt(((a-b)/2)**2+c*c);
  }
  return {x:x,sx:A,sy:B,sxy:C,vm:vm,s1:s1,it:it};
}
function peak(M,r){let v=0,ix=0;for(const k of M.pk)if(r.vm[k]>v){v=r.vm[k];ix=k;}return {v:v,k:ix};}
const TESTS=[
 ['Kirsch plate with a hole: stress concentration',()=>{const r=fe(win({a:1,b:1,rc:1,W:20,n:32,m:32,lam:0,q:1}),1000,.3,0);return [r.sy[0],3,.02];}],
 ['Cantilever tip deflection vs. beam theory',()=>{const M=beam({L:100,h:10,P:1,ny:8}),r=fe(M,2e5,.3,0);return [-r.x[2*(M.ny/2*M.N+M.nx)+1],1e6/(3*2e5*1000/12),.06];}],
 ['Cantilever mid-span bending stress',()=>{const M=beam({L:100,h:10,P:1,ny:8}),r=fe(M,2e5,.3,0);return [r.sx[M.ny*M.N+M.nx/2],3,.05];}],
 ['Thick cylinder bore hoop stress (Lamé)',()=>{const r=fe(cyl({a:50,b:100,n:32,q:20}),2e5,.3,1);return [r.sy[0],20*12500/7500,.025];}]
];

if(typeof module!=='undefined')module.exports={mesh:mesh,win:win,beam:beam,cyl:cyl,fe:fe,peak:peak,TESTS:TESTS,MATS:MATS};
