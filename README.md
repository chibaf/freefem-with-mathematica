# freefem-with-mathematica

<HTML>

<HEAD>
<TITLE>Working with Mathematica</TITLE>
<META http-equiv="Content-Type" content="text/html; charset=shift_jis">
<META NAME="keyword" lang="en" CONTENT="Working with Mathematica">
<META NAME="author" content="Fumihiro Chiba">
<META NAME="description" content="Freefem++: Working with Mathematica">
</HEAD>

<BODY><!-- BGCOLOR="#FFFFFF" background="http://homepage.mac.com/chibaf/back02.JPG">-->

<i><H1>Working with Mathematica</H1></i>
<h2>Problem: PDE</h2>
<p>We consider the following problem:</p>
<ul>
<p>(-Delta + 1)u=f in Omega, Delta=d^2/dx^2+d^2/dy^2</p>
<p>u=g on the boundary of Omega</p>
</ul>
<img width="616" height="616" alt="image" src="https://github.com/user-attachments/assets/17e8a251-5226-4b7b-ac4d-018f4e3ed85c" />

<img width="616" height="616" alt="image" src="https://github.com/user-attachments/assets/6e0361aa-6617-43b3-b137-40aa12d44ece" />

<h2>Generating mesh and matrices by Frefem++</h2>
<p>We choose the following f and g: f=-10*sin(x)sin(y), g=0.</p>
<h3>Freefem code</h3>
<pre>
//   define borders
real rad=1.0;
int n1;
func f=-10*sin(x)*sin(y);
func g=0;
n1=64;
border a(t=0,2*pi){x=rad*(1-0.5*sin(2*t))*cos(t);y=rad*(1-0.25*cos(3*t))*sin(t);label=1;};

//   define domain and mesh
mesh Omega=buildmesh(a(n1));
savemesh(Omega,"./wwm.msh");
plot(Omega,ps="./Omega.eps");

//generate fem space
fespace Vh(Omega,P1);
Vh fh=f;

//define problem
varf mat1(u,v)=
  int2d(Omega)(dx(u)*dx(v)+dy(u)*dy(v)+u*v)+on(a,u=g);
varf mat2(u,v)= int2d(Omega)(u*v);

// define matrix and functions
matrix&lt;real&gt; a1=mat1(Vh,Vh);
matrix&lt;real&gt; b1=mat2(Vh,Vh);

//  write the sparse matrix A
 {
  ofstream file("./mat1.dat");
  file &lt;&lt; a1 &lt;&lt; endl;
 }

//  write the sparse matrix B
 {
  ofstream file("./mat2.dat");
  file &lt;&lt; b1 &lt;&lt; endl;
 }

//  write the nonhomogenous trem f
{ 
ofstream file("./f.dat"); 
file &lt;&lt; fh[]; 
} 
</pre>

<h2>Mathematica: Reading mesh and matrix, Solving equation and writing solution</h2>
<a href="./mcode"><h3>Mathematica code</h3></a>

<h2>Displaying solution by Freefem++ and medit</h2>
<p><img src="./ss.png" align=Top alt="The Solution of Equation plotted by Freefem++" width="400" height="400"><br>The Solution of Equation plotted by Freefem++</p>
<p><img src="./s_me.png" align=Top alt="The Solution of Equation plotted by medit" width="400" height="400"><br>The Solution of Equation plotted by medit</p>
<h3>Freefem++ Code</h3>
<pre>
string dr,fn;
dr="./";

mesh Omega = readmesh(dr+"wwm.msh");
//plot(Omega,ps=dr+"ss.eps");
fespace femp1(Omega,P1);
femp1 f;

fn=dr+"s.dat";
{
ifstream file(fn);
file &gt;&gt; f[] ;
}
plot(f,fill=1,ps=dr+"ss.eps",wait=1);
savemesh(Omega,"./s",[x,y,f]);
</pre>
<br>



</BODY>

</HTML>
