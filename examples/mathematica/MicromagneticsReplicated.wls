#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
   SPDX-FileCopyrightText: 2025 Andreas Buchheit <buchheit@num.uni-sb.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Micromagnetics Replicated*)


(* ::Subsection:: *)
(*Setup and Auxiliary Functions*)


(* ::Text:: *)
(*Set to 0\:2070 = 1 in every kernel.*)


Unprotect[Power];
Power[0|0.|0*I|0.*I,0|0.|0*I|0.*I]:=1;
Power[z_,0]/;Abs[z]==0:=1;
Protect[Power];
ParallelEvaluate[
Unprotect[Power];
Power[0|0.|0*I|0.*I,0|0.|0*I|0.*I]:=1;
Power[z_,0]/;Abs[z]==0:=1;
Protect[Power];
];


(* ::Text:: *)
(*Function calculating the minimum of absolute and relative error.*)


Error[ex_,approx_]=Min[{Abs[ex-approx],Abs[ex-approx]/Abs[ex]}];
SetAttributes[Error,Listable]


(* ::Text:: *)
(*Load MaTeX for figures.*)


<<MaTeX`
SetOptions[MaTeX,"Preamble"->{"\\usepackage{amsmath}","\\usepackage{bm}","\\usepackage{dsfont}","\\usepackage{xcolor}"}];
SetOptions[MaTeX,FontSize->18];


(* ::Text:: *)
(*Access colors of scientific color scheme.*)


Colors=(("DefaultPlotStyle"/.(Method/. Charting`ResolvePlotTheme["Scientific",ListLinePlot]))/. Directive[x_,__]:>x);


(* ::Subsection:: *)
(*Crandall and Bessel Functions*)


(* ::Text:: *)
(*Compiled 1D Bessel function for the evaluation in the series and recursive region of the algorithm.*)


KCompiled = Compile[{
	{s, _Real}
	,{x, _Real}
	,{y, _Real}
	}
	,Module[{result=0.,num1=0.,num2=0.,num3=0.,num=0.,den1=0.,den2=0.,den3=0.,den=0.}
	,(* Series around the origin *)
		If[x+y<1.5
			,result=x^s Gamma[-s,x]
			;result+=Total[
			Table[
			x^(s+j) Gamma[-s-j,x](-y)^j/j!
			,{j,20}
			]
		]
		(* Recursive algorithm *)
		,num3 = 1.
		;den2 = Exp[x+y]
		;den3 =(x+s+1-y)den2
		;Do[
			num=((x+s+1+2(it-1)-y)num3+(2y-s-(it-1))num2-y num1)/it
			;den=((x+s+1+2(it-1)-y)den3+(2y-s-(it-1))den2-y den1)/it
			;num1 = num2
			;num2=num3
			;num3=num
			;den1 = den2
			;den2=den3
			;den3=den
			,{it,2,100,1}
		]
		;result=num/den
		]
		;Return[result]
	],CompilationTarget->"C",RuntimeOptions->"Speed"
];


(* ::Text:: *)
(*Combined incomplete Bessel function.*)


GBessel[\[Nu]_,k_,r_]:=With[{s=N[-\[Nu]/2],x=N[Pi k . k], y=N[Pi r . r]}
	,Piecewise[{
		      {
			1/s
			,x+y<10^(-40)
		},{
			Gamma[s,0,y]/y^s
			,x<10^(-40)
		},{
			x^s Gamma[-s,x]
			,y<10^(-40)
		},{
			0
			,2 Sqrt[\[Pi]]x^(-(\[Nu]+1)/4)  y^((\[Nu]-1)/4)Exp[(\[Nu]^2/16-2x y)/Sqrt[x y]]<10^-16
		},{
			2(x/y)^(s/2) BesselK[-s,2Sqrt[x y]]-KCompiled[-s,y,x]
			,x+1/5<y
		},{
			KCompiled[s,x,y]
		,True
		}
	}]
]


(* ::Text:: *)
(*Crandall functions.*)


GLower[\[Nu]_,z_]:=Piecewise[{
	{ 2/\[Nu],z . z<10^(-40)}
	,{Gamma[\[Nu]/2,0,Pi*z . z]/(Pi*z . z)^(\[Nu]/2),Norm[z]>=10^(-64)}
}];
GUpper[\[Nu]_,z_]:=Piecewise[{
	{ -2/\[Nu],Norm[z]<10^(-40)}
	,{Gamma[\[Nu]/2,Pi*z . z]/(Pi*z . z)^(\[Nu]/2),Norm[z]>=10^(-64)&&Norm[z]<4}
	,{0,True}
}];


(* ::Text:: *)
(*Derivatives of Crandall functions and the incomplete Bessel functions with respect to r.*)


p[z_,\[Alpha]_,\[Beta]_]:=Apply[Times,(-\[Pi])^(\[Alpha]-\[Beta])Binomial[\[Alpha],\[Beta]](\[Alpha]-\[Beta])!/(\[Alpha]-2\[Beta])!(2z)^(\[Alpha]-2\[Beta])]
GUpperDer[\[Nu]_,r_,\[Alpha]_]:=Total[Flatten[Table[
	p[r,\[Alpha],\[Beta]]GUpper[\[Nu]+2Total[\[Alpha]-\[Beta]],r]
	,{\[Beta],Tuples[Range[ConstantArray[0,Length[\[Alpha]]],\[Alpha]/2]]}
]]]
GLowerDer[\[Nu]_,r_,\[Alpha]_]:=Total[Flatten[Table[
	p[r,\[Alpha],\[Beta]]GLower[\[Nu]+2Total[\[Alpha]-\[Beta]],r]
	,{\[Beta],Tuples[Range[ConstantArray[0,Length[\[Alpha]]],\[Alpha]/2]]}
]]]
GBesselDer[\[Nu]_,k_,r_,\[Alpha]_]:=Total[Flatten[Table[
	p[r,\[Alpha],\[Beta]]GBessel[\[Nu]-2Total[\[Alpha]-\[Beta]],k,r]
	,{\[Beta],Tuples[Range[ConstantArray[0,Length[\[Alpha]]],\[Alpha]/2]]}
]]]


(* ::Input:: *)
(**)


(* ::Subsection:: *)
(*Far Field Zeta Function*)


(* ::Text:: *)
(*Pad a list of tuples with zeros, such that they have dimension newDim. Used to go from the n-dim lattice to the d dimensional point set L.*)


extendTuples[tuples_List,newDim_]:=PadRight[#,newDim,0]&/@tuples;


(* ::Text:: *)
(*Far field sum over the real space.*)


sumRealFarField[\[Nu]_,A_,r_,\[Alpha]_,nNear_,sumLim_]:=
Module[{d,n,zRangeFar,zRangeNear,res}
	,d=Length[r]
	;n=Length[A]
	;(*Drop zero term*)
	zRangeFar=extendTuples[A . #&/@Select[Tuples[Range[-sumLim,sumLim],n],Max[Abs[#]]>(nNear+10^(-15))&],d]
	;zRangeNear=extendTuples[A . #&/@Tuples[Range[-nNear,nNear],n],d]
	;res=Pi^(\[Nu]/2)/Gamma[\[Nu]/2]Total[ GUpperDer[\[Nu],r-#,\[Alpha]]& /@zRangeFar]-Pi^(\[Nu]/2)/Gamma[\[Nu]/2]*Total[GLowerDer[\[Nu],r-#,\[Alpha]]& /@zRangeNear]
	;Return[res]
]


(* ::Text:: *)
(*Far field sum over the reciprocal lattice.*)


sumFourier[\[Nu]_,A_,r_,\[Alpha]_,sumLim_]:=Module[{d,n,AStar,kRange,\[Alpha]Par,\[Alpha]Orth,rPar,rOrth,res}
	,d=Length[r]
	;n=Length[A]
	;(*Projections*)
	\[Alpha]Par=Take[\[Alpha],n]
	;\[Alpha]Orth=Take[\[Alpha],{n+1,d}]
	;rPar=Take[r,n]
	;rOrth=Take[r,{n+1,d}]
	;(*Reciprocal lattice matrix*)
	AStar=Transpose[Inverse[A]]
	;kRange=AStar . #&/@Tuples[Range[-sumLim,sumLim],n]
	;res=Pi^(\[Nu]/2)/Gamma[\[Nu]/2]/Abs[Det[A]]Total[Apply[Times,(2Pi*I*#)^\[Alpha]Par]Exp[2Pi*I*# . rPar]GBesselDer[n-\[Nu],#,rOrth,\[Alpha]Orth]& /@kRange]
	;Return[res]
]


(* ::Text:: *)
(*Far field zeta function.*)


ZetaLFarField[\[Nu]_,A_,r_,\[Alpha]_,nNear_,sumLim_:4]:=sumRealFarField[\[Nu],A,r,\[Alpha],nNear,sumLim]+sumFourier[\[Nu],A,r,\[Alpha],sumLim]


(* ::Text:: *)
(*Correction term of the far field zeta function.*)


ZetaLFarFieldCorrDer[\[Nu]_,A_,r_,c_,m_,nNear_,\[Alpha]Max_:2]:=Module[{d,\[Alpha]Range,res},
	d=Length[r];
	\[Alpha]Range=Select[Tuples[Range[0,\[Alpha]Max],d],Total[#]<=\[Alpha]Max&];
	res=Total[2^d/Apply[Times,(2(#+1))!]ZetaLFarField[\[Nu],A,r,2#+m,nNear]*Apply[Times,c^(2#)]&/@\[Alpha]Range];
	Return[res];
]


(* ::Subsection:: *)
(*Direct Sum over the Short Range Double Integral*)


(* ::Text:: *)
(*The integrand functions needed in the algorithm. [m] specifies the number of derivatives taken (up to 2).*)


psi[0][r_,t_,c_]=1/2 t (((-2 E^(-((\[Pi] r^2)/t^2))+E^(-((\[Pi] (-2 c+r)^2)/t^2))+E^(-((\[Pi] (2 c+r)^2)/t^2))) t)/\[Pi]+(2 c-r) Erf[(Sqrt[\[Pi]] (2 c-r))/t]-2 r Erf[(Sqrt[\[Pi]] r)/t]+(2 c+r) Erf[(Sqrt[\[Pi]] (2 c+r))/t]);


psi[1][r_,t_,c_]=1/2 t (-2 Erf[(Sqrt[\[Pi]] r)/t]+Erf[(Sqrt[\[Pi]] (-2 c+r))/t]+Erf[(Sqrt[\[Pi]] (2 c+r))/t]);


psi[2][r_,t_,c_]=E^((8 c \[Pi] r)/t^2-(\[Pi] (2 c+r)^2)/t^2)-2 E^((4 c \[Pi] (c+r))/t^2-(\[Pi] (2 c+r)^2)/t^2)+E^(-((\[Pi] (2 c+r)^2)/t^2));


(* ::Text:: *)
(*Double integral of the direct sum.*)


SummandDoubleIntegral[\[Nu]_,r_,c_,m_]:=Module[{d,res}
	,d=Length[r]
	;res=(1/Pi)^(-\[Nu]/2)/Gamma[\[Nu]/2]/Product[c[[i]]^2,{i,1,d}]2NIntegrate[
	t^(-\[Nu]-1)*Product[psi[m[[i]]][r[[i]],t,c[[i]]/2],{i,1,d}]
	,{t,0,Infinity}
	,Exclusions->{0}
	,WorkingPrecision->20
	]
	;Return[res]
]


(* ::Text:: *)
(*Summand function in real space.*)


UNearCube[\[Nu]_,A_,r_,c_,m_,nNear_]:=Module[{d,n,zRangeNear,res}
	,d=Length[r]
	;n=Length[A]
	;zRangeNear=extendTuples[A . #&/@Tuples[Range[-nNear,nNear],n],d]
	;res=Total[ParallelTable[SummandDoubleIntegral[\[Nu],r-z,c,m],{z,zRangeNear}]]
	;Return[res]
]


(* ::Subsection:: *)
(*Potential Function*)


(* ::Text:: *)
(*Potential function for a n - dimensional lattice with lattice matrix A  in d = length (r)  dimensional  space  of  cuboids[-c_i/2, ci/2] .*)
(**)
(*Arguments*)
(*	A : lattice matrix; columns  are  primitive  lattice  vectors  (n*n).*)
(*	r : evaluation point in \:211d^d (vector of length d).*)
(*	c : side  lengths  of  the  cuboid  in  each  coordinate  (vector length d).*)
(*	c : vector describing cuboid = prod_i [-c_i/2, c_i/2]  (vector of length d).*)
(*	m : derivative  vector  (vector of length d).*)
(*	nNear : nonnegative  integer. All  lattice  points  within  A . {-nNear, \[Ellipsis], nNear}^n  are  included  by  direct  summation.*)
(*	\[Alpha]Max : maximum derivative order of the zeta-function correction term  (controls exponential convergence of the remainder).*)
(**)
(*Error  Estimate: For  distance  r  to  the  nearest  singularity  (componentwise distance r_j), the  truncation  error  behaves  like Error = O ( prod_j (c_j/r_j)^(2 \[Alpha] + 2)) with  \[Alpha] = \[Alpha]Max . Thus increasing \[Alpha]Max reduces the error exponentially in \[Alpha] .*)


UCube[\[Nu]_,A_,r_,c_,m_,nNear_,\[Alpha]Max_]:=
	Module[{\[Nu]Num,ANum,rNum,cNum,mNum,nNearNum,\[Alpha]MaxNum,res},
	{\[Nu]Num,ANum,rNum,cNum,mNum,nNearNum,\[Alpha]MaxNum}=N[{\[Nu],A,r,c,m,nNear,\[Alpha]Max}];
	res=ZetaLFarFieldCorrDer[\[Nu],A,r,c,m,nNear,\[Alpha]Max]+UNearCube[\[Nu],A,r,c,m,nNear];
	Return[res];
]


(* ::Subsection::Closed:: *)
(*Demagnetization Factor*)


(* ::Text:: *)
(*Demagnetization factor in z direction Dz.*)
(**)
(*Arguments:*)
(*	NN : number of divisions of the three-dimensional cuboids per dimension.*)
(*	nNear : Truncation term for the direct summation.*)
(*	\[Alpha]Max : maximum derivative order of the zeta-function correction term  (controls exponential convergence of the remainder).*)


(* ::Input:: *)
(*Dz[NN_,nNear_,\[Alpha]Max_]:=1-1/(4Pi)(1/NN^3)(-UCube[1,{{1,0},{0,1}},{0,0,0},{1/NN,1/NN,1/NN},{2,0,0},nNear,\[Alpha]Max]-UCube[1,{{1,0},{0,1}},{0,0,0},{1/NN,1/NN,1/NN},{0,2,0},nNear,\[Alpha]Max]);*)


(* ::Text:: *)
(*Asymptotic behavior that yields high-precision values for large NN.*)


(* ::Input:: *)
(*DzAsymptotic[NN_]:=N[1/3+DirichletBeta[3/2]Zeta[3/2]/(Pi*NN^3),60]*)


(* ::Section:: *)
(*Asymptotic Expansion Error*)


(* ::Text:: *)
(*Calculate Asymptotic behavior for a range of subdivision values NN.*)


NNRange=Range[5,100,5];
DzAsymptoticTab=ParallelTable[DzAsymptotic[NN],{NN,NNRange}];


(* ::Text:: *)
(*Calculate high-precision demagnetization values. Takes less than 40 seconds on a standard laptop with 12 kernels.*)


{evalTime,DzTab}=AbsoluteTiming[Table[Dz[NN,3,2],{NN,NNRange}]]


(* ::Text:: *)
(*Show Error of asymptotic expansion in comparison with high-precision values.*)


pDeviation=ListLogLogPlot[
	Transpose[{NNRange,Error[DzTab,DzAsymptoticTab]}]
	,PlotRange->{{4.35,115},{10^-5,2 10^-16}}
	,PlotTheme->"Scientific"
	,PlotStyle->{Colors[[2]],PointSize[0.025]}
	,GridLines->{{5,10,50,100},10^Range[-15,0,3]}
	,FrameLabel->{MaTeX["N"],MaTeX["|D_z-D_z^{\\mathrm{asym}}|"]}
	,FrameTicks->{
		{
			Join[
			{10^#,MaTeX[If[#>=-1,ToString[10^#],"10^{"<>ToString[#]<>"}"]],{0.014,0}}&/@Range[-15,0,3]
			,{#,,{0.008,0}}&/@ Flatten[Table[i 10^j,{i,{1}},{j,-20,1,1}]]
			]
			,Join[
			{10^#,,{0.014,0}}&/@Range[-16,0,2]
			,{#,,{0.008,0}}&/@ Flatten[Table[i 10^j,{i,{1}},{j,-20,1,1}]]
			]
		}
		,{
			Join[
			{#,MaTeX[ToString@#],{0.014,0}}&/@{5,10,50,100}
			,{#,,{0.008,0}}&/@Join[Range[5,10,1],Range[10,50,10],Range[50,100,20]]
			]
			,Join[
			{#,,{0.014,0}}&/@{5,10,50,100}
			,{#,,{0.008,0}}&/@Join[Range[5,10,1],Range[10,50,10],Range[50,100,20]]
			]
		}
	}
	,ImageSize->400
]
