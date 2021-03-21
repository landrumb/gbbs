(window.webpackJsonp=window.webpackJsonp||[]).push([[17],{120:function(a,e,t){"use strict";t.d(e,"a",(function(){return i})),t.d(e,"b",(function(){return h}));var n=t(0),s=t.n(n);function m(a,e,t){return e in a?Object.defineProperty(a,e,{value:t,enumerable:!0,configurable:!0,writable:!0}):a[e]=t,a}function p(a,e){var t=Object.keys(a);if(Object.getOwnPropertySymbols){var n=Object.getOwnPropertySymbols(a);e&&(n=n.filter((function(e){return Object.getOwnPropertyDescriptor(a,e).enumerable}))),t.push.apply(t,n)}return t}function r(a){for(var e=1;e<arguments.length;e++){var t=null!=arguments[e]?arguments[e]:{};e%2?p(Object(t),!0).forEach((function(e){m(a,e,t[e])})):Object.getOwnPropertyDescriptors?Object.defineProperties(a,Object.getOwnPropertyDescriptors(t)):p(Object(t)).forEach((function(e){Object.defineProperty(a,e,Object.getOwnPropertyDescriptor(t,e))}))}return a}function c(a,e){if(null==a)return{};var t,n,s=function(a,e){if(null==a)return{};var t,n,s={},m=Object.keys(a);for(n=0;n<m.length;n++)t=m[n],e.indexOf(t)>=0||(s[t]=a[t]);return s}(a,e);if(Object.getOwnPropertySymbols){var m=Object.getOwnPropertySymbols(a);for(n=0;n<m.length;n++)t=m[n],e.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(a,t)&&(s[t]=a[t])}return s}var b=s.a.createContext({}),l=function(a){var e=s.a.useContext(b),t=e;return a&&(t="function"==typeof a?a(e):r(r({},e),a)),t},i=function(a){var e=l(a.components);return s.a.createElement(b.Provider,{value:e},a.children)},o={inlineCode:"code",wrapper:function(a){var e=a.children;return s.a.createElement(s.a.Fragment,{},e)}},N=s.a.forwardRef((function(a,e){var t=a.components,n=a.mdxType,m=a.originalType,p=a.parentName,b=c(a,["components","mdxType","originalType","parentName"]),i=l(t),N=n,h=i["".concat(p,".").concat(N)]||i[N]||o[N]||m;return t?s.a.createElement(h,r(r({ref:e},b),{},{components:t})):s.a.createElement(h,r({ref:e},b))}));function h(a,e){var t=arguments,n=e&&e.mdxType;if("string"==typeof a||n){var m=t.length,p=new Array(m);p[0]=N;var r={};for(var c in e)hasOwnProperty.call(e,c)&&(r[c]=e[c]);r.originalType=a,r.mdxType="string"==typeof a?a:n,p[1]=r;for(var b=2;b<m;b++)p[b]=t[b];return s.a.createElement.apply(null,p)}return s.a.createElement.apply(null,t)}N.displayName="MDXCreateElement"},87:function(a,e,t){"use strict";t.r(e),t.d(e,"frontMatter",(function(){return p})),t.d(e,"metadata",(function(){return r})),t.d(e,"toc",(function(){return c})),t.d(e,"default",(function(){return l}));var n=t(3),s=t(7),m=(t(0),t(120)),p={id:"general_weight_sssp",title:"General-Weight SSSP (Bellman-Ford)"},r={unversionedId:"benchmarks/sssp/general_weight_sssp",id:"benchmarks/sssp/general_weight_sssp",isDocsHomePage:!1,title:"General-Weight SSSP (Bellman-Ford)",description:"Problem Specification",source:"@site/docs/benchmarks/sssp/general_weight_sssp.md",slug:"/benchmarks/sssp/general_weight_sssp",permalink:"/gbbs/docs/benchmarks/sssp/general_weight_sssp",version:"current",sidebar:"docs",previous:{title:"Positive-Weight SSSP (Delta Stepping)",permalink:"/gbbs/docs/benchmarks/sssp/positive_weight_sssp"},next:{title:"Single-Source Widest Path",permalink:"/gbbs/docs/benchmarks/sssp/ss_widest_path"}},c=[{value:"Problem Specification",id:"problem-specification",children:[]},{value:"Algorithm Implementations",id:"algorithm-implementations",children:[]},{value:"Cost Bounds",id:"cost-bounds",children:[]},{value:"Compiling and Running",id:"compiling-and-running",children:[]},{value:"References",id:"references",children:[]}],b={toc:c};function l(a){var e=a.components,t=Object(s.a)(a,["components"]);return Object(m.b)("wrapper",Object(n.a)({},b,t,{components:e,mdxType:"MDXLayout"}),Object(m.b)("h2",{id:"problem-specification"},"Problem Specification"),Object(m.b)("h4",{id:"input"},"Input"),Object(m.b)("p",null,Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"G"),Object(m.b)("mo",{parentName:"mrow"},"="),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(m.b)("mi",{parentName:"mrow"},"V"),Object(m.b)("mo",{parentName:"mrow",separator:"true"},","),Object(m.b)("mi",{parentName:"mrow"},"E"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},")")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"G=(V, E)")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"G"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}}),Object(m.b)("span",{parentName:"span",className:"mrel"},"="),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}})),Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(m.b)("span",{parentName:"span",className:"mopen"},"("),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.22222em"}},"V"),Object(m.b)("span",{parentName:"span",className:"mpunct"},","),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.16666666666666666em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.05764em"}},"E"),Object(m.b)("span",{parentName:"span",className:"mclose"},")"))))),", an unweighted graph, and a source, ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"s"),Object(m.b)("mo",{parentName:"mrow"},"\u2208"),Object(m.b)("mi",{parentName:"mrow"},"V")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"s \\in V")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.5782em",verticalAlign:"-0.0391em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"s"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}}),Object(m.b)("span",{parentName:"span",className:"mrel"},"\u2208"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}})),Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.22222em"}},"V"))))),". The input\ngraph can either be undirected or directed."),Object(m.b)("h4",{id:"output"},"Output"),Object(m.b)("p",null,"Output: ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"D")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"D")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.02778em"}},"D"))))),", a mapping where ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"D"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"["),Object(m.b)("mi",{parentName:"mrow"},"v"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"]")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"D[v]")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.02778em"}},"D"),Object(m.b)("span",{parentName:"span",className:"mopen"},"["),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.03588em"}},"v"),Object(m.b)("span",{parentName:"span",className:"mclose"},"]")))))," is the shortest path distance from\n",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"s")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"s")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.43056em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"s")))))," to ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"v")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"v")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.43056em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.03588em"}},"v")))))," in ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"G")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"G")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"G")))))," and ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow",mathvariant:"normal"},"\u221e")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"\\infty")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.43056em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord"},"\u221e")))))," if ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"v")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"v")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.43056em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.03588em"}},"v")))))," is unreachable.  If the graph\ncontains any negative-weight cycles reachable from ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"s")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"s")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.43056em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"s"))))),", the vertices\nof these negative-weight cycles and vertices reachable from them must\nhave a distance of ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mo",{parentName:"mrow"},"\u2212"),Object(m.b)("mi",{parentName:"mrow",mathvariant:"normal"},"\u221e")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"-\\infty")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.66666em",verticalAlign:"-0.08333em"}}),Object(m.b)("span",{parentName:"span",className:"mord"},"\u2212"),Object(m.b)("span",{parentName:"span",className:"mord"},"\u221e"))))),"."),Object(m.b)("h2",{id:"algorithm-implementations"},"Algorithm Implementations"),Object(m.b)("p",null,"The code for our implemenation is available\n",Object(m.b)("a",{parentName:"p",href:"https://github.com/ldhulipala/gbbs/tree/master/benchmarks/GeneralWeightSSSP/BellmanFord"},"here"),".\nWe provide more details about our implementation in ","[1]","."),Object(m.b)("h2",{id:"cost-bounds"},"Cost Bounds"),Object(m.b)("p",null,"The algorithm runs in ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"O"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(m.b)("mrow",{parentName:"mrow"},Object(m.b)("mi",{parentName:"mrow",mathvariant:"sans-serif"},"D"),Object(m.b)("mi",{parentName:"mrow",mathvariant:"sans-serif"},"i"),Object(m.b)("mi",{parentName:"mrow",mathvariant:"sans-serif"},"a"),Object(m.b)("mi",{parentName:"mrow",mathvariant:"sans-serif"},"m")),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(m.b)("mi",{parentName:"mrow"},"G"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},")"),Object(m.b)("mi",{parentName:"mrow"},"m"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},")")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"O(\\mathsf{Diam}(G)m)")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.02778em"}},"O"),Object(m.b)("span",{parentName:"span",className:"mopen"},"("),Object(m.b)("span",{parentName:"span",className:"mord"},Object(m.b)("span",{parentName:"span",className:"mord mathsf"},"D"),Object(m.b)("span",{parentName:"span",className:"mord mathsf"},"i"),Object(m.b)("span",{parentName:"span",className:"mord mathsf"},"a"),Object(m.b)("span",{parentName:"span",className:"mord mathsf"},"m")),Object(m.b)("span",{parentName:"span",className:"mopen"},"("),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"G"),Object(m.b)("span",{parentName:"span",className:"mclose"},")"),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"m"),Object(m.b)("span",{parentName:"span",className:"mclose"},")")))))," work and ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"O"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(m.b)("mrow",{parentName:"mrow"},Object(m.b)("mi",{parentName:"mrow",mathvariant:"sans-serif"},"D"),Object(m.b)("mi",{parentName:"mrow",mathvariant:"sans-serif"},"i"),Object(m.b)("mi",{parentName:"mrow",mathvariant:"sans-serif"},"a"),Object(m.b)("mi",{parentName:"mrow",mathvariant:"sans-serif"},"m")),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(m.b)("mi",{parentName:"mrow"},"G"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},")"),Object(m.b)("mi",{parentName:"mrow"},"log"),Object(m.b)("mo",{parentName:"mrow"},"\u2061"),Object(m.b)("mi",{parentName:"mrow"},"n"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},")")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"O(\\mathsf{Diam}(G) \\log n)")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.02778em"}},"O"),Object(m.b)("span",{parentName:"span",className:"mopen"},"("),Object(m.b)("span",{parentName:"span",className:"mord"},Object(m.b)("span",{parentName:"span",className:"mord mathsf"},"D"),Object(m.b)("span",{parentName:"span",className:"mord mathsf"},"i"),Object(m.b)("span",{parentName:"span",className:"mord mathsf"},"a"),Object(m.b)("span",{parentName:"span",className:"mord mathsf"},"m")),Object(m.b)("span",{parentName:"span",className:"mopen"},"("),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"G"),Object(m.b)("span",{parentName:"span",className:"mclose"},")"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.16666666666666666em"}}),Object(m.b)("span",{parentName:"span",className:"mop"},"lo",Object(m.b)("span",{parentName:"span",style:{marginRight:"0.01389em"}},"g")),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.16666666666666666em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"n"),Object(m.b)("span",{parentName:"span",className:"mclose"},")"))))),"\ndepth. Please ","[1]"," for details."),Object(m.b)("h2",{id:"compiling-and-running"},"Compiling and Running"),Object(m.b)("p",null,"The benchmark can be compiled by running:"),Object(m.b)("pre",null,Object(m.b)("code",{parentName:"pre"},"bazel build -c opt //benchmarks/GeneralWeightSSSP/BellmanFord:BellmanFord_main\n")),Object(m.b)("p",null,"It can then be run on an input graph in the ",Object(m.b)("em",{parentName:"p"},"uncompressed format")," as follows:"),Object(m.b)("pre",null,Object(m.b)("code",{parentName:"pre"},"numactl -i all ./bazel-bin/benchmarks/GeneralWeightSSSP/BellmanFord/BellmanFord_main -s -m -src 1 inputs/rMatGraph_J_5_100\n")),Object(m.b)("p",null,"It can then be run on an input graph in the ",Object(m.b)("em",{parentName:"p"},"compressed format")," as follows:"),Object(m.b)("pre",null,Object(m.b)("code",{parentName:"pre"},"numactl -i all ./bazel-bin/benchmarks/GeneralWeightSSSP/BellmanFord/BellmanFord_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda\n")),Object(m.b)("h2",{id:"references"},"References"),Object(m.b)("p",null,"[1]"," Laxman Dhulipala, Guy Blelloch, and Julian Shun",Object(m.b)("br",null),"\n",Object(m.b)("a",{parentName:"p",href:"https://ldhulipala.github.io/papers/gbbs_topc.pdf"},Object(m.b)("em",{parentName:"a"},"Theoretically Efficient Parallel Graph Algorithms Can Be Fast and Scalable")),Object(m.b)("br",null),"\nProceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 393-404, 2018. ",Object(m.b)("br",null)))}l.isMDXComponent=!0}}]);