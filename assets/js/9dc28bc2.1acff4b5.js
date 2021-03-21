(window.webpackJsonp=window.webpackJsonp||[]).push([[24],{120:function(e,a,t){"use strict";t.d(a,"a",(function(){return l})),t.d(a,"b",(function(){return O}));var n=t(0),m=t.n(n);function s(e,a,t){return a in e?Object.defineProperty(e,a,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[a]=t,e}function p(e,a){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var n=Object.getOwnPropertySymbols(e);a&&(n=n.filter((function(a){return Object.getOwnPropertyDescriptor(e,a).enumerable}))),t.push.apply(t,n)}return t}function c(e){for(var a=1;a<arguments.length;a++){var t=null!=arguments[a]?arguments[a]:{};a%2?p(Object(t),!0).forEach((function(a){s(e,a,t[a])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):p(Object(t)).forEach((function(a){Object.defineProperty(e,a,Object.getOwnPropertyDescriptor(t,a))}))}return e}function r(e,a){if(null==e)return{};var t,n,m=function(e,a){if(null==e)return{};var t,n,m={},s=Object.keys(e);for(n=0;n<s.length;n++)t=s[n],a.indexOf(t)>=0||(m[t]=e[t]);return m}(e,a);if(Object.getOwnPropertySymbols){var s=Object.getOwnPropertySymbols(e);for(n=0;n<s.length;n++)t=s[n],a.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(m[t]=e[t])}return m}var b=m.a.createContext({}),i=function(e){var a=m.a.useContext(b),t=a;return e&&(t="function"==typeof e?e(a):c(c({},a),e)),t},l=function(e){var a=i(e.components);return m.a.createElement(b.Provider,{value:a},e.children)},o={inlineCode:"code",wrapper:function(e){var a=e.children;return m.a.createElement(m.a.Fragment,{},a)}},N=m.a.forwardRef((function(e,a){var t=e.components,n=e.mdxType,s=e.originalType,p=e.parentName,b=r(e,["components","mdxType","originalType","parentName"]),l=i(t),N=n,O=l["".concat(p,".").concat(N)]||l[N]||o[N]||s;return t?m.a.createElement(O,c(c({ref:a},b),{},{components:t})):m.a.createElement(O,c({ref:a},b))}));function O(e,a){var t=arguments,n=a&&a.mdxType;if("string"==typeof e||n){var s=t.length,p=new Array(s);p[0]=N;var c={};for(var r in a)hasOwnProperty.call(a,r)&&(c[r]=a[r]);c.originalType=e,c.mdxType="string"==typeof e?e:n,p[1]=c;for(var b=2;b<s;b++)p[b]=t[b];return m.a.createElement.apply(null,p)}return m.a.createElement.apply(null,t)}N.displayName="MDXCreateElement"},94:function(e,a,t){"use strict";t.r(a),t.d(a,"frontMatter",(function(){return p})),t.d(a,"metadata",(function(){return c})),t.d(a,"toc",(function(){return r})),t.d(a,"default",(function(){return i}));var n=t(3),m=t(7),s=(t(0),t(120)),p={id:"connectivity",title:"Connectivity"},c={unversionedId:"benchmarks/connectivity/connectivity",id:"benchmarks/connectivity/connectivity",isDocsHomePage:!1,title:"Connectivity",description:"Problem Specification",source:"@site/docs/benchmarks/connectivity/connectivity.md",slug:"/benchmarks/connectivity/connectivity",permalink:"/gbbs/docs/benchmarks/connectivity/connectivity",version:"current",sidebar:"docs",previous:{title:"Biconnectivity",permalink:"/gbbs/docs/benchmarks/connectivity/biconnectivity"},next:{title:"Low Diameter Decomposition",permalink:"/gbbs/docs/benchmarks/connectivity/low_diameter_decomposition"}},r=[{value:"Problem Specification",id:"problem-specification",children:[]},{value:"Algorithm Implementations",id:"algorithm-implementations",children:[]},{value:"Cost Bounds",id:"cost-bounds",children:[]},{value:"Compiling and Running",id:"compiling-and-running",children:[]},{value:"References",id:"references",children:[]}],b={toc:r};function i(e){var a=e.components,t=Object(m.a)(e,["components"]);return Object(s.b)("wrapper",Object(n.a)({},b,t,{components:a,mdxType:"MDXLayout"}),Object(s.b)("h2",{id:"problem-specification"},"Problem Specification"),Object(s.b)("h4",{id:"input"},"Input"),Object(s.b)("p",null,Object(s.b)("span",{parentName:"p",className:"math math-inline"},Object(s.b)("span",{parentName:"span",className:"katex"},Object(s.b)("span",{parentName:"span",className:"katex-mathml"},Object(s.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(s.b)("semantics",{parentName:"math"},Object(s.b)("mrow",{parentName:"semantics"},Object(s.b)("mi",{parentName:"mrow"},"G"),Object(s.b)("mo",{parentName:"mrow"},"="),Object(s.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(s.b)("mi",{parentName:"mrow"},"V"),Object(s.b)("mo",{parentName:"mrow",separator:"true"},","),Object(s.b)("mi",{parentName:"mrow"},"E"),Object(s.b)("mo",{parentName:"mrow",stretchy:"false"},")")),Object(s.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"G=(V, E)")))),Object(s.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(s.b)("span",{parentName:"span",className:"base"},Object(s.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(s.b)("span",{parentName:"span",className:"mord mathdefault"},"G"),Object(s.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}}),Object(s.b)("span",{parentName:"span",className:"mrel"},"="),Object(s.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}})),Object(s.b)("span",{parentName:"span",className:"base"},Object(s.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(s.b)("span",{parentName:"span",className:"mopen"},"("),Object(s.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.22222em"}},"V"),Object(s.b)("span",{parentName:"span",className:"mpunct"},","),Object(s.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.16666666666666666em"}}),Object(s.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.05764em"}},"E"),Object(s.b)("span",{parentName:"span",className:"mclose"},")"))))),", an undirected graph on ",Object(s.b)("span",{parentName:"p",className:"math math-inline"},Object(s.b)("span",{parentName:"span",className:"katex"},Object(s.b)("span",{parentName:"span",className:"katex-mathml"},Object(s.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(s.b)("semantics",{parentName:"math"},Object(s.b)("mrow",{parentName:"semantics"},Object(s.b)("mi",{parentName:"mrow"},"n")),Object(s.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"n")))),Object(s.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(s.b)("span",{parentName:"span",className:"base"},Object(s.b)("span",{parentName:"span",className:"strut",style:{height:"0.43056em",verticalAlign:"0em"}}),Object(s.b)("span",{parentName:"span",className:"mord mathdefault"},"n")))))," vertices."),Object(s.b)("h4",{id:"output"},"Output"),Object(s.b)("p",null,Object(s.b)("span",{parentName:"p",className:"math math-inline"},Object(s.b)("span",{parentName:"span",className:"katex"},Object(s.b)("span",{parentName:"span",className:"katex-mathml"},Object(s.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(s.b)("semantics",{parentName:"math"},Object(s.b)("mrow",{parentName:"semantics"},Object(s.b)("mi",{parentName:"mrow"},"C")),Object(s.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"C")))),Object(s.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(s.b)("span",{parentName:"span",className:"base"},Object(s.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(s.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.07153em"}},"C"))))),", a ",Object(s.b)("a",{parentName:"p",href:"/docs/benchmarks/definitions"},"mapping")," where ",Object(s.b)("span",{parentName:"p",className:"math math-inline"},Object(s.b)("span",{parentName:"span",className:"katex"},Object(s.b)("span",{parentName:"span",className:"katex-mathml"},Object(s.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(s.b)("semantics",{parentName:"math"},Object(s.b)("mrow",{parentName:"semantics"},Object(s.b)("mi",{parentName:"mrow"},"C"),Object(s.b)("mo",{parentName:"mrow",stretchy:"false"},"["),Object(s.b)("mi",{parentName:"mrow"},"v"),Object(s.b)("mo",{parentName:"mrow",stretchy:"false"},"]")),Object(s.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"C[v]")))),Object(s.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(s.b)("span",{parentName:"span",className:"base"},Object(s.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(s.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.07153em"}},"C"),Object(s.b)("span",{parentName:"span",className:"mopen"},"["),Object(s.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.03588em"}},"v"),Object(s.b)("span",{parentName:"span",className:"mclose"},"]")))))," is a unique id\nbetween ",Object(s.b)("span",{parentName:"p",className:"math math-inline"},Object(s.b)("span",{parentName:"span",className:"katex"},Object(s.b)("span",{parentName:"span",className:"katex-mathml"},Object(s.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(s.b)("semantics",{parentName:"math"},Object(s.b)("mrow",{parentName:"semantics"},Object(s.b)("mo",{parentName:"mrow",stretchy:"false"},"["),Object(s.b)("mn",{parentName:"mrow"},"0"),Object(s.b)("mo",{parentName:"mrow",separator:"true"},","),Object(s.b)("mi",{parentName:"mrow"},"n"),Object(s.b)("mo",{parentName:"mrow",stretchy:"false"},")")),Object(s.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"[0, n)")))),Object(s.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(s.b)("span",{parentName:"span",className:"base"},Object(s.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(s.b)("span",{parentName:"span",className:"mopen"},"["),Object(s.b)("span",{parentName:"span",className:"mord"},"0"),Object(s.b)("span",{parentName:"span",className:"mpunct"},","),Object(s.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.16666666666666666em"}}),Object(s.b)("span",{parentName:"span",className:"mord mathdefault"},"n"),Object(s.b)("span",{parentName:"span",className:"mclose"},")")))))," representing the component of ",Object(s.b)("span",{parentName:"p",className:"math math-inline"},Object(s.b)("span",{parentName:"span",className:"katex"},Object(s.b)("span",{parentName:"span",className:"katex-mathml"},Object(s.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(s.b)("semantics",{parentName:"math"},Object(s.b)("mrow",{parentName:"semantics"},Object(s.b)("mi",{parentName:"mrow"},"v")),Object(s.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"v")))),Object(s.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(s.b)("span",{parentName:"span",className:"base"},Object(s.b)("span",{parentName:"span",className:"strut",style:{height:"0.43056em",verticalAlign:"0em"}}),Object(s.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.03588em"}},"v")))))," s.t. ",Object(s.b)("span",{parentName:"p",className:"math math-inline"},Object(s.b)("span",{parentName:"span",className:"katex"},Object(s.b)("span",{parentName:"span",className:"katex-mathml"},Object(s.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(s.b)("semantics",{parentName:"math"},Object(s.b)("mrow",{parentName:"semantics"},Object(s.b)("mi",{parentName:"mrow"},"C"),Object(s.b)("mo",{parentName:"mrow",stretchy:"false"},"["),Object(s.b)("mi",{parentName:"mrow"},"u"),Object(s.b)("mo",{parentName:"mrow",stretchy:"false"},"]"),Object(s.b)("mo",{parentName:"mrow"},"="),Object(s.b)("mi",{parentName:"mrow"},"C"),Object(s.b)("mo",{parentName:"mrow",stretchy:"false"},"["),Object(s.b)("mi",{parentName:"mrow"},"v"),Object(s.b)("mo",{parentName:"mrow",stretchy:"false"},"]")),Object(s.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"C[u] = C[v]")))),Object(s.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(s.b)("span",{parentName:"span",className:"base"},Object(s.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(s.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.07153em"}},"C"),Object(s.b)("span",{parentName:"span",className:"mopen"},"["),Object(s.b)("span",{parentName:"span",className:"mord mathdefault"},"u"),Object(s.b)("span",{parentName:"span",className:"mclose"},"]"),Object(s.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}}),Object(s.b)("span",{parentName:"span",className:"mrel"},"="),Object(s.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}})),Object(s.b)("span",{parentName:"span",className:"base"},Object(s.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(s.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.07153em"}},"C"),Object(s.b)("span",{parentName:"span",className:"mopen"},"["),Object(s.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.03588em"}},"v"),Object(s.b)("span",{parentName:"span",className:"mclose"},"]"))))),"\nif and only if ",Object(s.b)("span",{parentName:"p",className:"math math-inline"},Object(s.b)("span",{parentName:"span",className:"katex"},Object(s.b)("span",{parentName:"span",className:"katex-mathml"},Object(s.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(s.b)("semantics",{parentName:"math"},Object(s.b)("mrow",{parentName:"semantics"},Object(s.b)("mi",{parentName:"mrow"},"u")),Object(s.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"u")))),Object(s.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(s.b)("span",{parentName:"span",className:"base"},Object(s.b)("span",{parentName:"span",className:"strut",style:{height:"0.43056em",verticalAlign:"0em"}}),Object(s.b)("span",{parentName:"span",className:"mord mathdefault"},"u")))))," and ",Object(s.b)("span",{parentName:"p",className:"math math-inline"},Object(s.b)("span",{parentName:"span",className:"katex"},Object(s.b)("span",{parentName:"span",className:"katex-mathml"},Object(s.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(s.b)("semantics",{parentName:"math"},Object(s.b)("mrow",{parentName:"semantics"},Object(s.b)("mi",{parentName:"mrow"},"v")),Object(s.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"v")))),Object(s.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(s.b)("span",{parentName:"span",className:"base"},Object(s.b)("span",{parentName:"span",className:"strut",style:{height:"0.43056em",verticalAlign:"0em"}}),Object(s.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.03588em"}},"v")))))," are in the same connected component in ",Object(s.b)("span",{parentName:"p",className:"math math-inline"},Object(s.b)("span",{parentName:"span",className:"katex"},Object(s.b)("span",{parentName:"span",className:"katex-mathml"},Object(s.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(s.b)("semantics",{parentName:"math"},Object(s.b)("mrow",{parentName:"semantics"},Object(s.b)("mi",{parentName:"mrow"},"G")),Object(s.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"G")))),Object(s.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(s.b)("span",{parentName:"span",className:"base"},Object(s.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(s.b)("span",{parentName:"span",className:"mord mathdefault"},"G"))))),"."),Object(s.b)("h2",{id:"algorithm-implementations"},"Algorithm Implementations"),Object(s.b)("p",null,"We provide multiple implementations of connectivity in GBBS. The\nprimary implementation is based on the ",Object(s.b)("a",{parentName:"p",href:"low_diameter_decomposition"},"low-diameter\ndecomposition")," based algorithm from Shun\net al. ","[1]","."),Object(s.b)("p",null,"The code for the primary implementation is available\n",Object(s.b)("a",{parentName:"p",href:"https://github.com/ldhulipala/gbbs/tree/master/benchmarks/Connectivity/WorkEfficientSDB"},"here"),"."),Object(s.b)("h2",{id:"cost-bounds"},"Cost Bounds"),Object(s.b)("p",null,"The algorithm runs in ",Object(s.b)("span",{parentName:"p",className:"math math-inline"},Object(s.b)("span",{parentName:"span",className:"katex"},Object(s.b)("span",{parentName:"span",className:"katex-mathml"},Object(s.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(s.b)("semantics",{parentName:"math"},Object(s.b)("mrow",{parentName:"semantics"},Object(s.b)("mi",{parentName:"mrow"},"O"),Object(s.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(s.b)("mi",{parentName:"mrow"},"n"),Object(s.b)("mo",{parentName:"mrow"},"+"),Object(s.b)("mi",{parentName:"mrow"},"m"),Object(s.b)("mo",{parentName:"mrow",stretchy:"false"},")")),Object(s.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"O(n + m)")))),Object(s.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(s.b)("span",{parentName:"span",className:"base"},Object(s.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(s.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.02778em"}},"O"),Object(s.b)("span",{parentName:"span",className:"mopen"},"("),Object(s.b)("span",{parentName:"span",className:"mord mathdefault"},"n"),Object(s.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2222222222222222em"}}),Object(s.b)("span",{parentName:"span",className:"mbin"},"+"),Object(s.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2222222222222222em"}})),Object(s.b)("span",{parentName:"span",className:"base"},Object(s.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(s.b)("span",{parentName:"span",className:"mord mathdefault"},"m"),Object(s.b)("span",{parentName:"span",className:"mclose"},")")))))," expected work and ",Object(s.b)("span",{parentName:"p",className:"math math-inline"},Object(s.b)("span",{parentName:"span",className:"katex"},Object(s.b)("span",{parentName:"span",className:"katex-mathml"},Object(s.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(s.b)("semantics",{parentName:"math"},Object(s.b)("mrow",{parentName:"semantics"},Object(s.b)("mi",{parentName:"mrow"},"O"),Object(s.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(s.b)("msup",{parentName:"mrow"},Object(s.b)("mo",{parentName:"msup"},Object(s.b)("mi",{parentName:"mo"},"log"),Object(s.b)("mo",{parentName:"mo"},"\u2061")),Object(s.b)("mn",{parentName:"msup"},"3")),Object(s.b)("mi",{parentName:"mrow"},"n"),Object(s.b)("mo",{parentName:"mrow",stretchy:"false"},")")),Object(s.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"O(\\log^{3} n)")))),Object(s.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(s.b)("span",{parentName:"span",className:"base"},Object(s.b)("span",{parentName:"span",className:"strut",style:{height:"1.148448em",verticalAlign:"-0.25em"}}),Object(s.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.02778em"}},"O"),Object(s.b)("span",{parentName:"span",className:"mopen"},"("),Object(s.b)("span",{parentName:"span",className:"mop"},Object(s.b)("span",{parentName:"span",className:"mop"},"lo",Object(s.b)("span",{parentName:"span",style:{marginRight:"0.01389em"}},"g")),Object(s.b)("span",{parentName:"span",className:"msupsub"},Object(s.b)("span",{parentName:"span",className:"vlist-t"},Object(s.b)("span",{parentName:"span",className:"vlist-r"},Object(s.b)("span",{parentName:"span",className:"vlist",style:{height:"0.8984479999999999em"}},Object(s.b)("span",{parentName:"span",style:{top:"-3.1473400000000002em",marginRight:"0.05em"}},Object(s.b)("span",{parentName:"span",className:"pstrut",style:{height:"2.7em"}}),Object(s.b)("span",{parentName:"span",className:"sizing reset-size6 size3 mtight"},Object(s.b)("span",{parentName:"span",className:"mord mtight"},Object(s.b)("span",{parentName:"span",className:"mord mtight"},"3"))))))))),Object(s.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.16666666666666666em"}}),Object(s.b)("span",{parentName:"span",className:"mord mathdefault"},"n"),Object(s.b)("span",{parentName:"span",className:"mclose"},")"))))),"\ndepth w.h.p., and the proof can be found in the Shun et al. paper ","[1]","."),Object(s.b)("h2",{id:"compiling-and-running"},"Compiling and Running"),Object(s.b)("p",null,"The benchmark can be compiled by running:"),Object(s.b)("pre",null,Object(s.b)("code",{parentName:"pre"},"bazel build -c opt //benchmarks/Connectivity/WorkEfficientSDB14/...\n")),Object(s.b)("p",null,"It can then be run on a test input graph in the ",Object(s.b)("em",{parentName:"p"},"uncompressed format")," as follows:"),Object(s.b)("pre",null,Object(s.b)("code",{parentName:"pre"},"numactl -i all ./bazel-bin/benchmarks/Connectivity/WorkEfficientSDB14/Connectivity_main -s -m -src 1 inputs/rMatGraph_J_5_100\n")),Object(s.b)("p",null,"It can then be run on a test input graph in the ",Object(s.b)("em",{parentName:"p"},"compressed format")," as follows:"),Object(s.b)("pre",null,Object(s.b)("code",{parentName:"pre"},"numactl -i all ./bazel-bin/benchmarks/Connectivity/WorkEfficientSDB14/Connectivity_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda\n")),Object(s.b)("h2",{id:"references"},"References"),Object(s.b)("p",null,"[1]"," Julian Shun, Laxman Dhulipala, and Guy Blelloch",Object(s.b)("br",null),"\n",Object(s.b)("em",{parentName:"p"},"A Simple and Practical Linear-Work Parallel Algorithm for Connectivity"),Object(s.b)("br",null),"\nProceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 143-153, 2014."))}i.isMDXComponent=!0}}]);